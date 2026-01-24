/** 
 * @file      jpegultrascan.cpp 
 * @brief     JPEG lossless recompressor that tries all scan possibilities to minimize size
 * @author    Aaron Kaluszka
 * @version   2.0.0
 * @date      2026-01-23
 * @copyright Copyright 2015-2026 Aaron Kaluszka
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ***********************************************/

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <limits>
#include <set>
#include <iomanip>
#include <memory>
#include <atomic>
#include <mutex>
#include <thread>
#include <chrono>
#include <boost/program_options.hpp>

#ifdef USE_MOZJPEG_API
#include "jpeglib.h"
#define jpegli_std_error jpeg_std_error
#define jpegli_create_compress jpeg_create_compress
#define jpegli_mem_dest jpeg_mem_dest
#define jpegli_copy_critical_parameters jpeg_copy_critical_parameters
#define jpegli_write_coefficients jpeg_write_coefficients
#define jpegli_finish_compress jpeg_finish_compress
#define jpegli_destroy_compress jpeg_destroy_compress
#define jpegli_write_marker jpeg_write_marker
#define jpegli_create_decompress jpeg_create_decompress
#define jpegli_mem_src jpeg_mem_src
#define jpegli_save_markers jpeg_save_markers
#define jpegli_read_header jpeg_read_header
#define jpegli_read_coefficients jpeg_read_coefficients
#define jpegli_finish_decompress jpeg_finish_decompress
#define jpegli_destroy_decompress jpeg_destroy_decompress
#else
#include "lib/jpegli/encode.h"
#include "lib/jpegli/decode.h"
#endif
#ifndef ENABLE_TURBO
#define ENABLE_TURBO 0
#endif
#ifndef ENABLE_MOZJPEG
#define ENABLE_MOZJPEG 0
#endif

#include "scan_types.h"
#if ENABLE_TURBO
#include "turbo_encode.h"
#endif
#if ENABLE_MOZJPEG
#include "mozjpeg_encode.h"
#endif

#define JPEG_APP0       0xE0    /* APP0 marker code  (JFIF) */
#define JPEG_APP14      0xEE    /* APP14 marker code  (Adobe) */
#define JPEG_COM        0xFE    /* COM marker code */


static std::string components_key(const std::vector<int>& components) {
  std::ostringstream oss;
  for (size_t i = 0; i < components.size(); ++i) {
    if (i) oss << ' ';
    oss << components[i];
  }
  return oss.str();
}

static std::string scan_key(const ScanLine& scan) {
  std::ostringstream oss;
  oss << components_key(scan.components) << ": " << scan.Ss << ' ' << scan.Se
      << ' ' << scan.Ah << ' ' << scan.Al;
  return oss.str();
}

static std::string scan_label(const ScanLine& scan) {
  std::ostringstream oss;
  oss << components_key(scan.components) << ": " << scan.Ss << ' ' << scan.Se
      << ' ' << scan.Ah << ' ' << scan.Al;
  return oss.str();
}

static std::string component_only_label(const std::vector<int>& components) {
  std::ostringstream oss;
  oss << components_key(components) << ": 0 " << (DCTSIZE2 - 1) << " 0 0";
  return oss.str();
}

struct SegmentInfo {
  size_t offset = 0;
  size_t end_offset = 0;
  int marker = 0;
};


static std::vector<SegmentInfo> collect_segments(const std::vector<uint8_t>& data) {
  std::vector<SegmentInfo> segments;
  size_t i = 0;
  while (i + 1 < data.size()) {
    if (data[i] != 0xFF) {
      ++i;
      continue;
    }
    size_t j = i;
    while (j < data.size() && data[j] == 0xFF) ++j;
    if (j >= data.size()) break;
    int marker = data[j];
    if (marker == 0x00) {
      i = j;
      continue;
    }
    if (marker == 0xD8 || marker == 0xD9) {
      segments.push_back({i, j + 1, marker});
      i = j + 1;
      continue;
    }
    if (j + 2 >= data.size()) break;
    size_t length = (static_cast<size_t>(data[j + 1]) << 8) + data[j + 2];
    segments.push_back({i, j + 1 + length, marker});
    i = j + 1 + length;
  }
  return segments;
}

static const char* marker_name(int marker) {
  switch (marker) {
    case 0xD8: return "SOI";
    case 0xD9: return "EOI";
    case 0xDB: return "DQT";
    case 0xC4: return "DHT";
    case 0xDA: return "SOS";
    case 0xC0: return "SOF0";
    case 0xC1: return "SOF1";
    case 0xC2: return "SOF2";
    case 0xC9: return "SOF9";
    case 0xFE: return "COM";
    default:
      if (marker >= 0xE0 && marker <= 0xEF) return "APP";
      return "MARK";
  }
}

static std::vector<uint8_t> read_file_bytes(const std::string& path) {
  std::ifstream file(path, std::ios::binary);
  if (!file) return {};
  file.seekg(0, std::ios::end);
  std::streamoff size = file.tellg();
  file.seekg(0, std::ios::beg);
  std::vector<uint8_t> data(size);
  if (size > 0) file.read(reinterpret_cast<char*>(data.data()), size);
  return data;
}

struct ScanSizes {
  std::vector<size_t> sizes;
  std::vector<size_t> dht_bytes;
  size_t total_dht_bytes = 0;
  size_t total_dht_markers = 0;
  size_t pre_sos_dht_bytes = 0;
};

static ScanSizes scan_sizes(const std::vector<uint8_t>& data) {
  ScanSizes result;
  size_t scanstart = 0;
  size_t nextoff = 0;
  bool seen_first_sos = false;
  for (size_t i = 0; i + 1 < data.size(); ++i) {
    if (data[i] != 0xFF) continue;
    size_t j = i;
    while (j < data.size() && data[j] == 0xFF) ++j;
    if (j >= data.size()) break;
    uint8_t type = data[j];
    if (type == 0) {
      i = j;
      continue;
    }
    if (j < nextoff + 2) {
      i = j;
      continue;
    }
    size_t tagstart = i < nextoff ? nextoff : i;
    size_t len = 0;
    if (type == 0xD8 || type == 0xD9) {
      len = 0;
    } else if (j + 2 < data.size()) {
      len = (static_cast<size_t>(data[j + 1]) << 8) + data[j + 2];
      len += j - tagstart;
    } else {
      break;
    }
    nextoff = tagstart + len;
    if (scanstart) {
      // End current scan at first marker after entropy data.
      result.sizes.back() += tagstart - scanstart;
      scanstart = 0;
      if (type == 0xD9) break;
    }

    if (type == 0xC4) {
      result.total_dht_bytes += len;
      result.total_dht_markers += 1;
      if (!seen_first_sos) {
        result.pre_sos_dht_bytes += len;
      }
    }

    if (type == 0xDA) {
      seen_first_sos = true;
      result.sizes.push_back(len);
      result.dht_bytes.push_back(0);
      scanstart = tagstart + len;
    }
    i = j;
  }
  return result;
}

static ScanLine make_scan(const std::vector<int>& components, int Ss, int Se, int Ah, int Al) {
  ScanLine scan;
  scan.components = components;
  std::sort(scan.components.begin(), scan.components.end());
  scan.components.erase(std::unique(scan.components.begin(), scan.components.end()), scan.components.end());
  scan.Ss = Ss;
  scan.Se = Se;
  scan.Ah = Ah;
  scan.Al = Al;
  return scan;
}

static jpeg_scan_info to_jpeg_scan_info(const ScanLine& scan) {
  jpeg_scan_info info{};
  info.comps_in_scan = scan.components.size();
  for (size_t i = 0; i < scan.components.size(); ++i) info.component_index[i] = scan.components[i];
  info.Ss = scan.Ss;
  info.Se = scan.Se;
  info.Ah = scan.Ah;
  info.Al = scan.Al;
  return info;
}

enum class StripMode {
  All,
  Comments,
  None,
};

static const size_t kInvalidSize = std::numeric_limits<size_t>::max() / 4;
static size_t add_sizes(size_t a, size_t b) {
  if (a >= kInvalidSize || b >= kInvalidSize) return kInvalidSize;
  if (a > kInvalidSize - b) return kInvalidSize;
  return a + b;
}

static size_t compute_scan_set_cost(const std::vector<ScanLine>& scans);

struct scan_set {
  size_t size = 0;
  std::vector<ScanLine> all_scans;
  scan_set operator+(const scan_set& r) {
    scan_set ss = *this;
    ss.all_scans.insert(std::end(ss.all_scans), std::cbegin(r.all_scans), std::cend(r.all_scans));
    ss.size = compute_scan_set_cost(ss.all_scans);
    return ss;
  }
  scan_set& operator+=(const scan_set& r) {
    all_scans.insert(std::end(all_scans), std::cbegin(r.all_scans), std::cend(r.all_scans));
    size = compute_scan_set_cost(all_scans);
    return *this;
  }
};

void queue_scans(std::vector<std::vector<ScanLine>>& scan_queue, const std::vector<ScanLine>& scans) {
  scan_queue.push_back(scans);
}

void queue_scans_multi(std::vector<std::vector<ScanLine>>& scan_queue, const std::vector<ScanLine>& scans,
                       size_t max_per_file) {
  std::vector<ScanLine> pending = scans;
  std::vector<ScanLine> dc;
  std::vector<ScanLine> ac;
  size_t i = 0;
  while (i < pending.size()) {
    if (dc.size() + ac.size() > max_per_file) {
      do {
        if (ac.empty()) break;
        ac.pop_back();
        if (i > 0) --i;
      } while (!ac.empty() && ac.back().Al != 0);
      if (!ac.empty()) {
        std::vector<ScanLine> part = dc;
        part.insert(part.end(), ac.begin(), ac.end());
        scan_queue.push_back(std::move(part));
        ac.clear();
      }
    }
    if (pending[i].Ss == 0 && pending[i].Se == 0) {
      dc.push_back(pending[i]);
      pending.erase(pending.begin() + static_cast<long>(i));
      continue;
    } else {
      ac.push_back(pending[i]);
    }
    ++i;
  }
  if (!ac.empty()) {
    std::vector<ScanLine> part = dc;
    part.insert(part.end(), ac.begin(), ac.end());
    scan_queue.push_back(std::move(part));
  }
}

std::vector<ScanLine> generate_bit_scans(const std::vector<int>& plane, int a, int b, int num_bits) {
  std::vector<ScanLine> scans;
  scans.push_back(make_scan(plane, a, b, 0, num_bits));
  for (auto i = num_bits; i > 0; --i) scans.push_back(make_scan(plane, a, b, i, i - 1));
  return scans;
}

std::vector<ScanLine> generate_bit_scans(int plane, int a, int b, int num_bits) {
  return generate_bit_scans(std::vector<int>{plane}, a, b, num_bits);
}

void generate_dc_scans(const std::vector<int>& plane, const std::vector<int>& all_components, int max_bits,
                       std::vector<std::vector<ScanLine>>& scan_queue) {
  std::vector<int> dummy_components;
  if (plane != all_components) {
    for (int c : all_components) {
      if (std::find(plane.begin(), plane.end(), c) == plane.end()) dummy_components.push_back(c);
    }
  }
  for (auto i = 0; i <= max_bits; ++i) {
    auto scans = generate_bit_scans(plane, 0, 0, i);
    if (!dummy_components.empty()) {
      scans.push_back(make_scan(dummy_components, 0, 0, 0, 9));
    }
    queue_scans(scan_queue, scans);
  }
}

void generate_ac_scans(int plane, const std::vector<int>& all_components, int max_bits,
                       std::vector<std::vector<ScanLine>>& scan_queue) {
  const size_t max_per_file = 60;
  ScanLine dummy = make_scan(all_components, 0, 0, 0, 9);
  for (auto num_coeffs = 0; num_coeffs < DCTSIZE2 - 1; ++num_coeffs) {
    auto offset_max = std::min(num_coeffs, DCTSIZE2 - 2 - num_coeffs);
    for (auto offset = 0; offset <= offset_max; ++offset) {
      for (auto num_bits = 0; num_bits <= max_bits; ++num_bits) {
        std::vector<ScanLine> scans;
        scans.push_back(dummy);
        if (offset > 0) {
          auto new_scans = generate_bit_scans(plane, 1, offset, num_bits);
          scans.insert(std::end(scans), std::cbegin(new_scans), std::cend(new_scans));
        }
        auto i = offset + 1;
        for (; i + num_coeffs < DCTSIZE2 - 1; i += num_coeffs + 1) {
          auto new_scans = generate_bit_scans(plane, i, i + num_coeffs, num_bits);
          scans.insert(std::end(scans), std::cbegin(new_scans), std::cend(new_scans));
        }
        if (i < DCTSIZE2) {
          auto new_scans = generate_bit_scans(plane, i, DCTSIZE2 - 1, num_bits);
          scans.insert(std::end(scans), std::cbegin(new_scans), std::cend(new_scans));
        }
        if (scans.size() > max_per_file) {
          queue_scans_multi(scan_queue, scans, max_per_file);
        } else {
          queue_scans(scan_queue, scans);
        }
      }
    }
  }
}

static std::unordered_map<std::string, scan_set> optimal;
static std::unordered_map<std::string, scan_set> optimal2;
static std::unordered_map<std::string, size_t> sizes;
static std::unordered_set<std::string> component_only_sizes;
static std::unordered_map<std::string, scan_set> dc_cache;
static scan_set ac_all;
static scan_set best_scans;
static int g_max_bits = 3;
static std::vector<int> g_all_components;
static std::unique_ptr<std::ofstream> g_scan_dump;
static std::unique_ptr<std::ofstream> g_scan_set_dump;
static std::unique_ptr<std::ofstream> g_symbol_dump;
static std::unique_ptr<std::ofstream> g_fail_scan_dump;
static bool g_dump_dht = false;
static bool g_dump_segments = false;
static bool g_score_dht = true;
static int g_threads = 1;
static std::atomic<uint64_t> g_lock_wait_ns{0};
static std::atomic<uint64_t> g_lock_hold_ns{0};
static std::atomic<uint64_t> g_lock_acquires{0};
static std::mutex g_optimal_mutex;
static std::mutex g_optimal2_mutex;
static std::atomic<uint64_t> g_worker_busy_ns{0};
static std::atomic<uint64_t> g_do_scans_wall_ns{0};
static std::atomic<uint64_t> g_do_scans_capacity_ns{0};
static std::atomic<uint64_t> g_do_scans_calls{0};

enum class EncoderBackend {
  Jpegli,
  Turbo,
  Mozjpeg,
};

static EncoderBackend g_encoder_backend = EncoderBackend::Jpegli;

struct DhtInfo {
  std::string summary;
  size_t total_bytes = 0;
};

struct ScanDhtInfo {
  std::string hash;
  size_t bytes = 0;
  bool valid = false;
};

struct HuffmanTable {
  int table_class = 0;
  int table_id = 0;
  std::array<uint8_t, 16> bits{};
  std::array<uint8_t, 256> huffval{};
  int symbols = 0;
};

static std::unordered_map<std::string, ScanDhtInfo> scan_dht_info;

struct JpegInput {
  jpeg_decompress_struct dinfo{};
  jpeg_error_mgr jerr{};
  std::vector<uint8_t> data;
  jvirt_barray_ptr* coef_arrays = nullptr;
};

static bool prepare_input(JpegInput& input, const std::vector<uint8_t>& data, StripMode strip_mode) {
  input.data = data;
  input.dinfo.err = jpegli_std_error(&input.jerr);
  jpegli_create_decompress(&input.dinfo);
  jpegli_mem_src(&input.dinfo, input.data.data(), input.data.size());
  if (strip_mode != StripMode::None) {
    jpegli_save_markers(&input.dinfo, JPEG_COM, 0xFFFF);
    for (int app = 0; app <= 15; ++app) {
      jpegli_save_markers(&input.dinfo, JPEG_APP0 + app, 0xFFFF);
    }
  }
  if (jpegli_read_header(&input.dinfo, TRUE) != JPEG_HEADER_OK) {
    jpegli_destroy_decompress(&input.dinfo);
    return false;
  }
  input.coef_arrays = jpegli_read_coefficients(&input.dinfo);
  if (!input.coef_arrays) {
    jpegli_destroy_decompress(&input.dinfo);
    return false;
  }
  return true;
}

static void destroy_input(JpegInput& input) {
  jpegli_finish_decompress(&input.dinfo);
  jpegli_destroy_decompress(&input.dinfo);
}

struct SymbolStats {
  std::array<uint64_t, 12> dc_symbols{};
  std::array<uint64_t, 256> ac_symbols{};
  uint64_t dc_bits = 0;
  uint64_t dc_refine_bits = 0;
  uint64_t ac_bits = 0;
  uint64_t ac_refine_bits = 0;
  uint64_t ac_new_sign_bits = 0;
  uint64_t ac_eob_symbols = 0;
  uint64_t ac_zrl_symbols = 0;
};

static constexpr int kZigZag[64] = {
  0, 1, 8, 16, 9, 2, 3, 10,
  17, 24, 32, 25, 18, 11, 4, 5,
  12, 19, 26, 33, 40, 48, 41, 34,
  27, 20, 13, 6, 7, 14, 21, 28,
  35, 42, 49, 56, 57, 50, 43, 36,
  29, 22, 15, 23, 30, 37, 44, 51,
  58, 59, 52, 45, 38, 31, 39, 46,
  53, 60, 61, 54, 47, 55, 62, 63,
};

static int num_bits(int value) {
  if (value == 0) return 0;
  int bits = 0;
  unsigned int v = static_cast<unsigned int>(value);
  while (v) {
    ++bits;
    v >>= 1;
  }
  return bits;
}

static void emit_eobrun(SymbolStats& stats, int eobrun) {
  if (eobrun <= 0) return;
  int rbits = num_bits(eobrun);
  stats.ac_symbols[rbits] += 1;
  stats.ac_bits += static_cast<uint64_t>(rbits);
  stats.ac_eob_symbols += 1;
}

struct CompState {
  JDIMENSION width = 0;
  JDIMENSION height = 0;
  std::vector<uint8_t> significant;
};

static std::vector<SymbolStats> compute_symbol_stats(JpegInput& input,
                                                     const std::vector<ScanLine>& scans) {
  std::vector<SymbolStats> stats;
  stats.reserve(scans.size());
  const int num_components = input.dinfo.num_components;
  std::vector<int> last_dc(num_components, 0);
  std::vector<CompState> comp_states(num_components);
  for (int ci = 0; ci < num_components; ++ci) {
    const auto& comp = input.dinfo.comp_info[ci];
    CompState state;
    state.width = comp.width_in_blocks;
    state.height = comp.height_in_blocks;
    state.significant.resize(static_cast<size_t>(state.width) * state.height * DCTSIZE2, 0);
    comp_states[ci] = std::move(state);
  }

  for (const auto& scan : scans) {
    SymbolStats stat;
    const bool is_dc = (scan.Ss == 0 && scan.Se == 0);
    if (is_dc) {
      for (JDIMENSION mcu_row = 0; mcu_row < input.dinfo.MCU_rows_in_scan; ++mcu_row) {
        for (JDIMENSION mcu_col = 0; mcu_col < input.dinfo.MCUs_per_row; ++mcu_col) {
          for (int comp_index : scan.components) {
            auto& comp = input.dinfo.comp_info[comp_index];
            for (int v = 0; v < comp.MCU_height; ++v) {
              JDIMENSION block_row = mcu_row * comp.MCU_height + v;
              if (block_row >= comp.height_in_blocks) continue;
              JBLOCKARRAY buffer = (*input.dinfo.mem->access_virt_barray)(
                  reinterpret_cast<j_common_ptr>(&input.dinfo),
                  input.coef_arrays[comp_index], block_row, 1, FALSE);
              for (int h = 0; h < comp.MCU_width; ++h) {
                JDIMENSION block_col = mcu_col * comp.MCU_width + h;
                if (block_col >= comp.width_in_blocks) continue;
                JCOEFPTR block = buffer[0][block_col];
                int dc = static_cast<int>(block[0]) >> scan.Al;
                if (scan.Ah == 0) {
                  int diff = dc - last_dc[comp_index];
                  last_dc[comp_index] = dc;
                  int size = num_bits(std::abs(diff));
                  if (size >= static_cast<int>(stat.dc_symbols.size())) size = stat.dc_symbols.size() - 1;
                  stat.dc_symbols[size] += 1;
                  stat.dc_bits += static_cast<uint64_t>(size);
                } else {
                  stat.dc_refine_bits += 1;
                }
              }
            }
          }
        }
      }
      stats.push_back(stat);
      continue;
    }

    int eobrun = 0;
    for (JDIMENSION mcu_row = 0; mcu_row < input.dinfo.MCU_rows_in_scan; ++mcu_row) {
      for (JDIMENSION mcu_col = 0; mcu_col < input.dinfo.MCUs_per_row; ++mcu_col) {
        for (int comp_index : scan.components) {
          auto& comp = input.dinfo.comp_info[comp_index];
          auto& state = comp_states[comp_index];
          for (int v = 0; v < comp.MCU_height; ++v) {
            JDIMENSION block_row = mcu_row * comp.MCU_height + v;
            if (block_row >= comp.height_in_blocks) continue;
            JBLOCKARRAY buffer = (*input.dinfo.mem->access_virt_barray)(
                reinterpret_cast<j_common_ptr>(&input.dinfo),
                input.coef_arrays[comp_index], block_row, 1, FALSE);
            for (int h = 0; h < comp.MCU_width; ++h) {
              JDIMENSION block_col = mcu_col * comp.MCU_width + h;
              if (block_col >= comp.width_in_blocks) continue;
              JCOEFPTR block = buffer[0][block_col];
              size_t block_index = static_cast<size_t>(block_row) * state.width + block_col;
              uint8_t* sig = state.significant.data() + block_index * DCTSIZE2;

              if (scan.Ah == 0) {
                int run = 0;
                int last_nonzero = -1;
                for (int k = scan.Ss; k <= scan.Se; ++k) {
                  int coef = block[kZigZag[k]];
                  int absval = std::abs(coef) >> scan.Al;
                  if (absval == 0) {
                    run++;
                    continue;
                  }
                  if (eobrun > 0) {
                    emit_eobrun(stat, eobrun);
                    eobrun = 0;
                  }
                  last_nonzero = k;
                  while (run > 15) {
                    stat.ac_symbols[0xF0] += 1;
                    stat.ac_zrl_symbols += 1;
                    run -= 16;
                  }
                  int size = num_bits(absval);
                  int symbol = (run << 4) | size;
                  stat.ac_symbols[symbol] += 1;
                  stat.ac_bits += static_cast<uint64_t>(size);
                  sig[kZigZag[k]] = 1;
                  run = 0;
                }
                if (last_nonzero < scan.Se) {
                  eobrun++;
                  if (eobrun == 0x7FFF) {
                    emit_eobrun(stat, eobrun);
                    eobrun = 0;
                  }
                }
              } else {
                int run = 0;
                int eob = -1;
                int absvals[DCTSIZE2];
                for (int k = scan.Ss; k <= scan.Se; ++k) {
                  int coef = block[kZigZag[k]];
                  int absval = std::abs(coef) >> scan.Al;
                  absvals[k] = absval;
                  if (absval == 1) eob = k;
                }
                int br = 0;
                for (int k = scan.Ss; k <= scan.Se; ++k) {
                  int absval = absvals[k];
                  if (absval == 0) {
                    run++;
                    continue;
                  }
                  if (absval > 1 || sig[kZigZag[k]]) {
                    stat.ac_refine_bits += 1;
                    br++;
                    continue;
                  }
                  if (eobrun > 0) {
                    emit_eobrun(stat, eobrun);
                    eobrun = 0;
                  }
                  while (run > 15 && k <= eob) {
                    stat.ac_symbols[0xF0] += 1;
                    stat.ac_zrl_symbols += 1;
                    run -= 16;
                    br = 0;
                  }
                  int symbol = (run << 4) | 1;
                  stat.ac_symbols[symbol] += 1;
                  stat.ac_bits += 1;
                  stat.ac_new_sign_bits += 1;
                  sig[kZigZag[k]] = 1;
                  run = 0;
                  br = 0;
                }
                if (run > 0 || br > 0) {
                  eobrun++;
                  if (eobrun == 0x7FFF) {
                    emit_eobrun(stat, eobrun);
                    eobrun = 0;
                  }
                }
              }
            }
          }
        }
      }
    }
    if (eobrun > 0) emit_eobrun(stat, eobrun);
    stats.push_back(stat);
  }

  return stats;
}

scan_set choose_best(const scan_set& left, const scan_set& right) {
  if (!left.size && !right.size) return left;
  if (!right.size || (left.size && left.size < right.size)) return left;
  if (!left.size || (right.size && right.size < left.size) || left.all_scans.size() > right.all_scans.size()) return right;
  return left;
}

scan_set try_bit_splits(const std::vector<int>& plane, int a, int b, int c, int d);
scan_set try_splits(const std::vector<int>& plane, int a, int b, int max_bits);
static size_t copy_markers(j_decompress_ptr srcinfo, j_compress_ptr dstinfo, StripMode strip_mode,
                           const std::set<int>& strip_apps);

static uint64_t hash_bytes(const uint8_t* data, size_t len, uint64_t seed) {
  uint64_t h = seed ? seed : 1469598103934665603ull;
  for (size_t i = 0; i < len; ++i) {
    h ^= data[i];
    h *= 1099511628211ull;
  }
  return h;
}

static DhtInfo collect_dht_info(j_compress_ptr cinfo) {
  DhtInfo info;
  std::ostringstream oss;
  bool first = true;
  for (int i = 0; i < NUM_HUFF_TBLS; ++i) {
    if (cinfo->dc_huff_tbl_ptrs[i]) {
      const JHUFF_TBL* tbl = cinfo->dc_huff_tbl_ptrs[i];
      uint64_t h = hash_bytes(tbl->bits + 1, 16, 0);
      h = hash_bytes(tbl->huffval, 256, h);
      if (!first) oss << ' ';
      oss << "DC" << i << "=" << std::hex << h << std::dec;
      first = false;
      size_t symbols = 0;
      for (int k = 1; k <= 16; ++k) symbols += tbl->bits[k];
      info.total_bytes += 17 + symbols;
    }
    if (cinfo->ac_huff_tbl_ptrs[i]) {
      const JHUFF_TBL* tbl = cinfo->ac_huff_tbl_ptrs[i];
      uint64_t h = hash_bytes(tbl->bits + 1, 16, 0);
      h = hash_bytes(tbl->huffval, 256, h);
      if (!first) oss << ' ';
      oss << "AC" << i << "=" << std::hex << h << std::dec;
      first = false;
      size_t symbols = 0;
      for (int k = 1; k <= 16; ++k) symbols += tbl->bits[k];
      info.total_bytes += 17 + symbols;
    }
  }
  info.summary = oss.str();
  return info;
}

static bool encode_with_scans_jpegli_bytes(JpegInput& input, const std::vector<ScanLine>& scans,
                                           bool arithmetic, StripMode strip_mode,
                                           const std::set<int>& strip_apps,
                                           std::vector<uint8_t>* out_bytes,
                                           DhtInfo* dht_info) {
  jpeg_compress_struct cinfo{};
  jpeg_error_mgr cerr{};
  cinfo.err = jpegli_std_error(&cerr);
  jpegli_create_compress(&cinfo);

  unsigned char* outbuffer = nullptr;
  unsigned long outsize = 0;
  jpegli_mem_dest(&cinfo, &outbuffer, &outsize);
  jpegli_copy_critical_parameters(&input.dinfo, &cinfo);
  cinfo.optimize_coding = !arithmetic;
  cinfo.arith_code = arithmetic;

  std::vector<jpeg_scan_info> scan_infos;
  if (!scans.empty()) {
    scan_infos.reserve(scans.size());
    for (const auto& scan : scans) scan_infos.push_back(to_jpeg_scan_info(scan));
    cinfo.scan_info = scan_infos.data();
    cinfo.num_scans = scan_infos.size();
  }

  jpegli_write_coefficients(&cinfo, input.coef_arrays);
  if (strip_mode != StripMode::None) {
    copy_markers(&input.dinfo, &cinfo, strip_mode, strip_apps);
  }
  jpegli_finish_compress(&cinfo);
  if (dht_info) {
    *dht_info = collect_dht_info(&cinfo);
  }

  if (out_bytes) out_bytes->assign(outbuffer, outbuffer + outsize);
  jpegli_destroy_compress(&cinfo);
  free(outbuffer);
  return true;
}

static std::vector<std::vector<uint8_t>> extract_dht_segments(const std::vector<uint8_t>& data) {
  std::vector<std::vector<uint8_t>> segments;
  size_t i = 0;
  while (i + 1 < data.size()) {
    if (data[i] != 0xFF) {
      ++i;
      continue;
    }
    size_t j = i;
    while (j < data.size() && data[j] == 0xFF) ++j;
    if (j >= data.size()) break;
    uint8_t marker = data[j];
    size_t tag_start = i;
    size_t cursor = j + 1;
    if (marker == 0xD8 || marker == 0xD9) {
      i = cursor;
      continue;
    }
    if (marker == 0x01 || (marker >= 0xD0 && marker <= 0xD7)) {
      i = cursor;
      continue;
    }
    if (cursor + 1 >= data.size()) break;
    uint16_t len = static_cast<uint16_t>((data[cursor] << 8) | data[cursor + 1]);
    size_t seg_end = cursor + len;
    if (seg_end > data.size()) break;
    if (marker == 0xC4) {
      segments.emplace_back(data.begin() + tag_start, data.begin() + seg_end);
    }
    i = seg_end;
  }
  return segments;
}

static std::vector<HuffmanTable> parse_dht_tables_from_bytes(const std::vector<uint8_t>& data) {
  std::vector<HuffmanTable> tables;
  size_t i = 0;
  bool in_scan = false;
  while (i + 1 < data.size()) {
    if (!in_scan) {
      if (data[i] != 0xFF) {
        ++i;
        continue;
      }
      size_t j = i;
      while (j < data.size() && data[j] == 0xFF) ++j;
      if (j >= data.size()) break;
      uint8_t marker = data[j];
      size_t cursor = j + 1;
      if (marker == 0xD8) {
        i = cursor;
        continue;
      }
      if (marker == 0xD9) break;
      if (marker == 0xDA) {
        if (cursor + 1 >= data.size()) break;
        uint16_t len = static_cast<uint16_t>((data[cursor] << 8) | data[cursor + 1]);
        size_t seg_end = cursor + len;
        if (seg_end > data.size()) break;
        i = seg_end;
        in_scan = true;
        continue;
      }
      if (marker == 0x01 || (marker >= 0xD0 && marker <= 0xD7)) {
        i = cursor;
        continue;
      }
      if (cursor + 1 >= data.size()) break;
      uint16_t len = static_cast<uint16_t>((data[cursor] << 8) | data[cursor + 1]);
      size_t seg_end = cursor + len;
      if (seg_end > data.size()) break;
      if (marker == 0xC4) {
        size_t k = cursor + 2;
        while (k < seg_end) {
          uint8_t htinfo = data[k++];
          int tc = (htinfo >> 4) & 0xF;
          int tid = htinfo & 0xF;
          if (k + 16 > seg_end) break;
          HuffmanTable table;
          table.table_class = tc;
          table.table_id = tid;
          for (int b = 0; b < 16; ++b) {
            table.bits[static_cast<size_t>(b)] = data[k + b];
          }
          k += 16;
          int count = 0;
          for (int b = 0; b < 16; ++b) count += table.bits[static_cast<size_t>(b)];
          table.symbols = count;
          if (k + static_cast<size_t>(count) > seg_end) break;
          for (int b = 0; b < count; ++b) {
            table.huffval[static_cast<size_t>(b)] = data[k + b];
          }
          k += static_cast<size_t>(count);
          tables.push_back(table);
        }
      }
      i = seg_end;
    } else {
      size_t j = i;
      while (j + 1 < data.size()) {
        if (data[j] == 0xFF && data[j + 1] != 0x00) break;
        ++j;
      }
      i = j;
      in_scan = false;
    }
  }
  return tables;
}

static DhtInfo compute_dht_info_from_tables(const std::vector<HuffmanTable>& tables) {
  DhtInfo info;
  std::ostringstream oss;
  bool first = true;
  for (int tc = 0; tc <= 1; ++tc) {
    for (int tid = 0; tid < NUM_HUFF_TBLS; ++tid) {
      const HuffmanTable* table = nullptr;
      for (const auto& t : tables) {
        if (t.table_class == tc && t.table_id == tid) {
          table = &t;
          break;
        }
      }
      if (!table) continue;
      uint64_t h = hash_bytes(table->bits.data(), 16, 0);
      h = hash_bytes(table->huffval.data(), static_cast<size_t>(table->symbols), h);
      if (!first) oss << ' ';
      oss << (tc == 0 ? "DC" : "AC") << tid << "=" << std::hex << h << std::dec;
      first = false;
      info.total_bytes += 17 + static_cast<size_t>(table->symbols);
    }
  }
  info.summary = oss.str();
  return info;
}

static std::vector<ScanDhtInfo> scan_dht_states(const std::vector<uint8_t>& data) {
  std::vector<ScanDhtInfo> results;
  std::vector<HuffmanTable> tables;
  tables.reserve(8);
  size_t i = 0;
  bool in_scan = false;
  while (i + 1 < data.size()) {
    if (!in_scan) {
      if (data[i] != 0xFF) {
        ++i;
        continue;
      }
      size_t j = i;
      while (j < data.size() && data[j] == 0xFF) ++j;
      if (j >= data.size()) break;
      uint8_t marker = data[j];
      size_t cursor = j + 1;
      if (marker == 0xD8) {
        i = cursor;
        continue;
      }
      if (marker == 0xD9) break;
      if (marker == 0xDA) {
        if (cursor + 1 >= data.size()) break;
        uint16_t len = static_cast<uint16_t>((data[cursor] << 8) | data[cursor + 1]);
        size_t seg_end = cursor + len;
        if (seg_end > data.size()) break;
        ScanDhtInfo info;
        DhtInfo summary = compute_dht_info_from_tables(tables);
        info.hash = summary.summary;
        info.bytes = summary.total_bytes;
        info.valid = !info.hash.empty();
        results.push_back(info);
        i = seg_end;
        in_scan = true;
        continue;
      }
      if (marker == 0x01 || (marker >= 0xD0 && marker <= 0xD7)) {
        i = cursor;
        continue;
      }
      if (cursor + 1 >= data.size()) break;
      uint16_t len = static_cast<uint16_t>((data[cursor] << 8) | data[cursor + 1]);
      size_t seg_end = cursor + len;
      if (seg_end > data.size()) break;
      if (marker == 0xC4) {
        size_t k = cursor + 2;
        while (k < seg_end) {
          uint8_t htinfo = data[k++];
          int tc = (htinfo >> 4) & 0xF;
          int tid = htinfo & 0xF;
          if (k + 16 > seg_end) break;
          HuffmanTable table;
          table.table_class = tc;
          table.table_id = tid;
          for (int b = 0; b < 16; ++b) {
            table.bits[static_cast<size_t>(b)] = data[k + b];
          }
          k += 16;
          int count = 0;
          for (int b = 0; b < 16; ++b) count += table.bits[static_cast<size_t>(b)];
          table.symbols = count;
          if (k + static_cast<size_t>(count) > seg_end) break;
          for (int b = 0; b < count; ++b) {
            table.huffval[static_cast<size_t>(b)] = data[k + b];
          }
          k += static_cast<size_t>(count);
          bool replaced = false;
          for (auto& existing : tables) {
            if (existing.table_class == tc && existing.table_id == tid) {
              existing = table;
              replaced = true;
              break;
            }
          }
          if (!replaced) tables.push_back(table);
        }
      }
      i = seg_end;
    } else {
      size_t j = i;
      while (j + 1 < data.size()) {
        if (data[j] == 0xFF && data[j + 1] != 0x00) break;
        ++j;
      }
      i = j;
      in_scan = false;
    }
  }
  return results;
}


static bool encode_with_scans(JpegInput& input, const std::vector<ScanLine>& scans, bool arithmetic,
                              StripMode strip_mode, const std::set<int>& strip_apps,
                              std::vector<uint8_t>* out_bytes, ScanSizes* out_scan_sizes,
                              DhtInfo* dht_info) {
  if (g_encoder_backend == EncoderBackend::Turbo) {
#if ENABLE_TURBO
    std::vector<uint8_t> data;
    if (!encode_with_scans_turbo(input.data, scans, arithmetic, &data)) {
      return false;
    }
    if (out_bytes) *out_bytes = std::move(data);
    if (out_scan_sizes && out_bytes) *out_scan_sizes = scan_sizes(*out_bytes);
    if (dht_info) {
      dht_info->summary = "turbo";
      dht_info->total_bytes = 0;
    }
    return true;
#else
    std::cerr << "Turbo encoder not enabled in this build.\n";
    return false;
#endif
  }

  if (g_encoder_backend == EncoderBackend::Mozjpeg) {
#if ENABLE_MOZJPEG
    if (scans.size() > 64) {
      return false;
    }
    std::vector<uint8_t> data;
    if (!encode_with_scans_mozjpeg_from_coeffs(&input.dinfo, input.coef_arrays,
                                               scans, arithmetic, false, &data)) {
      return false;
    }
    if (out_bytes) *out_bytes = std::move(data);
    if (out_scan_sizes && out_bytes) *out_scan_sizes = scan_sizes(*out_bytes);
    if (dht_info) {
      dht_info->summary = "mozjpeg";
      dht_info->total_bytes = 0;
    }
    return true;
#else
    std::cerr << "Mozjpeg encoder not enabled in this build.\n";
    return false;
#endif
  }

  std::vector<uint8_t> data;
  DhtInfo info;
  if (!encode_with_scans_jpegli_bytes(input, scans, arithmetic, strip_mode, strip_apps,
                                      &data, dht_info ? &info : nullptr)) {
    return false;
  }
  if (dht_info) {
    *dht_info = info;
  }
  if (out_bytes) *out_bytes = std::move(data);
  if (out_scan_sizes && out_bytes) *out_scan_sizes = scan_sizes(*out_bytes);
  return true;
}

static size_t count_ac_contexts(const std::vector<ScanLine>& scans) {
  size_t contexts = 0;
  for (const auto& scan : scans) {
    if (scan.Se > 0) {
      contexts += scan.components.size();
    }
  }
  return contexts;
}

static bool apply_scan_to_mask(const ScanLine& scan, std::vector<std::vector<uint16_t>>& comp_mask) {
  for (int ci : scan.components) {
    if (scan.Ss != 0 && comp_mask[ci][0] == 0) return false;
    for (int k = scan.Ss; k <= scan.Se; ++k) {
      if (comp_mask[ci][k] == 0) {
        if (scan.Ah != 0) return false;
        comp_mask[ci][k] = static_cast<uint16_t>((0xffffu << scan.Al) & 0xffffu);
      } else {
        if (comp_mask[ci][k] != static_cast<uint16_t>((0xffffu << scan.Ah) & 0xffffu) ||
            scan.Al != scan.Ah - 1) {
          return false;
        }
        comp_mask[ci][k] = static_cast<uint16_t>(comp_mask[ci][k] | (1u << scan.Al));
      }
    }
  }
  return true;
}

static bool complete_scan_script(const std::vector<ScanLine>& base_scans, int num_components,
                                 std::vector<ScanLine>& full_scans) {
  std::vector<std::vector<uint16_t>> comp_mask(
      num_components, std::vector<uint16_t>(DCTSIZE2, 0));
  full_scans = base_scans;
  for (const auto& scan : base_scans) {
    if (!apply_scan_to_mask(scan, comp_mask)) return false;
  }

  auto expected_mask = [](int level) -> uint16_t {
    return static_cast<uint16_t>((0xffffu << level) & 0xffffu);
  };

  for (int c = 0; c < num_components; ++c) {
    // DC completion
    if (comp_mask[c][0] == 0) {
      ScanLine scan = make_scan(std::vector<int>{c}, 0, 0, 0, 0);
      full_scans.push_back(scan);
      comp_mask[c][0] = 0xffffu;
    } else {
      for (int level = 15; level > 0; --level) {
        if (comp_mask[c][0] == expected_mask(level)) {
          ScanLine scan = make_scan(std::vector<int>{c}, 0, 0, level, level - 1);
          full_scans.push_back(scan);
          comp_mask[c][0] = expected_mask(level - 1);
        }
      }
    }

    // AC completion: base scans for missing coefficients
    int start = -1;
    for (int k = 1; k < DCTSIZE2; ++k) {
      if (comp_mask[c][k] == 0) {
        if (start < 0) start = k;
      } else if (start >= 0) {
        ScanLine scan = make_scan(std::vector<int>{c}, start, k - 1, 0, 0);
        full_scans.push_back(scan);
        for (int j = start; j <= k - 1; ++j) comp_mask[c][j] = 0xffffu;
        start = -1;
      }
    }
    if (start >= 0) {
      ScanLine scan = make_scan(std::vector<int>{c}, start, DCTSIZE2 - 1, 0, 0);
      full_scans.push_back(scan);
      for (int j = start; j <= DCTSIZE2 - 1; ++j) comp_mask[c][j] = 0xffffu;
    }

    // AC completion: refinement scans for partial coefficients
    for (int level = 15; level > 0; --level) {
      int run_start = -1;
      for (int k = 1; k < DCTSIZE2; ++k) {
        if (comp_mask[c][k] == expected_mask(level)) {
          if (run_start < 0) run_start = k;
        } else if (run_start >= 0) {
          ScanLine scan = make_scan(std::vector<int>{c}, run_start, k - 1, level, level - 1);
          full_scans.push_back(scan);
          for (int j = run_start; j <= k - 1; ++j) comp_mask[c][j] = expected_mask(level - 1);
          run_start = -1;
        }
      }
      if (run_start >= 0) {
        ScanLine scan = make_scan(std::vector<int>{c}, run_start, DCTSIZE2 - 1, level, level - 1);
        full_scans.push_back(scan);
        for (int j = run_start; j <= DCTSIZE2 - 1; ++j) comp_mask[c][j] = expected_mask(level - 1);
      }
    }
  }

  for (int c = 0; c < num_components; ++c) {
    for (int k = 0; k < DCTSIZE2; ++k) {
      if (comp_mask[c][k] != 0xffffu) return false;
    }
  }
  return true;
}

static void dump_scan_set(const std::vector<ScanLine>& scans, const ScanSizes& sizes_local,
                          const DhtInfo* dht_info, const std::vector<SegmentInfo>* segments) {
  if (!g_scan_dump) return;
  if (dht_info) {
    (*g_scan_dump) << "#DHT " << dht_info->summary
                   << " bytes=" << dht_info->total_bytes
                   << " stream_dht_bytes=" << sizes_local.total_dht_bytes
                   << " stream_dht_markers=" << sizes_local.total_dht_markers
                   << " pre_sos_dht_bytes=" << sizes_local.pre_sos_dht_bytes
                   << "\n";
  }
  if (g_score_dht) {
    (*g_scan_dump) << "#COLS scan\tsize\tdht_bytes\tdht_hash\tdht_table_bytes\n";
  }
  if (segments) {
    (*g_scan_dump) << "#SEG";
    for (const auto& seg : *segments) {
      const char* name = marker_name(seg.marker);
      (*g_scan_dump) << " " << name << "@0x" << std::hex << seg.offset
                     << "..0x" << seg.end_offset << std::dec;
    }
    (*g_scan_dump) << "\n";
  }
  for (size_t i = 0; i < scans.size(); ++i) {
    size_t dht_bytes = i < sizes_local.dht_bytes.size() ? sizes_local.dht_bytes[i] : 0;
    (*g_scan_dump) << scan_label(scans[i]) << "\t" << sizes_local.sizes[i] << "\t" << dht_bytes;
    if (g_score_dht) {
      std::string key = scan_key(scans[i]);
      auto it = scan_dht_info.find(key);
      if (it != scan_dht_info.end() && it->second.valid) {
        (*g_scan_dump) << "\t" << it->second.hash << "\t" << it->second.bytes;
      }
    }
    (*g_scan_dump) << "\n";
  }
  (*g_scan_dump) << "\n";
}

static void dump_symbol_set(const std::vector<ScanLine>& scans,
                            const std::vector<SymbolStats>& stats) {
  if (!g_symbol_dump) return;
  (*g_symbol_dump) << "#SYM\n";
  size_t count = std::min(scans.size(), stats.size());
  for (size_t i = 0; i < count; ++i) {
    const auto& st = stats[i];
    (*g_symbol_dump) << scan_label(scans[i])
                     << "\tDCbits=" << st.dc_bits
                     << " DCref=" << st.dc_refine_bits
                     << " ACbits=" << st.ac_bits
                     << " ACref=" << st.ac_refine_bits
                     << " ACnew=" << st.ac_new_sign_bits
                     << " EOB=" << st.ac_eob_symbols
                     << " ZRL=" << st.ac_zrl_symbols;

    bool first = true;
    for (size_t k = 0; k < st.dc_symbols.size(); ++k) {
      if (st.dc_symbols[k] == 0) continue;
      (*g_symbol_dump) << (first ? " DCsym{" : ",");
      (*g_symbol_dump) << k << "=" << st.dc_symbols[k];
      first = false;
    }
    if (!first) (*g_symbol_dump) << "}";

    first = true;
    for (size_t k = 0; k < st.ac_symbols.size(); ++k) {
      if (st.ac_symbols[k] == 0) continue;
      (*g_symbol_dump) << (first ? " ACsym{" : ",");
      (*g_symbol_dump) << "0x" << std::hex << std::setw(2) << std::setfill('0')
                       << k << std::dec << "=" << st.ac_symbols[k];
      first = false;
    }
    if (!first) (*g_symbol_dump) << "}";

    (*g_symbol_dump) << "\n";
  }
  (*g_symbol_dump) << "\n";
}

static void dump_scan_set_lines(const std::vector<ScanLine>& scans) {
  if (!g_scan_set_dump) return;
  for (const auto& scan : scans) {
    (*g_scan_set_dump) << scan_label(scan) << "\n";
  }
  (*g_scan_set_dump) << "\n";
}

static void dump_fail_scan_set(const std::vector<ScanLine>& scans, const std::string& reason) {
  if (!g_fail_scan_dump) return;
  (*g_fail_scan_dump) << "#FAIL " << reason << "\n";
  for (const auto& scan : scans) {
    (*g_fail_scan_dump) << scan_label(scan) << "\n";
  }
  (*g_fail_scan_dump) << "\n";
}

static void dump_component_only_set(const std::vector<std::vector<int>>& groups) {
  if (!g_scan_set_dump) return;
  for (const auto& group : groups) {
    (*g_scan_set_dump) << component_only_label(group) << "\n";
  }
  (*g_scan_set_dump) << "\n";
}

static void do_scans(JpegInput& input, std::vector<std::vector<ScanLine>>& scan_queue, bool arithmetic,
                     StripMode strip_mode, const std::set<int>& strip_apps, bool verbose, bool quiet) {
  g_lock_wait_ns.store(0, std::memory_order_relaxed);
  g_lock_hold_ns.store(0, std::memory_order_relaxed);
  g_lock_acquires.store(0, std::memory_order_relaxed);
  const auto batch_start = std::chrono::steady_clock::now();
  g_do_scans_calls.fetch_add(1, std::memory_order_relaxed);
  auto with_lock = [&](std::mutex* lock, const auto& fn) {
    if (!lock) {
      fn();
      return;
    }
    auto wait_start = std::chrono::steady_clock::now();
    lock->lock();
    auto lock_start = std::chrono::steady_clock::now();
    g_lock_wait_ns.fetch_add(
        static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(lock_start - wait_start).count()),
        std::memory_order_relaxed);
    g_lock_acquires.fetch_add(1, std::memory_order_relaxed);
    fn();
    auto lock_end = std::chrono::steady_clock::now();
    g_lock_hold_ns.fetch_add(
        static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(lock_end - lock_start).count()),
        std::memory_order_relaxed);
    lock->unlock();
  };

  auto process_one = [&](JpegInput& local_input, const std::vector<ScanLine>& scans,
                         std::mutex* lock) {
    with_lock(lock, [&] { dump_scan_set_lines(scans); });
    std::vector<ScanLine> full_scans;
    if (g_encoder_backend == EncoderBackend::Mozjpeg) {
      full_scans = scans;
    } else {
      if (!complete_scan_script(scans, local_input.dinfo.num_components, full_scans)) {
        with_lock(lock, [&] {
          for (const auto& scan : scans) {
            std::string key = scan_key(scan);
            if (!sizes.contains(key)) sizes[key] = kInvalidSize;
          }
          if (verbose) {
            std::cerr << "\nSkipping invalid scan script:\n";
            for (const auto& scan : scans) std::cerr << scan_label(scan) << '\n';
          }
        });
        return;
      }
    }
    if (count_ac_contexts(full_scans) > 252) {
      with_lock(lock, [&] {
        for (const auto& scan : scans) {
          std::string key = scan_key(scan);
          if (!sizes.contains(key)) sizes[key] = kInvalidSize;
        }
        if (verbose) {
          std::cerr << "\nSkipping scan set with too many AC scans:\n";
          for (const auto& scan : scans) std::cerr << scan_label(scan) << '\n';
        }
        dump_fail_scan_set(scans, "too-many-ac");
      });
      return;
    }
    std::vector<SymbolStats> symbol_stats;
    if (g_symbol_dump) {
      symbol_stats = compute_symbol_stats(local_input, full_scans);
    }
    std::vector<uint8_t> data;
    ScanSizes sizes_local;
    DhtInfo dht_info;
    DhtInfo* dht_ptr = g_dump_dht ? &dht_info : nullptr;
    if (!encode_with_scans(local_input, full_scans, arithmetic, strip_mode, strip_apps,
                           &data, &sizes_local, dht_ptr)) {
      with_lock(lock, [&] { dump_fail_scan_set(scans, "encode-failed"); });
      return;
    }
    if (sizes_local.sizes.size() < scans.size()) {
      with_lock(lock, [&] {
        std::cerr << "\nIncorrect encoding for scans:\n";
        for (const auto& scan : scans) std::cerr << scan_label(scan) << '\n';
        dump_fail_scan_set(scans, "sizes-mismatch");
      });
      return;
    }
    std::vector<ScanDhtInfo> dht_states;
    if (g_score_dht) {
      dht_states = scan_dht_states(data);
      if (dht_states.size() < scans.size()) {
        with_lock(lock, [&] {
          std::cerr << "\nMissing DHT state for scans:\n";
          for (const auto& scan : scans) std::cerr << scan_label(scan) << '\n';
          dump_fail_scan_set(scans, "missing-dht-state");
        });
        return;
      }
    }
    std::vector<SegmentInfo> segments;
    std::vector<SegmentInfo>* seg_ptr = nullptr;
    if (g_dump_segments) {
      segments = collect_segments(data);
      seg_ptr = &segments;
    }
    with_lock(lock, [&] {
      if (g_score_dht) {
        for (size_t i = 0; i < scans.size(); ++i) {
          std::string key = scan_key(scans[i]);
          if (!scan_dht_info.contains(key)) {
            scan_dht_info[key] = dht_states[i];
          }
        }
      }
      dump_scan_set(scans, sizes_local, dht_ptr, seg_ptr);
      if (g_symbol_dump) {
        dump_symbol_set(scans, symbol_stats);
      }
      if (!verbose && !quiet) {
        std::cout << '.';
      }
      for (size_t i = 0; i < scans.size(); ++i) {
        std::string key = scan_key(scans[i]);
        if (!sizes.contains(key)) sizes[key] = sizes_local.sizes[i];
      }
    });
  };

  if (g_threads <= 1 || scan_queue.size() < 2) {
    for (const auto& scans : scan_queue) {
      process_one(input, scans, nullptr);
    }
    scan_queue.clear();
    const auto batch_end = std::chrono::steady_clock::now();
    const auto wall_ns = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(batch_end - batch_start).count());
    g_do_scans_wall_ns.fetch_add(wall_ns, std::memory_order_relaxed);
    g_do_scans_capacity_ns.fetch_add(wall_ns, std::memory_order_relaxed);
    g_worker_busy_ns.fetch_add(wall_ns, std::memory_order_relaxed);
    return;
  }

  const size_t total = scan_queue.size();
  const int workers = std::min<int>(g_threads, static_cast<int>(total));
  std::atomic<size_t> index{0};
  std::mutex lock;
  std::atomic<uint64_t> batch_busy_ns{0};
  std::vector<std::thread> threads;
  threads.reserve(workers);
  for (int t = 0; t < workers; ++t) {
    threads.emplace_back([&] {
      JpegInput local_input;
      if (!prepare_input(local_input, input.data, strip_mode)) {
        return;
      }
      uint64_t local_busy = 0;
      while (true) {
        size_t i = index.fetch_add(1);
        if (i >= total) break;
        auto work_start = std::chrono::steady_clock::now();
        process_one(local_input, scan_queue[i], &lock);
        auto work_end = std::chrono::steady_clock::now();
        local_busy += static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(work_end - work_start).count());
      }
      batch_busy_ns.fetch_add(local_busy, std::memory_order_relaxed);
      destroy_input(local_input);
    });
  }
  for (auto& th : threads) th.join();
  scan_queue.clear();

  if (!quiet) {
    auto batch_end = std::chrono::steady_clock::now();
    const auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(batch_end - batch_start).count();
    const auto wait_ms = g_lock_wait_ns.load(std::memory_order_relaxed) / 1'000'000;
    const auto hold_ms = g_lock_hold_ns.load(std::memory_order_relaxed) / 1'000'000;
    const auto acquires = g_lock_acquires.load(std::memory_order_relaxed);
    const auto busy_ms = batch_busy_ns.load(std::memory_order_relaxed) / 1'000'000;
    const auto capacity_ms = static_cast<uint64_t>(total_ms) * static_cast<uint64_t>(workers);
    const double util = capacity_ms ? (static_cast<double>(busy_ms) / static_cast<double>(capacity_ms)) * 100.0 : 0.0;
    std::cout << "\nLock stats: wait=" << wait_ms << "ms hold=" << hold_ms
              << "ms acquires=" << acquires << " batch=" << total_ms << "ms"
              << " busy=" << busy_ms << "ms util=" << std::fixed << std::setprecision(1)
              << util << "%\n";

    const auto wall_ns = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(batch_end - batch_start).count());
    g_do_scans_wall_ns.fetch_add(wall_ns, std::memory_order_relaxed);
    g_do_scans_capacity_ns.fetch_add(static_cast<uint64_t>(workers) * wall_ns, std::memory_order_relaxed);
    g_worker_busy_ns.fetch_add(batch_busy_ns.load(std::memory_order_relaxed), std::memory_order_relaxed);
  }
}

static scan_set get_dc_optimal(JpegInput& input, const std::vector<int>& plane, bool arithmetic,
                               StripMode strip_mode, const std::set<int>& strip_apps, bool verbose, bool quiet) {
  std::string key = components_key(plane);
  if (!dc_cache.contains(key)) {
    std::vector<std::vector<ScanLine>> scan_queue;
    generate_dc_scans(plane, g_all_components, g_max_bits, scan_queue);
    do_scans(input, scan_queue, arithmetic, strip_mode, strip_apps, verbose, quiet);
    dc_cache[key] = try_splits(plane, 0, 0, g_max_bits);
  }
  return dc_cache[key];
}

static void partitionscans(JpegInput& input, size_t offset, const std::vector<std::vector<int>>& partition,
                           bool arithmetic, StripMode strip_mode, const std::set<int>& strip_apps, bool verbose, bool quiet) {
  scan_set dcall_scans;
  for (const auto& group : partition) {
    auto dc_scans = get_dc_optimal(input, group, arithmetic, strip_mode, strip_apps, verbose, quiet);
    dcall_scans += dc_scans;
  }
  best_scans = choose_best(dcall_scans + ac_all, best_scans);

  scan_set all_scans;
  std::vector<std::vector<int>> groups = partition;
  for (const auto& group : partition) {
    ScanLine base = make_scan(group, 0, DCTSIZE2 - 1, 0, 0);
    std::string key = scan_key(base);
    if (!sizes.contains(key)) {
      std::vector<std::vector<ScanLine>> scan_queue;
      generate_dc_scans(group, g_all_components, g_max_bits, scan_queue);
      do_scans(input, scan_queue, arithmetic, strip_mode, strip_apps, verbose, quiet);
    }
    all_scans.all_scans.push_back(base);
    all_scans.size = add_sizes(all_scans.size, sizes.contains(key) ? sizes[key] : kInvalidSize);
  }
  best_scans = choose_best(all_scans, best_scans);

  bool need_component_only = false;
  for (const auto& group : groups) {
    std::string key = components_key(group);
    if (!component_only_sizes.contains(key)) {
      need_component_only = true;
      break;
    }
  }
  if (need_component_only) {
    dump_component_only_set(groups);
    for (const auto& group : groups) {
      component_only_sizes.insert(components_key(group));
    }
  }

  for (size_t i = offset; i < partition.size(); ++i) {
    for (size_t j = i + 1; j < partition.size(); ++j) {
      std::vector<std::vector<int>> merged = partition;
      merged[i].insert(merged[i].end(), merged[j].begin(), merged[j].end());
      merged.erase(merged.begin() + j);
      partitionscans(input, i, merged, arithmetic, strip_mode, strip_apps, verbose, quiet);
    }
  }
}

static std::vector<ScanLine> order_scans(const std::vector<ScanLine>& scans) {
  std::vector<ScanLine> sorted = scans;
  std::sort(sorted.begin(), sorted.end(), [](const ScanLine& a, const ScanLine& b) {
    if (a.Se != b.Se) return a.Se < b.Se;
    return components_key(a.components) < components_key(b.components);
  });

  size_t n = sorted.size();
  std::vector<std::vector<size_t>> prev(n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      if (i == j) continue;
      const auto& ai = sorted[i];
      const auto& bj = sorted[j];
      if (components_key(ai.components) == components_key(bj.components)) {
        bool overlap = (ai.Ss <= bj.Ss && bj.Ss <= ai.Se) || (bj.Ss <= ai.Ss && ai.Ss <= bj.Se);
        if (overlap && ai.Al > bj.Al) {
          prev[j].push_back(i);
          continue;
        }
      }
      if (components_key(bj.components) > components_key(ai.components) && bj.Se >= ai.Se) {
        prev[j].push_back(i);
      }
    }
  }

  std::vector<ScanLine> ordered;
  std::vector<bool> used(n, false);
  while (ordered.size() < n) {
    bool progressed = false;
    for (size_t i = 0; i < n; ++i) {
      if (used[i]) continue;
      bool ready = true;
      for (size_t p : prev[i]) {
        if (!used[p]) {
          ready = false;
          break;
        }
      }
      if (!ready) continue;
      ordered.push_back(sorted[i]);
      used[i] = true;
      progressed = true;
    }
    if (!progressed) break;
  }
  return ordered.empty() ? sorted : ordered;
}

static size_t compute_scan_set_cost(const std::vector<ScanLine>& scans) {
  size_t total = 0;
  if (scans.empty()) return total;
  if (!g_score_dht) {
    for (const auto& scan : scans) {
      std::string key = scan_key(scan);
      if (!sizes.contains(key)) return kInvalidSize;
      total = add_sizes(total, sizes[key]);
    }
    return total;
  }
  auto ordered = order_scans(scans);
  std::string prev_hash;
  for (const auto& scan : ordered) {
    std::string key = scan_key(scan);
    if (!sizes.contains(key)) return kInvalidSize;
    total = add_sizes(total, sizes[key]);
    auto it = scan_dht_info.find(key);
    if (it == scan_dht_info.end() || !it->second.valid) return kInvalidSize;
    if (it->second.hash != prev_hash) {
      total = add_sizes(total, it->second.bytes);
      prev_hash = it->second.hash;
    }
  }
  return total;
}

static bool get_optimal(const std::string& key, scan_set* out) {
  std::lock_guard<std::mutex> guard(g_optimal_mutex);
  auto it = optimal.find(key);
  if (it == optimal.end()) return false;
  *out = it->second;
  return true;
}

static scan_set set_optimal_if_absent(const std::string& key, const scan_set& value) {
  std::lock_guard<std::mutex> guard(g_optimal_mutex);
  auto [it, inserted] = optimal.emplace(key, value);
  return it->second;
}

static bool get_optimal2(const std::string& key, scan_set* out) {
  std::lock_guard<std::mutex> guard(g_optimal2_mutex);
  auto it = optimal2.find(key);
  if (it == optimal2.end()) return false;
  *out = it->second;
  return true;
}

static scan_set set_optimal2_if_absent(const std::string& key, const scan_set& value) {
  std::lock_guard<std::mutex> guard(g_optimal2_mutex);
  auto [it, inserted] = optimal2.emplace(key, value);
  return it->second;
}

scan_set try_bit_splits(const std::vector<int>& plane, int a, int b, int c, int d) {
  ScanLine min_scan = make_scan(plane, a, b, c, d);
  std::string key = scan_key(min_scan);
  scan_set cached;
  if (get_optimal2(key, &cached)) {
    return cached;
  }
  scan_set min_scans;
  min_scans.all_scans.push_back(min_scan);
  min_scans.size = compute_scan_set_cost(min_scans.all_scans);
  for (auto i = a; i < b; ++i) {
    auto a_scans = try_bit_splits(plane, a, i, c, d);
    auto b_scans = try_bit_splits(plane, i + 1, b, c, d);
    min_scans = choose_best(a_scans + b_scans, min_scans);
  }
  return set_optimal2_if_absent(key, min_scans);
}

scan_set try_splits(const std::vector<int>& plane, int a, int b, int max_bits) {
  std::ostringstream oss;
  oss << components_key(plane) << ": " << a << ' ' << b << " split";
  std::string key = oss.str();
  scan_set cached;
  if (get_optimal(key, &cached)) {
    return cached;
  }
  scan_set min_scans;
  min_scans.size = 0;
  for (auto i = 0; i <= max_bits; ++i) {
    auto all_scans = try_bit_splits(plane, a, b, 0, i);
    for (auto j = i; j > 0; --j) {
      all_scans += try_bit_splits(plane, a, b, j, j - 1);
    }
    min_scans = choose_best(all_scans, min_scans);
  }
  for (auto i = a; i < b; ++i) {
    auto a_scans = try_splits(plane, a, i, max_bits);
    auto b_scans = try_splits(plane, i + 1, b, max_bits);
    min_scans = choose_best(a_scans + b_scans, min_scans);
  }
  return set_optimal_if_absent(key, min_scans);
}

static size_t copy_markers(j_decompress_ptr srcinfo, j_compress_ptr dstinfo, StripMode strip_mode,
                           const std::set<int>& strip_apps) {
  size_t size = 0;
  for (jpeg_saved_marker_ptr marker = srcinfo->marker_list; marker; marker = marker->next) {
    if (marker->marker == JPEG_COM) {
      if (strip_mode == StripMode::None) continue;
      if (strip_mode != StripMode::Comments && strip_mode != StripMode::All) continue;
    } else if (marker->marker >= JPEG_APP0 && marker->marker <= JPEG_APP0 + 15) {
      int app_id = marker->marker - JPEG_APP0;
      if (strip_apps.contains(app_id)) continue;
      if (strip_mode != StripMode::All) continue;
      if (dstinfo->write_JFIF_header &&
          marker->marker == JPEG_APP0 &&
          marker->data_length >= 5 &&
          marker->data[0] == 'J' &&
          marker->data[1] == 'F' &&
          marker->data[2] == 'I' &&
          marker->data[3] == 'F' &&
          marker->data[4] == 0) {
        continue;
      }
      if (dstinfo->write_Adobe_marker &&
          marker->marker == JPEG_APP14 &&
          marker->data_length >= 5 &&
          marker->data[0] == 'A' &&
          marker->data[1] == 'd' &&
          marker->data[2] == 'o' &&
          marker->data[3] == 'b' &&
          marker->data[4] == 'e') {
        continue;
      }
    } else if (strip_mode != StripMode::All) {
      continue;
    }
    jpegli_write_marker(dstinfo, marker->marker, marker->data, marker->data_length);
    if (marker->marker == JPEG_COM ||
        (marker->marker >= JPEG_APP0 && marker->marker <= JPEG_APP0 + 15)) {
      size += marker->data_length;
    }
  }
  return size;
}

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
  const auto program_start = std::chrono::steady_clock::now();
  po::options_description desc("Switches");
  desc.add_options()
    ("s,s", po::value<std::string>()->implicit_value("none")->default_value("all"), "Strip extra markers")
    ("i,i", po::value<int>()->implicit_value(1)->default_value(0), "Allow multiple planes per DC scan, which may be incompatible with some software")
    ("j,j", po::value<std::vector<int>>()->implicit_value(std::vector<int>{0}, "0")->composing(), "Strip APP0 or APPn segment")
    ("a,a", "Use arithmetic coding (unsupported by most software)")
    ("b,b", po::value<int>()->default_value(3), "Maximum number of bit splits to test (default: 3)")
    ("t,t", po::value<int>()->implicit_value(0)->default_value(1),
     "Number of simultaneous processes (default: 1, -t uses all cores)")
    ("q,q", "Suppress all output")
    ("v,v", "Verbose output")
    ("h,h", "Help")
    ("?", "Help")
    ("dump-scan-sizes", po::value<std::string>(), "Dump per-scan sizes to file")
    ("dump-dht", "Include DHT summary in scan sizes dump")
    ("dump-segments", "Include segment sequence (type@offset) in scan sizes dump")
    ("dump-symbols", po::value<std::string>(), "Dump per-scan symbol stats to file")
    ("dump-scan-sets", po::value<std::string>(), "Dump scan sets to file")
    ("dump-fail-scans", po::value<std::string>(), "Dump scan sets that fail encoding")
    ("encoder", po::value<std::string>()->default_value("jpegli"), "Encoding backend: jpegli, turbo, or mozjpeg")
    ("I", po::value<std::string>(), "Input filename")
    ("O", po::value<std::string>(), "Output filename");

  po::positional_options_description pd;
  pd.add("I", 1).add("O", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pd).run(), vm);
  po::notify(vm);

  bool help = vm.count("h") || vm.count("?");
  if (help || !vm.count("I") || !vm.count("O")) {
    std::cout << "usage: jpegultrascan [switches] inputfile outputfile\n"
              << "JPEG lossless recompressor that tries all scan possibilities to minimize size\n"
              << "Switches:\n"
              << "  -s       Strip all extra markers\n"
              << "  -s some  Strip only non-comment markers\n"
              << "  -i       Allow multiple planes per DC scan, which may be incompatible with some software\n"
              << "  -i -1    Disallow single scan DC components, which may be incompatible with some software\n"
              << "  -j       Strip APP0 segment for 18-byte savings (generates non-compliant JPEG that may be incompatible with some software)\n"
              << "  -j 0-14  Strip APPn segment; specify multiple times to strip different segments\n"
              << "  -a       Use arithmetic coding (unsupported by most software)\n"
              << "  -b 0-13  Maximum number of bit splits to test (default: 3)\n"
              << "  -t [N]   Number of simultaneous processes (default: 1)\n"
              << "  -q       Suppress all output\n"
              << "  -v       Verbose output\n"
              << "  --dump-fail-scans FILE     Dump scan sets that fail encoding\n";
    return help ? 0 : 1;
  }

  const std::string fin = vm["I"].as<std::string>();
  const std::string fout = vm["O"].as<std::string>();
  const bool verbose = vm.count("v") > 0;
  const bool quiet = vm.count("q") > 0;
  const bool arithmetic = vm.count("a") > 0;
  const int incompat = vm["i"].as<int>();
  int max_bits = vm["b"].as<int>();
  if (vm.count("t")) {
    int requested = vm["t"].as<int>();
    if (requested <= 0) {
      unsigned detected = std::thread::hardware_concurrency();
      g_threads = detected ? static_cast<int>(detected) : 1;
    } else {
      g_threads = std::max(1, requested);
    }
  } else {
    unsigned detected = std::thread::hardware_concurrency();
    g_threads = detected > 0 ? static_cast<int>(detected) : 1;
  }

  std::string strip_opt = vm["s"].as<std::string>();
  StripMode strip_mode = StripMode::All;
  if (strip_opt == "some") strip_mode = StripMode::Comments;
  else if (strip_opt == "none") strip_mode = StripMode::None;

  std::set<int> strip_apps;
  if (vm.count("j")) {
    for (int app : vm["j"].as<std::vector<int>>()) {
      if (app >= 0 && app <= 14) strip_apps.insert(app);
    }
  }

  if (vm.count("dump-scan-sizes")) {
    g_scan_dump = std::make_unique<std::ofstream>(vm["dump-scan-sizes"].as<std::string>());
    if (!*g_scan_dump) {
      std::cerr << "Couldn't write scan sizes to " << vm["dump-scan-sizes"].as<std::string>() << "\n";
      return 1;
    }
  }
  g_dump_dht = vm.count("dump-dht") > 0;
  g_dump_segments = vm.count("dump-segments") > 0;
  if (vm.count("dump-symbols")) {
    g_symbol_dump = std::make_unique<std::ofstream>(vm["dump-symbols"].as<std::string>());
    if (!*g_symbol_dump) {
      std::cerr << "Couldn't write scan symbols to " << vm["dump-symbols"].as<std::string>() << "\n";
      return 1;
    }
  }
  if (vm.count("dump-scan-sets")) {
    g_scan_set_dump = std::make_unique<std::ofstream>(vm["dump-scan-sets"].as<std::string>());
    if (!*g_scan_set_dump) {
      std::cerr << "Couldn't write scan sets to " << vm["dump-scan-sets"].as<std::string>() << "\n";
      return 1;
    }
  }
  if (vm.count("dump-fail-scans")) {
    g_fail_scan_dump = std::make_unique<std::ofstream>(vm["dump-fail-scans"].as<std::string>());
    if (!*g_fail_scan_dump) {
      std::cerr << "Couldn't write fail scans to " << vm["dump-fail-scans"].as<std::string>() << "\n";
      return 1;
    }
  }

  std::string encoder = vm["encoder"].as<std::string>();
  if (encoder == "jpegli") {
    g_encoder_backend = EncoderBackend::Jpegli;
  } else if (encoder == "turbo") {
    g_encoder_backend = EncoderBackend::Turbo;
  } else if (encoder == "mozjpeg") {
    g_encoder_backend = EncoderBackend::Mozjpeg;
  } else {
    std::cerr << "Unknown encoder backend: " << encoder
              << " (expected jpegli, turbo, or mozjpeg)\n";
    return 1;
  }

  max_bits = std::clamp(max_bits, 0, 13);

  JpegInput input;
  input.data = read_file_bytes(fin);
  if (input.data.empty()) {
    std::cerr << "Couldn't read file " << fin << "\n";
    return 1;
  }

  input.dinfo.err = jpegli_std_error(&input.jerr);
  jpegli_create_decompress(&input.dinfo);
  jpegli_mem_src(&input.dinfo, input.data.data(), input.data.size());

  if (strip_mode != StripMode::None) {
    jpegli_save_markers(&input.dinfo, JPEG_COM, 0xFFFF);
    for (int app = 0; app <= 15; ++app) {
      jpegli_save_markers(&input.dinfo, JPEG_APP0 + app, 0xFFFF);
    }
  }

  jpegli_read_header(&input.dinfo, TRUE);
  input.coef_arrays = jpegli_read_coefficients(&input.dinfo);

  if (max_bits > 10 && input.dinfo.data_precision == 8) max_bits = 10;
  g_max_bits = max_bits;

  g_all_components.clear();
  for (int i = 0; i < input.dinfo.num_components; ++i) g_all_components.push_back(i);

  sizes.clear();
  optimal.clear();
  optimal2.clear();
  dc_cache.clear();
  scan_dht_info.clear();
  component_only_sizes.clear();
  ac_all = scan_set{};
  best_scans = scan_set{};
  best_scans.size = kInvalidSize;

  std::vector<scan_set> ac_scans(input.dinfo.num_components);
  std::vector<std::vector<std::vector<ScanLine>>> ac_scan_queues(input.dinfo.num_components);
  std::mutex io_lock;

  if (g_threads > 1 && input.dinfo.num_components > 1) {
    std::vector<std::thread> plane_threads;
    plane_threads.reserve(input.dinfo.num_components);
    for (int c = 0; c < input.dinfo.num_components; ++c) {
      plane_threads.emplace_back([&, c] {
        if (!quiet) {
          std::lock_guard<std::mutex> guard(io_lock);
          std::cout << "Calculating sizes of AC scans, plane " << c << "\n";
        }
        generate_ac_scans(c, g_all_components, g_max_bits, ac_scan_queues[c]);
      });
    }
    for (auto& th : plane_threads) th.join();
  } else {
    for (int c = 0; c < input.dinfo.num_components; ++c) {
      if (!quiet) std::cout << "Calculating sizes of AC scans, plane " << c << "\n";
      generate_ac_scans(c, g_all_components, g_max_bits, ac_scan_queues[c]);
    }
  }

  for (int c = 0; c < input.dinfo.num_components; ++c) {
    do_scans(input, ac_scan_queues[c], arithmetic, strip_mode, strip_apps, verbose, quiet);
  }

  if (g_threads > 1 && input.dinfo.num_components > 1) {
    std::vector<std::thread> split_threads;
    split_threads.reserve(input.dinfo.num_components);
    for (int c = 0; c < input.dinfo.num_components; ++c) {
      split_threads.emplace_back([&, c] {
        if (!quiet) {
          std::lock_guard<std::mutex> guard(io_lock);
          std::cout << "\nSearching for best AC scan, plane " << c << "\n";
        }
        ac_scans[c] = try_splits(std::vector<int>{c}, 1, DCTSIZE2 - 1, g_max_bits);
        if (verbose) {
          std::lock_guard<std::mutex> guard(io_lock);
          std::cout << "Best AC scan, plane " << c << ":\n";
          for (const auto& scan : ac_scans[c].all_scans) std::cout << scan_label(scan) << "\n";
          std::cout << "Size: " << ac_scans[c].size << "\n\n";
        }
      });
    }
    for (auto& th : split_threads) th.join();
  } else {
    for (int c = 0; c < input.dinfo.num_components; ++c) {
      if (!quiet) std::cout << "\nSearching for best AC scan, plane " << c << "\n";
      ac_scans[c] = try_splits(std::vector<int>{c}, 1, DCTSIZE2 - 1, g_max_bits);
      if (verbose) {
        std::cout << "Best AC scan, plane " << c << ":\n";
        for (const auto& scan : ac_scans[c].all_scans) std::cout << scan_label(scan) << "\n";
        std::cout << "Size: " << ac_scans[c].size << "\n\n";
      }
    }
  }

  for (const auto& set : ac_scans) ac_all += set;

  if (!quiet) std::cout << "Calculating sizes and searching for best DC scan\n";
  if (incompat > 0 && input.dinfo.num_components > 1) {
    std::vector<std::vector<int>> partition;
    for (int c = 0; c < input.dinfo.num_components; ++c) partition.push_back({c});
    partitionscans(input, 0, partition, arithmetic, strip_mode, strip_apps, verbose, quiet);
  } else {
    if (input.dinfo.num_components > 1 && incompat >= 0) {
      std::vector<std::vector<int>> partition;
      for (int c = 0; c < input.dinfo.num_components; ++c) partition.push_back({c});
      partitionscans(input, input.dinfo.num_components + 1, partition, arithmetic, strip_mode, strip_apps, verbose, quiet);
    }
    std::vector<std::vector<int>> partition = {g_all_components};
    partitionscans(input, 2, partition, arithmetic, strip_mode, strip_apps, verbose, quiet);
  }

  auto ordered_scans = order_scans(best_scans.all_scans);
  if (!quiet) {
    std::cout << "\nBest scans:\n";
    for (const auto& scan : ordered_scans) std::cout << scan_label(scan) << "\n";
  }

  std::vector<ScanLine> final_scans;
  if (!complete_scan_script(ordered_scans, input.dinfo.num_components, final_scans) ||
      count_ac_contexts(final_scans) > 252) {
    final_scans.clear();
    final_scans.push_back(make_scan(g_all_components, 0, DCTSIZE2 - 1, 0, 0));
  }

  std::vector<uint8_t> out_data;
  std::string out_encoder = "jpegli";
  std::vector<uint8_t> jpegli_data;
  if (!encode_with_scans_jpegli_bytes(input, final_scans, arithmetic, strip_mode, strip_apps,
                                      &jpegli_data, nullptr)) {
    jpegli_data.clear();
  }
  out_data = jpegli_data;
  if (ENABLE_TURBO) {
    std::vector<uint8_t> turbo_data;
    if (encode_with_scans_turbo(input.data, final_scans, arithmetic, &turbo_data)) {
      if (out_data.empty() || (!turbo_data.empty() && turbo_data.size() < out_data.size())) {
        out_data = std::move(turbo_data);
        out_encoder = "turbo";
      }
    }
  }
  if (!quiet) {
    std::cout << "Chosen encoder: " << out_encoder << "\n";
  }

  size_t insize = input.data.size();
  size_t outsize = out_data.size();
  if (!quiet) {
    std::cout << "Size: " << outsize << " bytes\n";
    std::cout << "Change: " << std::fixed << std::setprecision(6)
              << (static_cast<double>(outsize) / static_cast<double>(insize) - 1.0) * 100.0 << "%\n";
  }

  const std::vector<uint8_t>& to_write = (outsize && outsize <= insize) ? out_data : input.data;
  std::ofstream out(fout, std::ios::binary);
  if (!out) {
    std::cerr << "Couldn't write file " << fout << "\n";
    return 1;
  }
  out.write(reinterpret_cast<const char*>(to_write.data()), to_write.size());

  destroy_input(input);

  if (!quiet) {
    const auto program_end = std::chrono::steady_clock::now();
    const auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(program_end - program_start).count();
    const auto do_scans_ms = g_do_scans_wall_ns.load(std::memory_order_relaxed) / 1'000'000;
    const auto capacity_ms = g_do_scans_capacity_ns.load(std::memory_order_relaxed) / 1'000'000;
    const auto busy_ms = g_worker_busy_ns.load(std::memory_order_relaxed) / 1'000'000;
    const double util = capacity_ms ? (static_cast<double>(busy_ms) / static_cast<double>(capacity_ms)) * 100.0 : 0.0;
    const auto other_ms = total_ms - do_scans_ms;
    std::cout << "Timing: total=" << total_ms << "ms do_scans=" << do_scans_ms
              << "ms other=" << other_ms << "ms util=" << std::fixed
              << std::setprecision(1) << util << "%\n";
  }
  return 0;
}
