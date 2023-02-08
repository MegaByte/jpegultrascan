/** 
 * @file      jpegultrascan.cpp 
 * @brief     JPEG lossless recompressor that tries all scan possibilities to minimize size
 * @author    Aaron Kaluszka
 * @version   2.0.0
 * @date      2022-12-22
 * @copyright Copyright 2015-2023 Aaron Kaluszka
 * \n \n
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * \n \n
 *     http://www.apache.org/licenses/LICENSE-2.0
 * \n \n
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
#include <boost/program_options.hpp>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>


#include "mozjpeg/jinclude.h"
#include "mozjpeg/jpeglib.h"
#include "mozjpeg/jpegint.h"
#include "mozjpeg/transupp.h"
#include "mozjpeg/jcmaster.h"

const uint16_t MAX_HT = 65535;
const uint16_t MAX_COEFFS = 63;

#define JPEG_APP1       0xE1    /* APP1 marker code  (EXIF) */
#define JPEG_APP14       0xEE    /* APP14 marker code  (Adobe) */
#define EXIF_HEADER_LENGTH 6

/* Map Exif orientation values to correct JXFORM_CODE */
static JXFORM_CODE orient_jxform[9] = {
  JXFORM_NONE,
  JXFORM_NONE,
  JXFORM_FLIP_H,
  JXFORM_ROT_180,
  JXFORM_FLIP_V,
  JXFORM_TRANSPOSE,
  JXFORM_ROT_90,
  JXFORM_TRANSVERSE,
  JXFORM_ROT_270
};


struct component {
  uint8_t component;
  uint8_t scale_v: 4;
  uint8_t scale_h: 4;
  uint8_t qt: 4;
};

struct scan_component {
  uint8_t component;
  uint8_t HTDC: 4;
  uint8_t HTAC: 4;
};

struct scan {
  uint8_t plane_start;
  uint8_t plane_end;
  uint8_t bit_start;
  uint8_t bit_end;
};

class JPEG {
public:
  uint16_t width = 0;
  uint16_t height = 0;
  uint8_t bits_sample = 8;

  std::vector<uint8_t*> segments;
  std::vector<std::vector<std::vector<int16_t> > > fdimage; 
  std::vector<component> components;
  std::vector<std::vector<uint8_t> > QT;
  std::vector<std::array<uint8_t, MAX_HT> > HT;

  void dqt_decode(uint8_t* i) {
    uint16_t size = ((uint16_t)(*(i + 1)) << 8) + (*(i + 2));
    uint16_t extended = (*(i + 3)) >> 4;
    uint16_t id = (*(i + 3)) & 0xF;
    if (QT.size() <= id) {
      QT.resize(id + 1);
      QT[id].resize(size);
      // copy data
    }
    std::cout << "Size " << size << std::endl;
    std::cout << "Ex:" << extended << std::endl;
    std::cout << "ID:" << id << std::endl;
  }

  void dht_decode(uint8_t* i) {
    uint16_t size = ((uint16_t)(*(i + 1)) << 8)+(*(i + 2));
    i += 3;
    uint16_t number = *i >> 4;
    uint16_t acdc = (*i & 0xF) >> 3;
    if (HT.size() <= number) {
      HT.resize(number + 1);
      for (size_t i = 0; i < 65536; ++i) {
        HT[number][i] = 0xFF;
      }
    }
    std::cout << "Size " << size << std::endl;
    std::cout << "Number " << number << std::endl;
    std::cout << (acdc ? "AC" : "DC") << std::endl;
    std::cout << "Symbol count";
    uint16_t total = 0;
    ++i;
    for (size_t j = 0; j < 16; ++j) {
      uint16_t count = (uint16_t)*(i + j);
      std::cout << ' ' << count;
    }
    std::cout << std::endl;
    i += 16;
    uint16_t code = 0;
    for (size_t j = 0; j < 16; ++j) {
      uint16_t jcount = *(i - 16 + j);
      for (size_t k = 0; k < jcount; ++k) {
        ++code;
        uint8_t value = *(i + total + k);
        std::cout << code << ' ' << (uint16_t)value << std::endl;
        HT[number][code] = value;
      }
      total += jcount;
      code <<= 1;
    }

    std::cout << std::endl;
  }

  void sof_decode(uint8_t* i) {
    uint16_t size = ((uint16_t)(*(i + 1)) << 8)+(*(i + 2));
    i += 3;
    bits_sample = *i;
    width = ((uint16_t)(*(i + 1)) << 8) + *(i + 2);
    height = ((uint16_t)(*(i + 3)) << 8) + *(i + 4);
    uint16_t num_components = *(i + 5);
    std::cout << "Size " << size << std::endl;
    std::cout << "Bits/sample " << bits_sample << std::endl;
    std::cout << "Width " << width << std::endl;
    std::cout << "Height " << height << std::endl;
    std::cout << "Components " << num_components << std::endl;
    components.resize(num_components);
    fdimage.resize(num_components);
    for (size_t j = 0; j < num_components; ++j) {
      components[j].component = *(i + 6 + j * 3);
      components[j].scale_v = (*(i + 6 + j * 3 + 1) >> 4);
      components[j].scale_h = (*(i + 6 + j * 3 + 1) & 0xF);
      components[j].qt = *(i + 6 + j * 3 + 2);
      std::cout << "Component " << components[j].component << std::endl;
      std::cout << "SF V " << components[j].scale_v << std::endl;
      std::cout << "SF H " << components[j].scale_h << std::endl;
      std::cout << "QT " << components[j].qt << std::endl;
      fdimage[j].resize(MAX_COEFFS + 1);
      for (size_t k = 0; k <= MAX_COEFFS; ++k) {
        fdimage[j][k].assign(width * height / (MAX_COEFFS + 1) / components[j].scale_v / components[j].scale_h, 0);
      }
    }
  }

  uint8_t* sos_decode(uint8_t* i) {
    uint16_t size = ((uint16_t)(*(i + 1)) << 8) + (*(i + 2));
    i += 3;
    uint8_t num_components = *i;
    std::cout << "Size " << size << std::endl;
    std::cout << "Components " << uint16_t(num_components) << std::endl;
    std::vector<scan_component> scan_components(num_components);
    for (size_t j = 0; j < num_components; ++j) {
      scan_components[j].component = *(i + 1 + j * 2);
      scan_components[j].HTAC = (*(i + 1 + j * 2 + 1) >> 4);
      scan_components[j].HTDC = (*(i + 1 + j * 2 + 1) & 0xF);
      std::cout << "Component " << scan_components[j].component << std::endl;
      std::cout << "HT AC " << scan_components[j].HTAC << std::endl;
      std::cout << "HT DC " << scan_components[j].HTDC << std::endl;
    }
    i += 1 + num_components * 2;
    uint8_t start_of_selection = *i;
    uint8_t end_of_selection =  *(i + 1);
    uint8_t bit_position =  *(i + 2);
    std::cout << "Start of selection " << uint16_t(start_of_selection) << std::endl;
    std::cout << "End of selection " << uint16_t(end_of_selection) << std::endl;
    std::cout << "Bit position " << uint16_t(bit_position) << std::endl;
    uint8_t num_planes = (end_of_selection - start_of_selection + 1);
    size_t blocksize = size_t(num_planes) * num_components;

    uint8_t bit_start = 0;
    uint8_t bit_end = 0;
    bool skip_next = false;
  // scan decode
    i += 3;
    uint16_t index;
    uint8_t cur_component = 0;
    uint8_t cur_plane = start_of_selection;
    int16_t DC = 0;
    size_t k = 0;
    size_t block = 0;
    while (1) {
      ++bit_end;
      if (bit_end == 1) {
        if (*i == 0xFF) {
          if (*(i + 1) != 0) return i;
        }
        index = *i & ((2 << bit_start) - 1);
      }
      else if (bit_end == 9 || bit_end == 17) {
        if (*i == 0xFF) {
          if (*(i + 1) != 0) return i;
          ++i;
        }
        index = (index << 8) | *(++i);
      }

      uint16_t indexX = index >> (8 - (bit_end & 7));
      uint8_t val = HT[cur_plane == 0 ? scan_components[cur_component].HTDC : scan_components[cur_component].HTAC][indexX];
      if (val == 0xFF ) {
        // not a valid bit sequence, get another bit
        continue;
      }
      bit_end &= 7;
      bit_start = bit_end;
      if (val == 0) {
        // next block
        k += k % blocksize;
        cur_component = (k % (num_planes * num_components)) / num_planes;
        cur_plane = (k % (num_planes * num_components)) % num_planes;
        std::cout << "next" << std::endl;
        continue;
      }
      uint8_t zeros;
      uint8_t nextsize;
      if (cur_plane == 0) {
        nextsize = val;
      } else {
        zeros = val >> 4;
        nextsize = val & 0xF;
        k += zeros;
        if (nextsize == 0) ++k;
      }
      cur_component = (k % (num_planes * num_components)) / num_planes;
      cur_plane = (k % (num_planes * num_components)) % num_planes;

      bit_end += nextsize;
      
      index = *i & ((2 << bit_start) - 1);

      if (bit_end >= 8) {
        ++i;
        if (*i == 0xFF) {
          if (*(i + 1) != 0) return i;
        }
        index = (index << 8) | *(++i);
      }
      if (bit_end >= 16) {
        ++i;
        if (*i == 0xFF) {
          if (*(i + 1) != 0) return i;
        }
        index = (index << 8) | *(++i);
      }
      block = k / (num_planes * num_components);
      std::cout << k << ' ' << indexX << ' ' << int16_t(val) << ' ' << uint16_t(block) << ' ' << uint16_t(cur_component) << ' ' << uint16_t(cur_plane + start_of_selection) << ' ' << blocksize << std::endl;
      if (cur_plane == 0) { //DC
        int16_t dc_value = index >> (nextsize - 1) ? index : ~(int16_t)index;
        fdimage[scan_components[cur_component].component][0][block] |= dc_value + DC;
        DC = fdimage[scan_components[cur_component].component][0][block];
      } else {
        fdimage[scan_components[cur_component].component][cur_plane + start_of_selection][block] |= index >> (nextsize - 1) ? index : ~(int16_t)index;
      }
      ++k;
      cur_component = (k % (num_planes * num_components)) / num_planes;
      cur_plane = (k % (num_planes * num_components)) % num_planes;
    }
    
    return i;
  }

#define JPEG_SOI 0xD8
#define JPEG_DQT 0xDB
#define JPEG_APP0 0xE0
#define JPEG_SOF0 0xC0
#define JPEG_SOF1 0xC1
#define JPEG_SOF2 0xC2
#define JPEG_SOF9 0xC9
#define JPEG_DHT 0xC4
#define JPEG_SOS 0xDA
#define JPEG_EOI 0xD9
#define JPEG_COM 0xFE
#define JPEG_TAG 0xFF



  void decode(uint8_t* data, uint8_t* data_end) {
    for (uint8_t* i = data; i < data_end; ++i) {
      if (*i == JPEG_TAG) {
        if (*(i + 1) == 0) {
          ++i;
          continue;
        }
        segments.push_back(i);

        std::cout << (uint64_t)(i - data) << ' ';
        switch (*(++i)) {
          case JPEG_SOI:
            std::cout << "SOI";
          break;
          case JPEG_APP0:
            std::cout << "APP0";
          break;
          case JPEG_APP1:
            std::cout << "APP1";
          break;
          case JPEG_DQT:
            std::cout << "DQT" << std::endl;
            dqt_decode(i);
          break;
          case JPEG_SOF0:
          case JPEG_SOF1:
          case JPEG_SOF2:
          case JPEG_SOF9:
            std::cout << "SOF";
            std::cout << ((uint16_t)(*i) & 3) << std::endl;
            sof_decode(i);
          break;
          case JPEG_DHT:
            std::cout << "DHT" << std::endl;
            dht_decode(i);
          break;
          case JPEG_SOS:
            std::cout << "SOS" << std::endl;
            i = sos_decode(i);
            continue;
          break;
          case JPEG_EOI:
            std::cout << "EOI";
          break;
          case JPEG_COM:
            std::cout << "COM";
          break;
          default:
            std::cout << (uint16_t)(*i);
          break;
        }
        if (*i == JPEG_SOI || *i == JPEG_EOI || *i == JPEG_SOS) {
          std::cout << std::endl;
          continue;
        }
        uint16_t size = ((uint16_t)(*(i + 1))<<8)+(*(i + 2));
  //      std::cout << ' ' << size;
        i += size;
        std::cout << std::endl;
      }
    }
  }

  void optimize_huffman() {
// just copy for now
  };

  void encode(bool app0, bool app1, bool com) {
    optimize_huffman();



    // SOI
    // APP0
    // APP1
    // DQT
    // SOF
    // DHT
    // SOS
    // COM
    // EOI
  }

};

struct scan_set {
  uint32_t size;
  std::vector<jpeg_scan_info> all_scans;
  scan_set operator+(const scan_set& r) {
    scan_set ss = *this;
    ss.size += r.size;
    ss.all_scans.insert(std::end(ss.all_scans), std::cbegin(r.all_scans), std::cend(r.all_scans));
    return ss;
  }
  scan_set& operator+=(const scan_set& r) {
    size += r.size;
    all_scans.insert(std::end(all_scans), std::cbegin(r.all_scans), std::cend(r.all_scans));
    return *this;
  }
};

void queue_scans(std::vector<jpeg_scan_info>& scan_queue, const std::vector<jpeg_scan_info>& scans) {
  scan_queue.insert(std::end(scan_queue), std::cbegin(scans), std::cend(scans));
}

void init_jpeg_scan_info(jpeg_scan_info& info, const std::vector<int>& components, int Ss, int Se, int Ah, int Al) {
  info.comps_in_scan = components.size();
  for (size_t i = 0; i < components.size(); ++i) info.component_index[i] = components[i];
  info.Ss = Ss;
  info.Se = Se;
  info.Ah = Ah;
  info.Al = Al;
}

void init_jpeg_scan_info(jpeg_scan_info& info, int component, int Ss, int Se, int Ah, int Al) {
  info.comps_in_scan = 1;
  info.component_index[0] = component;
  info.Ss = Ss;
  info.Se = Se;
  info.Ah = Ah;
  info.Al = Al;
}

std::vector<jpeg_scan_info> generate_bit_scans(std::vector<int>& plane, int a, int b, int num_bits) {
  std::vector<jpeg_scan_info> scans(num_bits + 1);
  init_jpeg_scan_info(scans[0], plane, a, b, 0, num_bits);
  for (auto i = num_bits; i > 0; --i) init_jpeg_scan_info(scans[num_bits - i + 1], plane, a, b, i, i - 1);
  return scans;
}

std::vector<jpeg_scan_info> generate_bit_scans(int plane, int a, int b, int num_bits) {
  std::vector<jpeg_scan_info> scans(num_bits + 1);
  init_jpeg_scan_info(scans[0], plane, a, b, 0, num_bits);
  for (auto i = num_bits; i > 0; --i) init_jpeg_scan_info(scans[num_bits - i + 1], plane, a, b, i, i - 1);
  return scans;
}

void generate_dc_scans(std::vector<int>& plane, int max_bits, std::vector<jpeg_scan_info>& scan_queue) {
  for (auto i = 0; i <= max_bits; ++i) {
    auto scans = generate_bit_scans(plane, 0, 0, i);
    queue_scans(scan_queue, scans);
  }
}

void generate_ac_scans(int plane, int max_bits, std::vector<jpeg_scan_info>& scan_queue) {
  for (auto num_coeffs = 0; num_coeffs < DCTSIZE2; ++num_coeffs) {
    auto offset_max = std::min(num_coeffs + 1, DCTSIZE2 - 1 - num_coeffs);
    for (auto offset = 0; offset <= offset_max; ++offset) {
      for (auto num_bits = 0; num_bits <= max_bits; ++num_bits) {
        std::vector<jpeg_scan_info> scans;
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
        queue_scans(scan_queue, scans);
      }
    }
  }
}

static std::unordered_map<std::string, scan_set> optimal;
static std::unordered_map<std::string, scan_set> optimal2;
static std::unordered_map<std::string, int> sizes;

scan_set choose_best(const scan_set& left, const scan_set& right) {
  if (!left.size && !right.size) return left;
  if (!right.size || (left.size && left.size < right.size)) return left;
  if (!left.size || (right.size && right.size < left.size) || left.all_scans.size() > right.all_scans.size()) return right;
  return left;
}

scan_set try_bit_splits(std::vector<int>& plane, int a, int b, int c, int d) {
  std::stringstream key;
  for (auto i = 0; i < plane.size(); ++i) key << plane[i] << ' ';
  key << a << ' ' << b << ' ' << c << ' ' << d;
  if (!optimal2.contains(key.str())) {
    scan_set min_scans;
    jpeg_scan_info scan;
    init_jpeg_scan_info(scan, plane, a, b, c, d);
    min_scans.all_scans.push_back(scan);
    min_scans.size = sizes[key.str()];
    for (auto i = a; i < b; ++i) {
      auto a_scans = try_bit_splits(plane, a, i, c, d);
      auto b_scans = try_bit_splits(plane, i + 1, b, c, d);
      min_scans = choose_best(a_scans + b_scans, min_scans);
    }
    optimal2[key.str()] = min_scans;
  }
  return optimal2[key.str()];
}

scan_set try_splits(std::vector<int>& plane, int a, int b, int max_bits) {
  std::stringstream key;
  for (auto i = 0; i < plane.size(); ++i) key << plane[i] << ' ';
  key << a << ' ' << b;
  if (!optimal.contains(key.str())) {
    scan_set min_scans{0};
    for (auto i = 0; i < max_bits; ++i) {
      auto all_scans = try_bit_splits(plane, a, b, 0, i);
      for (auto j = i; i > 0; --j) {
        all_scans += try_bit_splits(plane, a, b, i, i - 1);
      }
      min_scans = choose_best(all_scans, min_scans);
    }
    for (auto i = a; i < b; ++i) {
      auto a_scans = try_splits(plane, a, i, max_bits);
      auto b_scans = try_splits(plane, i, b, max_bits);
      min_scans = choose_best(a_scans + b_scans, min_scans);
    }
  }
  return optimal[key.str()];
}

static size_t jcopy_markers_execute_s (j_decompress_ptr srcinfo, j_compress_ptr dstinfo)
{
  size_t size = 0;
  for (jpeg_saved_marker_ptr marker = srcinfo->marker_list; marker; marker = marker->next) {
    if (dstinfo->write_JFIF_header &&
        marker->marker == JPEG_APP1 &&
        marker->data_length >= 5 &&
        (marker->data[0]) == 'J' &&
        (marker->data[1]) == 'F' &&
        (marker->data[2]) == 'I' &&
        (marker->data[3]) == 'F' &&
        (marker->data[4]) == 0)
      continue;                 /* reject duplicate JFIF */
    if (dstinfo->write_Adobe_marker &&
        marker->marker == JPEG_APP14 &&
        marker->data_length >= 5 &&
        (marker->data[0]) == 'A' &&
        (marker->data[1]) == 'd' &&
        (marker->data[2]) == 'o' &&
        (marker->data[3]) == 'b' &&
        (marker->data[4]) == 'e')
      continue;                 /* reject duplicate Adobe */
    jpeg_write_marker(dstinfo, marker->marker, marker->data, marker->data_length);
    if (marker->marker == JPEG_COM || (marker->marker >= JPEG_APP0 && marker->marker <= JPEG_APP0 + 15)){
      size += marker->data_length;
    }
  }
  return size;
}

/* Get Exif image orientation. Copied from adjust_exif_parameters. */
LOCAL(unsigned int)
get_exif_orientation (JOCTET* data, unsigned int length) {
  boolean is_motorola; /* Flag for byte order */
  unsigned int number_of_tags, tagnum;
  unsigned int offset;
  JOCTET *word[2], *dword[4];
  if (length < 12) return 0; /* Length of an IFD entry */

  /* Discover byte order */
  if (data[0] == 0x49 && data[1] == 0x49) {
    word[0] = data + 1;
    word[1] = data + 0;
    dword[0] = data + 3;
    dword[1] = data + 2;
    dword[2] = data + 1;
    dword[3] = data + 0;
  } else if (data[0] == 0x4D && data[1] == 0x4D) {
    word[0] = data + 0;
    word[1] = data + 1;
    dword[0] = data + 0;
    dword[1] = data + 1;
    dword[2] = data + 2;
    dword[3] = data + 3;
  } else return 0;

  /* Check Tag Mark */
  if (word[0][2] != 0 || word[1][2] != 0x2A) return 0;

  /* Get first IFD offset (offset to IFD0) */
  if (dword[0][4] != 0 || dword[1][4] != 0) return 0;
  offset = dword[2][4];
  offset <<= 8;
  offset += dword[3][4];
  if (offset > length - 2) return 0; /* check end of data segment */

  /* Get the number of directory entries contained in this IFD */
  number_of_tags = (word[0][offset] << 8) + word[1][offset];
  if (number_of_tags == 0) return 0;
  offset += 2;

  /* Search for Orientation Tag in this IFD */
  do {
    if (offset > length - 12) return 0; /* check end of data segment */
    const int ORIENTATION_TAG = 0x0112;
    if ((word[0][offset] << 8) + word[1][offset] == ORIENTATION_TAG) return word[0][offset + 8];
    offset += 12;
  } while (--number_of_tags);
  return 0;
}

namespace po = boost::program_options;
static const char* progname;

int print_error(const char* format, ...) {
    int result;
    va_list args;
    va_start(args, format);
    fprintf(stderr, "%s: ", progname);
    result = vprintf(format, args);
    va_end(args);
    return result;
}

METHODDEF(void)
output_message (j_common_ptr cinfo)
{
  char buffer[JMSG_LENGTH_MAX];

  /* Create the message */
  (*cinfo->err->format_message) (cinfo, buffer);

  /* Send it to stderr, adding a newline */
  print_error("%s\n", cinfo->err->addon_message_table[0], buffer);
}

long long filesize (const char* file) {
    struct stat stats;
    return stat(file, &stats) == 0 ? stats.st_size : -1;
}

void select_scans(j_compress_ptr cinfo, int next_scan_number) {
  std::cout << "SELECTION" << std::endl;
  my_master_ptr master = (my_master_ptr) cinfo->master;
  for (int i = 0; i < cinfo->num_scans; ++i) {
    std::cout << master->scan_size[i] << ' ';
  }
  std::cout << std::endl;
}

void finish_pass(j_compress_ptr cinfo) {
  my_master_ptr master = (my_master_ptr) cinfo->master;
  std::cout << "FINISHING " << master->pass_type << std::endl;
  /* The entropy coder always needs an end-of-pass call,
   * either to analyze statistics or to flush its output buffer.
   */
  (*cinfo->entropy->finish_pass) (cinfo);

  /* Update state for next pass */
  switch (master->pass_type) {
  case main_pass:
    /* next pass is either output of scan 0 (after optimization)
     * or output of scan 1 (if no optimization).
     */
    if (cinfo->master->trellis_quant)
      master->pass_type = trellis_pass;
    else {
    master->pass_type = output_pass;
    if (! cinfo->optimize_coding)
      master->scan_number++;
    }
    break;
  case huff_opt_pass:
    /* next pass is always output of current scan */
    master->pass_type = (master->pass_number < master->pass_number_scan_opt_base-1) ? trellis_pass : output_pass;
    break;
  case output_pass:
    /* next pass is either optimization or output of next scan */
    if (cinfo->optimize_coding)
      master->pass_type = huff_opt_pass;

      (*cinfo->dest->term_destination)(cinfo);
      cinfo->dest = master->saved_dest;
      select_scans(cinfo, master->scan_number + 1);

    master->scan_number++;
    break;
  case trellis_pass:
    if (cinfo->optimize_coding)
      master->pass_type = huff_opt_pass;
    else
      master->pass_type = (master->pass_number < master->pass_number_scan_opt_base-1) ? trellis_pass : output_pass;
      
    if ((master->pass_number + 1) %
        (cinfo->num_components * (cinfo->master->use_scans_in_trellis ? 4 : 2)) == 0 &&
        cinfo->master->trellis_q_opt) {
      int i, j;

      for (i = 0; i < NUM_QUANT_TBLS; i++) {
        for (j = 1; j < DCTSIZE2; j++) {
          if (cinfo->master->norm_coef[i][j] != 0.0) {
            int q = (int)(cinfo->master->norm_src[i][j] /
                          cinfo->master->norm_coef[i][j] + 0.5);
            if (q > 254) q = 254;
            if (q < 1) q = 1;
            cinfo->quant_tbl_ptrs[i]->quantval[j] = q;
          }
        }
      }
    }
    break;
  }

  master->pass_number++;
}

int main(int argc, char* argv[]) {
  bool strip = true;
  int autorotate = 2;
  const char* fin;
  const char* fout;

  po::options_description desc("Switches");
  desc.add_options()
    ("s", po::value<std::string>()->default_value("all"), "Strip extra markers")
    ("i", po::value<int>()->default_value(1), "Allow multiple planes per DC scan, which may be incompatible with some software")
    ("j", "Strip APP0 segment for 18-byte savings (generates non-compliant JPEG that may be incompatible with some software)")
    ("a", "Use arithmetic coding (unsupported by most software)")
    ("b", po::value<int>()->default_value(3), "Maximum number of bit splits to test (default: 3)")
    ("t", po::value<int>()->default_value(8), "Number of simultaneous processes (default: 8 if specified, 1 otherwise)")
    ("q", "Suppress all output")
    ("v", "Verbose output")
    ("h", "Help")
    ("?", "Help")
    ("I", po::value<std::string>(), "Input filename")
    ("O", po::value<std::string>(), "Output filename")
    ;

  po::positional_options_description pd;
  pd.add("I", 1).add("O", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pd).run(), vm);
  po::notify(vm);

  std::cout << "usage: jpegultrascan.pl [switches] inputfile outputfile" << std::endl
    << "JPEG lossless recompressor that tries all scan possibilities to minimize size" << std::endl;

  progname = argv[0];
  fin = vm["I"].as<std::string>().data();
  fout = vm["O"].as<std::string>().data();
  bool verbose = vm.count("v");
  bool quiet = vm.count("q");
  bool help = vm.count("h") || vm.count("?");
  bool arithmetic = vm.count("a");
  bool strip_app0 = vm.count("j");
  int max_bits = vm["b"].as<int>();

  int fd;
  struct stat sb;
  fd = open(argv[1], O_RDWR); 
  if (fd < 0) {
    fprintf(stderr, "Could not open %s\n", argv[1]);
    return 1;
  }
  fstat(fd, &sb);
  uint8_t* data;
  uint8_t* data_end;
  data = (uint8_t*) mmap(NULL, sb.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
  if (data == MAP_FAILED) exit(1);
 // file.rindex = file.windex = file.lindex = file.data_start;
  data_end = data + sb.st_size;
  JPEG j;
  j.decode(data, data_end);
  j.encode(false, false, false);
  return 0;
}

/*
  struct jpeg_decompress_struct srcinfo;
  struct jpeg_compress_struct dstinfo;
  struct jpeg_error_mgr jsrcerr, jdsterr;
  jpeg_transform_info transformoption; // image transformation options
  FILE* fp;
  unsigned char *outbuffer = 0;
  unsigned long outsize = 0;
  size_t extrasize = 0;
  unsigned char copy_exif = 0;

  // Initialize the JPEG decompression object with default error handling.
  srcinfo.err = jpeg_std_error(&jsrcerr);
  srcinfo.err->output_message = output_message;
  const char* addon = fin;
  srcinfo.err->addon_message_table = &addon;
  jpeg_create_decompress(&srcinfo);

  // Initialize the JPEG compression object with default error handling.
  dstinfo.err = jpeg_std_error(&jdsterr);
  jpeg_create_compress(&dstinfo);
 
  // Open the input file.
  if (!(fp = fopen(fin, "rb"))) {
    print_error("can't open %s for reading\n", progname, fin);
    return 2;
  }

  long long insize = filesize(fin);
  if (insize < 0){
    print_error("can't read from %s\n", progname, fin);
    return 2;
  }

  unsigned char* inbuffer = (unsigned char*) malloc(insize);
  if (!inbuffer) {
    print_error("memory allocation failure\n");
    exit(1);
  }

  if (fread(inbuffer, 1, insize, fp) < (size_t) insize) {
    print_error("can't read from %s\n", progname, fin);
  }
  fclose(fp);

  jpeg_mem_src(&srcinfo, inbuffer, insize);

  // Enable saving of extra markers that we want to copy
  if (!strip) {
    jcopy_markers_setup(&srcinfo, JCOPYOPT_ALL);
  } else if (autorotate > 0) {
    jpeg_save_markers(&srcinfo, JPEG_APP1, 0xFFFF);
  }

  // Read file header
  jpeg_read_header(&srcinfo, 1);

  // Determine orientation from Exif
  transformoption.transform = JXFORM_NONE;

  const unsigned char EXIF_HEADER[EXIF_HEADER_LENGTH] = {'E', 'x', 'i', 'f', '\0', '\0'};
  if (autorotate > 0 && srcinfo.marker_list != NULL &&
      srcinfo.marker_list->marker == JPEG_APP1 &&
      srcinfo.marker_list->data_length >= EXIF_HEADER_LENGTH &&
      memcmp(srcinfo.marker_list->data, EXIF_HEADER, EXIF_HEADER_LENGTH) == 0) {
    unsigned int orientation = get_exif_orientation(srcinfo.marker_list->data + EXIF_HEADER_LENGTH, srcinfo.marker_list->data_length - EXIF_HEADER_LENGTH);
    // Setup transform options for auto-rotate
    if (orientation > 1 && orientation <= 8) {
      transformoption.transform = orient_jxform[orientation];
      transformoption.perfect = autorotate > 1;
      transformoption.trim = TRUE;
      transformoption.force_grayscale = FALSE;
      transformoption.crop = FALSE;
      transformoption.slow_hflip = FALSE;
      // If perfect requested but not possible, show warning and do not transform
      if (!jtransform_request_workspace(&srcinfo, &transformoption)) {
        print_error("%s can't be transformed perfectly\n", progname, fin);
        transformoption.transform = JXFORM_NONE;
        copy_exif = 1;
      }
    }
  }

  // Read source file as DCT coefficients
  jvirt_barray_ptr * src_coef_arrays = jpeg_read_coefficients(&srcinfo);

  // Initialize destination compression parameters from source values
  jpeg_copy_critical_parameters(&srcinfo, &dstinfo);

  // Adjust destination parameters if required by transform options;
  //  also find out which set of coefficient arrays will hold the output.
  jvirt_barray_ptr * dst_coef_arrays = src_coef_arrays;
  if (transformoption.transform != JXFORM_NONE) {
    dst_coef_arrays = jtransform_adjust_parameters(&srcinfo, &dstinfo, src_coef_arrays, &transformoption);
  }
  dstinfo.write_JFIF_header = strip_app0;

  // Adjust default compression parameters by re-parsing the options
  dstinfo.optimize_coding = !arithmetic;
  dstinfo.arith_code = arithmetic;

  // Specify data destination for compression
  jpeg_mem_dest(&dstinfo, &outbuffer, &outsize);

  std::vector<jpeg_scan_info> scan_queue;

  std::cout << "Calculating sizes and searching for best DC scan" << std::endl;

  std::vector components;
  for (size_t i = 0; i < dstinfo.num_components; ++i) components.push_back(i);
  generate_dc_scans(components, max_bits, scan_queue);
  for (size_t i = 0; i < dstinfo.num_components; ++i) {
    components.clear();
    components.push_back(i);
    generate_dc_scans(components, max_bits, scan_queue);
  }


//  if ($incompat > 0 && $planes > 1) {
//    partitionscans 0, 0 .. $planes;
//  } else {
//    if $planes > 0 && $incompat >= 0;
//      partitionscans $planes + 1, 0 .. $planes
//    partitionscans 2, join ' ', 0 .. $planes;
//  }

  for (size_t i = 0; i < dstinfo.num_components; ++i) {
    std::cout << "Calculating sizes of AC scans, plane " << i << std::endl;
    generate_ac_scans(i, max_bits, scan_queue);

    // doscans
    // std::cout << "Searching for best AC scan, plane " << i << std::endl;
    // auto ac_all_scans = try_splits(i, 1, DCTSIZE2 - 1);
    // if (verbose) std::cout << "Best AC scan, plane " << i << ":" << std::endl;
    // ac_all_scans[i].all_scans[i]
    //std::cout << "Size: " << ac_all_scans.size << std::endl;
  }

  jpeg_c_set_bool_param(&dstinfo, JBOOLEAN_OPTIMIZE_SCANS, TRUE);

  // Start compressor (note no image data is actually written here)
  jpeg_write_coefficients(&dstinfo, dst_coef_arrays);

  // Write list of all scans

  auto scanptr = (jpeg_scan_info *) (dstinfo.mem->alloc_small) ((j_common_ptr) &dstinfo, JPOOL_IMAGE, scan_queue.size() * sizeof(jpeg_scan_info));
  memcpy(scanptr, scan_queue.data(), scan_queue.size() * sizeof(jpeg_scan_info));

for (size_t i = 0; i < scan_queue.size(); ++i) {
      std::cout << scanptr[i].Ss << ' ' << scanptr[i].Se << ' ' << scanptr[i].Ah << ' ' << scanptr[i].Al << std::endl;
}

  dstinfo.scan_info = scanptr;
  dstinfo.num_scans = scan_queue.size();
  jpeg_c_set_bool_param(&dstinfo, JBOOLEAN_OPTIMIZE_SCANS, FALSE);

  // Copy to the output file any extra markers that we want to preserve
  if (!strip || copy_exif) {
    extrasize = jcopy_markers_execute_s(&srcinfo, &dstinfo);
  }

  // Execute image transformation, if any
  if (transformoption.transform != JXFORM_NONE) {
    jtransform_execute_transformation(&srcinfo, &dstinfo, src_coef_arrays, &transformoption);
  }

  dstinfo.master->finish_pass = finish_pass;
  // Finish compression and release memory
  jpeg_finish_compress(&dstinfo);
  free(inbuffer);

  bool x = insize < outsize;

  if (outsize < insize){
    // Open the output file.
    if (!(fp = fopen(fout, "wb"))) {
      print_error("can't open %s for writing\n", progname, fout);
      free(outbuffer);
      return 2;
    }

    // Write new file.
    if (fwrite(outbuffer, 1, outsize, fp) < outsize) {
      print_error("can't write to %s\n", progname, fout);
    }
    fclose(fp);
  }

  jpeg_destroy_compress(&dstinfo);
  jpeg_finish_decompress(&srcinfo);
  jpeg_destroy_decompress(&srcinfo);
  free(outbuffer);
//  (*stripped_outsize) = x ? insize : outsize - extrasize;
  return x;
}
*/
