#pragma once

#include <cstdint>
#include <vector>

#include "scan_types.h"
#include "jpeglib.h"

bool encode_with_scans_mozjpeg(const std::vector<uint8_t>& input_data,
                               const std::vector<ScanLine>& scans,
                               bool arithmetic,
                               std::vector<uint8_t>* out_bytes);

bool encode_with_scans_mozjpeg_from_coeffs(j_decompress_ptr dinfo,
                                           jvirt_barray_ptr* coef_arrays,
                                           const std::vector<ScanLine>& scans,
                                           bool arithmetic,
                                           bool minimal_scan,
                                           std::vector<uint8_t>* out_bytes);
