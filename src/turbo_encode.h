#pragma once

#include <cstdint>
#include <vector>

#include "scan_types.h"

bool encode_with_scans_turbo(const std::vector<uint8_t>& input_data,
                             const std::vector<ScanLine>& scans,
                             bool arithmetic,
                             std::vector<uint8_t>* out_bytes);
