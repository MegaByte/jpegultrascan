#pragma once

#include <vector>

struct ScanLine {
  std::vector<int> components;
  int Ss = 0;
  int Se = 0;
  int Ah = 0;
  int Al = 0;
};
