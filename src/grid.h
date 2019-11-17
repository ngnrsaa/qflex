#ifndef GRID__H
#define GRID__H

#include <algorithm>
#include <fstream>
#include <vector>

#include "errors.h"
#include "global.h"

namespace qflex {

struct QflexGrid {
  int I{0}, J{0};
  std::vector<std::vector<int>> qubits_off;
  void clear();
  void load(std::istream& istream);
  void load(std::istream&& istream);
  void load(const std::string& filename);
};

}  // namespace qflex

#endif
