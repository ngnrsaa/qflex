#ifndef GRID__H
#define GRID__H

#include <algorithm>
#include <fstream>
#include <vector>

#include "global.h"

namespace qflex {

struct QflexGrid {
  std::size_t I{0}, J{0};
  std::vector<std::vector<std::size_t>> qubits_off;
  void clear();
  void load(std::istream& istream);
  void load(std::istream&& istream);
  void load(const std::string& filename);
};

}  // namespace qflex

#endif
