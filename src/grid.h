#ifndef GRID__H
#define GRID__H

#include <algorithm>
#include <fstream>
#include <vector>

namespace qflex {

struct QflexGrid {
  std::size_t I{0}, J{0};
  #warning FIX HERE
  std::vector<std::vector<int>> qubits_off;
  void load(std::istream& istream);
  void load(std::istream&& istream);
  void load(const std::string& filename);
};

}  // namespace qflex

#endif
