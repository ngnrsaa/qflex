#include "grid.h"

namespace qflex {

void QflexGrid::clear() {
  this->I = 0;
  this->J = 0;
  qubits_off.clear();
}
void QflexGrid::load(std::istream& istream) { this->load(std::move(istream)); }
void QflexGrid::load(std::istream&& istream) {
  // Clear this grid
  this->clear();

  std::string line;
  while (std::getline(istream, line))
    if (std::size(line) and line[0] != '#') {
      // String unnecessary chars
      line.erase(
          std::remove_if(std::begin(line), std::end(line),
                         [](auto&& x) { return not(x == '0' || x == '1'); }),
          std::end(line));

      // Continue if line is empty
      if (std::empty(line)) continue;

      // Get number of columns
      if (J != 0 and J != std::size(line))
        throw std::string("Grid size is inconsistent");
      else
        J = std::size(line);

      // Get off qubits
      for (int q = 0; q < J; ++q)
        if (line[q] == '0') qubits_off.push_back({I, q});

      // Update number of rows
      ++I;
    }
}

void QflexGrid::load(const std::string& filename) {
  if (auto in = std::ifstream(filename); in.good())
    this->load(in);
  else
    throw std::string("Cannot open grid file: ") + filename;
}

}  // namespace qflex
