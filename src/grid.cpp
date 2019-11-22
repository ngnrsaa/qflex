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
  std::size_t line_counter = 0;
  while (std::getline(istream, line))
    if (std::size(line) and line[0] != '#') {
      line_counter++;

      // String unnecessary chars
      line.erase(
          std::remove_if(std::begin(line), std::end(line),
                         [](auto&& x) { return not(x == '0' || x == '1'); }),
          std::end(line));

      // Continue if line is empty
      if (std::empty(line)) continue;

      // Get number of columns
      if (J != 0 and J != std::size(line))
        throw ERROR_MSG("Grid size at line ", line_counter,  " is inconsistent with a width of ", std::size(line), " instead of ", J, ".");
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
    throw ERROR_MSG("Cannot open grid file: ", filename, ".");
}

}  // namespace qflex
