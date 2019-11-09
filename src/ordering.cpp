#include "ordering.h"

namespace qflex {

void QflexOrdering::clear() {
  this->instructions.clear();
}

void QflexOrdering::load(std::istream& istream) {
  this->load(std::move(istream));
}
void QflexOrdering::load(std::istream&& istream) {

  // Strip everthing from a line that doesn't follow the right format
  auto strip_line = [](std::string line) {
    // Remove everything after '#'
    line = std::regex_replace(line, std::regex("#.*"), "");

    // Remove any special character
    line = std::regex_replace(line, std::regex("[^)(\\s\\ta-zA-Z0-9_.,-]"), "");

    // Convert tabs to spaces
    line = std::regex_replace(line, std::regex("[\\t]"), " ");

    // Remove multiple spaces
    line = std::regex_replace(line, std::regex("[\\s]{2,}"), " ");

    // Remove last space
    line = std::regex_replace(line, std::regex("\\s+$"), "");

    // Remove spaces between parentheses
    line = std::regex_replace(line, std::regex("\\s+(?=[^()]*\\))"), "");

    // After stripping, line should follow the format
    // cut (1,0) 39 45
    // expand A 10
    // ...

    return line;
  };

  // Clear this ordering
  this->clear();

  std::string line;
  while (std::getline(istream, line))
    if (std::size(line = strip_line(line))) this->instructions.push_back(line);
}
void QflexOrdering::load(const std::string& filename) {
  if (auto in = std::ifstream(filename); in.good())
    this->load(in);
  else
    throw std::string("Cannot open ordering file: ") + filename;
}

}  // namespace qflex
