#include "ordering.h"

namespace qflex {

void QflexOrdering::load(std::istream& istream) { this->load(std::move(istream)); }
void QflexOrdering::load(std::istream&& istream) {}

void QflexOrdering::load(const std::string& filename) {
  if (auto in = std::ifstream(filename); in.good())
    this->load(in);
  else
    throw std::string("Cannot open ordering file: ") + filename;
}

}
