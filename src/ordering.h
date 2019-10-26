#ifndef ORDERING__H
#define ORDERING__H

#include <fstream>

namespace qflex {

struct QflexOrdering {
  void load(std::istream& istream);
  void load(std::istream&& istream);
  void load(const std::string& filename);
};

}

#endif
