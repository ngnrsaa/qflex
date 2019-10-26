#ifndef ORDERING__H
#define ORDERING__H

#include <fstream>
#include <vector>
#include <string>
#include <regex>

namespace qflex {

struct QflexOrdering {
  std::vector<std::string> instructions;
  void load(std::istream& istream);
  void load(std::istream&& istream);
  void load(const std::string& filename);
};

}

#endif
