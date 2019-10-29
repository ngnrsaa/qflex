#ifndef ORDERING__H
#define ORDERING__H

#include <fstream>
#include <regex>
#include <string>
#include <vector>

namespace qflex {

struct QflexOrdering {
  std::vector<std::string> instructions;
  void load(std::istream& istream);
  void load(std::istream&& istream);
  void load(const std::string& filename);
};

}  // namespace qflex

#endif
