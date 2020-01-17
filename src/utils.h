#ifndef UTILS__H
#define UTILS__H

#include <sstream>
#include <stdexcept>
#include <string>

namespace qflex::utils {

template <typename T>
std::string concat(const T& x) {
  std::stringstream ss;
  ss << x;
  return ss.str();
}

template <typename T, typename... Q>
std::string concat(const T& x, const Q&... y) {
  std::stringstream ss;
  ss << x;
  return ss.str() + concat(y...);
}

}  // namespace qflex::utils

#endif
