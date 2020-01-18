#ifndef UTILS__H
#define UTILS__H

#include <cmath>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>

#include "errors.h"

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

// Convert more readable strings into 'memory' bytes.
std::size_t from_readable_memory_string(std::string memory);

// Convert 'memory' bytes into a more readable string.
std::string readable_memory_string(double memory);

}  // namespace qflex::utils

#endif
