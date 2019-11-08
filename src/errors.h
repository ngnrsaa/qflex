#ifndef ERROR__H
#define ERROR__H

#include <sstream>
#include <stdexcept>
#include <string>

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

#define ERROR_MSG(...)                                                   \
  concat("ERROR (" __FILE__ ": ", __FUNCTION__, ":", __LINE__, ") --> ", \
         __VA_ARGS__)

#endif
