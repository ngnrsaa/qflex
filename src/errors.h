#ifndef ERROR__H
#define ERROR__H

#include <stdexcept>
#include <sstream>
#include <string>

template<typename T, typename... Q> std::string concat(const T& x, const Q&...y) {
  std::stringstream ss;
  ss << x;
  if constexpr (sizeof...(y)) return ss.str() + concat(y...);
  else return ss.str();
}

#define ERROR_MSG(...) concat("ERROR (" __FILE__ ": ", __FUNCTION__, ":", __LINE__, ") --> ", __VA_ARGS__)

#endif
