#ifndef ERROR__H
#define ERROR__H

#include "global.h"
#include "utils.h"

#define ERROR_MSG(...)                                                       \
  qflex::utils::concat("ERROR (" __FILE__ ": ", __FUNCTION__, ":", __LINE__, \
                       ") --> ", __VA_ARGS__)

template <typename... T>
std::string WARN_MSG(const T&... x) {
  std::time_t t = std::time(nullptr);
  if (char mbstr[100];
      std::strftime(mbstr, sizeof(mbstr), "%c", std::localtime(&t)))
    return qflex::utils::concat("[", mbstr, "] ", x...);
  else
    return "";
}

#endif
