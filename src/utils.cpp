#ifndef UTILS__CPP
#define UTILS__CPP

#include "utils.h"

namespace qflex::utils {

// Convert 'memory' bytes into a more readable string.
std::string readable_memory_string(double memory) {
  std::string suffix[] = {" B", " kB", " MB", " GB"};
  std::size_t scale = 0;
  while (memory >= (1 << 10)) {
    ++scale;
    memory /= (1 << 10);
  }
  return concat(static_cast<std::size_t>(memory*100)/100., suffix[scale]);
}

}  // namespace qflex::utils

#endif
