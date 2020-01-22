#ifndef MEMORY__H
#define MEMORY__H

#include <signal.h>
#include <unistd.h>

#include <array>
#include <fstream>
#include <iostream>

#include "utils.h"

namespace qflex::memory {

inline std::array<std::string, 2> get_memory_usage() noexcept {
  if (auto in = std::ifstream("/proc/self/status"); in.good()) {
    // Skip to the desired line
    std::string line;
    for (std::size_t i = 0; i < 17; ++i) std::getline(in, line);

    // Strip everything which is not a number at the beginning of line
    line = line.substr(line.find_first_of("0123456789"), std::size(line));

    // Strip everything after the first space and get memory peak
    const std::string m_peak = utils::readable_memory_string(
        std::stod(line.substr(0, line.find(' '))) * 1024);

    // Get next line
    std::getline(in, line);

    // Strip everything which is not a number at the beginning of line
    line = line.substr(line.find_first_of("0123456789"), std::size(line));

    // Strip everything after the first space and get memory usage
    const std::string m_used = utils::readable_memory_string(
        std::stod(line.substr(0, line.find(' '))) * 1024);

    return {m_used, m_peak};

  } else
    return {std::string("n/a"), std::string("n/a")};
}

inline void print_memory_usage(int _ = 0) noexcept {
  const auto [vm_used, vm_peak] = get_memory_usage();
  std::cerr << "Memory usage: " << vm_used << " (Peak: " << vm_peak << ")"
            << std::endl;
}

}  // namespace qflex::memory

#endif
