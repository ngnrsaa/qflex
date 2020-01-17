#ifndef MEMORY__H
#define MEMORY__H

#include <signal.h>
#include <unistd.h>

#include <fstream>
#include <iostream>

namespace qflex::memory {

inline std::string get_peak_memory_usage() noexcept {
  if (auto in = std::ifstream("/proc/self/status"); in.good()) {
    // Skip to the desired line
    std::string line;
    for (std::size_t i = 0; i < 17; ++i) std::getline(in, line);

    // Strip everything which is not a number at the beginning of line
    return line.substr(line.find_first_of("0123456789"), std::size(line));

  } else
    return std::string("n/a");
}

inline void print_peak_memory_usage(int _ = 0) noexcept {
  static std::string old_vm_size = get_peak_memory_usage();
  if (const std::string vm_size = get_peak_memory_usage();
      old_vm_size != vm_size) {
    std::cerr << "Memory usage: " << vm_size << std::endl;
    old_vm_size = vm_size;
  }
}

}  // namespace qflex::memory

#endif
