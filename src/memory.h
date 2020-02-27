#ifndef MEMORY__H
#define MEMORY__H

#include <signal.h>
#include <unistd.h>

#include <array>
#include <ctime>
#include <iostream>

#include "utils.h"

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
// Not tested
//#include <windows.h>
//#include <psapi.h>
#elif __APPLE__
#include <mach/mach.h>
#include <sys/resource.h>
#elif __unix__
#include <sys/resource.h>
#include <sys/time.h>

#include <fstream>
#endif

namespace qflex::memory {

/**
 * Get memory usage.
 * @return a std::array<std::string, 2>, with a[0] = "used resident memory" and
 * a[1] = "peak resident memory".
 */
inline std::array<std::string, 2> get_memory_usage() noexcept {
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)

  // Not tested
  // PROCESS_MEMORY_COUNTERS memCounter;
  // if (GetProcessMemoryInfo(GetCurrentProcess(), &memCounter, sizeof
  // memCounter))
  //  return {utils::readable_memory_string(info.WorkingSetSize),
  //  utils::readable_memory_string(info.PeakWorkingSetSize)};

#elif __APPLE__

  struct mach_task_basic_info info;
  mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info,
                &count) == KERN_SUCCESS)
    return {utils::readable_memory_string(info.resident_size),
            utils::readable_memory_string(info.resident_size_max)};

#elif __unix__

  if (auto in = std::ifstream("/proc/self/status"); in.good()) {
    // Get the value of a string in the format "label: value kB"
    auto get_value = [](std::string line) {
      // Strip everything which is not a number at the beginning of line
      line = line.substr(line.find_first_of("0123456789"), std::size(line));

      // Strip everything after the first space and get value
      return utils::readable_memory_string(
          std::stod(line.substr(0, line.find(' '))) * 1024);
    };

    // Read file and get values
    std::string line, m_used, m_peak;
    while (std::getline(in, line)) {
      if (line.find("VmHWM", 0) != std::string::npos)
        m_peak = get_value(line);
      else if (line.find("VmRSS", 0) != std::string::npos)
        m_used = get_value(line);
    }

    return {m_used, m_peak};
  }

#endif

  return {std::string("n/a"), std::string("n/a")};
}

/**
 * Print memory usage.
 * @param unused_signal_number is included to match the signature required by
 * signal.
 */
inline void print_memory_usage(int unused_signal_number = 0) noexcept {
  const auto [vm_used, vm_peak] = get_memory_usage();

  std::cerr << WARN_MSG("Memory usage: ", vm_used, " (Peak: ", vm_peak, ")")
            << std::endl;

  // Reset alarm
  if (global::track_memory_seconds > 0) alarm(global::track_memory_seconds);
}

}  // namespace qflex::memory

#endif
