#ifndef UTILS__CPP
#define UTILS__CPP

#include "utils.h"

namespace qflex::utils {

// Return date and time in the representation: Thu Aug 23 14:55:02 2001
std::string get_date_time() {
  std::time_t t = std::time(nullptr);
  if (char mbstr[100];
      std::strftime(mbstr, sizeof(mbstr), "%c", std::localtime(&t)))
    return mbstr;
  else
    return "[time-readout error]";
}

// Convert more readable strings into 'memory' bytes.
std::size_t from_readable_memory_string(std::string memory) {
  // Remove any special character
  memory =
      std::regex_replace(memory, std::regex("[^)(\\s\\ta-zA-Z0-9_.,-]"), "");

  // Convert tabs to spaces
  memory = std::regex_replace(memory, std::regex("[\\t]"), " ");

  // Remove all spaces
  memory = std::regex_replace(memory, std::regex(" "), "");

  // Convert all letters to uppercase
  for (auto &c : memory) c = std::toupper(c);

  std::size_t m =
      std::stoll(memory.substr(0, memory.find_first_not_of("0123456789")));

  {
    // Find the first non-number and get the suffix
    const std::string suffix = memory.substr(
        std::min(memory.find_first_not_of("0123456789"), std::size(memory)),
        std::size(memory));

    if (suffix == "GB")
      m *= std::pow(1024, 3);
    else if (suffix == "MB")
      m *= std::pow(1024, 2);
    else if (suffix == "KB")
      m *= std::pow(1024, 1);
    else if (suffix == "B")
      m *= std::pow(1024, 0);
    else if (suffix == "")
      m *= std::pow(1024, 0);
    else
      throw ERROR_MSG("Suffix ", suffix, " not recognized.");
  }

  return m;
}

// Convert 'memory' bytes into a more readable string.
std::string readable_memory_string(double memory) {
  std::string suffix[] = {" B", " kB", " MB", " GB"};
  std::size_t scale = 0;
  while (memory >= (1 << 10)) {
    ++scale;
    memory /= (1 << 10);
  }
  return concat(static_cast<std::size_t>(memory * 100) / 100., suffix[scale]);
}

}  // namespace qflex::utils

#endif
