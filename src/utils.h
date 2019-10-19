#ifndef UTILS__H
#define UTILS__H

#include <iostream>
#include <regex>
#include <unordered_map>
#include <vector>

namespace std {

template <typename T, typename U>
struct hash<std::pair<T, U>> {
  std::size_t operator()(const std::pair<T, U>& p) const {
    return std::hash<T>()(p.first) ^ (std::hash<U>()(p.second) << 1);
  }
};

}  // namespace std

namespace qflex {

std::size_t compute_depth(std::istream&& istream);

}

#endif
