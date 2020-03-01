#ifndef GLOBAL__H
#define GLOBAL__H

#include <cstddef>

namespace qflex::global {

// Default verbose level
inline int verbose = 0;

// Max allowed memory (default: 1GB)
inline std::size_t memory_limit = 1L << 30;

// Interval to track memory (default: 0)
inline std::size_t track_memory_seconds = 0;

}  // namespace qflex::global

#endif
