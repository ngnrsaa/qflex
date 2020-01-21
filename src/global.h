#ifndef GLOBAL__H
#define GLOBAL__H

namespace qflex::global {

inline int verbose = 0;

// Default limit of one gigabyte.
inline std::size_t memory_limit = 1L << 30;

// Interval to track memory
inline std::size_t track_memory;

}  // namespace qflex::global

#endif
