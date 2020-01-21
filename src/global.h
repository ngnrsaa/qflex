#ifndef GLOBAL__H
#define GLOBAL__H

namespace qflex::global {

// Default verbose level
inline int verbose = 0;

// Max allowed memory (default: 1GB)
inline std::size_t memory_limit = 1L << 30;

}  // namespace qflex::global

#endif
