#ifndef STOPWATCH__H
#define STOPWATCH__H

#include <chrono>

namespace qflex::utils {

using nanoseconds = std::chrono::nanoseconds;
using microseconds = std::chrono::microseconds;
using seconds = std::chrono::seconds;
using minutes = std::chrono::minutes;
using hours = std::chrono::hours;

template <typename _clock = std::chrono::steady_clock>
struct Stopwatch {
  using clock = _clock;

  /**
   * Initialize Stopwatch
   */
  Stopwatch() { start(); }

  Stopwatch(const Stopwatch &) = delete;
  Stopwatch(Stopwatch &&) = delete;

  Stopwatch &operator=(const Stopwatch &) = delete;
  Stopwatch &operator=(Stopwatch &&) = delete;

  /**
   * Start Stopwatch
   */
  void start() { _start = _latest = clock::now(); }

  /**
   * Stop Stopwatch
   */
  void stop() { _latest = clock::now(); }

  /**
   * Compute the difference in time between _start and _latest. If _latest <=
   * _start, then it compute the difference in time between _start and
   * clock::now()
   * @return a difference in time in number (std::size_t) of units
   */
  template <typename unit>
  std::size_t time_passed() const {
    return std::chrono::duration_cast<unit>(
               (_latest > _start ? _latest : clock::now()) - _start)
        .count();
  }

 private:
  std::chrono::time_point<clock> _start, _latest;
};

}  // namespace qflex::utils

#endif
