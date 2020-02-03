#ifndef STOPWATCH__H
#define STOPWATCH__H

#include <chrono>

namespace qflex::utils {

using nanoseconds = std::chrono::nanoseconds;
using microseconds = std::chrono::microseconds;
using milliseconds = std::chrono::milliseconds;
using seconds = std::chrono::seconds;
using minutes = std::chrono::minutes;
using hours = std::chrono::hours;

struct Stopwatch {
  using clock = std::chrono::steady_clock;

  /**
   * Initialize Stopwatch
   */
  Stopwatch() {}

  Stopwatch(const Stopwatch &) = delete;
  Stopwatch(Stopwatch &&) = delete;

  Stopwatch &operator=(const Stopwatch &) = delete;
  Stopwatch &operator=(Stopwatch &&) = delete;

  /**
   * Check if the stopwatch is running
   * @return true if stopwatch is running, false otherwise.
   */
  bool running() const { return _running; }

  /**
   * Start Stopwatch
   */
  void start();

  /**
   * Stop Stopwatch
   */
  void stop();

  /**
   * Reset Stopwatch
   */
  void reset();

  /**
   * Compute the difference in time between clock::now() and the last _split.
   * @return a difference in time in number (std::size_t) of units
   */
  template <typename unit>
  std::size_t split();

  /**
   * Compute the difference in time between _start and _latest. If _latest <=
   * _start, then it compute the difference in time between _start and
   * clock::now()
   * @return a difference in time in number (std::size_t) of units
   */
  template <typename unit>
  std::size_t time_passed() const;

 private:
  bool _started{false};
  bool _running{false};
  std::chrono::time_point<clock> _start, _latest, _split;
};

}  // namespace qflex::utils

#endif
