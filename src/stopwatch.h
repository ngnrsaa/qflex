#ifndef STOPWATCH__H
#define STOPWATCH__H

#include <chrono>

#include "errors.h"

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
   * Check if stopwatch is running
   * @return true if stopwatch is running, false otherwise.
   */
  bool is_running() const;

  /**
   * Check if stopwatch has started
   * @return true id stopwatch has already started, false otherwise.
   */
  bool has_started() const;

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
   * _split is updated everytime this function is called.
   * @return a difference in time in number (std::size_t) of units
   */
  template <typename unit>
  std::size_t split() {
    if (!has_started()) throw ERROR_MSG("Stopwatch never started.");

    if (!is_running()) throw ERROR_MSG("Stopwatch is not running.");

    auto _now = clock::now();
    auto _delta_time = std::chrono::duration_cast<unit>(_now - _split);
    _split = _now;
    return _delta_time.count();
  }

  /**
   * Compute the difference in time between _start and _latest. If _latest <=
   * _start, then it compute the difference in time between _start and
   * clock::now()
   * @return a difference in time in number (std::size_t) of units
   */
  template <typename unit>
  std::size_t time_passed() const {
    if (!has_started()) throw ERROR_MSG("Stopwatch never started.");

    if (is_running()) throw ERROR_MSG("Stopwatch is still running.");

    return std::chrono::duration_cast<unit>(
               (_latest > _start ? _latest : clock::now()) - _start)
        .count();
  }

 private:
  bool _started{false};
  bool _running{false};
  std::chrono::time_point<clock> _start, _latest, _split;
};

}  // namespace qflex::utils

#endif
