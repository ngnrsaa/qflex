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

template <typename _clock = std::chrono::steady_clock>
struct Stopwatch {
  using clock = _clock;

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
  void start() {
    if (_started) throw ERROR_MSG("Stopwatch already started.");

    if (_running) throw ERROR_MSG("Stopwatch already running.");

    _start = _latest = _split = clock::now();
    _started = _running = true;
  }

  /**
   * Stop Stopwatch
   */
  void stop() {
    if (!_started) throw ERROR_MSG("Stopwatch never started.");

    if (!_running) throw ERROR_MSG("Stopwatch is not running.");

    _latest = clock::now();
    _running = false;
  }

  /**
   * Reset Stopwatch
   */
  void reset() { _started = _running = false; }

  /**
   * Compute the difference in time between clock::now() and the last _split.
   * @return a difference in time in number (std::size_t) of units
   */
  template <typename unit>
  std::size_t split() {
    if (!_started) throw ERROR_MSG("Stopwatch never started.");

    if (!_running) throw ERROR_MSG("Stopwatch is not running.");

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
    if (!_started) throw ERROR_MSG("Stopwatch never started.");

    if (_running) throw ERROR_MSG("Stopwatch is still running.");

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
