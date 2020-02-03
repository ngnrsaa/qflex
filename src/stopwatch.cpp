#include "stopwatch.h"

#include "errors.h"

namespace qflex::utils {

void Stopwatch::start() {
  if (_started) throw ERROR_MSG("Stopwatch already started.");
  if (_running) throw ERROR_MSG("Stopwatch already running.");
  _start = _latest = _split = clock::now();
  _started = _running = true;
}

void Stopwatch::stop() {
  if (!_started) throw ERROR_MSG("Stopwatch never started.");
  if (!_running) throw ERROR_MSG("Stopwatch is not running.");
  _latest = clock::now();
  _running = false;
}

void Stopwatch::reset() { _started = _running = false; }

template <typename unit>
std::size_t Stopwatch::split() {
  if (!_started) throw ERROR_MSG("Stopwatch never started.");

  if (!_running) throw ERROR_MSG("Stopwatch is not running.");

  auto _now = clock::now();
  auto _delta_time = std::chrono::duration_cast<unit>(_now - _split);
  _split = _now;
  return _delta_time.count();
}

template <typename unit>
std::size_t Stopwatch::time_passed() const {
  if (!_started) throw ERROR_MSG("Stopwatch never started.");

  if (_running) throw ERROR_MSG("Stopwatch is still running.");

  return std::chrono::duration_cast<unit>(
             (_latest > _start ? _latest : clock::now()) - _start)
      .count();
}

}  // namespace qflex::utils
