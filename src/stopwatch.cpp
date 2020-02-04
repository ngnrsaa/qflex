#include "stopwatch.h"

#include "errors.h"

namespace qflex::utils {

bool Stopwatch::is_running() const { return _running; }

bool Stopwatch::has_started() const { return _started; }

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

}  // namespace qflex::utils
