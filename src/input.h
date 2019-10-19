#ifndef INPUT__H
#define INPUT__H

#include <fstream>
#include <string>

#include "circuit.h"
#include "grid.h"

namespace qflex {

struct QflexInput {
  int K;
  std::istream* circuit_data;
  std::istream* ordering_data;
  QflexCircuit circuit;
  QflexGrid grid;
  std::string initial_state;
  std::string final_state;
  bool enable_timing_logs = false;
};

}

#endif
