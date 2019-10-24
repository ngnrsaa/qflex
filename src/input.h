#ifndef INPUT__H
#define INPUT__H

#include "circuit.h"

namespace qflex {

struct QflexGrid {
  int I{0}, J{0};
  std::vector<std::vector<int>> qubits_off;
  void load(std::istream& istream);
  void load(const std::string& filename);
};

struct QflexInput {
  std::istream* ordering_data;
  QflexCircuit circuit;
  QflexGrid grid;
  std::string initial_state;
  std::string final_state;
  bool enable_timing_logs = false;
};

}

#endif
