#ifndef INPUT__H
#define INPUT__H

#include "ordering.h"
#include "circuit.h"
#include "grid.h"

namespace qflex {

struct QflexInput {
  QflexOrdering ordering;
  QflexCircuit circuit;
  QflexGrid grid;
  std::string initial_state;
  std::string final_state;
  bool enable_timing_logs = false;
};

}  // namespace qflex

#endif
