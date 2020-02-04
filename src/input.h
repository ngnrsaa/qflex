#ifndef INPUT__H
#define INPUT__H

#include "circuit.h"
#include "global.h"
#include "grid.h"
#include "ordering.h"

namespace qflex {

struct QflexInput {
  QflexOrdering ordering;
  QflexCircuit circuit;
  QflexGrid grid;
  std::vector<std::string> initial_states;
  std::vector<std::string> final_states;
};

}  // namespace qflex

#endif
