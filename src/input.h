#ifndef INPUT__H
#define INPUT__H

#include "global.h"
#include "circuit.h"
#include "grid.h"
#include "ordering.h"

namespace qflex {

struct QflexInput {
  QflexOrdering ordering;
  QflexCircuit circuit;
  QflexGrid grid;
  std::string initial_state;
  std::string final_state;
};

}  // namespace qflex

#endif
