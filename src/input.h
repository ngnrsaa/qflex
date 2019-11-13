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
  std::string initial_state;
  std::string final_state;
};

}  // namespace qflex

#endif
