#ifndef EVALUATE_CIRCUIT_
#define EVALUATE_CIRCUIT_

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "contraction_utils.h"
#include "mkl_tensor.h"
#include "read_circuit.h"

namespace qflex {

struct QflexInput {
  int I, J, K;
  double fidelity;
  std::string circuit_filename;
  std::string ordering_filename;
  std::string grid_filename;
  std::string initial_state;
  std::string final_state_A;
  bool enable_timing_logs = false;
};

// TODO(martinop): replace all "read from file" methods with stream-passing.
// Reads in grid layout from a file, which should be formatted as an I x J grid
// of zeroes (for "off" qubits) and ones (for "on" qubits).
std::vector<std::vector<int>> read_grid_layout_from_file(
    int I, int J, std::string grid_filename);

// Determines the final qubit positions and output states for a given ordering
void get_output_states(const ContractionOrdering& ordering,
                       std::vector<std::vector<int>>* final_qubits,
                       std::vector<std::string>* output_states);

/**
 * Evaluates a circuit and returns the final amplitudes of each state resulting
 * from the provided contraction ordering.
 *
 * Usage:
 * $ ./qflex.x I J K fidelity circuit_filename ordering_filename \
 *       grid_filename [initial_state] [final_state_A]
 *
 * @param input args required to specify a circuit for evaluation.
 * @return vector of <state bitstring, amplitude> pairs for each output state.
 */
std::vector<std::pair<std::string, std::complex<double>>> EvaluateCircuit(
    QflexInput* input);

}  // namespace qflex

#endif  // EVALUATE_CIRCUIT_
