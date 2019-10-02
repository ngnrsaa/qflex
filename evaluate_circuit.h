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
#include "read_circuit.h"
#include "tensor.h"

namespace qflex {

struct QflexInput {
  int grid_height;
  int grid_width;
  int super_cycles;

  //deprecated?
  double fidelity;

  bool enable_timing_logs = false;

  std::istream* circuit_data;
  std::istream* ordering_data;
  std::istream* grid_data;
  std::string initial_state;
  std::string final_state_A;
};

/**
 * Reads in grid layout from a file, which should be formatted as an grid_height x grid_width grid
 * of zeroes (for "off" qubits) and ones (for "on" qubits).\
 * @param grid_data std::istream* containing grid layout stored as a string.
 * @param grid_height int with the first spatial dimension of the grid of qubits.
 * @param grid_width int with the second spatial dimension of the grid of qubits.
 * @return a list of coordinates for "off" qubits in the grid_height x grid_width grid provided.
 */
std::vector<std::vector<int>> read_grid_layout_from_stream(
    std::istream* grid_data, int grid_height, int grid_width);

/**
 * Determines the final qubit positions and output states for a given ordering.
 * @param ordering std::list<ContractionOperation> to parse output states from
 * @param final_qubits vector of coordinates for qubits with terminal cuts, to
 * be populated by this method.
 * @param output_states vector of output states for the given contraction
 * ordering, to be populated by this method.
 */
void get_output_states(const std::list<ContractionOperation>& ordering,
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
 * States for qubits with terminal cuts are listed at the end onf the state
 * bitstring, in the order of their terminal cuts.
 */
std::vector<std::pair<std::string, std::complex<double>>> EvaluateCircuit(
    QflexInput* input);

}  // namespace qflex

#endif  // EVALUATE_CIRCUIT_
