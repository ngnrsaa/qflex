/**
 * @file evaluate_circuit.h
 * @see https://github.com/benjaminvillalonga/optimized_parallel_QC_with_TN
 *
 * @author Benjamin Villalonga (main contributor), Bron Nelson, Sergio Boixo and
 * Salvatore Mandra
 * @date Created: August 2018
 * @date Modified: August 2018
 *
 * @copyright: Copyright © 2019, United States Government, as represented
 * by the Administrator of the National Aeronautics and Space Administration.
 * All rights reserved.
 * @licence: Apache License, Version 2.0
 */

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

struct QflexGrid {
  int I{0}, J{0};
  std::vector<std::vector<int>> qubits_off;
  void load(std::istream& istream);
  void load(const std::string& filename);
};

struct QflexInput {
  int K;
  std::istream* circuit_data;
  std::istream* ordering_data;
  QflexGrid grid;
  std::string initial_state;
  std::string final_state_A;
  bool enable_timing_logs = false;
};

/**
 * Reads in grid layout from a file, which should be formatted as an I x J grid
 * of zeroes (for "off" qubits) and ones (for "on" qubits).\
 * @param grid_data std::istream* containing grid layout stored as a string.
 * @param I int with the first spatial dimension of the grid of qubits.
 * @param J int with the second spatial dimension of the grid of qubits.
 * @return a list of coordinates for "off" qubits in the I x J grid provided.
 */
std::vector<std::vector<int>> read_grid_layout_from_stream(
    std::istream* grid_data, int I, int J);

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
