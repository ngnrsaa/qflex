/**
 * @file evaluate_circuit.h
 *
 * @author Benjamin Villalonga (main contributor), Bron Nelson, Sergio Boixo and
 * Salvatore Mandra
 * @contributors: The qFlex Developers (see CONTRIBUTORS.md)
 * @date Created: August 2018
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
#include "global.h"
#include "input.h"
#include "read_circuit.h"
#include "tensor.h"

namespace qflex {

/**
 * Reads in grid layout from a file, which should be formatted as an I x J grid
 * of zeroes (for "off" qubits) and ones (for "on" qubits).\
 * @param grid_data std::istream* containing grid layout stored as a string.
 * @param I std::size_t with the first spatial dimension of the grid of qubits.
 * @param J std::size_t with the second spatial dimension of the grid of qubits.
 * @return a list of coordinates for "off" qubits in the I x J grid provided.
 */
std::vector<std::vector<std::size_t>> read_grid_layout_from_stream(
    std::istream* grid_data, std::size_t I, std::size_t J);

struct QflexFinalQubits {
  std::vector<std::vector<std::size_t>> qubits;
  std::vector<std::size_t> output_pos_map;
  std::vector<std::vector<std::size_t>> output_values_map;
};

QflexFinalQubits get_final_qubits(
    const QflexGrid& grid, const std::list<ContractionOperation>& ordering);

/**
 * Determines the final qubit positions and output states for a given ordering.
 * @param input QflexInput generated from the command line.
 * @param ordering std::list<ContractionOperation> to parse output states from
 * @param final_qubits vector of coordinates for qubits with terminal cuts, to
 * be populated by this method.
 * @param output_states vector of output states for the given contraction
 * ordering, to be populated by this method.
 */
std::vector<std::string> get_output_states(
    const std::string& base_state, const QflexFinalQubits& final_qubits);

void apply_terminal_cuts(
    const QflexGrid& grid, const QflexFinalQubits& final_qubits,
    std::vector<std::vector<std::vector<Tensor>>>* tensor_grid_3D_ptr);

void apply_delta_output(
    const QflexGrid& grid, const std::string& final_state,
    const QflexFinalQubits& final_qubits,
    const std::vector<std::vector<std::vector<Tensor>>>& tensor_grid_3D,
    std::vector<std::vector<Tensor>>* tensor_grid_prt,
    std::vector<s_type>* scratch_2D_ptr);

/**
 * Evaluates a circuit and returns the final amplitudes of each state resulting
 * from the provided contraction ordering.
 *
 * Usage:
 * $ ./qflex.x I J K fidelity circuit_filename ordering_filename \
 *       grid_filename [initial_state] [final_state]
 *
 * @param input args required to specify a circuit for evaluation.
 * @return vector of <state bitstring, amplitude> pairs for each output state.
 * States for qubits with terminal cuts are listed at the end of the state
 * bitstring, in the order of their terminal cuts.
 */
std::vector<std::tuple<std::string, std::string, std::complex<double>>>
EvaluateCircuit(const QflexInput& input);

}  // namespace qflex

#endif  // EVALUATE_CIRCUIT_
