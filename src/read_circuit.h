/**
 * @file read_circuit.h
 * Helper functions to read quantum circuits from a file.
 * @see https://github.com/benjaminvillalonga/optimized_parallel_QC_with_TN
 *
 * @author Benjamin Villalonga (main contributor), Bron Nelson, Sergio Boixo and
 * Salvatore Mandra
 * @contributors: The qFlex Developers (see CONTRIBUTORS.md)
 * @date Created: September 2018
 * @date Modified: October 2019
 *
 * @copyright: Copyright © 2019, United States Government, as represented
 * by the Administrator of the National Aeronautics and Space Administration.
 * All rights reserved.
 * @licence: Apache License, Version 2.0
 */

#ifndef READ_CIRCUIT_
#define READ_CIRCUIT_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <ctime>
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "contraction_utils.h"
#include "circuit.h"
#include "tensor.h"

namespace std {

template <typename T, typename U>
struct hash<std::pair<T, U>> {
  std::size_t operator()(const std::pair<T, U>& p) const {
    return std::hash<T>()(p.first) ^ (std::hash<U>()(p.second) << 1);
  }
};

}  // namespace std

namespace qflex {

// TODO: Use math library calls for these. Also use constexpr.
const double _SQRT_2 = 1.41421356237309504880168872;
const double _INV_SQRT_2 = 1. / _SQRT_2;
const double _PI = 3.14159265358979323846264338;
const int SUPER_CYCLE_DEPTH = 8;
const int DIM = 2;

/**
 * Read circuit from stream and fill in a 2D grid of vectors of tensors.
 * @param qflex::QflexCircuit containing circuit information.
 * @param I int with the first spatial dimension of the grid of qubits.
 * @param J int with the second spatial dimension of the grid of qubits.
 * @param initial_conf string with 0s and 1s with the input configuration of
 * the circuit.
 * @param final_conf string with 0s and 1s with the output configuration on B.
 * @param final_qubit_region vector<vector<int>> with the coords. of the
 * qubits in qubits with terminal cuts.
 * @param off vector<vector<int>> with the coords. of the qubits turned off.
 * @param grid_of_tensors referenced to a vector<vector<vector<Tensor>>> with
 * tensors (gates) at each position of the grid.
 * @param scratch pointer to s_type array with scratch space for all operations
 * performed in this function.
 */
void circuit_data_to_tensor_network(
    const QflexCircuit &circuit, int I, int J,
    const std::string initial_conf, const std::string final_conf,
    const std::optional<std::vector<std::vector<int>>>& final_qubit_region,
    const std::optional<std::vector<std::vector<int>>>& off,
    std::vector<std::vector<std::vector<Tensor>>>& grid_of_tensors,
    s_type* scratch);

/**
 * Contracts a 2D grid of vectors of tensors onto a 2D grid of tensors,
 * contracting in the time (third) direction (i.e. flattening the tensor
 * network), and renaming the indices accordingly.
 * @param grid_of_tensors reference to a vector<vector<vector<Tensor>>> with the 2D
 * grid of vectors of tensors. The typical names for the indices in a grid are
 * assumed.
 * @param grid_of_tensors_2D reference to a vector<vector<Tensor>> where the
 * 2D grid of tensors will be stored. The typical names for the indices will
 * be used.
 * @param final_qubit_region optional<vector<vector<int>>> with the coords.
 * of the qubits in qubits with terminal cuts.
 * @param off optional<vector<vector<int>>> with the coords. of the qubits
 * turned off.
 * @param ordering std::list<ContractionOperation> providing the steps required
 * to contract the tensor grid.
 * @param scratch pointer to s_type array with enough space for all scratch
 * work.
*/
void flatten_grid_of_tensors(
    std::vector<std::vector<std::vector<Tensor>>>& grid_of_tensors,
    std::vector<std::vector<Tensor>>& grid_of_tensors_2D,
    std::optional<std::vector<std::vector<int>>> final_qubit_region,
    std::optional<std::vector<std::vector<int>>> off,
    const std::list<ContractionOperation>& ordering, s_type* scratch);

/**
 * Read circuit from file and fill vector of tensors (of gates), vector with
 * names of the input indices of the tensors and vector with the output indices
 * of the tensors.
 * @param filename string with the name of the circuit file.
 * @param I int with the number of qubits.
 * @param gates reference to a vector<Tensor> to be filled with the gates.
 * @param inputs reference to a vector<vector<string>> to be filled with the
 * input indexes of the gates.
 * @param outputs reference to a vector<vector<string>> to be filled with the
 * output indexes of the gates.
 * @param scratch pointer to s_type array with scratch space for operations
 * performed in this function.
 *
 */
void read_wave_function_evolution(
    std::string filename, int I, std::vector<Tensor>& gates,
    std::vector<std::vector<std::string>>& inputs,
    std::vector<std::vector<std::string>>& outputs, s_type* scratch);

}  // namespace qflex

#endif  // READ_CIRCUIT_
