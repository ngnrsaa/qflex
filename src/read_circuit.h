/**
 * @file read_circuit.h
 * Helper functions to read quantum circuits from a file.
 *
 * @author Benjamin Villalonga (main contributor), Bron Nelson, Sergio Boixo and
 * Salvatore Mandra
 * @contributors: The qFlex Developers (see CONTRIBUTORS.md)
 * @date Created: September 2018
 *
 * @copyright: Copyright © 2019, United States Government, as represented
 * by the Administrator of the National Aeronautics and Space Administration.
 * All rights reserved.
 * @licence: Apache License, Version 2.0
 */

#ifndef READ_CIRCUIT_
#define READ_CIRCUIT_

#include <algorithm>
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

#include "circuit.h"
#include "contraction_utils.h"
#include "global.h"
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

const std::size_t SUPER_CYCLE_DEPTH = 8;
const std::size_t DIM = 2;

/**
 * Return gate as an array.
 * @param name of the gate.
 * @param parameters for the gate.
 * @return an s_type array corresponding to the gate.
 */
std::vector<s_type> gate_array(const std::string& gate_name,
                               const std::vector<double>& params);

/**
 * Read circuit from stream and fill in a 2D grid of vectors of tensors.
 * @param qflex::QflexCircuit containing circuit information.
 * @param I std::size_t with the first spatial dimension of the grid of qubits.
 * @param J std::size_t with the second spatial dimension of the grid of qubits.
 * @param initial_conf string with 0s and 1s with the input configuration of
 * the circuit.
 * @param final_conf string with 0s and 1s with the output configuration on B.
 * @param final_qubit_region vector<vector<std::size_t>> with the coords. of the
 * qubits in qubits with terminal cuts.
 * @param off vector<vector<std::size_t>> with the coords. of the qubits turned
 * off.
 * @param grid_of_tensors referenced to a vector<vector<vector<Tensor>>> with
 * tensors (gates) at each position of the grid.
 * @param scratch pointer to s_type array with scratch space for all operations
 * performed in this function.
 */
void circuit_data_to_tensor_network(
    const QflexCircuit& circuit, std::size_t I, std::size_t J,
    const std::string& initial_conf,
    const std::optional<std::vector<std::vector<std::size_t>>>& off,
    std::vector<std::vector<std::vector<Tensor>>>& grid_of_tensors,
    s_type* scratch);

/**
 * Contracts a 2D grid of vectors of tensors onto a 2D grid of tensors,
 * contracting in the time (third) direction (i.e. flattening the tensor
 * network), and renaming the indices accordingly.
 * @param grid_of_tensors reference to a vector<vector<vector<Tensor>>> with the
 * 2D grid of vectors of tensors. The typical names for the indices in a grid
 * are assumed.
 * @param grid_of_tensors_2D reference to a vector<vector<Tensor>> where the
 * 2D grid of tensors will be stored. The typical names for the indices will
 * be used.
 * @param final_qubit_region optional<vector<vector<std::size_t>>> with the
 * coords. of the qubits in qubits with terminal cuts.
 * @param off optional<vector<vector<std::size_t>>> with the coords. of the
 * qubits turned off.
 * @param ordering std::list<ContractionOperation> providing the steps required
 * to contract the tensor grid.
 * @param scratch pointer to s_type array with enough space for all scratch
 * work.
 */
std::vector<std::vector<Tensor>> flatten_grid_of_tensors(
    std::vector<std::vector<std::vector<Tensor>>>& grid_of_tensors,
    const std::optional<std::vector<std::vector<std::size_t>>>& off,
    s_type* scratch);

void reorder_grid_of_tensors(
    std::vector<std::vector<Tensor>>* grid_of_tensors_2D_ptr,
    const std::optional<std::vector<std::vector<std::size_t>>>&
        final_qubit_region,
    const std::optional<std::vector<std::vector<std::size_t>>>& off,
    const std::list<ContractionOperation>& ordering, s_type* scratch);

/**
 * Read circuit from file and fill vector of tensors (of gates), vector with
 * names of the input indices of the tensors and vector with the output indices
 * of the tensors.
 * @param filename string with the name of the circuit file.
 * @param I std::size_t with the number of qubits.
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
    std::string filename, std::size_t I, std::vector<Tensor>& gates,
    std::vector<std::vector<std::string>>& inputs,
    std::vector<std::vector<std::string>>& outputs, s_type* scratch);

}  // namespace qflex

#endif  // READ_CIRCUIT_
