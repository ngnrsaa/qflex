/**
 * @file read_circuit.h
 * Helper functions to read quantum circuits from a file.
 * @see https://github.com/benjaminvillalonga/optimized_parallel_QC_with_TN
 *
 * @author Benjamin Villalonga
 * @date Created: September 2018
 * @date Modified: February 2019
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
#include <vector>

#include "contraction_utils.h"
#include "mkl_tensor.h"

namespace qflex {

// TODO: Use math library calls for these. Also use constexpr.
const double _SQRT_2 = 1.41421356237309504880168872;
const double _INV_SQRT_2 = 1. / _SQRT_2;
const double _PI = 3.14159265358979323846264338;
const int SUPER_CYCLE_DEPTH = 8;
const int DIM = 2;

// Methods in this namespace should only be used in read_circuit.cpp and
// read_circuit_test.cpp.
namespace internal {

/**
 * Helper method for google_circuit_file_to_grid_of_tensors. External users
 * should call that method instead.
 * @param circuit_data std::istream containing circuit as a string.
 * @param I int with the first spatial dimension of the grid of qubits.
 * @param J int with the second spatial dimension of the grid of qubits.
 * @param K int with depth of the grid of tensors (depth_of_circuit-2)/8.
 * @param initial_conf string with 0s and 1s with the input configuration of
 * the circuit.
 * @param final_conf_B string with 0s and 1s with the output configuration on B.
 * @param A vector<vector<int>> with the coords. of the qubits in A.
 * @param off vector<vector<int>> with the coords. of the qubits turned off.
 * @param grid_of_tensors referenced to a vector<vector<vector<MKLTensor>>> with
 * tensors at each position of the grid.
 * @param scratch pointer to s_type array with scratch space for all operations
 * performed in this function.
 */
void circuit_data_to_grid_of_tensors(
    std::istream* circuit_data, int I, int J, int K,
    const std::string initial_conf, const std::string final_conf_B,
    const std::optional<std::vector<std::vector<int>>>& A,
    const std::optional<std::vector<std::vector<int>>>& off,
    std::vector<std::vector<std::vector<MKLTensor>>>& grid_of_tensors,
    s_type* scratch);

}  // namespace internal

/**
 * Read circuit from file and fill in a 3D grid of tensors with the indices
 * labelled as "(i1,j1,k1),(i2,j2,k2)", where i1<=i2, j1<=j2, and k1<=k2.
 * I*J must be equal to the number of qubits; K must be equal to
 * (depth_of_circuit-2)/8; initial_conf and final_conf must have the length
 * equal to the number of qubits.
 * @param filename string with the name of the circuit file.
 * @param I int with the first spatial dimension of the grid of qubits.
 * @param J int with the second spatial dimension of the grid of qubits.
 * @param K int with depth of the grid of tensors (depth_of_circuit-2)/8.
 * @param initial_conf string with 0s and 1s with the input configuration of
 * the circuit.
 * @param final_conf_B string with 0s and 1s with the output configuration on B.
 * @param A vector<vector<int>> with the coords. of the qubits in A.
 * @param off vector<vector<int>> with the coords. of the qubits turned off.
 * @param grid_of_tensors referenced to a vector<vector<vector<MKLTensor>>>
 * with tensors at * each position of the grid.
 * @param scratch pointer to s_type array with scratch space for all operations
 * performed in this function.
 */
void google_circuit_file_to_grid_of_tensors(
    std::string filename, int I, int J, int K, const std::string initial_conf,
    const std::string final_conf_B,
    const std::optional<std::vector<std::vector<int>>>& A,
    const std::optional<std::vector<std::vector<int>>>& off,
    std::vector<std::vector<std::vector<MKLTensor>>>& grid_of_tensors,
    s_type* scratch);

/**
 * Contracts a 3D grid of tensors onto a 2D grid of tensors, contracting
 * in the time (third) direction, and renaming the indices accordingly.
 * @param grid_of_tensors_3D reference to a
 * vector<vector<vector<MKLTensor>>> with the 3D grid of tensors. It must be a
 * grid dimensionswise. The typical names for the indices in a grid is assumed.
 * @param grid_of_tensors_2D reference to a vector<vector<MKLTensor>> where the
 * 2D grid of tensors will be stored. The typical names for the indices will
 * be used.
 * @param A optional<vector<vector<int>>> with the coords. of the qubits in A.
 * @param off optional<vector<vector<int>>> with the coords. of the qubits
 * turned off.
 * @param ordering ContractionOrdering providing the steps required to contract
 * the tensor grid.
 * @param scratch pointer to s_type array with enough space for all scratch
 * work.
 */
void grid_of_tensors_3D_to_2D(
    std::vector<std::vector<std::vector<MKLTensor>>>& grid_of_tensors_3D,
    std::vector<std::vector<MKLTensor>>& grid_of_tensors_2D,
    std::optional<std::vector<std::vector<int>>> A,
    std::optional<std::vector<std::vector<int>>> off,
    const ContractionOrdering& ordering, s_type* scratch);

/**
 * Read circuit from file and fill vector of tensors (of gates), vector with
 * names of the input indices of the tensors and vector with the output indices
 * of the tensors.
 * @param filename string with the name of the circuit file.
 * @param I int with the number of qubits.
 * @param gates reference to a vector<MKLTensor> to be filled with the gates.
 * @param inputs reference to a vector<vector<string>> to be filled with the
 * input indexes of the gates.
 * @param outputs reference to a vector<vector<string>> to be filled with the
 * output indexes of the gates.
 * @param scratch pointer to s_type array with scratch space for operations
 * performed in this function.
 *
 */
void read_wave_function_evolution(
    std::string filename, int I, std::vector<MKLTensor>& gates,
    std::vector<std::vector<std::string>>& inputs,
    std::vector<std::vector<std::string>>& outputs, s_type* scratch);

}  // namespace qflex

#endif  // READ_CIRCUIT_
