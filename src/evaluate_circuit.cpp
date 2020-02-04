/**
 * @file evaluate_circuit.cpp
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

#include "evaluate_circuit.h"

#include "stopwatch.h"

namespace qflex {

// Gets the position in the output state vector of the qubit at tensor_pos.
std::size_t find_output_pos(const QflexInput* input,
                            std::vector<std::size_t> tensor_pos) {
  if (input == nullptr) {
    throw ERROR_MSG("Input must be non-null.");
  }
  std::size_t pos = (tensor_pos[0] * input->grid.J) + tensor_pos[1];
  for (const auto off_pos : input->grid.qubits_off) {
    if (off_pos[0] < tensor_pos[0]) {
      --pos;
    } else if (off_pos[0] == tensor_pos[0] && off_pos[1] < tensor_pos[1]) {
      --pos;
    }
  }
  return pos;
}

std::string get_output_states(
    const QflexInput* input, const std::list<ContractionOperation>& ordering,
    std::vector<std::vector<std::size_t>>* final_qubits,
    std::vector<std::string>* output_states) {
  if (input == nullptr) {
    throw ERROR_MSG("Input must be non-null.");
  }
  if (final_qubits == nullptr) {
    throw ERROR_MSG("Final qubits must be non-null.");
  }
  if (output_states == nullptr) {
    throw ERROR_MSG("Output states must be non-null");
  }
  std::vector<std::size_t> output_pos_map;
  std::vector<std::vector<std::size_t>> output_values_map;
  std::string base_state = input->final_state;
  // If the final state isn't provided, it should be all zeroes except for
  // qubits with terminal cuts (which should have 'x').
  bool final_state_unspecified = false;
  if (input->final_state.empty()) {
    final_state_unspecified = true;
    base_state = std::string(input->initial_state.length(), '0');
  }
  for (const auto& op : ordering) {
    // TODO(martinop): update to use the new operation.
    if (op.op_type != ContractionOperation::CUT) continue;
    // Any qubit with a terminal cut is in the final region.
    if (op.cut.tensors.size() != 1) continue;
    const std::size_t output_pos = find_output_pos(input, op.cut.tensors[0]);
    const auto tensor_pos = op.cut.tensors[0];
    if (final_state_unspecified) {
      base_state[output_pos] = 'x';
    }
    // TODO(martinop): reconsider requiring 'x' for cut indices.
    output_pos_map.push_back(output_pos);
    if (op.cut.values.empty()) {
      output_values_map.push_back({0, 1});
    } else {
      output_values_map.push_back(op.cut.values);
    }
    final_qubits->push_back(tensor_pos);
  }
  // Construct the full set of output state strings.
  std::vector<std::string> temp_output_states;
  output_states->push_back(base_state);
  for (std::size_t i = 0; i < final_qubits->size(); ++i) {
    const std::size_t pos = output_pos_map[i];
    for (const std::string& state : *output_states) {
      for (const std::size_t val : output_values_map[i]) {
        std::string partial_state = state;
        partial_state[pos] = std::to_string(val).at(0);
        temp_output_states.push_back(partial_state);
      }
    }
    *output_states = temp_output_states;
    temp_output_states.clear();
  }
  // Verify that output states have no leftover "x" after replacement.
  for (std::size_t i = 0; i < output_states->at(0).length(); ++i) {
    char c = output_states->at(0)[i];
    if (c != '0' && c != '1') {
      throw ERROR_MSG("Final state has non-binary character ", c, " at index ",
                      i, "despite having no terminal cut there.");
    }
  }
  return base_state;
}

std::vector<std::pair<std::string, std::complex<double>>> EvaluateCircuit(
    QflexInput* input) {
  if (input == nullptr) {
    throw ERROR_MSG("Input must be non-null.");
  }

  utils::Stopwatch stopwatch;

  // Start stopwatch
  if (global::verbose > 0) stopwatch.start();

  // Create the ordering for this tensor contraction from file.
  std::list<ContractionOperation> ordering;

  // Parse ordering
  ordering_data_to_contraction_ordering(*input, &ordering);

  if (global::verbose > 0)
    std::cerr << WARN_MSG("Time spent making contraction ordering: ",
                          stopwatch.split<utils::milliseconds>() / 1000., "s")
              << std::endl;

  std::size_t init_length =
      input->grid.I * input->grid.J - input->grid.qubits_off.size();
  if (input->initial_state.empty()) {
    input->initial_state = std::string(init_length, '0');
  }

  // Get a list of qubits and output states for the final region, and set the
  // final_state if one wasn't provided.
  std::vector<std::vector<std::size_t>> final_qubits;
  std::vector<std::string> output_states;

  try {
    input->final_state =
        get_output_states(input, ordering, &final_qubits, &output_states);
  } catch (const std::string& err_msg) {
    throw ERROR_MSG("Failed to call get_output_states(). Error:\n\t[", err_msg,
                    "]");
  }

  // Declaring and then filling 2D grid of tensors.
  std::vector<std::vector<Tensor>> tensor_grid(
      input->grid.I, std::vector<Tensor>(input->grid.J));

  if (global::verbose > 0)
    std::cerr << WARN_MSG("Time spent allocating 2D grid of tensors: ",
                          stopwatch.split<utils::milliseconds>() / 1000., "s")
              << std::endl;

  // Scope so that scratch space and the 3D grid of tensors are destructed.
  {
    // Scratch space for creating 3D tensor network. The largest single-gate
    // tensor we currently support is rank 4.
    s_type scratch_3D[16];

    // Creating 3D grid of tensors from file.
    std::vector<std::vector<std::vector<Tensor>>> tensor_grid_3D;
    circuit_data_to_tensor_network(input->circuit, input->grid.I, input->grid.J,
                                   input->initial_state, input->final_state,
                                   final_qubits, input->grid.qubits_off,
                                   tensor_grid_3D, scratch_3D);

    if (global::verbose > 0)
      std::cerr << WARN_MSG(
                       "Time spent creating 3D grid of tensors from file: ",
                       stopwatch.split<utils::milliseconds>() / 1000., "s")
                << std::endl;

    std::size_t max_size = 0;
    for (std::size_t i = 0; i < tensor_grid_3D.size(); ++i) {
      for (std::size_t j = 0; j < tensor_grid_3D[i].size(); ++j) {
        std::unordered_map<std::string, std::size_t> index_dim;
        for (const auto tensor : tensor_grid_3D[i][j]) {
          for (const auto& [index, dim] : tensor.get_index_to_dimension()) {
            if (index_dim.find(index) == index_dim.end()) {
              index_dim[index] = dim;
            } else {
              // Index is shared between adjacent tensors; remove it.
              index_dim.erase(index);
            }
          }
        }
        std::size_t tensor_size = 1;
        for (const auto& [index, dim] : index_dim) {
          tensor_size *= dim;
        }
        if (tensor_size > max_size) {
          max_size = tensor_size;
        }
      }
    }

    if (global::verbose > 0)
      std::cerr << WARN_MSG(
                       "Time spent to determine maximum size for tensors: ",
                       stopwatch.split<utils::milliseconds>() / 1000., "s")
                << std::endl;

    // Scratch space for contracting 3D to 2D grid. This must have enough space
    // to hold the largest single-qubit tensor in the 2D grid.
    std::vector<s_type> scratch_2D(max_size);

    if (global::verbose > 0)
      std::cerr << WARN_MSG("Time spent allocating scratch space for 2D grid: ",
                            stopwatch.split<utils::milliseconds>() / 1000., "s")
                << std::endl;

    // Contract 3D grid onto 2D grid of tensors, as usual.
    flatten_grid_of_tensors(tensor_grid_3D, tensor_grid, final_qubits,
                            input->grid.qubits_off, ordering,
                            scratch_2D.data());

    if (global::verbose > 0)
      std::cerr << WARN_MSG(
                       "Time spent creating 2D grid of tensors from 3D one: ",
                       stopwatch.split<utils::milliseconds>() / 1000., "s")
                << std::endl;
  }

  // Perform tensor grid contraction.
  std::vector<std::complex<double>> amplitudes(output_states.size());
  std::vector<std::pair<std::string, std::complex<double>>> result;
  try {
    ContractGrid(ordering, &tensor_grid, &amplitudes);
  } catch (const std::string& err_msg) {
    throw ERROR_MSG("Failed to call ContractGrid(). Error:\n\t[", err_msg, "]");
  }
  for (std::size_t c = 0; c < amplitudes.size(); ++c) {
    result.push_back(std::make_pair(output_states[c], amplitudes[c]));
  }

  if (global::verbose > 0)
    std::cerr << WARN_MSG("Time spent contracting tensors: ",
                          stopwatch.split<utils::milliseconds>() / 1000., "s")
              << std::endl;

  // Final time
  if (global::verbose > 0) {
    stopwatch.stop();
    std::cerr << WARN_MSG("Total time: ",
                          stopwatch.time_passed<utils::milliseconds>() / 1000.,
                          "s")
              << std::endl;
  }

  return result;
}

}  // namespace qflex
