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

#include <set>

namespace qflex {

// Gets the position in the output state vector of the qubit at tensor_pos.
std::size_t find_output_pos(const QflexGrid& grid,
                            std::vector<std::size_t> tensor_pos) {
  std::size_t pos = (tensor_pos[0] * grid.J) + tensor_pos[1];
  for (const auto off_pos : grid.qubits_off) {
    if (off_pos[0] < tensor_pos[0]) {
      --pos;
    } else if (off_pos[0] == tensor_pos[0] && off_pos[1] < tensor_pos[1]) {
      --pos;
    }
  }
  return pos;
}

QflexFinalQubits get_final_qubits(
    const QflexGrid& grid, const std::list<ContractionOperation>& ordering) {
  std::vector<std::vector<std::size_t>> final_qubits;
  std::vector<std::size_t> output_pos_map;
  std::vector<std::vector<std::size_t>> output_values_map;

  // If the final state isn't provided, it should be all zeroes except for
  // qubits with terminal cuts (which should have 'x').
  for (const auto& op : ordering) {
    // TODO(martinop): update to use the new operation.
    if (op.op_type != ContractionOperation::CUT) continue;

    // Any qubit with a terminal cut is in the final region.
    if (op.cut.tensors.size() != 1) continue;
    const std::size_t output_pos = find_output_pos(grid, op.cut.tensors[0]);
    const auto tensor_pos = op.cut.tensors[0];

    // TODO(martinop): reconsider requiring 'x' for cut indices.
    output_pos_map.push_back(output_pos);
    if (op.cut.values.empty()) {
      output_values_map.push_back({0, 1});
    } else {
      output_values_map.push_back(op.cut.values);
    }
    final_qubits.push_back(tensor_pos);
  }

  return {final_qubits, output_pos_map, output_values_map};
}

std::vector<std::string> get_output_states(
    const std::string& base_state, const QflexFinalQubits& final_qubits) {
  std::vector<std::string> output_states;

  // Construct the full set of output state strings.
  std::vector<std::string> temp_output_states;
  output_states.push_back(base_state);
  for (std::size_t i = 0; i < std::size(final_qubits.qubits); ++i) {
    const std::size_t pos = final_qubits.output_pos_map[i];
    for (const std::string& state : output_states) {
      for (const std::size_t val : final_qubits.output_values_map[i]) {
        std::string partial_state = state;
        partial_state[pos] = std::to_string(val).at(0);
        temp_output_states.push_back(partial_state);
      }
    }
    output_states = temp_output_states;
    temp_output_states.clear();
  }
  // Verify that output states have no leftover "x" after replacement.
  for (std::size_t i = 0; i < output_states.at(0).length(); ++i) {
    char c = output_states.at(0)[i];
    if (c != '0' && c != '1') {
      throw ERROR_MSG("Final state has non-binary character ", c, " at index ",
                      i, "despite having no terminal cut there.");
    }
  }

  return output_states;
}

std::vector<std::pair<std::string, std::complex<double>>> EvaluateCircuit(
    const QflexInput& input) {
  std::chrono::steady_clock::time_point t_output_0;
  std::chrono::steady_clock::time_point t0, t1;

  // Print difference in time
  auto print_diff_time = [&t0, &t1](const std::string& msg) {
    const double delta_time_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            (t1 = std::chrono::steady_clock::now()) - t0)
            .count() /
        1000.;
    std::cerr << WARN_MSG(msg, delta_time_ms, "s") << std::endl;
    t0 = std::move(t1);
  };

  // First initialization
  if (global::verbose > 0) {
    // Timing variables.
    t_output_0 = std::chrono::steady_clock::now();
    t0 = std::chrono::steady_clock::now();
  }

  // Create the ordering for this tensor contraction from file.
  std::list<ContractionOperation> ordering;

  // Parse ordering
  ordering_data_to_contraction_ordering(input, &ordering);

  if (global::verbose > 0)
    print_diff_time("Time spent making contraction ordering: ");

  // Get a list of qubits and output states for the final region.
  const auto final_qubits = get_final_qubits(input.grid, ordering);

  std::vector<std::vector<std::vector<Tensor>>> tensor_grid_3D;
  std::vector<std::vector<Tensor>> tensor_grid_orig;
  std::vector<s_type> scratch_2D;

  // Scope so that scratch space and the 3D grid of tensors are destructed.
  {
    // Scratch space for creating 3D tensor network. The largest single-gate
    // tensor we currently support is rank 4.
    s_type scratch_3D[16];

    // Creating 3D grid of tensors from file.
    circuit_data_to_tensor_network(input.circuit, input.grid.I, input.grid.J,
                                   input.initial_state, input.grid.qubits_off,
                                   tensor_grid_3D, scratch_3D);

    if (global::verbose > 0)
      print_diff_time("Time spent creating 3D grid of tensors from file: ");

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
        for (const auto& [index, dim] : index_dim)
          if (tensor_size *= dim; max_size < tensor_size)
            max_size = tensor_size;
      }
    }

    // Scratch space for contracting 3D to 2D grid. This must have enough space
    // to hold the largest single-qubit tensor in the 2D grid.
    scratch_2D.resize(max_size);

    if (global::verbose > 0)
      print_diff_time("Time spent allocating scratch space for 2D grid: ");

    // Rename last index when in final_qubit_region to
    // "(i,j),(o)".
    {
      std::size_t idx = 0;
      for (std::size_t i = 0; i < input.grid.I; ++i)
        for (std::size_t j = 0; j < input.grid.J; ++j) {
          // Skip if (i,j) is off
          if (find_grid_coord_in_list(input.grid.qubits_off, i, j)) continue;

          std::string last_name = utils::concat(
              "(", i, ",", j, "),(", std::size(tensor_grid_3D[i][j]) - 1, ")");

          if (find_grid_coord_in_list(final_qubits.qubits, i, j)) {
            std::string output_name = utils::concat("(", i, ",", j, "),(o)");

            try {
              tensor_grid_3D[i][j].back().rename_index(last_name, output_name);
            } catch (const std::string& err_msg) {
              throw ERROR_MSG("Failed to call rename_index(). Error:\n\t[",
                              err_msg, "]");
            }
          }

          ++idx;
        }
    }

    // Contract 3D grid onto 2D grid of tensors, as usual. At this point,
    // tensor_grid should be independent of the given final_state
    tensor_grid_orig = flatten_grid_of_tensors(
        tensor_grid_3D, input.grid.qubits_off, scratch_2D.data());

    if (global::verbose > 0)
      print_diff_time("Time spent creating 2D grid of tensors from 3D one: ");
  }

  // Perform contraction
  std::vector<std::pair<std::string, std::complex<double>>> result;

  for (std::size_t k = std::size(input.final_states); k--;) {
    const std::string final_state = input.final_states[k];

    std::vector<std::vector<Tensor>> tensor_grid;
    tensor_grid = k > 0 ? tensor_grid_orig : std::move(tensor_grid_orig);

    if (global::verbose > 0) print_diff_time("Time spent copying: ");

    // Insert deltas to last layer on qubits that are in not in
    // final_qubit_region.
    {
      std::size_t idx = 0;
      for (std::size_t i = 0; i < input.grid.I; ++i)
        for (std::size_t j = 0; j < input.grid.J; ++j) {
          // Skip if (i,j) is off
          if (find_grid_coord_in_list(input.grid.qubits_off, i, j)) continue;

          std::string last_name = utils::concat(
              "(", i, ",", j, "),(", std::size(tensor_grid_3D[i][j]) - 1, ")");

          if (!find_grid_coord_in_list(final_qubits.qubits, i, j)) {
            std::string delta_gate =
                (final_state[idx] == '0') ? "delta_0" : "delta_1";

            Tensor& A = tensor_grid[i][j];
            Tensor B = Tensor({last_name}, {2}, gate_array(delta_gate, {}));
            Tensor C({""}, {result_size(A, B)});
            try {
              multiply(A, B, C, scratch_2D.data());
            } catch (const std::string& err_msg) {
              throw ERROR_MSG("Failed to call multiply(). Error:\n\t[", err_msg,
                              "]");
            }

            // Move the result of A*B --> A
            A = std::move(C);
          }

          ++idx;
        }
    }

    // Reorder the 2D grid
    reorder_grid_of_tensors(tensor_grid, final_qubits.qubits,
                            input.grid.qubits_off, ordering, scratch_2D.data());

    if (global::verbose > 0)
      print_diff_time("Time spent applying final state: ");

    // Perform tensor grid contraction.
    const auto output_states = get_output_states(final_state, final_qubits);
    std::vector<std::complex<double>> amplitudes(output_states.size());
    try {
      ContractGrid(ordering, &tensor_grid, &amplitudes);
    } catch (const std::string& err_msg) {
      throw ERROR_MSG("Failed to call ContractGrid(). Error:\n\t[", err_msg,
                      "]");
    }
    for (std::size_t c = 0; c < amplitudes.size(); ++c) {
      result.push_back(std::make_pair(output_states[c], amplitudes[c]));
    }

    // Contracting time
    if (global::verbose > 0) print_diff_time("Time spent in contracting: ");
  }

  // Final time
  if (global::verbose > 0) {
    t1 = std::chrono::steady_clock::now();
    const auto time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(t1 -
                                                                  t_output_0);
    std::cerr << WARN_MSG("Total time: ", time_span.count(), "s") << std::endl;
  }

  return result;
}

}  // namespace qflex
