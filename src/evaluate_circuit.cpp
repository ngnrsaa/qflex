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

namespace qflex {

void QflexGrid::load(std::istream& istream) {
  I = J = 0;
  std::string line;
  while (std::getline(istream, line))
    if (std::size(line) and line[0] != '#') {
      // String unnecessary chars
      line.erase(
          std::remove_if(std::begin(line), std::end(line),
                         [](auto&& x) { return not(x == '0' || x == '1'); }),
          std::end(line));

      // Continue if line is empty
      if (std::empty(line)) continue;

      // Get number of columns
      if (J != 0 and J != std::size(line))
        throw std::string("Grid size is inconsistent");
      else
        J = std::size(line);

      // Get off qubits
      for (int q = 0; q < J; ++q)
        if (line[q] == '0') qubits_off.push_back({I, q});

      // Update number of rows
      ++I;
    }
};

void QflexGrid::load(const std::string& filename) {
  if (auto in = std::ifstream(filename); in.good())
    this->load(in);
  else
    throw std::string("Cannot open grid file: ") + filename;
};

// Gets the position in the output state vector of the qubit at tensor_pos.
int find_output_pos(const QflexInput* input, std::vector<int> tensor_pos) {
  int pos = (tensor_pos[0] * input->grid.J) + tensor_pos[1];
  for (const auto off_pos : input->grid.qubits_off) {
    if (off_pos[0] < tensor_pos[0]) {
      --pos;
    } else if (off_pos[0] == tensor_pos[0] && off_pos[1] < tensor_pos[1]) {
      --pos;
    }
  }
  return pos;
}

std::string get_output_states(const QflexInput* input,
                              const std::list<ContractionOperation>& ordering,
                              std::vector<std::vector<int>>* final_qubits,
                              std::vector<std::string>* output_states) {
  if (final_qubits == nullptr) {
    std::cout << "Final qubits must be non-null." << std::endl;
    assert(final_qubits != nullptr);
  }
  std::vector<int> output_pos_map;
  std::vector<std::vector<int>> output_values_map;
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
    const int pos = find_output_pos(input, op.cut.tensors[0]);
    const auto tensor_pos = op.cut.tensors[0];
    if (final_state_unspecified) {
      base_state[pos] = 'x';
    }
    // TODO(martinop): reconsider requiring 'x' for cut indices.
    output_pos_map.push_back(pos);
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
  for (int i = 0; i < final_qubits->size(); ++i) {
    const int pos = output_pos_map[i];
    for (const std::string& state : *output_states) {
      for (const int val : output_values_map[i]) {
        std::string partial_state = state;
        partial_state[pos] = std::to_string(val).at(0);
        temp_output_states.push_back(partial_state);
      }
    }
    *output_states = temp_output_states;
    temp_output_states.clear();
  }
  // Verify that output states have no leftover "x" after replacement.
  for (int i = 0; i < output_states->at(0).length(); ++i) {
    char c = output_states->at(0)[i];
    if (c != '0' && c != '1') {
      std::cout << "Final state has non-binary character " << c << " at index "
                << i << "despite having no terminal cut there.";
      assert(c == '0' || c == '1');
    }
  }
  return base_state;
}

std::vector<std::pair<std::string, std::complex<double>>> EvaluateCircuit(
    QflexInput* input) {
  if (input == nullptr) {
    std::cout << "Input must be non-null." << std::endl;
    assert(input != nullptr);
  }
  // Set precision for the printed floats.
  std::cout.precision(12);

  // Timing variables.
  std::chrono::high_resolution_clock::time_point t_output_0, t_output_1;
  t_output_0 = std::chrono::high_resolution_clock::now();
  std::chrono::high_resolution_clock::time_point t0, t1;
  std::chrono::duration<double> time_span;

  // Reading input.
  const int super_dim = (int)pow(DIM, input->circuit.depth);

  // Create the ordering for this tensor contraction from file.
  t0 = std::chrono::high_resolution_clock::now();
  std::list<ContractionOperation> ordering;
  ordering_data_to_contraction_ordering(input->ordering_data, input->grid.I,
                                        input->grid.J, input->grid.qubits_off,
                                        &ordering);
  t1 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
  if (input->enable_timing_logs) {
    std::cout << "Time spent making contraction ordering: " << time_span.count()
              << "s" << std::endl;
  }

  int init_length =
      input->grid.I * input->grid.J - input->grid.qubits_off.size();
  if (input->initial_state.empty()) {
    input->initial_state = std::string(init_length, '0');
  }

  // Get a list of qubits and output states for the final region, and set the
  // final_state if one wasn't provided.
  std::vector<std::vector<int>> final_qubits;
  std::vector<std::string> output_states;
  input->final_state =
      get_output_states(input, ordering, &final_qubits, &output_states);

  // Scratch space to be reused for operations.
  t0 = std::chrono::high_resolution_clock::now();
  // This scratch space is used for reading circuit and building tensor
  // network. At most we need super_dim * 4 for square grids, and times 2
  // when qubits are cut on the output index.
  s_type* scratch = new s_type[(int)pow(super_dim, 4) * 2];
  t1 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
  if (input->enable_timing_logs) {
    std::cout << "Time spent reading allocating scratch space: "
              << time_span.count() << "s" << std::endl;
  }

  // Declaring and then filling 2D grid of tensors.
  std::vector<std::vector<Tensor>> tensor_grid(input->grid.I);
  for (int i = 0; i < input->grid.I; ++i) {
    tensor_grid[i] = std::vector<Tensor>(input->grid.J);
  }
  // Scope so that the 3D grid of tensors is destructed.
  {
    // Creating 3D grid of tensors from file.
    t0 = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<std::vector<Tensor>>> tensor_grid_3D;
    circuit_data_to_tensor_network(input->circuit, input->grid.I, input->grid.J,
                                   input->initial_state, input->final_state,
                                   final_qubits, input->grid.qubits_off,
                                   tensor_grid_3D, scratch);

    t1 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    if (input->enable_timing_logs) {
      std::cout << "Time spent creating 3D grid of tensors from file: "
                << time_span.count() << "s" << std::endl;
    }

    // Contract 3D grid onto 2D grid of tensors, as usual.
    t0 = std::chrono::high_resolution_clock::now();
    flatten_grid_of_tensors(tensor_grid_3D, tensor_grid, final_qubits,
                            input->grid.qubits_off, ordering, scratch);
    t1 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    if (input->enable_timing_logs) {
      std::cout << "Time spent creating 2D grid of tensors from 3D one: "
                << time_span.count() << "s" << std::endl;
    }

    // Freeing scratch data: delete and NULL.
    delete[] scratch;
    scratch = NULL;
  }

  // Perform tensor grid contraction.
  std::vector<std::complex<double>> amplitudes(output_states.size());
  std::vector<std::pair<std::string, std::complex<double>>> result;
  ContractGrid(ordering, &tensor_grid, &amplitudes);
  for (int c = 0; c < amplitudes.size(); ++c) {
    result.push_back(std::make_pair(output_states[c], amplitudes[c]));
  }

  // Final time
  t_output_1 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
      t_output_1 - t_output_0);
  std::cout << "Total time: " << time_span.count() << "s" << std::endl;

  return result;
}

}  // namespace qflex
