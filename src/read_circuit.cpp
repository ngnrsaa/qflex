/**
 * @file read_circuit.cpp
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

#include "read_circuit.h"

#include <memory>

#include "errors.h"

namespace qflex {
// TODO: when building tensor networks and pulling gates from gate_array(s),
// stop relying on dimensions given by hand. They should all rely on dimensions
// either stored in gate_arrays or on DIM, which is 2 as a global variable.

const std::unordered_map<std::string, std::vector<s_type>> _GATES_DATA(
    {// Deltas.
     {"delta_0", std::vector<s_type>({1.0, 0.0})},
     {"delta_1", std::vector<s_type>({0.0, 1.0})},
     // For one-qubit gates, the first index is input an second is output.
     {"h", std::vector<s_type>({static_cast<s_type::value_type>(M_SQRT1_2),
                                static_cast<s_type::value_type>(M_SQRT1_2),
                                static_cast<s_type::value_type>(M_SQRT1_2),
                                -static_cast<s_type::value_type>(M_SQRT1_2)})},
     {"hz_1_2",
      std::vector<s_type>({{0.5, 0.5},
                           {static_cast<s_type::value_type>(M_SQRT1_2), 0},
                           {0., -static_cast<s_type::value_type>(M_SQRT1_2)},
                           {0.5, 0.5}})},
     {"t", std::vector<s_type>({1.0,
                                0.,
                                0.,
                                {static_cast<s_type::value_type>(M_SQRT1_2),
                                 static_cast<s_type::value_type>(M_SQRT1_2)}})},
     {"x_1_2",
      std::vector<s_type>({{0.5, 0.5}, {0.5, -0.5}, {0.5, -0.5}, {0.5, 0.5}})},
     {"y_1_2",
      std::vector<s_type>({{0.5, 0.5}, {0.5, 0.5}, {-0.5, -0.5}, {0.5, 0.5}})},
     // For cz, both q1 and q2 get indices in the order (input, virtual,
     // output).
     // {"cz_q1", vector<s_type>({1.,0.,0.,0.,0.,0.,0.,1.})},
     // {"cz_q2", vector<s_type>({1.,0.,1.,0.,0.,1.,0.,-1.})}});
     // Use more "balanced" SVD. Otherwise tensors are very sparse.
     {"cz_q1",
      std::vector<s_type>({-0.3446133714352872, 0., 1.1381806476131544, 0., 0.,
                           -1.1381806476131544, 0., -0.3446133714352872})},
     {"cz_q2",
      std::vector<s_type>({-1.0484937059720079, 0., 0.5611368023131075, 0., 0.,
                           0.5611368023131075, 0., 1.0484937059720079})},
     // For the non-decomposed cz, the convention is (in1, in2, out1, out2).
     {"cz", std::vector<s_type>({1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 1., 0.,
                                 0., 0., 0., -1.})},

     // For cx, both q1 and q2 get indices in the order (input, virtual,
     // output).
     //{"cx_q1", std::vector<s_type>({1.,0.,0.,0.,0.,0.,0.,1.})},
     //{"cx_q2", std::vector<s_type>({1.,0.,0.,1.,0.,1.,1.,0.})},
     // Use more "balanced" SVD. Otherwise tensors are very sparse.
     {"cx_q1",
      std::vector<s_type>({0.8408964152537143, 0., 0.8408964152537143, 0., 0.,
                           -0.8408964152537143, 0., 0.8408964152537143})},
     {"cx_q2", std::vector<s_type>({0.5946035575013604, -0.5946035575013604,
                                    0.5946035575013604, 0.5946035575013604,
                                    -0.5946035575013604, 0.5946035575013604,
                                    0.5946035575013604, 0.5946035575013604})},
     // For the non-decomposed cx, the convention is (in1, in2, out1, out2).
     {"cx", std::vector<s_type>({1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 1.,
                                 0., 0., 1., 0.})}});

// TODO: make params optional.
std::vector<s_type> gate_array(const std::string& gate_name,
                               const std::vector<double>& params) {
  if (gate_name == "rz") {
    const double angle_rads = M_PI * params[0];
    const s_type phase = std::exp(s_type(0., angle_rads / 2));
    return std::vector<s_type>({conj(phase), 0., 0., phase});
  }

  bool valid_gate = _GATES_DATA.find(gate_name) != _GATES_DATA.end();
  if (!valid_gate) {
    throw ERROR_MSG("Invalid gate name provided: ", gate_name);
  }
  return _GATES_DATA.at(gate_name);
}

/**
 * Get fSim decomposed tensors given angles theta and phi, in units of PI.
 * @param double theta_rads with the first angle of the fSim gate.
 * @param double phi_rads with the second angle of the fSim gate.
 * @return vector<vector<stype>> with the two tensors of the fSim gate in
 * row major storage and order (input, virtual, output) for indices on both
 * tensors.
 */
std::vector<std::vector<s_type>> fSim(s_type::value_type theta_rads,
                                      s_type::value_type phi_rads) {
  s_type c1, c2, c3, c4, d;
  c1 = s_type(0.5) * (std::exp(s_type(0.0, -1 * phi_rads / 2)) +
                      s_type(std::cos(theta_rads), 0.0));
  c2 = s_type(0.5) * (std::exp(s_type(0.0, -1 * phi_rads / 2)) -
                      s_type(std::cos(theta_rads), 0.0));
  c3 = s_type(0.0, -1 / 2.) * s_type(std::sin(theta_rads), 0.0);
  c4 = c3;
  s_type s1, s2, s3, s4;
  s1 = std::sqrt(c1);
  s2 = std::sqrt(c2);
  s3 = std::sqrt(c3);
  s4 = std::sqrt(c4);
  d = std::exp(s_type(0.0, phi_rads / 4.));
  std::vector<s_type> q1_tensor_array({d * s1,
                                       {0.0, 0.0},
                                       d * s2,
                                       {0.0, 0.0},
                                       {0.0, 0.0},
                                       s3,
                                       {0.0, 0.0},
                                       s_type(0.0, 1.0) * s4,
                                       {0.0, 0.0},
                                       conj(d) * s1,
                                       {0.0, 0.0},
                                       -conj(d) * s2,
                                       s3,
                                       {0.0, 0.0},
                                       s_type(0.0, -1.0) * s4,
                                       {0.0, 0.0}});

  std::vector<s_type> q2_tensor_array(q1_tensor_array);
  return std::vector<std::vector<s_type>>({q1_tensor_array, q2_tensor_array});
}

// TODO: Refactor this so that all gate arrays are handled similarly regardless
// of how many qubits they operate on and whether they're parametrized. Put gate
// array handling into its own class.
// TODO: make params optional.
// TODO: implement c-phase, cpf(phi).
std::tuple<std::vector<s_type>, std::vector<s_type>, std::vector<std::size_t>>
gate_arrays(const std::string& gate_name, const std::vector<double>& params) {
  if (gate_name == "cz") {
    return std::tuple<std::vector<s_type>, std::vector<s_type>,
                      std::vector<std::size_t>>(
        gate_array("cz_q1", params), gate_array("cz_q2", params), {2, 2, 2});
  } else if (gate_name == "cx") {
    return std::tuple<std::vector<s_type>, std::vector<s_type>,
                      std::vector<std::size_t>>(
        gate_array("cx_q1", params), gate_array("cx_q2", params), {2, 2, 2});
  } else if (gate_name == "fsim") {
    const double theta_rads = M_PI * params[0];
    const double phi_rads = M_PI * params[1];
    std::vector<std::vector<s_type>> ret_val = fSim(theta_rads, phi_rads);
    return std::tuple<std::vector<s_type>, std::vector<s_type>,
                      std::vector<std::size_t>>(ret_val[0], ret_val[1],
                                                {2, 4, 2});
  }
  throw ERROR_MSG("Invalid gate name provided: ", gate_name);
}

/**
 * Helper method for grid_of_tensors_3D_to_2D which generates a comparator
 * function for sorting the indices of a tensor at grid position "local" into
 * the order prescribed by "ordering".
 *
 * The resulting comparator takes two grid positions adjacent to "local" as
 * input ("lhs" and "rhs") and determines the relative order in which the
 * tensors at each position are contracted with the tensor at "local". The
 * comparator returns:
 *   - true if the lhs tensor is contracted with the local tensor before the
 *     rhs tensor is contracted with the local tensor
 *   - false if the rhs tensor is contracted with the local tensor before the
 *     lhs tensor is contracted with the local tensor
 *   - error if either the lhs or rhs tensor is never contracted with the local
 *     tensor (this usually indicates an issue in the provided ordering)
 *
 * Generating and applying this comparator for all grid tensors ensures that
 * contraction (using the order given by "ordering") will never require further
 * reordering of tensor indices.
 * @param ordering std::list<ContractionOperation> providing the steps required
 * to contract the tensor grid.
 * @param local vector<std::size_t> of the target qubit coordinates.
 * @return function for use as a comparator in std::sort.
 */
std::function<bool(std::vector<std::size_t>, std::vector<std::size_t>)>
order_func(const std::list<ContractionOperation>& ordering,
           std::vector<std::size_t> local) {
  return [&ordering, local](const std::vector<std::size_t> lhs,
                            const std::vector<std::size_t> rhs) {
    // Cuts are projected early and should be ordered first.
    std::vector<std::vector<std::size_t>> lhs_pair, rhs_pair;
    if (local[0] < lhs[0] || local[1] < lhs[1]) {
      lhs_pair = {local, lhs};
    } else {
      lhs_pair = {lhs, local};
    }
    if (local[0] < rhs[0] || local[1] < rhs[1]) {
      rhs_pair = {local, rhs};
    } else {
      rhs_pair = {rhs, local};
    }
    for (const auto& op : ordering) {
      if (op.op_type != ContractionOperation::CUT) continue;
      if (lhs_pair == op.cut.tensors) return true;
      if (rhs_pair == op.cut.tensors) return false;
    }

    std::string lpatch = "null";
    std::string rpatch = "null";
    std::optional<std::size_t> lpos;
    std::optional<std::size_t> rpos;
    std::size_t op_num = 0;
    for (const auto& op : ordering) {
      if (op.op_type != ContractionOperation::EXPAND) continue;
      if (lhs == op.expand.tensor) {
        lpatch = op.expand.id;
        lpos = op_num;
        break;
      }
      op_num++;
    }
    if (not lpos.has_value()) {
      throw ERROR_MSG("Left hand side of pair not found: (", local[0], ',',
                      local[1], "),(", lhs[0], ',', lhs[1], ").");
    }

    op_num = 0;
    for (const auto& op : ordering) {
      if (op.op_type != ContractionOperation::EXPAND) continue;
      if (rhs == op.expand.tensor) {
        rpatch = op.expand.id;
        rpos = op_num;
        break;
      }
      op_num++;
    }
    if (not rpos.has_value()) {
      throw ERROR_MSG("Right hand side of pair not found: (", local[0], ',',
                      local[1], "),(", rhs[0], ',', rhs[1], ").");
    }
    if (lpatch == rpatch) {
      // lhs and rhs are in the same patch.
      return lpos < rpos;
    }

    std::string local_patch;
    for (const auto& op : ordering) {
      if (op.op_type != ContractionOperation::EXPAND) continue;
      if (local == op.expand.tensor) {
        local_patch = op.expand.id;
        break;
      }
    }
    if (lpatch == local_patch) {
      // rhs won't be connected until patches merge.
      return true;
    }
    if (rpatch == local_patch) {
      // lhs won't be connected until patches merge.
      return false;
    }
    // Both lhs and rhs are in different patches from local_patch; find out
    // which merges with the local patch first.
    for (const auto& op : ordering) {
      if (op.op_type != ContractionOperation::MERGE) continue;
      if (local_patch == op.merge.source_id) {
        if (lpatch == op.merge.target_id) return true;
        if (rpatch == op.merge.target_id) return false;
        local_patch = op.merge.target_id;
      } else if (local_patch == op.merge.target_id) {
        if (lpatch == op.merge.source_id) return true;
        if (rpatch == op.merge.source_id) return false;
        // local_patch is already the target ID.
      } else if (lpatch == op.merge.source_id) {
        lpatch = op.merge.target_id;
      } else if (rpatch == op.merge.source_id) {
        rpatch = op.merge.target_id;
      }
    }
    // Error in comparison - likely issue in contraction ordering.
    throw ERROR_MSG("Failed to compare (", lhs[0], ',', lhs[1], ") and (",
                    rhs[0], ',', rhs[1], ") for local (", local[0], ',',
                    local[1], ").");
  };
}

// TODO: add tests for this function. Compactify code. Use index_name()
// function for all index names used here. Use smart pointers where possible.
// Add depth functionality to read a circuit up to a certain cycle.
void circuit_data_to_tensor_network(
    const QflexCircuit& circuit, std::size_t I, std::size_t J,
    const std::string initial_conf, const std::string final_conf,
    const std::optional<std::vector<std::vector<std::size_t>>>&
        final_qubit_region,
    const std::optional<std::vector<std::vector<std::size_t>>>& off,
    std::vector<std::vector<std::vector<Tensor>>>& grid_of_tensors,
    s_type* scratch) {
  if (scratch == nullptr) {
    throw ERROR_MSG("Scratch must be non-null.");
  }
  // Useful for plugging into the tensor network:
  std::vector<std::size_t> i_j_1, i_j_2;
  // Calculated from input.
  const std::size_t grid_size = I * J;
  const std::size_t off_size = off.has_value() ? off.value().size() : 0;
  const std::size_t num_active_qubits_from_grid = grid_size - off_size;

  if (circuit.num_active_qubits != num_active_qubits_from_grid) {
    throw ERROR_MSG(
        "The number of active qubits read from the file: ",
        circuit.num_active_qubits,
        ", does not match the number of active qubits read from the grid: ",
        num_active_qubits_from_grid, ".");
  }

  // Check for the length of initial_conf and final_conf.
  {
    // std::size_t off_size = off.has_value() ? off.value().size() : 0;
    if (initial_conf.size() != num_active_qubits_from_grid) {
      throw ERROR_MSG("Size of initial_conf: ", initial_conf.size(),
                      ", must be equal to the number of qubits: ",
                      num_active_qubits_from_grid, ".");
    }
    if (final_conf.size() != initial_conf.size()) {
      throw ERROR_MSG("Size of final_conf: ", final_conf.size(),
                      ", must be equal to size of initial_conf: ",
                      initial_conf.size(), ".");
    }
  }

  // Creating grid variables.
  grid_of_tensors = std::vector<std::vector<std::vector<Tensor>>>(I);
  std::vector<std::vector<std::size_t>> grid_of_counters(I);
  std::unordered_map<std::string, std::size_t> link_counters;
  for (std::size_t i = 0; i < I; ++i) {
    grid_of_tensors[i] = std::vector<std::vector<Tensor>>(J);
    grid_of_counters[i] = std::vector<std::size_t>(J);
    for (std::size_t j = 0; j < J; ++j) {
      grid_of_tensors[i][j] = std::vector<Tensor>();
      grid_of_counters[i][j] = 0;
    }
  }

  // Insert deltas to first layer.
  std::size_t idx = 0;
  for (std::size_t q = 0; q < grid_size; ++q) {
    std::vector<std::size_t> i_j = get_qubit_coords(q, J);
    std::size_t i = i_j[0], j = i_j[1];
    if (find_grid_coord_in_list(off, i, j)) {
      continue;
    }
    std::string delta_gate = (initial_conf[idx] == '0') ? "delta_0" : "delta_1";
    std::string output_name = "(" + std::to_string(i_j[0]) + "," +
                              std::to_string(i_j[1]) + "),(" +
                              std::to_string(grid_of_counters[i][j]) + ")";
    try {
      grid_of_tensors[i][j].push_back(
          Tensor({output_name}, {2}, gate_array(delta_gate, {})));
    } catch (const std::string& err_msg) {
      throw ERROR_MSG("Failed to call Tensor(). Error:\n\t[", err_msg, "]");
    }
    idx += 1;
  }

  // Check if off contains valid qubits
  if (off.has_value())
    for (const auto& w : off.value()) {
      const auto& x = w[0];
      const auto& y = w[1];
      if (x >= I or y >= J)
        throw ERROR_MSG("Off qubit '(", x, ", ", y, ")' is outside the grid.");
    }

  std::unordered_set<std::size_t> used_qubits;
  std::size_t last_cycle{0};
  for (auto gate : circuit.gates) {
    // gate.name = gate name
    // gate.cycle = gate cycle
    // gate.qubits = vector of qubits
    // gate.params = vector of params

    if (last_cycle != gate.cycle) {
      last_cycle = gate.cycle;
      used_qubits.clear();
    }

    // Check that qubits haven't been used in cycle yet
    for (const auto& q : gate.qubits)
      if (used_qubits.find(q) != std::end(used_qubits))
        throw ERROR_MSG("Qubit ", q,
                        " has been used twice in the same cycle: ", gate.cycle);
      else
        used_qubits.insert(q);

    // One qubit gate
    if (std::size_t num_qubits = std::size(gate.qubits); num_qubits == 1) {
      // Get qubit
      std::size_t q1 = gate.qubits[0];
      i_j_1 = get_qubit_coords(q1, J);

      // Check that position is an active qubit
      bool qubit_off = find_grid_coord_in_list(off, i_j_1[0], i_j_1[1]);
      if (qubit_off)
        throw ERROR_MSG("Qubit '", q1, "' in gate '", gate.raw,
                        "' must correspond to an active qubit.");

      std::string input_name =
          "(" + std::to_string(i_j_1[0]) + "," + std::to_string(i_j_1[1]) +
          "),(" + std::to_string(grid_of_counters[i_j_1[0]][i_j_1[1]]) + ")";
      ++grid_of_counters[i_j_1[0]][i_j_1[1]];
      std::string output_name =
          "(" + std::to_string(i_j_1[0]) + "," + std::to_string(i_j_1[1]) +
          "),(" + std::to_string(grid_of_counters[i_j_1[0]][i_j_1[1]]) + ")";
      try {
        grid_of_tensors[i_j_1[0]][i_j_1[1]].push_back(
            Tensor({input_name, output_name}, {2, 2},
                   gate_array(gate.name, gate.params)));
      } catch (const std::string& err_msg) {
        throw ERROR_MSG("Failed to call Tensor(). Error:\n\t[", err_msg, "]");
      }

      // Two qubit gate
    } else if (num_qubits == 2) {
      // Get qubits
      std::size_t q1 = gate.qubits[0];
      std::size_t q2 = gate.qubits[1];
      i_j_1 = get_qubit_coords(q1, J);
      i_j_2 = get_qubit_coords(q2, J);

      // Check that positions are active qubits
      bool first_qubit_off = find_grid_coord_in_list(off, i_j_1[0], i_j_1[1]);
      bool second_qubit_off = find_grid_coord_in_list(off, i_j_2[0], i_j_2[1]);

      if (first_qubit_off)
        throw ERROR_MSG("Qubit '", q1, "' in gate '", gate.raw,
                        "' must correspond to an active qubit.");
      if (second_qubit_off)
        throw ERROR_MSG("Qubit '", q2, "' in gate '", gate.raw,
                        "' must correspond to an active qubit.");

      {
        std::size_t x_dist =
            std::max(i_j_1[0], i_j_2[0]) - std::min(i_j_1[0], i_j_2[0]);
        std::size_t y_dist =
            std::max(i_j_1[1], i_j_2[1]) - std::min(i_j_1[1], i_j_2[1]);
        if (x_dist + y_dist != 1)
          throw ERROR_MSG("Qubits ", q1, " and ", q2,
                          " are not nearest neighbors.");
      }

      std::vector<s_type> gate_q1;
      std::vector<s_type> gate_q2;
      std::vector<std::size_t> dimensions;
      tie(gate_q1, gate_q2, dimensions) = gate_arrays(gate.name, gate.params);
      std::string link_name;
      try {
        link_name = index_name(i_j_1, i_j_2);
      } catch (const std::string& err_msg) {
        throw ERROR_MSG("Failed to call index_name(). Error:\n\t[", err_msg,
                        "]");
      }
      // If no error is caught, link_name will be initialized.
      link_counters[link_name]++;
      std::size_t counter = link_counters[link_name];
      std::string virtual_name =
          "(" + std::to_string(std::min(i_j_1[0], i_j_2[0])) + "," +
          std::to_string(std::min(i_j_1[1], i_j_2[1])) + "," +
          std::to_string(counter) + "),(" +
          std::to_string(std::max(i_j_1[0], i_j_2[0])) + "," +
          std::to_string(std::max(i_j_1[1], i_j_2[1])) + "," +
          std::to_string(counter) + ")";
      std::string input_name_1 =
          "(" + std::to_string(i_j_1[0]) + "," + std::to_string(i_j_1[1]) +
          "),(" + std::to_string(grid_of_counters[i_j_1[0]][i_j_1[1]]) + ")";
      std::string input_name_2 =
          "(" + std::to_string(i_j_2[0]) + "," + std::to_string(i_j_2[1]) +
          "),(" + std::to_string(grid_of_counters[i_j_2[0]][i_j_2[1]]) + ")";
      ++grid_of_counters[i_j_1[0]][i_j_1[1]];
      ++grid_of_counters[i_j_2[0]][i_j_2[1]];
      std::string output_name_1 =
          "(" + std::to_string(i_j_1[0]) + "," + std::to_string(i_j_1[1]) +
          "),(" + std::to_string(grid_of_counters[i_j_1[0]][i_j_1[1]]) + ")";
      std::string output_name_2 =
          "(" + std::to_string(i_j_2[0]) + "," + std::to_string(i_j_2[1]) +
          "),(" + std::to_string(grid_of_counters[i_j_2[0]][i_j_2[1]]) + ")";
      try {
        grid_of_tensors[i_j_1[0]][i_j_1[1]].push_back(Tensor(
            {input_name_1, virtual_name, output_name_1}, dimensions, gate_q1));
      } catch (const std::string& err_msg) {
        throw ERROR_MSG("Failed to call Tensor(). Error:\n\t[", err_msg, "]");
      }
      try {
        grid_of_tensors[i_j_2[0]][i_j_2[1]].push_back(Tensor(
            {input_name_2, virtual_name, output_name_2}, dimensions, gate_q2));
      } catch (const std::string& err_msg) {
        throw ERROR_MSG("Failed to call Tensor(). Error:\n\t[", err_msg, "]");
      }
    } else
      throw ERROR_MSG("k-qubit gates with k > 2 not yet implemented.");
  }

  // Insert deltas to last layer on qubits that are in not in
  // final_qubit_region. Rename last index when in final_qubit_region to
  // "(i,j),(o)".
  idx = 0;
  for (std::size_t q = 0; q < grid_size; ++q) {
    std::vector<std::size_t> i_j = get_qubit_coords(q, J);
    std::size_t i = i_j[0], j = i_j[1];
    if (find_grid_coord_in_list(off, i, j)) {
      continue;
    }
    std::string last_name = "(" + std::to_string(i_j[0]) + "," +
                            std::to_string(i_j[1]) + "),(" +
                            std::to_string(grid_of_counters[i][j]) + ")";
    if (find_grid_coord_in_list(final_qubit_region, i, j)) {
      std::string output_name =
          "(" + std::to_string(i_j[0]) + "," + std::to_string(i_j[1]) + "),(o)";
      try {
        grid_of_tensors[i][j].back().rename_index(last_name, output_name);
      } catch (const std::string& err_msg) {
        throw ERROR_MSG("Failed to call rename_index(). Error:\n\t[", err_msg,
                        "]");
      }
    } else {
      std::string delta_gate = (final_conf[idx] == '0') ? "delta_0" : "delta_1";
      try {
        grid_of_tensors[i][j].push_back(
            Tensor({last_name}, {2}, gate_array(delta_gate, {})));
      } catch (const std::string& err_msg) {
        throw ERROR_MSG("Failed to call Tensor(). Error:\n\t[", err_msg, "]");
      }
    }
    idx += 1;
  }

  // Be proper about pointers.
  scratch = NULL;
}

// TODO: add tests for this function. Compactify code. Use smart pointers
// where possible. Improved contraction procedure.
void flatten_grid_of_tensors(
    std::vector<std::vector<std::vector<Tensor>>>& grid_of_tensors,
    std::vector<std::vector<Tensor>>& grid_of_tensors_2D,
    const std::optional<std::vector<std::vector<std::size_t>>>&
        final_qubit_region,
    const std::optional<std::vector<std::vector<std::size_t>>>& off,
    const std::list<ContractionOperation>& ordering, s_type* scratch) {
  if (scratch == nullptr) {
    throw ERROR_MSG("Scratch must be non-null.");
  }

  for (std::size_t i = 0, I = std::size(grid_of_tensors); i < I; ++i) {
    for (std::size_t j = 0, J = std::size(grid_of_tensors[i]); j < J; ++j) {
      // Skip if (i,j) doesn't belong to layout
      if (find_grid_coord_in_list(off, i, j)) continue;

      // Let's check that the K direction is not empty at this point
      if (std::empty(grid_of_tensors[i][j]))
        throw ERROR_MSG("Time-direction cannot be empty");

      // Get first tensor in time direction
      Tensor A = grid_of_tensors[i][j][0];

      for (std::size_t k = 1, K = grid_of_tensors[i][j].size(); k < K; ++k) {
        Tensor& B = grid_of_tensors[i][j][k];
        Tensor C({""}, {result_size(A, B)});
        try {
          multiply(A, B, C, scratch);
        } catch (const std::string& err_msg) {
          throw ERROR_MSG("Failed to call multiply(). Error:\n\t[", err_msg,
                          "]");
        }

        // Move the result of A*B --> A
        A = std::move(C);
      }

      // Move the final tensor to 2D grid of tensors
      grid_of_tensors_2D[i][j] = std::move(A);
    }
  }

  for (std::size_t i = 0, I = std::size(grid_of_tensors); i < I; ++i) {
    for (std::size_t j = 0, J = std::size(grid_of_tensors[i]); j < J; ++j) {
      if (find_grid_coord_in_list(off, i, j)) {
        continue;
      }

      // Positions of connected active qubits.
      std::vector<std::size_t> local({i, j});
      std::vector<std::vector<std::size_t>> pairs;

      if (i < I - 1 && !find_grid_coord_in_list(off, i + 1, j)) {
        pairs.push_back({i + 1, j});
      }
      if (j < J - 1 && !find_grid_coord_in_list(off, i, j + 1)) {
        pairs.push_back({i, j + 1});
      }

      for (const auto& pair : pairs) {
        std::string between_name;
        try {
          between_name = index_name(local, pair);
        } catch (const std::string& err_msg) {
          throw ERROR_MSG("Failed to call index_name(). Error:\n\t[", err_msg,
                          "]");
        }
        // If no error is caught, between_name will be initialized.
        bundle_between(grid_of_tensors_2D[i][j],
                       grid_of_tensors_2D[pair[0]][pair[1]], between_name,
                       scratch);
      }
    }
  }

  // Reorder.
  for (std::size_t i = 0, I = std::size(grid_of_tensors); i < I; ++i) {
    for (std::size_t j = 0, J = std::size(grid_of_tensors[i]); j < J; ++j) {
      if (find_grid_coord_in_list(off, i, j)) {
        continue;
      }

      std::vector<std::string> ordered_indices_2D;

      // Positions of connected active qubits.
      std::vector<std::vector<std::size_t>> pairs;

      if (i > 0 && !find_grid_coord_in_list(off, i - 1, j)) {
        pairs.push_back({i - 1, j});
      }
      if (j > 0 && !find_grid_coord_in_list(off, i, j - 1)) {
        pairs.push_back({i, j - 1});
      }
      if (i < I - 1 && !find_grid_coord_in_list(off, i + 1, j)) {
        pairs.push_back({i + 1, j});
      }
      if (j < J - 1 && !find_grid_coord_in_list(off, i, j + 1)) {
        pairs.push_back({i, j + 1});
      }

      // If this qubit is in the final region, bundling must be adjusted.
      // std::size_t fr_buffer = 0;
      if (find_grid_coord_in_list(final_qubit_region, i, j)) {
        try {
          ordered_indices_2D.push_back(index_name({i, j}, {}));
        } catch (const std::string& err_msg) {
          throw ERROR_MSG("Failed to call index_name(). Error:\n\t[", err_msg,
                          "]");
        }
      }

      std::vector<std::size_t> local = {i, j};
      auto order_fn = order_func(ordering, local);
      std::sort(pairs.begin(), pairs.end(), order_fn);

      for (const auto& pair : pairs) {
        std::vector<std::size_t> q1, q2;
        if (pair[0] < i || pair[1] < j) {
          q1 = pair;
          q2 = local;
        } else {
          q1 = local;
          q2 = pair;
        }
        try {
          ordered_indices_2D.push_back(
              index_name({q1[0], q1[1]}, {q2[0], q2[1]}));
        } catch (const std::string& err_msg) {
          throw ERROR_MSG("Failed to call index_name(). Error:\n\t[", err_msg,
                          "]");
        }
      }

      // Reorder.
      try {
        grid_of_tensors_2D[i][j].reorder(ordered_indices_2D, scratch);
      } catch (const std::string& err_msg) {
        throw ERROR_MSG("Failed to call reorder(). Error:\n\t[", err_msg, "]");
      }
    }
  }
}

// This function is currently not being called.
// TODO: Decide whether or not to deprecate function, also needs to be tested.
void read_wave_function_evolution(
    std::string filename, std::size_t I, std::vector<Tensor>& gates,
    std::vector<std::vector<std::string>>& inputs,
    std::vector<std::vector<std::string>>& outputs, s_type* scratch) {
  if (scratch == nullptr) {
    throw ERROR_MSG("Scratch must be non-null.");
  }
  // Open file.
  auto io = std::ifstream(filename);
  if (io.bad()) {
    throw ERROR_MSG("Cannot open file: ", filename);
  }

  // Gotten from the file.
  std::size_t num_qubits, cycle, q1;
  std::optional<std::size_t> q2;
  std::string gate;
  // Useful for plugging into the tensor network:
  std::vector<std::size_t> i_j_1, i_j_2;

  // The first element should be the number of qubits
  io >> num_qubits;

  // Check for the number of qubits.
  if (num_qubits != I) {
    throw ERROR_MSG(
        "I: ", I, " must be equal to the number of qubits: ", num_qubits, ".");
  }

  std::string line;
  while (getline(io, line))
    if (line.size() && line[0] != '#') {  // Avoid comments
      std::stringstream ss(line);
      // The first element is the cycle
      ss >> cycle;
      // The second element is the gate
      ss >> gate;
      // Get the first position
      ss >> q1;
      // Get the second position in the case
      if (gate == "cz" || gate == "cx" || gate.rfind("fsim", 0) == 0) ss >> *q2;

      // Fill in one-qubit gates.
      if (q2.has_value()) {
        std::string input_index = std::to_string(q1) + ",i";
        std::string output_index = std::to_string(q1) + ",o";
        try {
          gates.push_back(Tensor({input_index, output_index}, {DIM, DIM},
                                 gate_array(gate, {})));
        } catch (const std::string& err_msg) {
          throw ERROR_MSG("Failed to call Tensor(). Error:\n\t[", err_msg, "]");
        }
        inputs.push_back({input_index});
        outputs.push_back({output_index});
      } else {
        std::string input_index1 = std::to_string(q1) + ",i";
        std::string output_index1 = std::to_string(q1) + ",o";
        std::string input_index2 = std::to_string(q2.value()) + ",i";
        std::string output_index2 = std::to_string(q2.value()) + ",o";
        inputs.push_back({input_index1, input_index2});
        outputs.push_back({output_index1, output_index2});
        try {
          gates.push_back(
              Tensor({input_index1, input_index2, output_index1, output_index2},
                     {DIM, DIM, DIM, DIM}, gate_array(gate, {})));
        } catch (const std::string& err_msg) {
          throw ERROR_MSG("Failed to call Tensor(). Error:\n\t[", err_msg, "]");
        }
      }
    }

  // Be proper about pointers.
  scratch = NULL;
}

}  // namespace qflex
