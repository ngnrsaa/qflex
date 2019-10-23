/**
 * @file read_circuit.cpp
 * Helper functions to read quantum circuits from a file.
 * @see https://github.com/benjaminvillalonga/optimized_parallel_QC_with_TN
 *
 * @author Benjamin Villalonga (main contributor), Bron Nelson, Sergio Boixo and
 * Salvatore Mandra
 * @date Created: September 2018
 * @date Modified: February 2019
 *
 * @copyright: Copyright © 2019, United States Government, as represented
 * by the Administrator of the National Aeronautics and Space Administration.
 * All rights reserved.
 * @licence: Apache License, Version 2.0
 */

#include <memory>
#include "read_circuit.h"

namespace qflex {

std::size_t compute_depth(std::istream&& istream) {
  auto is_number = [](const std::string& token) {
    try {
      std::stol(token);
    } catch (...) {
      return false;
    }
    return true;
  };

  auto is_integer = [&is_number](const std::string& token) {
    return is_number(token) and std::stol(token) == std::stod(token);
  };

  auto strip_line = [](std::string line) {
    // Remove everything after '#'
    line = std::regex_replace(line, std::regex("#.*"), "");

    // Remove any special character
    line = std::regex_replace(line, std::regex("[^)(\\s\\ta-zA-Z0-9_.,-]"), "");

    // Convert tabs to spaces
    line = std::regex_replace(line, std::regex("[\\t]"), " ");

    // Remove multiple spaces
    line = std::regex_replace(line, std::regex("[\\s]{2,}"), " ");

    // Remove last space
    line = std::regex_replace(line, std::regex("\\s+$"), "");

    // Remove any space before '('
    line = std::regex_replace(line, std::regex("[\\s]+[(]"), "(");

    // Remove spaces between parentheses
    line = std::regex_replace(line, std::regex("\\s+(?=[^()]*\\))"), "");

    return line;
  };

  auto tokenize = [](const std::string& line,
                     const std::string& regex_expr = "[^\\s]+") {
    std::vector<std::string> tokens;
    auto word_regex = std::regex(regex_expr);
    for (auto w =
             std::sregex_iterator(std::begin(line), std::end(line), word_regex);
         w != std::sregex_iterator(); ++w)
      tokens.push_back(w->str());
    return tokens;
  };

  std::size_t line_counter{0}, last_cycle_number{0};
  std::string line;

  auto error_msg = [&line, &line_counter](const std::string& msg) {
    std::string err_msg =
        "[" + std::to_string(line_counter + 1) + ": " + line + "] " + msg;
    std::cerr << err_msg << std::endl;
    return err_msg;
  };

  std::size_t depth{0};
  std::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t> layers;

  while (std::getline(istream, line)) {
    if (std::size(line = strip_line(line))) {
      // Check that there are only one '(' and one ')'
      if (std::count(std::begin(line), std::end(line), '(') > 1 or
          std::count(std::begin(line), std::end(line), ')') > 1)
        throw error_msg("Wrong format.");

      // Tokenize the line
      auto tokens = tokenize(line);

      // Enforce first line to be a single number which correspond to the number
      // of qubits
      if (line_counter == 0) {
        if (std::size(tokens) != 1 or std::stol(tokens[0]) <= 0)
          throw error_msg(
              "First line in circuit must be the number of active qubits.");
      } else {
        // Check the correct number of tokens
        if (std::size(tokens) < 3)
          throw error_msg(
              "Gate must be specified as: cycle gate_name[(p1[,p2,...])] q1 "
              "[q2, ...]");

        // Check the first token is actually a number
        if (not is_integer(tokens[0]))
          throw error_msg("First token must be a valid cycle number.");

        std::size_t cycle = std::stol(tokens[0]);

        // Check that cycle number is monotonically increasing
        if (cycle < last_cycle_number)
          throw error_msg("Cycle number can only increase.");

        // Add all the qubits
        if (std::size_t num_qubits = std::size(tokens) - 2; num_qubits == 2) {
          std::size_t q1 = std::stol(tokens[2]);
          std::size_t q2 = std::stol(tokens[3]);
          if (q1 > q2) std::swap(q1, q2);

          if (auto new_depth = ++layers[{q1, q2}]; new_depth > depth)
            depth = new_depth;
          else if (num_qubits > 2)
            throw error_msg("k-qubit gates are not supported.");
        }
      }
    }

    // Increment line counter
    ++line_counter;
  }

  return depth;
}

const std::unordered_map<std::string, std::vector<s_type>> _GATES_DATA(
    {// Deltas.
     {"delta_0", std::vector<s_type>({1.0, 0.0})},
     {"delta_1", std::vector<s_type>({0.0, 1.0})},
     // For one-qubit gates, the first index is input an second is output.
     {"h", std::vector<s_type>({_INV_SQRT_2, _INV_SQRT_2, _INV_SQRT_2,
                                -_INV_SQRT_2})},
     {"hz_1_2",
      std::vector<s_type>(
          {{0.5, 0.5}, {_INV_SQRT_2, 0}, {0., -_INV_SQRT_2}, {0.5, 0.5}})},
     {"t", std::vector<s_type>({1.0, 0., 0., {_INV_SQRT_2, _INV_SQRT_2}})},
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

std::vector<s_type> gate_array(const std::string& gate_name) {
  static const std::regex rz_regex("rz\\((.*)\\)");
  std::smatch match;
  if (std::regex_match(gate_name, match, rz_regex) && match.size() > 1) {
    const double angle_rads = _PI * stod(match.str(1));
    const s_type phase = exp(s_type(0., angle_rads / 2));
    return std::vector<s_type>({conj(phase), 0., 0., phase});
  }

  bool valid_gate = _GATES_DATA.find(gate_name) != _GATES_DATA.end();
  if (!valid_gate) {
    std::cout << "Invalid gate name provided: " << gate_name << std::endl;
    assert(false);
  }
  return _GATES_DATA.at(gate_name);
}

/**
 * Returns vector of vectors of s_type with the Schmidt decomposition of the
 * fSim gate.
 * @param theta s_type::value_type with the angle $theta$. Modulo $2\pi$ is
 * taken.
 * @param phi s_type::value_type with the angle $phi$. Modulo $2\pi$ is taken.
 * @param scratch pointer to s_type array with scratch space for all operations
 * performed in this function.
 * @return vector<vector<s_type>> with three elements: first the vector with
 * entries of the first qubit, second the vector with entries of the second
 * qubit, and third the vector with the singular values (informational only).
 */
std::vector<std::vector<s_type>> fSim(s_type::value_type theta,
                                      s_type::value_type phi, s_type* scratch) {
  if (scratch == nullptr) {
    std::cout << "Scratch must be non-null." << std::endl;
    assert(scratch != nullptr);
  }

  static_assert(std::is_floating_point<typename s_type::value_type>::value);

  std::vector<s_type> coeffs(
      {{0.0, -std::sin(theta) / 2},
       {0.0, -std::sin(theta) / 2},
       {(std::cos(-theta / 2) - std::cos(theta)) / 2, std::sin(-theta / 2) / 2},
       {(std::cos(-theta / 2) + std::cos(theta)) / 2,
        std::sin(-theta / 2) / 2}});

  std::vector<double> norm_coeffs(coeffs.size());
  for (int i = 0; i < coeffs.size(); ++i) {
    norm_coeffs[i] = abs(coeffs[i]);
  }

  std::vector<std::vector<s_type>> q1_matrices(
      {{{0., 0.}, {1., 0.}, {1., 0.}, {0., 0.}},
       {{0., 0.}, {0., 1.}, {0., -1.}, {0., 0.}},
       {{std::cos(phi / 4), std::sin(phi / 4)},
        {0., 0.},
        {0., 0.},
        {-std::cos(-phi / 4), -std::sin(-phi / 4)}},
       {{std::cos(phi / 4), std::sin(phi / 4)},
        {0., 0.},
        {0., 0.},
        {std::cos(-phi / 4), std::sin(-phi / 4)}}});
  std::vector<std::vector<s_type>> q2_matrices(
      {{{0., 0.}, {1., 0.}, {1., 0.}, {0., 0.}},
       {{0., 0.}, {0., 1.}, {0., -1.}, {0., 0.}},
       {{std::cos(phi / 4), std::sin(phi / 4)},
        {0., 0.},
        {0., 0.},
        {-std::cos(-phi / 4), -std::sin(-phi / 4)}},
       {{std::cos(phi / 4), std::sin(phi / 4)},
        {0., 0.},
        {0., 0.},
        {std::cos(-phi / 4), std::sin(-phi / 4)}}});

  for (int i = 0; i < coeffs.size(); ++i) {
    for (int j = 0; j < q1_matrices[i].size(); ++j) {
      q1_matrices[i][j] *= coeffs[i];
      q2_matrices[i][j] *= coeffs[i];
    }
  }

  struct cnmm {
    s_type c;
    double n;
    std::vector<s_type> m1;
    std::vector<s_type> m2;
    cnmm(s_type c_, double n_, std::vector<s_type> m1_, std::vector<s_type> m2_)
        : c(c_), n(n_), m1(m1_), m2(m2_) {}
    bool operator<(const cnmm& other) const { return n < other.n; }
  };

  std::vector<cnmm> my_cnmm;
  for (int i = 0; i < coeffs.size(); ++i) {
    my_cnmm.emplace_back(coeffs[i], norm_coeffs[i], q1_matrices[i],
                         q2_matrices[i]);
  }

  sort(my_cnmm.begin(), my_cnmm.end());
  reverse(my_cnmm.begin(), my_cnmm.end());

  std::vector<s_type> vec_q1_tensor, vec_q2_tensor;
  for (auto v : my_cnmm) {
    for (auto w : v.m1) vec_q1_tensor.emplace_back(w);
    for (auto w : v.m2) vec_q2_tensor.emplace_back(w);
  }

  Tensor q1_tensor({"v", "q1i", "q2i"}, {4, 2, 2}, vec_q1_tensor);
  Tensor q2_tensor({"v", "q1i", "q2i"}, {4, 2, 2}, vec_q2_tensor);
  q1_tensor.reorder({"q1i", "v", "q2i"}, scratch);
  q2_tensor.reorder({"q1i", "v", "q2i"}, scratch);
  std::vector<s_type> q1_reordered_tensor(q1_tensor.size());
  std::vector<s_type> q2_reordered_tensor(q2_tensor.size());
  for (int i = 0; i < q1_tensor.size(); ++i) {
    q1_reordered_tensor[i] = *(q1_tensor.data() + i);
    q2_reordered_tensor[i] = *(q2_tensor.data() + i);
  }

  sort(norm_coeffs.begin(), norm_coeffs.end());
  reverse(norm_coeffs.begin(), norm_coeffs.end());

  // TODO: Convert norm_coeffs to vector<s_type> and return it.
  std::vector<std::vector<s_type>> ret_val(
      {q1_reordered_tensor, q2_reordered_tensor});

  return ret_val;
}

// TODO: Refactor this so that all gate arrays are handled similarly regardless
// of how many qubits they operate on and whether they're parametrized. Put gate
// array handling into its own class.
std::tuple<std::vector<s_type>, std::vector<s_type>, std::vector<size_t>>
gate_arrays(const std::string& gate_name, s_type* scratch) {
  if (scratch == nullptr) {
    std::cout << "Scratch must be non-null." << std::endl;
    assert(scratch != nullptr);
  }
  static const std::regex fsim_regex("fsim\\((.*),(.*)\\)");
  std::smatch match;
  if (gate_name == "cz") {
    return std::tuple<std::vector<s_type>, std::vector<s_type>,
                      std::vector<size_t>>(gate_array("cz_q1"),
                                           gate_array("cz_q2"), {2, 2, 2});
  } else if (gate_name == "cx") {
    return std::tuple<std::vector<s_type>, std::vector<s_type>,
                      std::vector<size_t>>(gate_array("cx_q1"),
                                           gate_array("cx_q2"), {2, 2, 2});
  } else if (std::regex_match(gate_name, match, fsim_regex) &&
             match.size() > 2) {
    const double theta_rads = _PI * stod(match.str(1));
    const double phi_rads = _PI * stod(match.str(2));
    std::vector<std::vector<s_type>> ret_val =
        fSim(theta_rads, phi_rads, scratch);
    return std::tuple<std::vector<s_type>, std::vector<s_type>,
                      std::vector<size_t>>(ret_val[0], ret_val[1], {4, 2, 2});
  }
  std::cout << "Invalid gate name provided: " << gate_name << std::endl;
  assert(false);
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
 * @param local vector<int> of the target qubit coordinates.
 * @return function for use as a comparator in std::sort.
 */
std::function<bool(std::vector<int>, std::vector<int>)> order_func(
    const std::list<ContractionOperation>& ordering, std::vector<int> local) {
  return [&ordering, local](const std::vector<int> lhs,
                            const std::vector<int> rhs) {
    // Cuts are projected early and should be ordered first.
    std::vector<std::vector<int>> lhs_pair, rhs_pair;
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
    int lpos = -1;
    int rpos = -1;
    int op_num = 0;
    for (const auto& op : ordering) {
      if (op.op_type != ContractionOperation::EXPAND) continue;
      if (lhs == op.expand.tensor) {
        lpatch = op.expand.id;
        lpos = op_num;
        break;
      }
      op_num++;
    }
    if (lpos == -1) {
      char error[200];
      snprintf(error, sizeof(error),
               "Left hand side of pair not found: (%d,%d),(%d,%d)", local[0],
               local[1], lhs[0], lhs[1]);
      std::cout << error << std::endl;
      std::cout << "Halting reordering." << std::endl;
      assert(false);
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
    if (rpos == -1) {
      char error[200];
      snprintf(error, sizeof(error),
               "Right hand side of pair not found: (%d,%d),(%d,%d)", local[0],
               local[1], rhs[0], rhs[1]);
      std::cout << error << std::endl;
      std::cout << "Halting reordering." << std::endl;
      assert(false);
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
    char error[200];
    snprintf(error, sizeof(error),
             "Failed to compare (%d,%d) and (%d,%d) for local (%d,%d).", lhs[0],
             lhs[1], rhs[0], rhs[1], local[0], local[1]);
    std::cout << error << std::endl;
    std::cout << "Halting reordering." << std::endl;
    assert(false);
  };
}

// This function is currently not being called.
// TODO: Decide whether or not to deprecate function, also needs to be tested.
// TODO(martinop): remove "final_qubit_region" argument?
// This can be derived from the 'x' states in final_conf.
void circuit_data_to_grid_of_tensors(
    std::istream* circuit_data, int I, int J, int K,
    const std::string initial_conf, const std::string final_conf,
    const std::optional<std::vector<std::vector<int>>>& final_qubit_region,
    const std::optional<std::vector<std::vector<int>>>& off,
    std::vector<std::vector<std::vector<Tensor>>>& grid_of_tensors,
    s_type* scratch) {
  if (circuit_data == nullptr) {
    std::cout << "Circuit data stream must be non-null." << std::endl;
    assert(circuit_data != nullptr);
  }
  if (scratch == nullptr) {
    std::cout << "Scratch must be non-null." << std::endl;
    assert(scratch != nullptr);
  }
  // Gotten from the file.
  int circuit_data_num_qubits, cycle, q1, q2;
  std::string gate;
  // Useful for plugging into the tensor network:
  std::vector<int> i_j_1, i_j_2;
  int super_cycle;
  // Calculated from input.
  const int grid_size = I * J;
  const int off_size = off.has_value() ? off.value().size() : 0;
  const int num_active_qubits_from_grid = grid_size - off_size;

  // The first element is required to be the number of active qubits.
  *(circuit_data) >> circuit_data_num_qubits;
  if (circuit_data_num_qubits == 0) {
    std::cout
        << "First line in circuit file must be the number of active qubits."
        << std::endl;
    assert(circuit_data_num_qubits != 0);
  }
  if (circuit_data_num_qubits != num_active_qubits_from_grid) {
    std::cout
        << "The number of active qubits read from the file: "
        << circuit_data_num_qubits
        << ", does not match the number of active qubits read from the grid: "
        << num_active_qubits_from_grid << "." << std::endl;
    assert(circuit_data_num_qubits == num_active_qubits_from_grid);
  }

  // Assert for the length of initial_conf and final_conf.
  {
    // size_t off_size = off.has_value() ? off.value().size() : 0;
    if (initial_conf.size() != num_active_qubits_from_grid) {
      std::cout << "Size of initial_conf: " << initial_conf.size()
                << ", must be equal to the number of qubits: "
                << num_active_qubits_from_grid << "." << std::endl;
      assert(initial_conf.size() == num_active_qubits_from_grid);
    }
    if (final_conf.size() != initial_conf.size()) {
      std::cout << "Size of final_conf: " << final_conf.size()
                << ", must be equal to size of initial_conf: "
                << initial_conf.size() << "." << std::endl;
      assert(final_conf.size() == initial_conf.size());
    }
  }

  // Creating grid variables.
  std::vector<std::vector<std::vector<std::vector<Tensor>>>>
      grid_of_groups_of_tensors(I);
  grid_of_tensors = std::vector<std::vector<std::vector<Tensor>>>(I);
  std::vector<std::vector<std::vector<int>>> counter_group(I);
  for (int i = 0; i < I; ++i) {
    grid_of_groups_of_tensors[i] =
        std::vector<std::vector<std::vector<Tensor>>>(J);
    grid_of_tensors[i] = std::vector<std::vector<Tensor>>(J);
    counter_group[i] = std::vector<std::vector<int>>(J);
    for (int j = 0; j < J; ++j) {
      grid_of_groups_of_tensors[i][j] = std::vector<std::vector<Tensor>>(K);
      grid_of_tensors[i][j] = std::vector<Tensor>(K);
      counter_group[i][j] = std::vector<int>(K, 0);
      for (int k = 0; k < K; ++k) {
        grid_of_groups_of_tensors[i][j][k] = std::vector<Tensor>();
      }
    }
  }

  // Insert deltas and Hadamards to first layer.
  int idx = 0;
  for (int q = 0; q < grid_size; ++q) {
    std::vector<int> i_j = get_qubit_coords(q, J);
    int i = i_j[0], j = i_j[1];
    if (find_grid_coord_in_list(off, i, j)) {
      continue;
    }
    std::string delta_gate = (initial_conf[idx] == '0') ? "delta_0" : "delta_1";
    grid_of_groups_of_tensors[i][j][0].push_back(
        Tensor({"th"}, {2}, gate_array(delta_gate)));
    grid_of_groups_of_tensors[i][j][0].push_back(
        Tensor({"th", "t0"}, {2, 2}, gate_array("h")));
    idx += 1;
  }

  std::string line;
  int line_counter = 0;
  int cycle_holder = 0;
  std::unordered_set<int> used_qubits;
  // Read one line at a time from the circuit, skipping comments.
  while ((++line_counter, getline(*circuit_data, line))) {
    if (line.size() && line[0] != '#') {
      std::stringstream ss(line);
      // The first element is the cycle
      ss >> cycle;
      if (cycle != cycle_holder) {
        cycle_holder = cycle;
        used_qubits.clear();
      }
      // The second element is the gate
      ss >> gate;
      // Get the first position
      // TODO: Stop relying on the assumption that the gate and its parameters
      // can be read as one token without spaces. This is (mostly) fine for
      // "rz(0.5)", but will fail for, e.g., "fsim(0.25, -0.5)".
      ss >> q1;
      // Get the second position if needed
      // TODO: Two-qubit gates should be encapsulated better.
      if (gate == "cz" || gate == "cx" || gate.rfind("fsim", 0) == 0) {
        ss >> q2;
      } else {
        q2 = -1;
      }

      // Check that q1 hasn't already been used in this cycle.
      std::unordered_set<int>::const_iterator q1_used = used_qubits.find(q1);
      if (q1_used != used_qubits.end()) {
        std::cout << "The qubit " << q1 << " in '" << line_counter << ": "
                  << line << "' has already been used in this cycle."
                  << std::endl;
        assert(q1_used == used_qubits.end());
      } else {
        used_qubits.insert(q1);
      }
      // Check that q2 hasn't already been used in this cycle when applicable.
      if (q2 != -1) {
        std::unordered_set<int>::const_iterator q2_used = used_qubits.find(q2);
        if (q2_used != used_qubits.end()) {
          std::cout << "The qubit " << q2 << " in '" << line_counter << ": "
                    << line << "' has already been used in this cycle. "
                    << std::endl;
          assert(q2_used == used_qubits.end());
        } else {
          used_qubits.insert(q2);
        }
      }

      // Get i, j and super_cycle
      i_j_1 = get_qubit_coords(q1, J);
      if (q2 >= 0) {
        i_j_2 = get_qubit_coords(q2, J);
      }
      super_cycle = (cycle - 1) / SUPER_CYCLE_DEPTH;

      // Fill in one-qubit gates.
      if (q2 < 0 && cycle > 0 && cycle <= SUPER_CYCLE_DEPTH * K) {
        // Check that position is an active qubit
        bool qubit_off = find_grid_coord_in_list(off, i_j_1[0], i_j_1[1]);
        if (qubit_off) {
          std::cout << "The qubit in '" << line << "' references (" << i_j_1[0]
                    << ", " << i_j_1[1]
                    << ") which must be coordinates of an active qubit."
                    << std::endl;
          assert(!qubit_off);
        }
        std::string input_index =
            "t" +
            std::to_string(counter_group[i_j_1[0]][i_j_1[1]][super_cycle]);
        std::string output_index =
            "t" +
            std::to_string(counter_group[i_j_1[0]][i_j_1[1]][super_cycle] + 1);
        ++counter_group[i_j_1[0]][i_j_1[1]][super_cycle];
        grid_of_groups_of_tensors[i_j_1[0]][i_j_1[1]][super_cycle].push_back(
            Tensor({input_index, output_index}, {2, 2}, gate_array(gate)));
      }
      // Fill in two-qubit gates.
      if (q2 >= 0 && cycle > 0 && cycle <= SUPER_CYCLE_DEPTH * K) {
        // Check that positions are active qubits
        bool first_qubit_off = find_grid_coord_in_list(off, i_j_1[0], i_j_1[1]);
        bool second_qubit_off =
            find_grid_coord_in_list(off, i_j_2[0], i_j_2[1]);
        if (first_qubit_off) {
          std::cout << "The first qubit of '" << line << "' references ("
                    << i_j_1[0] << ", " << i_j_1[1]
                    << ") which must be coordinates of an active qubit."
                    << std::endl;
          assert(!first_qubit_off);
        }
        if (second_qubit_off) {
          std::cout << "The second qubit of '" << line << "' references ("
                    << i_j_2[0] << ", " << i_j_2[1]
                    << ") which must be coordinates of an active qubit."
                    << std::endl;
          assert(!second_qubit_off);
        }
        std::vector<s_type> gate_q1;
        std::vector<s_type> gate_q2;
        std::vector<size_t> dimensions;
        tie(gate_q1, gate_q2, dimensions) = gate_arrays(gate, scratch);
        std::string input_index_1 =
            "t" +
            std::to_string(counter_group[i_j_1[0]][i_j_1[1]][super_cycle]);
        std::string output_index_1 =
            "t" +
            std::to_string(counter_group[i_j_1[0]][i_j_1[1]][super_cycle] + 1);
        std::string virtual_index =
            index_name({i_j_1[0], i_j_1[1], super_cycle},
                       {i_j_2[0], i_j_2[1], super_cycle});
        std::string input_index_2 =
            "t" +
            std::to_string(counter_group[i_j_2[0]][i_j_2[1]][super_cycle]);
        std::string output_index_2 =
            "t" +
            std::to_string(counter_group[i_j_2[0]][i_j_2[1]][super_cycle] + 1);
        ++counter_group[i_j_1[0]][i_j_1[1]][super_cycle];
        ++counter_group[i_j_2[0]][i_j_2[1]][super_cycle];
        grid_of_groups_of_tensors[i_j_1[0]][i_j_1[1]][super_cycle].push_back(
            Tensor({input_index_1, virtual_index, output_index_1}, dimensions,
                   gate_q1));
        grid_of_groups_of_tensors[i_j_2[0]][i_j_2[1]][super_cycle].push_back(
            Tensor({input_index_2, virtual_index, output_index_2}, dimensions,
                   gate_q2));
      }
    }
  }
  // Insert Hadamards and deltas to last layer.
  idx = -1;
  for (int q = 0; q < grid_size; ++q) {
    std::vector<int> i_j = get_qubit_coords(q, J);
    int i = i_j[0], j = i_j[1];
    int k = K - 1;
    if (find_grid_coord_in_list(off, i, j)) {
      continue;
    }
    idx += 1;
    std::string last_index = "t" + std::to_string(counter_group[i][j][k]);
    grid_of_groups_of_tensors[i][j][k].push_back(
        Tensor({"th", last_index}, {2, 2}, gate_array("h")));
    if (find_grid_coord_in_list(final_qubit_region, i, j)) {
      continue;
    }
    std::string delta_gate = (final_conf[idx] == '0') ? "delta_0" : "delta_1";
    grid_of_groups_of_tensors[i][j][k].push_back(
        Tensor({"th"}, {2}, gate_array(delta_gate)));
  }

  // Contracting each group of gates into a single tensor.
  for (int i = 0; i < I; ++i)
    for (int j = 0; j < J; ++j)
      for (int k = 0; k < K; ++k) {
        if (find_grid_coord_in_list(off, i, j)) {
          continue;
        }

        std::vector<Tensor>& group = grid_of_groups_of_tensors[i][j][k];
        std::vector<Tensor> group_containers(
            SUPER_CYCLE_DEPTH + 2,  // +2 for d and H.
            Tensor({""}, {(int)pow(DIM, 6)}));
        group_containers[0] = group[0];
        int t = 1;
        for (t = 1; t < group.size(); ++t) {
          multiply(group_containers[t - 1], group[t], group_containers[t],
                   scratch);
        }
        grid_of_tensors[i][j][k] = group_containers[t - 1];
      }

  // Rename "t..." indices.
  for (int i = 0; i < I; ++i)
    for (int j = 0; j < J; ++j)
      for (int k = 0; k < K; ++k) {
        if (find_grid_coord_in_list(off, i, j)) {
          continue;
        }
        if (k > 0) {
          std::string new_first_index = index_name({i, j, k - 1}, {i, j, k});
          grid_of_tensors[i][j][k].rename_index("t0", new_first_index);
        }
        if (k < K - 1) {
          std::string last_index = "t" + std::to_string(counter_group[i][j][k]);
          std::string new_last_index = index_name({i, j, k}, {i, j, k + 1});
          grid_of_tensors[i][j][k].rename_index(last_index, new_last_index);
        }
        if (k == K - 1 && find_grid_coord_in_list(final_qubit_region, i, j)) {
          std::string last_index = "th";
          std::string new_last_index = index_name({i, j}, {});
          grid_of_tensors[i][j][k].rename_index(last_index, new_last_index);
        }
      }

  // Be proper about pointers.
  scratch = NULL;
}

// This function is currently not being called.
// TODO: Decide whether or not to deprecate function, also needs to be tested.
void grid_of_tensors_3D_to_2D( std::vector<std::vector<std::vector<Tensor>>>& grid_of_tensors_3D,
    std::vector<std::vector<Tensor>>& grid_of_tensors_2D,
    std::optional<std::vector<std::vector<int>>> final_qubit_region,
    std::optional<std::vector<std::vector<int>>> off,
    const std::list<ContractionOperation>& ordering, s_type* scratch) {
  if (scratch == nullptr) {
    std::cout << "Scratch must be non-null." << std::endl;
    assert(scratch != nullptr);
  }
  // Get dimensions and super_dim = DIM^k.
  const int I = grid_of_tensors_3D.size();
  const int J = grid_of_tensors_3D[0].size();
  const int K = grid_of_tensors_3D[0][0].size();
  const int super_dim = (int)pow(DIM, K);

  // Contract vertically and fill grid_of_tensors_2D.
  for (int i = 0; i < I; ++i) {
    for (int j = 0; j < J; ++j) {
      if (find_grid_coord_in_list(off, i, j)) {
        continue;
      }

      size_t container_dim = (int)pow(super_dim, 4);
      if (find_grid_coord_in_list(final_qubit_region, i, j)) {
        container_dim *= DIM;
      }
      std::vector<Tensor> group_containers =
          std::vector<Tensor>(2, Tensor({""}, {container_dim}));
      Tensor* source_container = &group_containers[0];
      Tensor* target_container = &group_containers[1];

      if (K == 1) {
        grid_of_tensors_2D[i][j] = grid_of_tensors_3D[i][j][0];
      } else {
        multiply(grid_of_tensors_3D[i][j][0], grid_of_tensors_3D[i][j][1],
                 *source_container, scratch);
        for (int k = 1; k < K - 1; ++k) {
          multiply(*source_container, grid_of_tensors_3D[i][j][k + 1],
                   *target_container, scratch);
          Tensor* swap_container = source_container;
          source_container = target_container;
          target_container = swap_container;
        }
        grid_of_tensors_2D[i][j] = *source_container;
      }
    }
  }

  // Reorder and bundle.
  for (int i = 0; i < I; ++i) {
    for (int j = 0; j < J; ++j) {
      if (find_grid_coord_in_list(off, i, j)) {
        continue;
      }

      std::vector<std::string> ordered_indices_3D;
      std::vector<std::string> indices_2D;

      // Positions of connected active qubits.
      std::vector<std::vector<int>> pairs;

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
      int fr_buffer = 0;
      if (find_grid_coord_in_list(final_qubit_region, i, j)) {
        ordered_indices_3D.push_back(index_name({i, j}, {}));
        fr_buffer = 1;
      }

      std::vector<int> local = {i, j};
      auto order_fn = order_func(ordering, local);
      std::sort(pairs.begin(), pairs.end(), order_fn);

      for (const auto& pair : pairs) {
        std::vector<int> q1, q2;
        if (pair[0] < i || pair[1] < j) {
          q1 = pair;
          q2 = local;
        } else {
          q1 = local;
          q2 = pair;
        }
        for (int k = 0; k < K; ++k) {
          ordered_indices_3D.push_back(
              index_name({q1[0], q1[1], k}, {q2[0], q2[1], k}));
        }
        indices_2D.push_back(index_name(q1, q2));
      }

      // Reorder.
      grid_of_tensors_2D[i][j].reorder(ordered_indices_3D, scratch);

      // Bundle.
      int max_idx;
      max_idx = indices_2D.size();
      for (int idx_num = 0; idx_num < max_idx; ++idx_num) {
        std::vector<std::string> indices_to_bundle(
            ordered_indices_3D.begin() + idx_num * K + fr_buffer,
            ordered_indices_3D.begin() + (idx_num + 1) * K + fr_buffer);
        grid_of_tensors_2D[i][j].bundle(indices_to_bundle, indices_2D[idx_num]);
      }
    }
  }

  // Be proper about pointers.
  scratch = NULL;
}

// TODO: add tests for this function. Compactify code. Use index_name()
// function for all index names used here. Use smart pointers where possible.
// Add depth functionality to read a circuit up to a certain cycle.
void circuit_data_to_tensor_network(
    std::istream* circuit_data, int I, int J,
    const std::string initial_conf, const std::string final_conf,
    const std::optional<std::vector<std::vector<int>>>& final_qubit_region,
    const std::optional<std::vector<std::vector<int>>>& off,
    std::vector<std::vector<std::vector<Tensor>>>& grid_of_tensors,
    s_type* scratch) {
  if (circuit_data == nullptr) {
    std::cout << "Circuit data stream must be non-null." << std::endl;
    assert(circuit_data != nullptr);
  }
  if (scratch == nullptr) {
    std::cout << "Scratch must be non-null." << std::endl;
    assert(scratch != nullptr);
  }
  // Gotten from the file.
  int circuit_data_num_qubits, cycle, q1, q2;
  std::string gate;
  // Useful for plugging into the tensor network:
  std::vector<int> i_j_1, i_j_2;
  // Calculated from input.
  const int grid_size = I * J;
  const int off_size = off.has_value() ? off.value().size() : 0;
  const int num_active_qubits_from_grid = grid_size - off_size;

  // The first element is required to be the number of active qubits.
  *(circuit_data) >> circuit_data_num_qubits;
  if (circuit_data_num_qubits == 0) {
    std::cout
        << "First line in circuit file must be the number of active qubits."
        << std::endl;
    assert(circuit_data_num_qubits != 0);
  }
  if (circuit_data_num_qubits != num_active_qubits_from_grid) {
    std::cout
        << "The number of active qubits read from the file: "
        << circuit_data_num_qubits
        << ", does not match the number of active qubits read from the grid: "
        << num_active_qubits_from_grid << "." << std::endl;
    assert(circuit_data_num_qubits == num_active_qubits_from_grid);
  }

  // Assert for the length of initial_conf and final_conf.
  {
    // size_t off_size = off.has_value() ? off.value().size() : 0;
    if (initial_conf.size() != num_active_qubits_from_grid) {
      std::cout << "Size of initial_conf: " << initial_conf.size()
                << ", must be equal to the number of qubits: "
                << num_active_qubits_from_grid << "." << std::endl;
      assert(initial_conf.size() == num_active_qubits_from_grid);
    }
    if (final_conf.size() != initial_conf.size()) {
      std::cout << "Size of final_conf: " << final_conf.size()
                << ", must be equal to size of initial_conf: "
                << initial_conf.size() << "." << std::endl;
      assert(final_conf.size() == initial_conf.size());
    }
  }

  // Creating grid variables.
  grid_of_tensors = std::vector<std::vector<std::vector<Tensor>>>(I);
  std::vector<std::vector<int>> grid_of_counters(I);
  std::unordered_map<std::string, int> link_counters;
  for (int i = 0; i < I; ++i) {
    grid_of_tensors[i] = std::vector<std::vector<Tensor>>(J);
    grid_of_counters[i] = std::vector<int>(J);
    for (int j = 0; j < J; ++j) {
      grid_of_tensors[i][j] = std::vector<Tensor>();
      grid_of_counters[i][j] = 0;
    }
  }

  // Insert deltas to first layer.
  int idx = 0;
  for (int q = 0; q < grid_size; ++q) {
    std::vector<int> i_j = get_qubit_coords(q, J);
    int i = i_j[0], j = i_j[1];
    if (find_grid_coord_in_list(off, i, j)) {
      continue;
    }
    std::string delta_gate = (initial_conf[idx] == '0') ? "delta_0" : "delta_1";
    std::string output_name = "(" + std::to_string(i_j[0]) + ","
        + std::to_string(i_j[1]) + "),("
        + std::to_string(grid_of_counters[i][j]) + ")";
    grid_of_tensors[i][j].push_back(
        Tensor({output_name}, {2}, gate_array(delta_gate)));
    idx += 1;
  }

  std::string line;
  int line_counter = 0;
  int cycle_holder = 0;
  std::unordered_set<int> used_qubits;
  // Read one line at a time from the circuit, skipping comments.
  while (getline(*circuit_data, line)) {
    if (line.size() && line[0] != '#') {
      std::stringstream ss(line);
      // The first element is the cycle
      ss >> cycle;
      if (cycle != cycle_holder) {
        cycle_holder = cycle;
        used_qubits.clear();
      }
      // The second element is the gate
      ss >> gate;
      // Get the first position
      // TODO: Stop relying on the assumption that the gate and its parameters
      // can be read as one token without spaces. This is (mostly) fine for
      // "rz(0.5)", but will fail for, e.g., "fsim(0.25, -0.5)".
      ss >> q1;
      // Get the second position in the case
      // TODO: Two-qubit gates should be encapsulated better.
      if (gate == "cz" || gate == "cx" || gate.rfind("fsim", 0) == 0) {
        ss >> q2;
      } else {
        q2 = -1;
      }


      // Check that q1 hasn't already been used in this cycle.
      std::unordered_set<int>::const_iterator q1_used = used_qubits.find(q1);
      if (q1_used != used_qubits.end()) {
        std::cout << "The qubit " << q1 << " in '" << line_counter << ": "
                  << line << "' has already been used in this cycle."
                  << std::endl;
        assert(q1_used == used_qubits.end());
      } else {
        used_qubits.insert(q1);
      }
      // Check that q2 hasn't already been used in this cycle when applicable.
      if (q2 != -1) {
        std::unordered_set<int>::const_iterator q2_used = used_qubits.find(q2);
        if (q2_used != used_qubits.end()) {
          std::cout << "The qubit " << q2 << " in '" << line_counter << ": "
                    << line << "' has already been used in this cycle. "
                    << std::endl;
          assert(q2_used == used_qubits.end());
        } else {
          used_qubits.insert(q2);
        }
      }

      // Get i and j
      i_j_1 = get_qubit_coords(q1, J);
      if (q2 >= 0) {
        i_j_2 = get_qubit_coords(q2, J);
      }

      // Fill in one-qubit gates.
      if (q2 < 0) {
        // Check that position is an active qubit
        bool qubit_off = find_grid_coord_in_list(off, i_j_1[0], i_j_1[1]);
        if (qubit_off) {
          std::cout << "The qubit in '" << line << "' references (" << i_j_1[0]
                    << ", " << i_j_1[1]
                    << ") which must be coordinates of an active qubit."
                    << std::endl;
          assert(!qubit_off);
        }
        std::string input_name = "(" + std::to_string(i_j_1[0]) + ","
            + std::to_string(i_j_1[1]) + "),("
            + std::to_string(grid_of_counters[i_j_1[0]][i_j_1[1]]) + ")";
        ++grid_of_counters[i_j_1[0]][i_j_1[1]];
        std::string output_name = "(" + std::to_string(i_j_1[0]) + ","
            + std::to_string(i_j_1[1]) + "),("
            + std::to_string(grid_of_counters[i_j_1[0]][i_j_1[1]]) + ")";
        grid_of_tensors[i_j_1[0]][i_j_1[1]].push_back(
            Tensor({input_name, output_name}, {2, 2}, gate_array(gate)));
      }
      // Fill in two-qubit gates.
      if (q2 >= 0) {
        // Check that positions are active qubits
        bool first_qubit_off = find_grid_coord_in_list(off, i_j_1[0], i_j_1[1]);
        bool second_qubit_off =
            find_grid_coord_in_list(off, i_j_2[0], i_j_2[1]);
        if (first_qubit_off) {
          std::cout << "The first qubit of '" << line << "' references ("
                    << i_j_1[0] << ", " << i_j_1[1]
                    << ") which must be coordinates of an active qubit."
                    << std::endl;
          assert(!first_qubit_off);
        }
        if (second_qubit_off) {
          std::cout << "The second qubit of '" << line << "' references ("
                    << i_j_2[0] << ", " << i_j_2[1]
                    << ") which must be coordinates of an active qubit."
                    << std::endl;
          assert(!second_qubit_off);
        }
        bool nearest_neighbors = (std::abs(i_j_1[0] - i_j_2[0])
            + std::abs(i_j_1[1] - i_j_2[1])) == 1;
        if (!nearest_neighbors) {
          std::cout << "Qubits in '" << line << "' are not nearest neighbors."
                    << std::endl;
          assert(nearest_neighbors);
        }
        std::vector<s_type> gate_q1;
        std::vector<s_type> gate_q2;
        std::vector<size_t> dimensions;
        tie(gate_q1, gate_q2, dimensions) = gate_arrays(gate, scratch);
        std::string link_name = index_name(i_j_1, i_j_2);
        link_counters[link_name]++;
        int counter = link_counters[link_name];
        std::string virtual_name = "("
            + std::to_string(std::min(i_j_1[0], i_j_2[0])) + ","
            + std::to_string(std::min(i_j_1[1], i_j_2[1])) + ","
            + std::to_string(counter) + "),("
            + std::to_string(std::max(i_j_1[0], i_j_2[0])) + ","
            + std::to_string(std::max(i_j_1[1], i_j_2[1])) + ","
            + std::to_string(counter) + ")";
        std::string input_name_1 = "(" + std::to_string(i_j_1[0]) + ","
            + std::to_string(i_j_1[1]) + "),("
            + std::to_string(grid_of_counters[i_j_1[0]][i_j_1[1]]) + ")";
        std::string input_name_2 = "(" + std::to_string(i_j_2[0]) + ","
            + std::to_string(i_j_2[1]) + "),("
            + std::to_string(grid_of_counters[i_j_2[0]][i_j_2[1]]) + ")";
        ++grid_of_counters[i_j_1[0]][i_j_1[1]];
        ++grid_of_counters[i_j_2[0]][i_j_2[1]];
        std::string output_name_1 = "(" + std::to_string(i_j_1[0]) + ","
            + std::to_string(i_j_1[1]) + "),("
            + std::to_string(grid_of_counters[i_j_1[0]][i_j_1[1]]) + ")";
        std::string output_name_2 = "(" + std::to_string(i_j_2[0]) + ","
            + std::to_string(i_j_2[1]) + "),("
            + std::to_string(grid_of_counters[i_j_2[0]][i_j_2[1]]) + ")";
        grid_of_tensors[i_j_1[0]][i_j_1[1]].push_back(
            Tensor({input_name_1, virtual_name, output_name_1}, dimensions,
                   gate_q1));
        grid_of_tensors[i_j_2[0]][i_j_2[1]].push_back(
            Tensor({input_name_2, virtual_name, output_name_2}, dimensions,
                   gate_q2));
      }
    }
  }

  // Insert deltas to last layer on qubits that are in not in
  // final_qubit_region. Rename last index when in final_qubit_region to
  // "(i,j),(o)".
  idx = -1;
  for (int q = 0; q < grid_size; ++q) {
    std::vector<int> i_j = get_qubit_coords(q, J);
    int i = i_j[0], j = i_j[1];
    if (find_grid_coord_in_list(off, i, j)) {
      continue;
    }
    idx += 1;
    std::string last_name = "(" + std::to_string(i_j[0]) + ","
        + std::to_string(i_j[1]) + "),("
        + std::to_string(grid_of_counters[i][j]) + ")";
    if (find_grid_coord_in_list(final_qubit_region, i, j)) {
      std::string output_name = "(" + std::to_string(i_j[0]) + ","
          + std::to_string(i_j[1]) + "),(o)";
      grid_of_tensors[i][j].back().rename_index(last_name, output_name);
    } else {
      std::string delta_gate = (final_conf[idx] == '0') ? "delta_0"
          : "delta_1";
      grid_of_tensors[i][j].push_back(
          Tensor({last_name}, {2}, gate_array(delta_gate)));
    }
  }

  // Be proper about pointers.
  scratch = NULL; 
}

// TODO: add tests for this function. Optimize contraction procedure; it uses
// far more memory than needed currently. Compactify code. Use smart pointers
// where possible.
void flatten_grid_of_tensors(
    std::vector<std::vector<std::vector<Tensor>>>& grid_of_tensors,
    std::vector<std::vector<Tensor>>& grid_of_tensors_2D,
    std::optional<std::vector<std::vector<int>>> final_qubit_region,
    std::optional<std::vector<std::vector<int>>> off,
    const std::list<ContractionOperation>& ordering, s_type* scratch) {

  if (scratch == nullptr) {
    std::cout << "Scratch must be non-null." << std::endl;
    assert(scratch != nullptr);
  }

  // Contract vertically and fill grid_of_tensors_2D.
  int I = grid_of_tensors.size();
  for (int i = 0; i < I; ++i) {
    int J = grid_of_tensors[i].size();
    for (int j = 0; j < J; ++j) {
      if (find_grid_coord_in_list(off, i, j))
        continue;
      int K = grid_of_tensors[i][j].size();
      std::vector<Tensor> column_of_tensors(K);
      column_of_tensors[0] = Tensor(grid_of_tensors[i][j][0]);
      for (int k = 0; k < K-1; ++k) {
        Tensor A(column_of_tensors[k]);
        Tensor B(grid_of_tensors[i][j][k + 1]);
        size_t result_dimension = result_size(A, B);
        Tensor C({""}, {result_dimension});
        multiply(A, B, C, scratch);
        column_of_tensors[k + 1] = Tensor(C);
      }
      grid_of_tensors_2D[i][j] = Tensor(column_of_tensors.back());
    }
  }

  for (int i = 0; i < I; ++i) {
    int J = grid_of_tensors[i].size();
    for (int j = 0; j < J; ++j) {
      if (find_grid_coord_in_list(off, i, j)) {
        continue;
      }

      // Positions of connected active qubits.
      std::vector<int> local({i, j});
      std::vector<std::vector<int>> pairs;

      if (i < I - 1 && !find_grid_coord_in_list(off, i + 1, j)) {
        pairs.push_back({i + 1, j});
      }
      if (j < J - 1 && !find_grid_coord_in_list(off, i, j + 1)) {
        pairs.push_back({i, j + 1});
      }

      for (const auto& pair : pairs) {
        std::string between_name = index_name(local, pair);
        bundle_between(grid_of_tensors_2D[i][j],
            grid_of_tensors_2D[pair[0]][pair[1]], between_name, scratch);
      }
    }
  }

  // Reorder.
  for (int i = 0; i < I; ++i) {
    int J = grid_of_tensors[i].size();
    for (int j = 0; j < J; ++j) {
      if (find_grid_coord_in_list(off, i, j)) {
        continue;
      }

      std::vector<std::string> ordered_indices_2D;

      // Positions of connected active qubits.
      std::vector<std::vector<int>> pairs;

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
      int fr_buffer = 0;
      if (find_grid_coord_in_list(final_qubit_region, i, j)) {
        ordered_indices_2D.push_back(index_name({i, j}, {}));
        fr_buffer = 1;
      }

      std::vector<int> local = {i, j};
      auto order_fn = order_func(ordering, local);
      std::sort(pairs.begin(), pairs.end(), order_fn);

      for (const auto& pair : pairs) {
        std::vector<int> q1, q2;
        if (pair[0] < i || pair[1] < j) {
          q1 = pair;
          q2 = local;
        } else {
          q1 = local;
          q2 = pair;
        }
        ordered_indices_2D.push_back(
            index_name({q1[0], q1[1]}, {q2[0], q2[1]}));
      }

      // Reorder.
      grid_of_tensors_2D[i][j].reorder(ordered_indices_2D, scratch);
    }
  }
}

// This function is currently not being called.
// TODO: Decide whether or not to deprecate function, also needs to be tested.
void read_wave_function_evolution(
    std::string filename, int I, std::vector<Tensor>& gates,
    std::vector<std::vector<std::string>>& inputs,
    std::vector<std::vector<std::string>>& outputs, s_type* scratch) {
  if (scratch == nullptr) {
    std::cout << "Scratch must be non-null." << std::endl;
    assert(scratch != nullptr);
  }
  // Open file.
  auto io = std::ifstream(filename);
  if (io.bad()) {
    std::cout << "Cannot open file: " << filename << std::endl;
    assert(io.good());
  }

  // Gotten from the file.
  int num_qubits, cycle, q1, q2;
  std::string gate;
  // Useful for plugging into the tensor network:
  std::vector<int> i_j_1, i_j_2;
  int super_cycle;

  // The first element should be the number of qubits
  io >> num_qubits;

  // Assert for the number of qubits.
  if (num_qubits != I) {
    std::cout << "I: " << I
              << " must be equal to the number of qubits: " << num_qubits
              << std::endl;
    assert(num_qubits == I);
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
      if (gate == "cz" || gate == "cx" || gate.rfind("fsim", 0) == 0)
        ss >> q2;
      else
        q2 = -1;

      // Fill in one-qubit gates.
      if (q2 < 0) {
        std::string input_index = std::to_string(q1) + ",i";
        std::string output_index = std::to_string(q1) + ",o";
        gates.push_back(
            Tensor({input_index, output_index}, {DIM, DIM}, gate_array(gate)));
        inputs.push_back({input_index});
        outputs.push_back({output_index});
      }
      if (q2 >= 0) {
        std::string input_index1 = std::to_string(q1) + ",i";
        std::string output_index1 = std::to_string(q1) + ",o";
        std::string input_index2 = std::to_string(q2) + ",i";
        std::string output_index2 = std::to_string(q2) + ",o";
        inputs.push_back({input_index1, input_index2});
        outputs.push_back({output_index1, output_index2});
        gates.push_back(
            Tensor({input_index1, input_index2, output_index1, output_index2},
                   {DIM, DIM, DIM, DIM}, gate_array(gate)));
      }
    }

  // Be proper about pointers.
  scratch = NULL;
}

}  // namespace qflex
