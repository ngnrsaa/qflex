/**
 * @file read_circuit.cpp
 * Helper functions to read quantum circuits from a file.
 * @see https://github.com/benjaminvillalonga/optimized_parallel_QC_with_TN
 *
 * @author Benjamin Villalonga
 * @date Created: September 2018
 * @date Modified: February 2019
 */

#include "read_circuit.h"

namespace {

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
                                 0., 0., 0., -1.})}});

std::vector<s_type> gate_array(const std::string& gate_name) {
  static const std::regex rz_regex("rz\\((.*)\\)");
  std::smatch match;
  if (std::regex_match(gate_name, match, rz_regex) && match.size() > 1) {
    const double angle_rads = _PI * stod(match.str(1));
    const s_type phase = exp(s_type(0., angle_rads / 2));
    return std::vector<s_type>({conj(phase), 0., 0., phase});
  }

  // TODO: Better error-checking if the gate name isn't recognised.
  return _GATES_DATA.at(gate_name);
}

/**
 * Returns vector of vectors of s_type with the Schmidt decomposition of the
 * fSim gate.
 * @param theta double with the angle $theta$. Modulo $2\pi$ is taken.
 * @param phi double with the angle $phi$. Modulo $2\pi$ is taken.
 * @param scratch pointer to s_type array with scratch space for all operations
 * performed in this function.
 * @return vector<vector<s_type>> with three elements: first the vector with
 * entries of the first qubit, second the vector with entries of the second
 * qubit, and third the vector with the singular values (informational only).
 */
std::vector<std::vector<s_type>> fSim(double theta, double phi,
                                      s_type* scratch) {
  std::vector<s_type> coeffs(
      {{0.0, -0.5 * sin(theta)},
       {0.0, -0.5 * sin(theta)},
       {0.5 * (cos(-theta / 2.) - cos(theta)), 0.5 * sin(-theta / 2.)},
       {0.5 * (cos(-theta / 2.) + cos(theta)), 0.5 * sin(-theta / 2.)}});

  std::vector<double> norm_coeffs(coeffs.size());
  for (int i = 0; i < coeffs.size(); ++i) {
    norm_coeffs[i] = abs(coeffs[i]);
  }

  std::vector<std::vector<s_type>> q1_matrices(
      {{{0., 0.}, {1., 0.}, {1., 0.}, {0., 0.}},
       {{0., 0.}, {0., 1.}, {0., -1.}, {0., 0.}},
       {{cos(phi / 4.), sin(phi / 4.)},
        {0., 0.},
        {0., 0.},
        {-cos(-phi / 4.), -sin(-phi / 4.)}},
       {{cos(phi / 4.), sin(phi / 4.)},
        {0., 0.},
        {0., 0.},
        {cos(-phi / 4.), sin(-phi / 4.)}}});
  std::vector<std::vector<s_type>> q2_matrices(
      {{{0., 0.}, {1., 0.}, {1., 0.}, {0., 0.}},
       {{0., 0.}, {0., 1.}, {0., -1.}, {0., 0.}},
       {{cos(phi / 4.), sin(phi / 4.)},
        {0., 0.},
        {0., 0.},
        {-cos(-phi / 4.), -sin(-phi / 4.)}},
       {{cos(phi / 4.), sin(phi / 4.)},
        {0., 0.},
        {0., 0.},
        {cos(-phi / 4.), sin(-phi / 4.)}}});

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

  std::vector<s_type> q1_tensor, q2_tensor;
  for (auto v : my_cnmm) {
    for (auto w : v.m1) q1_tensor.emplace_back(w);
    for (auto w : v.m2) q2_tensor.emplace_back(w);
  }

  MKLTensor q1_mkltensor({"v", "q1i", "q2i"}, {4, 2, 2}, q1_tensor);
  MKLTensor q2_mkltensor({"v", "q1i", "q2i"}, {4, 2, 2}, q2_tensor);
  q1_mkltensor.reorder({"q1i", "v", "q2i"}, scratch);
  q2_mkltensor.reorder({"q1i", "v", "q2i"}, scratch);
  std::vector<s_type> q1_reordered_tensor(q1_mkltensor.size());
  std::vector<s_type> q2_reordered_tensor(q2_mkltensor.size());
  for (int i = 0; i < q1_mkltensor.size(); ++i) {
    q1_reordered_tensor[i] = *(q1_mkltensor.data() + i);
    q2_reordered_tensor[i] = *(q2_mkltensor.data() + i);
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
  static const std::regex fsim_regex("fsim\\((.*),(.*)\\)");
  std::smatch match;
  if (gate_name == "cz") {
    return std::tuple<std::vector<s_type>, std::vector<s_type>,
                      std::vector<size_t>>(gate_array("cz_q1"),
                                           gate_array("cz_q2"), {2, 2, 2});
  } else if (std::regex_match(gate_name, match, fsim_regex) &&
             match.size() > 2) {
    const double theta_rads = _PI * stod(match.str(1));
    const double phi_rads = _PI * stod(match.str(2));
    std::vector<std::vector<s_type>> ret_val =
        fSim(theta_rads, phi_rads, scratch);
    return std::tuple<std::vector<s_type>, std::vector<s_type>,
                      std::vector<size_t>>(ret_val[0], ret_val[1], {4, 2, 2});
  }
}

/**
 * Returns spatial coordinates i and j on the grid given qubit number q.
 * @param q int with the qubit number.
 * @param J int with the second spatial dimension of the grid of qubits.
 * @return int with the spatial coordinates i and j of qubit q on the grid.
 */
std::vector<int> _q_to_i_j(int q, int J) {
  int i = q / J;
  return std::vector<int>({i, q - i * J});
}

std::string tensor_name(std::vector<int> p1, std::vector<int> p2) {
  char buffer[64];
  if (p1.size() == 2 && p2.size() == 2) {
    // Two-qubit contraction.
    snprintf(buffer, sizeof(buffer), "(%d,%d),(%d,%d)", p1[0], p1[1], p2[0],
             p2[1]);
    return buffer;
  }
  if (p1.size() == 3 && p2.size() == 3) {
    // Single-qubit contraction, or virtual index.
    snprintf(buffer, sizeof(buffer), "(%d,%d,%d),(%d,%d,%d)", p1[0], p1[1],
             p1[2], p2[0], p2[1], p2[2]);
    return buffer;
  }
  // Final qubit output value assignment.
  if (p1.size() == 2 && p2.empty()) {
    snprintf(buffer, sizeof(buffer), "(%d,%d),(o)", p1[0], p1[1]);
    return buffer;
  }
  assert(false && "Failed to construct tensor name.");
  return "";
}

/**
 * Helper function to find a grid coordinate in a list of coordinates.
 */
bool find_grid_coord_in_list(
    const std::optional<std::vector<std::vector<int>>>& coord_list, const int i,
    const int j) {
  return coord_list.has_value() &&
         find(coord_list.value().begin(), coord_list.value().end(),
              std::vector<int>({i, j})) != coord_list.value().end();
}

/**
 * Helper method for grid_of_tensors_3D_to_2D which determines the order for
 * tensor indices around a target qubit based on contraction order and cuts.
 * @param ordering vector<vector<vector<int>>> listing the order in which
 * contractions are to be performed.
 * @param cuts vector<vector<vector<int>>> listing the cuts applied to the grid.
 * @param local vector<int> of the target qubit coordinates.
 * @return function for use as a comparator in std::sort.
 */
std::function<bool(std::vector<int>, std::vector<int>)> order_func(
    const std::vector<std::vector<std::vector<int>>>& ordering,
    const std::vector<std::vector<std::vector<int>>>& cuts,
    std::vector<int> local) {
  return [ordering, cuts, local](const std::vector<int> lhs,
                                 const std::vector<int> rhs) {
    // Cuts are projected early and should be ordered first.
    std::vector<std::vector<int>> lhs_pair, rhs_pair;
    if (local[0] < lhs[0] || local[1] < lhs[1]) {
      lhs_pair = {local, lhs};
    } else {
      lhs_pair = {lhs, local};
    }
    if (std::find(cuts.begin(), cuts.end(), lhs_pair) != cuts.end()) {
      return true;
    }
    if (local[0] < rhs[0] || local[1] < rhs[1]) {
      rhs_pair = {local, rhs};
    } else {
      rhs_pair = {rhs, local};
    }
    if (std::find(cuts.begin(), cuts.end(), rhs_pair) != cuts.end()) {
      return false;
    }
    std::vector<int> lpos, rpos;
    for (int i = 0; i < ordering.size(); i++) {
      const auto& outer = ordering[i];
      auto find_lpos = std::find(outer.begin(), outer.end(), lhs);
      if (find_lpos != outer.end()) {
        lpos = {i, std::distance(outer.begin(), find_lpos)};
        break;
      }
    }
    if (lpos.empty()) {
      char error[200];
      snprintf(error, sizeof(error), "Pair not found: (%d,%d),(%d,%d)",
               local[0], local[1], lhs[0], lhs[1]);
      std::cout << error << std::endl;
    }
    for (int i = 0; i < ordering.size(); i++) {
      const auto& outer = ordering[i];
      auto find_rpos = std::find(outer.begin(), outer.end(), rhs);
      if (find_rpos != outer.end()) {
        rpos = {i, std::distance(outer.begin(), find_rpos)};
        break;
      }
    }
    if (rpos.empty()) {
      char error[200];
      snprintf(error, sizeof(error), "Pair not found: (%d,%d),(%d,%d)",
               local[0], local[1], rhs[0], rhs[1]);
      std::cout << error << std::endl;
    }
    if (lpos[0] == rpos[0]) {
      // lhs and rhs are in the same patch.
      return lpos[1] < rpos[1];
    }

    int local_patch;
    for (int i = 0; i < ordering.size(); i++) {
      const auto& outer = ordering[i];
      if (std::find(outer.begin(), outer.end(), local) != outer.end()) {
        local_patch = i;
        break;
      }
    }
    if (lpos[0] == local_patch) {
      // rhs won't be connected until patches merge.
      return true;
    }
    if (rpos[0] == local_patch) {
      // lhs won't be connected until patches merge.
      return false;
    }
    // Both lhs and rhs are in different patches from local_patch.
    return lhs[0] < rhs[0];
  };
}

}  // namespace

void circuit_data_to_grid_of_tensors(
    std::istream* circuit_data, int I, int J, int K,
    const std::string initial_conf, const std::string final_conf_B,
    const std::optional<std::vector<std::vector<int>>>& A,
    const std::optional<std::vector<std::vector<int>>>& off,
    std::vector<std::vector<std::vector<MKLTensor>>>& grid_of_tensors,
    s_type* scratch) {
  // Gotten from the file.
  int num_qubits, cycle, q1, q2;
  std::string gate;
  // Useful for plugging into the tensor network:
  std::vector<int> i_j_1, i_j_2;
  int super_cycle;

  // The first element should be the number of qubits
  *(circuit_data) >> num_qubits;
  num_qubits = I * J;

  // Assert for the number of qubits.
  assert(num_qubits == I * J && "I*J must be equal to the number of qubits.");
  // Assert for the length of initial_conf and final_conf_B.
  {
    size_t off_size = off.has_value() ? off.value().size() : 0;
    size_t A_size = A.has_value() ? A.value().size() : 0;
    assert(initial_conf.size() == num_qubits - off_size &&
           "initial_conf must be of size equal to the number of qubits.");
    assert(final_conf_B.size() == num_qubits - off_size - A_size &&
           "final_conf_B must be of size equal to the number of qubits.");
  }

  // Creating grid variables.
  std::vector<std::vector<std::vector<std::vector<MKLTensor>>>>
      grid_of_groups_of_tensors(I);
  grid_of_tensors = std::vector<std::vector<std::vector<MKLTensor>>>(I);
  std::vector<std::vector<std::vector<int>>> counter_group(I);
  for (int i = 0; i < I; ++i) {
    grid_of_groups_of_tensors[i] =
        std::vector<std::vector<std::vector<MKLTensor>>>(J);
    grid_of_tensors[i] = std::vector<std::vector<MKLTensor>>(J);
    counter_group[i] = std::vector<std::vector<int>>(J);
    for (int j = 0; j < J; ++j) {
      grid_of_groups_of_tensors[i][j] = std::vector<std::vector<MKLTensor>>(K);
      grid_of_tensors[i][j] = std::vector<MKLTensor>(K);
      counter_group[i][j] = std::vector<int>(K, 0);
      for (int k = 0; k < K; ++k) {
        grid_of_groups_of_tensors[i][j][k] = std::vector<MKLTensor>();
      }
    }
  }

  // Insert deltas and Hadamards to first layer.
  int idx = 0;
  for (int q = 0; q < num_qubits; ++q) {
    std::vector<int> i_j = _q_to_i_j(q, J);
    int i = i_j[0], j = i_j[1];
    if (find_grid_coord_in_list(off, i, j)) {
      continue;
    }
    std::string delta_gate = (initial_conf[idx] == '0') ? "delta_0" : "delta_1";
    grid_of_groups_of_tensors[i][j][0].push_back(
        MKLTensor({"th"}, {2}, gate_array(delta_gate)));
    grid_of_groups_of_tensors[i][j][0].push_back(
        MKLTensor({"th", "t0"}, {2, 2}, gate_array("h")));
    idx += 1;
  }

  std::string line;
  // Read one line at a time from the circuit, skipping comments.
  while (getline(*circuit_data, line))
    if (line.size() && line[0] != '#') {
      std::stringstream ss(line);
      // The first element is the cycle
      ss >> cycle;
      // The second element is the gate
      ss >> gate;
      // Get the first position
      // TODO: Stop relying on the assumption that the gate and its parameters
      // can be read as one token without spaces. This is (mostly) fine for
      // "rz(0.5)", but will fail for, e.g., "fsim(0.25, -0.5)".
      ss >> q1;
      // Get the second position in the case
      // TODO: Two-qubit gates should be encapsulated better.
      if (gate == "cz" || gate.rfind("fsim", 0) == 0) {
        ss >> q2;
      } else {
        q2 = -1;
      }

      // Get i, j and super_cycle
      i_j_1 = _q_to_i_j(q1, J);
      if (q2 >= 0) {
        i_j_2 = _q_to_i_j(q2, J);
      }
      super_cycle = (cycle - 1) / SUPER_CYCLE_DEPTH;

      // Fill in one-qubit gates.
      if (q2 < 0 && cycle > 0 && cycle <= SUPER_CYCLE_DEPTH * K) {
        if (find_grid_coord_in_list(off, i_j_1[0], i_j_1[1])) {
          continue;
        }
        std::string input_index =
            "t" +
            std::to_string(counter_group[i_j_1[0]][i_j_1[1]][super_cycle]);
        std::string output_index =
            "t" +
            std::to_string(counter_group[i_j_1[0]][i_j_1[1]][super_cycle] + 1);
        ++counter_group[i_j_1[0]][i_j_1[1]][super_cycle];
        grid_of_groups_of_tensors[i_j_1[0]][i_j_1[1]][super_cycle].push_back(
            MKLTensor({input_index, output_index}, {2, 2}, gate_array(gate)));
      }
      if (q2 >= 0 && cycle > 0 && cycle <= SUPER_CYCLE_DEPTH * K) {
        if (find_grid_coord_in_list(off, i_j_1[0], i_j_1[1]) ||
            find_grid_coord_in_list(off, i_j_2[0], i_j_2[1])) {
          continue;
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
            tensor_name({i_j_1[0], i_j_1[1], super_cycle},
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
            MKLTensor({input_index_1, virtual_index, output_index_1},
                      dimensions, gate_q1));
        grid_of_groups_of_tensors[i_j_2[0]][i_j_2[1]][super_cycle].push_back(
            MKLTensor({input_index_2, virtual_index, output_index_2},
                      dimensions, gate_q2));
      }
    }
  // Insert Hadamards and deltas to last layer.
  idx = 0;
  for (int q = 0; q < num_qubits; ++q) {
    std::vector<int> i_j = _q_to_i_j(q, J);
    int i = i_j[0], j = i_j[1];
    int k = K - 1;
    if (find_grid_coord_in_list(off, i, j)) {
      continue;
    }
    std::string last_index = "t" + std::to_string(counter_group[i][j][k]);
    grid_of_groups_of_tensors[i][j][k].push_back(
        MKLTensor({"th", last_index}, {2, 2}, gate_array("h")));
    if (find_grid_coord_in_list(A, i, j)) {
      continue;
    }
    std::string delta_gate = (final_conf_B[idx] == '0') ? "delta_0" : "delta_1";
    grid_of_groups_of_tensors[i][j][k].push_back(
        MKLTensor({"th"}, {2}, gate_array(delta_gate)));
    idx += 1;  // Move in B only.
  }

  // Contracting each group of gates into a single tensor.
  for (int i = 0; i < I; ++i)
    for (int j = 0; j < J; ++j)
      for (int k = 0; k < K; ++k) {
        if (find_grid_coord_in_list(off, i, j)) {
          continue;
        }

        std::vector<MKLTensor>& group = grid_of_groups_of_tensors[i][j][k];
        std::vector<MKLTensor> group_containers(
            SUPER_CYCLE_DEPTH + 2,  // +2 for d and H.
            MKLTensor({""}, {(int)pow(DIM, 6)}));
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
          std::string new_first_index = tensor_name({i, j, k - 1}, {i, j, k});
          grid_of_tensors[i][j][k].rename_index("t0", new_first_index);
        }
        if (k < K - 1) {
          std::string last_index = "t" + std::to_string(counter_group[i][j][k]);
          std::string new_last_index = tensor_name({i, j, k}, {i, j, k + 1});
          grid_of_tensors[i][j][k].rename_index(last_index, new_last_index);
        }
        if (k == K - 1 && find_grid_coord_in_list(A, i, j)) {
          std::string last_index = "th";
          std::string new_last_index = tensor_name({i, j}, {});
          grid_of_tensors[i][j][k].rename_index(last_index, new_last_index);
        }
      }

  // Be proper about pointers.
  scratch = NULL;
}

void google_circuit_file_to_grid_of_tensors(
    std::string filename, int I, int J, int K, const std::string initial_conf,
    const std::string final_conf_B,
    const std::optional<std::vector<std::vector<int>>>& A,
    const std::optional<std::vector<std::vector<int>>>& off,
    std::vector<std::vector<std::vector<MKLTensor>>>& grid_of_tensors,
    s_type* scratch) {
  // Open file.
  auto io = std::ifstream(filename);
  assert(io.good() && "Cannot open file.");
  circuit_data_to_grid_of_tensors(&io, I, J, K, initial_conf, final_conf_B, A,
                                  off, grid_of_tensors, scratch);
}

void grid_of_tensors_3D_to_2D(
    std::vector<std::vector<std::vector<MKLTensor>>>& grid_of_tensors_3D,
    std::vector<std::vector<MKLTensor>>& grid_of_tensors_2D,
    std::optional<std::vector<std::vector<int>>> A,
    std::optional<std::vector<std::vector<int>>> off,
    const std::vector<std::vector<std::vector<int>>>& ordering,
    const std::vector<std::vector<std::vector<int>>>& cuts,
    s_type* scratch) {
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
      if (find_grid_coord_in_list(A, i, j)) {
        container_dim *= DIM;
      }
      std::vector<MKLTensor> group_containers =
          std::vector<MKLTensor>(2, MKLTensor({""}, {container_dim}));
      MKLTensor* source_container = &group_containers[0];
      MKLTensor* target_container = &group_containers[1];

      if (K == 1) {
        grid_of_tensors_2D[i][j] = grid_of_tensors_3D[i][j][0];
      } else {
        multiply(grid_of_tensors_3D[i][j][0], grid_of_tensors_3D[i][j][1],
                 *source_container, scratch);
        for (int k = 1; k < K - 1; ++k) {
          multiply(*source_container, grid_of_tensors_3D[i][j][k + 1],
                   *target_container, scratch);
          MKLTensor* swap_container = source_container;
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
      std::string index_name;

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
      if (find_grid_coord_in_list(A, i, j)) {
        ordered_indices_3D.push_back(tensor_name({i, j}, {}));
        fr_buffer = 1;
      }

      std::vector<int> local = {i, j};
      auto order_fn = order_func(ordering, cuts, local);
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
              tensor_name({q1[0], q1[1], k}, {q2[0], q2[1], k}));
        }
        indices_2D.push_back(tensor_name(q1, q2));
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

void read_wave_function_evolution(
    std::string filename, int I, std::vector<MKLTensor>& gates,
    std::vector<std::vector<std::string>>& inputs,
    std::vector<std::vector<std::string>>& outputs, s_type* scratch) {
  // Open file.
  auto io = std::ifstream(filename);
  assert(io.good() && "Cannot open file.");

  // Gotten from the file.
  int num_qubits, cycle, q1, q2;
  std::string gate;
  // Useful for plugging into the tensor network:
  std::vector<int> i_j_1, i_j_2;
  int super_cycle;

  // The first element should be the number of qubits
  io >> num_qubits;

  // Assert for the number of qubits.
  assert(num_qubits == I && "I must be equal to the number of qubits.");

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
      if (gate == "cz")
        ss >> q2;
      else
        q2 = -1;

      // Fill in one-qubit gates.
      if (q2 < 0) {
        std::string input_index = std::to_string(q1) + ",i";
        std::string output_index = std::to_string(q1) + ",o";
        gates.push_back(MKLTensor({input_index, output_index}, {DIM, DIM},
                                  gate_array(gate)));
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
        gates.push_back(MKLTensor(
            {input_index1, input_index2, output_index1, output_index2},
            {DIM, DIM, DIM, DIM}, gate_array(gate)));
      }
    }

  // Be proper about pointers.
  scratch = NULL;
}
