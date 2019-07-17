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
#include "mkl_tensor.h"
#include "read_circuit.h"

using ::qflex::ContractGrid;
using ::qflex::ContractionOperation;
using ::qflex::ContractionOrdering;
using ::qflex::CutIndex;
using ::qflex::DIM;
using ::qflex::ExpandPatch;
using ::qflex::MergePatches;
using ::qflex::MKLTensor;
using ::qflex::s_type;

// TODO(martinop): move these methods to a library and add tests.
std::vector<std::vector<int>> read_grid_layout_from_file(
    int I, int J, std::string grid_filename) {
  auto io = std::ifstream(grid_filename);
  assert(io.good() && "Cannot open grid file.");
  std::vector<std::vector<int>> qubits_off;
  bool on;
  for (int i = 0; i < I; ++i) {
    for (int j = 0; j < J; ++j) {
      io >> on;
      if (on == 0) qubits_off.push_back({i, j});
    }
  }
  return qubits_off;
}

std::vector<std::vector<int>> get_final_qubits_from_ordering(
    const ContractionOrdering& ordering) {
  std::vector<std::vector<int>> final_qubits;
  for (const auto& op : ordering) {
    if (op->op_type != ContractionOperation::CUT) continue;
    const auto* cut = dynamic_cast<const CutIndex*>(op.get());
    // Any qubit with a terminal cut is in the final region.
    // TODO(martinop): update to use the new operation.
    if (cut->tensors.size() == 1) {
      final_qubits.push_back(cut->tensors[0]);
    }
  }
  return final_qubits;
}

// Input: ./qflex.x I J K fidelity circuit_filename ordering_filename \
//            grid_filename [initial_conf] [final_conf]
//
// Example:
// $ ./qflex.x 11 12 2 0.005 ./circuits/ben_11_16_0.txt \
//       ./ordering/bristlecone_48.txt ./grid/bristlecone_48.txt
int main(int argc, char** argv) {
  // Set precision for the printed floats.
  std::cout.precision(12);

  // Timing variables.
  std::chrono::high_resolution_clock::time_point t_output_0, t_output_1;
  t_output_0 = std::chrono::high_resolution_clock::now();
  std::chrono::high_resolution_clock::time_point t0, t1;
  std::chrono::duration<double> time_span;

  // Reading input.
  t0 = std::chrono::high_resolution_clock::now();
  if (argc < 8) throw std::logic_error("ERROR: Not enough arguments.");
  int current_arg = 1;
  const int I = atoi(argv[current_arg++]);
  const int J = atoi(argv[current_arg++]);
  const int K = atoi(argv[current_arg++]);
  double fidelity = atof(argv[current_arg++]);
  const int super_dim = (int)pow(DIM, K);
  const std::string circuit_filename = std::string(argv[current_arg++]);
  const std::string ordering_filename = std::string(argv[current_arg++]);
  const std::string grid_filename = std::string(argv[current_arg++]);
  auto qubits_off = read_grid_layout_from_file(I, J, grid_filename);
  int init_length = I * J - qubits_off.size();
  int final_region_size = 8;
  std::string initial_conf(init_length, '0');
  std::string final_conf_B(init_length - final_region_size, '0');
  std::vector<std::string> final_conf_A(1, std::string(final_region_size, '0'));
  if (argc > current_arg) initial_conf = std::string(argv[current_arg++]);
  if (argc > current_arg) final_conf_B = std::string(argv[current_arg++]);
  if (argc > current_arg) {
    final_conf_A = std::vector<std::string>(argc - current_arg);
    for (int s = 0; s < final_conf_A.size(); ++s)
      final_conf_A[s] = std::string(argv[s + current_arg]);
  }
  const int num_Cs = final_conf_A.size();
  int num_qubits = I * J;
  t1 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
  // std::cout << "Time spent reading input: "
  //     << time_span.count()
  //     << "s\n\n";

  // Create the ordering for this tensor contraction from file.
  ContractionOrdering ordering;
  google_ordering_file_to_contraction_ordering(ordering_filename, I, J,
                                               qubits_off, &ordering);
  // Get a list of qubits in the final region.
  auto final_qubits = get_final_qubits_from_ordering(ordering);

  // Scratch space to be reused for operations.
  t0 = std::chrono::high_resolution_clock::now();
  s_type* scratch = new s_type[(int)pow(super_dim, 7)];
  t1 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
  // std::cout << "Time spent reading allocating scratch space: "
  //     << time_span.count()
  //     << "s\n\n";

  // Declaring and then filling 2D grid of tensors.
  std::vector<std::vector<MKLTensor>> tensor_grid(I);
  for (int i = 0; i < I; ++i) {
    tensor_grid[i] = std::vector<MKLTensor>(J);
  }
  // Scope so that the 3D grid of tensors is destructed.
  {
    // Creating 3D grid of tensors from file.
    t0 = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<std::vector<MKLTensor>>> tensor_grid_3D;
    google_circuit_file_to_grid_of_tensors(
        circuit_filename, I, J, K, initial_conf, final_conf_B, final_qubits,
        qubits_off, tensor_grid_3D, scratch);
    t1 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    // std::cout << "Time spent creating 3D grid of tensors from file: "
    //     << time_span.count()
    //     << "s\n\n";

    // Contract 3D grid onto 2D grid of tensors, as usual.
    t0 = std::chrono::high_resolution_clock::now();
    grid_of_tensors_3D_to_2D(tensor_grid_3D, tensor_grid, final_qubits,
                             qubits_off, ordering, scratch);
    t1 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    // std::cout << "Time spent creating 2D grid of tensors from 3D one: "
    //     << time_span.count()
    //     << "s\n\n";

    // Freeing scratch data: delete and NULL.
    delete[] scratch;
    scratch = NULL;
  }

  // Perform tensor grid contraction.
  std::vector<std::complex<double>> amplitudes(num_Cs);
  ContractGrid(ordering, &tensor_grid, &amplitudes);

  // Printing output
  for (int c = 0; c < num_Cs; ++c) {
    std::cout << initial_conf << " ";
    std::cout << final_conf_B << " ";
    std::cout << final_conf_A[c] << " ";
    std::cout << std::real(amplitudes[c]) << " " << std::imag(amplitudes[c]);
    std::cout << std::endl;
  }

  // Final time
  t_output_1 = std::chrono::high_resolution_clock::now();
  time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
      t_output_1 - t_output_0);
  std::cout << time_span.count() << " s\n\n";

  return 0;
}
