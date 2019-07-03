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

// Input: I J K fidelity filename initial_conf (optional) final_conf (optional)
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
  if (argc < 6) throw std::logic_error("ERROR: Not enough arguments.");
  const int I = atoi(argv[1]);
  const int J = atoi(argv[2]);
  const int K = atoi(argv[3]);
  double fidelity = atof(argv[4]);
  const int super_dim = (int)pow(DIM, K);
  const std::string filename = std::string(argv[5]);
  std::string initial_conf(70, '0'), final_conf_B(62, '0');
  std::vector<std::string> final_conf_A(1, std::string(8, '0'));
  if (argc > 6) initial_conf = std::string(argv[6]);
  if (argc > 7) final_conf_B = std::string(argv[7]);
  if (argc > 8) {
    final_conf_A = std::vector<std::string>(argc - 8);
    for (int s = 0; s < final_conf_A.size(); ++s)
      final_conf_A[s] = std::string(argv[s + 8]);
  }
  const int num_Cs = final_conf_A.size();
  int num_qubits = I * J;
  t1 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
  // cout << "Time spent reading input: "
  //     << time_span.count()
  //     << "s\n\n";

  // clang-format off
  // List of qubits to remove (off) and qubits in A.
  std::vector<std::vector<int>> qubits_off({
    {0,0},{0,1},{0,2},{0,3},{0,4}            ,{0,7},{0,8},{0,9},{0,10},{0,11},
    {1,0},{1,1},{1,2},{1,3}                        ,{1,8},{1,9},{1,10},{1,11},
    {2,0},{2,1},{2,2}                                    ,{2,9},{2,10},{2,11},
    {3,0},{3,1}                                                ,{3,10},{3,11},
    {4,0}                                                             ,{4,11},
    {5,0}                                                             ,{5,11},
    {6,0}                                                             ,{6,11},
    {7,0},{7,1}                                                ,{7,10},{7,11},
    {8,0},{8,1},{8,2},                                    {8,9},{8,10},{8,11},
    {9,0},{9,1},{9,2},{9,3},                        {9,8},{9,9},{9,10},{9,11},
    {10,0},{10,1},{10,2},{10,3},{10,4},              {10,7},{10,8},{10,9},{10,10},{10,11}
        });
  // clang-format on
  std::vector<std::vector<int>> qubits_A(
      {{3, 9}, {4, 9}, {4, 10}, {5, 9}, {5, 10}, {6, 9}, {6, 10}, {7, 9}});

  // Construct the ordering for this tensor contraction.
  ContractionOrdering ordering;
  const std::vector<std::vector<std::vector<int>>> cuts_a = {
      {{4, 1}, {5, 1}}, {{4, 2}, {5, 2}}, {{4, 3}, {5, 3}}};
  for (const auto& cut : cuts_a) {
    ordering.emplace_back(new CutIndex(cut, {0}));
  }
  // This cut has different values and must be created separately.
  ordering.emplace_back(new CutIndex({{4, 4}, {5, 4}}, {0, 1, 2}));
  const std::vector<std::vector<int>> order_a = {
      {5, 1}, {6, 1},  {5, 2},  {6, 2}, {7, 2}, {5, 3}, {6, 3}, {7, 3}, {8, 3},
      {5, 4}, {6, 4},  {7, 4},  {8, 4}, {9, 4}, {5, 5}, {6, 5}, {7, 5}, {8, 5},
      {9, 5}, {10, 5}, {10, 6}, {9, 6}, {8, 6}, {7, 6}, {6, 6}, {5, 6}, {9, 7},
      {8, 7}, {7, 7},  {6, 7},  {5, 7}, {8, 8}, {7, 8}, {6, 8}, {5, 8}};
  for (const auto& coord : order_a) {
    ordering.emplace_back(new ExpandPatch("A", coord));
  }
  const std::vector<std::vector<int>> order_b = {
      {4, 1}, {4, 2}, {3, 2}, {4, 3}, {3, 3}, {2, 3}, {4, 4}, {3, 4}, {2, 4},
      {1, 4}, {4, 5}, {3, 5}, {2, 5}, {1, 5}, {0, 5}, {0, 6}, {1, 6}, {2, 6},
      {3, 6}, {4, 6}, {1, 7}, {2, 7}, {3, 7}, {4, 7}, {2, 8}, {3, 8}, {4, 8}};
  for (const auto& coord : order_b) {
    ordering.emplace_back(new ExpandPatch("B", coord));
  }
  ordering.emplace_back(new MergePatches("A", "B"));
  // Add a terminal cut for every qubit in the final region.
  for (const auto& cut : qubits_A) {
    ordering.emplace_back(new CutIndex({cut}, {0}));
  }
  const std::vector<std::vector<int>> order_c = {
      {4, 10}, {5, 10}, {6, 10}, {4, 9}, {3, 9}, {5, 9}, {6, 9}, {7, 9}};
  for (const auto& coord : order_c) {
    ordering.emplace_back(new ExpandPatch("C", coord));
  }
  ordering.emplace_back(new MergePatches("B", "C"));

  // Scratch space to be reused for operations.
  t0 = std::chrono::high_resolution_clock::now();
  s_type* scratch = new s_type[(int)pow(super_dim, 7)];
  t1 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
  // cout << "Time spent reading allocating scratch space: "
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
    google_circuit_file_to_grid_of_tensors(filename, I, J, K, initial_conf,
                                           final_conf_B, qubits_A, qubits_off,
                                           tensor_grid_3D, scratch);
    t1 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    // cout << "Time spent creating 3D grid of tensors from file: "
    //     << time_span.count()
    //     << "s\n\n";

    // Contract 3D grid onto 2D grid of tensors, as usual.
    t0 = std::chrono::high_resolution_clock::now();
    grid_of_tensors_3D_to_2D(tensor_grid_3D, tensor_grid, qubits_A, qubits_off,
                             ordering, scratch);
    t1 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    // cout << "Time spent creating 2D grid of tensors from 3D one: "
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
