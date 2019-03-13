#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <complex>
#include <ctime>
#include <chrono>

#include <iostream>
#include <fstream>
#include <sstream>

#include "read_circuit_full_talsh.h"

#include "talshxx.hpp"
#include "talsh_wrapper.h"

#include <omp.h>

using namespace std;
using namespace chrono;


// Input: I J K fidelity filename initial_conf (optional) final_conf (optional)
int main(int argc, char **argv) {

  // Reading input.
  if (argc<6) throw logic_error("ERROR: Not enough arguments.");
  const int I = atoi(argv[1]);
  const int J = atoi(argv[2]);
  const int K = atoi(argv[3]);
  size_t super_dim = (size_t)pow(DIM,K);
  double fidelity = atof(argv[4]);
  const string filename = string(argv[5]);
  string initial_conf(51, '0'), final_conf_B(45, '0');
  vector<string> final_conf_A(1, string(6, '0'));
  if (argc>6)
    initial_conf = string(argv[6]);
  if (argc>7)
    final_conf_B = string(argv[7]);
  if (argc>8)
  {
    final_conf_A = vector<string>(argc-8);
    for (int s=0; s<final_conf_A.size(); ++s)
      final_conf_A[s] = string(argv[s+8]);
  }
  const int num_Cs = final_conf_A.size();
  int num_qubits = I*J;

  // List of qubits to remove (off) and qubits in A.
  vector<vector<int>> qubits_off({
    {0,0},{0,1},{0,2},{0,3},{0,4}            ,{0,7},{0,8},{0,9},{0,10},{0,11},
    {1,0},{1,1},{1,2},{1,3}                        ,{1,8},{1,9},{1,10},{1,11},
    {2,0},{2,1},{2,2},{2,3}                              ,{2,9},{2,10},{2,11},
    {3,0},{3,1}                                                ,{3,10},{3,11},
    {4,0}                                                      ,{4,10},{4,11},
    {5,0}                                                ,{5,9},{5,10},{5,11},
    {6,0}                                          ,{6,8},{6,9},{6,10},{6,11},
    {7,0},{7,1}                              ,{7,7},{7,8},{7,9},{7,10},{7,11},
    {8,0},{8,1},{8,2}                  ,{8,6},{8,7},{8,8},{8,9},{8,10},{8,11},
    {9,0},{9,1},{9,2},{9,3},{9,4},{9,5},{9,6},{9,7},{9,8},{9,9},{9,10},{9,11},
    {10,0},{10,1},{10,2},{10,3},{10,4},{10,5},{10,6},{10,7},{10,8},{10,9},{10,10},{10,11}
    }); 
  vector<vector<int>> qubits_A({       {0,5},{0,6},
                                 {1,4},{1,5},{1,6},{1,7} });


  // Vector of amplitudes with size num_Cs
  vector<s_type> amplitudes;

  // Initialize TALSH
  talsh::initialize();
  {
    int errc;

    // Declaring and then filling 2D grid of tensors.
    vector<vector<shared_ptr<talsh::Tensor>>> tensor_grid(I);
    for (int i=0; i<I; ++i)
    {
      tensor_grid[i] = vector<shared_ptr<talsh::Tensor>>(J);
    }
    google_circuit_file_to_grid_of_tensors(filename, I, J, initial_conf,
                 final_conf_B, qubits_A, qubits_off, tensor_grid);

  }
  // Shut down TALSH
  talsh::shutdown();


  return 0;
} 
