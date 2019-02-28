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

#include "read_circuit_no_mkl.h"

#include "talshxx.hpp"
#include "talsh_wrapper.h"

#include <omp.h>

using namespace std;
using namespace chrono;


// Input: I J K fidelity filename initial_conf (optional) final_conf (optional)
int main(int argc, char **argv) {

  // Reading input.
  if (argc<5) throw logic_error("ERROR: Not enough arguments.");
  const int I = atoi(argv[1]);
  const int J = atoi(argv[2]);
  double fidelity = atof(argv[3]);
  const string filename = string(argv[4]);
  string initial_conf(70, '0'), final_conf_B(62, '0');
  vector<string> final_conf_A(1, string(8, '0'));
  if (argc>5)
    initial_conf = string(argv[5]);
  if (argc>6)
    final_conf_B = string(argv[6]);
  if (argc>7)
  {
    final_conf_A = vector<string>(argc-7);
    for (int s=0; s<final_conf_A.size(); ++s)
      final_conf_A[s] = string(argv[s+7]);
  }
  const int num_Cs = final_conf_A.size();
  int num_qubits = I*J;

  // List of qubits to remove (off) and qubits in A.
  vector<vector<int>> qubits_off({
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
  vector<vector<int>> qubits_A({{3,9},{4,9},{4,10},{5,9},{5,10},
                                {6,9},{6,10},{7,9}});


  // Initialize TALSH
  talsh::initialize();
  {
    // Declaring and then filling 2D grid of tensors.
    vector<vector<talsh::Tensor *>> tensor_data_grid(I);
    for (int i=0; i<I; ++i)
    {
      tensor_data_grid[i] = vector<talsh::Tensor *>(J);
    }
    google_circuit_file_to_grid_of_tensors(filename, I, J, initial_conf,
                 final_conf_B, qubits_A, qubits_off, tensor_data_grid);
  }
  // Shut down TALSH
  talsh::shutdown();

  

  return 0;
} 
