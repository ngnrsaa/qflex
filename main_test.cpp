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

  // Set precision for the printed floats.
  cout.precision(12);

  // Reading input.
  if (argc<6) throw logic_error("ERROR: Not enough arguments.");
  const int I = atoi(argv[1]);
  const int J = atoi(argv[2]);
  const int K = atoi(argv[3]);
  size_t super_dim = (size_t)pow(DIM,K);
  double fidelity = atof(argv[4]);
  const string filename = string(argv[5]);
  string initial_conf(9, '0'), final_conf_B(9, '0');
  vector<string> final_conf_A(1, string(0, '0'));
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
    {0,0},{0,1},{0,2},{0,3},{0,4},{0,5},{0,6},{0,7},{0,8},{0,9},{0,10},{0,11},
    {1,0},{1,1},{1,2},{1,3},{1,4},{1,5},{1,6},{1,7},{1,8},{1,9},{1,10},{1,11},
    {2,0},{2,1},{2,2},{2,3},{2,4},{2,5},{2,6},{2,7},{2,8},{2,9},{2,10},{2,11},
    {3,0},{3,1},{3,2},{3,3},{3,4},{3,5},{3,6},{3,7},{3,8},{3,9},{3,10},{3,11},
    {4,0},{4,1},{4,2},{4,3},{4,4},{4,5},{4,6},{4,7},{4,8},{4,9},{4,10},{4,11},
    {5,0},{5,1},{5,2},{5,3},{5,4},{5,5},{5,6},{5,7},{5,8},{5,9},{5,10},{5,11},
    {6,0},{6,1},{6,2}                  ,{6,6},{6,7},{6,8},{6,9},{6,10},{6,11},
    {7,0},{7,1},{7,2}                  ,{7,6},{7,7},{7,8},{7,9},{7,10},{7,11},
    {8,0},{8,1},{8,2}                  ,{8,6},{8,7},{8,8},{8,9},{8,10},{8,11},
    {9,0},{9,1},{9,2},{9,3},{9,4},{9,5},{9,6},{9,7},{9,8},{9,9},{9,10},{9,11},
    {10,0},{10,1},{10,2},{10,3},{10,4},{10,5},{10,6},{10,7},{10,8},{10,9},{10,10},{10,11}
    }); 
  vector<vector<int>> qubits_A({});


  // Initialize TALSH
  unsigned long size(2000000000);
  talsh::initialize(&size);
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

    // Allocate talsh::Tensors involved in the contraction
    vector<int> dims_1(1, super_dim);
    vector<int> dims_2(2, super_dim);
    vector<int> dims_3(3, super_dim);
    vector<int> dims_4(4, super_dim);
    size_t vol_1 = (size_t)pow(super_dim,1);
    size_t vol_2 = (size_t)pow(super_dim,2);
    size_t vol_3 = (size_t)pow(super_dim,3);
    size_t vol_4 = (size_t)pow(super_dim,4);
    talsh::Tensor H_1_legs_a(dims_1, s_type(0.0));
    talsh::Tensor H_2_legs_a(dims_2, s_type(0.0));
    talsh::Tensor H_2_legs_b(dims_2, s_type(0.0));
    talsh::Tensor H_3_legs_a(dims_3, s_type(0.0));
    talsh::Tensor H_3_legs_b(dims_3, s_type(0.0));
    talsh::Tensor H_4_legs_a(dims_4, s_type(0.0));
    talsh::Tensor H_4_legs_b(dims_4, s_type(0.0));
    // Finally, scalar S
    talsh::Tensor S({}, s_type(0.0));

    // Print tensors
    /*
    s_type const * ptr_T;
    tensor_grid[6][5]->getDataAccessHostConst(&ptr_T);
    int size_T = tensor_grid[6][5]->getVolume();
    for (int t=0; t<size_T; ++t)
      cout << (ptr_T[t]) << " ";
    cout << "\n\n";
    ptr_T = nullptr;
    */

    // Start contracting.
    TensContraction tc("D(a,c,d)+=L(b,c,d)*R(a,b)", &H_3_legs_a,
                        tensor_grid[6][4].get(), tensor_grid[6][3].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,d)+=L(c,d)*R(a,b,c)", &H_3_legs_b,
                        tensor_grid[6][5].get(), &H_3_legs_a);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(d,e,b,c)+=L(a,d,e)*R(a,b,c)", &H_4_legs_a,
                        tensor_grid[7][3].get(), &H_3_legs_b);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,e,f,d)+=L(c,b,e,f)*R(a,b,c,d)", &H_4_legs_b,
                        tensor_grid[7][4].get(), &H_4_legs_a);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,e)+=L(d,c,e)*R(a,b,c,d)", &H_3_legs_a,
                        tensor_grid[7][5].get(), &H_4_legs_b);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(d,b,c)+=L(a,d)*R(a,b,c)", &H_3_legs_b,
                        tensor_grid[8][3].get(), &H_3_legs_a);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(d,c)+=L(b,a,d)*R(a,b,c)", &H_2_legs_a,
                        tensor_grid[8][4].get(), &H_3_legs_b);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D()+=L(b,a)*R(a,b)", &S,
                        tensor_grid[8][5].get(), &H_2_legs_a);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));

    // Store result data_S
    s_type const * ptr_S;
    S.getDataAccessHostConst(&ptr_S);
    cout << (ptr_S[0]) << endl;
    ptr_S = nullptr;
    /*
    */

  }
  // Shut down TALSH
  talsh::shutdown();


  return 0;
} 
