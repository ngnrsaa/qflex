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
  if (argc<6) throw logic_error("ERROR: Not enough arguments.");
  const int I = atoi(argv[1]);
  const int J = atoi(argv[2]);
  const int K = atoi(argv[3]);
  size_t super_dim = (size_t)pow(DIM,K);
  double fidelity = atof(argv[4]);
  const string filename = string(argv[5]);
  string initial_conf(52, '0'), final_conf_B(46, '0');
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
    {2,0},{2,1},{2,2}                                    ,{2,9},{2,10},{2,11},
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


  // Initialize TALSH
  talsh::initialize();
  {
    int errc;

    // Declaring and then filling 2D grid of tensors.
    vector<vector<talsh::Tensor *>> tensor_grid(I);
    for (int i=0; i<I; ++i)
    {
      tensor_grid[i] = vector<talsh::Tensor *>(J);
    }
    vector<s_type *> delete_this;
    vector<size_t> volumes;
    google_circuit_file_to_grid_of_tensors(filename, I, J, initial_conf,
                 final_conf_B, qubits_A, qubits_off, tensor_grid, delete_this,
                 volumes);

    // Allocate talsh::Tensors involved in the contraction
    vector<int> dims_2(2, super_dim);
    vector<int> dims_3(3, super_dim);
    vector<int> dims_4(4, super_dim);
    vector<int> dims_5(5, super_dim);
    vector<int> dims_6(6, super_dim);
    vector<int> dims_7(7, super_dim);
    size_t vol_2 = (size_t)pow(super_dim,2);
    size_t vol_3 = (size_t)pow(super_dim,3);
    size_t vol_4 = (size_t)pow(super_dim,4);
    size_t vol_5 = (size_t)pow(super_dim,5);
    size_t vol_6 = (size_t)pow(super_dim,6);
    size_t vol_7 = (size_t)pow(super_dim,7);
    s_type * data_2_legs_a = (s_type *) malloc(vol_2*sizeof(s_type));
    errc = talsh::pinHostMemory(data_2_legs_a,vol_2*sizeof(s_type));
    assert(errc==0);
    s_type * data_2_legs_b = (s_type *) malloc(vol_2*sizeof(s_type));
    errc = talsh::pinHostMemory(data_2_legs_b,vol_2*sizeof(s_type));
    assert(errc==0);
    s_type * data_3_legs_a = (s_type *) malloc(vol_3*sizeof(s_type));
    errc = talsh::pinHostMemory(data_3_legs_a,vol_2*sizeof(s_type));
    assert(errc==0);
    s_type * data_3_legs_b = (s_type *) malloc(vol_3*sizeof(s_type));
    errc = talsh::pinHostMemory(data_3_legs_b,vol_2*sizeof(s_type));
    assert(errc==0);
    s_type * data_4_legs_a = (s_type *) malloc(vol_4*sizeof(s_type));
    errc = talsh::pinHostMemory(data_4_legs_a,vol_2*sizeof(s_type));
    assert(errc==0);
    s_type * data_4_legs_b = (s_type *) malloc(vol_4*sizeof(s_type));
    errc = talsh::pinHostMemory(data_4_legs_b,vol_2*sizeof(s_type));
    assert(errc==0);
    s_type * data_5_legs_a = (s_type *) malloc(vol_5*sizeof(s_type));
    errc = talsh::pinHostMemory(data_5_legs_a,vol_2*sizeof(s_type));
    assert(errc==0);
    s_type * data_5_legs_b = (s_type *) malloc(vol_5*sizeof(s_type));
    errc = talsh::pinHostMemory(data_5_legs_b,vol_2*sizeof(s_type));
    assert(errc==0);
    s_type * data_6_legs_a = (s_type *) malloc(vol_6*sizeof(s_type));
    errc = talsh::pinHostMemory(data_6_legs_a,vol_2*sizeof(s_type));
    assert(errc==0);
    s_type * data_6_legs_b = (s_type *) malloc(vol_6*sizeof(s_type));
    errc = talsh::pinHostMemory(data_6_legs_b,vol_2*sizeof(s_type));
    assert(errc==0);
    s_type * data_7_legs_a = (s_type *) malloc(vol_7*sizeof(s_type));
    errc = talsh::pinHostMemory(data_7_legs_a,vol_2*sizeof(s_type));
    assert(errc==0);
    s_type * data_7_legs_b = (s_type *) malloc(vol_7*sizeof(s_type));
    errc = talsh::pinHostMemory(data_7_legs_b,vol_2*sizeof(s_type));
    assert(errc==0);
    delete_this.push_back(data_2_legs_a);
    delete_this.push_back(data_2_legs_b);
    delete_this.push_back(data_3_legs_a);
    delete_this.push_back(data_3_legs_b);
    delete_this.push_back(data_4_legs_a);
    delete_this.push_back(data_4_legs_b);
    delete_this.push_back(data_5_legs_a);
    delete_this.push_back(data_5_legs_b);
    delete_this.push_back(data_6_legs_a);
    delete_this.push_back(data_6_legs_b);
    delete_this.push_back(data_7_legs_a);
    delete_this.push_back(data_7_legs_b);
    volumes.push_back(vol_2);
    volumes.push_back(vol_2);
    volumes.push_back(vol_3);
    volumes.push_back(vol_3);
    volumes.push_back(vol_4);
    volumes.push_back(vol_4);
    volumes.push_back(vol_5);
    volumes.push_back(vol_5);
    volumes.push_back(vol_6);
    volumes.push_back(vol_6);
    volumes.push_back(vol_7);
    volumes.push_back(vol_7);
    talsh::Tensor H_2_legs_a(dims_2, data_2_legs_a);
    talsh::Tensor H_2_legs_b(dims_2, data_2_legs_b);
    talsh::Tensor H_3_legs_a(dims_3, data_3_legs_a);
    talsh::Tensor H_3_legs_b(dims_3, data_3_legs_b);
    talsh::Tensor H_4_legs_a(dims_4, data_4_legs_a);
    talsh::Tensor H_4_legs_b(dims_4, data_4_legs_b);
    talsh::Tensor H_5_legs_a(dims_5, data_5_legs_a);
    talsh::Tensor H_5_legs_b(dims_5, data_5_legs_b);
    talsh::Tensor H_6_legs_a(dims_6, data_6_legs_a);
    talsh::Tensor H_6_legs_b(dims_6, data_6_legs_b);
    talsh::Tensor H_7_legs_a(dims_7, data_7_legs_a);
    talsh::Tensor H_7_legs_b(dims_7, data_7_legs_b);


    // Start contracting
    // Level 1
    TensContraction tc("D(c,d,b)+=L(a,b)*R(a,c,d)",
                        &H_3_legs_a, tensor_grid[4][1], tensor_grid[5][1]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(d,b,c)+=L(a,b,c)*R(a,d)",
                        &H_3_legs_b, &H_3_legs_a, tensor_grid[6][1]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Level 2
    tc = TensContraction("D(f,a,b,d,e)+=L(a,b,c)*R(c,d,e,f)",
                        &H_5_legs_a, &H_3_legs_b, tensor_grid[4][2]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,d,f)+=L(a,b,c,d,e)*R(e,f)",
                        &H_5_legs_b, &H_5_legs_a, tensor_grid[3][2]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,f,g,d,e)+=L(a,b,c,d,e)*R(c,b,f,g)",
                        &H_5_legs_a, &H_5_legs_b, tensor_grid[5][2]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(f,g,c,d,e)+=L(a,b,c,d,e)*R(b,a,f,g)",
                        &H_5_legs_b, &H_5_legs_a, tensor_grid[6][2]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(f,b,c,d,e)+=L(a,b,c,d,e)*R(a,f)",
                        &H_5_legs_a, &H_5_legs_b, tensor_grid[7][2]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));

    // Level 3
    // Up to here
    /*
    */


    // Delete all pointers
    for (auto & v : delete_this)
    {
      errc = talsh::unpinHostMemory(v); assert(errc == 0);
      delete [] v; v = nullptr;
    }
    data_2_legs_a = nullptr;
    data_2_legs_b = nullptr;
    data_3_legs_a = nullptr;
    data_3_legs_b = nullptr;
    data_4_legs_a = nullptr;
    data_4_legs_b = nullptr;
    data_5_legs_a = nullptr;
    data_5_legs_b = nullptr;
    data_6_legs_a = nullptr;
    data_6_legs_b = nullptr;
    data_7_legs_a = nullptr;
    data_7_legs_b = nullptr;
  }
  // Shut down TALSH
  talsh::shutdown();

  

  return 0;
} 
