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
    // First, tensors for region C. Done by hand right now. Change in future.
    vector<talsh::Tensor *> Cs;
    s_type * data_c05 = (s_type *) malloc(vol_2*sizeof(s_type));
    errc = talsh::pinHostMemory(data_c05,vol_2*sizeof(s_type));
    assert(errc==0);
    s_type * data_c06 = (s_type *) malloc(vol_2*sizeof(s_type));
    errc = talsh::pinHostMemory(data_c06,vol_2*sizeof(s_type));
    assert(errc==0);
    s_type * data_c14 = (s_type *) malloc(vol_2*sizeof(s_type));
    errc = talsh::pinHostMemory(data_c14,vol_2*sizeof(s_type));
    assert(errc==0);
    s_type * data_c15 = (s_type *) malloc(vol_4*sizeof(s_type));
    errc = talsh::pinHostMemory(data_c15,vol_4*sizeof(s_type));
    assert(errc==0);
    s_type * data_c16 = (s_type *) malloc(vol_4*sizeof(s_type));
    errc = talsh::pinHostMemory(data_c16,vol_4*sizeof(s_type));
    assert(errc==0);
    s_type * data_c17 = (s_type *) malloc(vol_2*sizeof(s_type));
    errc = talsh::pinHostMemory(data_c17,vol_2*sizeof(s_type));
    assert(errc==0);
    delete_this.push_back(data_c05);
    delete_this.push_back(data_c06);
    delete_this.push_back(data_c14);
    delete_this.push_back(data_c15);
    delete_this.push_back(data_c16);
    delete_this.push_back(data_c17);
    volumes.push_back(vol_2);
    volumes.push_back(vol_2);
    volumes.push_back(vol_2);
    volumes.push_back(vol_4);
    volumes.push_back(vol_4);
    volumes.push_back(vol_2);
    Cs.push_back(new talsh::Tensor(dims_2, data_c05));
    Cs.push_back(new talsh::Tensor(dims_2, data_c06));
    Cs.push_back(new talsh::Tensor(dims_2, data_c14));
    Cs.push_back(new talsh::Tensor(dims_4, data_c15));
    Cs.push_back(new talsh::Tensor(dims_4, data_c16));
    Cs.push_back(new talsh::Tensor(dims_2, data_c17));
    // Second, helper tensors
    s_type * data_2_legs_a = (s_type *) malloc(vol_2*sizeof(s_type));
    errc = talsh::pinHostMemory(data_2_legs_a,vol_2*sizeof(s_type));
    assert(errc==0);
    s_type * data_2_legs_b = (s_type *) malloc(vol_2*sizeof(s_type));
    errc = talsh::pinHostMemory(data_2_legs_b,vol_2*sizeof(s_type));
    assert(errc==0);
    s_type * data_3_legs_a = (s_type *) malloc(vol_3*sizeof(s_type));
    errc = talsh::pinHostMemory(data_3_legs_a,vol_3*sizeof(s_type));
    assert(errc==0);
    s_type * data_3_legs_b = (s_type *) malloc(vol_3*sizeof(s_type));
    errc = talsh::pinHostMemory(data_3_legs_b,vol_3*sizeof(s_type));
    assert(errc==0);
    s_type * data_4_legs_a = (s_type *) malloc(vol_4*sizeof(s_type));
    errc = talsh::pinHostMemory(data_4_legs_a,vol_4*sizeof(s_type));
    assert(errc==0);
    s_type * data_4_legs_b = (s_type *) malloc(vol_4*sizeof(s_type));
    errc = talsh::pinHostMemory(data_4_legs_b,vol_4*sizeof(s_type));
    assert(errc==0);
    s_type * data_4_legs_c = (s_type *) malloc(vol_4*sizeof(s_type));
    errc = talsh::pinHostMemory(data_4_legs_c,vol_4*sizeof(s_type));
    assert(errc==0);
    s_type * data_5_legs_a = (s_type *) malloc(vol_5*sizeof(s_type));
    errc = talsh::pinHostMemory(data_5_legs_a,vol_5*sizeof(s_type));
    assert(errc==0);
    s_type * data_5_legs_b = (s_type *) malloc(vol_5*sizeof(s_type));
    errc = talsh::pinHostMemory(data_5_legs_b,vol_5*sizeof(s_type));
    assert(errc==0);
    s_type * data_6_legs_a = (s_type *) malloc(vol_6*sizeof(s_type));
    errc = talsh::pinHostMemory(data_6_legs_a,vol_6*sizeof(s_type));
    assert(errc==0);
    s_type * data_6_legs_b = (s_type *) malloc(vol_6*sizeof(s_type));
    errc = talsh::pinHostMemory(data_6_legs_b,vol_6*sizeof(s_type));
    assert(errc==0);
    s_type * data_7_legs_a = (s_type *) malloc(vol_7*sizeof(s_type));
    errc = talsh::pinHostMemory(data_7_legs_a,vol_7*sizeof(s_type));
    assert(errc==0);
    s_type * data_7_legs_b = (s_type *) malloc(vol_7*sizeof(s_type));
    errc = talsh::pinHostMemory(data_7_legs_b,vol_7*sizeof(s_type));
    assert(errc==0);
    delete_this.push_back(data_2_legs_a);
    delete_this.push_back(data_2_legs_b);
    delete_this.push_back(data_3_legs_a);
    delete_this.push_back(data_3_legs_b);
    delete_this.push_back(data_4_legs_a);
    delete_this.push_back(data_4_legs_b);
    delete_this.push_back(data_4_legs_c);
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
    talsh::Tensor H_4_legs_c(dims_4, data_4_legs_c);
    talsh::Tensor H_5_legs_a(dims_5, data_5_legs_a);
    talsh::Tensor H_5_legs_b(dims_5, data_5_legs_b);
    talsh::Tensor H_6_legs_a(dims_6, data_6_legs_a);
    talsh::Tensor H_6_legs_b(dims_6, data_6_legs_b);
    talsh::Tensor H_7_legs_a(dims_7, data_7_legs_a);
    talsh::Tensor H_7_legs_b(dims_7, data_7_legs_b);
    // Finally, scalar S
    s_type * data_S = (s_type *) malloc(1*sizeof(s_type));
    errc = talsh::pinHostMemory(data_S,1*sizeof(s_type));
    assert(errc==0);
    delete_this.push_back(data_S);
    volumes.push_back(1);
    talsh::Tensor S({}, data_S);


    // Start contracting.
    // It's working, but I haven't respected the ordering criterion!
    // Level 1
    TensContraction tcg("D(c,d,b)+=L(a,b)*R(a,c,d)",
                        &H_3_legs_a, tensor_grid[4][1], tensor_grid[5][1]);
    errc = tcg.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tcg.sync(DEV_HOST,0));
    TensContraction tc("D(d,b,c)+=L(a,b,c)*R(a,d)",
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
    tc = TensContraction("D(a,b,c,d,f,g)+=L(a,b,c,d,e)*R(e,f,g)",
                        &H_6_legs_a, &H_5_legs_a, tensor_grid[3][3]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,g,h,f)+=L(a,b,c,d,e,f)*R(e,d,g,h)",
                        &H_6_legs_b, &H_6_legs_a, tensor_grid[4][3]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,g,h,e,f)+=L(a,b,c,d,e,f)*R(d,c,g,h)",
                        &H_6_legs_a, &H_6_legs_b, tensor_grid[5][3]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,g,h,d,e,f)+=L(a,b,c,d,e,f)*R(c,b,g,h)",
                        &H_6_legs_b, &H_6_legs_a, tensor_grid[6][3]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(g,h,c,d,e,f)+=L(a,b,c,d,e,f)*R(b,a,g,h)",
                        &H_6_legs_a, &H_6_legs_b, tensor_grid[7][3]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(g,b,c,d,e,f)+=L(a,b,c,d,e,f)*R(a,g)",
                        &H_6_legs_b, &H_6_legs_a, tensor_grid[8][3]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Level 4
    tc = TensContraction("D(h,g,b,c,d,e,f)+=L(a,b,c,d,e,f)*R(g,a,h)",
                        &H_7_legs_a, &H_6_legs_b, tensor_grid[8][4]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,i,h,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,c,b,i)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[7][4]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,i,h,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,d,c,i)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[6][4]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,i,h,f,g)+=L(a,b,c,d,e,f,g)*R(h,e,d,i)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[5][4]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,d,i,h,g)+=L(a,b,c,d,e,f,g)*R(h,f,e,i)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[4][4]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,d,e,i,h)+=L(a,b,c,d,e,f,g)*R(h,g,f,i)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[3][4]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Level 5
    tc = TensContraction("D(h,b,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,a)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[8][5]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(i,h,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,b,a,i)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[7][5]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,i,h,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,c,b,i)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[6][5]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,i,h,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,d,c,i)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[5][5]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,i,h,f,g)+=L(a,b,c,d,e,f,g)*R(h,e,d,i)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[4][5]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,d,i,h,g)+=L(a,b,c,d,e,f,g)*R(h,f,e,i)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[3][5]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Level 6
    tc = TensContraction("D(h,b,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,a)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[7][6]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(i,h,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,b,a,i)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[6][6]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,i,h,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,c,b,i)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[5][6]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,i,h,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,d,c,i)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[4][6]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,i,h,f,g)+=L(a,b,c,d,e,f,g)*R(h,e,d,i)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[3][6]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Level 7
    tc = TensContraction("D(h,b,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,a)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[6][7]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(i,h,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,b,a,i)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[5][7]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,i,h,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,c,b,i)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[4][7]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,i,h,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,d,c,i)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[3][7]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Corner
    tc = TensContraction("D(h,b,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,a)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[5][8]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(i,h,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,b,a,i)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[4][8]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(h,b,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,a)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[4][9]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(h,b,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,a)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[3][9]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(h,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,c,b,a)",
                        &H_5_legs_a, &H_7_legs_a, tensor_grid[3][8]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Row 2
    tc = TensContraction("D(f,b,c,d,e)+=L(a,b,c,d,e)*R(f,a)",
                        &H_5_legs_b, &H_5_legs_a, tensor_grid[2][8]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(f,g,c,d,e)+=L(a,b,c,d,e)*R(f,g,b,a)",
                        &H_5_legs_a, &H_5_legs_b, tensor_grid[2][7]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,f,g,d,e)+=L(a,b,c,d,e)*R(f,g,c,b)",
                        &H_5_legs_b, &H_5_legs_a, tensor_grid[2][6]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,f,g,e)+=L(a,b,c,d,e)*R(f,g,d,c)",
                        &H_5_legs_a, &H_5_legs_b, tensor_grid[2][5]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,f)+=L(a,b,c,d,e)*R(f,e,d)",
                        &H_4_legs_c, &H_5_legs_a, tensor_grid[2][4]);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Finished contracting A

    // Now contract C
    for (int c=0; c<num_Cs; ++c)
    {
      // Build Cs
      for (int t=0; t<qubits_A.size(); ++t)
      {
        int i = qubits_A[t][0]; int j = qubits_A[t][1];
        string delta_gate = (final_conf_A[c][t]=='0')?"delta_0":"delta_1";
        vector<s_type> gate_vector(_GATES_DATA.at(delta_gate));
        s_type * gate_data = (s_type *) malloc(gate_vector.size()*sizeof(s_type));
        errc = talsh::pinHostMemory(gate_data,gate_vector.size()*sizeof(s_type));
        assert(errc==0);
        for (size_t p=0; p<gate_vector.size(); ++p)
          *(gate_data+p) = gate_vector[p];
        vector<int> dims_delta({DIM});
        talsh::Tensor delta(dims_delta, gate_data);
        delete_this.push_back(gate_data);
        volumes.push_back(gate_vector.size());
        gate_data = nullptr;
        if (t==3 || t==4)
        {
          tc = TensContraction("D(b,c,d,e)+=L(a)*R(a,b,c,d,e)",
                              Cs[t], &delta, tensor_grid[i][j]);
          errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
          assert(tc.sync(DEV_HOST,0));
        }
        else
        {
          tc = TensContraction("D(b,c)+=L(a)*R(a,b,c)",
                              Cs[t], &delta, tensor_grid[i][j]);
          errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
          assert(tc.sync(DEV_HOST,0));
        }
      }
      // Contract C
      tc = TensContraction("D(a,c)+=L(a,b)*R(b,c)",
                          &H_2_legs_a, Cs[0], Cs[1]);
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(c,d,e,b)+=L(a,b)*R(a,c,d,e)",
                          &H_4_legs_a, &H_2_legs_a, Cs[3]);
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(e,b,c,d)+=L(a,b,c,d)*R(e,a)",
                          &H_4_legs_b, &H_4_legs_a, Cs[2]);
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,e,f)+=L(a,b,c,d)*R(d,c,e,f)",
                          &H_4_legs_a, &H_4_legs_b, Cs[4]);
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,c,e)+=L(a,b,c,d)*R(d,e)",
                          &H_4_legs_b, &H_4_legs_a, Cs[5]);
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      // Contract A and C
      tc = TensContraction("D()+=L(a,b,c,d)*R(d,c,b,a)",
                          &S, &H_4_legs_c, &H_4_legs_b);
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      // Store result data_S
      amplitudes.push_back(*data_S);
    }


    // Delete all pointers
    for (auto & v : delete_this)
    {
      errc = talsh::unpinHostMemory(v); assert(errc == 0);
      delete [] v; v = nullptr;
    }
    data_c05 = nullptr;
    data_c06 = nullptr;
    data_c14 = nullptr;
    data_c15 = nullptr;
    data_c16 = nullptr;
    data_c17 = nullptr;
    data_2_legs_a = nullptr;
    data_2_legs_b = nullptr;
    data_3_legs_a = nullptr;
    data_3_legs_b = nullptr;
    data_4_legs_a = nullptr;
    data_4_legs_b = nullptr;
    data_4_legs_c = nullptr;
    data_5_legs_a = nullptr;
    data_5_legs_b = nullptr;
    data_6_legs_a = nullptr;
    data_6_legs_b = nullptr;
    data_7_legs_a = nullptr;
    data_7_legs_b = nullptr;
    data_S = nullptr;
    // Delete all tensors
    for (int i=0; i<I; ++i) for (int j=0; j<J; ++j)
    {
      if (tensor_grid[i][j]!=nullptr)
      {
        delete tensor_grid[i][j];
        tensor_grid[i][j] = nullptr;
      }
    }
    for (int t=0; t<Cs.size(); ++t)
    {
      delete Cs[t];
      Cs[t] = nullptr;
    }
  }
  // Shut down TALSH
  talsh::shutdown();

  // Printing output
  for (int c=0; c<num_Cs; ++c)
  {
    cout << initial_conf << " ";
    cout << final_conf_B << " ";
    cout << final_conf_A[c] << " ";
    cout << real(amplitudes[c]) << " "
         << imag(amplitudes[c]) << " ";
    cout << "\n";
  }


  return 0;
} 
