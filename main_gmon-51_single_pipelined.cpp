/*

  Copyright © 2019, United States Government, as represented by the Administrator
  of the National Aeronautics and Space Administration. All rights reserved.
  
  The Flexible Quantum Circuit Simulator (qFlex)  platform is licensed under the
  Apache License, Version 2.0 (the "License"); you may not use this file except in
  compliance with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0. 
  
  Unless required by applicable law or agreed to in writing, software distributed
  under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
  CONDITIONS OF ANY KIND, either express or implied. See the License for the
  specific language governing permissions and limitations under the License.

*/


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

#include <omp.h>

#include "read_circuit_full_talsh.h"

#include "talshxx.hpp"
#include "talsh_wrapper.h"

using namespace std;
using namespace chrono;


// Input: I J K fidelity filename initial_conf (optional) final_conf (optional)
int main(int argc, char **argv) {

  // Set precision for the printed floats.
  cout.precision(12);

  // Timing variables.
  high_resolution_clock::time_point t_output_0, t_output_1;
  t_output_0 = high_resolution_clock::now();
  high_resolution_clock::time_point t0, t1;
  duration<double> time_span;

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


  int num_threads = 4;
  omp_set_num_threads(num_threads);

  // Vector of amplitudes with size num_Cs
  vector<vector<s_type>> amplitudes(num_threads);

  // Initialize TALSH
  unsigned long size(12000000000);
  talsh::initialize(&size);
  {

    // Declaring and then filling 2D grid of tensors.
    vector<vector<vector<shared_ptr<talsh::Tensor>>>> tensor_grid(num_threads);
    for (int t=0; t<num_threads; ++t)
    {
      tensor_grid[t] = vector<vector<shared_ptr<talsh::Tensor>>>(I);
      for (int i=0; i<I; ++i)
      {
        tensor_grid[t][i] = vector<shared_ptr<talsh::Tensor>>(J);
      }
      google_circuit_file_to_grid_of_tensors(filename, I, J, initial_conf,
                   final_conf_B, qubits_A, qubits_off, tensor_grid[t]);
    }


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
    vector<vector<shared_ptr<talsh::Tensor>>> Cs(num_threads); 
    for (int t=0; t<num_threads; ++t)
    {
      Cs[t].push_back(shared_ptr<talsh::Tensor>(
                          new talsh::Tensor(dims_2, s_type(0.0))));
      Cs[t].push_back(shared_ptr<talsh::Tensor>(
                          new talsh::Tensor(dims_2, s_type(0.0))));
      Cs[t].push_back(shared_ptr<talsh::Tensor>(
                          new talsh::Tensor(dims_2, s_type(0.0))));
      Cs[t].push_back(shared_ptr<talsh::Tensor>(
                          new talsh::Tensor(dims_4, s_type(0.0))));
      Cs[t].push_back(shared_ptr<talsh::Tensor>(
                          new talsh::Tensor(dims_4, s_type(0.0))));
      Cs[t].push_back(shared_ptr<talsh::Tensor>(
                          new talsh::Tensor(dims_2, s_type(0.0))));
    }

    // Second, helper tensors.
    vector<talsh::Tensor> H_2_legs_a;
    vector<talsh::Tensor> H_2_legs_b;
    vector<talsh::Tensor> H_3_legs_a;
    vector<talsh::Tensor> H_3_legs_b;
    vector<talsh::Tensor> H_4_legs_a;
    vector<talsh::Tensor> H_4_legs_b;
    vector<talsh::Tensor> H_4_legs_c;
    vector<talsh::Tensor> H_5_legs_a;
    vector<talsh::Tensor> H_5_legs_b;
    vector<talsh::Tensor> H_6_legs_a;
    vector<talsh::Tensor> H_6_legs_b;
    vector<talsh::Tensor> H_7_legs_a;
    vector<talsh::Tensor> H_7_legs_b;
    for (int t=0; t<num_threads; ++t)
    {
      H_2_legs_a.push_back(talsh::Tensor(dims_2, s_type(0.0)));
      H_2_legs_b.push_back(talsh::Tensor(dims_2, s_type(0.0)));
      H_3_legs_a.push_back(talsh::Tensor(dims_3, s_type(0.0)));
      H_3_legs_b.push_back(talsh::Tensor(dims_3, s_type(0.0)));
      H_4_legs_a.push_back(talsh::Tensor(dims_4, s_type(0.0)));
      H_4_legs_b.push_back(talsh::Tensor(dims_4, s_type(0.0)));
      H_4_legs_c.push_back(talsh::Tensor(dims_4, s_type(0.0)));
      H_5_legs_a.push_back(talsh::Tensor(dims_5, s_type(0.0)));
      H_5_legs_b.push_back(talsh::Tensor(dims_5, s_type(0.0)));
      H_6_legs_a.push_back(talsh::Tensor(dims_6, s_type(0.0)));
      H_6_legs_b.push_back(talsh::Tensor(dims_6, s_type(0.0)));
      H_7_legs_a.push_back(talsh::Tensor(dims_7, s_type(0.0)));
      H_7_legs_b.push_back(talsh::Tensor(dims_7, s_type(0.0)));
    }

    // Finally, scalar S
    vector<talsh::Tensor> S;
    for (int t=0; t<num_threads; ++t)
    {
      S.push_back(talsh::Tensor({}, s_type(0.0)));
    }

    for (int t=0; t<num_threads; ++t)
    {
      int errc;
      //int t = omp_get_thread_num();
      cout << t << endl;
      // Start contracting.
      // It's working, but I haven't respected the ordering criterion!
      // Level 1 (starting with 1)
      TensContraction tc("D(c,d,b)+=L(a,b)*R(a,c,d)", &H_3_legs_a[t],
                          tensor_grid[t][4][1].get(), tensor_grid[t][5][1].get());
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(d,b,c)+=L(a,b,c)*R(a,d)", &H_3_legs_b[t],
                            &H_3_legs_a[t], tensor_grid[t][6][1].get());
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));

      // Level 2 (starting with 4)
      tc = TensContraction("D(a,b,e,f,d)+=L(a,b,c)*R(d,c,e,f)", &H_5_legs_a[t],
                            &H_3_legs_b[t], tensor_grid[t][4][2].get());
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,c,d,f)+=L(a,b,c,d,e)*R(e,f)", &H_5_legs_b[t],
                            &H_5_legs_a[t], tensor_grid[t][3][2].get());
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,f,g,d,e)+=L(a,b,c,d,e)*R(c,b,f,g)", &H_5_legs_a[t],
                            &H_5_legs_b[t], tensor_grid[t][5][2].get());
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(f,g,c,d,e)+=L(a,b,c,d,e)*R(b,a,f,g)", &H_5_legs_b[t],
                            &H_5_legs_a[t], tensor_grid[t][6][2].get());
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(f,b,c,d,e)+=L(a,b,c,d,e)*R(a,f)", &H_5_legs_a[t],
                            &H_5_legs_b[t], tensor_grid[t][7][2].get());
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));

      // Level 3 (starting with 9)
      tc = TensContraction("D(a,b,c,d,f,g)+=L(a,b,c,d,e)*R(e,f,g)", &H_6_legs_a[t],
                            &H_5_legs_a[t], tensor_grid[t][3][3].get());
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,c,g,h,f)+=L(a,b,c,d,e,f)*R(e,d,g,h)",
                            &H_6_legs_b[t], &H_6_legs_a[t], tensor_grid[t][4][3].get());
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,g,h,e,f)+=L(a,b,c,d,e,f)*R(d,c,g,h)",
                            &H_6_legs_a[t], &H_6_legs_b[t], tensor_grid[t][5][3].get());
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,g,h,d,e,f)+=L(a,b,c,d,e,f)*R(c,b,g,h)",
                            &H_6_legs_b[t], &H_6_legs_a[t], tensor_grid[t][6][3].get());
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(g,h,c,d,e,f)+=L(a,b,c,d,e,f)*R(b,a,g,h)",
                            &H_6_legs_a[t], &H_6_legs_b[t], tensor_grid[t][7][3].get());
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(g,b,c,d,e,f)+=L(a,b,c,d,e,f)*R(a,g)",
                            &H_6_legs_b[t], &H_6_legs_a[t], tensor_grid[t][8][3].get());
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
    }


    // Start using pipelines for GPU contractions
    vector<vector<TensContraction>> conts_gpu;
    int c = 0;

    // First, push all contractions
    // Level 4 (starting with 15)
    t0 = high_resolution_clock::now();
    conts_gpu.push_back(vector<TensContraction>());
    for (int t=0; t<num_threads; ++t)
    {
      conts_gpu[c].emplace_back(
        TensContraction("D(h,g,b,c,d,e,f)+=L(a,b,c,d,e,f)*R(g,a,h)",
        &(H_7_legs_a[t]), &(H_6_legs_b[t]), tensor_grid[t][8][4].get()));
    }
    ++c;
    conts_gpu.push_back(vector<TensContraction>());
    for (int t=0; t<num_threads; ++t)
    {
      conts_gpu[c].emplace_back(
        TensContraction("D(a,i,h,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,c,b,i)",
        &H_7_legs_b[t], &H_7_legs_a[t], tensor_grid[t][7][4].get()));
    }
    ++c;
    conts_gpu.push_back(vector<TensContraction>());
    for (int t=0; t<num_threads; ++t)
    {
      conts_gpu[c].emplace_back(
        TensContraction("D(a,b,i,h,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,d,c,i)",
        &H_7_legs_a[t], &H_7_legs_b[t], tensor_grid[t][6][4].get()));
    }
    ++c;
    conts_gpu.push_back(vector<TensContraction>());
    for (int t=0; t<num_threads; ++t)
    {
      conts_gpu[c].emplace_back(
        TensContraction("D(a,b,c,i,h,f,g)+=L(a,b,c,d,e,f,g)*R(h,e,d,i)",
        &H_7_legs_b[t], &H_7_legs_a[t], tensor_grid[t][5][4].get()));
    }
    ++c;
    conts_gpu.push_back(vector<TensContraction>());
    for (int t=0; t<num_threads; ++t)
    {
      conts_gpu[c].emplace_back(
        TensContraction("D(a,b,c,d,i,h,g)+=L(a,b,c,d,e,f,g)*R(h,f,e,i)",
        &H_7_legs_a[t], &H_7_legs_b[t], tensor_grid[t][4][4].get()));
    }
    ++c;
    conts_gpu.push_back(vector<TensContraction>());
    for (int t=0; t<num_threads; ++t)
    {
      conts_gpu[c].emplace_back(
        TensContraction("D(a,b,c,d,e,i,h)+=L(a,b,c,d,e,f,g)*R(h,g,f,i)",
        &H_7_legs_b[t], &H_7_legs_a[t], tensor_grid[t][3][4].get()));
    }
    ++c;

    /*
    // Level 5 (starting with 21)
    tc = TensContraction("D(h,b,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,a)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[8][5].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(i,h,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,b,a,i)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[7][5].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,i,h,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,c,b,i)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[6][5].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,i,h,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,d,c,i)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[5][5].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,i,h,f,g)+=L(a,b,c,d,e,f,g)*R(h,e,d,i)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[4][5].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,d,i,h,g)+=L(a,b,c,d,e,f,g)*R(h,f,e,i)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[3][5].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));

    // Level 6 (starting with 27)
    tc = TensContraction("D(h,b,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,a)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[7][6].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(i,h,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,b,a,i)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[6][6].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,i,h,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,c,b,i)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[5][6].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,i,h,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,d,c,i)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[4][6].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,i,h,f,g)+=L(a,b,c,d,e,f,g)*R(h,e,d,i)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[3][6].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));

    // Level 7 (starting with 32)
    tc = TensContraction("D(h,b,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,a)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[6][7].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(i,h,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,b,a,i)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[5][7].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,i,h,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,c,b,i)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[4][7].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,i,h,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,d,c,i)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[3][7].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));

    // Corner (starting with 36)
    tc = TensContraction("D(h,b,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,a)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[5][8].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(i,h,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,b,a,i)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[4][8].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(h,b,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,a)",
                        &H_7_legs_b, &H_7_legs_a, tensor_grid[4][9].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(h,b,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,a)",
                        &H_7_legs_a, &H_7_legs_b, tensor_grid[3][9].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(h,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,c,b,a)",
                        &H_5_legs_a, &H_7_legs_a, tensor_grid[3][8].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    t1 = high_resolution_clock::now();
    time_span = duration_cast<duration<double>>(t1 - t0);
    cout << "Time spent on large contractions: "
         << time_span.count()
         << "s\n\n";
    */

    // Second, actually contract
    for (int d=0; d<c; ++d)
    {
      unsigned int remains = conts_gpu[d].size();
      auto its = conts_gpu[d].begin();
      // Need to be independent contractions. Otherwise, it could scehdule
      // two of them with (C = A*B and D = C*E), and there could be problems.
      while (remains > 0)
      {
        while (its != conts_gpu[d].end())
        {
          int errc = its->execute(DEV_NVIDIA_GPU,0); // Tries schedule execution
          assert(errc==TALSH_SUCCESS || errc==TRY_LATER); // Checks correctness
          if (errc==TRY_LATER) break; // If not scheduled, then do sth different
          ++its; // If scheduled, go for the next one
        }
        for (auto itc=conts_gpu[d].begin(); itc!=its; ++itc) // Check scheduled cs
        {
          int sts;
          if (itc->ready(&sts,DEV_HOST,0)) // If done (and move memory)
          {
            if (sts==TALSH_TASK_COMPLETED) // If this status, then --remains
            {
              --remains;
              cout << " Left " << remains << " tensor contractions" << endl;
            }
          }
        }
      }
    }
    /*
    // Row 2 (starting with 41)
    tc = TensContraction("D(f,b,c,d,e)+=L(a,b,c,d,e)*R(f,a)",
                        &H_5_legs_b, &H_5_legs_a, tensor_grid[2][8].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(f,g,c,d,e)+=L(a,b,c,d,e)*R(f,g,b,a)",
                        &H_5_legs_a, &H_5_legs_b, tensor_grid[2][7].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,f,g,d,e)+=L(a,b,c,d,e)*R(f,g,c,b)",
                        &H_5_legs_b, &H_5_legs_a, tensor_grid[2][6].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,f,g,e)+=L(a,b,c,d,e)*R(f,g,d,c)",
                        &H_5_legs_a, &H_5_legs_b, tensor_grid[2][5].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,f)+=L(a,b,c,d,e)*R(f,e,d)",
                        &H_4_legs_c, &H_5_legs_a, tensor_grid[2][4].get());
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
        vector<int> dims_delta({DIM});
        talsh::Tensor delta(dims_delta, gate_vector);
        if (t==3 || t==4)
        {
          tc = TensContraction("D(b,c,d,e)+=L(a)*R(a,b,c,d,e)",
                              Cs[t].get(), &delta, tensor_grid[i][j].get());
          errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
          assert(tc.sync(DEV_HOST,0));
        }
        else
        {
          tc = TensContraction("D(b,c)+=L(a)*R(a,b,c)",
                              Cs[t].get(), &delta, tensor_grid[i][j].get());
          errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
          assert(tc.sync(DEV_HOST,0));
        }
      }
      // Contract C
      tc = TensContraction("D(a,c)+=L(a,b)*R(b,c)",
                          &H_2_legs_a, Cs[0].get(), Cs[1].get());
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(c,d,e,b)+=L(a,b)*R(a,c,d,e)",
                          &H_4_legs_a, &H_2_legs_a, Cs[3].get());
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(e,b,c,d)+=L(a,b,c,d)*R(e,a)",
                          &H_4_legs_b, &H_4_legs_a, Cs[2].get());
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,e,f)+=L(a,b,c,d)*R(d,c,e,f)",
                          &H_4_legs_a, &H_4_legs_b, Cs[4].get());
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,c,e)+=L(a,b,c,d)*R(d,e)",
                          &H_4_legs_b, &H_4_legs_a, Cs[5].get());
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      // Contract A and C
      tc = TensContraction("D()+=L(a,b,c,d)*R(d,c,b,a)",
                          &S, &H_4_legs_c, &H_4_legs_b);
      errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      // Store result data_S
      s_type const * ptr_S;
      S.getDataAccessHostConst(&ptr_S);
      amplitudes[omp_get_thread_num()].push_back(ptr_S[0]);
      cout << ptr_S[0] << endl;
    }
    }
    */

  }
  // Shut down TALSH
  talsh::shutdown();

  /*
  // Printing output
  for (int t=0; t<num_threads; ++t)
  {
    for (int c=0; c<num_Cs; ++c)
    {
      cout << initial_conf << " ";
      cout << final_conf_B << " ";
      cout << final_conf_A[c] << " ";
      cout << real(amplitudes[t][c]) << " "
           << imag(amplitudes[t][c]) << " ";
      cout << "\n";
    }
  }
  */

  return 0;
} 
