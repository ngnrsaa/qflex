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

#include <vector>
#include <string>
#include <cmath>
#include <complex>
#include <ctime>
#include <chrono>
#include <cassert>

#include <iostream>
#include <fstream>
#include <sstream>

#include "contraction.h"
#include "read_circuit_full_talsh.h"

// Time
#include <ctime>
#include <chrono>
using namespace chrono;


///////////////////////////// CLASS FUNCTIONS /////////////////////////////////

Contraction::Contraction(string input_string, int _num_args, int _num_amps)
{
  stringstream ss(input_string);
  ss >> I;
  ss >> J;
  ss >> K;
  ss >> fidelity;
  ss >> filename;

  num_args = _num_args;
  num_amps = _num_amps;

  super_dim = (size_t)pow(DIM,K);
  num_Cs = num_args-2;
  amplitudes.resize(num_amps*num_Cs);
  int num_qubits = I*J;
  // List of qubits to remove (off) and qubits in A.
  qubits_off = vector<vector<int>>({
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
  qubits_A = vector<vector<int>>({       {0,5},{0,6},
                                  {1,4},{1,5},{1,6},{1,7} });

   
  int errc;

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
  vector<shared_ptr<talsh::Tensor>> Cs; 
  Cs.push_back(shared_ptr<talsh::Tensor>(
                      new talsh::Tensor(dims_2, s_type(0.0))));
  Cs.push_back(shared_ptr<talsh::Tensor>(
                      new talsh::Tensor(dims_2, s_type(0.0))));
  Cs.push_back(shared_ptr<talsh::Tensor>(
                      new talsh::Tensor(dims_2, s_type(0.0))));
  Cs.push_back(shared_ptr<talsh::Tensor>(
                      new talsh::Tensor(dims_4, s_type(0.0))));
  Cs.push_back(shared_ptr<talsh::Tensor>(
                      new talsh::Tensor(dims_4, s_type(0.0))));
  Cs.push_back(shared_ptr<talsh::Tensor>(
                      new talsh::Tensor(dims_2, s_type(0.0))));

  // Second, helper tensors.
  H_2_legs_a =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_2, s_type(0.0)));
  H_2_legs_b =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_2, s_type(0.0)));
  H_3_legs_a =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_3, s_type(0.0)));
  H_3_legs_b =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_3, s_type(0.0)));
  H_4_legs_a =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_4, s_type(0.0)));
  H_4_legs_b =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_4, s_type(0.0)));
  H_4_legs_c =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_4, s_type(0.0)));
  H_5_legs_a =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_5, s_type(0.0)));
  H_5_legs_b =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_5, s_type(0.0)));
  H_6_legs_a =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_6, s_type(0.0)));
  H_6_legs_b =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_6, s_type(0.0)));
  H_7_legs_a =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_7, s_type(0.0)));
  H_7_legs_b =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_7, s_type(0.0)));

  // Finally, scalar S
  S = shared_ptr<talsh::Tensor>(new talsh::Tensor({}, s_type(0.0)));

  cout << "Created contraction with input string " << input_string << endl;
}


// This function defines all the details of the contraction for each circuit
// size.
void Contraction::contract(string input_string)
{
  int errc;


  string initial_conf, final_conf_B;
  vector<string> final_conf_A(num_Cs, "");

  // Read string with a stringstream
  stringstream ss(input_string);
  ss >> initial_conf;
  ss >> final_conf_B;
  string conf_A;
  for (int s=0; s<num_Cs; ++s)
  {
    ss >> conf_A;
    final_conf_A[s] = string(conf_A);
  }

  // Declaring and then filling 2D grid of tensors.
  vector<vector<shared_ptr<talsh::Tensor>>> tensor_grid(I);
  for (int i=0; i<I; ++i)
  {
    tensor_grid[i] = vector<shared_ptr<talsh::Tensor>>(J);
  }
  // Fix this so that it is only read once from file! Leave legs up open
  google_circuit_file_to_grid_of_tensors(filename, I, J, initial_conf,
               final_conf_B, qubits_A, qubits_off, tensor_grid);

  // Start actual contraction.
  // It's working, but I haven't respected the ordering criterion!
  // Level 1 (starting with 1)
  TensContraction tc("D(c,d,b)+=L(a,b)*R(a,c,d)", H_3_legs_a.get(),
                      tensor_grid[4][1].get(), tensor_grid[5][1].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(d,b,c)+=L(a,b,c)*R(a,d)", H_3_legs_b.get(),
                        H_3_legs_a.get(), tensor_grid[6][1].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));

  // Level 2 (starting with 4)
  tc = TensContraction("D(a,b,e,f,d)+=L(a,b,c)*R(d,c,e,f)", H_5_legs_a.get(),
                        H_3_legs_b.get(), tensor_grid[4][2].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,b,c,d,f)+=L(a,b,c,d,e)*R(e,f)", H_5_legs_b.get(),
                        H_5_legs_a.get(), tensor_grid[3][2].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,f,g,d,e)+=L(a,b,c,d,e)*R(c,b,f,g)",
                  H_5_legs_a.get(), H_5_legs_b.get(), tensor_grid[5][2].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(f,g,c,d,e)+=L(a,b,c,d,e)*R(b,a,f,g)",
                  H_5_legs_b.get(), H_5_legs_a.get(), tensor_grid[6][2].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(f,b,c,d,e)+=L(a,b,c,d,e)*R(a,f)", H_5_legs_a.get(),
                        H_5_legs_b.get(), tensor_grid[7][2].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));

  // Level 3 (starting with 9)
  tc = TensContraction("D(a,b,c,d,f,g)+=L(a,b,c,d,e)*R(e,f,g)",
                  H_6_legs_a.get(), H_5_legs_a.get(), tensor_grid[3][3].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,b,c,g,h,f)+=L(a,b,c,d,e,f)*R(e,d,g,h)",
                  H_6_legs_b.get(), H_6_legs_a.get(), tensor_grid[4][3].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,b,g,h,e,f)+=L(a,b,c,d,e,f)*R(d,c,g,h)",
                  H_6_legs_a.get(), H_6_legs_b.get(), tensor_grid[5][3].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,g,h,d,e,f)+=L(a,b,c,d,e,f)*R(c,b,g,h)",
                  H_6_legs_b.get(), H_6_legs_a.get(), tensor_grid[6][3].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(g,h,c,d,e,f)+=L(a,b,c,d,e,f)*R(b,a,g,h)",
                  H_6_legs_a.get(), H_6_legs_b.get(), tensor_grid[7][3].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(g,b,c,d,e,f)+=L(a,b,c,d,e,f)*R(a,g)",
                  H_6_legs_b.get(), H_6_legs_a.get(), tensor_grid[8][3].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));


  // Start using pipelines for GPU contractions
  vector<TensContraction> conts_gpu;

  // First, emplace all contractions

  // Level 4 (starting with 15)
  conts_gpu.emplace_back(
    TensContraction("D(h,g,b,c,d,e,f)+=L(a,b,c,d,e,f)*R(g,a,h)",
    H_7_legs_a.get(), H_6_legs_b.get(), tensor_grid[8][4].get()));
  conts_gpu.emplace_back(
    TensContraction("D(a,i,h,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,c,b,i)",
    H_7_legs_b.get(), H_7_legs_a.get(), tensor_grid[7][4].get()));
  conts_gpu.emplace_back(
    TensContraction("D(a,b,i,h,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,d,c,i)",
    H_7_legs_a.get(), H_7_legs_b.get(), tensor_grid[6][4].get()));
  conts_gpu.emplace_back(
    TensContraction("D(a,b,c,i,h,f,g)+=L(a,b,c,d,e,f,g)*R(h,e,d,i)",
    H_7_legs_b.get(), H_7_legs_a.get(), tensor_grid[5][4].get()));
  conts_gpu.emplace_back(
    TensContraction("D(a,b,c,d,i,h,g)+=L(a,b,c,d,e,f,g)*R(h,f,e,i)",
    H_7_legs_a.get(), H_7_legs_b.get(), tensor_grid[4][4].get()));
  conts_gpu.emplace_back(
    TensContraction("D(a,b,c,d,e,i,h)+=L(a,b,c,d,e,f,g)*R(h,g,f,i)",
    H_7_legs_b.get(), H_7_legs_a.get(), tensor_grid[3][4].get()));


  // Level 5 (starting with 21)
  conts_gpu.emplace_back(
    TensContraction("D(h,b,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,a)",
    H_7_legs_a.get(), H_7_legs_b.get(), tensor_grid[8][5].get()));
  conts_gpu.emplace_back(
    TensContraction("D(i,h,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,b,a,i)",
    H_7_legs_b.get(), H_7_legs_a.get(), tensor_grid[7][5].get()));
  conts_gpu.emplace_back(
    TensContraction("D(a,i,h,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,c,b,i)",
    H_7_legs_a.get(), H_7_legs_b.get(), tensor_grid[6][5].get()));
  conts_gpu.emplace_back(
    TensContraction("D(a,b,i,h,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,d,c,i)",
    H_7_legs_b.get(), H_7_legs_a.get(), tensor_grid[5][5].get()));
  conts_gpu.emplace_back(
    TensContraction("D(a,b,c,i,h,f,g)+=L(a,b,c,d,e,f,g)*R(h,e,d,i)",
    H_7_legs_a.get(), H_7_legs_b.get(), tensor_grid[4][5].get()));
  conts_gpu.emplace_back(
    TensContraction("D(a,b,c,d,i,h,g)+=L(a,b,c,d,e,f,g)*R(h,f,e,i)",
    H_7_legs_b.get(), H_7_legs_a.get(), tensor_grid[3][5].get()));

  // Level 6 (starting with 27)
  conts_gpu.emplace_back(
    TensContraction("D(h,b,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,a)",
    H_7_legs_a.get(), H_7_legs_b.get(), tensor_grid[7][6].get()));
  conts_gpu.emplace_back(
    TensContraction("D(i,h,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,b,a,i)",
    H_7_legs_b.get(), H_7_legs_a.get(), tensor_grid[6][6].get()));
  conts_gpu.emplace_back(
    TensContraction("D(a,i,h,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,c,b,i)",
    H_7_legs_a.get(), H_7_legs_b.get(), tensor_grid[5][6].get()));
  conts_gpu.emplace_back(
    TensContraction("D(a,b,i,h,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,d,c,i)",
    H_7_legs_b.get(), H_7_legs_a.get(), tensor_grid[4][6].get()));
  conts_gpu.emplace_back(
    TensContraction("D(a,b,c,i,h,f,g)+=L(a,b,c,d,e,f,g)*R(h,e,d,i)",
    H_7_legs_a.get(), H_7_legs_b.get(), tensor_grid[3][6].get()));

  // Level 7 (starting with 32)
  conts_gpu.emplace_back(
    TensContraction("D(h,b,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,a)",
    H_7_legs_b.get(), H_7_legs_a.get(), tensor_grid[6][7].get()));
  conts_gpu.emplace_back(
    TensContraction("D(i,h,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,b,a,i)",
    H_7_legs_a.get(), H_7_legs_b.get(), tensor_grid[5][7].get()));
  conts_gpu.emplace_back(
    TensContraction("D(a,i,h,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,c,b,i)",
    H_7_legs_b.get(), H_7_legs_a.get(), tensor_grid[4][7].get()));
  conts_gpu.emplace_back(
    TensContraction("D(a,b,i,h,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,d,c,i)",
    H_7_legs_a.get(), H_7_legs_b.get(), tensor_grid[3][7].get()));

  // Corner (starting with 36)
  conts_gpu.emplace_back(
    TensContraction("D(h,b,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,a)",
    H_7_legs_b.get(), H_7_legs_a.get(), tensor_grid[5][8].get()));
  conts_gpu.emplace_back(
    TensContraction("D(i,h,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,b,a,i)",
    H_7_legs_a.get(), H_7_legs_b.get(), tensor_grid[4][8].get()));
  conts_gpu.emplace_back(
    TensContraction("D(h,b,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,a)",
    H_7_legs_b.get(), H_7_legs_a.get(), tensor_grid[4][9].get()));
  conts_gpu.emplace_back(
    TensContraction("D(h,b,c,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,a)",
    H_7_legs_a.get(), H_7_legs_b.get(), tensor_grid[3][9].get()));
  conts_gpu.emplace_back(
    TensContraction("D(h,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,c,b,a)",
    H_5_legs_a.get(), H_7_legs_a.get(), tensor_grid[3][8].get()));


  // Second, actually contract
  for (auto & tcg : conts_gpu)
  {
    errc = tcg.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tcg.sync(DEV_NVIDIA_GPU,0));
  }
  cout << "Done GPU, corner!" << endl;

  // Corner (starting with 36) (only last one, because of sync on DEV_HOST)
  tc = TensContraction("D(h,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,c,b,a)",
                  H_5_legs_a.get(), H_7_legs_a.get(), tensor_grid[3][8].get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));


  // Row 2 (starting with 41)
  tc = TensContraction("D(f,b,c,d,e)+=L(a,b,c,d,e)*R(f,a)",
                  H_5_legs_b.get(), H_5_legs_a.get(), tensor_grid[2][8].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(f,g,c,d,e)+=L(a,b,c,d,e)*R(f,g,b,a)",
                  H_5_legs_a.get(), H_5_legs_b.get(), tensor_grid[2][7].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,f,g,d,e)+=L(a,b,c,d,e)*R(f,g,c,b)",
                  H_5_legs_b.get(), H_5_legs_a.get(), tensor_grid[2][6].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,b,f,g,e)+=L(a,b,c,d,e)*R(f,g,d,c)",
                  H_5_legs_a.get(), H_5_legs_b.get(), tensor_grid[2][5].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,b,c,f)+=L(a,b,c,d,e)*R(f,e,d)",
                  H_4_legs_c.get(), H_5_legs_a.get(), tensor_grid[2][4].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  // Finished contracting B


  /*
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
                        H_2_legs_a.get(), Cs[0].get(), Cs[1].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(c,d,e,b)+=L(a,b)*R(a,c,d,e)",
                        H_4_legs_a.get(), H_2_legs_a.get(), Cs[3].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(e,b,c,d)+=L(a,b,c,d)*R(e,a)",
                        H_4_legs_b.get(), H_4_legs_a.get(), Cs[2].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,e,f)+=L(a,b,c,d)*R(d,c,e,f)",
                        H_4_legs_a.get(), H_4_legs_b.get(), Cs[4].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,e)+=L(a,b,c,d)*R(d,e)",
                        H_4_legs_b.get(), H_4_legs_a.get(), Cs[5].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Contract A and C
    tc = TensContraction("D()+=L(a,b,c,d)*R(d,c,b,a)",
                        S.get(), H_4_legs_c.get(), H_4_legs_b.get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Store result data_S
    s_type const * ptr_S;
    S->getDataAccessHostConst(&ptr_S);
    amplitudes[c] = ptr_S[0];
  }
  */
}


vector<s_type> const & Contraction::get_amplitudes() const
{
  return amplitudes;
}

/////////////////////////// EXTERNAL FUNCTIONS ////////////////////////////////
