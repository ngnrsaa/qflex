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

#include "contraction_7x7x40.h"
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
  qubits_off = vector<vector<int>>({ }); 
  qubits_A = vector<vector<int>>({ {5,0},{5,1},{5,2},
                                   {6,0},{6,1},{6,2} });

   
  int errc;

  // Allocate talsh::Tensors involved in the contraction
  vector<int> dims_1(1, super_dim);
  vector<int> dims_2(2, super_dim);
  vector<int> dims_3(3, super_dim);
  vector<int> dims_4(4, super_dim);
  vector<int> dims_5(5, super_dim);
  vector<int> dims_6(6, super_dim);

  // First, tensors for region C. Done by hand right now. Change in future.
  Cs.push_back(shared_ptr<talsh::Tensor>(
                      new talsh::Tensor(dims_3, s_type(0.0))));
  Cs.push_back(shared_ptr<talsh::Tensor>(
                      new talsh::Tensor(dims_4, s_type(0.0))));
  Cs.push_back(shared_ptr<talsh::Tensor>(
                      new talsh::Tensor(dims_4, s_type(0.0))));
  Cs.push_back(shared_ptr<talsh::Tensor>(
                      new talsh::Tensor(dims_2, s_type(0.0))));
  Cs.push_back(shared_ptr<talsh::Tensor>(
                      new talsh::Tensor(dims_3, s_type(0.0))));
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
  H_5_legs_a =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_5, s_type(0.0)));
  H_5_legs_b =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_5, s_type(0.0)));
  H_6_legs_a =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_6, s_type(0.0)));
  H_6_legs_b =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_6, s_type(0.0)));
  H_6_legs_c =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_6, s_type(0.0)));

  
  // Region tensors
  AB = shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_6, s_type(0.0)));
  // A doesn't need a tensor, because 6_a and 6_b are not used anywhere else.

  // Also, need tensors to hold slices.
  // Dimensions of the ones that are open are moved one position back.
  // First position should be the open leg, by convention.
  vector<int> dims_S26(3, super_dim); dims_S26[2] = 1;
  S26 = shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_S26, s_type(0.0)));
  vector<int> dims_S36(3, super_dim); dims_S36[0] = 1;
  S36 = shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_S36, s_type(0.0)));
  vector<int> dims_S62(4, super_dim); dims_S62[0] = DIM; dims_S62[3] = 1;
  S62 = shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_S62, s_type(0.0)));
  vector<int> dims_S63(3, super_dim); dims_S63[1] = 1;
  vector<int> dims_P62(3, super_dim); dims_P62[0] = DIM;
  S63 = shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_S63, s_type(0.0)));

  // Finally, scalar S
  S = shared_ptr<talsh::Tensor>(new talsh::Tensor({}, s_type(0.0)));

  // Cut configurations (hard coded for the moment)
  slice_confs = vector<vector<int>>({ {0, 0},
                                      {0, 1}  });

  cout << "Created contraction with input string " << input_string << endl;
}


// This function might soon build all tensors in region C (if strings are
// repeated across batches).
void Contraction::load_circuit(string input_string)
{
  // Timing variables.
  high_resolution_clock::time_point t0, t1;
  duration<double> span; 
  int errc;

  // Input and output bit-strings
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
  t0 = high_resolution_clock::now();
  open_tensor_grid = vector<vector<shared_ptr<talsh::Tensor>>>(I);
  for (int i=0; i<I; ++i)
  {
    open_tensor_grid[i] = vector<shared_ptr<talsh::Tensor>>(J);
  }
  google_circuit_file_to_open_grid_of_tensors(filename, I, J, initial_conf,
                                              qubits_off, open_tensor_grid);
  t1 = high_resolution_clock::now();
  span = duration_cast<duration<double>>(t1 - t0);
  cout << "Time spent loading circuit is " << span.count() << endl;
}


// This function defines all the details of the contraction for each circuit
// size.
// Allocation of tensors might go out too.
void Contraction::contract(string input_string)
{
  // Timing variables.
  high_resolution_clock::time_point t0, t1;
  duration<double> span;

  int errc;

  // Set amplitudes to 0, in order to start adding paths.
  for (auto & v : amplitudes)
    v = s_type({0.0,0.0});

  // Input and output bit-strings
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

  t0 = high_resolution_clock::now();
  vector<vector<shared_ptr<talsh::Tensor>>> tensor_grid(I);
  for (int i=0; i<I; ++i)
  {
    tensor_grid[i] = vector<shared_ptr<talsh::Tensor>>(J);
  }
  close_circuit(I, J, final_conf_B, qubits_A, qubits_off,
                open_tensor_grid, tensor_grid);
  t1 = high_resolution_clock::now();
  span = duration_cast<duration<double>>(t1 - t0);
  cout << "Time spent closing circuit is " << span.count() << endl;
  

  t0 = high_resolution_clock::now();
  // Start actual contraction.
  // I'll follow L = inherited tensor, R = new tensor.
  bool done;
  
  // For the reshaping of the slicing tensors.
  vector<int> dims_2(2, super_dim);

  // A
  // No XL
  TensContraction tc("D(a,c,d)+=L(a,b)*R(b,c,d)", H_3_legs_a.get(),
                      tensor_grid[0][0].get(), tensor_grid[0][1].get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,b,d,e)+=L(a,b,c)*R(c,d,e)", H_4_legs_a.get(),
                      H_3_legs_a.get(), tensor_grid[0][2].get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(e,f,b,c,d)+=L(a,b,c,d)*R(a,e,f)", H_5_legs_a.get(),
                      H_4_legs_a.get(), tensor_grid[1][0].get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,f,g,d,e)+=L(a,b,c,d,e)*R(c,b,f,g)",
                H_5_legs_b.get(), H_5_legs_a.get(), tensor_grid[1][1].get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,b,f,g,e)+=L(a,b,c,d,e)*R(d,c,f,g)",
                H_5_legs_a.get(), H_5_legs_b.get(), tensor_grid[1][2].get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  // XL
  errc = H_6_legs_a->contractAccumulateXL(nullptr,
          "D(f,g,b,c,d,e)+=L(a,b,c,d,e)*R(a,f,g)", *H_5_legs_a,
          *tensor_grid[2][0], DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}));
  done = H_6_legs_a->sync();
  errc = H_6_legs_b->contractAccumulateXL(nullptr,
          "D(a,g,h,d,e,f)+=L(a,b,c,d,e,f)*R(c,b,g,h)", *H_6_legs_a,
          *tensor_grid[2][1], DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}));
  done = H_6_legs_b->sync();
  errc = H_6_legs_c->contractAccumulateXL(nullptr,
          "D(a,b,g,h,e,f)+=L(a,b,c,d,e,f)*R(d,c,g,h)", *H_6_legs_b,
          *tensor_grid[2][1], DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}));
  done = H_6_legs_c->sync();

  // B
  // Slice tensor_grid[2][6] => S26 => P26.
  vector<int> dims_S26(3, super_dim); dims_S26[2] = 1;
  S26->reshape(dims_S26);
  int i0 = 0;
  talsh::TensorTask task_hl;
  errc = tensor_grid[2][6]->extractSlice(
                      &task_hl,*S26,vector<int>{0,0,i0},DEV_HOST,0);
  done = tensor_grid[2][6]->sync();
  task_hl.clean();
  S26->reshape(dims_2);
  // Contract
  // No XL
  tc = TensContraction("D(a,c,d)+=L(a,b)*R(b,c,d)", H_3_legs_a.get(),
                      tensor_grid[0][6].get(), tensor_grid[1][6].get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,b,d)+=L(a,b,c)*R(c,d)", H_3_legs_b.get(),
                      H_3_legs_a.get(), S26.get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(d,e,b,c)+=L(a,b,c)*R(d,e,a)", H_4_legs_a.get(),
                      H_3_legs_b.get(), tensor_grid[0][5].get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,e,f,d)+=L(a,b,c,d)*R(b,e,f,c)", H_4_legs_b.get(),
                      H_4_legs_a.get(), tensor_grid[1][5].get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,b,e,f)+=L(a,b,c,d)*R(c,e,f,d)", H_4_legs_a.get(),
                      H_4_legs_b.get(), tensor_grid[2][5].get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(e,f,b,c,d)+=L(a,b,c,d)*R(e,f,a)", H_5_legs_a.get(),
                      H_4_legs_a.get(), tensor_grid[0][4].get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,f,g,d,e)+=L(a,b,c,d,e)*R(b,f,g,c)",
              H_5_legs_b.get(), H_5_legs_a.get(), tensor_grid[1][4].get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,b,f,g,e)+=L(a,b,c,d,e)*R(c,f,g,d)",
              H_5_legs_a.get(), H_5_legs_b.get(), tensor_grid[2][4].get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  // XL
  errc = H_6_legs_a->contractAccumulateXL(nullptr,
          "D(f,g,b,c,d,e)+=L(a,b,c,d,e)*R(f,g,a)", *H_5_legs_a,
          *tensor_grid[0][3], DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}));
  done = H_6_legs_a->sync();
  errc = H_6_legs_b->contractAccumulateXL(nullptr,
          "D(a,g,h,d,e,f)+=L(a,b,c,d,e,f)*R(b,g,h,c)", *H_6_legs_a,
          *tensor_grid[1][3], DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}));
  done = H_6_legs_b->sync();
  errc = H_6_legs_a->contractAccumulateXL(nullptr,
          "D(a,b,g,h,e,f)+=L(a,b,c,d,e,f)*R(c,g,h,d)", *H_6_legs_b,
          *tensor_grid[2][3], DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}));
  done = H_6_legs_a->sync();

  // Contract AB: H_6_legs_c and H_6_legs_b onto AB.
  errc = AB->contractAccumulateXL(nullptr,
          "D(a,b,c,g,h,i)+=L(a,b,c,d,e,f)*R(f,e,d,g,h,i)", *H_6_legs_c,
          *H_6_legs_b, DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}));
  done = H_6_legs_a->sync();

  // E
  // Slice tensor_grid[3][6] => S36
  vector<int> dims_S36(3, super_dim); dims_S36[0] = 1;
  S36->reshape(dims_S36);
  errc = tensor_grid[3][6]->extractSlice(
                      &task_hl,*S36,vector<int>{i0,0,0},DEV_HOST,0);
  done = tensor_grid[3][6]->sync();
  task_hl.clean();
  S36->reshape(dims_2);
  // Contract
  // No XL
  tc = TensContraction("D(c,d,b)+=L(a,b)*R(c,d,a)", H_3_legs_a.get(),
                      tensor_grid[6][6].get(), tensor_grid[5][6].get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(d,e,b,c)+=L(a,b,c)*R(d,e,a)", H_4_legs_a.get(),
                      H_3_legs_a.get(), tensor_grid[4][6].get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(e,b,c,d)+=L(a,b,c,d)*R(e,a)", H_4_legs_b.get(),
                      H_4_legs_a.get(), S36.get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,b,c,e,f)+=L(a,b,c,d)*R(e,f,d)", H_5_legs_a.get(),
                      H_4_legs_b.get(), tensor_grid[6][5].get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,b,f,g,e)+=L(a,b,c,d,e)*R(f,g,d,c)",
                H_5_legs_b.get(), H_5_legs_a.get(), tensor_grid[5][5].get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,f,g,d,e)+=L(a,b,c,d,e)*R(f,g,c,b)",
                H_5_legs_a.get(), H_5_legs_b.get(), tensor_grid[4][5].get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(f,g,c,d,e)+=L(a,b,c,d,e)*R(f,g,b,a)",
                H_5_legs_b.get(), H_5_legs_a.get(), tensor_grid[3][5].get());
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  // XL
  errc = H_6_legs_a->contractAccumulateXL(nullptr,
          "D(a,b,c,d,f,g)+=L(a,b,c,d,e)*R(f,g,e)", *H_5_legs_b,
          *tensor_grid[6][4], DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}));
  done = H_6_legs_a->sync();
  errc = H_6_legs_b->contractAccumulateXL(nullptr,
          "D(a,b,c,g,h,f)+=L(a,b,c,d,e,f)*R(g,h,e,d)", *H_6_legs_a,
          *tensor_grid[5][4], DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}));
  done = H_6_legs_b->sync();
  errc = H_6_legs_a->contractAccumulateXL(nullptr,
          "D(a,b,g,h,e,f)+=L(a,b,c,d,e,f)*R(g,h,d,c)", *H_6_legs_b,
          *tensor_grid[4][4], DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}));
  done = H_6_legs_a->sync();
  errc = H_6_legs_b->contractAccumulateXL(nullptr,
          "D(a,g,h,d,e,f)+=L(a,b,c,d,e,f)*R(g,h,c,b)", *H_6_legs_a,
          *tensor_grid[4][4], DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}));
  done = H_6_legs_b->sync();
  // H_6_legs_b is pD!

  // Begin inner loop
  for (int i1=0; i1<6; ++i1)
  {
    // Slice tensor_grid[6][2] => S62
    vector<int> dims_S62(4, super_dim); dims_S62[0] = DIM; dims_S62[3] = 1;
    S62->reshape(dims_S62);
    errc = tensor_grid[6][2]->extractSlice(
                        &task_hl,*S62,vector<int>{0,0,0,i1},DEV_HOST,0);
    done = tensor_grid[6][2]->sync();
    task_hl.clean();
    vector<int> dims_P62(3, super_dim); dims_S62[0] = DIM;
    S62->reshape(dims_P62);
    // Slice tensor_grid[6][3] => S63
    vector<int> dims_S63(3, super_dim); dims_S63[1] = 1;
    S63->reshape(dims_S63);
    errc = tensor_grid[6][3]->extractSlice(
                        &task_hl,*S63,vector<int>{0,i1,0},DEV_HOST,0);
    done = tensor_grid[6][3]->sync();
    task_hl.clean();
    S63->reshape(dims_2);
    
    // Finish contracting D
    tc = TensContraction("D(c,d,b,e)+=L(a,b)*R(c,d,a,e)",
                  H_4_legs_a.get(), S63.get(), tensor_grid[5][3].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    errc = H_6_legs_a->contractAccumulateXL(nullptr,
          "D(a,b,c,d,g,h)+=L(a,b,c,d,e,f)*R(g,h,f,e)", *H_6_legs_a,
          *H_4_legs_a, DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}));
    done = H_6_legs_a->sync();
    errc = H_6_legs_b->contractAccumulateXL(nullptr,
          "D(a,b,c,g,h,f)+=L(a,b,c,d,e,f)*R(g,h,e,d)", *H_6_legs_a,
          *tensor_grid[4][3], DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}));
    done = H_6_legs_b->sync();
    errc = H_6_legs_a->contractAccumulateXL(nullptr,
          "D(a,b,g,h,e,f)+=L(a,b,c,d,e,f)*R(g,h,d,c)", *H_6_legs_b,
          *tensor_grid[3][3], DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}));
    done = H_6_legs_a->sync();

    // Contract ABD: AB and H_6_legs_a onto H_6_legs_c
    errc = H_6_legs_c->contractAccumulateXL(nullptr,
          "D(a,b,c,g,h,i)+=L(a,b,c,d,e,f)*R(f,e,d,g,h,i)", *AB,
          *H_6_legs_a, DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}));
    done = H_6_legs_c->sync();

    // Start contracting ABCD
    errc = H_6_legs_a->contractAccumulateXL(nullptr,
          "D(a,b,g,h,e,f)+=L(a,b,c,d,e,f)*R(c,g,h,d)", *H_6_legs_c,
          *tensor_grid[3][2], DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}));
    done = H_6_legs_a->sync();
    errc = H_6_legs_b->contractAccumulateXL(nullptr,
          "D(a,g,h,d,e,f)+=L(a,b,c,d,e,f)*R(b,g,h,c)", *H_6_legs_a,
          *tensor_grid[3][1], DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}));
    done = H_6_legs_b->sync();
    errc = H_5_legs_b->contractAccumulateXL(nullptr,
          "D(g,c,d,e,f)+=L(a,b,c,d,e,f)*R(a,g,b)", *H_6_legs_a,
          *tensor_grid[3][0], DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}));
    done = H_5_legs_b->sync();
    // No XL
    tc = TensContraction("D(a,b,f,g,e)+=L(a,b,c,d,e)*R(c,f,g,d)",
                  H_5_legs_a.get(), H_5_legs_b.get(), tensor_grid[4][2].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,f,g,d,e)+=L(a,b,c,d,e)*R(b,f,g,c)",
                  H_5_legs_b.get(), H_5_legs_a.get(), tensor_grid[4][1].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(f,c,d,e)+=L(a,b,c,d,e)*R(a,f,b)",
                  H_4_legs_a.get(), H_5_legs_b.get(), tensor_grid[4][0].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));

    /*
    t0 = high_resolution_clock::now();
    // Tensors live in the GPU and stay there!
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
        if (t==0)
        {
          tc = TensContraction("D(b,c)+=L(a)*R(a,b,c)",
                              Cs[t].get(), &delta, S06.get());
        } else if (t==1) {
          tc = TensContraction("D(b,c,d,e)+=L(a)*R(a,b,c,d,e)",
                              Cs[t].get(), &delta, S16.get());
        } else if (t==2) {
          tc = TensContraction("D(b,c)+=L(a)*R(a,b,c)",
                              Cs[t].get(), &delta, tensor_grid[i][j].get());
        } else if (t==3) {
          tc = TensContraction("D(b,c,d,e)+=L(a)*R(a,b,c,d,e)",
                              Cs[t].get(), &delta, S26.get());
        } else if (t==4) {
          tc = TensContraction("D(b,c,d,e)+=L(a)*R(a,b,c,d,e)",
                              Cs[t].get(), &delta, tensor_grid[i][j].get());
        } else if (t==5) {
          tc = TensContraction("D(b,c)+=L(a)*R(a,b,c)",
                              Cs[t].get(), &delta, tensor_grid[i][j].get());
        }
        errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
        assert(tc.sync(DEV_HOST,0));
      }
      // Contract C
      tc = TensContraction("D(a,c,d,e)+=L(a,b)*R(b,c,d,e)",
                      H_2_legs_pp_a.get(), Cs[0].get(), Cs[1].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,c,e)+=L(a,b,c,d)*R(d,e)",
                      H_2_legs_pp_b.get(), H_2_legs_pp_a.get(), Cs[2].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,e,f,g,d)+=L(a,b,c,d)*R(c,e,f,g)",
                      H_3_legs_ppp_a.get(), H_2_legs_pp_b.get(), Cs[3].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,c,d,g,h)+=L(a,b,c,d,e,f)*R(f,e,g,h)",
                      H_3_legs_ppp_b.get(), H_3_legs_ppp_a.get(), Cs[4].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,c,d,e,g)+=L(a,b,c,d,e,f)*R(f,g)",
                      H_3_legs_ppp_a.get(), H_3_legs_ppp_b.get(), Cs[5].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      // Contract A and C
      tc = TensContraction("D()+=L(a,b,c,d,e,f)*R(a,b,c,f,e,d)",
                          S.get(), H_3_legs_ppp_c.get(), H_3_legs_ppp_a.get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      // Store result data_S
      s_type const * ptr_S;
      S->getDataAccessHostConst(&ptr_S);
      amplitudes[c] += ptr_S[0];
    }
    t1 = high_resolution_clock::now();
    span = duration_cast<duration<double>>(t1 - t0);
    cout << "Time contracting all Cs is " << span.count() << endl;
    */
  }


  
  // UP TO HERE!
  /*

  // Begin iteration over slice_confs
  for (int s=0; s<slice_confs.size(); ++s)
  {
    int i0 = slice_confs[s][0];
    int i1 = slice_confs[s][0];
    int i2 = slice_confs[s][0];

    // Slice tensors cut
    t0 = high_resolution_clock::now();
    bool done;
    talsh::TensorTask task_hl;
    errc = tensor_grid[0][5]->extractSlice(
                        &task_hl,*S05,std::vector<int>{i0,0},DEV_HOST,0);
    done = tensor_grid[0][5]->sync();
    task_hl.clean();
    errc = tensor_grid[0][6]->extractSlice(
                        &task_hl,*S06,std::vector<int>{i0,0,0},DEV_HOST,0);
    done = tensor_grid[0][6]->sync();
    task_hl.clean();
    errc = tensor_grid[1][5]->extractSlice(
                        &task_hl,*S15,std::vector<int>{i1,0,0,0},DEV_HOST,0);
    done = tensor_grid[1][5]->sync();
    task_hl.clean();
    errc = tensor_grid[1][6]->extractSlice(
                        &task_hl,*S16,std::vector<int>{i1,0,0,0,0},DEV_HOST,0);
    done = tensor_grid[1][6]->sync();
    task_hl.clean();
    errc = tensor_grid[2][5]->extractSlice(
                        &task_hl,*S25,std::vector<int>{i2,0,0,0},DEV_HOST,0);
    done = tensor_grid[2][5]->sync();
    task_hl.clean();
    errc = tensor_grid[2][6]->extractSlice(
                        &task_hl,*S26,std::vector<int>{i2,0,0,0,0},DEV_HOST,0);
    done = tensor_grid[2][6]->sync();
    task_hl.clean();
    t1 = high_resolution_clock::now();
    span = duration_cast<duration<double>>(t1 - t0);
    cout << "Time spent slicing tensors is " << span.count() << endl;


    // Tensor A
    tc = TensContraction("D(b,e,c,d)+=L(a,b)*R(a,c,d,e)", H_2_legs_pp_a.get(),
                        S05.get(), S15.get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,e,d)+=L(a,b,c,d)*R(e,c)", H_2_legs_pp_b.get(),
                        H_2_legs_pp_a.get(), tensor_grid[1][4].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,g,c,e,f)+=L(a,b,c,d)*R(d,e,f,g)",
            H_3_legs_ppp_a.get(), H_2_legs_pp_b.get(), S25.get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,g,h,f)+=L(a,b,c,d,e,f)*R(d,g,h,e)",
            H_3_legs_ppp_b.get(), H_3_legs_ppp_a.get(), tensor_grid[2][4].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,g,e,f)+=L(a,b,c,d,e,f)*R(g,d)",
            H_3_legs_ppp_a.get(), H_3_legs_ppp_b.get(), tensor_grid[2][3].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,d,e,g,h,i)+=L(a,b,c,d,e,f)*R(f,g,h,i)",
            H_5_legs_ppp_a.get(), H_3_legs_ppp_a.get(), tensor_grid[3][5].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,d,i,j,g,h)+=L(a,b,c,d,e,f,g,h)*R(e,i,j,f)",
            H_5_legs_ppp_b.get(), H_5_legs_ppp_a.get(), tensor_grid[3][4].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,i,j,f,g,h)+=L(a,b,c,d,e,f,g,h)*R(d,i,j,e)",
            H_5_legs_ppp_a.get(), H_5_legs_ppp_b.get(), tensor_grid[3][3].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,d,e,f,i,j,k,h)+=L(a,b,c,d,e,f,g,h)*R(g,i,j,k)",
            H_7_legs_ppp_a.get(), H_5_legs_ppp_a.get(), tensor_grid[4][5].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
            "D(a,b,c,d,e,k,l,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(f,k,l,g)",
            H_7_legs_ppp_b.get(), H_7_legs_ppp_a.get(), tensor_grid[4][4].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
            "D(a,b,c,d,k,l,g,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(e,k,l,f)",
            H_7_legs_ppp_a.get(), H_7_legs_ppp_b.get(), tensor_grid[4][3].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // A_corn
    tc = TensContraction(
            "D(a,b,c,k,l,f,g,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(k,l,e,d)",
            H_7_legs_ppp_b.get(), H_7_legs_ppp_a.get(), A_corn.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Continue
    tc = TensContraction(
            "D(a,b,c,d,k,l,g,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(f,e,k,l)",
            H_7_legs_ppp_a.get(), H_7_legs_ppp_b.get(), tensor_grid[5][3].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
            "D(a,b,c,d,e,k,l,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(g,f,k,l)",
            H_7_legs_ppp_b.get(), H_7_legs_ppp_a.get(), tensor_grid[5][4].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
            "D(a,b,c,d,e,f,k,l,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(h,g,k,l)",
            H_7_legs_ppp_a.get(), H_7_legs_ppp_b.get(), tensor_grid[5][5].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Pair of tensors by A_corn
    tc = TensContraction(
            "D(a,b,c,k,l,f,g,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(e,d,k,l)",
            H_7_legs_ppp_b.get(), H_7_legs_ppp_a.get(), A_pair.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Continue
    tc = TensContraction(
            "D(a,b,c,d,k,l,g,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(f,e,k,l)",
            H_7_legs_ppp_a.get(), H_7_legs_ppp_b.get(), tensor_grid[6][4].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
            "D(a,b,c,d,e,k,l,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(g,f,k,l)",
            H_7_legs_ppp_b.get(), H_7_legs_ppp_a.get(), tensor_grid[6][5].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
            "D(a,b,c,k,l,f,g,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(e,d,k,l)",
            H_7_legs_ppp_a.get(), H_7_legs_ppp_b.get(), tensor_grid[7][4].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
            "D(a,b,c,d,k,l,g,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(f,e,k,l)",
            H_7_legs_ppp_b.get(), H_7_legs_ppp_a.get(), tensor_grid[7][5].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // D_corn
    tc = TensContraction(
            "D(a,b,c,k,l,f,g,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(k,l,e,d)",
            H_7_legs_ppp_a.get(), H_7_legs_ppp_b.get(), D_corn.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Continue
    tc = TensContraction(
            "D(a,b,c,d,l,k,g,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(k,f,e,l)",
            H_7_legs_ppp_b.get(), H_7_legs_ppp_a.get(), tensor_grid[7][6].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
            "D(a,b,c,l,k,f,g,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(k,e,d,l)",
            H_7_legs_ppp_a.get(), H_7_legs_ppp_b.get(), tensor_grid[7][7].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
            "D(a,b,c,d,e,l,k,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(k,g,f,l)",
            H_7_legs_ppp_b.get(), H_7_legs_ppp_a.get(), tensor_grid[6][6].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
            "D(a,b,c,d,l,k,g,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(k,f,e,l)",
            H_7_legs_ppp_a.get(), H_7_legs_ppp_b.get(), tensor_grid[6][7].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Pair of tensors by B_corn
    tc = TensContraction(
            "D(a,b,c,l,k,f,g,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(k,e,d,l)",
            H_7_legs_ppp_b.get(), H_7_legs_ppp_a.get(), B_pair.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Continue
    tc = TensContraction(
            "D(a,b,c,d,e,f,l,k,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(k,h,g,l)",
            H_7_legs_ppp_a.get(), H_7_legs_ppp_b.get(), tensor_grid[5][6].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
            "D(a,b,c,d,e,l,k,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(k,g,f,l)",
            H_7_legs_ppp_b.get(), H_7_legs_ppp_a.get(), tensor_grid[5][7].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
            "D(a,b,c,d,l,k,g,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(k,f,e,l)",
            H_7_legs_ppp_a.get(), H_7_legs_ppp_b.get(), tensor_grid[5][8].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // B_corn
    tc = TensContraction(
            "D(a,b,c,k,l,f,g,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(k,l,e,d)",
            H_7_legs_ppp_b.get(), H_7_legs_ppp_a.get(), B_corn.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Continue
    tc = TensContraction(
            "D(a,b,c,d,e,f,g,l,k,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(k,i,h,l)",
            H_7_legs_ppp_a.get(), H_7_legs_ppp_b.get(), tensor_grid[4][6].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
            "D(a,b,c,d,e,f,l,k,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(k,h,g,l)",
            H_7_legs_ppp_b.get(), H_7_legs_ppp_a.get(), tensor_grid[4][7].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
            "D(a,b,c,d,k,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(k,g,f,e)",
            H_5_legs_ppp_a.get(), H_7_legs_ppp_b.get(), tensor_grid[4][8].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
            "D(a,b,c,d,e,f,j,i)+=L(a,b,c,d,e,f,g,h)*R(i,h,g,j)",
            H_5_legs_ppp_b.get(), H_5_legs_ppp_a.get(), tensor_grid[3][6].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
            "D(a,b,c,d,e,j,i,h)+=L(a,b,c,d,e,f,g,h)*R(i,g,f,j)",
            H_5_legs_ppp_a.get(), H_5_legs_ppp_b.get(), tensor_grid[3][7].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
            "D(a,b,c,i,g,h)+=L(a,b,c,d,e,f,g,h)*R(i,f,e,d)",
            H_3_legs_ppp_c.get(), H_5_legs_ppp_a.get(), tensor_grid[3][8].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Done with A!



    t0 = high_resolution_clock::now();
    // Tensors live in the GPU and stay there!
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
        if (t==0)
        {
          tc = TensContraction("D(b,c)+=L(a)*R(a,b,c)",
                              Cs[t].get(), &delta, S06.get());
        } else if (t==1) {
          tc = TensContraction("D(b,c,d,e)+=L(a)*R(a,b,c,d,e)",
                              Cs[t].get(), &delta, S16.get());
        } else if (t==2) {
          tc = TensContraction("D(b,c)+=L(a)*R(a,b,c)",
                              Cs[t].get(), &delta, tensor_grid[i][j].get());
        } else if (t==3) {
          tc = TensContraction("D(b,c,d,e)+=L(a)*R(a,b,c,d,e)",
                              Cs[t].get(), &delta, S26.get());
        } else if (t==4) {
          tc = TensContraction("D(b,c,d,e)+=L(a)*R(a,b,c,d,e)",
                              Cs[t].get(), &delta, tensor_grid[i][j].get());
        } else if (t==5) {
          tc = TensContraction("D(b,c)+=L(a)*R(a,b,c)",
                              Cs[t].get(), &delta, tensor_grid[i][j].get());
        }
        errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
        assert(tc.sync(DEV_HOST,0));
      }
      // Contract C
      tc = TensContraction("D(a,c,d,e)+=L(a,b)*R(b,c,d,e)",
                      H_2_legs_pp_a.get(), Cs[0].get(), Cs[1].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,c,e)+=L(a,b,c,d)*R(d,e)",
                      H_2_legs_pp_b.get(), H_2_legs_pp_a.get(), Cs[2].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,e,f,g,d)+=L(a,b,c,d)*R(c,e,f,g)",
                      H_3_legs_ppp_a.get(), H_2_legs_pp_b.get(), Cs[3].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,c,d,g,h)+=L(a,b,c,d,e,f)*R(f,e,g,h)",
                      H_3_legs_ppp_b.get(), H_3_legs_ppp_a.get(), Cs[4].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,c,d,e,g)+=L(a,b,c,d,e,f)*R(f,g)",
                      H_3_legs_ppp_a.get(), H_3_legs_ppp_b.get(), Cs[5].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      // Contract A and C
      tc = TensContraction("D()+=L(a,b,c,d,e,f)*R(a,b,c,f,e,d)",
                          S.get(), H_3_legs_ppp_c.get(), H_3_legs_ppp_a.get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      // Store result data_S
      s_type const * ptr_S;
      S->getDataAccessHostConst(&ptr_S);
      amplitudes[c] += ptr_S[0];
    }
    t1 = high_resolution_clock::now();
    span = duration_cast<duration<double>>(t1 - t0);
    cout << "Time contracting all Cs is " << span.count() << endl;
  }
  */
}


vector<s_type> const & Contraction::get_amplitudes() const
{
  return amplitudes;
}

/////////////////////////// EXTERNAL FUNCTIONS ////////////////////////////////
// NONE
