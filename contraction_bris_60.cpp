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

#include "contraction_bris_60.h"
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
    {2,0},{2,1},{2,2}                                    ,{2,9},{2,10},{2,11},
    {3,0},{3,1}                                                ,{3,10},{3,11},
    {4,0}                                                             ,{4,11},
    {5,0}                                                             ,{5,11},
    {6,0},{6,1}                                                ,{6,10},{6,11},
    {7,0},{7,1},{7,2}                                    ,{7,9},{7,10},{7,11},
    {8,0},{8,1},{8,2},{8,3}                        ,{8,8},{8,9},{8,10},{8,11},
    {9,0},{9,1},{9,2},{9,3},{9,4}            ,{9,7},{9,8},{9,9},{9,10},{9,11},
    {10,0},{10,1},{10,2},{10,3},{10,4},{10,5},{10,6},{10,7},{10,8},{10,9},{10,10},{10,11}
    }); 
  qubits_A = vector<vector<int>>({ {0,6},
                                   {1,6},{1,7},
                                   {2,6},{2,7},{2,8} });

   
  int errc;

  // Allocate talsh::Tensors involved in the contraction
  vector<int> dims_1(1, super_dim);
  vector<int> dims_2(2, super_dim);
  vector<int> dims_3(3, super_dim);
  vector<int> dims_4(4, super_dim);
  vector<int> dims_5(5, super_dim);
  vector<int> dims_6(6, super_dim);
  vector<int> dims_7(7, super_dim);

  // First, tensors for region C. Done by hand right now. Change in future.
  vector<int> dims_C0(2,super_dim); dims_C0[0] = 1;
  Cs.push_back(shared_ptr<talsh::Tensor>(
                      new talsh::Tensor(dims_C0, s_type(0.0))));
  vector<int> dims_C1(4,super_dim); dims_C1[1] = 1;
  Cs.push_back(shared_ptr<talsh::Tensor>(
                      new talsh::Tensor(dims_C1, s_type(0.0))));
  vector<int> dims_C2(2,super_dim);
  Cs.push_back(shared_ptr<talsh::Tensor>(
                      new talsh::Tensor(dims_C2, s_type(0.0))));
  vector<int> dims_C3(4,super_dim); dims_C3[1] = 1;
  Cs.push_back(shared_ptr<talsh::Tensor>(
                      new talsh::Tensor(dims_C3, s_type(0.0))));
  vector<int> dims_C4(4,super_dim);
  Cs.push_back(shared_ptr<talsh::Tensor>(
                      new talsh::Tensor(dims_C4, s_type(0.0))));
  vector<int> dims_C5(2,super_dim);
  Cs.push_back(shared_ptr<talsh::Tensor>(
                      new talsh::Tensor(dims_C5, s_type(0.0))));

  // Second, helper tensors.
  H_2_legs_a =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_2, s_type(0.0)));
  H_2_legs_b =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_2, s_type(0.0)));
  H_4_legs_a =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_4, s_type(0.0)));
  H_4_legs_b =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_4, s_type(0.0)));

  vector<int> dims_2_pp(4,super_dim); dims_2_pp[0] = 1; dims_2_pp[1] = 1;
  H_2_legs_pp_a =
        shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_2_pp, s_type(0.0)));
  H_2_legs_pp_b =
        shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_2_pp, s_type(0.0)));
  vector<int> dims_3_ppp(6,super_dim);
              dims_3_ppp[0] = 1; dims_3_ppp[1] = 1; dims_3_ppp[2] = 1;
  H_3_legs_ppp_a =
        shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_3_ppp, s_type(0.0)));
  H_3_legs_ppp_b =
        shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_3_ppp, s_type(0.0)));
  H_3_legs_ppp_c =
        shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_3_ppp, s_type(0.0)));
  vector<int> dims_5_ppp(8,super_dim);
              dims_5_ppp[0] = 1; dims_5_ppp[1] = 1; dims_5_ppp[2] = 1;
  H_5_legs_ppp_a =
        shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_5_ppp, s_type(0.0)));
  H_5_legs_ppp_b =
        shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_5_ppp, s_type(0.0)));
  vector<int> dims_7_ppp(10,super_dim);
              dims_7_ppp[0] = 1; dims_7_ppp[1] = 1; dims_7_ppp[2] = 1;
  H_7_legs_ppp_a =
        shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_7_ppp, s_type(0.0)));
  H_7_legs_ppp_b =
        shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_7_ppp, s_type(0.0)));
  

  // Corner tensors
  A_pair = shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_4, s_type(0.0)));
  B_pair = shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_4, s_type(0.0)));
  A_corn = shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_4, s_type(0.0)));
  B_corn = shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_4, s_type(0.0)));
  D_corn = shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_4, s_type(0.0)));

  // Also, need tensors to hold slices.
  // Dimensions of the ones that are open are moved one position back.
  // First position should be the open leg, by convention.
  vector<int> dims_S05(2, super_dim); dims_S05[1] = 1;
  vector<int> dims_S06(3, super_dim); dims_S06[1] = 1; dims_S06[0] = DIM;
  vector<int> dims_S15(4, super_dim); dims_S15[3] = 1;
  vector<int> dims_S16(5, super_dim); dims_S16[2] = 1; dims_S16[0] = DIM;
  vector<int> dims_S25(4, super_dim); dims_S25[3] = 1;
  vector<int> dims_S26(5, super_dim); dims_S26[2] = 1; dims_S26[0] = DIM;

  S05 = shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_S05, s_type(0.0)));
  S06 = shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_S06, s_type(0.0)));
  S15 = shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_S15, s_type(0.0)));
  S16 = shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_S16, s_type(0.0)));
  S25 = shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_S25, s_type(0.0)));
  S26 = shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_S26, s_type(0.0)));

  // Finally, scalar S
  S = shared_ptr<talsh::Tensor>(new talsh::Tensor({}, s_type(0.0)));

  // Cut configurations (hard coded for the moment)
  slice_confs = vector<vector<int>>({ {0, 0, 0},
                                      {0, 0, 1},
                                      {0, 0, 2},
                                      {0, 0, 3},
                                      {0, 0, 4},
                                      {0, 0, 5},
                                      {0, 0, 6},
                                      {0, 0, 7},
                                      {0, 0, 8},
                                      {0, 0, 9},
                                      {0, 0, 10},
                                      {0, 0, 12},
                                      {0, 0, 14},
                                      {0, 0, 15},
                                      {0, 1, 0},
                                      {0, 1, 1},
                                      {0, 1, 2},
                                      {0, 1, 3},
                                      {0, 1, 4},
                                      {0, 1, 5}  });

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
  // Fix this so that it is only read once from file! Leave legs up open
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

  // Build corners
  // A_corn
  TensContraction tc("D(b,c)+=L(a,b)*R(a,c)", H_2_legs_a.get(),
                      tensor_grid[4][1].get(), tensor_grid[5][1].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(b,d,e,c)+=L(a,b)*R(c,a,d,e)", H_4_legs_a.get(),
                      H_2_legs_a.get(), tensor_grid[4][2].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,b,c,e)+=L(a,b,c,d)*R(d,e)", H_4_legs_b.get(),
                      H_4_legs_a.get(), tensor_grid[3][2].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(e,f,c,d)+=L(a,b,c,d)*R(b,a,e,f)", H_4_legs_a.get(),
                      H_4_legs_b.get(), tensor_grid[5][2].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(e,b,c,d)+=L(a,b,c,d)*R(a,e)", A_corn.get(),
                      H_4_legs_a.get(), tensor_grid[6][2].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));

  // B_corn
  tc = TensContraction("D(a,c)+=L(a,b)*R(c,b)", H_2_legs_a.get(),
                      tensor_grid[4][10].get(), tensor_grid[5][10].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(c,d,e,b)+=L(a,b)*R(c,d,e,a)", H_4_legs_a.get(),
                      H_2_legs_a.get(), tensor_grid[4][9].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(e,b,c,d)+=L(a,b,c,d)*R(e,a)", H_4_legs_b.get(),
                      H_4_legs_a.get(), tensor_grid[3][9].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,b,e,f)+=L(a,b,c,d)*R(c,e,f,d)", H_4_legs_a.get(),
                      H_4_legs_b.get(), tensor_grid[5][9].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,b,c,e)+=L(a,b,c,d)*R(d,e)", B_corn.get(),
                      H_4_legs_a.get(), tensor_grid[6][9].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));

  // D_corn
  tc = TensContraction("D(a,c)+=L(a,b)*R(c,b)", H_2_legs_a.get(),
                      tensor_grid[9][5].get(), tensor_grid[9][6].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,e,c,d)+=L(a,b)*R(c,d,b,e)", H_4_legs_a.get(),
                      H_2_legs_a.get(), tensor_grid[8][5].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(a,b,c,e)+=L(a,b,c,d)*R(d,e)", H_4_legs_b.get(),
                      H_4_legs_a.get(), tensor_grid[8][4].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(f,e,c,d)+=L(a,b,c,d)*R(e,b,a,f)", H_4_legs_a.get(),
                      H_4_legs_b.get(), tensor_grid[8][6].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  tc = TensContraction("D(e,b,c,d)+=L(a,b,c,d)*R(e,a)", D_corn.get(),
                      H_4_legs_a.get(), tensor_grid[8][7].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  t1 = high_resolution_clock::now();
  // Corners built!
  // Pair of tensors by A_corn
  tc = TensContraction(
          "D(a,b,e,d)+=L(a,b,c,d)*R(c,e)",
          A_pair.get(), tensor_grid[6][3].get(), tensor_grid[7][3].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  // Pair of tensors by B_corn
  tc = TensContraction(
          "D(c,d,b,e)+=L(a,b)*R(c,d,a,e)",
          B_pair.get(), tensor_grid[7][8].get(), tensor_grid[6][8].get());
  errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  // Pairs built!
  

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
}


vector<s_type> const & Contraction::get_amplitudes() const
{
  return amplitudes;
}

/////////////////////////// EXTERNAL FUNCTIONS ////////////////////////////////
// NONE
