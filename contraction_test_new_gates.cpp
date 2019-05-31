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

#include "contraction_test_new_gates.h"
#include "read_circuit_full_talsh.h"

// Time
#include <ctime>
#include <chrono>
using namespace chrono;


///////////////////////////// CLASS FUNCTIONS /////////////////////////////////

Contraction::Contraction(string input_string, int _num_args, int _num_amps)
{
  cout << "Calling constructor" << endl << flush;

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
  amplitudes.resize(num_amps*max(1,num_Cs));
  int num_qubits = I*J;
  // List of qubits to remove (off) and qubits in A.
  qubits_off = vector<vector<int>>({ }); 
  qubits_A = vector<vector<int>>({ {0,1} });

   
  int errc;

  // Allocate talsh::Tensors involved in the contraction
  vector<int> dims_1(1, super_dim);
  vector<int> dims_2(2, super_dim);
  vector<int> dims_3(3, super_dim);


  // Second, helper tensors.
  H_2_legs_a =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_2, s_type(0.0)));
  H_2_legs_b =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_2, s_type(0.0)));
  H_3_legs_a =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_3, s_type(0.0)));
  H_3_legs_b =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_3, s_type(0.0)));
  

  // Finally, scalar S
  S = shared_ptr<talsh::Tensor>(new talsh::Tensor({}, s_type(0.0)));


  // HARD CODED norm_factor
  norm_factor = s_type(1.2);
  global_norm_factor = s_type(1.);
  for (int q=0; q<I*J-qubits_off.size(); ++q)
    global_norm_factor *= norm_factor;

  //cout << "Created contraction with input string " << input_string << endl;
}


// This function might soon build all tensors in region C (if strings are
// repeated across batches).
void Contraction::load_circuit(string input_string)
{
  cout << "Calling load_circuit" << endl << flush;
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


  // If add renormalization, do it here!
  renormalize_circuit(I, J, qubits_off, open_tensor_grid, norm_factor);


  // Initialize Cs here, because for general gates we don't know the
  // dimension of each C tensor until here.
  // First: get dimensions of each C. 
  for (int t=0; t<qubits_A.size(); ++t)
  {
    int i = qubits_A[t][0]; int j = qubits_A[t][1];
    vector<int> dims_T = dims_from_tensor(open_tensor_grid[i][j].get());
    vector<int> dims_C(dims_T.size()-1);
    for (int p=0; p<dims_C.size(); ++p)
    {
      dims_C[p] = dims_T[p+1];
    }
    Cs.push_back(shared_ptr<talsh::Tensor>(
                        new talsh::Tensor(dims_C, s_type(0.0))));
  }


  t1 = high_resolution_clock::now();
  span = duration_cast<duration<double>>(t1 - t0);
  //cout << "Time spent loading circuit is " << span.count() << endl;
}


// This function defines all the details of the contraction for each circuit
// size.
// Allocation of tensors might go out too.
void Contraction::contract(string input_string)
{
  cout << "Calling contract" << endl << flush;
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
  //cout << "Time spent closing circuit is " << span.count() << endl;
  

  t0 = high_resolution_clock::now();
  // Start actual contraction.
  // I'll follow L = inherited tensor, R = new tensor.
  bool done;
  
  // For the reshaping of the slicing tensors.
  vector<int> dims_2(2, super_dim);


  // Contraction
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
      TensContraction tc("D(b,c)+=L(a)*R(a,b,c)",
                          Cs[t].get(), &delta, tensor_grid[i][j].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
    }
    // Contract A and C
    TensContraction tc("D()+=L(a)*R(a)", S.get(),
                        tensor_grid[0][0].get(), Cs[0].get());
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Store result data_S
    s_type const * ptr_S;
    S->getDataAccessHostConst(&ptr_S);
    amplitudes[c] += ptr_S[0] / global_norm_factor;
  }
  t1 = high_resolution_clock::now();
  span = duration_cast<duration<double>>(t1 - t0);

}


vector<s_type> const & Contraction::get_amplitudes() const
{
  cout << "Calling get_amplitudes" << endl << flush;
  return amplitudes;
}

double Contraction::get_time_largest_contraction() const
{
  cout << "Calling get_time_largest_contraction" << endl << flush;
  return time_largest_contraction;
}

/////////////////////////// EXTERNAL FUNCTIONS ////////////////////////////////
// NONE