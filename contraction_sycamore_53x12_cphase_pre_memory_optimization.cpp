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

#include "contraction_sycamore_53x12_cphase.h"
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
       {0,0},{0,1},{0,2},{0,3},{0,4},            {0,7},{0,8},{0,9},
       {1,0},{1,1},{1,2},{1,3},                        {1,8},{1,9},
       {2,0},{2,1},{2,2},{2,3},                              {2,9},
       {3,0},{3,1},                                                
       {4,0},                                                      
                                                             {5,9},
       {6,0},                                          {6,8},{6,9},
       {7,0},{7,1},                              {7,7},{7,8},{7,9},
       {8,0},{8,1},{8,2},                  {8,6},{8,7},{8,8},{8,9},
       {9,0},{9,1},{9,2},{9,3},      {9,5},{9,6},{9,7},{9,8},{9,9}
                                    }); 
  qubits_A = vector<vector<int>>({       {0,5},{0,6},
                                   {1,4},{1,5},{1,6},{1,7} });

   
  int errc;


  // Scalar S
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
  //cout << "Loading circuit..." << endl << flush;
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
  // Then storing dimensions.
  t0 = high_resolution_clock::now();
  open_tensor_grid = vector<vector<shared_ptr<talsh::Tensor>>>(I);
  for (int i=0; i<I; ++i)
  {
    open_tensor_grid[i] = vector<shared_ptr<talsh::Tensor>>(J);
  }
  google_circuit_file_to_open_grid_of_tensors(filename, I, J, initial_conf,
                                              qubits_off, open_tensor_grid);
  // Dimensions
  dims_grid = vector<vector<vector<int>>>(I);
  for (int i=0; i<I; ++i)
  {
    dims_grid[i] = vector<vector<int>>(J);
                                                 
    for (int j=0; j<J; ++j)
    {
      if (find(qubits_off.begin(),qubits_off.end(),
          vector<int>({i,j}))!=qubits_off.end())
      {
        continue;
      }
      vector<int> dims_T = dims_from_tensor(open_tensor_grid[i][j].get());
      dims_grid[i][j] = vector<int>(dims_T.size()-1);
      for (int p=0; p<dims_grid[i][j].size(); ++p)
      {
        dims_grid[i][j][p] = dims_T[p+1];
      }
    }
  }


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
  //cout << "Closing circuit ... " << endl << flush;
  close_circuit(I, J, final_conf_B, qubits_A, qubits_off,
                open_tensor_grid, tensor_grid);
  t1 = high_resolution_clock::now();
  span = duration_cast<duration<double>>(t1 - t0);
  //cout << "Time spent closing circuit is " << span.count() << endl;
  

  t0 = high_resolution_clock::now();
  // Start actual contraction.
  // I'll follow L = inherited tensor, R = new tensor.
  bool done;
  

  // Start contraction
  // Pre-cut
  // Corners ((5,0) and (9,4))
  unique_ptr<talsh::Tensor> corner_51(new talsh::Tensor({dims_grid[5][1][0],
                                                         dims_grid[5][1][2],
                                                         dims_grid[5][1][3]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
  errc = corner_51->contractAccumulate(nullptr,
        "D(b,c,d)+=L(a)*R(b,a,c,d)",
        *tensor_grid[5][0].get(),
        *tensor_grid[5][1].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = corner_51->sync();
  unique_ptr<talsh::Tensor> corner_84(new talsh::Tensor({dims_grid[8][4][0],
                                                         dims_grid[8][4][1],
                                                         dims_grid[8][4][3]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
  errc = corner_84->contractAccumulate(nullptr,
        "D(a,b,d)+=L(a,b,c,d)*R(c)",
        *tensor_grid[8][4].get(),
        *tensor_grid[9][4].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = corner_84->sync();

  // 4,1 & 4,2 ; 6,1 & 6,2 ; 6,6 & 6,7 ; 5,7 & 5,8 ; 2,7 & 2,8
  unique_ptr<talsh::Tensor> pair_41_42(new talsh::Tensor({dims_grid[4][1][0],
                                                         dims_grid[4][2][2],
                                                         dims_grid[4][2][3],
                                                         dims_grid[4][2][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
  errc = pair_41_42->contractAccumulate(nullptr,
        "D(a,d,e,c)+=L(a,b)*R(c,b,d,e)",
        *tensor_grid[4][1].get(),
        *tensor_grid[4][2].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = pair_41_42->sync();
  unique_ptr<talsh::Tensor> pair_61_62(new talsh::Tensor({dims_grid[6][2][2],
                                                         dims_grid[6][2][3],
                                                         dims_grid[6][2][0],
                                                         dims_grid[6][1][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
  errc = pair_61_62->contractAccumulate(nullptr,
        "D(d,e,c,a)+=L(a,b)*R(c,b,d,e)",
        *tensor_grid[6][1].get(),
        *tensor_grid[6][2].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = pair_61_62->sync();
  unique_ptr<talsh::Tensor> pair_66_67(new talsh::Tensor({dims_grid[6][7][0],
                                                         dims_grid[6][6][0],
                                                         dims_grid[6][6][1],
                                                         dims_grid[6][6][2]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
  errc = pair_66_67->contractAccumulate(nullptr,
        "D(e,a,b,c)+=L(a,b,c,d)*R(e,d)",
        *tensor_grid[6][6].get(),
        *tensor_grid[6][7].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = pair_66_67->sync();
  unique_ptr<talsh::Tensor> pair_57_58(new talsh::Tensor({dims_grid[5][8][0],
                                                         dims_grid[5][7][0],
                                                         dims_grid[5][7][1],
                                                         dims_grid[5][7][2]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
  errc = pair_57_58->contractAccumulate(nullptr,
        "D(e,a,b,c)+=L(a,b,c,d)*R(e,d)",
        *tensor_grid[5][7].get(),
        *tensor_grid[5][8].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = pair_57_58->sync();
  unique_ptr<talsh::Tensor> pair_27_28(new talsh::Tensor({dims_grid[2][7][0],
                                                         dims_grid[2][7][1],
                                                         dims_grid[2][7][2],
                                                         dims_grid[2][8][1]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
  errc = pair_27_28->contractAccumulate(nullptr,
        "D(a,b,c,e)+=L(a,b,c,d)*R(d,e)",
        *tensor_grid[2][7].get(),
        *tensor_grid[2][8].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = pair_27_28->sync();




  // Tensor C
  // C loop!
  // Temporarily calling C K.
  // First indexes are the open ones.
  unique_ptr<talsh::Tensor> ka(new talsh::Tensor({DIM, DIM,
                                    dims_grid[0][5][0],
                                    dims_grid[0][6][1]},
                                    talsh::TensorData<s_type>::kind,
                                    talsh_tens_no_init));
  errc = ka->contractAccumulate(nullptr,
        "D(o,p,a,c)+=L(o,a,b)*R(p,b,c)",
        *tensor_grid[0][5].get(),
        *tensor_grid[0][6].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = ka->sync();
  unique_ptr<talsh::Tensor> kb(new talsh::Tensor({DIM, DIM,
                                    dims_grid[1][4][0],
                                    dims_grid[1][5][2],
                                    dims_grid[1][5][3],
                                    dims_grid[1][5][0]},
                                    talsh::TensorData<s_type>::kind,
                                    talsh_tens_no_init));
  errc = kb->contractAccumulate(nullptr,
        "D(o,p,a,d,e,c)+=L(o,a,b)*R(p,c,b,d,e)",
        *tensor_grid[1][4].get(),
        *tensor_grid[1][5].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = kb->sync();
  unique_ptr<talsh::Tensor> kc(new talsh::Tensor({DIM, DIM, DIM, DIM,
                                    dims_grid[1][4][0],
                                    dims_grid[1][5][2],
                                    dims_grid[1][5][3],
                                    dims_grid[0][6][1]},
                                    talsh::TensorData<s_type>::kind,
                                    talsh_tens_no_init));
  errc = kc->contractAccumulate(nullptr,
        "D(q,r,o,p,a,b,c,e)+=L(o,p,a,b,c,d)*R(q,r,d,e)",
        *kb.get(),
        *ka.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = kc->sync();
  ka.reset();
  kb.reset();
  unique_ptr<talsh::Tensor> kd(new talsh::Tensor({DIM, DIM,
                                    dims_grid[1][6][0],
                                    dims_grid[1][6][1],
                                    dims_grid[1][6][2],
                                    dims_grid[1][7][1]},
                                    talsh::TensorData<s_type>::kind,
                                    talsh_tens_no_init));
  errc = kd->contractAccumulate(nullptr,
        "D(o,p,a,b,c,e)+=L(o,a,b,c,d)*R(p,d,e)",
        *tensor_grid[1][6].get(),
        *tensor_grid[1][7].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = kd->sync();
  unique_ptr<talsh::Tensor> K(new talsh::Tensor({DIM, DIM, DIM, DIM, DIM, DIM,
                                    dims_grid[1][4][0],
                                    dims_grid[1][5][2],
                                    dims_grid[1][6][2],
                                    dims_grid[1][7][1]},
                                    talsh::TensorData<s_type>::kind,
                                    talsh_tens_no_init));
  errc = K->contractAccumulate(nullptr,
        "D(o,p,q,r,s,t,a,b,e,f)+=L(o,p,q,r,a,b,c,d)*R(s,t,d,c,e,f)",
        *kc.get(),
        *kd.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = K->sync();
  kc.reset();
  kd.reset();
  /*
    unique_ptr<talsh::Tensor> C(new talsh::Tensor({dims_grid[1][4][0],
                                      dims_grid[1][5][2],
                                      dims_grid[1][6][2],
                                      dims_grid[1][7][1]},
                                      talsh::TensorData<s_type>::kind,
                                      talsh_tens_no_init));
    errc = C->contractAccumulate(nullptr,
          "D(a,b,e,f)+=L(a,b,c,d)*R(d,c,e,f)",
          *cc.get(),
          *cd.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = C->sync();
    cc.reset();
    cd.reset();
    // Contract onto scalar!!
    errc = S->contractAccumulate(nullptr,
          "D()+=L(a,b,c,d)*R(a,b,c,d)",
          *F.get(),
          *C.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = S->sync();
    C.reset();
    // Store result data_S
    s_type const * ptr_S;
    S->getDataAccessHostConst(&ptr_S);
    //amplitudes[c] += ptr_S[0];
    amplitudes[c] += ptr_S[0] / global_norm_factor;
    ptr_S = nullptr;
  }
  // Tensor C - UP TO HERE
  */







  // Tensor B
  unique_ptr<talsh::Tensor> ba(new talsh::Tensor({dims_grid[8][3][0],
                                                  dims_grid[8][4][0],
                                                  dims_grid[8][4][3]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
  errc = ba->contractAccumulate(nullptr,
        "D(a,c,d)+=L(a,b)*R(c,b,d)",
        *tensor_grid[8][3].get(),
        *corner_84.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = ba->sync();
  unique_ptr<talsh::Tensor> bb(new talsh::Tensor({dims_grid[8][3][0],
                                                  dims_grid[8][4][0],
                                                  dims_grid[8][5][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
  errc = bb->contractAccumulate(nullptr,
        "D(a,b,d)+=L(a,b,c)*R(d,c)",
        *ba.get(),
        *tensor_grid[8][5].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = bb->sync();
  ba.reset();
  // Group of 2 at the end of row 7
  unique_ptr<talsh::Tensor> bc(new talsh::Tensor({dims_grid[7][5][2],
                                                  dims_grid[7][5][1],
                                                  dims_grid[7][5][0],
                                                  dims_grid[7][6][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
  errc = bc->contractAccumulate(nullptr,
        "D(c,b,a,e)+=L(a,b,c,d)*R(e,d)",
        *tensor_grid[7][5].get(),
        *tensor_grid[7][6].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = bc->sync();

  // Up to here just use contractAccumulate

  unique_ptr<talsh::Tensor> bd(new talsh::Tensor({dims_grid[8][3][0],
                                                  dims_grid[8][4][0],
                                                  dims_grid[7][5][1],
                                                  dims_grid[7][5][0],
                                                  dims_grid[7][6][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
  errc = bd->contractAccumulateXL(nullptr,
        "D(a,b,d,e,f)+=L(a,b,c)*R(c,d,e,f)",
        *bb.get(),
        *bc.get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
  done = bd->sync();
  bb.reset();
  bc.reset();
  unique_ptr<talsh::Tensor> be(new talsh::Tensor({dims_grid[8][3][0],
                                                  dims_grid[7][4][1],
                                                  dims_grid[7][4][0],
                                                  dims_grid[7][5][0],
                                                  dims_grid[7][6][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
  // Tensor 7,4 (middle of row 7)
  errc = be->contractAccumulateXL(nullptr,
        "D(a,g,f,d,e)+=L(a,b,c,d,e)*R(f,g,b,c)",
        *bd.get(),
        *tensor_grid[7][4].get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
  done = be->sync();
  bd.reset();
  // Group of 2 at the beginning of row 7
  unique_ptr<talsh::Tensor> bf(new talsh::Tensor({dims_grid[7][2][0],
                                                  dims_grid[7][3][2],
                                                  dims_grid[7][3][3],
                                                  dims_grid[7][3][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
  errc = bf->contractAccumulate(nullptr,
        "D(a,d,e,c)+=L(a,b)*R(c,b,d,e)",
        *tensor_grid[7][2].get(),
        *tensor_grid[7][3].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = bf->sync();
  unique_ptr<talsh::Tensor> B(new talsh::Tensor({dims_grid[7][2][0],
                                                  dims_grid[7][3][0],
                                                  dims_grid[7][4][0],
                                                  dims_grid[7][5][0],
                                                  dims_grid[7][6][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
  errc = B->contractAccumulate(nullptr,
        "D(f,g,c,d,e)+=L(a,b,c,d,e)*R(f,a,b,g)",
        *be.get(),
        *bf.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = B->sync();
  be.reset();
  bf.reset();
  // bg is B

  // Tensor D comes here
  unique_ptr<talsh::Tensor> da(new talsh::Tensor({dims_grid[3][9][0],
                                                  dims_grid[4][9][1]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
  errc = da->contractAccumulate(nullptr,
        "D(a,c)+=L(a,b)*R(b,c)",
        *tensor_grid[3][9].get(),
        *tensor_grid[4][9].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = da->sync();
  unique_ptr<talsh::Tensor> db(new talsh::Tensor({dims_grid[3][8][0],
                                                  dims_grid[3][8][1],
                                                  dims_grid[3][8][2],
                                                  dims_grid[4][9][1]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
  errc = db->contractAccumulate(nullptr,
        "D(c,d,e,b)+=L(a,b)*R(c,d,e,a)",
        *da.get(),
        *tensor_grid[3][8].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = db->sync();
  da.reset();
  unique_ptr<talsh::Tensor> D(new talsh::Tensor({dims_grid[3][8][0],
                                                  dims_grid[3][8][1],
                                                  dims_grid[4][8][1],
                                                  dims_grid[4][8][2]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
  errc = D->contractAccumulate(nullptr,
        "D(a,b,e,f)+=L(a,b,c,d)*R(c,e,f,d)",
        *db.get(),
        *tensor_grid[4][8].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = D->sync();
  db.reset();
  // dc is D
  

  // Start first and second cut loops
  talsh::TensorTask task_hl;
  int first_cut_dim(dims_grid[3][3][2]);
  int second_cut_dim(dims_grid[4][3][3]);
  vector<int> first_cut_offsets({0});
  vector<int> second_cut_offsets({0});
  /*
  vector<int> first_cut_offsets;
  vector<int> second_cut_offsets;
  for (int p=0; p<first_cut_dim; ++p)
  {
    first_cut_offsets.push_back(p);
  }
  for (int p=0; p<second_cut_dim; ++p)
  {
    second_cut_offsets.push_back(p);
  }
  */
  int first_slice_size(1);
  int second_slice_size(32);
  for (auto i0 : first_cut_offsets) for (auto i1 : second_cut_offsets)
  {
    cout << "i0 = " << i0 << endl << flush;
    // First slice on [3][3] and [3][4].
    unique_ptr<talsh::Tensor> first_3_3(new talsh::Tensor({dims_grid[3][3][0],
                                                         dims_grid[3][3][1],
                                                         first_slice_size},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
    unique_ptr<talsh::Tensor> first_3_4(new talsh::Tensor({dims_grid[3][4][0],
                                                         first_slice_size,
                                                         dims_grid[3][4][2],
                                                         dims_grid[3][4][3]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
    errc = tensor_grid[3][3]->extractSlice(
                  &task_hl,*first_3_3,vector<int>{0,0,i0},DEV_HOST,0);
    done = tensor_grid[3][3]->sync();
    task_hl.clean();
    errc = tensor_grid[3][4]->extractSlice(
                  &task_hl,*first_3_4,vector<int>{0,i0,0,0},DEV_HOST,0);
    done = tensor_grid[3][4]->sync();
    task_hl.clean();
    // Second slice on [4][3] and [4][4].
    unique_ptr<talsh::Tensor> second_4_3(new talsh::Tensor({dims_grid[4][3][0],
                                                         dims_grid[4][3][1],
                                                         dims_grid[4][3][2],
                                                         second_slice_size},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
    unique_ptr<talsh::Tensor> second_4_4(new talsh::Tensor({dims_grid[4][4][0],
                                                         second_slice_size,
                                                         dims_grid[4][4][2],
                                                         dims_grid[4][4][3]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
    errc = tensor_grid[4][3]->extractSlice(
                  &task_hl,*second_4_3,vector<int>{0,0,0,i1},DEV_HOST,0);
    done = tensor_grid[4][3]->sync();
    task_hl.clean();
    errc = tensor_grid[4][4]->extractSlice(
                  &task_hl,*second_4_4,vector<int>{0,i1,0,0},DEV_HOST,0);
    done = tensor_grid[4][4]->sync();
    task_hl.clean();

    // Continue contraction
    // Finish A
    unique_ptr<talsh::Tensor> aa(new talsh::Tensor({dims_grid[3][2][0],
                                                    dims_grid[3][3][1],
                                                    first_slice_size},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
    errc = aa->contractAccumulate(nullptr,
          "D(a,c,d)+=L(a,b)*R(b,c,d)",
          *tensor_grid[3][2].get(),
          *first_3_3.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = aa->sync();
    unique_ptr<talsh::Tensor> ab(new talsh::Tensor({dims_grid[4][1][0],
                                                    dims_grid[4][2][2],
                                                    dims_grid[4][2][3],
                                                    dims_grid[3][3][1],
                                                    first_slice_size},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
    errc = ab->contractAccumulate(nullptr,
          "D(a,b,c,e,f)+=L(a,b,c,d)*R(d,e,f)",
          *pair_41_42.get(),
          *aa.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = ab->sync();
    aa.reset();
    unique_ptr<talsh::Tensor> ac(new talsh::Tensor({dims_grid[4][1][0],
                                                    dims_grid[4][2][2],
                                                    dims_grid[4][3][2],
                                                    second_slice_size,
                                                    first_slice_size},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
    errc = ac->contractAccumulate(nullptr,
          "D(a,b,f,g,e)+=L(a,b,c,d,e)*R(d,c,f,g)",
          *ab.get(),
          *second_4_3.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = ac->sync();
    ab.reset();
    unique_ptr<talsh::Tensor> ad(new talsh::Tensor({dims_grid[5][1][2],
                                                    dims_grid[5][1][3],
                                                    dims_grid[4][2][2],
                                                    dims_grid[4][3][2],
                                                    second_slice_size,
                                                    first_slice_size},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
    errc = ad->contractAccumulate(nullptr,
          "D(f,g,b,c,d,e)+=L(a,b,c,d,e)*R(a,f,g)",
          *ac.get(),
          *corner_51.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = ad->sync();
    ac.reset();
    unique_ptr<talsh::Tensor> ae(new talsh::Tensor({dims_grid[5][1][2],
                                                    dims_grid[5][2][2],
                                                    dims_grid[5][2][3],
                                                    dims_grid[4][3][2],
                                                    second_slice_size,
                                                    first_slice_size},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
    errc = ae->contractAccumulate(nullptr,
          "D(a,g,h,d,e,f)+=L(a,b,c,d,e,f)*R(c,b,g,h)",
          *ad.get(),
          *tensor_grid[5][2].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = ae->sync();
    ad.reset();
    unique_ptr<talsh::Tensor> af(new talsh::Tensor({dims_grid[5][1][2],
                                                    dims_grid[5][2][2],
                                                    dims_grid[5][3][2],
                                                    dims_grid[5][3][3],
                                                    second_slice_size,
                                                    first_slice_size},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
    errc = af->contractAccumulate(nullptr,
          "D(a,b,g,h,e,f)+=L(a,b,c,d,e,f)*R(d,c,g,h)",
          *ae.get(),
          *tensor_grid[5][3].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = af->sync();
    ae.reset();
    unique_ptr<talsh::Tensor> ag(new talsh::Tensor({dims_grid[6][2][2],
                                                    dims_grid[6][2][3],
                                                    dims_grid[5][3][2],
                                                    dims_grid[5][3][3],
                                                    second_slice_size,
                                                    first_slice_size},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
    errc = ag->contractAccumulate(nullptr,
          "D(g,h,c,d,e,f)+=L(a,b,c,d,e,f)*R(g,h,b,a)",
          *af.get(),
          *pair_61_62.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = ag->sync();
    af.reset();
    unique_ptr<talsh::Tensor> A(new talsh::Tensor({dims_grid[6][2][2],
                                                    dims_grid[6][3][2],
                                                    dims_grid[6][3][3],
                                                    dims_grid[5][3][3],
                                                    second_slice_size,
                                                    first_slice_size},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
    errc = A->contractAccumulate(nullptr,
          "D(a,g,h,d,e,f)+=L(a,b,c,d,e,f)*R(c,b,g,h)",
          *ag.get(),
          *tensor_grid[6][3].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = A->sync();
    ag.reset();
    // Done with A. ah is A.

    // Contract A with B into AB.
    unique_ptr<talsh::Tensor> AB(new talsh::Tensor({first_slice_size,
                                                    second_slice_size,
                                                    dims_grid[5][3][3],
                                                    dims_grid[6][3][3],
                                                    dims_grid[7][4][0],
                                                    dims_grid[7][5][0],
                                                    dims_grid[7][6][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
    errc = AB->contractAccumulateXL(nullptr,
          "D(f,e,d,c,g,h,i)+=L(a,b,c,d,e,f)*R(a,b,g,h,i)",
          *A.get(),
          *B.get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
    done = AB->sync();
    A.reset();

    // Contract AB with qubits in E, one by one.
    unique_ptr<talsh::Tensor> ea(new talsh::Tensor({first_slice_size,
                                                    second_slice_size,
                                                    dims_grid[5][3][3],
                                                    dims_grid[6][4][0],
                                                    dims_grid[6][4][3],
                                                    dims_grid[7][5][0],
                                                    dims_grid[7][6][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
    t0 = high_resolution_clock::now();
    errc = ea->contractAccumulateXL(nullptr,
          "D(a,b,c,h,i,f,g)+=L(a,b,c,d,e,f,g)*R(h,d,e,i)",
          *AB.get(),
          *tensor_grid[6][4].get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
    done = ea->sync();
    t1 = high_resolution_clock::now();
    span = duration_cast<duration<double>>(t1 - t0);
    //time_largest_contraction = double(span.count());
    cout << double(span.count()) << "s" << endl;
    AB.reset();
    unique_ptr<talsh::Tensor> eb(new talsh::Tensor({first_slice_size,
                                                    second_slice_size,
                                                    dims_grid[5][3][3],
                                                    dims_grid[6][4][0],
                                                    dims_grid[6][5][0],
                                                    dims_grid[6][5][3],
                                                    dims_grid[7][6][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));

    vector<int> dims_ea({first_slice_size*
                        second_slice_size*
                        dims_grid[5][3][3]*
                        dims_grid[6][4][0],
                        dims_grid[6][4][3],
                        dims_grid[7][5][0],
                        dims_grid[7][6][0]});
    ea->reshape(dims_ea);
    vector<int> dims_eb({first_slice_size*
                        second_slice_size*
                        dims_grid[5][3][3]*
                        dims_grid[6][4][0],
                        dims_grid[6][5][0],
                        dims_grid[6][5][3],
                        dims_grid[7][6][0]}); 
    eb->reshape(dims_eb);
    t0 = high_resolution_clock::now();
    errc = eb->contractAccumulateXL(nullptr,
          /*
          "D(a,b,c,d,h,i,g)+=L(a,b,c,d,e,f,g)*R(h,e,f,i)",
          */
          "D(a,h,i,g)+=L(a,e,f,g)*R(h,e,f,i)",
          *ea.get(),
          *tensor_grid[6][5].get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
    done = eb->sync();
    t1 = high_resolution_clock::now();
    dims_eb = vector<int>({first_slice_size*
                        second_slice_size*
                        dims_grid[5][3][3]*
                        dims_grid[6][4][0]*
                        dims_grid[6][5][0],
                        dims_grid[6][5][3],
                        dims_grid[7][6][0]}); 
    eb->reshape(dims_eb);

    span = duration_cast<duration<double>>(t1 - t0);
    time_largest_contraction = double(span.count());
    //cout << double(span.count()) << "s" << endl;
    ea.reset();
    unique_ptr<talsh::Tensor> ec(new talsh::Tensor({first_slice_size,
                                                    second_slice_size,
                                                    dims_grid[5][3][3],
                                                    dims_grid[6][4][0],
                                                    dims_grid[6][5][0],
                                                    dims_grid[6][6][0],
                                                    dims_grid[6][7][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));

    vector<int> dims_ec({first_slice_size*
                        second_slice_size*
                        dims_grid[5][3][3]*
                        dims_grid[6][4][0]*
                        dims_grid[6][5][0],
                        dims_grid[6][6][0],
                        dims_grid[6][7][0]});
    ec->reshape(dims_ec);
    errc = ec->contractAccumulateXL(nullptr,
          /*
          "D(a,b,c,d,e,i,h)+=L(a,b,c,d,e,f,g)*R(h,i,f,g)",
          */
          "D(a,i,h)+=L(a,f,g)*R(h,i,f,g)",
          *eb.get(),
          *pair_66_67.get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
    done = ec->sync();
    dims_ec = vector<int>({first_slice_size,
                        second_slice_size,
                        dims_grid[5][3][3],
                        dims_grid[6][4][0],
                        dims_grid[6][5][0],
                        dims_grid[6][6][0],
                        dims_grid[6][7][0]});
    ec->reshape(dims_ec);
    eb.reset();

    // Row 5
    unique_ptr<talsh::Tensor> ed(new talsh::Tensor({first_slice_size,
                                                    second_slice_size,
                                                    dims_grid[5][4][0],
                                                    dims_grid[5][4][3],
                                                    dims_grid[6][5][0],
                                                    dims_grid[6][6][0],
                                                    dims_grid[6][7][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
    errc = ed->contractAccumulateXL(nullptr,
          "D(a,b,h,i,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,c,d,i)",
          *ec.get(),
          *tensor_grid[5][4].get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
    done = ed->sync();
    ec.reset();
    unique_ptr<talsh::Tensor> ee(new talsh::Tensor({first_slice_size,
                                                    second_slice_size,
                                                    dims_grid[5][4][0],
                                                    dims_grid[5][5][0],
                                                    dims_grid[5][5][3],
                                                    dims_grid[6][6][0],
                                                    dims_grid[6][7][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
    errc = ee->contractAccumulateXL(nullptr,
          "D(a,b,c,h,i,f,g)+=L(a,b,c,d,e,f,g)*R(h,d,e,i)",
          *ed.get(),
          *tensor_grid[5][5].get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
    done = ee->sync();
    ed.reset();
    unique_ptr<talsh::Tensor> ef(new talsh::Tensor({first_slice_size,
                                                    second_slice_size,
                                                    dims_grid[5][4][0],
                                                    dims_grid[5][5][0],
                                                    dims_grid[5][6][0],
                                                    dims_grid[5][6][3],
                                                    dims_grid[6][7][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
    vector<int> dims_ee({first_slice_size*
                    second_slice_size*
                    dims_grid[5][4][0]*
                    dims_grid[5][5][0],
                    dims_grid[5][5][3],
                    dims_grid[6][6][0],
                    dims_grid[6][7][0]});
    ee->reshape(dims_ee);
    vector<int> dims_ef({first_slice_size*
                        second_slice_size*
                        dims_grid[5][4][0]*
                        dims_grid[5][5][0],
                        dims_grid[5][6][0],
                        dims_grid[5][6][3],
                        dims_grid[6][7][0]});
    ef->reshape(dims_ef);
    errc = ef->contractAccumulateXL(nullptr,
          /*
          "D(a,b,c,d,h,i,g)+=L(a,b,c,d,e,f,g)*R(h,e,f,i)",
          */
          "D(a,h,i,g)+=L(a,e,f,g)*R(h,e,f,i)",
          *ee.get(),
          *tensor_grid[5][6].get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
    done = ef->sync();
    dims_ef = vector<int>({first_slice_size*
                        second_slice_size*
                        dims_grid[5][4][0]*
                        dims_grid[5][5][0]*
                        dims_grid[5][6][0],
                        dims_grid[5][6][3],
                        dims_grid[6][7][0]});
    ef->reshape(dims_ef);
    ee.reset();
    unique_ptr<talsh::Tensor> E(new talsh::Tensor({first_slice_size,
                                                    second_slice_size,
                                                    dims_grid[5][4][0],
                                                    dims_grid[5][5][0],
                                                    dims_grid[5][6][0],
                                                    dims_grid[5][7][0],
                                                    dims_grid[5][8][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
    vector<int> dims_E({first_slice_size*
                        second_slice_size*
                        dims_grid[5][4][0]*
                        dims_grid[5][5][0]*
                        dims_grid[5][6][0],
                        dims_grid[5][7][0],
                        dims_grid[5][8][0]});
    E->reshape(dims_E);
    errc = E->contractAccumulateXL(nullptr,
          /*
          "D(a,b,c,d,e,i,h)+=L(a,b,c,d,e,f,g)*R(h,i,f,g)",
          */
          "D(a,i,h)+=L(a,f,g)*R(h,i,f,g)",
          *ef.get(),
          *pair_57_58.get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
    done = E->sync();
    dims_E = vector<int>({first_slice_size,
                          second_slice_size,
                          dims_grid[5][4][0],
                          dims_grid[5][5][0],
                          dims_grid[5][6][0],
                          dims_grid[5][7][0],
                          dims_grid[5][8][0]});
    E->reshape(dims_E);
    ef.reset();

    // Start third cut loop
    int third_cut_dim(dims_grid[3][4][2]);
    //vector<int> third_cut_offsets({0});
    vector<int> third_cut_offsets({0});
    /*
    vector<int> third_cut_offsets;
    for (int p=0; p<third_cut_dim; ++p)
    {
      third_cut_offsets.push_back(p);
    }
    */
    int third_slice_size(32);
    for (auto i2 : third_cut_offsets)
    {
      cout << "i2 = " << i2 << endl << flush;
      // Time i2 loop.
      t0 = high_resolution_clock::now();
      unique_ptr<talsh::Tensor> third_3_4(new talsh::Tensor({dims_grid[3][4][0],
                                                           first_slice_size,
                                                           third_slice_size,
                                                           dims_grid[3][4][3]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
      unique_ptr<talsh::Tensor> third_4_4(new talsh::Tensor({third_slice_size,
                                                           second_slice_size,
                                                           dims_grid[4][4][2],
                                                           dims_grid[4][4][3]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
      errc = first_3_4->extractSlice(
                    &task_hl,*third_3_4,vector<int>{0,0,i2,0},DEV_HOST,0);
      done = first_3_4->sync();
      task_hl.clean();
      errc = second_4_4->extractSlice(
                    &task_hl,*third_4_4,vector<int>{i2,0,0,0},DEV_HOST,0);
      done = second_4_4->sync();
      task_hl.clean();

      // Continue contraction of E into F, after third cut. E is stored.
      unique_ptr<talsh::Tensor> fa(new talsh::Tensor({first_slice_size,
                                                      third_slice_size,
                                                      dims_grid[4][4][3],
                                                      dims_grid[5][5][0],
                                                      dims_grid[5][6][0],
                                                      dims_grid[5][7][0],
                                                      dims_grid[5][8][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
      errc = fa->contractAccumulateXL(nullptr,
            "D(a,h,i,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,b,c,i)",
            *E.get(),
            *third_4_4.get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
      done = fa->sync();
      unique_ptr<talsh::Tensor> fb(new talsh::Tensor({first_slice_size,
                                                      third_slice_size,
                                                      dims_grid[4][5][0],
                                                      dims_grid[4][5][3],
                                                      dims_grid[5][6][0],
                                                      dims_grid[5][7][0],
                                                      dims_grid[5][8][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
      errc = fb->contractAccumulateXL(nullptr,
            "D(a,b,h,i,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,c,d,i)",
            *fa.get(),
            *tensor_grid[4][5].get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
      done = fb->sync();
      fa.reset();
      unique_ptr<talsh::Tensor> fc(new talsh::Tensor({first_slice_size,
                                                      third_slice_size,
                                                      dims_grid[4][5][0],
                                                      dims_grid[4][6][0],
                                                      dims_grid[4][6][3],
                                                      dims_grid[5][7][0],
                                                      dims_grid[5][8][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
      errc = fc->contractAccumulateXL(nullptr,
            "D(a,b,c,h,i,f,g)+=L(a,b,c,d,e,f,g)*R(h,d,e,i)",
            *fb.get(),
            *tensor_grid[4][6].get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
      done = fc->sync();
      fb.reset();
      unique_ptr<talsh::Tensor> fd(new talsh::Tensor({first_slice_size,
                                                      third_slice_size,
                                                      dims_grid[4][5][0],
                                                      dims_grid[4][6][0],
                                                      dims_grid[4][7][0],
                                                      dims_grid[4][7][3],
                                                      dims_grid[5][8][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
      vector<int> dims_fc({first_slice_size*
                          third_slice_size*
                          dims_grid[4][5][0]*
                          dims_grid[4][6][0],
                          dims_grid[4][6][3],
                          dims_grid[5][7][0],
                          dims_grid[5][8][0]});
      fc->reshape(dims_fc);
      vector<int> dims_fd({first_slice_size*
                          third_slice_size*
                          dims_grid[4][5][0]*
                          dims_grid[4][6][0],
                          dims_grid[4][7][0],
                          dims_grid[4][7][3],
                          dims_grid[5][8][0]});
      fd->reshape(dims_fd);
      errc = fd->contractAccumulateXL(nullptr,
            /*
            "D(a,b,c,d,h,i,g)+=L(a,b,c,d,e,f,g)*R(h,e,f,i)",
            */
            "D(a,h,i,g)+=L(a,e,f,g)*R(h,e,f,i)",
            *fc.get(),
            *tensor_grid[4][7].get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
      done = fd->sync();
      dims_fd = vector<int>({first_slice_size*
                          third_slice_size*
                          dims_grid[4][5][0]*
                          dims_grid[4][6][0]*
                          dims_grid[4][7][0],
                          dims_grid[4][7][3],
                          dims_grid[5][8][0]});
      fd->reshape(dims_fd);
      fc.reset();
      // Now contract fd with D!
      cout << "Before D" << endl << flush;
      unique_ptr<talsh::Tensor> fe(new talsh::Tensor({first_slice_size,
                                                      third_slice_size,
                                                      dims_grid[4][5][0],
                                                      dims_grid[4][6][0],
                                                      dims_grid[4][7][0],
                                                      dims_grid[3][8][1],
                                                      dims_grid[3][8][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
      vector<int> dims_fe({first_slice_size*
                          third_slice_size*
                          dims_grid[4][5][0]*
                          dims_grid[4][6][0]*
                          dims_grid[4][7][0],
                          dims_grid[3][8][1],
                          dims_grid[3][8][0]});
      fe->reshape(dims_fe);
      errc = fe->contractAccumulateXL(nullptr,
            /*
            "D(a,b,c,d,e,i,h)+=L(a,b,c,d,e,f,g)*R(h,i,f,g)",
            */
            "D(a,i,h)+=L(a,f,g)*R(h,i,f,g)",
            *fd.get(),
            *D.get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
      done = fe->sync();
      dims_fe = vector<int>({first_slice_size*
                          third_slice_size*
                          dims_grid[4][5][0]*
                          dims_grid[4][6][0],
                          dims_grid[4][7][0],
                          dims_grid[3][8][1],
                          dims_grid[3][8][0]});
      fe->reshape(dims_fe);
      fd.reset();
      // Continue with row 3.
      cout << "After D" << endl << flush;
      unique_ptr<talsh::Tensor> ff(new talsh::Tensor({first_slice_size,
                                                      third_slice_size,
                                                      dims_grid[4][5][0],
                                                      dims_grid[4][6][0],
                                                      dims_grid[3][7][1],
                                                      dims_grid[3][7][0],
                                                      dims_grid[3][8][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
      vector<int> dims_ff({first_slice_size*
                      third_slice_size*
                      dims_grid[4][5][0]*
                      dims_grid[4][6][0],
                      dims_grid[3][7][1],
                      dims_grid[3][7][0],
                      dims_grid[3][8][0]});
      ff->reshape(dims_ff);
      errc = ff->contractAccumulateXL(nullptr,
            /*
            "D(a,b,c,d,i,h,g)+=L(a,b,c,d,e,f,g)*R(h,i,e,f)",
            */
            "D(a,i,h,g)+=L(a,e,f,g)*R(h,i,e,f)",
            *fe.get(),
            *tensor_grid[3][7].get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
      done = ff->sync();
      dims_ff = vector<int>({first_slice_size,
                            third_slice_size,
                            dims_grid[4][5][0],
                            dims_grid[4][6][0],
                            dims_grid[3][7][1],
                            dims_grid[3][7][0],
                            dims_grid[3][8][0]});
      ff->reshape(dims_ff);
      fe.reset();
      unique_ptr<talsh::Tensor> fg(new talsh::Tensor({first_slice_size,
                                                      third_slice_size,
                                                      dims_grid[4][5][0],
                                                      dims_grid[3][6][1],
                                                      dims_grid[3][6][0],
                                                      dims_grid[3][7][0],
                                                      dims_grid[3][8][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
      errc = fg->contractAccumulateXL(nullptr,
            "D(a,b,c,i,h,f,g)+=L(a,b,c,d,e,f,g)*R(h,i,d,e)",
            *ff.get(),
            *tensor_grid[3][6].get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
      done = fg->sync();
      ff.reset();
      unique_ptr<talsh::Tensor> fh(new talsh::Tensor({first_slice_size,
                                                      third_slice_size,
                                                      dims_grid[3][5][1],
                                                      dims_grid[3][5][0],
                                                      dims_grid[3][6][0],
                                                      dims_grid[3][7][0],
                                                      dims_grid[3][8][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
      errc = fh->contractAccumulateXL(nullptr,
            "D(a,b,i,h,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,i,c,d)",
            *fg.get(),
            *tensor_grid[3][5].get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
      done = fh->sync();
      fg.reset();
      unique_ptr<talsh::Tensor> fi(new talsh::Tensor({dims_grid[3][4][0],
                                                      dims_grid[3][5][0],
                                                      dims_grid[3][6][0],
                                                      dims_grid[3][7][0],
                                                      dims_grid[3][8][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
      errc = fi->contractAccumulateXL(nullptr,
            "D(h,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,a,b,c)",
            *fh.get(),
            *third_3_4.get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
      done = fi->sync();
      fh.reset();
      // Row 2 with F sequence
      unique_ptr<talsh::Tensor> fj(new talsh::Tensor({dims_grid[3][4][0],
                                                      dims_grid[3][5][0],
                                                      dims_grid[3][6][0],
                                                      dims_grid[2][7][1],
                                                      dims_grid[2][7][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
      errc = fj->contractAccumulate(nullptr,
            "D(a,b,c,g,f)+=L(a,b,c,d,e)*R(f,g,d,e)",
            *fi.get(),
            *pair_27_28.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
      done = fj->sync();
      fi.reset();
      unique_ptr<talsh::Tensor> fk(new talsh::Tensor({dims_grid[3][4][0],
                                                      dims_grid[3][5][0],
                                                      dims_grid[2][6][1],
                                                      dims_grid[2][6][0],
                                                      dims_grid[2][7][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
      errc = fk->contractAccumulateXL(nullptr,
            "D(a,b,g,f,e)+=L(a,b,c,d,e)*R(f,g,c,d)",
            *fj.get(),
            *tensor_grid[2][6].get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
      done = fk->sync();
      fj.reset();
      unique_ptr<talsh::Tensor> fl(new talsh::Tensor({dims_grid[3][4][0],
                                                      dims_grid[2][5][1],
                                                      dims_grid[2][5][0],
                                                      dims_grid[2][6][0],
                                                      dims_grid[2][7][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
      errc = fl->contractAccumulateXL(nullptr,
            "D(a,g,f,d,e)+=L(a,b,c,d,e)*R(f,g,b,c)",
            *fk.get(),
            *tensor_grid[2][5].get(), DEV_NVIDIA_GPU, DEV_DEFAULT, s_type(1.0), false);
      done = fl->sync();
      fk.reset();
      unique_ptr<talsh::Tensor> F(new talsh::Tensor({dims_grid[2][4][0],
                                                      dims_grid[2][5][0],
                                                      dims_grid[2][6][0],
                                                      dims_grid[2][7][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
      errc = F->contractAccumulate(nullptr,
            "D(f,c,d,e)+=L(a,b,c,d,e)*R(f,a,b)",
            *fl.get(),
            *tensor_grid[2][4].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
      done = F->sync();
      fl.reset();
      // fm is F. Done with F. Start fast sampling on corner.
   

      /*
      // C loop!
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
          string contraction_string;
          if (t==0 || t==1 || t==2 || t==5)
          {
            contraction_string = "D(b,c)+=L(a)*R(a,b,c)";
          } else {
            contraction_string = "D(b,c,d,e)+=L(a)*R(a,b,c,d,e)";
          }
          errc = Cs[t]->contractAccumulate(nullptr, contraction_string,
                delta,
                *tensor_grid[i][j].get(),
                DEV_NVIDIA_GPU, 0, s_type(1.0), false);
          done = Cs[t]->sync();
        }

        // Contract C
        unique_ptr<talsh::Tensor> ca(new talsh::Tensor({dims_grid[0][5][0],
                                          dims_grid[0][6][1]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
        errc = ca->contractAccumulate(nullptr,
              "D(a,c)+=L(a,b)*R(b,c)",
              *Cs[0].get(),
              *Cs[1].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
        done = ca->sync();
        unique_ptr<talsh::Tensor> cb(new talsh::Tensor({dims_grid[1][4][0],
                                          dims_grid[1][5][2],
                                          dims_grid[1][5][3],
                                          dims_grid[1][5][0]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
        errc = cb->contractAccumulate(nullptr,
              "D(a,d,e,c)+=L(a,b)*R(c,b,d,e)",
              *Cs[2].get(),
              *Cs[3].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
        done = cb->sync();
        unique_ptr<talsh::Tensor> cc(new talsh::Tensor({dims_grid[1][4][0],
                                          dims_grid[1][5][2],
                                          dims_grid[1][5][3],
                                          dims_grid[0][6][1]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
        errc = cc->contractAccumulate(nullptr,
              "D(a,b,c,e)+=L(a,b,c,d)*R(d,e)",
              *cb.get(),
              *ca.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
        done = cc->sync();
        ca.reset();
        cb.reset();
        unique_ptr<talsh::Tensor> cd(new talsh::Tensor({dims_grid[1][6][0],
                                          dims_grid[1][6][1],
                                          dims_grid[1][6][2],
                                          dims_grid[1][7][1]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
        errc = cd->contractAccumulate(nullptr,
              "D(a,b,c,e)+=L(a,b,c,d)*R(d,e)",
              *Cs[4].get(),
              *Cs[5].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
        done = cd->sync();
        unique_ptr<talsh::Tensor> C(new talsh::Tensor({dims_grid[1][4][0],
                                          dims_grid[1][5][2],
                                          dims_grid[1][6][2],
                                          dims_grid[1][7][1]},
                                          talsh::TensorData<s_type>::kind,
                                          talsh_tens_no_init));
        errc = C->contractAccumulate(nullptr,
              "D(a,b,e,f)+=L(a,b,c,d)*R(d,c,e,f)",
              *cc.get(),
              *cd.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
        done = C->sync();
        cc.reset();
        cd.reset();
        // Contract onto scalar!!
        errc = S->contractAccumulate(nullptr,
              "D()+=L(a,b,c,d)*R(a,b,c,d)",
              *F.get(),
              *C.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
        done = S->sync();
        C.reset();
        // Store result data_S
        s_type const * ptr_S;
        S->getDataAccessHostConst(&ptr_S);
        amplitudes[c] += ptr_S[0] / global_norm_factor;
        ptr_S = nullptr;
      }
      */



      unique_ptr<talsh::Tensor> S_vector(new talsh::Tensor({
                                        DIM, DIM, DIM, DIM, DIM, DIM},
                                        talsh::TensorData<s_type>::kind,
                                        talsh_tens_no_init));
      errc = S_vector->contractAccumulate(nullptr,
            "D(o,p,q,r,s,t)+=L(a,b,c,d)*R(o,p,q,r,s,t,a,b,c,d)",
            *F.get(),
            *K.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
      done = S_vector->sync();
      s_type const * ptr_S_vector;
      S_vector->getDataAccessHostConst(&ptr_S_vector);
      for (int c=0; c<num_Cs; ++c) 
      {
        amplitudes[c] += ptr_S_vector[c] / global_norm_factor;
      }



      
      t1 = high_resolution_clock::now();
      span = duration_cast<duration<double>>(t1 - t0);
      cout << "Time in i2 loop: " << span.count() << endl;
      F.reset();
    }
    
  }

}


vector<s_type> const & Contraction::get_amplitudes() const
{
  return amplitudes;
}

double Contraction::get_time_largest_contraction() const
{
  return time_largest_contraction;
}

/////////////////////////// EXTERNAL FUNCTIONS ////////////////////////////////
// NONE
