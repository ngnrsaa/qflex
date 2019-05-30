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

#include "contraction_sycamore_51x10.h"
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
       {5,0},                                                {5,9},
       {6,0},                                          {6,8},{6,9},
       {7,0},{7,1},                              {7,7},{7,8},{7,9},
       {8,0},{8,1},{8,2},                  {8,6},{8,7},{8,8},{8,9},
       {9,0},{9,1},{9,2},{9,3},{9,4},{9,5},{9,6},{9,7},{9,8},{9,9}
                                    }); 
  qubits_A = vector<vector<int>>({       {0,5},{0,6},
                                   {1,4},{1,5},{1,6},{1,7} });

   
  int errc;


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
  // Column 1
  unique_ptr<talsh::Tensor> aa(new talsh::Tensor({dims_grid[5][1][1],
                                                  dims_grid[5][1][2],
                                                  dims_grid[4][1][1]},
                                                  s_type(0.0)));
  errc = aa->contractAccumulate(nullptr,
        "D(c,d,b)+=L(a,b)*R(a,c,d)",
        *tensor_grid[4][1].get(),
        *tensor_grid[5][1].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = aa->sync();
  unique_ptr<talsh::Tensor> ab(new talsh::Tensor({dims_grid[6][1][1],
                                                  dims_grid[5][1][2],
                                                  dims_grid[4][1][1]},
                                                  s_type(0.0)));
  errc = ab->contractAccumulate(nullptr,
        "D(d,b,c)+=L(a,b,c)*R(a,d)",
        *aa.get(),
        *tensor_grid[6][1].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = ab->sync();
  aa.reset();
  // Column 2
  // First group of 2
  unique_ptr<talsh::Tensor> ac(new talsh::Tensor({dims_grid[4][2][1],
                                                  dims_grid[4][2][2],
                                                  dims_grid[4][2][3],
                                                  dims_grid[3][2][1]},
                                                  s_type(0.0)));
  errc = ac->contractAccumulate(nullptr,
        "D(b,c,d,e)+=L(a,b,c,d)*R(a,e)",
        *tensor_grid[4][2].get(),
        *tensor_grid[3][2].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = ac->sync();
  unique_ptr<talsh::Tensor> ad(new talsh::Tensor({dims_grid[6][1][1],
                                                  dims_grid[5][1][2],
                                                  dims_grid[4][2][2],
                                                  dims_grid[4][2][3],
                                                  dims_grid[3][2][1]},
                                                  s_type(0.0)));
  errc = ad->contractAccumulate(nullptr,
        "D(a,b,d,e,f)+=L(a,b,c)*R(c,d,e,f)",
        *ab.get(),
        *ac.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = ad->sync();
  ab.reset();
  ac.reset();
  unique_ptr<talsh::Tensor> ae(new talsh::Tensor({dims_grid[6][1][1],
                                                  dims_grid[5][2][2],
                                                  dims_grid[5][2][3],
                                                  dims_grid[4][2][3],
                                                  dims_grid[3][2][1]},
                                                  s_type(0.0)));
  errc = ae->contractAccumulate(nullptr,
        "D(a,f,g,d,e)+=L(a,b,c,d,e)*R(c,b,f,g)",
        *ad.get(),
        *tensor_grid[5][2].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = ae->sync();
  ad.reset();
  unique_ptr<talsh::Tensor> af(new talsh::Tensor({dims_grid[6][2][0],
                                                  dims_grid[6][2][1],
                                                  dims_grid[7][2][1],
                                                  dims_grid[6][2][3]},
                                                  s_type(0.0)));
  errc = af->contractAccumulate(nullptr,
        "D(c,d,b,e)+=L(a,b)*R(c,d,a,e)",
        *tensor_grid[7][2].get(),
        *tensor_grid[6][2].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = af->sync();
  unique_ptr<talsh::Tensor> ag(new talsh::Tensor({dims_grid[7][2][1],
                                                  dims_grid[6][2][3],
                                                  dims_grid[5][2][3],
                                                  dims_grid[4][2][3],
                                                  dims_grid[3][2][1]},
                                                  s_type(0.0)));
  errc = ag->contractAccumulate(nullptr,
        "D(f,g,c,d,e)+=L(a,b,c,d,e)*R(b,a,f,g)",
        *ae.get(),
        *af.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = ae->sync();
  ae.reset();
  af.reset();
  // ag is the pre-first-cut tensor.

  // Tensor B
  unique_ptr<talsh::Tensor> ba(new talsh::Tensor({dims_grid[8][3][0],
                                                  dims_grid[8][4][0],
                                                  dims_grid[8][4][2]},
                                                  s_type(0.0)));
  errc = ba->contractAccumulate(nullptr,
        "D(a,c,d)+=L(a,b)*R(c,b,d)",
        *tensor_grid[8][3].get(),
        *tensor_grid[8][4].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = ba->sync();
  unique_ptr<talsh::Tensor> bb(new talsh::Tensor({dims_grid[8][3][0],
                                                  dims_grid[8][4][0],
                                                  dims_grid[8][5][0]},
                                                  s_type(0.0)));
  errc = bb->contractAccumulate(nullptr,
        "D(a,b,d)+=L(a,b,c)*R(d,c)",
        *ba.get(),
        *tensor_grid[8][5].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = bb->sync();
  ba.reset();
  unique_ptr<talsh::Tensor> bc(new talsh::Tensor({dims_grid[7][5][2],
                                                  dims_grid[7][5][1],
                                                  dims_grid[7][5][0],
                                                  dims_grid[7][6][0]},
                                                  s_type(0.0)));
  errc = bc->contractAccumulate(nullptr,
        "D(c,b,a,e)+=L(a,b,c,d)*R(e,d)",
        *tensor_grid[7][5].get(),
        *tensor_grid[7][6].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = bc->sync();
  unique_ptr<talsh::Tensor> bd(new talsh::Tensor({dims_grid[8][3][0],
                                                  dims_grid[8][4][0],
                                                  dims_grid[7][5][1],
                                                  dims_grid[7][5][0],
                                                  dims_grid[7][6][0]},
                                                  s_type(0.0)));
  errc = bd->contractAccumulate(nullptr,
        "D(a,b,d,e,f)+=L(a,b,c)*R(c,d,e,f)",
        *bb.get(),
        *bc.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = bd->sync();
  bb.reset();
  bc.reset();
  unique_ptr<talsh::Tensor> be(new talsh::Tensor({dims_grid[8][3][0],
                                                  dims_grid[7][4][1],
                                                  dims_grid[7][4][0],
                                                  dims_grid[7][5][0],
                                                  dims_grid[7][6][0]},
                                                  s_type(0.0)));
  errc = be->contractAccumulate(nullptr,
        "D(a,g,f,d,e)+=L(a,b,c,d,e)*R(f,g,b,c)",
        *bd.get(),
        *tensor_grid[7][4].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = be->sync();
  bd.reset();
  unique_ptr<talsh::Tensor> bf(new talsh::Tensor({dims_grid[6][6][2],
                                                  dims_grid[6][6][1],
                                                  dims_grid[6][6][0],
                                                  dims_grid[6][7][0]},
                                                  s_type(0.0)));
  errc = bf->contractAccumulate(nullptr,
        "D(c,b,a,e)+=L(a,b,c,d)*R(e,d)",
        *tensor_grid[6][6].get(),
        *tensor_grid[6][7].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = bf->sync();
  unique_ptr<talsh::Tensor> bg(new talsh::Tensor({dims_grid[8][3][0],
                                                  dims_grid[7][4][1],
                                                  dims_grid[7][4][0],
                                                  dims_grid[7][5][0],
                                                  dims_grid[6][6][1],
                                                  dims_grid[6][6][0],
                                                  dims_grid[6][7][0]},
                                                  s_type(0.0)));
  errc = bg->contractAccumulate(nullptr,
        "D(a,b,c,d,f,g,h)+=L(a,b,c,d,e)*R(e,f,g,h)",
        *be.get(),
        *bf.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = bg->sync();
  be.reset();
  bf.reset();
  unique_ptr<talsh::Tensor> bh(new talsh::Tensor({dims_grid[8][3][0],
                                                  dims_grid[7][4][1],
                                                  dims_grid[7][4][0],
                                                  dims_grid[6][5][1],
                                                  dims_grid[6][5][0],
                                                  dims_grid[6][6][0],
                                                  dims_grid[6][7][0]},
                                                  s_type(0.0)));
  errc = bh->contractAccumulate(nullptr,
        "D(a,b,c,i,h,f,g)+=L(a,b,c,d,e,f,g)*R(h,i,d,e)",
        *bg.get(),
        *tensor_grid[6][5].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = bh->sync();
  /*
  errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
  assert(tc.sync(DEV_HOST,0));
  */
  bg.reset();
  unique_ptr<talsh::Tensor> B(new talsh::Tensor({dims_grid[8][3][0],
                                                  dims_grid[7][4][1],
                                                  dims_grid[6][4][1],
                                                  dims_grid[6][4][0],
                                                  dims_grid[6][5][0],
                                                  dims_grid[6][6][0],
                                                  dims_grid[6][7][0]},
                                                  s_type(0.0)));
  errc = B->contractAccumulate(nullptr,
        "D(a,b,i,h,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,i,c,d)",
        *bh.get(),
        *tensor_grid[6][4].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = B->sync();
  bh.reset();
  // bi is B

  // Tensor D comes here
  unique_ptr<talsh::Tensor> da(new talsh::Tensor({dims_grid[4][9][1],
                                                  dims_grid[3][9][0]},
                                                  s_type(0.0)));
  errc = da->contractAccumulate(nullptr,
        "D(b,c)+=L(a,b)*R(c,a)",
        *tensor_grid[4][9].get(),
        *tensor_grid[3][9].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = da->sync();
  unique_ptr<talsh::Tensor> db(new talsh::Tensor({dims_grid[5][8][1],
                                                  dims_grid[4][8][1],
                                                  dims_grid[4][8][0],
                                                  dims_grid[4][8][3]},
                                                  s_type(0.0)));
  errc = db->contractAccumulate(nullptr,
        "D(b,d,c,e)+=L(a,b)*R(c,d,a,e)",
        *tensor_grid[5][8].get(),
        *tensor_grid[4][8].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = db->sync();
  unique_ptr<talsh::Tensor> dc(new talsh::Tensor({dims_grid[5][8][1],
                                                  dims_grid[4][8][1],
                                                  dims_grid[4][8][0],
                                                  dims_grid[3][9][0]},
                                                  s_type(0.0)));
  errc = dc->contractAccumulate(nullptr,
        "D(a,b,c,e)+=L(a,b,c,d)*R(d,e)",
        *db.get(),
        *da.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = dc->sync();
  db.reset();
  da.reset();
  unique_ptr<talsh::Tensor> dd(new talsh::Tensor({dims_grid[2][8][0],
                                                  dims_grid[3][8][1],
                                                  dims_grid[3][8][2],
                                                  dims_grid[3][8][3]},
                                                  s_type(0.0)));
  errc = dd->contractAccumulate(nullptr,
        "D(e,b,c,d)+=L(a,b,c,d)*R(e,a)",
        *tensor_grid[3][8].get(),
        *tensor_grid[2][8].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = dd->sync();
  unique_ptr<talsh::Tensor> D(new talsh::Tensor({dims_grid[5][8][1],
                                                  dims_grid[4][8][1],
                                                  dims_grid[3][8][1],
                                                  dims_grid[2][8][0]},
                                                  s_type(0.0)));
  errc = D->contractAccumulate(nullptr,
        "D(a,b,f,e)+=L(a,b,c,d)*R(e,f,c,d)",
        *dc.get(),
        *dd.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = D->sync();
  dc.reset();
  dd.reset();
  // de is D
  

  // Start first cut loop
  talsh::TensorTask task_hl;
  vector<int> first_cut_offsets({0});
  int first_cut_dim(dims_grid[3][3][2]);
  int first_slice_size(1);
  for (auto i0 : first_cut_offsets)
  {
    // First slice on [3][3] and [4][3].
    unique_ptr<talsh::Tensor> first_3_3(new talsh::Tensor({dims_grid[3][3][0],
                                                         dims_grid[3][3][1],
                                                         first_slice_size},
                                                         s_type(0.0)));
    unique_ptr<talsh::Tensor> first_4_3(new talsh::Tensor({dims_grid[4][3][0],
                                                         first_slice_size,
                                                         dims_grid[4][3][2],
                                                         dims_grid[4][3][3]},
                                                         s_type(0.0)));
    errc = tensor_grid[3][3]->extractSlice(
                  &task_hl,*first_3_3,vector<int>{0,0,i0},DEV_HOST,0);
    done = tensor_grid[3][3]->sync();
    task_hl.clean();
    errc = tensor_grid[4][3]->extractSlice(
                  &task_hl,*first_4_3,vector<int>{0,i0,0,0},DEV_HOST,0);
    done = tensor_grid[4][3]->sync();
    task_hl.clean();

    // Continue contraction
    // Column 3
    unique_ptr<talsh::Tensor> ah(new talsh::Tensor({dims_grid[7][2][1],
                                                    dims_grid[6][2][3],
                                                    dims_grid[5][2][3],
                                                    dims_grid[4][2][3],
                                                    dims_grid[3][3][1],
                                                    first_slice_size},
                                                    s_type(0.0)));
    errc = ah->contractAccumulate(nullptr,
          "D(a,b,c,d,f,g)+=L(a,b,c,d,e)*R(e,f,g)",
          *ag.get(),
          *first_3_3.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = ah->sync();
    ag.reset();
    unique_ptr<talsh::Tensor> ai(new talsh::Tensor({dims_grid[7][2][1],
                                                    dims_grid[6][2][3],
                                                    dims_grid[5][2][3],
                                                    dims_grid[4][3][2],
                                                    dims_grid[4][3][3],
                                                    first_slice_size},
                                                    s_type(0.0)));
    errc = ai->contractAccumulate(nullptr,
          "D(a,b,c,g,h,f)+=L(a,b,c,d,e,f)*R(e,d,g,h)",
          *ah.get(),
          *tensor_grid[4][3].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = ai->sync();
    ah.reset();
    
  }
  /*
  // Start outer loop (in practice it will take only one value
  vector<int> outer_values({3});
  //vector<int> outer_values({0,1,2,3});
  for (auto i0 : outer_values)
  //for (int i0=0; i0<1; ++i0)
  {

    // A
    // Group of 3
    TensContraction tc("D(d,e,c,a)+=L(a,b,c)*R(b,d,e)", H_4_legs_a.get(),
                        tensor_grid[1][0].get(), tensor_grid[2][0].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,c,e)+=L(a,b,c,d)*R(d,e)", H_4_legs_b.get(),
                        H_4_legs_a.get(), tensor_grid[0][0].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Group of 6
    // Group of 2
    tc = TensContraction("D(a,b,d,e)+=L(a,b,c)*R(c,d,e)", H_4_legs_a.get(),
                        tensor_grid[0][1].get(), tensor_grid[0][2].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Group of 4
    tc = TensContraction("D(a,b,c,f,g,e)+=L(a,b,c,d)*R(e,d,f,g)",
          H_6_legs_a.get(), tensor_grid[1][1].get(), tensor_grid[1][2].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,g,h,i,d,e,f)+=L(a,b,c,d,e,f)*R(c,g,h,i)",
                  H_8_legs_a.get(), H_6_legs_a.get(), tensor_grid[2][1].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
          "D(a,b,c,d,i,j,g,h)+=L(a,b,c,d,e,f,g,h)*R(f,e,i,j)",
          H_8_legs_b.get(), H_8_legs_a.get(), tensor_grid[2][2].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // 2 + 4
    tc = TensContraction(
          "D(a,e,f,g,h,i,j,d)+=L(a,b,c,d)*R(b,e,f,g,h,i,j,c)",
          H_8_legs_a.get(), H_4_legs_a.get(), H_8_legs_b.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // 3 + 6 qubits
    tc = TensContraction(
          "D(a,e,f,g,h,i)+=L(a,b,c,d)*R(d,c,b,e,f,g,h,i)",
          H_6_legs_c.get(), H_4_legs_b.get(), H_8_legs_a.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // 6_c is upper-left 9 qubits. All others are clean.
    // Group of 6
    // Group of 2
    tc = TensContraction("D(a,b,d,e)+=L(a,b,c)*R(c,d,e)", H_4_legs_a.get(),
                        tensor_grid[0][3].get(), tensor_grid[0][4].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Group of 4
    tc = TensContraction("D(a,b,c,f,g,e)+=L(a,b,c,d)*R(e,d,f,g)",
          H_6_legs_a.get(), tensor_grid[1][3].get(), tensor_grid[1][4].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,g,h,i,d,e,f)+=L(a,b,c,d,e,f)*R(c,g,h,i)",
                  H_8_legs_a.get(), H_6_legs_a.get(), tensor_grid[2][3].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
          "D(a,b,c,d,i,j,g,h)+=L(a,b,c,d,e,f,g,h)*R(f,e,i,j)",
          H_8_legs_b.get(), H_8_legs_a.get(), tensor_grid[2][4].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // 2 + 4
    tc = TensContraction(
          "D(a,e,f,g,h,i,j,d)+=L(a,b,c,d)*R(b,e,f,g,h,i,j,c)",
          H_8_legs_a.get(), H_4_legs_a.get(), H_8_legs_b.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // 6 + 8 -> 8
    tc = TensContraction(
          "D(a,b,c,g,h,i,j,k)+=L(a,b,c,d,e,f)*R(f,e,d,g,h,i,j,k)",
          H_8_legs_c.get(), H_6_legs_c.get(), H_8_legs_a.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));


    // 8_c is upper-left upper part. All other are clean.
    // Group of 6
    // Group of 2
    tc = TensContraction("D(d,e,c,a)+=L(a,b,c)*R(b,d,e)", H_4_legs_a.get(),
                        tensor_grid[3][0].get(), tensor_grid[4][0].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Group of 4
    tc = TensContraction("D(a,b,c,f,g,e)+=L(a,b,c,d)*R(e,d,f,g)",
          H_6_legs_a.get(), tensor_grid[3][1].get(), tensor_grid[3][2].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,g,h,i,d,e,f)+=L(a,b,c,d,e,f)*R(c,g,h,i)",
                  H_8_legs_a.get(), H_6_legs_a.get(), tensor_grid[4][1].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
          "D(a,b,c,d,j,i,g,h)+=L(a,b,c,d,e,f,g,h)*R(f,e,i,j)",
          H_8_legs_b.get(), H_8_legs_a.get(), tensor_grid[4][2].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // 2 + 4
    tc = TensContraction(
          "D(a,f,g,h,i,j,e,d)+=L(a,b,c,d)*R(e,c,b,f,g,h,i,j)",
          H_8_legs_a.get(), H_4_legs_a.get(), H_8_legs_b.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));

    // 15 + 6 -> 21
    errc = H_10_legs_a->contractAccumulateXL(nullptr,
          "D(i,j,k,l,m,d,e,f,g,h)+=L(a,b,c,d,e,f,g,h)*R(i,j,k,l,m,c,b,a)",
          *H_8_legs_c.get(),
          *H_8_legs_a.get(), DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}), false);
    done = H_10_legs_a->sync();
    // 10_a is all A but the last corner of 4. Everything else is clean.
    // Group of 4
    tc = TensContraction("D(a,b,c,f,g,e)+=L(a,b,c,d)*R(e,d,f,g)",
          H_6_legs_a.get(), tensor_grid[3][3].get(), tensor_grid[3][4].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,g,h,i,d,e,f)+=L(a,b,c,d,e,f)*R(c,g,h,i)",
                  H_8_legs_a.get(), H_6_legs_a.get(), tensor_grid[4][3].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
          "D(a,b,c,d,i,j,g,h)+=L(a,b,c,d,e,f,g,h)*R(f,e,i,j)",
          H_8_legs_b.get(), H_8_legs_a.get(), tensor_grid[4][4].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Finish A!
    errc = H_10_legs_c->contractAccumulateXL(nullptr,
          "D(a,b,c,k,l,m,n,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(f,e,d,k,l,m,n,g)",
          *H_10_legs_a.get(),
          *H_8_legs_b.get(), DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}), false);
    done = H_10_legs_c->sync();
    // Finished A!
    // 10_c is A. Everything else is clean.

    // B
    // Slice
    vector<int> dims_S_4_10(3, super_dim); dims_S_4_10[2] = 1;
    S_4_10->reshape(dims_S_4_10);
    talsh::TensorTask task_hl;
    errc = tensor_grid[4][10]->extractSlice(
                        &task_hl,*S_4_10,vector<int>{0,0,i0},DEV_HOST,0);
    done = tensor_grid[4][10]->sync();
    task_hl.clean();
    S_4_10->reshape(dims_2);
    // Group of 6
    // Group of 2
    tc = TensContraction("D(a,b,d)+=L(a,b,c)*R(c,d)", H_3_legs_a.get(),
                        tensor_grid[0][9].get(), tensor_grid[0][10].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Group of 4
    tc = TensContraction("D(a,b,c,f,e)+=L(a,b,c,d)*R(e,d,f)",
          H_5_legs_a.get(), tensor_grid[1][9].get(), tensor_grid[1][10].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,f,g,h,d,e)+=L(a,b,c,d,e)*R(c,f,g,h)",
                  H_7_legs_a.get(), H_5_legs_a.get(), tensor_grid[2][9].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
          "D(a,b,c,d,h,g)+=L(a,b,c,d,e,f,g)*R(f,e,h)",
          H_6_legs_a.get(), H_7_legs_a.get(), tensor_grid[2][10].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // 2+4 -> 6
    tc = TensContraction(
          "D(a,d,e,f,g)+=L(a,b,c)*R(b,d,e,f,g,c)",
          H_5_legs_c.get(), H_3_legs_a.get(), H_6_legs_a.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // 5_c is upper-right corner group of 6.
    // Second group of 4
    tc = TensContraction("D(a,b,c,f,e)+=L(a,b,c,d)*R(e,d,f)",
          H_5_legs_a.get(), tensor_grid[3][9].get(), tensor_grid[3][10].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,f,g,h,d,e)+=L(a,b,c,d,e)*R(c,f,g,h)",
                  H_7_legs_a.get(), H_5_legs_a.get(), tensor_grid[4][9].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
          "D(a,b,c,d,g)+=L(a,b,c,d,e,f,g)*R(f,e)",
          H_5_legs_a.get(), H_7_legs_a.get(), S_4_10.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Put together top-right corner with group of 4 bordering the cut
    tc = TensContraction(
          "D(a,b,c,f,g,h)+=L(a,b,c,d,e)*R(d,f,g,h,e)",
          H_6_legs_c.get(), H_5_legs_c.get(), H_5_legs_a.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // 6_c is the top-right right part of B. 10_c is taken, but all others
    // are clean
    // Next 2+4 from the right (to the left)
    // Group of 6
    // Group of 2
    tc = TensContraction("D(a,b,d,e)+=L(a,b,c)*R(c,d,e)", H_4_legs_a.get(),
                        tensor_grid[0][7].get(), tensor_grid[0][8].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Group of 4
    tc = TensContraction("D(a,b,c,f,g,e)+=L(a,b,c,d)*R(e,d,f,g)",
          H_6_legs_a.get(), tensor_grid[1][7].get(), tensor_grid[1][8].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,g,h,i,d,e,f)+=L(a,b,c,d,e,f)*R(c,g,h,i)",
                  H_8_legs_a.get(), H_6_legs_a.get(), tensor_grid[2][7].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
          "D(a,b,c,d,i,j,g,h)+=L(a,b,c,d,e,f,g,h)*R(f,e,i,j)",
          H_8_legs_b.get(), H_8_legs_a.get(), tensor_grid[2][8].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // 2 + 4
    tc = TensContraction(
          "D(a,e,f,g,h,i,j,d)+=L(a,b,c,d)*R(b,e,f,g,h,i,j,c)",
          H_8_legs_a.get(), H_4_legs_a.get(), H_8_legs_b.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Put this group of 6 together witht he previous group of 8
    tc = TensContraction(
          "D(a,b,c,d,e,i,j,k)+=L(a,b,c,d,e,f,g,h)*R(h,g,f,i,j,k)",
          H_8_legs_c.get(), H_8_legs_a.get(), H_6_legs_c.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // 8_c has what I have called alpha (top-right corner 14 qubits).
    // 10_c taken, others clean.
  
    // Next group of 4 (close to the cut)
    // Group of 4
    tc = TensContraction("D(a,b,c,f,g,e)+=L(a,b,c,d)*R(e,d,f,g)",
          H_6_legs_a.get(), tensor_grid[3][7].get(), tensor_grid[3][8].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,g,h,i,d,e,f)+=L(a,b,c,d,e,f)*R(c,g,h,i)",
                  H_8_legs_a.get(), H_6_legs_a.get(), tensor_grid[4][7].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
          "D(a,b,c,d,i,j,g,h)+=L(a,b,c,d,e,f,g,h)*R(f,e,i,j)",
          H_8_legs_b.get(), H_8_legs_a.get(), tensor_grid[4][8].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Put the previous corner of 14 together with this gorup of 4
    tc = TensContraction(
          "D(a,b,c,i,j,k,l,h)+=L(a,b,c,d,e,f,g,h)*R(d,i,j,k,l,g,f,e)",
          H_8_legs_d.get(), H_8_legs_c.get(), H_8_legs_b.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // 8_d is the top right corner of 20 qubits. 10_c taken, rest clean.
    // Group of 6
    // Group of 2
    tc = TensContraction("D(a,b,d,e)+=L(a,b,c)*R(c,d,e)", H_4_legs_a.get(),
                        tensor_grid[0][5].get(), tensor_grid[0][6].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Group of 4
    tc = TensContraction("D(a,b,c,f,g,e)+=L(a,b,c,d)*R(e,d,f,g)",
          H_6_legs_a.get(), tensor_grid[1][5].get(), tensor_grid[1][6].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,g,h,i,d,e,f)+=L(a,b,c,d,e,f)*R(c,g,h,i)",
                  H_8_legs_a.get(), H_6_legs_a.get(), tensor_grid[2][5].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
          "D(a,b,c,d,i,j,g,h)+=L(a,b,c,d,e,f,g,h)*R(f,e,i,j)",
          H_8_legs_b.get(), H_8_legs_a.get(), tensor_grid[2][6].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // 2 + 4
    tc = TensContraction(
          "D(a,e,f,g,h,i,j,d)+=L(a,b,c,d)*R(b,e,f,g,h,i,j,c)",
          H_8_legs_a.get(), H_4_legs_a.get(), H_8_legs_b.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Put it together with the 20 corner qubits
    errc = H_10_legs_a->contractAccumulateXL(nullptr,
          "D(a,b,c,d,e,i,j,k,l,m)+=L(a,b,c,d,e,f,g,h)*R(h,g,f,i,j,k,l,m)",
          *H_8_legs_a.get(),
          *H_8_legs_d.get(), DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}), false);
    done = H_10_legs_a->sync();
    // Final group of 4 for B
    // Group of 4
    tc = TensContraction("D(a,b,c,f,g,e)+=L(a,b,c,d)*R(e,d,f,g)",
          H_6_legs_a.get(), tensor_grid[3][5].get(), tensor_grid[3][6].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,g,h,i,d,e,f)+=L(a,b,c,d,e,f)*R(c,g,h,i)",
                  H_8_legs_a.get(), H_6_legs_a.get(), tensor_grid[4][5].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
          "D(a,b,c,d,i,j,g,h)+=L(a,b,c,d,e,f,g,h)*R(f,e,i,j)",
          H_8_legs_b.get(), H_8_legs_a.get(), tensor_grid[4][6].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Put it together with the 26 corner qubits
    errc = H_10_legs_b->contractAccumulateXL(nullptr,
          "D(a,b,c,k,l,m,n,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(d,k,l,m,n,g,f,e)",
          *H_10_legs_a.get(),
          *H_8_legs_b.get(), DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}), false);
    done = H_10_legs_b->sync();
    // 10_b is B!

    // Build AB (10_c & 10_b)
    // Time it!
    t0 = high_resolution_clock::now();
    errc = AB->contractAccumulateXL(nullptr,
      "D(a,b,c,d,e,k,l,m,n,o)+=L(a,b,c,d,e,f,g,h,i,j)*R(j,i,h,g,f,k,l,m,n,o)",
      *H_10_legs_c.get(),
      *H_10_legs_b.get(), DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}), false);
    done = AB->sync();
    t1 = high_resolution_clock::now();
    span = duration_cast<duration<double>>(t1 - t0);
    time_largest_contraction = double(span.count());


    // E
    // All clean but AB!
    // Slice
    vector<int> dims_S_5_10(3, super_dim); dims_S_5_10[0] = 1;
    S_5_10->reshape(dims_S_5_10);
    errc = tensor_grid[5][10]->extractSlice(
                        &task_hl,*S_5_10,vector<int>{i0,0,0},DEV_HOST,0);
    done = tensor_grid[5][10]->sync();
    task_hl.clean();
    S_5_10->reshape(dims_2);
    // Group of 4: bottom-right corner
    tc = TensContraction("D(a,b,c,f,e)+=L(a,b,c,d)*R(e,d,f)",
        H_5_legs_a.get(), tensor_grid[9][9].get(), tensor_grid[9][10].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,f,g,d,e)+=L(a,b,c,d,e)*R(c,f,g)",
                H_6_legs_a.get(), H_5_legs_a.get(), tensor_grid[10][9].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
          "D(a,b,c,f)+=L(a,b,c,d,e,f)*R(e,d)",
          H_4_legs_b.get(), H_6_legs_a.get(), tensor_grid[10][10].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // 4_b taken
    // Next group of 4 (going upwards)
    tc = TensContraction("D(a,b,c,f,e)+=L(a,b,c,d)*R(e,d,f)",
          H_5_legs_a.get(), tensor_grid[7][9].get(), tensor_grid[7][10].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,f,g,h,d,e)+=L(a,b,c,d,e)*R(c,f,g,h)",
                  H_7_legs_a.get(), H_5_legs_a.get(), tensor_grid[8][9].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
          "D(a,b,c,d,h,g)+=L(a,b,c,d,e,f,g)*R(f,e,h)",
          H_6_legs_a.get(), H_7_legs_a.get(), tensor_grid[8][10].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // 4+4 -> 8 qubits on bottom-right corner
    tc = TensContraction(
          "D(f,a,b,c,g,h)+=L(a,b,c,d,e,f)*R(d,g,h,e)",
          H_6_legs_c.get(), H_6_legs_a.get(), H_4_legs_b.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // 6_c is bottom-right 8 qubits. All others but AB clean!
    // Group of 4 by the cut
    tc = TensContraction("D(a,b,c,e)+=L(a,b,c,d)*R(d,e)",
          H_4_legs_a.get(), tensor_grid[5][9].get(), S_5_10.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,e,f,g,d)+=L(a,b,c,d)*R(c,e,f,g)",
                  H_6_legs_a.get(), H_4_legs_a.get(), tensor_grid[6][9].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
          "D(a,b,c,d,g)+=L(a,b,c,d,e,f)*R(f,e,g)",
          H_5_legs_a.get(), H_6_legs_a.get(), tensor_grid[6][10].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // 8 + 4 -> bottom-right 12 qubits
    tc = TensContraction("D(g,h,i,c,d,e,f)+=L(a,b,c,d,e,f)*R(g,h,i,b,a)",
          H_7_legs_b.get(), H_6_legs_c.get(), H_5_legs_a.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // 7_b has the bottom-right right 12 qubits! All others clean but AB
    // Next group of 4 to the left of the bottom-left corner
    tc = TensContraction("D(a,b,c,f,g,e)+=L(a,b,c,d)*R(e,d,f,g)",
          H_6_legs_a.get(), tensor_grid[9][7].get(), tensor_grid[9][8].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,h,i,d,e,f)+=L(a,b,c,d,e,f)*R(c,h,i)",
                H_7_legs_a.get(), H_6_legs_a.get(), tensor_grid[10][7].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
          "D(a,b,c,h,f,g)+=L(a,b,c,d,e,f,g)*R(e,d,h)",
          H_6_legs_a.get(), H_7_legs_a.get(), tensor_grid[10][8].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // This 4 with the bottom-right right 12 qubits -> 16
    tc = TensContraction(
          "D(a,b,c,d,e,k,h,i,j)+=L(a,b,c,d,e,f,g)*R(h,i,j,g,f,k)",
          H_9_legs_a.get(), H_7_legs_b.get(), H_6_legs_a.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Next group of 4 (left-right of the bottom-right corner)
    tc = TensContraction("D(a,b,c,f,g,e)+=L(a,b,c,d)*R(e,d,f,g)",
          H_6_legs_a.get(), tensor_grid[7][7].get(), tensor_grid[7][8].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,g,h,i,d,e,f)+=L(a,b,c,d,e,f)*R(c,g,h,i)",
                  H_8_legs_a.get(), H_6_legs_a.get(), tensor_grid[8][7].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
          "D(a,b,c,d,i,j,g,h)+=L(a,b,c,d,e,f,g,h)*R(f,e,i,j)",
          H_8_legs_b.get(), H_8_legs_a.get(), tensor_grid[8][8].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // 4+16 -> 20
    tc = TensContraction(
          "D(a,b,c,m,j,k,l,h,i)+=L(a,b,c,d,e,f,g,h,i)*R(j,k,l,g,f,e,d,m)",
          H_9_legs_b.get(), H_9_legs_a.get(), H_8_legs_b.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Last group of 4 for pE
    tc = TensContraction("D(a,b,c,f,g,e)+=L(a,b,c,d)*R(e,d,f,g)",
          H_6_legs_a.get(), tensor_grid[5][7].get(), tensor_grid[5][8].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction("D(a,b,g,h,i,d,e,f)+=L(a,b,c,d,e,f)*R(c,g,h,i)",
                  H_8_legs_a.get(), H_6_legs_a.get(), tensor_grid[6][7].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    tc = TensContraction(
          "D(a,b,c,d,i,j,g,h)+=L(a,b,c,d,e,f,g,h)*R(f,e,i,j)",
          H_8_legs_b.get(), H_8_legs_a.get(), tensor_grid[6][8].get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Finish pE (24 qubits)
    tc = TensContraction(
          "D(a,m,j,k,l,f,g,h,i)+=L(a,b,c,d,e,f,g,h,i)*R(j,k,l,e,d,c,b,m)",
          pE.get(), H_9_legs_b.get(), H_8_legs_b.get());
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    // Built pE!


    
    // Finish building E, then ABE, and start completeng C!
    // Only AB and pE are taken. All others are clean!
    vector<int> inner_values({2});
    //vector<int> inner_values({0,1,2,3});
    for (auto i1 : inner_values)
    //for (int i1=0; i1<1; ++i1)
    {
      // Slices
      vector<int> dims_S_10_5(3, super_dim); dims_S_10_5[1] = 1;
      S_10_5->reshape(dims_S_10_5);
      errc = tensor_grid[10][5]->extractSlice(
                        &task_hl,*S_10_5,vector<int>{0,i1,0},DEV_HOST,0);
      done = tensor_grid[10][5]->sync();
      task_hl.clean();
      S_10_5->reshape(dims_2);
      //
      vector<int> dims_S_10_4(3, super_dim); dims_S_10_4[2] = 1;
      S_10_4->reshape(dims_S_10_4);
      errc = tensor_grid[10][4]->extractSlice(
                        &task_hl,*S_10_4,vector<int>{0,0,i1},DEV_HOST,0);
      done = tensor_grid[10][4]->sync();
      task_hl.clean();
      S_10_4->reshape(dims_2);
      // Group of 4 by the cut
      tc = TensContraction("D(a,b,c,f,g,e)+=L(a,b,c,d)*R(e,d,f,g)",
          H_6_legs_a.get(), tensor_grid[9][5].get(), tensor_grid[9][6].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST));
      tc = TensContraction("D(a,b,g,d,e,f)+=L(a,b,c,d,e,f)*R(c,g)",
                    H_6_legs_b.get(), H_6_legs_a.get(), S_10_5.get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST));
      tc = TensContraction(
            "D(a,b,g,e,f)+=L(a,b,c,d,e,f)*R(d,c,g)",
            H_5_legs_a.get(), H_6_legs_b.get(), tensor_grid[10][6].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST));
      // Put it together with pE
      errc = H_10_legs_a->contractAccumulateXL(nullptr,
          "D(a,b,c,d,e,f,g,l,j,k)+=L(a,b,c,d,e,f,g,h,i)*R(j,k,i,h,l)",
          *pE,
          *H_5_legs_a, DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}), false);
      done = H_10_legs_a->sync();
      assert(tc.sync(DEV_HOST,0));
      // Second to last group of 4 for E
      tc = TensContraction("D(a,b,c,f,g,e)+=L(a,b,c,d)*R(e,d,f,g)",
          H_6_legs_a.get(), tensor_grid[7][5].get(), tensor_grid[7][6].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,g,h,i,d,e,f)+=L(a,b,c,d,e,f)*R(c,g,h,i)",
                  H_8_legs_a.get(), H_6_legs_a.get(), tensor_grid[8][5].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction(
            "D(a,b,c,d,i,j,g,h)+=L(a,b,c,d,e,f,g,h)*R(f,e,i,j)",
            H_8_legs_b.get(), H_8_legs_a.get(), tensor_grid[8][6].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      // Put it together: 28+4 -> 32 qubits
      errc = H_10_legs_b->contractAccumulateXL(nullptr,
          "D(a,b,c,d,e,n,k,l,m,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(k,l,m,i,h,g,f,n)",
          *H_10_legs_a,
          *H_8_legs_b, DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}), false);
      done = H_10_legs_b->sync();
      assert(tc.sync(DEV_HOST,0));
      // Last group of 4 for E!
      tc = TensContraction("D(a,b,c,f,g,e)+=L(a,b,c,d)*R(e,d,f,g)",
          H_6_legs_a.get(), tensor_grid[5][5].get(), tensor_grid[5][6].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,g,h,i,d,e,f)+=L(a,b,c,d,e,f)*R(c,g,h,i)",
                  H_8_legs_a.get(), H_6_legs_a.get(), tensor_grid[6][5].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction(
            "D(a,b,c,d,i,j,g,h)+=L(a,b,c,d,e,f,g,h)*R(f,e,i,j)",
            H_8_legs_b.get(), H_8_legs_a.get(), tensor_grid[6][6].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      // Finish E!
      errc = H_10_legs_a->contractAccumulateXL(nullptr,
          "D(a,b,c,n,k,l,m,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(k,l,m,g,f,e,d,n)",
          *H_10_legs_b,
          *H_8_legs_b, DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}), false);
      done = H_10_legs_a->sync();
      assert(tc.sync(DEV_HOST,0));
      // Finished E!

      // AB + 10_a -> ABD (10_b)
      errc = H_10_legs_b->contractAccumulateXL(nullptr,
      "D(a,b,c,d,e,k,l,m,n,o)+=L(a,b,c,d,e,f,g,h,i,j)*R(j,i,h,g,f,k,l,m,n,o)",
        *AB,
        *H_10_legs_a, DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}), false);
      done = H_10_legs_b->sync();
      assert(tc.sync(DEV_HOST,0));

      // Start C
      // Inner most group of 4
      tc = TensContraction("D(a,b,c,f,g,e)+=L(a,b,c,d)*R(e,d,f,g)",
          H_6_legs_a.get(), tensor_grid[5][3].get(), tensor_grid[5][4].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,g,h,i,d,e,f)+=L(a,b,c,d,e,f)*R(c,g,h,i)",
                  H_8_legs_a.get(), H_6_legs_a.get(), tensor_grid[6][3].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction(
            "D(a,b,c,d,i,j,g,h)+=L(a,b,c,d,e,f,g,h)*R(f,e,i,j)",
            H_8_legs_b.get(), H_8_legs_a.get(), tensor_grid[6][4].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      // Contract with ABD
      errc = H_10_legs_a->contractAccumulateXL(nullptr,
          "D(a,b,c,k,l,m,n,h,i,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(d,k,l,m,n,g,f,e)",
          *H_10_legs_b,
          *H_8_legs_b, DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}), false);
      done = H_10_legs_a->sync();
      assert(tc.sync(DEV_HOST,0));
      // Next group of 4 downwards
      tc = TensContraction("D(a,b,c,f,g,e)+=L(a,b,c,d)*R(e,d,f,g)",
          H_6_legs_a.get(), tensor_grid[7][3].get(), tensor_grid[7][4].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,g,h,i,d,e,f)+=L(a,b,c,d,e,f)*R(c,g,h,i)",
                  H_8_legs_a.get(), H_6_legs_a.get(), tensor_grid[8][3].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction(
            "D(a,b,c,d,i,j,g,h)+=L(a,b,c,d,e,f,g,h)*R(f,e,i,j)",
            H_8_legs_b.get(), H_8_legs_a.get(), tensor_grid[8][4].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      // Contract with the rest
      errc = H_10_legs_b->contractAccumulateXL(nullptr,
          "D(a,b,c,d,e,k,l,m,n,j)+=L(a,b,c,d,e,f,g,h,i,j)*R(f,k,l,m,n,i,h,g)",
          *H_10_legs_a,
          *H_8_legs_b, DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}), false);
      done = H_10_legs_b->sync();
      assert(tc.sync(DEV_HOST,0));
      // Group of 4 by the cut
      tc = TensContraction("D(a,b,c,f,g,e)+=L(a,b,c,d)*R(e,d,f,g)",
          H_6_legs_a.get(), tensor_grid[9][3].get(), tensor_grid[9][4].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,g,h,d,e,f)+=L(a,b,c,d,e,f)*R(c,g,h)",
                H_7_legs_a.get(), H_6_legs_a.get(), tensor_grid[10][3].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction(
            "D(a,b,c,f,g)+=L(a,b,c,d,e,f,g)*R(e,d)",
            H_5_legs_a.get(), H_7_legs_a.get(), S_10_4.get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      // Contract with the rest (start reducing size)
      errc = H_9_legs_a->contractAccumulateXL(nullptr,
          "D(a,b,c,d,e,f,g,k,l)+=L(a,b,c,d,e,f,g,h,i,j)*R(h,k,l,j,i)",
          *H_10_legs_b,
          *H_5_legs_a, DEV_NVIDIA_GPU, 0, s_type({1.0,0.0}), false);
      done = H_9_legs_a->sync();
      assert(tc.sync(DEV_HOST,0));
      // Next group of 4 (3 to go)
      tc = TensContraction("D(a,b,c,f,g,e)+=L(a,b,c,d)*R(e,d,f,g)",
          H_6_legs_a.get(), tensor_grid[5][1].get(), tensor_grid[5][2].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,g,h,i,d,e,f)+=L(a,b,c,d,e,f)*R(c,g,h,i)",
                  H_8_legs_a.get(), H_6_legs_a.get(), tensor_grid[6][1].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction(
            "D(a,b,c,d,i,j,g,h)+=L(a,b,c,d,e,f,g,h)*R(f,e,i,j)",
            H_8_legs_b.get(), H_8_legs_a.get(), tensor_grid[6][2].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      // Contract with the rest
      tc = TensContraction(
            "D(a,j,k,l,m,f,g,h,i)+=L(a,b,c,d,e,f,g,h,i)*R(b,j,k,l,m,e,d,c)",
            H_9_legs_b.get(), H_9_legs_a.get(), H_8_legs_b.get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      // Next group of 4 (2 to go)
      tc = TensContraction("D(a,b,c,f,g,e)+=L(a,b,c,d)*R(e,d,f,g)",
          H_6_legs_a.get(), tensor_grid[7][1].get(), tensor_grid[7][2].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,g,h,i,d,e,f)+=L(a,b,c,d,e,f)*R(c,g,h,i)",
                  H_8_legs_a.get(), H_6_legs_a.get(), tensor_grid[8][1].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction(
            "D(a,b,c,d,i,j,g,h)+=L(a,b,c,d,e,f,g,h)*R(f,e,i,j)",
            H_8_legs_b.get(), H_8_legs_a.get(), tensor_grid[8][2].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      // Contract with the rest
      tc = TensContraction(
            "D(a,b,c,j,k,l,m,h,i)+=L(a,b,c,d,e,f,g,h,i)*R(d,j,k,l,m,g,f,e)",
            H_9_legs_a.get(), H_9_legs_b.get(), H_8_legs_b.get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      // Last group of 4
      tc = TensContraction("D(a,b,c,f,g,e)+=L(a,b,c,d)*R(e,d,f,g)",
          H_6_legs_a.get(), tensor_grid[9][1].get(), tensor_grid[9][2].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction("D(a,b,g,h,d,e,f)+=L(a,b,c,d,e,f)*R(c,g,h)",
                H_7_legs_a.get(), H_6_legs_a.get(), tensor_grid[10][1].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      tc = TensContraction(
            "D(a,b,c,h,f,g)+=L(a,b,c,d,e,f,g)*R(e,d,h)",
            H_6_legs_a.get(), H_7_legs_a.get(), tensor_grid[10][2].get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      // Contract with the rest
      tc = TensContraction(
            "D(a,b,c,d,e,j,k)+=L(a,b,c,d,e,f,g,h,i)*R(f,j,k,i,h,g)",
            H_7_legs_b.get(), H_9_legs_a.get(), H_6_legs_a.get());
      errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(tc.sync(DEV_HOST,0));
      // 7_b is circuit-C. All others are free


      // C loop!
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
          if (t!=5)
          {
            tc = TensContraction("D(b,c,d)+=L(a)*R(a,b,c,d)",
                                Cs[t].get(), &delta, tensor_grid[i][j].get());
          } else {
            tc = TensContraction("D(b,c)+=L(a)*R(a,b,c)",
                                Cs[t].get(), &delta, tensor_grid[i][j].get());
          }
          errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
          assert(tc.sync(DEV_HOST,0));
        }
        // Contract C
        tc = TensContraction("D(c,d,b)+=L(a,b)*R(c,a,d)",
                      H_3_legs_a.get(), Cs[5].get(), Cs[4].get());
        errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
        assert(tc.sync(DEV_HOST,0));
        tc = TensContraction("D(d,e,b,c)+=L(a,b,c)*R(d,a,e)",
                      H_4_legs_a.get(), H_3_legs_a.get(), Cs[3].get());
        errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
        assert(tc.sync(DEV_HOST,0));
        tc = TensContraction("D(e,f,b,c,d)+=L(a,b,c,d)*R(e,a,f)",
                      H_5_legs_a.get(), H_4_legs_a.get(), Cs[2].get());
        errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
        assert(tc.sync(DEV_HOST,0));
        tc = TensContraction("D(f,g,b,c,d,e)+=L(a,b,c,d,e)*R(f,a,g)",
                      H_6_legs_a.get(), H_5_legs_a.get(), Cs[1].get());
        errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
        assert(tc.sync(DEV_HOST,0));
        tc = TensContraction("D(g,h,b,c,d,e,f)+=L(a,b,c,d,e,f)*R(g,a,h)",
                      H_7_legs_a.get(), H_6_legs_a.get(), Cs[0].get());
        errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
        assert(tc.sync(DEV_HOST,0));
        // Contract ABC (7_b) and C (7_a)
        tc = TensContraction("D()+=L(a,b,c,d,e,f,g)*R(a,b,c,d,e,f,g)",
                            S.get(), H_7_legs_b.get(), H_7_legs_a.get());
        errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
        assert(tc.sync(DEV_HOST,0));
        // Store result data_S
        s_type const * ptr_S;
        S->getDataAccessHostConst(&ptr_S);
        //amplitudes[c] += ptr_S[0];
        amplitudes[c] += ptr_S[0] / global_norm_factor;
        ptr_S = nullptr;
      }
      t1 = high_resolution_clock::now();
      span = duration_cast<duration<double>>(t1 - t0);
      //cout << "Time contracting all Cs is " << span.count() << endl;
    }

  }
  */

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
