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

#include "contraction_sycamore_53x10_cphase.h"
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
  // Corners ((5,0) and (9,4))
  unique_ptr<talsh::Tensor> corner_51(new talsh::Tensor({dims_grid[5][1][0],
                                                         dims_grid[5][1][2],
                                                         dims_grid[5][1][3]},
                                                         s_type(0.0)));
  errc = corner_51->contractAccumulate(nullptr,
        "D(b,c,d)+=L(a)*R(b,a,c,d)",
        *tensor_grid[5][0].get(),
        *tensor_grid[5][1].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = corner_51->sync();
  unique_ptr<talsh::Tensor> corner_84(new talsh::Tensor({dims_grid[8][4][0],
                                                         dims_grid[8][4][1],
                                                         dims_grid[8][4][3]},
                                                         s_type(0.0)));
  errc = corner_84->contractAccumulate(nullptr,
        "D(a,b,d)+=L(a,b,c,d)*R(c)",
        *tensor_grid[8][4].get(),
        *tensor_grid[9][4].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = corner_84->sync();
  // Column 1
  unique_ptr<talsh::Tensor> aa(new talsh::Tensor({dims_grid[5][1][2],
                                                  dims_grid[5][1][3],
                                                  dims_grid[4][1][1]},
                                                  s_type(0.0)));
  errc = aa->contractAccumulate(nullptr,
        "D(c,d,b)+=L(a,b)*R(a,c,d)",
        *tensor_grid[4][1].get(),
        *corner_51.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
  done = aa->sync();
  unique_ptr<talsh::Tensor> ab(new talsh::Tensor({dims_grid[6][1][1],
                                                  dims_grid[5][1][3],
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
                                                  dims_grid[5][1][3],
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
                                                  dims_grid[8][4][3]},
                                                  s_type(0.0)));
  errc = ba->contractAccumulate(nullptr,
        "D(a,c,d)+=L(a,b)*R(c,b,d)",
        *tensor_grid[8][3].get(),
        *corner_84.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
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
  // Group of 2 at the end of row 7
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
  int first_cut_dim(dims_grid[3][3][2]);
  //vector<int> first_cut_offsets({0});
  vector<int> first_cut_offsets;
  for (int p=0; p<first_cut_dim; ++p)
  {
    first_cut_offsets.push_back(p);
  }
  int first_slice_size(1);
  for (auto i0 : first_cut_offsets)
    cout << i0 << " ";
  cout << "\n" << flush;
  for (auto i0 : first_cut_offsets)
  {
    // First slice on [3][3] and [4][3].
    unique_ptr<talsh::Tensor> first_3_3(new talsh::Tensor({dims_grid[3][3][0],
                                                         dims_grid[3][3][1],
                                                         first_slice_size},
                                                         s_type(0.0)));
    unique_ptr<talsh::Tensor> first_3_4(new talsh::Tensor({dims_grid[3][4][0],
                                                         first_slice_size,
                                                         dims_grid[3][4][2],
                                                         dims_grid[3][4][3]},
                                                         s_type(0.0)));
    errc = tensor_grid[3][3]->extractSlice(
                  &task_hl,*first_3_3,vector<int>{0,0,i0},DEV_HOST,0);
    done = tensor_grid[3][3]->sync();
    task_hl.clean();
    errc = tensor_grid[3][4]->extractSlice(
                  &task_hl,*first_3_4,vector<int>{0,i0,0,0},DEV_HOST,0);
    done = tensor_grid[3][4]->sync();
    task_hl.clean();

    // Continue contraction
    // Finish A
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
    unique_ptr<talsh::Tensor> aj(new talsh::Tensor({dims_grid[7][2][1],
                                                    dims_grid[6][2][3],
                                                    dims_grid[5][3][2],
                                                    dims_grid[5][3][3],
                                                    dims_grid[4][3][3],
                                                    first_slice_size},
                                                    s_type(0.0)));
    errc = aj->contractAccumulate(nullptr,
          "D(a,b,g,h,e,f)+=L(a,b,c,d,e,f)*R(d,c,g,h)",
          *ai.get(),
          *tensor_grid[5][3].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = aj->sync();
    ai.reset();
    unique_ptr<talsh::Tensor> ak(new talsh::Tensor({dims_grid[7][2][1],
                                                    dims_grid[6][3][2],
                                                    dims_grid[6][3][3],
                                                    dims_grid[5][3][3],
                                                    dims_grid[4][3][3],
                                                    first_slice_size},
                                                    s_type(0.0)));
    errc = ak->contractAccumulate(nullptr,
          "D(a,g,h,d,e,f)+=L(a,b,c,d,e,f)*R(c,b,g,h)",
          *aj.get(),
          *tensor_grid[6][3].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = ak->sync();
    aj.reset();
    unique_ptr<talsh::Tensor> A(new talsh::Tensor({dims_grid[7][3][2],
                                                    dims_grid[7][3][3],
                                                    dims_grid[6][3][3],
                                                    dims_grid[5][3][3],
                                                    dims_grid[4][3][3],
                                                    first_slice_size},
                                                    s_type(0.0)));
    errc = A->contractAccumulate(nullptr,
          "D(g,h,c,d,e,f)+=L(a,b,c,d,e,f)*R(b,a,g,h)",
          *ak.get(),
          *tensor_grid[7][3].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = A->sync();
    ak.reset();
    // al is A

    // Contract A with B
    unique_ptr<talsh::Tensor> AB(new talsh::Tensor({first_slice_size,
                                                    dims_grid[4][3][3],
                                                    dims_grid[5][3][3],
                                                    dims_grid[6][4][0],
                                                    dims_grid[6][5][0],
                                                    dims_grid[6][6][0],
                                                    dims_grid[6][7][0]},
                                                    s_type(0.0)));
    errc = AB->contractAccumulate(nullptr,
          "D(f,e,d,g,h,i,j)+=L(a,b,c,d,e,f)*R(a,b,c,g,h,i,j)",
          *A.get(),
          *B.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = AB->sync();
    A.reset();

    // Contract AB with further tensors onto region E
    // First row of E (row 5)
    unique_ptr<talsh::Tensor> ea(new talsh::Tensor({first_slice_size,
                                                    dims_grid[4][3][3],
                                                    dims_grid[5][4][0],
                                                    dims_grid[5][4][3],
                                                    dims_grid[6][5][0],
                                                    dims_grid[6][6][0],
                                                    dims_grid[6][7][0]},
                                                    s_type(0.0)));
    errc = ea->contractAccumulate(nullptr,
          "D(a,b,h,i,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,c,d,i)",
          *AB.get(),
          *tensor_grid[5][4].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = ea->sync();
    AB.reset();
    unique_ptr<talsh::Tensor> eb(new talsh::Tensor({first_slice_size,
                                                    dims_grid[4][3][3],
                                                    dims_grid[5][4][0],
                                                    dims_grid[5][5][0],
                                                    dims_grid[5][5][3],
                                                    dims_grid[6][6][0],
                                                    dims_grid[6][7][0]},
                                                    s_type(0.0)));
    errc = eb->contractAccumulate(nullptr,
          "D(a,b,c,h,i,f,g)+=L(a,b,c,d,e,f,g)*R(h,d,e,i)",
          *ea.get(),
          *tensor_grid[5][5].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = eb->sync();
    ea.reset();
    unique_ptr<talsh::Tensor> ec(new talsh::Tensor({first_slice_size,
                                                    dims_grid[4][3][3],
                                                    dims_grid[5][4][0],
                                                    dims_grid[5][5][0],
                                                    dims_grid[5][6][0],
                                                    dims_grid[5][6][3],
                                                    dims_grid[6][7][0]},
                                                    s_type(0.0)));
    errc = ec->contractAccumulate(nullptr,
          "D(a,b,c,d,h,i,g)+=L(a,b,c,d,e,f,g)*R(h,e,f,i)",
          *eb.get(),
          *tensor_grid[5][6].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = ec->sync();
    eb.reset();
    unique_ptr<talsh::Tensor> ed(new talsh::Tensor({first_slice_size,
                                                    dims_grid[4][3][3],
                                                    dims_grid[5][4][0],
                                                    dims_grid[5][5][0],
                                                    dims_grid[5][6][0],
                                                    dims_grid[5][7][0],
                                                    dims_grid[5][7][3]},
                                                    s_type(0.0)));
    errc = ed->contractAccumulate(nullptr,
          "D(a,b,c,d,e,h,i)+=L(a,b,c,d,e,f,g)*R(h,f,g,i)",
          *ec.get(),
          *tensor_grid[5][7].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = ed->sync();
    ec.reset();
    // Second row of E (row 4)
    unique_ptr<talsh::Tensor> ee(new talsh::Tensor({first_slice_size,
                                                    dims_grid[4][4][0],
                                                    dims_grid[4][4][3],
                                                    dims_grid[5][5][0],
                                                    dims_grid[5][6][0],
                                                    dims_grid[5][7][0],
                                                    dims_grid[5][7][3]},
                                                    s_type(0.0)));
    errc = ee->contractAccumulate(nullptr,
          "D(a,h,i,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,b,c,i)",
          *ed.get(),
          *tensor_grid[4][4].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = ee->sync();
    ed.reset();
    unique_ptr<talsh::Tensor> ef(new talsh::Tensor({first_slice_size,
                                                    dims_grid[4][4][0],
                                                    dims_grid[4][5][0],
                                                    dims_grid[4][5][3],
                                                    dims_grid[5][6][0],
                                                    dims_grid[5][7][0],
                                                    dims_grid[5][7][3]},
                                                    s_type(0.0)));
    errc = ef->contractAccumulate(nullptr,
          "D(a,b,h,i,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,c,d,i)",
          *ee.get(),
          *tensor_grid[4][5].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = ef->sync();
    ee.reset();
    unique_ptr<talsh::Tensor> eg(new talsh::Tensor({first_slice_size,
                                                    dims_grid[4][4][0],
                                                    dims_grid[4][5][0],
                                                    dims_grid[4][6][0],
                                                    dims_grid[4][6][3],
                                                    dims_grid[5][7][0],
                                                    dims_grid[5][7][3]},
                                                    s_type(0.0)));
    errc = eg->contractAccumulate(nullptr,
          "D(a,b,c,h,i,f,g)+=L(a,b,c,d,e,f,g)*R(h,d,e,i)",
          *ef.get(),
          *tensor_grid[4][6].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = eg->sync();
    ef.reset();
    unique_ptr<talsh::Tensor> eh(new talsh::Tensor({first_slice_size,
                                                    dims_grid[4][4][0],
                                                    dims_grid[4][5][0],
                                                    dims_grid[4][6][0],
                                                    dims_grid[4][7][0],
                                                    dims_grid[4][7][3],
                                                    dims_grid[5][7][3]},
                                                    s_type(0.0)));
    errc = eh->contractAccumulate(nullptr,
          "D(a,b,c,d,h,i,g)+=L(a,b,c,d,e,f,g)*R(h,e,f,i)",
          *eg.get(),
          *tensor_grid[4][7].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = eh->sync();
    eg.reset();

    // Contract partial E with D
    unique_ptr<talsh::Tensor> ei(new talsh::Tensor({first_slice_size,
                                                    dims_grid[4][4][0],
                                                    dims_grid[4][5][0],
                                                    dims_grid[4][6][0],
                                                    dims_grid[4][7][0],
                                                    dims_grid[3][8][1],
                                                    dims_grid[2][8][0]},
                                                    s_type(0.0)));
    errc = ei->contractAccumulate(nullptr,
          "D(a,b,c,d,e,h,i)+=L(a,b,c,d,e,f,g)*R(g,f,h,i)",
          *eh.get(),
          *D.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = ei->sync();
    eh.reset();

    // Continue contracting
    // Third row of E (row 3)
    unique_ptr<talsh::Tensor> ej(new talsh::Tensor({first_slice_size,
                                                    dims_grid[4][4][0],
                                                    dims_grid[4][5][0],
                                                    dims_grid[4][6][0],
                                                    dims_grid[3][7][1],
                                                    dims_grid[3][7][0],
                                                    dims_grid[2][8][0]},
                                                    s_type(0.0)));
    errc = ej->contractAccumulate(nullptr,
          "D(a,b,c,d,i,h,g)+=L(a,b,c,d,e,f,g)*R(h,i,e,f)",
          *ei.get(),
          *tensor_grid[3][7].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = ej->sync();
    ei.reset();
    unique_ptr<talsh::Tensor> ek(new talsh::Tensor({first_slice_size,
                                                    dims_grid[4][4][0],
                                                    dims_grid[4][5][0],
                                                    dims_grid[3][6][1],
                                                    dims_grid[3][6][0],
                                                    dims_grid[3][7][0],
                                                    dims_grid[2][8][0]},
                                                    s_type(0.0)));
    errc = ek->contractAccumulate(nullptr,
          "D(a,b,c,i,h,f,g)+=L(a,b,c,d,e,f,g)*R(h,i,d,e)",
          *ej.get(),
          *tensor_grid[3][6].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = ek->sync();
    ej.reset();
    unique_ptr<talsh::Tensor> el(new talsh::Tensor({first_slice_size,
                                                    dims_grid[4][4][0],
                                                    dims_grid[3][5][1],
                                                    dims_grid[3][5][0],
                                                    dims_grid[3][6][0],
                                                    dims_grid[3][7][0],
                                                    dims_grid[2][8][0]},
                                                    s_type(0.0)));
    errc = el->contractAccumulate(nullptr,
          "D(a,b,i,h,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,i,c,d)",
          *ek.get(),
          *tensor_grid[3][5].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = el->sync();
    ek.reset();
    unique_ptr<talsh::Tensor> em(new talsh::Tensor({dims_grid[3][4][0],
                                                    dims_grid[3][5][0],
                                                    dims_grid[3][6][0],
                                                    dims_grid[3][7][0],
                                                    dims_grid[2][8][0]},
                                                    s_type(0.0)));
    errc = em->contractAccumulate(nullptr,
          "D(h,d,e,f,g)+=L(a,b,c,d,e,f,g)*R(h,a,b,c)",
          *el.get(),
          *first_3_4.get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = em->sync();
    el.reset();

    // Fourth (and last) row of E (row 2)
    unique_ptr<talsh::Tensor> en(new talsh::Tensor({dims_grid[3][4][0],
                                                    dims_grid[3][5][0],
                                                    dims_grid[3][6][0],
                                                    dims_grid[2][7][1],
                                                    dims_grid[2][7][0]},
                                                    s_type(0.0)));
    errc = en->contractAccumulate(nullptr,
          "D(a,b,c,g,f)+=L(a,b,c,d,e)*R(f,g,d,e)",
          *em.get(),
          *tensor_grid[2][7].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = en->sync();
    em.reset();
    unique_ptr<talsh::Tensor> eo(new talsh::Tensor({dims_grid[3][4][0],
                                                    dims_grid[3][5][0],
                                                    dims_grid[2][6][1],
                                                    dims_grid[2][6][0],
                                                    dims_grid[2][7][0]},
                                                    s_type(0.0)));
    errc = eo->contractAccumulate(nullptr,
          "D(a,b,g,f,e)+=L(a,b,c,d,e)*R(f,g,c,d)",
          *en.get(),
          *tensor_grid[2][6].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = eo->sync();
    en.reset();
    unique_ptr<talsh::Tensor> ep(new talsh::Tensor({dims_grid[3][4][0],
                                                    dims_grid[2][5][1],
                                                    dims_grid[2][5][0],
                                                    dims_grid[2][6][0],
                                                    dims_grid[2][7][0]},
                                                    s_type(0.0)));
    errc = ep->contractAccumulate(nullptr,
          "D(a,g,f,d,e)+=L(a,b,c,d,e)*R(f,g,b,c)",
          *eo.get(),
          *tensor_grid[2][5].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = ep->sync();
    eo.reset();
    unique_ptr<talsh::Tensor> eq(new talsh::Tensor({dims_grid[2][4][0],
                                                    dims_grid[2][5][0],
                                                    dims_grid[2][6][0],
                                                    dims_grid[2][7][0]},
                                                    s_type(0.0)));
    errc = eq->contractAccumulate(nullptr,
          "D(f,c,d,e)+=L(a,b,c,d,e)*R(f,a,b)",
          *ep.get(),
          *tensor_grid[2][4].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
    done = eq->sync();
    ep.reset();
    // eq is "everything but C"


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
        string contraction_string;
        if (t==0 || t==1 || t==2 || t==5)
        {
          contraction_string = "D(b,c)+=L(a)*R(a,b,c)";
        } else {
          contraction_string = "D(b,c,d,e)+=L(a)*R(a,b,c,d,e)";
        }
        errc = Cs[t]->contractAccumulate(nullptr,
              contraction_string,
              delta,
              *tensor_grid[i][j].get(),
              DEV_NVIDIA_GPU, 0, s_type(1.0), false);
        done = Cs[t]->sync();
      }
      // Contract C
      unique_ptr<talsh::Tensor> ca(new talsh::Tensor({dims_grid[0][5][0],
                                                      dims_grid[0][6][1]},
                                                      s_type(0.0)));
      errc = ca->contractAccumulate(nullptr,
            "D(a,c)+=L(a,b)*R(b,c)",
            *Cs[0].get(),
            *Cs[1].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
      done = ca->sync();
      unique_ptr<talsh::Tensor> cb(new talsh::Tensor({dims_grid[1][4][0],
                                                      dims_grid[1][5][2],
                                                      dims_grid[1][5][3],
                                                      dims_grid[1][5][0]},
                                                      s_type(0.0)));
      errc = cb->contractAccumulate(nullptr,
            "D(a,d,e,c)+=L(a,b)*R(c,b,d,e)",
            *Cs[2].get(),
            *Cs[3].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
      done = cb->sync();
      unique_ptr<talsh::Tensor> cc(new talsh::Tensor({dims_grid[1][4][0],
                                                      dims_grid[1][5][2],
                                                      dims_grid[1][5][3],
                                                      dims_grid[0][6][1]},
                                                      s_type(0.0)));
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
                                                      s_type(0.0)));
      errc = cd->contractAccumulate(nullptr,
            "D(a,b,c,e)+=L(a,b,c,d)*R(d,e)",
            *Cs[4].get(),
            *Cs[5].get(), DEV_NVIDIA_GPU, 0, s_type(1.0), false);
      done = cd->sync();
      unique_ptr<talsh::Tensor> C(new talsh::Tensor({dims_grid[1][4][0],
                                                      dims_grid[1][5][2],
                                                      dims_grid[1][6][2],
                                                      dims_grid[1][7][1]},
                                                      s_type(0.0)));
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
            *eq.get(),
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
    t1 = high_resolution_clock::now();
    span = duration_cast<duration<double>>(t1 - t0);
    //cout << "Time contracting all Cs is " << span.count() << endl;
    eq.reset();
    
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
