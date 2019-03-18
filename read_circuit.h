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


/**
* @file read_circuit.h
* Helper functions to read quantum circuits from a file.
* @see https://github.com/benjaminvillalonga/optimized_parallel_QC_with_TN
*
* @author Benjamin Villalonga
* @date Created: September 2018
* @date Modified: September 2018
*/

#ifndef READ_CIRCUIT_
#define READ_CIRCUIT_

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <complex>
#include <ctime>
#include <chrono>
#include <cassert>

#include <iostream>
#include <fstream>
#include <sstream>

#include "mkl_tensor.h"

using namespace std;
using namespace chrono;

const double _SQRT_2 = 1.41421356237309504880168872;
const double _INV_SQRT_2 = 1./_SQRT_2;
const int SUPER_CYCLE_DEPTH = 8;
const int DIM = 2;

const unordered_map<string,vector<s_type>> _GATES_DATA({
  // Deltas.
  {"delta_0", vector<s_type>({1.0, 0.0})},
  {"delta_1", vector<s_type>({0.0, 1.0})},
  // For one-qubit gates, the first index is input an second is output.
  {"h", vector<s_type>({_INV_SQRT_2,_INV_SQRT_2,
                        _INV_SQRT_2,-_INV_SQRT_2})},
  {"t", vector<s_type>({1.0,0.,0.,{_INV_SQRT_2,_INV_SQRT_2}})},
  {"x_1_2", vector<s_type>({{0.5,0.5},
                            {0.5,-0.5},
                            {0.5,-0.5},
                            {0.5,0.5}})},
  {"y_1_2", vector<s_type>({{0.5,0.5},
                            {0.5,0.5},
                            {-0.5,-0.5},
                            {0.5,0.5}})},
  // For cz, both q1 and q2 get indices in the order (input, virtual, output).
  //{"cz_q1", vector<s_type>({1.,0.,0.,0.,0.,0.,0.,1.})},
  //{"cz_q2", vector<s_type>({1.,0.,1.,0.,0.,1.,0.,-1.})}});
  // Use more "balanced" SVD. Otherwise tensors are very sparse.
  {"cz_q1", vector<s_type>({-0.3446133714352872, 0.,
                            1.1381806476131544, 0.,
                            0., -1.1381806476131544,
                            0., -0.3446133714352872})},
  {"cz_q2", vector<s_type>({-1.0484937059720079, 0.,
                            0.5611368023131075, 0.,
                            0., 0.5611368023131075,
                            0., 1.0484937059720079})},
  // For the non-decomposed cz, the convention is (in1, in2, out1, out2).
  {"cz", vector<s_type>({1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,-1.})}
  });


/**
* Returns vector of vectors of s_type with the Schmidt decomposition of the
* fSim gate.
* @param theta double with the angle $theta$. Modulo $2\pi$ is taken.
* @param phi double with the angle $phi$. Modulo $2\pi$ is taken.
* @param scratch pointer to s_type array with scratch space for all operations
* performed in this function.
* @return vector<vector<s_type>> with three elements: first the vector with
* entries of the first qubit, second the vector with entries of the second
* qubit, and third the vector with the singular values (informational only).
*/
/*
vector<vector<s_type>> fSim(double theta, double phi, s_type * scratch)
{

  vector<s_type> coeffs({
                   {0.0,-0.5*sin(theta)},
                   {0.0,-0.5*sin(theta)},
                   {0.5*(cos(-theta/2.)-cos(theta)),0.5*sin(-theta/2.)},
                   {0.5*(cos(-theta/2.)+cos(theta)),0.5*sin(-theta/2.)}});

  vector<double> norm_coeffs(coeffs.size());
  for (int i=0; i<coeffs.size(); ++i)
  {
    norm_coeffs[i] = abs(coeffs[i]);
  }

  vector<vector<s_type>> q1_matrices({
    {{0.,0.},{1.,0.},{1.,0.},{0.,0.}},
    {{0.,0.},{0.,1.},{0.,-1.},{0.,0.}},
    {{cos(phi/4.),sin(phi/4.)},{0.,0.},{0.,0.},{-cos(-phi/4.),-sin(-phi/4.)}},
    {{cos(phi/4.),sin(phi/4.)},{0.,0.},{0.,0.},{cos(-phi/4.),sin(-phi/4.)}}});
  vector<vector<s_type>> q2_matrices({
    {{0.,0.},{1.,0.},{1.,0.},{0.,0.}},
    {{0.,0.},{0.,1.},{0.,-1.},{0.,0.}},
    {{cos(phi/4.),sin(phi/4.)},{0.,0.},{0.,0.},{-cos(-phi/4.),-sin(-phi/4.)}},
    {{cos(phi/4.),sin(phi/4.)},{0.,0.},{0.,0.},{cos(-phi/4.),sin(-phi/4.)}}});

  for (int i=0; i<coeffs.size(); ++i)
  {
    for (int j=0; j<q1_matrices[i].size(); ++j)
    {
      q1_matrices[i][j] *= coeffs[i];
      q2_matrices[i][j] *= coeffs[i];
    }
  }

  struct cnmm { 
    s_type c;
    double n;
    vector<s_type> m1;
    vector<s_type> m2;
    cnmm(s_type c_, double n_, vector<s_type> m1_, vector<s_type> m2_) :
         c(c_), n(n_), m1(m1_), m2(m2_) {}
    bool operator<( const cnmm & other ) const
    { return n < other.n; }
  };

  vector<cnmm> my_cnmm;
  for (int i=0; i<coeffs.size(); ++i)
  {
    my_cnmm.emplace_back(
          coeffs[i], norm_coeffs[i], q1_matrices[i], q2_matrices[i]);
  }

  sort(my_cnmm.begin(), my_cnmm.end());
  reverse(my_cnmm.begin(), my_cnmm.end());

  vector<s_type> q1_tensor, q2_tensor;
  for (auto v : my_cnmm)
  {
    for (auto w : v.m1)
      q1_tensor.emplace_back(w);
    for (auto w : v.m2)
      q2_tensor.emplace_back(w);
  }

  MKLTensor q1_mkltensor({"v","q1i","q2i"}, {4,2,2}, q1_tensor);
  MKLTensor q2_mkltensor({"v","q1i","q2i"}, {4,2,2}, q2_tensor);
  q1_mkltensor.reorder({"q1i","v","q2i"}, scratch);
  q2_mkltensor.reorder({"q1i","v","q2i"}, scratch);
  vector<s_type> q1_reordered_tensor(q1_mkltensor.size());
  vector<s_type> q2_reordered_tensor(q2_mkltensor.size());
  for (int i=0; i<q1_mkltensor.size(); ++i)
  {
    q1_reordered_tensor[i] = *(q1_mkltensor.data()+i);
    q2_reordered_tensor[i] = *(q2_mkltensor.data()+i);
  }

  sort(norm_coeffs.begin(), norm_coeffs.end());
  reverse(norm_coeffs.begin(), norm_coeffs.end());

  vector<vector<s_type>> ret_val({q1_reordered_tensor,
                                  q2_reordered_tensor,
                                  norm_coeffs});

  return ret_val;
}
*/


/**
* Returns spatial coordinates i and j on the grid given qubit number q.
* @param q int with the qubit number.
* @param J int with the second spatial dimension of the grid of qubits.
* @return int with the spatial coordinates i and j of qubit q on the grid.
*/
vector<int> _q_to_i_j(int q, int J)
{
  int i = q / J;
  return vector<int>({i, q-i*J});
}

/**
* Read circuit from file and fill in a 3D grid of tensors with the indices
* labelled as "(i1,j1,k1),(i2,j2,k2)", where i1<=i2, j1<=j2, and k1<=k2.
* I*J must be equal to the number of qubits; K must be equal to
* (depth_of_circuit-2)/8; initial_conf and final_conf must have the length
* equal to the number of qubits.
* @param filename string with the name of the circuit file.
* @param I int with the first spatial dimension of the grid of qubits.
* @param J int with the second spatial dimension of the grid of qubits.
* @param K int with depth of the grid of tensors (depth_of_circuit-2)/8.
* @param initial_conf string with 0s and 1s with the input configuration of
* the circuit.
* @param final_conf string with 0s and 1s with the output configuration of
* @param grid_of_tensors referenced to a vector<vector<vector<MKLTensor>>>
* with tensors at * each position of the grid.
* @param scratch pointer to s_type array with scratch space for all operations
* performed in this function.
* 
* the circuit.
*/
void google_circuit_file_to_grid_of_tensors(string filename, int I, int J,
        int K, string initial_conf, string final_conf,
        vector<vector<vector<MKLTensor>>> & grid_of_tensors, s_type * scratch)
{
  // Open file.
  auto io = ifstream(filename);
  assert(io.good() && "Cannot open file.");

  // Gotten from the file.
  int num_qubits, cycle, q1, q2;
  string gate;
  // Useful for plugging into the tensor network:
  vector<int> i_j_1, i_j_2;
  int super_cycle;

  // The first element should be the number of qubits
  io >> num_qubits;

  // Assert for the number of qubits.
  assert(num_qubits==I*J && "I*J must be equal to the number of qubits.");
  // Assert for the length of initial_conf and final_conf.
  {
    assert(initial_conf.size()==num_qubits
            && "initial_conf must be of size equal to the number of qubits.");
    assert(final_conf.size()==num_qubits
            && "final_conf must be of size equal to the number of qubits.");
  }

  // Creating grid variables.
  vector<vector<vector<vector<MKLTensor>>>> grid_of_groups_of_tensors(I);
  grid_of_tensors = vector<vector<vector<MKLTensor>>>(I);
  vector<vector<vector<int>>> counter_group(I);
  for (int i=0; i<I; ++i)
  {
    grid_of_groups_of_tensors[i] = vector<vector<vector<MKLTensor>>>(J);
    grid_of_tensors[i] = vector<vector<MKLTensor>>(J);
    counter_group[i] = vector<vector<int>>(J);
    for (int j=0; j<J; ++j)
    {
      grid_of_groups_of_tensors[i][j] = vector<vector<MKLTensor>>(K);
      grid_of_tensors[i][j] = vector<MKLTensor>(K);
      counter_group[i][j] = vector<int>(K, 0);
      for (int k=0; k<K; ++k)
      {
        grid_of_groups_of_tensors[i][j][k] = vector<MKLTensor>();
      }
    }
  }

  // Insert deltas and Hadamards to first layer.
  for (int q=0; q<num_qubits; ++q)
  {
    vector<int> i_j = _q_to_i_j(q, J);
    int i = i_j[0], j = i_j[1];
    string delta_gate = (initial_conf[q]=='0')?"delta_0":"delta_1";
    grid_of_groups_of_tensors[i][j][0].push_back(
                    MKLTensor({"th"}, {2}, _GATES_DATA.at(delta_gate)));
    grid_of_groups_of_tensors[i][j][0].push_back(
                    MKLTensor({"th","t0"}, {2,2}, _GATES_DATA.at("h")));
  }

  string line;
  while(getline(io, line)) if(line.size() && line[0] != '#') {// Avoid comments
    stringstream ss(line);
    // The first element is the cycle
    ss >> cycle;
    // The second element is the gate
    ss >> gate;
    // Get the first position
    ss >> q1;
    // Get the second position in the case
    if (gate=="cz") ss >> q2;
    else q2 = -1;

    // Get i, j and super_cycle
    i_j_1 = _q_to_i_j(q1, J);
    if (q2>=0)
    {
      i_j_2 = _q_to_i_j(q2, J);
    }
    super_cycle = (cycle - 1) / SUPER_CYCLE_DEPTH;

    // Fill in one-qubit gates.
    if (q2<0 && cycle>0 && cycle<=SUPER_CYCLE_DEPTH*K)
    {
      string input_index = "t"
              + to_string(counter_group[i_j_1[0]][i_j_1[1]][super_cycle]);
      string output_index = "t"
              + to_string(counter_group[i_j_1[0]][i_j_1[1]][super_cycle] + 1);
      ++counter_group[i_j_1[0]][i_j_1[1]][super_cycle];
      grid_of_groups_of_tensors[i_j_1[0]][i_j_1[1]][super_cycle].push_back(
        MKLTensor({input_index,output_index}, {2,2}, _GATES_DATA.at(gate)));
    }
    if (q2>=0 && cycle>0 && cycle<=SUPER_CYCLE_DEPTH*K)
    {
      string input_index_1 = "t"
              + to_string(counter_group[i_j_1[0]][i_j_1[1]][super_cycle]);
      string output_index_1 = "t"
              + to_string(counter_group[i_j_1[0]][i_j_1[1]][super_cycle] + 1);
      string virtual_index =
                      "("+to_string(i_j_1[0])+","+to_string(i_j_1[1])+","
                         +to_string(super_cycle)+"),("+to_string(i_j_2[0])+","
                         +to_string(i_j_2[1])+","+to_string(super_cycle)+")";
      string input_index_2 = "t"
              + to_string(counter_group[i_j_2[0]][i_j_2[1]][super_cycle]);
      string output_index_2 = "t"
              + to_string(counter_group[i_j_2[0]][i_j_2[1]][super_cycle] + 1);
      ++counter_group[i_j_1[0]][i_j_1[1]][super_cycle];
      ++counter_group[i_j_2[0]][i_j_2[1]][super_cycle];
      grid_of_groups_of_tensors[i_j_1[0]][i_j_1[1]][super_cycle].push_back(
                      MKLTensor({input_index_1,virtual_index,output_index_1},
                                {2,2,2}, _GATES_DATA.at("cz_q1")));
      grid_of_groups_of_tensors[i_j_2[0]][i_j_2[1]][super_cycle].push_back(
                      MKLTensor({input_index_2,virtual_index,output_index_2},
                                {2,2,2}, _GATES_DATA.at("cz_q2")));
    }
  }
  // Insert Hadamards and deltas to last layer.
  for (int q=0; q<num_qubits; ++q)
  {
    vector<int> i_j = _q_to_i_j(q, J);
    int i = i_j[0], j = i_j[1];
    int k = K-1;
    string delta_gate = (final_conf[q]=='0')?"delta_0":"delta_1";
    string last_index = "t"+to_string(counter_group[i][j][k]);
    grid_of_groups_of_tensors[i][j][k].push_back(
                    MKLTensor({"th",last_index}, {2,2}, _GATES_DATA.at("h")));
    grid_of_groups_of_tensors[i][j][k].push_back(
                    MKLTensor({"th"}, {2}, _GATES_DATA.at(delta_gate)));
  }

  // Contracting each group of gates into a single tensor.
  for (int i=0; i<I; ++i) for (int j=0; j<J; ++j) for (int k=0; k<K; ++k)
  {
    vector<MKLTensor> & group = grid_of_groups_of_tensors[i][j][k];
    vector<MKLTensor> group_containers(SUPER_CYCLE_DEPTH+2, // +2 for d and H.
                                       MKLTensor({""},{(int)pow(DIM,6)}));
    group_containers[0] = group[0];
    int t=1;
    for (t=1; t<group.size(); ++t)
    {
      multiply(group_containers[t-1], group[t], group_containers[t], scratch);
    }
    grid_of_tensors[i][j][k] = group_containers[t-1];
  }

  // Rename "t..." indices.
  for (int i=0; i<I; ++i) for (int j=0; j<J; ++j) for (int k=0; k<K; ++k)
  {
    if (k>0)
    {
      string new_first_index = "("+to_string(i)+","+to_string(j)+","
                               +to_string(k-1)+"),("+to_string(i)+","
                               +to_string(j)+","+to_string(k)+")";
      grid_of_tensors[i][j][k].rename_index("t0", new_first_index);
    }
    if (k<K-1)
    {
      string last_index = "t"+to_string(counter_group[i][j][k]);
      string new_last_index = "("+to_string(i)+","+to_string(j)+","
                              +to_string(k)+"),("+to_string(i)+","
                              +to_string(j)+","+to_string(k+1)+")";
      grid_of_tensors[i][j][k].rename_index(last_index, new_last_index);
    }
  }

  // Be proper about pointers.
  scratch = NULL;
}

/**
* Contracts a 3D grid of tensors onto a 2D grid of tensors, contracting
* in the time (thrid) direction, and renaming the indices accordingly.
* @param grid_of_tensors_3D reference to a
* vector<vector<vector<MKLTensor>>> with the 3D grid of tensors. It must be a
* grid dimensionswise. The typical names for the indices in a grid is assumed.
* @param grid_of_tensors_2D reference to a vector<vector<MKLTensor>> where the
* 2D grid of tensors will be stored. The typical names for the indices will
* be used.
* @param scratch pointer to s_type array with enough space for all scratch
* work.
*/
void grid_of_tensors_3D_to_2D(
      vector<vector<vector<MKLTensor>>> & grid_of_tensors_3D,
      vector<vector<MKLTensor>> & grid_of_tensors_2D, s_type * scratch)
{
  // Get dimensions and super_dim = DIM^k.
  const int I = grid_of_tensors_3D.size();
  const int J = grid_of_tensors_3D[0].size();
  const int K = grid_of_tensors_3D[0][0].size();
  const int super_dim = (int)pow(DIM,K);

  // Contract vertically and fill grid_of_tensors_2D.
  for (int i=0; i<I; ++i) for (int j=0; j<J; ++j)
  {
    vector<MKLTensor> group_containers(K-1,
                                    MKLTensor({""}, {(int)pow(super_dim,4)}));
    if (K==1)
    {
      group_containers[0] = grid_of_tensors_3D[i][j][0];
    }
    else
    {
      multiply(grid_of_tensors_3D[i][j][0], grid_of_tensors_3D[i][j][1],
               group_containers[0], scratch);
      for (int k=1; k<K-1; ++k)
      {
        multiply(group_containers[k-1], grid_of_tensors_3D[i][j][k+1],
                 group_containers[k], scratch);
      }
    }
    grid_of_tensors_2D[i][j] = group_containers[K-2];
  }

  // Reorder and bundle.
  for (int i=0; i<I; ++i) for (int j=0; j<J; ++j)
  {
    vector<string> ordered_indices_3D;
    vector<string> indices_2D;
    string index_name;
    if (i>0)
    {
      for (int k=0; k<K; ++k)
      {
        index_name = "("+to_string(i-1)+","+to_string(j)+","+to_string(k)+"),("+
                     to_string(i)+","+to_string(j)+","+to_string(k)+")";
        ordered_indices_3D.push_back(index_name);
      }
      index_name = "("+to_string(i-1)+","+to_string(j)+"),("
                   +to_string(i)+","+to_string(j)+")";
      indices_2D.push_back(index_name);
    }
    if (j>0)
    {
      for (int k=0; k<K; ++k)
      {
        index_name = "("+to_string(i)+","+to_string(j-1)+","+to_string(k)+"),("+
                     to_string(i)+","+to_string(j)+","+to_string(k)+")";
        ordered_indices_3D.push_back(index_name);
      }
      index_name = "("+to_string(i)+","+to_string(j-1)+"),("
                   +to_string(i)+","+to_string(j)+")";
      indices_2D.push_back(index_name);
    }
    if (i<I-1)
    {
      for (int k=0; k<K; ++k)
      {
        index_name = "("+to_string(i)+","+to_string(j)+","+to_string(k)+"),("+
                     to_string(i+1)+","+to_string(j)+","+to_string(k)+")";
        ordered_indices_3D.push_back(index_name);
      }
      index_name = "("+to_string(i)+","+to_string(j)+"),("
                   +to_string(i+1)+","+to_string(j)+")";
      indices_2D.push_back(index_name);
    }
    if (j<J-1)
    {
      for (int k=0; k<K; ++k)
      {
        index_name = "("+to_string(i)+","+to_string(j)+","+to_string(k)+"),("+
                     to_string(i)+","+to_string(j+1)+","+to_string(k)+")";
        ordered_indices_3D.push_back(index_name);
      }
      index_name = "("+to_string(i)+","+to_string(j)+"),("
                   +to_string(i)+","+to_string(j+1)+")";
      indices_2D.push_back(index_name);
    }
    // Reorder.
    grid_of_tensors_2D[i][j].reorder(ordered_indices_3D, scratch);

    // Bundle.
    for (int idx_num=0; idx_num<indices_2D.size(); ++idx_num)
    {
      vector<string> indices_to_bundle(ordered_indices_3D.begin()+idx_num*K,
                                     ordered_indices_3D.begin()+(idx_num+1)*K);
      grid_of_tensors_2D[i][j].bundle(indices_to_bundle, indices_2D[idx_num]);
    }
  }

  // Be proper about pointers.
  scratch = NULL;
}


/**
* Read circuit from file and fill in a 3D grid of tensors with the indices
* labelled as "(i1,j1,k1),(i2,j2,k2)", where i1<=i2, j1<=j2, and k1<=k2.
* I*J must be equal to the number of qubits; K must be equal to
* (depth_of_circuit-2)/8; initial_conf and final_conf must have the length
* equal to the number of qubits.
* @param filename string with the name of the circuit file.
* @param I int with the first spatial dimension of the grid of qubits.
* @param J int with the second spatial dimension of the grid of qubits.
* @param K int with depth of the grid of tensors (depth_of_circuit-2)/8.
* @param initial_conf string with 0s and 1s with the input configuration of
* the circuit.
* @param final_conf_B string with 0s and 1s with the output configuration on B.
* @param A vector<vector<int>> with the coords. of the qubits in A.
* @param off vector<vector<int>> with the coords. of the qubits turned off.
* @param grid_of_tensors referenced to a vector<vector<vector<MKLTensor>>>
* with tensors at * each position of the grid.
* @param scratch pointer to s_type array with scratch space for all operations
* performed in this function.
* 
* the circuit.
*/
void google_circuit_file_to_grid_of_tensors(string filename, int I, int J,
        int K, string initial_conf, string final_conf_B,
        vector<vector<int>> A, vector<vector<int>> off,
        vector<vector<vector<MKLTensor>>> & grid_of_tensors, s_type * scratch)
{
  // Open file.
  auto io = ifstream(filename);
  assert(io.good() && "Cannot open file.");

  // Gotten from the file.
  int num_qubits, cycle, q1, q2;
  string gate;
  // Useful for plugging into the tensor network:
  vector<int> i_j_1, i_j_2;
  int super_cycle;

  // The first element should be the number of qubits
  io >> num_qubits;
  num_qubits = I * J;

  // Assert for the number of qubits.
  assert(num_qubits==I*J && "I*J must be equal to the number of qubits.");
  // Assert for the length of initial_conf and final_conf.
  {
    assert(initial_conf.size()==num_qubits-off.size()
            && "initial_conf must be of size equal to the number of qubits.");
    assert(final_conf_B.size()==num_qubits-off.size()-A.size()
            && "final_conf_B must be of size equal to the number of qubits.");
  }

  // Creating grid variables.
  vector<vector<vector<vector<MKLTensor>>>> grid_of_groups_of_tensors(I);
  grid_of_tensors = vector<vector<vector<MKLTensor>>>(I);
  vector<vector<vector<int>>> counter_group(I);
  for (int i=0; i<I; ++i)
  {
    grid_of_groups_of_tensors[i] = vector<vector<vector<MKLTensor>>>(J);
    grid_of_tensors[i] = vector<vector<MKLTensor>>(J);
    counter_group[i] = vector<vector<int>>(J);
    for (int j=0; j<J; ++j)
    {
      grid_of_groups_of_tensors[i][j] = vector<vector<MKLTensor>>(K);
      grid_of_tensors[i][j] = vector<MKLTensor>(K);
      counter_group[i][j] = vector<int>(K, 0);
      for (int k=0; k<K; ++k)
      {
        grid_of_groups_of_tensors[i][j][k] = vector<MKLTensor>();
      }
    }
  }

  // Insert deltas and Hadamards to first layer.
  int idx = 0;
  for (int q=0; q<num_qubits; ++q)
  {
    vector<int> i_j = _q_to_i_j(q, J);
    int i = i_j[0], j = i_j[1];
    if (find(off.begin(),off.end(),
             vector<int>({i,j}))!=off.end())
    {
      continue;
    }
    string delta_gate = (initial_conf[idx]=='0')?"delta_0":"delta_1";
    grid_of_groups_of_tensors[i][j][0].push_back(
                    MKLTensor({"th"}, {2}, _GATES_DATA.at(delta_gate)));
    grid_of_groups_of_tensors[i][j][0].push_back(
                    MKLTensor({"th","t0"}, {2,2}, _GATES_DATA.at("h")));
    idx += 1;
  }

  string line;
  while(getline(io, line)) if(line.size() && line[0] != '#') {// Avoid comments
    stringstream ss(line);
    // The first element is the cycle
    ss >> cycle;
    // The second element is the gate
    ss >> gate;
    // Get the first position
    ss >> q1;
    // Get the second position in the case
    if (gate=="cz") ss >> q2;
    else q2 = -1;

    // Get i, j and super_cycle
    i_j_1 = _q_to_i_j(q1, J);
    if (q2>=0)
    {
      i_j_2 = _q_to_i_j(q2, J);
    }
    super_cycle = (cycle - 1) / SUPER_CYCLE_DEPTH;

    // Fill in one-qubit gates.
    if (q2<0 && cycle>0 && cycle<=SUPER_CYCLE_DEPTH*K)
    {
      if (find(off.begin(),off.end(),
               vector<int>({i_j_1[0],i_j_1[1]}))!=off.end())
      {
        continue;
      }
      string input_index = "t"
              + to_string(counter_group[i_j_1[0]][i_j_1[1]][super_cycle]);
      string output_index = "t"
              + to_string(counter_group[i_j_1[0]][i_j_1[1]][super_cycle] + 1);
      ++counter_group[i_j_1[0]][i_j_1[1]][super_cycle];
      grid_of_groups_of_tensors[i_j_1[0]][i_j_1[1]][super_cycle].push_back(
        MKLTensor({input_index,output_index}, {2,2}, _GATES_DATA.at(gate)));
    }
    if (q2>=0 && cycle>0 && cycle<=SUPER_CYCLE_DEPTH*K)
    {
      if (find(off.begin(),off.end(),
               vector<int>({i_j_1[0],i_j_1[1]}))!=off.end()
          ||
          find(off.begin(),off.end(),
               vector<int>({i_j_2[0],i_j_2[1]}))!=off.end())
      {
        continue;
      }
      string input_index_1 = "t"
              + to_string(counter_group[i_j_1[0]][i_j_1[1]][super_cycle]);
      string output_index_1 = "t"
              + to_string(counter_group[i_j_1[0]][i_j_1[1]][super_cycle] + 1);
      string virtual_index =
                      "("+to_string(i_j_1[0])+","+to_string(i_j_1[1])+","
                         +to_string(super_cycle)+"),("+to_string(i_j_2[0])+","
                         +to_string(i_j_2[1])+","+to_string(super_cycle)+")";
      string input_index_2 = "t"
              + to_string(counter_group[i_j_2[0]][i_j_2[1]][super_cycle]);
      string output_index_2 = "t"
              + to_string(counter_group[i_j_2[0]][i_j_2[1]][super_cycle] + 1);
      ++counter_group[i_j_1[0]][i_j_1[1]][super_cycle];
      ++counter_group[i_j_2[0]][i_j_2[1]][super_cycle];
      grid_of_groups_of_tensors[i_j_1[0]][i_j_1[1]][super_cycle].push_back(
                      MKLTensor({input_index_1,virtual_index,output_index_1},
                                {2,2,2}, _GATES_DATA.at("cz_q1")));
      grid_of_groups_of_tensors[i_j_2[0]][i_j_2[1]][super_cycle].push_back(
                      MKLTensor({input_index_2,virtual_index,output_index_2},
                                {2,2,2}, _GATES_DATA.at("cz_q2")));
    }
  }
  // Insert Hadamards and deltas to last layer.
  idx = 0;
  for (int q=0; q<num_qubits; ++q)
  {
    vector<int> i_j = _q_to_i_j(q, J);
    int i = i_j[0], j = i_j[1];
    int k = K-1;
    if (find(off.begin(),off.end(),
             vector<int>({i,j}))!=off.end())
    { continue; }
    string last_index = "t"+to_string(counter_group[i][j][k]);
    grid_of_groups_of_tensors[i][j][k].push_back(
                    MKLTensor({"th",last_index}, {2,2}, _GATES_DATA.at("h")));
    if (find(A.begin(),A.end(),vector<int>({i,j}))!=A.end())
    { continue; }
    string delta_gate = (final_conf_B[idx]=='0')?"delta_0":"delta_1";
    grid_of_groups_of_tensors[i][j][k].push_back(
                    MKLTensor({"th"}, {2}, _GATES_DATA.at(delta_gate)));
    idx += 1; // Move in B only.
  }

  // Contracting each group of gates into a single tensor.
  for (int i=0; i<I; ++i) for (int j=0; j<J; ++j) for (int k=0; k<K; ++k)
  {
    if (find(off.begin(),off.end(), vector<int>({i,j}))!=off.end())
    {
      continue;
    }

    vector<MKLTensor> & group = grid_of_groups_of_tensors[i][j][k];
    vector<MKLTensor> group_containers(SUPER_CYCLE_DEPTH+2, // +2 for d and H.
                                       MKLTensor({""},{(int)pow(DIM,6)}));
    group_containers[0] = group[0];
    int t=1;
    for (t=1; t<group.size(); ++t)
    {
      multiply(group_containers[t-1], group[t], group_containers[t], scratch);
    }
    grid_of_tensors[i][j][k] = group_containers[t-1];
  }

  // Rename "t..." indices.
  for (int i=0; i<I; ++i) for (int j=0; j<J; ++j) for (int k=0; k<K; ++k)
  {
    if (find(off.begin(),off.end(),vector<int>({i,j}))!=off.end())
    { continue; }
    if (k>0)
    {
      string new_first_index = "("+to_string(i)+","+to_string(j)+","
                               +to_string(k-1)+"),("+to_string(i)+","
                               +to_string(j)+","+to_string(k)+")";
      grid_of_tensors[i][j][k].rename_index("t0", new_first_index);
    }
    if (k<K-1)
    {
      string last_index = "t"+to_string(counter_group[i][j][k]);
      string new_last_index = "("+to_string(i)+","+to_string(j)+","
                              +to_string(k)+"),("+to_string(i)+","
                              +to_string(j)+","+to_string(k+1)+")";
      grid_of_tensors[i][j][k].rename_index(last_index, new_last_index);
    }
    if (k==K-1 && find(A.begin(),A.end(),vector<int>({i,j}))!=A.end())
    {
      string last_index = "th";
      string new_last_index = "("+to_string(i)+","+to_string(j)+"),(o)";
      grid_of_tensors[i][j][k].rename_index(last_index, new_last_index);
    }
  }

  // Be proper about pointers.
  scratch = NULL;
}


/**
* Contracts a 3D grid of tensors onto a 2D grid of tensors, contracting
* in the time (thrid) direction, and renaming the indices accordingly.
* @param grid_of_tensors_3D reference to a
* vector<vector<vector<MKLTensor>>> with the 3D grid of tensors. It must be a
* grid dimensionswise. The typical names for the indices in a grid is assumed.
* @param grid_of_tensors_2D reference to a vector<vector<MKLTensor>> where the
* 2D grid of tensors will be stored. The typical names for the indices will
* be used.
* @param A vector<vector<int>> with the coords. of the qubits in A.
* @param off vector<vector<int>> with the coords. of the qubits turned off.
* @param scratch pointer to s_type array with enough space for all scratch
* work.
*/
void grid_of_tensors_3D_to_2D(
      vector<vector<vector<MKLTensor>>> & grid_of_tensors_3D,
      vector<vector<MKLTensor>> & grid_of_tensors_2D,
      vector<vector<int>> A, vector<vector<int>> off, s_type * scratch)
{
  // Get dimensions and super_dim = DIM^k.
  const int I = grid_of_tensors_3D.size();
  const int J = grid_of_tensors_3D[0].size();
  const int K = grid_of_tensors_3D[0][0].size();
  const int super_dim = (int)pow(DIM,K);

  // Contract vertically and fill grid_of_tensors_2D.
  for (int i=0; i<I; ++i) for (int j=0; j<J; ++j)
  {
    if (find(off.begin(),off.end(),vector<int>({i,j}))!=off.end())
    { continue; }

    vector<MKLTensor> group_containers;
    if (find(A.begin(),A.end(),vector<int>({i,j}))!=A.end())
    {
      group_containers = vector<MKLTensor>(K-1,
                                MKLTensor({""}, {(int)pow(super_dim,4)*DIM}));
    }
    else
    {
      group_containers = vector<MKLTensor>(K-1,
                                MKLTensor({""}, {(int)pow(super_dim,4)}));
    }

    if (K==1)
    {
      group_containers[0] = grid_of_tensors_3D[i][j][0];
    }
    else
    {
      multiply(grid_of_tensors_3D[i][j][0], grid_of_tensors_3D[i][j][1],
               group_containers[0], scratch);
      for (int k=1; k<K-1; ++k)
      {
        multiply(group_containers[k-1], grid_of_tensors_3D[i][j][k+1],
                 group_containers[k], scratch);
      }
    }
    grid_of_tensors_2D[i][j] = group_containers[K-2];
  }

  // Reorder and bundle.
  for (int i=0; i<I; ++i) for (int j=0; j<J; ++j)
  {
    if (find(off.begin(),off.end(),vector<int>({i,j}))!=off.end())
    { continue; }

    vector<string> ordered_indices_3D;
    vector<string> indices_2D;
    string index_name;
    if (i>0 && find(off.begin(),off.end(),vector<int>({i-1,j}))==off.end())
    {
      for (int k=0; k<K; ++k)
      {
        index_name = "("+to_string(i-1)+","+to_string(j)+","+to_string(k)+"),("+
                     to_string(i)+","+to_string(j)+","+to_string(k)+")";
        ordered_indices_3D.push_back(index_name);
      }
      index_name = "("+to_string(i-1)+","+to_string(j)+"),("
                   +to_string(i)+","+to_string(j)+")";
      indices_2D.push_back(index_name);
    }
    if (j>0 && find(off.begin(),off.end(),vector<int>({i,j-1}))==off.end())
    {
      for (int k=0; k<K; ++k)
      {
        index_name = "("+to_string(i)+","+to_string(j-1)+","+to_string(k)+"),("+
                     to_string(i)+","+to_string(j)+","+to_string(k)+")";
        ordered_indices_3D.push_back(index_name);
      }
      index_name = "("+to_string(i)+","+to_string(j-1)+"),("
                   +to_string(i)+","+to_string(j)+")";
      indices_2D.push_back(index_name);
    }
    if (i<I-1 && find(off.begin(),off.end(),vector<int>({i+1,j}))==off.end())
    {
      for (int k=0; k<K; ++k)
      {
        index_name = "("+to_string(i)+","+to_string(j)+","+to_string(k)+"),("+
                     to_string(i+1)+","+to_string(j)+","+to_string(k)+")";
        ordered_indices_3D.push_back(index_name);
      }
      index_name = "("+to_string(i)+","+to_string(j)+"),("
                   +to_string(i+1)+","+to_string(j)+")";
      indices_2D.push_back(index_name);
    }
    if (j<J-1 && find(off.begin(),off.end(),vector<int>({i,j+1}))==off.end())
    {
      for (int k=0; k<K; ++k)
      {
        index_name = "("+to_string(i)+","+to_string(j)+","+to_string(k)+"),("+
                     to_string(i)+","+to_string(j+1)+","+to_string(k)+")";
        ordered_indices_3D.push_back(index_name);
      }
      index_name = "("+to_string(i)+","+to_string(j)+"),("
                   +to_string(i)+","+to_string(j+1)+")";
      indices_2D.push_back(index_name);
    }
    if (find(A.begin(),A.end(),vector<int>({i,j}))!=A.end())
    {
      index_name = "("+to_string(i)+","+to_string(j)+"),(o)";
      ordered_indices_3D.push_back(index_name);
    }

    // Reorder.
    grid_of_tensors_2D[i][j].reorder(ordered_indices_3D, scratch);

    // Bundle.
    int max_idx;
    max_idx = indices_2D.size();
    for (int idx_num=0; idx_num<max_idx; ++idx_num)
    {
      vector<string> indices_to_bundle(ordered_indices_3D.begin()+idx_num*K,
                                     ordered_indices_3D.begin()+(idx_num+1)*K);
      grid_of_tensors_2D[i][j].bundle(indices_to_bundle, indices_2D[idx_num]);
    }
  }

  // Be proper about pointers.
  scratch = NULL;
}


/**
* Read circuit from file and fill vector of tensors (of gates), vector with 
* names of the input indices of the tensors and vector with the output indices
* of the tensors.
* @param filename string with the name of the circuit file.
* @param I int with the number of qubits.
* @param gates reference to a vector<MKLTensor> to be filled with the gates.
* @param inputs reference to a vector<vector<string>> to be filled with the
* input indexes of the gates.
* @param outputs reference to a vector<vector<string>> to be filled with the
* output indexes of the gates.
* @param scratch pointer to s_type array with scratch space for operations
* performed in this function.
* 
*/
void read_wave_function_evolution(string filename, int I,
          vector<MKLTensor> & gates, vector<vector<string>> & inputs,
          vector<vector<string>> & outputs, s_type * scratch)
{
  // Open file.
  auto io = ifstream(filename);
  assert(io.good() && "Cannot open file.");

  // Gotten from the file.
  int num_qubits, cycle, q1, q2;
  string gate;
  // Useful for plugging into the tensor network:
  vector<int> i_j_1, i_j_2;
  int super_cycle;

  // The first element should be the number of qubits
  io >> num_qubits;

  // Assert for the number of qubits.
  assert(num_qubits==I && "I must be equal to the number of qubits.");

  string line;
  while(getline(io, line)) if(line.size() && line[0] != '#') {// Avoid comments
    stringstream ss(line);
    // The first element is the cycle
    ss >> cycle;
    // The second element is the gate
    ss >> gate;
    // Get the first position
    ss >> q1;
    // Get the second position in the case
    if (gate=="cz") ss >> q2;
    else q2 = -1;

    // Fill in one-qubit gates.
    if (q2<0)
    {
      string input_index = to_string(q1) + ",i";
      string output_index = to_string(q1) + ",o";
      gates.push_back(MKLTensor({input_index,output_index},{DIM,DIM},
                                _GATES_DATA.at(gate)));
      inputs.push_back({input_index});
      outputs.push_back({output_index});
    }
    if (q2>=0)
    {
      string input_index1 = to_string(q1) + ",i";
      string output_index1 = to_string(q1) + ",o";
      string input_index2 = to_string(q2) + ",i";
      string output_index2 = to_string(q2) + ",o";
      inputs.push_back({input_index1,input_index2});
      outputs.push_back({output_index1,output_index2});
      gates.push_back(MKLTensor(
                    {input_index1,input_index2,output_index1,output_index2},
                    {DIM,DIM,DIM,DIM},_GATES_DATA.at(gate)));
    }
  }

  // Be proper about pointers.
  scratch = NULL;
}

#endif
