/*

  Copyright © 2019, United States Government, as represented by the Administratora
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


#ifndef READ_CIRCUIT_
#define READ_CIRCUIT_

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <complex>
#include <ctime>
#include <chrono>
#include <cassert>
#include <unordered_map>

#include <iostream>
#include <fstream>
#include <sstream>

#include "talshxx.hpp"
#include "talsh_wrapper.h"

using namespace std;
using namespace chrono;

using s_type = complex<float>;

const double _SQRT_2 = 1.41421356237309504880168872;
const double _INV_SQRT_2 = 1./_SQRT_2;
const int DIM = 2;
const vector<string> _ALPHABET({"a","b","c","d","e","f","g","h","i","j",
                                "k","l","m","n","o","p","q","r","s","t",
                                "u","v","w","x","y","z",
                                "A","B","C","D","E","F","G","H","I","J",
                                "K","L","M","N","O","P","Q","R","S","T",
                                "U","V","W","X","Y","Z"});

// I've flipped the order of the indexes due to TALSH taking Fortran-like
// storage
const unordered_map<string,vector<s_type>> _GATES_DATA({
  // Deltas.
  {"I", vector<s_type>({1.0, 0.0, 0.0, 1.0})}, // Used to reorder (while mult.)
  {"delta_0", vector<s_type>({1.0, 0.0})},
  {"delta_1", vector<s_type>({0.0, 1.0})},
  // For one-qubit gates, the first index is output an second is input.
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
  // For cz, both q1 and q2 get indices in the order (output, virtual, input).
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
  // For the non-decomposed cz, the convention is (out2, out1, in2, in1).
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
* Return comma separated string with the string entries of the vector.
* @param pattern_vector vector<string> with string entries to concatenate.
* @return string with the concatenated entries of pattern_vector, separated by
* commas.
*/
string pattern_from_vector(vector<string> pattern_vector)
{
  string pattern_string("");
  for (auto v : pattern_vector)
    pattern_string += v + ",";
  pattern_string.erase(pattern_string.size()-1);
  return pattern_string;
}


/**
* Return volume of tensor from dimensions.
* @param dims vector<int> vector of dimensions.
* @return size_t with the volume given the set of dimensions.
*/
size_t volume_from_dims(vector<int> dims)
{
  size_t volume = 1;
  for (auto v : dims)
    volume *= v;
  return volume;
}


/**
* Return dims of a talsh::Tensor.
* @param T pointer to constant talsh::Tensor.
* @return vector<dim> with the dimensions of a talsh::Tensor.
*/
vector<int> dims_from_tensor(talsh::Tensor const * T)
{
  unsigned int num_dims;
  int const * dims_ptr = T->getDimExtents(num_dims);
  vector<int> dims;
  for (int i=0; i<num_dims; ++i)
    dims.push_back(dims_ptr[i]);
  dims_ptr = nullptr;
  return dims;
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
* @param initial_conf string with 0s and 1s with the input configuration of
* the circuit.
* @param final_conf_B string with 0s and 1s with the output configuration on B.
* @param A vector<vector<int>> with the coords. of the qubits in A.
* @param off vector<vector<int>> with the coords. of the qubits turned off.
* @param grid_of_tensors reference to a
* vector<vector<shared_ptr<talsh::Tensor>>>
* pointers a talsh::Tensor on each position of the grid.
* 
* the circuit.
*/
void google_circuit_file_to_grid_of_tensors(string filename, int I, int J,
        string initial_conf, string final_conf_B,
        vector<vector<int>> A, vector<vector<int>> off,
        vector<vector<shared_ptr<talsh::Tensor>>> & grid_of_tensors)
{

  int errc;

  // Open file.
  auto io = ifstream(filename);
  assert(io.good() && "Cannot open file.");

  // Gotten from the file.
  int num_qubits, cycle, q1, q2;
  string gate;
  // Useful for plugging into the tensor network:
  vector<int> i_j_1, i_j_2;

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
            && "final_conf_B must be of size equal to the number of qubits minus the number of qubits in A.");
  }

  // Create grid of vectors of tensors, of dims and of bonds
  vector<vector<vector<shared_ptr<talsh::Tensor>>>> grid_of_lists_of_tensors(I);
  vector<vector<vector<int>>> grid_of_lists_of_bonds(I);
  for (int i=0; i<I; ++i)
  {
    grid_of_lists_of_tensors[i] = vector<vector<shared_ptr<talsh::Tensor>>>(J);
    grid_of_lists_of_bonds[i] = vector<vector<int>>(J);
  }

  // Insert deltas to first layer.
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
    vector<s_type> gate_vector(_GATES_DATA.at(delta_gate));

    vector<int> dims({DIM});
    grid_of_lists_of_tensors[i][j].push_back(
          shared_ptr<talsh::Tensor>(new talsh::Tensor(dims, gate_vector)));
    grid_of_lists_of_bonds[i][j].push_back(-1);
    ++idx;
  }
  

  // Read gates from file
  string line;
  while(getline(io, line)) if(line.size() && line[0] != '#') // Avoid comments
  {
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

    // Get i, j
    i_j_1 = _q_to_i_j(q1, J);
    if (q2>=0)
    {
      i_j_2 = _q_to_i_j(q2, J);
    }

    // Push one-qubit gates
    if (q2<0)
    {
      if (find(off.begin(),off.end(),
               vector<int>({i_j_1[0],i_j_1[1]}))!=off.end())
      {
        continue;
      }
      vector<s_type> gate_vector(_GATES_DATA.at(gate));
      vector<int> dims({DIM,DIM});
      grid_of_lists_of_tensors[i_j_1[0]][i_j_1[1]].push_back(
        shared_ptr<talsh::Tensor>(new talsh::Tensor(dims, gate_vector)));
      grid_of_lists_of_bonds[i_j_1[0]][i_j_1[1]].push_back(-1);
    }

    // Push two-qubit gates
    if (q2>=0)
    {
      // If either qubit is in 'off'
      if (find(off.begin(),off.end(),
               vector<int>({i_j_1[0],i_j_1[1]}))!=off.end()
          ||
          find(off.begin(),off.end(),
               vector<int>({i_j_2[0],i_j_2[1]}))!=off.end())
      {
        continue;
      }
      // Gate 1
      vector<s_type> gate_vector_1(_GATES_DATA.at("cz_q1"));
      vector<int> dims({DIM,DIM,DIM});
      grid_of_lists_of_tensors[i_j_1[0]][i_j_1[1]].push_back(
        shared_ptr<talsh::Tensor>(new talsh::Tensor(dims, gate_vector_1)));
      vector<int> bond_type;
      if (i_j_1[0]<i_j_2[0]) bond_type = {2, 0};
      else if (i_j_1[1]<i_j_2[1]) bond_type = {3, 1};
      else if (i_j_1[0]>i_j_2[0]) bond_type = {0, 2};
      else if (i_j_1[1]>i_j_2[1]) bond_type = {1, 3};
      grid_of_lists_of_bonds[i_j_1[0]][i_j_1[1]].push_back(bond_type[0]);
      // Gate 2
      vector<s_type> gate_vector_2(_GATES_DATA.at("cz_q2"));
      grid_of_lists_of_tensors[i_j_2[0]][i_j_2[1]].push_back(
        shared_ptr<talsh::Tensor>(new talsh::Tensor(dims, gate_vector_2)));
      grid_of_lists_of_bonds[i_j_2[0]][i_j_2[1]].push_back(bond_type[1]);
    }

  }


  // Grid of tensors with out of order indexes, of out of order bond types,
  // of new to old bond positions, and of number of each bond type
  vector<vector<shared_ptr<talsh::Tensor>>> grid_of_ooo_tensors(I);
  vector<vector<vector<int>>> grid_of_ooo_bond_types(I);
  vector<vector<vector<int>>> grid_of_new_to_old_bond_positions(I);
  vector<vector<vector<int>>> grid_of_number_of_bond_types(I);
  for (int i=0; i<I; ++i)
  {
    grid_of_ooo_tensors[i] = vector<shared_ptr<talsh::Tensor>>(J);
    grid_of_ooo_bond_types[i] = vector<vector<int>>(J);
    grid_of_new_to_old_bond_positions[i] = vector<vector<int>>(J);
    grid_of_number_of_bond_types[i] = vector<vector<int>>(J, vector<int>(4,0));
  }


  // Contract each list of tensors onto a single one
  for (int q=0; q<num_qubits; ++q)
  {
    vector<int> i_j = _q_to_i_j(q, J);
    int i = i_j[0], j = i_j[1];
    if (find(off.begin(),off.end(), vector<int>({i,j}))!=off.end())
    { continue; }
    // Rank of the delta tensor is 1
    int current_rank(1), new_tensor_rank, new_rank;
    string pattern_D(""), pattern_L, pattern_R;
    vector<shared_ptr<talsh::Tensor>> list_of_tensors;
    // Assume there is at least 2 tensors per world-line
    for (int k=1; k<grid_of_lists_of_tensors[i][j].size(); ++k)
    {
      new_tensor_rank = grid_of_lists_of_tensors[i][j][k]->getRank();
      new_rank = current_rank + new_tensor_rank - 2;
      pattern_R = pattern_from_vector(vector<string>(_ALPHABET.cbegin(),
                                             _ALPHABET.cbegin()+current_rank));
      if (current_rank>1)
      {
        pattern_D = pattern_from_vector(vector<string>(_ALPHABET.cbegin()+1,
                                             _ALPHABET.cbegin()+current_rank));
      }
      else
        pattern_D = "";

      if (new_tensor_rank==1)
      {
        pattern_L = _ALPHABET[0];
      }
      else if (new_tensor_rank==2)
      {
        if (pattern_D=="")
          pattern_D = _ALPHABET[current_rank];
        else
          pattern_D = _ALPHABET[current_rank] + "," +  pattern_D;
        pattern_L = _ALPHABET[current_rank] + "," + _ALPHABET[0];
      }
      else if (new_tensor_rank==3)
      {
        if (pattern_D=="")
        {
          pattern_D = _ALPHABET[current_rank+1] + "," +
                      _ALPHABET[current_rank];
        }
        else
        {
          pattern_D = _ALPHABET[current_rank+1] + "," +
                      _ALPHABET[current_rank] + "," +  pattern_D;
        }
        pattern_L = _ALPHABET[current_rank+1] + "," +
                    _ALPHABET[current_rank] + "," + _ALPHABET[0];
      }
      else
        cout << "Some gate is not properly decomposed into tensors!" << endl;
      if (pattern_D[pattern_D.length()]==',')
        pattern_D.erase(pattern_D.length());
      string pattern = "D("+pattern_D+")+=L("+pattern_L+")*R("+pattern_R+")";

      shared_ptr<talsh::Tensor> L;
      shared_ptr<talsh::Tensor> R;
      shared_ptr<talsh::Tensor> D;
      // L is always the small one, and R the one that has been contracted
      // from the base.
      if (k==1)
      {
        L = grid_of_lists_of_tensors[i][j][1];
        R = grid_of_lists_of_tensors[i][j][0];
      } else {
        L = grid_of_lists_of_tensors[i][j][k];
        R = list_of_tensors[k-2];
      }
      vector<int> dims_L = dims_from_tensor(L.get());
      vector<int> dims_R = dims_from_tensor(R.get());
      vector<int> dims_D;
      for (int t=1; t<dims_L.size(); ++t)
        dims_D.push_back(dims_L[t]);
      for (int t=0; t<dims_R.size()-1; ++t)
        dims_D.push_back(dims_R[t]);
      if (k==grid_of_lists_of_tensors[i][j].size()-1)
      {
        grid_of_ooo_tensors[i][j] =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_D, s_type(0.0)));
        D = grid_of_ooo_tensors[i][j];
      } else {
        list_of_tensors.push_back(
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_D, s_type(0.0))));
        D = shared_ptr<talsh::Tensor>(list_of_tensors[k-1]);
      }
      TensContraction contraction(pattern, D.get(), L.get(), R.get());
      errc = contraction.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(contraction.sync(DEV_HOST,0));
      current_rank = new_rank;
    }
  }


  // Reorder data by multiplying times the identity
  // Create grid of new to old bond positions and grid of number of bonds of
  // each type
  for (int q=0; q<num_qubits; ++q)
  {
    vector<int> i_j = _q_to_i_j(q, J);
    int i = i_j[0], j = i_j[1];
    if (find(off.begin(),off.end(), vector<int>({i,j}))!=off.end())
    { continue; }
    for (auto v : grid_of_lists_of_bonds[i][j])
    { if (v>=0)
        grid_of_ooo_bond_types[i][j].push_back(v);
    }
    // Loop over bond types
    for (int b=0; b<4; ++b)
    {
      for (int p=0; p<grid_of_ooo_bond_types[i][j].size(); ++p)
      {
        int v = grid_of_ooo_bond_types[i][j][p];
        if (v==b)
        {
          grid_of_new_to_old_bond_positions[i][j].push_back(p);
          grid_of_number_of_bond_types[i][j][b] += 1;
        }
      }
    }
  }

  // Because of multiplhying every new gate as the L tensor, and having the
  // built tensor as R, the indexes are backwards. First, flip the entries
  // of the map, so that they take into account the fact that indexes entering
  // first were put last in the list (the 'old' part). Second, make sure to
  // have the right destination, so that indexes that I want to put first are
  // actually first.
  for (int q=0; q<num_qubits; ++q)
  {
    vector<int> i_j = _q_to_i_j(q, J);
    int i = i_j[0], j = i_j[1];
    if (find(off.begin(),off.end(), vector<int>({i,j}))!=off.end())
    { continue; }
    vector<int> helper;
    for (auto v : grid_of_new_to_old_bond_positions[i][j])
      helper.push_back(v);
    int size = helper.size();
    for (int p=0; p<size; ++p)
    {
      // 'new' = 'old'
      grid_of_new_to_old_bond_positions[i][j][p] = size-helper[p]-1;
    }
  }

  // Reorder tensors and bundle indexes
  idx = 0;
  for (int q=0; q<num_qubits; ++q)
  {
    vector<int> i_j = _q_to_i_j(q, J);
    int i = i_j[0], j = i_j[1];
    if (find(off.begin(),off.end(), vector<int>({i,j}))!=off.end())
    { continue; }

    // Reorder everything onto tensor T
    vector<int> ooo_dims = dims_from_tensor(grid_of_ooo_tensors[i][j].get());
    vector<int> dims_T(ooo_dims);
    // The real index is in the first position
    for (int p=0; p<dims_T.size()-1; ++p)
    {
      dims_T[p+1] = ooo_dims[grid_of_new_to_old_bond_positions[i][j][p]+1];
    }
    talsh::Tensor T(dims_T, s_type(0.0));
    vector<s_type> vector_I(_GATES_DATA.at("I"));
    vector<int> dims_I({DIM,DIM});
    talsh::Tensor I(dims_I, vector_I);
    vector<string> pattern_vector_D({_ALPHABET[0]});
    for (auto v : grid_of_new_to_old_bond_positions[i][j])
      pattern_vector_D.push_back(_ALPHABET[v+2]);
    string pattern_R(pattern_from_vector(vector<string>(_ALPHABET.cbegin()+2,
                                         _ALPHABET.cbegin()+dims_T.size()+1)));
    pattern_R = _ALPHABET[1] + "," + pattern_R;
    string pattern_D(pattern_from_vector(pattern_vector_D));
    string pattern_L(pattern_from_vector(vector<string>(_ALPHABET.cbegin(),
                                         _ALPHABET.cbegin()+2)));
    string pattern = "D("+pattern_D+")+=L("+pattern_L+")*R("+pattern_R+")";
    TensContraction contraction(pattern, &T, &I,
                                         grid_of_ooo_tensors[i][j].get());
    errc = contraction.execute(DEV_NVIDIA_GPU,0);
    assert(errc==TALSH_SUCCESS);
    assert(contraction.sync(DEV_HOST,0));

    // Bundle indexes: copy T to B, with the right bundled dimensions.
    // If the qubit is not contracted by a delta, don't do B, but the tensor
    // at the right position on the final grid.
    // If a qubit is contracted by delta, contract onto the 
    // the tensor at the right position of the final grid.
    vector<int> dim_per_super_bond;
    int super_idx = 1;
    for (int b=0; b<4; ++b)
    {
      int number_of_bonds = grid_of_number_of_bond_types[i][j][b];
      if (number_of_bonds == 0) { continue; }
      vector<int> dims(dims_T.cbegin()+super_idx,
                       dims_T.cbegin()+super_idx+number_of_bonds);
      dim_per_super_bond.push_back(volume_from_dims(dims));
      super_idx += number_of_bonds;
    }
    vector<int> super_dims({DIM});
    for (int p=0; p<dim_per_super_bond.size(); ++p)
      super_dims.push_back(dim_per_super_bond[p]);
    shared_ptr<talsh::Tensor> B(new talsh::Tensor(super_dims, s_type(0.0)));
    s_type const * data_T;
    s_type * data_B;
    T.getDataAccessHostConst(&data_T);
    B->getDataAccessHost(&data_B);
    for (size_t p=0; p<T.getVolume(); ++p)
      data_B[p] = data_T[p];
    data_T = nullptr; data_B = nullptr;
    if (find(A.begin(),A.end(), vector<int>({i,j}))!=A.end()) // No delta
    {
      grid_of_tensors[i][j] = B;
    } else { // Contract with delta
      string delta_gate = (final_conf_B[idx]=='0')?"delta_0":"delta_1";
      vector<s_type> gate_vector(_GATES_DATA.at(delta_gate));
      vector<int> dims_delta({DIM});
      talsh::Tensor delta(dims_delta, gate_vector);

      vector<int> dims_C(super_dims.cbegin()+1, super_dims.cend());
      shared_ptr<talsh::Tensor> C(new talsh::Tensor(dims_C, s_type(0.0)));

      pattern_D = pattern_from_vector(vector<string>(_ALPHABET.cbegin()+1,
                                       _ALPHABET.cbegin()+super_dims.size()));
      pattern_L = _ALPHABET[0];
      pattern_R = pattern_from_vector(vector<string>(_ALPHABET.cbegin(),
                                       _ALPHABET.cbegin()+super_dims.size()));
      string pattern_delta = "D("+pattern_D+")+=L("+pattern_L+
                             ")*R("+pattern_R+")";
      TensContraction contraction_delta(pattern_delta, C.get(),
                                                       &delta, B.get());
      errc = contraction_delta.execute(DEV_NVIDIA_GPU,0);
      assert(errc==TALSH_SUCCESS);
      assert(contraction_delta.sync(DEV_HOST,0));

      grid_of_tensors[i][j] = C;
      ++idx;
    }
  }

  // Close io file
  io.close();
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
* @param initial_conf string with 0s and 1s with the input configuration of
* the circuit.
* @param off vector<vector<int>> with the coords. of the qubits turned off.
* @param grid_of_tensors reference to a
* vector<vector<shared_ptr<talsh::Tensor>>>
* pointers a talsh::Tensor on each position of the grid.
* 
* the circuit.
*/
void google_circuit_file_to_open_grid_of_tensors(string filename, int I, int J,
        string initial_conf, vector<vector<int>> off,
        vector<vector<shared_ptr<talsh::Tensor>>> & grid_of_tensors)
{

  int errc;

  // Open file.
  auto io = ifstream(filename);
  assert(io.good() && "Cannot open file.");

  // Gotten from the file.
  int num_qubits, cycle, q1, q2;
  string gate;
  // Useful for plugging into the tensor network:
  vector<int> i_j_1, i_j_2;

  // The first element should be the number of qubits
  io >> num_qubits;
  num_qubits = I * J;

  // Assert for the number of qubits.
  assert(num_qubits==I*J && "I*J must be equal to the number of qubits.");
  // Assert for the length of initial_conf and final_conf.
  {
    assert(initial_conf.size()==num_qubits-off.size()
            && "initial_conf must be of size equal to the number of qubits.");
  }

  // Create grid of vectors of tensors, of dims and of bonds
  vector<vector<vector<shared_ptr<talsh::Tensor>>>> grid_of_lists_of_tensors(I);
  vector<vector<vector<int>>> grid_of_lists_of_bonds(I);
  for (int i=0; i<I; ++i)
  {
    grid_of_lists_of_tensors[i] = vector<vector<shared_ptr<talsh::Tensor>>>(J);
    grid_of_lists_of_bonds[i] = vector<vector<int>>(J);
  }

  // Insert deltas to first layer.
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
    vector<s_type> gate_vector(_GATES_DATA.at(delta_gate));

    vector<int> dims({DIM});
    grid_of_lists_of_tensors[i][j].push_back(
          shared_ptr<talsh::Tensor>(new talsh::Tensor(dims, gate_vector)));
    grid_of_lists_of_bonds[i][j].push_back(-1);
    ++idx;
  }
  

  // Read gates from file
  string line;
  while(getline(io, line)) if(line.size() && line[0] != '#') // Avoid comments
  {
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

    // Get i, j
    i_j_1 = _q_to_i_j(q1, J);
    if (q2>=0)
    {
      i_j_2 = _q_to_i_j(q2, J);
    }

    // Push one-qubit gates
    if (q2<0)
    {
      if (find(off.begin(),off.end(),
               vector<int>({i_j_1[0],i_j_1[1]}))!=off.end())
      {
        continue;
      }
      vector<s_type> gate_vector(_GATES_DATA.at(gate));
      vector<int> dims({DIM,DIM});
      grid_of_lists_of_tensors[i_j_1[0]][i_j_1[1]].push_back(
        shared_ptr<talsh::Tensor>(new talsh::Tensor(dims, gate_vector)));
      grid_of_lists_of_bonds[i_j_1[0]][i_j_1[1]].push_back(-1);
    }

    // Push two-qubit gates
    if (q2>=0)
    {
      // If either qubit is in 'off'
      if (find(off.begin(),off.end(),
               vector<int>({i_j_1[0],i_j_1[1]}))!=off.end()
          ||
          find(off.begin(),off.end(),
               vector<int>({i_j_2[0],i_j_2[1]}))!=off.end())
      {
        continue;
      }
      // Gate 1
      vector<s_type> gate_vector_1(_GATES_DATA.at("cz_q1"));
      vector<int> dims({DIM,DIM,DIM});
      grid_of_lists_of_tensors[i_j_1[0]][i_j_1[1]].push_back(
        shared_ptr<talsh::Tensor>(new talsh::Tensor(dims, gate_vector_1)));
      vector<int> bond_type;
      if (i_j_1[0]<i_j_2[0]) bond_type = {2, 0};
      else if (i_j_1[1]<i_j_2[1]) bond_type = {3, 1};
      else if (i_j_1[0]>i_j_2[0]) bond_type = {0, 2};
      else if (i_j_1[1]>i_j_2[1]) bond_type = {1, 3};
      grid_of_lists_of_bonds[i_j_1[0]][i_j_1[1]].push_back(bond_type[0]);
      // Gate 2
      vector<s_type> gate_vector_2(_GATES_DATA.at("cz_q2"));
      grid_of_lists_of_tensors[i_j_2[0]][i_j_2[1]].push_back(
        shared_ptr<talsh::Tensor>(new talsh::Tensor(dims, gate_vector_2)));
      grid_of_lists_of_bonds[i_j_2[0]][i_j_2[1]].push_back(bond_type[1]);
    }

  }


  // Grid of tensors with out of order indexes, of out of order bond types,
  // of new to old bond positions, and of number of each bond type
  vector<vector<shared_ptr<talsh::Tensor>>> grid_of_ooo_tensors(I);
  vector<vector<vector<int>>> grid_of_ooo_bond_types(I);
  vector<vector<vector<int>>> grid_of_new_to_old_bond_positions(I);
  vector<vector<vector<int>>> grid_of_number_of_bond_types(I);
  for (int i=0; i<I; ++i)
  {
    grid_of_ooo_tensors[i] = vector<shared_ptr<talsh::Tensor>>(J);
    grid_of_ooo_bond_types[i] = vector<vector<int>>(J);
    grid_of_new_to_old_bond_positions[i] = vector<vector<int>>(J);
    grid_of_number_of_bond_types[i] = vector<vector<int>>(J, vector<int>(4,0));
  }


  // Contract each list of tensors onto a single one
  for (int q=0; q<num_qubits; ++q)
  {
    vector<int> i_j = _q_to_i_j(q, J);
    int i = i_j[0], j = i_j[1];
    if (find(off.begin(),off.end(), vector<int>({i,j}))!=off.end())
    { continue; }
    // Rank of the delta tensor is 1
    int current_rank(1), new_tensor_rank, new_rank;
    string pattern_D(""), pattern_L, pattern_R;
    vector<shared_ptr<talsh::Tensor>> list_of_tensors;
    // Assume there is at least 2 tensors per world-line
    for (int k=1; k<grid_of_lists_of_tensors[i][j].size(); ++k)
    {
      new_tensor_rank = grid_of_lists_of_tensors[i][j][k]->getRank();
      new_rank = current_rank + new_tensor_rank - 2;
      pattern_R = pattern_from_vector(vector<string>(_ALPHABET.cbegin(),
                                             _ALPHABET.cbegin()+current_rank));
      if (current_rank>1)
      {
        pattern_D = pattern_from_vector(vector<string>(_ALPHABET.cbegin()+1,
                                             _ALPHABET.cbegin()+current_rank));
      }
      else
        pattern_D = "";

      if (new_tensor_rank==1)
      {
        pattern_L = _ALPHABET[0];
      }
      else if (new_tensor_rank==2)
      {
        if (pattern_D=="")
          pattern_D = _ALPHABET[current_rank];
        else
          pattern_D = _ALPHABET[current_rank] + "," +  pattern_D;
        pattern_L = _ALPHABET[current_rank] + "," + _ALPHABET[0];
      }
      else if (new_tensor_rank==3)
      {
        if (pattern_D=="")
        {
          pattern_D = _ALPHABET[current_rank+1] + "," +
                      _ALPHABET[current_rank];
        }
        else
        {
          pattern_D = _ALPHABET[current_rank+1] + "," +
                      _ALPHABET[current_rank] + "," +  pattern_D;
        }
        pattern_L = _ALPHABET[current_rank+1] + "," +
                    _ALPHABET[current_rank] + "," + _ALPHABET[0];
      }
      else
        cout << "Some gate is not properly decomposed into tensors!" << endl;
      if (pattern_D[pattern_D.length()]==',')
        pattern_D.erase(pattern_D.length());
      string pattern = "D("+pattern_D+")+=L("+pattern_L+")*R("+pattern_R+")";

      shared_ptr<talsh::Tensor> L;
      shared_ptr<talsh::Tensor> R;
      shared_ptr<talsh::Tensor> D;
      // L is always the small one, and R the one that has been contracted
      // from the base.
      if (k==1)
      {
        L = grid_of_lists_of_tensors[i][j][1];
        R = grid_of_lists_of_tensors[i][j][0];
      } else {
        L = grid_of_lists_of_tensors[i][j][k];
        R = list_of_tensors[k-2];
      }
      vector<int> dims_L = dims_from_tensor(L.get());
      vector<int> dims_R = dims_from_tensor(R.get());
      vector<int> dims_D;
      for (int t=1; t<dims_L.size(); ++t)
        dims_D.push_back(dims_L[t]);
      for (int t=0; t<dims_R.size()-1; ++t)
        dims_D.push_back(dims_R[t]);
      if (k==grid_of_lists_of_tensors[i][j].size()-1)
      {
        grid_of_ooo_tensors[i][j] =
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_D, s_type(0.0)));
        D = grid_of_ooo_tensors[i][j];
      } else {
        list_of_tensors.push_back(
            shared_ptr<talsh::Tensor>(new talsh::Tensor(dims_D, s_type(0.0))));
        D = shared_ptr<talsh::Tensor>(list_of_tensors[k-1]);
      }
      TensContraction contraction(pattern, D.get(), L.get(), R.get());
      errc = contraction.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
      assert(contraction.sync(DEV_HOST,0));
      current_rank = new_rank;
    }
  }


  // Reorder data by multiplying times the identity
  // Create grid of new to old bond positions and grid of number of bonds of
  // each type
  for (int q=0; q<num_qubits; ++q)
  {
    vector<int> i_j = _q_to_i_j(q, J);
    int i = i_j[0], j = i_j[1];
    if (find(off.begin(),off.end(), vector<int>({i,j}))!=off.end())
    { continue; }
    for (auto v : grid_of_lists_of_bonds[i][j])
    { if (v>=0)
        grid_of_ooo_bond_types[i][j].push_back(v);
    }
    // Loop over bond types
    for (int b=0; b<4; ++b)
    {
      for (int p=0; p<grid_of_ooo_bond_types[i][j].size(); ++p)
      {
        int v = grid_of_ooo_bond_types[i][j][p];
        if (v==b)
        {
          grid_of_new_to_old_bond_positions[i][j].push_back(p);
          grid_of_number_of_bond_types[i][j][b] += 1;
        }
      }
    }
  }

  // Because of multiplhying every new gate as the L tensor, and having the
  // built tensor as R, the indexes are backwards. First, flip the entries
  // of the map, so that they take into account the fact that indexes entering
  // first were put last in the list (the 'old' part). Second, make sure to
  // have the right destination, so that indexes that I want to put first are
  // actually first.
  for (int q=0; q<num_qubits; ++q)
  {
    vector<int> i_j = _q_to_i_j(q, J);
    int i = i_j[0], j = i_j[1];
    if (find(off.begin(),off.end(), vector<int>({i,j}))!=off.end())
    { continue; }
    vector<int> helper;
    for (auto v : grid_of_new_to_old_bond_positions[i][j])
      helper.push_back(v);
    int size = helper.size();
    for (int p=0; p<size; ++p)
    {
      // 'new' = 'old'
      grid_of_new_to_old_bond_positions[i][j][p] = size-helper[p]-1;
    }
  }

  // Reorder tensors and bundle indexes
  idx = 0;
  for (int q=0; q<num_qubits; ++q)
  {
    vector<int> i_j = _q_to_i_j(q, J);
    int i = i_j[0], j = i_j[1];
    if (find(off.begin(),off.end(), vector<int>({i,j}))!=off.end())
    { continue; }

    // Reorder everything onto tensor T
    vector<int> ooo_dims = dims_from_tensor(grid_of_ooo_tensors[i][j].get());
    vector<int> dims_T(ooo_dims);
    // The real index is in the first position
    for (int p=0; p<dims_T.size()-1; ++p)
    {
      dims_T[p+1] = ooo_dims[grid_of_new_to_old_bond_positions[i][j][p]+1];
    }
    talsh::Tensor T(dims_T, s_type(0.0));
    vector<s_type> vector_I(_GATES_DATA.at("I"));
    vector<int> dims_I({DIM,DIM});
    talsh::Tensor I(dims_I, vector_I);
    vector<string> pattern_vector_D({_ALPHABET[0]});
    for (auto v : grid_of_new_to_old_bond_positions[i][j])
      pattern_vector_D.push_back(_ALPHABET[v+2]);
    string pattern_R(pattern_from_vector(vector<string>(_ALPHABET.cbegin()+2,
                                         _ALPHABET.cbegin()+dims_T.size()+1)));
    pattern_R = _ALPHABET[1] + "," + pattern_R;
    string pattern_D(pattern_from_vector(pattern_vector_D));
    string pattern_L(pattern_from_vector(vector<string>(_ALPHABET.cbegin(),
                                         _ALPHABET.cbegin()+2)));
    string pattern = "D("+pattern_D+")+=L("+pattern_L+")*R("+pattern_R+")";
    TensContraction contraction(pattern, &T, &I,
                                         grid_of_ooo_tensors[i][j].get());
    errc = contraction.execute(DEV_NVIDIA_GPU,0);
    assert(errc==TALSH_SUCCESS);
    assert(contraction.sync(DEV_HOST,0));

    // Bundle indexes: copy T to B, with the right bundled dimensions.
    // If the qubit is not contracted by a delta, don't do B, but the tensor
    // at the right position on the final grid.
    // If a qubit is contracted by delta, contract onto the 
    // the tensor at the right position of the final grid.
    vector<int> dim_per_super_bond;
    int super_idx = 1;
    for (int b=0; b<4; ++b)
    {
      int number_of_bonds = grid_of_number_of_bond_types[i][j][b];
      if (number_of_bonds == 0) { continue; }
      vector<int> dims(dims_T.cbegin()+super_idx,
                       dims_T.cbegin()+super_idx+number_of_bonds);
      dim_per_super_bond.push_back(volume_from_dims(dims));
      super_idx += number_of_bonds;
    }
    vector<int> super_dims({DIM});
    for (int p=0; p<dim_per_super_bond.size(); ++p)
      super_dims.push_back(dim_per_super_bond[p]);
    shared_ptr<talsh::Tensor> B(new talsh::Tensor(super_dims, s_type(0.0)));
    s_type const * data_T;
    s_type * data_B;
    T.getDataAccessHostConst(&data_T);
    B->getDataAccessHost(&data_B);
    for (size_t p=0; p<T.getVolume(); ++p)
      data_B[p] = data_T[p];
    data_T = nullptr; data_B = nullptr;

    // grid_of_tensors points to the right tensor
    grid_of_tensors[i][j] = B;
  }

  // Close io file
  io.close();
}


/**
* Close circuit on B region.
* @param filename string with the name of the circuit file.
* @param I int with the first spatial dimension of the grid of qubits.
* @param J int with the second spatial dimension of the grid of qubits.
* @param final_conf_B string with 0s and 1s with the output configuration on B.
* @param A vector<vector<int>> with the coords. of the qubits in A.
* @param off vector<vector<int>> with the coords. of the qubits turned off.
* @param grid_of_tensors reference to a
* vector<vector<shared_ptr<talsh::Tensor>>>
* pointers a talsh::Tensor on each position of the grid.
* 
* the circuit.
*/
void close_circuit(int I, int J, string final_conf_B, vector<vector<int>> A,
            vector<vector<int>> off,
            vector<vector<shared_ptr<talsh::Tensor>>> & open_grid_of_tensors,
            vector<vector<shared_ptr<talsh::Tensor>>> & grid_of_tensors)
{
  
  int errc;

  // Reorder tensors and bundle indexes
  int num_qubits = I*J;
  int idx = 0;
  string pattern_D(""), pattern_L(""), pattern_R("");
  for (int q=0; q<num_qubits; ++q)
  {
    vector<int> i_j = _q_to_i_j(q, J);
    int i = i_j[0], j = i_j[1];
    if (find(off.begin(),off.end(), vector<int>({i,j}))!=off.end())
    { continue; }

    if (find(A.begin(),A.end(), vector<int>({i,j}))!=A.end()) // No delta
    {
      grid_of_tensors[i][j] = open_grid_of_tensors[i][j];
    } else { // Contract with delta
      vector<int> dims_open_T(dims_from_tensor(
                                          open_grid_of_tensors[i][j].get()));
      vector<int> dims_T(dims_open_T.cbegin()+1,dims_open_T.cend());
      shared_ptr<talsh::Tensor> T(new talsh::Tensor(dims_T, s_type(0.0)));

      string delta_gate = (final_conf_B[idx]=='0')?"delta_0":"delta_1";
      vector<s_type> gate_vector(_GATES_DATA.at(delta_gate));
      vector<int> dims_delta({DIM});
      talsh::Tensor delta(dims_delta, gate_vector);

      pattern_D = pattern_from_vector(vector<string>(_ALPHABET.cbegin()+1,
                                     _ALPHABET.cbegin()+dims_open_T.size()));
      pattern_L = _ALPHABET[0];
      pattern_R = pattern_from_vector(vector<string>(_ALPHABET.cbegin(),
                                     _ALPHABET.cbegin()+dims_open_T.size()));
      string pattern_delta = "D("+pattern_D+")+=L("+pattern_L+
                             ")*R("+pattern_R+")";
      TensContraction contraction_delta(pattern_delta, T.get(),
                                   &delta, open_grid_of_tensors[i][j].get());
      errc = contraction_delta.execute(DEV_NVIDIA_GPU,0);
      assert(errc==TALSH_SUCCESS);
      assert(contraction_delta.sync(DEV_HOST,0));

      grid_of_tensors[i][j] = T;
      ++idx;
    }
  }
 
}


#endif
