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
* @param grid_of_tensors reference to a vector<vector<talsh::Tensor *>>
* pointers a talsh::Tensor on each position of the grid.
* 
* the circuit.
*/
void google_circuit_file_to_grid_of_tensors(string filename, int I, int J,
        string initial_conf, string final_conf_B,
        vector<vector<int>> A, vector<vector<int>> off,
        vector<vector<talsh::Tensor *>> & grid_of_tensors)
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

  // Create grid of vectors of tensors, of dims and of bonds
  vector<vector<vector<talsh::Tensor>>> grid_of_lists_of_tensors(I);
  vector<vector<vector<vector<int>>>> grid_of_lists_of_dims(I);
  vector<vector<vector<int>>> grid_of_lists_of_bonds(I);
  for (int i=0; i<I; ++i)
  {
    grid_of_lists_of_tensors[i] = vector<vector<talsh::Tensor>>(J);
    grid_of_lists_of_dims[i] = vector<vector<vector<int>>>(J);
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
    s_type * gate_data = (s_type *) malloc(gate_vector.size()*sizeof(s_type));
    for (size_t p=0; p<gate_vector.size(); ++p)
      *(gate_data+p) = gate_vector[p];
    vector<int> dims({DIM});
    grid_of_lists_of_tensors[i][j].push_back(talsh::Tensor(dims, gate_data));
    gate_data = NULL;
    grid_of_lists_of_dims[i][j].push_back(dims);
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
      s_type * gate_data = (s_type *) malloc(gate_vector.size()*sizeof(s_type));
      for (size_t p=0; p<gate_vector.size(); ++p)
        *(gate_data+p) = gate_vector[p];
      vector<int> dims({DIM,DIM});
      grid_of_lists_of_tensors[i_j_1[0]][i_j_1[1]].push_back(
        talsh::Tensor(dims, gate_data));
      gate_data = NULL;
      grid_of_lists_of_dims[i_j_1[0]][i_j_1[1]].push_back(dims);
      grid_of_lists_of_bonds[i_j_1[0]][i_j_1[1]].push_back(-1);
    }

    // Push one-qubit gates
    if (q2>=0)
    {
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
      s_type * gate_data_1 = (s_type *) malloc(gate_vector_1.size()*sizeof(s_type));
      for (size_t p=0; p<gate_vector_1.size(); ++p)
        *(gate_data_1+p) = gate_vector_1[p];
      vector<int> dims({DIM,DIM,DIM});
      grid_of_lists_of_tensors[i_j_1[0]][i_j_1[1]].push_back(
        talsh::Tensor(dims, gate_data_1));
      gate_data_1 = NULL;
      grid_of_lists_of_dims[i_j_1[0]][i_j_1[1]].push_back(dims);
      vector<int> bond_type;
      if (i_j_1[0]<i_j_2[0]) bond_type = {2, 0};
      else if (i_j_1[1]<i_j_2[1]) bond_type = {3, 1};
      else if (i_j_1[0]>i_j_2[0]) bond_type = {0, 2};
      else if (i_j_1[1]>i_j_2[1]) bond_type = {1, 3};
      grid_of_lists_of_bonds[i_j_1[0]][i_j_1[1]].push_back(bond_type[0]);
      // Gate 2
      vector<s_type> gate_vector_2(_GATES_DATA.at("cz_q2"));
      s_type * gate_data_2 = (s_type *) malloc(gate_vector_2.size()*sizeof(s_type));
      for (size_t p=0; p<gate_vector_2.size(); ++p)
        *(gate_data_2+p) = gate_vector_2[p];
      grid_of_lists_of_tensors[i_j_2[0]][i_j_2[1]].push_back(
        talsh::Tensor(dims, gate_data_2));
      gate_data_2 = NULL;
      grid_of_lists_of_dims[i_j_2[0]][i_j_2[1]].push_back(dims);
      grid_of_lists_of_bonds[i_j_2[0]][i_j_2[1]].push_back(bond_type[1]);
    }

  }


  // Grid of tensors with out of order indexes, of final out of order dims,
  // of final out of order volumes, of out of order bond types,
  // of new to old bond positions, and of number of each bond type
  vector<vector<talsh::Tensor *>> grid_of_ooo_tensors(I);
  vector<vector<vector<int>>> grid_of_ooo_dims(I);
  vector<vector<size_t>> grid_of_ooo_volumes(I);
  vector<vector<vector<int>>> grid_of_ooo_bond_types(I);
  vector<vector<vector<int>>> grid_of_new_to_old_bond_positions(I);
  vector<vector<vector<int>>> grid_of_number_of_bond_types(I);
  for (int i=0; i<I; ++i)
  {
    grid_of_ooo_tensors[i] = vector<talsh::Tensor *>(J);
    grid_of_ooo_dims[i] = vector<vector<int>>(J);
    grid_of_ooo_volumes[i] = vector<size_t>(J);
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
    //cout << i << " " << j << endl;
    // Reuse pointer to data for all iterations
    vector<s_type *> data(2);
    data[0] = (s_type *) malloc((int)pow(2,20)*sizeof(s_type));
    data[1] = (s_type *) malloc((int)pow(2,20)*sizeof(s_type));
    vector<int> dims_S, dims_T;
    int idx = 0;
    // Assume there is at least 2 tensors per world-line
    for (int k=1; k<grid_of_lists_of_tensors[i][j].size(); ++k)
    {
      new_tensor_rank = grid_of_lists_of_tensors[i][j][k].getRank();
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
      dims_S = vector<int>(current_rank,DIM);
      dims_T = vector<int>(new_rank,DIM);
      talsh::Tensor S(dims_S, data[-idx+1]);
      talsh::Tensor T(dims_T, data[idx]);
      //cout << k << " " << pattern << endl;
      if (k==1)
      {
        TensContraction contraction(pattern, &T,
                                         &grid_of_lists_of_tensors[i][j][1],
                                         &grid_of_lists_of_tensors[i][j][0]);
        int errc = contraction.execute(DEV_HOST,0);
        assert(errc==TALSH_SUCCESS);
        assert(contraction.sync(DEV_HOST,0));
      }
      else
      {
        TensContraction contraction(pattern, &T,
                                         &grid_of_lists_of_tensors[i][j][k],
                                         &S);
        int errc = contraction.execute(DEV_HOST,0);
        assert(errc==TALSH_SUCCESS);
        assert(contraction.sync(DEV_HOST,0));
      }
      current_rank = new_rank;
      idx = -idx + 1;
    }
    grid_of_ooo_tensors[i][j] = new talsh::Tensor(dims_T, data[-idx+1]);
    grid_of_ooo_dims[i][j] = dims_T;
    grid_of_ooo_volumes[i][j] = volume_from_dims(dims_T);
    for (auto v : data)
      v = NULL;
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
    {
      if (v>=0)
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
  // Reorder tensors
  for (int q=0; q<num_qubits; ++q)
  {
    vector<int> i_j = _q_to_i_j(q, J);
    int i = i_j[0], j = i_j[1];
    if (find(off.begin(),off.end(), vector<int>({i,j}))!=off.end())
    { continue; }
    vector<int> dims_T(grid_of_ooo_dims[i][j].size());
    for (int p=0; p<dims_T.size(); ++p)
      dims_T[p] = grid_of_ooo_dims[i][j][grid_of_new_to_old_bond_positions[i][j][p]];
    s_type * data_T = (s_type *) malloc(grid_of_ooo_volumes[i][j]*sizeof(s_type));
    talsh::Tensor T(dims_T, data_T);
    vector<s_type> vector_I(_GATES_DATA.at("I"));
    s_type * data_I = (s_type *) malloc(vector_I.size()*sizeof(s_type));
    for (size_t p=0; p<vector_I.size(); ++p)
      *(data_I+p) = vector_I[p];
    vector<int> dims_I({DIM,DIM});
    talsh::Tensor I(dims_I, data_I);
    data_T = NULL;
    data_I = NULL;
    vector<string> pattern_vector_D({_ALPHABET[1]});
    for (auto v : grid_of_new_to_old_bond_positions[i][j])
      pattern_vector_D.push_back(_ALPHABET[v+2]);
    string pattern_L(pattern_from_vector(vector<string>(_ALPHABET.cbegin()+2,
                                         _ALPHABET.cbegin()+dims_T.size()+1)));
    pattern_L = _ALPHABET[0] + "," + pattern_L;
    string pattern_D(pattern_from_vector(pattern_vector_D));
    string pattern_R(pattern_from_vector(vector<string>(_ALPHABET.cbegin(),
                                         _ALPHABET.cbegin()+2)));
    string pattern = "D("+pattern_D+")+=L("+pattern_L+")*R("+pattern_R+")";
    cout << endl << pattern << endl;
    T.print();
    grid_of_ooo_tensors[i][j]->print();
    I.print();
    /*
    TensContraction contraction(pattern, &T, &(*grid_of_ooo_tensors[i][j]), &I);
    int errc = contraction.execute(DEV_HOST,0);
    assert(errc==TALSH_SUCCESS);
    assert(contraction.sync(DEV_HOST,0));
    */
  }
  

  /*
  // Insert deltas to last layer.
  idx = 0;
  for (int q=0; q<num_qubits; ++q)
  {
    vector<int> i_j = _q_to_i_j(q, J);
    int i = i_j[0], j = i_j[1];
    if (find(off.begin(),off.end(), vector<int>({i,j}))!=off.end())
    { continue; }
    if (find(A.begin(),A.end(),vector<int>({i,j}))!=A.end()) { continue; }
    string delta_gate = (final_conf_B[idx]=='0')?"delta_0":"delta_1";
    vector<s_type> gate_vector(_GATES_DATA.at(delta_gate));
    s_type * gate_data = (s_type *) malloc(gate_vector.size()*sizeof(s_type));
    for (size_t p=0; p<gate_vector.size(); ++p)
      *(gate_data+p) = gate_vector[p];
    vector<int> dims({DIM});
    grid_of_lists_of_tensors[i][j].push_back(talsh::Tensor(dims, gate_data));
    gate_data = NULL;
    grid_of_lists_of_dims[i][j].push_back(dims);
    grid_of_lists_of_bonds[i][j].push_back(-1);
    ++idx; // Move in B only.
  }
  */
  

}



#endif
