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
* @file mkl_tensor.cpp
* Implementation of the MKLTensor class.
* @see https://github.com/benjaminvillalonga/optimized_parallel_QC_with_TN
*
* @author Benjamin Villalonga
* @date Created: August 2018
* @date Modified: August 2018
*/

#include <algorithm>
#include <cmath>
#include <cassert>
#include <omp.h>

#include "mkl_tensor.h"

// Time
#include <ctime>
#include <chrono>
using namespace chrono;

/**
* Cache friendly size (for complex<float>) to move things around.
*/
#ifndef MAX_RIGHT_DIM
#define MAX_RIGHT_DIM 1024
#endif

/**
* Smallest size of cache friendly blocks (for complex<float>).
*/
#ifndef MIN_RIGHT_DIM
#define MIN_RIGHT_DIM 32
#endif

/**
* Global vector<string> with the alphabet.
*/
const vector<string> _ALPHABET({"a","b","c","d","e","f","g","h","i","j",
                                "k","l","m","n","o","p","q","r","s","t",
                                "u","v","w","x","y","z","A","B","C","D",
                                "E","F","G","H","I","J","K","L","M","N",
                                "O","P","Q","R","S","T","U","V","W","X",
                                "Y","Z"});

/**
* unordered_map<int,int> with the log2 of powers of 2 up to 2^30, in order to
* quickly switch to smart reordering and look up the logs.
*/
const unordered_map<int,int> _LOG_2({ { 2          , 1  },
                                      { 4          , 2  },
                                      { 8          , 3  },
                                      { 16         , 4  },
                                      { 32         , 5  },
                                      { 64         , 6  },
                                      { 128        , 7  },
                                      { 256        , 8  },
                                      { 512        , 9  },
                                      { 1024       , 10 },
                                      { 2048       , 11 },
                                      { 4096       , 12 },
                                      { 8192       , 13 },
                                      { 16384      , 14 },
                                      { 32768      , 15 },
                                      { 65536      , 16 },
                                      { 131072     , 17 },
                                      { 262144     , 18 },
                                      { 524288     , 19 },
                                      { 1048576    , 20 },
                                      { 2097152    , 21 },
                                      { 4194304    , 22 },
                                      { 8388608    , 23 },
                                      { 16777216   , 24 },
                                      { 33554432   , 25 },
                                      { 67108864   , 26 },
                                      { 134217728  , 27 },
                                      { 268435456  , 28 },
                                      { 536870912  , 29 },
                                      { 1073741824 , 30 },
                                      { 1073741824 , 30 } });

/**
* Global unordered_map<string,vector<int>> of reordering maps.
*/
unordered_map<string,vector<int>> _REORDER_MAPS;

/**
* Max size for dim_right in cache-friendly reordering.
*/
#define MAX_RIGHT_DIM 1024

/**
* Min size for dim_right in cache-friendly reordering.
*/
#define MIN_RIGHT_DIM 32

///////////////////////////// CLASS FUNCTIONS /////////////////////////////////

void MKLTensor::_init(const vector<string> & indices,
                      const vector<int> & dimensions)
{
  assert (indices.size()==dimensions.size()
          && "indices and dimensions should have equal size.");
  _indices = indices;
  _dimensions = dimensions;
  for (int i=0; i<_dimensions.size(); ++i)
    _index_to_dimension[indices[i]] = dimensions[i];
}

void MKLTensor::_clear()
{
  delete[] _data;
  _data = NULL;
}

void MKLTensor::_copy(const MKLTensor & other)
{
  if (_indices.empty())
  {
    _data = new s_type[other.size()];
  }
  _init(other.get_indices(), other.get_dimensions());
  #pragma omp parallel for schedule(static, MAX_RIGHT_DIM)
  for (int p=0; p<other.size(); ++p)
    *(_data+p) = *(other.data()+p);
}

MKLTensor::MKLTensor()
{
  _data = NULL;
}

MKLTensor::MKLTensor(vector<string> indices, vector<int> dimensions)
{
  _init(indices, dimensions);
  int total_dim = 1;
  for (int i=0; i<indices.size(); ++i)
    total_dim *= dimensions[i];
  _data = new s_type[total_dim];
}

MKLTensor::MKLTensor(vector<string> indices, vector<int> dimensions,
                     const vector<s_type> & data) :
                     MKLTensor(indices, dimensions)
{
  // Check that the data has the same length as this MKLTensor's size().
  int this_size = size();
  assert (this_size==data.size()
          && "The vector data has to match the size of the MKLTensor.");
  // Fill in the _data.
  for (int i=0; i<this_size; ++i)
    *(_data+i) = data[i];
}

MKLTensor::MKLTensor(vector<string> indices, vector<int> dimensions, s_type * data)
{
  _init(indices, dimensions);
  _data = data;
}

MKLTensor::MKLTensor(const MKLTensor & other)
{
  _copy(other);
}

MKLTensor::~MKLTensor()
{
  _clear();
}

const MKLTensor & MKLTensor::operator=(const MKLTensor & other)
{
  if (this != &other)
  {
    _copy(other);
  }
  return *this;
}

const vector<string> & MKLTensor::get_indices() const
{
  return _indices;
}

void MKLTensor::set_indices(const vector<string> & indices)
{
  _indices = indices;
}

const vector<int> & MKLTensor::get_dimensions() const
{
  return _dimensions;
}

void MKLTensor::set_dimensions(const vector<int> & dimensions)
{
  // Assert.
  if (_data)
  {
    int total_dim = 1;
    for (int i=0; i<dimensions.size(); ++i)
      total_dim *= dimensions[i];
    assert (size()>=total_dim
            && "The dimensions must match the size of the MKLTensor.");
  }
  _dimensions = dimensions;
}

void MKLTensor::set_indices_and_dimensions(const vector<string> & indices,
                                           const vector<int> & dimensions)
{
  // The following line takes care of the total size of the dimensions.
  set_dimensions(dimensions);
  _init(indices, dimensions);
}

const unordered_map<string,int> &
             MKLTensor::get_index_to_dimension() const
{
  return _index_to_dimension;
}

void MKLTensor::generate_index_to_dimension()
{
  for (int i=0; i<_indices.size(); ++i)
    _index_to_dimension[_indices[i]] = _dimensions[i];
}

int MKLTensor::size() const
{
  int total_dim = 1;
  for (int i=0; i<_dimensions.size(); ++i)
    total_dim *= _dimensions[i];
  return total_dim;
}

s_type * MKLTensor::data()
{
  return _data;
}

const s_type * MKLTensor::data() const
{
  return _data;
}

void MKLTensor::project(string index, int index_value,
                        MKLTensor & projection_tensor) const
{
  assert (index==_indices[0] && "index has to be indices[0].");
  assert ((index_value>=0 && index_value<_dimensions[0])
          && "index_value must be contained in [0, dimensions[0]).");
  s_type * projection_data = projection_tensor.data();
  int projection_size = projection_tensor.size();
  int projection_begin = projection_size * index_value;
  #pragma omp parallel for schedule(static, MAX_RIGHT_DIM)
  for (int p=0; p<projection_size; ++p)
    *(projection_data+p) = *(_data+projection_begin+p);
  vector<string> projection_indices(_indices.begin()+1, _indices.end());
  vector<int> projection_dimensions(_dimensions.begin()+1,
                                    _dimensions.end());
  projection_tensor.set_indices(projection_indices);
  projection_tensor.set_dimensions(projection_dimensions);
  projection_tensor.generate_index_to_dimension();
}

void MKLTensor::rename_index(string old_name, string new_name)
{
  auto it = find(_indices.begin(), _indices.end(), old_name);
  assert (it!=_indices.end() && "old_name has to be a valid index.");
  *it = new_name;
  _index_to_dimension[new_name] = _index_to_dimension[old_name];
  _index_to_dimension.erase(old_name);
}

void MKLTensor::bundle(vector<string> indices_to_bundle, string bundled_index)
{
  // Asserts.
  assert (_vector_s_in_vector_s(indices_to_bundle, _indices)
          && "indices_to_bundle has to be contained in indices.");
  vector<string> subtracted_indices(_vector_subtraction(_indices,
                                                        indices_to_bundle));
  vector<string> indices_to_bundled_original_order(
                            _vector_subtraction(_indices, subtracted_indices));
  assert (indices_to_bundled_original_order==indices_to_bundle
          && "indices_to_bundle must be in the order they appear in indices.");

  int bundled_dim = 1;
  for (int i=0; i<indices_to_bundle.size(); ++i)
  {
    bundled_dim *= _index_to_dimension[indices_to_bundle[i]];
    _index_to_dimension.erase(indices_to_bundle[i]);
  }
  _index_to_dimension[bundled_index] = bundled_dim;
  int bundled_idxpos = 0;
  for (int i=0; i<_indices.size(); ++i)
  {
    if (_string_in_vector(_indices[i], indices_to_bundle))
    {
      bundled_idxpos = i;
      break;
    }
  }
  vector<string> new_indices(subtracted_indices);
  new_indices.insert(new_indices.begin()+bundled_idxpos, bundled_index);
  vector<int> new_dimensions(new_indices.size());
  for (int i=0; i<new_dimensions.size(); ++i)
    new_dimensions[i] = _index_to_dimension[new_indices[i]];
  _indices = new_indices;
  _dimensions = new_dimensions;
}

void MKLTensor::_naive_reorder(vector<string> new_ordering,
                               s_type * scratch_copy)
{
  // Don't do anything if there is nothing to reorder.
  if (new_ordering==_indices)
    return;

  vector<string> old_ordering(_indices);
  vector<int> old_dimensions(_dimensions);
  int num_indices = old_ordering.size();
  int total_dim = size();

  // Create map_old_to_new_idxpos from old to new indices, and new_dimensions.
  vector<int> map_old_to_new_idxpos(num_indices);
  vector<int> new_dimensions(num_indices);
  for (int i=0; i<num_indices; ++i)
  {
    for (int j=0; j<num_indices; ++j)
    {
      if (old_ordering[i]==new_ordering[j])
      {
        map_old_to_new_idxpos[i] = j;
        new_dimensions[j] = old_dimensions[i];
        break;
      }
    }
  }

  // Create super dimensions (combined dimension of all to the right of i).
  vector<int> old_super_dimensions(num_indices);
  vector<int> new_super_dimensions(num_indices);
  old_super_dimensions[num_indices-1] = 1;
  new_super_dimensions[num_indices-1] = 1;
  for (int i=old_dimensions.size()-2; i>=0; --i)
  {
    old_super_dimensions[i] = old_super_dimensions[i+1] * old_dimensions[i+1];
    new_super_dimensions[i] = new_super_dimensions[i+1] * new_dimensions[i+1];
  }

  // Allocating small_map_old_to_new_position.
  vector<unsigned short int> small_map_old_to_new_position(MAX_RIGHT_DIM);

  // Start moving data around.
  // First copy all data into scratch.
  #pragma omp parallel for schedule(static, MAX_RIGHT_DIM)
  for (int p=0; p<total_dim; ++p)
    *(scratch_copy+p) = *(_data+p);

  // No combined efficient mapping from old to new positions with actual
  // copies in memory, all in small cache friendly (for old data, not new,
  // which could be very scattered) blocks.
  // Define i and j once for the whole iteration.
  int i, j;
  // Position old and new.
  int po = 0, pn;
  // Counter of the values of each indices in the iteration (old ordering).
  vector<int> old_counter(num_indices, 0);
  // offset is important when doing this in blocks, as it's indeed implemented.
  int offset = 0;
  // internal_po keeps track of interations within a block.
  // Blocks have size MAX_RIGHT_DIM.
  int internal_po = 0;
  // External loop loops over blocks.
  while(true)
  {
    // If end of entire opration, break.
    if (po==total_dim-1)
      break;

    internal_po = 0;
    // Each iteration of the while block goes through a new position.
    // Inside the while, j takes care of increasing indices properly.
    while(true)
    {
      po = 0;
      pn = 0;
      for (i=0; i<num_indices; ++i)
      {
        po += old_super_dimensions[i] * old_counter[i];
        pn += new_super_dimensions[map_old_to_new_idxpos[i]]*old_counter[i];
      }
      small_map_old_to_new_position[po-offset] = pn;
      for (j=num_indices-1; j>=0 ; --j)
      {
        if (++old_counter[j]<old_dimensions[j])
          break;
        else
          old_counter[j]=0;
      }
      // If end of block or end of entire operation, break.
      if ((++internal_po==MAX_RIGHT_DIM) || (po==total_dim-1))
        break;
      // If last index (0) was increased, then go back to fastest index.
      if (j<0)
        break;
    }
    // Copy data for this block, taking into account offset of small_map...
    for (int p=0; p<min(MAX_RIGHT_DIM,total_dim); ++p)
      *(_data+small_map_old_to_new_position[p]) = *(scratch_copy+offset+p);

    offset += MAX_RIGHT_DIM;
  }

  _init(new_ordering, new_dimensions);

  scratch_copy = NULL;
}

void MKLTensor::_fast_reorder(vector<string> new_ordering,
                              s_type * scratch_copy)
{
  // Create binary orderings.
  vector<string> old_ordering(_indices);
  vector<int> old_dimensions(_dimensions);
  int num_indices = old_ordering.size();
  int total_dim = 1;
  for (int i=0; i<num_indices; ++i)
    total_dim *= old_dimensions[i];
  // Create map_old_to_new_idxpos from old to new indices, and new_dimensions.
  vector<int> map_old_to_new_idxpos(num_indices);
  vector<int> new_dimensions(num_indices);
  for (int i=0; i<num_indices; ++i)
  {
    for (int j=0; j<num_indices; ++j)
    {
      if (old_ordering[i]==new_ordering[j])
      {
        map_old_to_new_idxpos[i] = j;
        new_dimensions[j] = old_dimensions[i];
        break;
      }
    }
  }
  // Create binary orderings:
  vector<int> old_logs(num_indices);
  for (int i=0; i<num_indices; ++i)
  {
    old_logs[i] = _LOG_2.at(old_dimensions[i]);
  }
  int num_binary_indices = _LOG_2.at(total_dim);
  // Create map from old letter to new group of letters.
  unordered_map<string,vector<string>> binary_groups;
  int alphabet_position = 0;
  for (int i=0; i<num_indices; ++i)
  {
    vector<string> group(old_logs[i]);
    for (int j=0; j<old_logs[i]; ++j)
    {
      group[j] = _ALPHABET[alphabet_position];
      ++alphabet_position;
    }
    binary_groups[old_ordering[i]] = group;
  }
  // Create old and new binary ordering in letters.
  vector<string> old_binary_ordering(num_binary_indices);
  vector<string> new_binary_ordering(num_binary_indices);
  int binary_position = 0;
  for (int i=0; i<num_indices; ++i)
  {
    string old_index = old_ordering[i];
    for (int j=0; j<binary_groups[old_index].size(); ++j)
    {
      old_binary_ordering[binary_position] = binary_groups[old_index][j];
      ++binary_position;
    }
  }
  binary_position = 0;
  for (int i=0; i<num_indices; ++i)
  {
    string new_index = new_ordering[i];
    for (int j=0; j<binary_groups[new_index].size(); ++j)
    {
      new_binary_ordering[binary_position] = binary_groups[new_index][j];
      ++binary_position;
    }
  }
  // Up to here, I have created old_binary_ordering and new_binary_ordering.

  // Change _indices and _dimensions, as well as _index_to_dimension.
  // This is common to all cases, special or default (worst case).
  _init(new_ordering, new_dimensions);

  // Now special cases, before the default L-R-L worst case.
  // Tensor doesn't have enough size to pass MAX_RIGHT_DIM => only one R.
  if (num_binary_indices<=_LOG_2.at(MAX_RIGHT_DIM))
  {
    _right_reorder(old_binary_ordering, new_binary_ordering,
                   num_binary_indices);
    scratch_copy = NULL;
    return;
  }
  // Reordering needs only one right move or one left move.
  // Left moves might benefit a lot from being applied on shorter strings,
  // up to L10. Computation times are L4>L5>L6>...>L10. I'll consider
  // all of these cases.
  {
    int Lr = _LOG_2.at(MAX_RIGHT_DIM);
    int Ll = new_binary_ordering.size() - Lr;
    int Rr = _LOG_2.at(MIN_RIGHT_DIM);
    int Rl = new_binary_ordering.size() - Rr;
    vector<string> Ll_old_indices(old_binary_ordering.begin(),
                                  old_binary_ordering.begin()+Ll);
    vector<string> Ll_new_indices(new_binary_ordering.begin(),
                                  new_binary_ordering.begin()+Ll);
    // Only one R10.
    if (Ll_old_indices==Ll_new_indices)
    {
      vector<string> Lr_old_indices(old_binary_ordering.begin()+Ll,
                                    old_binary_ordering.end());
      vector<string> Lr_new_indices(new_binary_ordering.begin()+Ll,
                                    new_binary_ordering.end());
      _right_reorder(Lr_old_indices, Lr_new_indices, Lr);
      scratch_copy = NULL;
      return;
    }
    // Only one L\nu move.
    for (int i=5; i>=-1; --i)
    {
      int extended_Rr = Rr + i;
      vector<string> Rr_old_indices(old_binary_ordering.end()-extended_Rr,
                                    old_binary_ordering.end());
      vector<string> Rr_new_indices(new_binary_ordering.end()-extended_Rr,
                                    new_binary_ordering.end());
      if (Rr_old_indices==Rr_new_indices)
      {
        vector<string> Rl_old_indices(old_binary_ordering.begin(),
                                      old_binary_ordering.end()-extended_Rr);
        vector<string> Rl_new_indices(new_binary_ordering.begin(),
                                      new_binary_ordering.end()-extended_Rr);
        _left_reorder(Rl_old_indices, Rl_new_indices, extended_Rr,
                      scratch_copy);
        scratch_copy = NULL;
        return;
      }
    }
  }
  

  // Worst case.
  {
    // There are two boundaries, L and R.
    // The worst case is the following. It can be optimized, in order to do
    // work early and maybe save the later steps. Think about that, but first
    // let's have something that already works:
    // 1) L5 All indices that are to the left of R and need to end up to its
    //    right are placed in the bucket.
    // 2) R10 All indices to the right of R are placed in their final ordering.
    // 3) L5 All indices to the left of L are placed in their final ordering.
    // Then hardcode special cases.
    // Add conditional to _left_reorder and _right_reorder, so that they don't
    // do anything when not needed.
    // Debug from here!
    int Lr = _LOG_2.at(MAX_RIGHT_DIM);
    int Ll = new_binary_ordering.size() - Lr;
    int Rr = _LOG_2.at(MIN_RIGHT_DIM);
    int Rl = new_binary_ordering.size() - Rr;
    // Helper vectors that can be reused.
    vector<string> Lr_indices(Lr), Ll_indices(Ll), Rr_indices(Rr),
                   Rl_indices(Rl);
    for (int i=0; i<Rr; ++i)
      Rr_indices[i] = new_binary_ordering[i+Rl];
    for (int i=0; i<Rl; ++i)
      Rl_indices[i] = old_binary_ordering[i];
    vector<string> Rr_new_in_Rl_old =
                    _vector_intersection(Rl_indices, Rr_indices);
    vector<string> Rl_old_not_in_Rr_new =
                    _vector_subtraction(Rl_indices, Rr_new_in_Rl_old);
    vector<string> Rl_first_step = _vector_concatenation(Rl_old_not_in_Rr_new,
                                                         Rr_new_in_Rl_old);
    vector<string> Rl_zeroth_step(Rl);
    for (int i=0; i<Rl; ++i)
      Rl_zeroth_step[i] = old_binary_ordering[i];
    _left_reorder(Rl_zeroth_step, Rl_first_step, Rr, scratch_copy);
    // Done with 1).
    // Let's go with 2).
    vector<string> Lr_first_step = _vector_concatenation(
      vector<string>(Rl_first_step.begin()+Ll, Rl_first_step.end()),
      vector<string>(old_binary_ordering.begin()+Rl,
                     old_binary_ordering.end()));
    Rr_indices = vector<string>(new_binary_ordering.begin()+Rl,
                                new_binary_ordering.end());
    vector<string> Lr_second_step = _vector_concatenation(
      _vector_subtraction(Lr_first_step, Rr_indices),
      vector<string>(Rr_indices));
    _right_reorder(Lr_first_step, Lr_second_step, Lr);
    // Done with 2).
    // Let's go with 3).
    vector<string> Rl_second_step = _vector_concatenation(
      vector<string>(Rl_first_step.begin(), Rl_first_step.begin()+Ll),
      vector<string>(Lr_second_step.begin(), Lr_second_step.begin()+Lr-Rr));
    vector<string> Rl_thrid_step(new_binary_ordering.begin(),
                                 new_binary_ordering.begin()+Rl);
    _left_reorder(Rl_second_step, Rl_thrid_step, Rr, scratch_copy);
    // done with 3).

    scratch_copy = NULL;
  }
}

// Assuming all indices are binary for old_ordering and new_ordering.
// old_ordering and new_ordering refer to the right.
void MKLTensor::_right_reorder(const vector<string> & old_ordering,
                               const vector<string> & new_ordering,
                               int num_indices_right)
{
  // Don't do anything if there is nothing to reorder.
  if (new_ordering==old_ordering)
    return;

  // Create dim, num_indices, map_old_to_new_idxpos from old to new indices,
  // old_dimensions, new_dimensions, and total_dim.
  int dim = 2;
  int num_indices = old_ordering.size();
  vector<int> map_old_to_new_idxpos(num_indices);
  vector<int> old_dimensions(num_indices, dim);
  vector<int> new_dimensions(num_indices, dim);
  int total_dim = 1;
  for (int i=0; i<num_indices; ++i)
    total_dim *= old_dimensions[i];
  for (int i=0; i<num_indices; ++i)
  {
    for (int j=0; j<num_indices; ++j)
    {
      if (old_ordering[i]==new_ordering[j])
      {
        map_old_to_new_idxpos[i] = j;
        new_dimensions[j] = old_dimensions[i];
        break;
      }
    }
  }

  // Create the map_old_to_new_position, or get a reference to it if it exists
  // on _REORDER_MAPS.
  string name = _reordering_to_string(map_old_to_new_idxpos, old_dimensions);
  if (_REORDER_MAPS.find(name)==_REORDER_MAPS.end())
  {
    _REORDER_MAPS[name] = vector<int>(total_dim);
    _generate_binary_reordering_map(map_old_to_new_idxpos,
                                    _REORDER_MAPS.at(name));
  }
  const vector<int> & map_old_to_new_position = _REORDER_MAPS.at(name);

  // With the map_old_to_new_position, we are ready to reorder within
  // small chuncks.
  int dim_right = total_dim;
  int dim_left = size() / dim_right; // Remember, it's all powers of 2, so OK.
  #pragma omp parallel
  {
    // For some reason, allocating these spaces and using them is about 2
    // times faster than bringing a pointer to a scratch space and using
    // different chunks of it.
    s_type * temp_data = new s_type[dim_right];
    #pragma omp for schedule(static)
    for (int pl=0; pl<dim_left; ++pl)
    {
      int current_thread = omp_get_thread_num();
      int offset = pl * dim_right;
      for (int pr=0; pr<dim_right; ++pr)
        *(temp_data+pr) = *(_data+offset+pr);
      for (int pr=0; pr<dim_right; ++pr)
        *(_data+offset+map_old_to_new_position[pr]) = *(temp_data+pr);
    }
    delete[] temp_data;
    temp_data = NULL;
  }
}

// Assuming all indices are binary for old_ordering and new_ordering.
// old_ordering and new_ordering refer to the left.
void MKLTensor::_left_reorder(const vector<string> & old_ordering,
                              const vector<string> & new_ordering,
                              int num_indices_right, s_type * scratch_copy)
{
  // Don't do anything if there is nothing to reorder.
  if (new_ordering==old_ordering)
    return;

  // Create dim, num_indices, map_old_to_new_idxpos from old to new indices,
  // old_dimensions, new_dimensions, and total_dim.
  int dim = 2;
  int num_indices = old_ordering.size();
  vector<int> map_old_to_new_idxpos(num_indices);
  vector<int> old_dimensions(num_indices, dim);
  vector<int> new_dimensions(num_indices, dim);
  int total_dim = 1;
  for (int i=0; i<num_indices; ++i)
    total_dim *= old_dimensions[i];
  for (int i=0; i<num_indices; ++i)
  {
    for (int j=0; j<num_indices; ++j)
    {
      if (old_ordering[i]==new_ordering[j])
      {
        map_old_to_new_idxpos[i] = j;
        new_dimensions[j] = old_dimensions[i];
        break;
      }
    }
  }
  // total_dim is the total dimension of the left!

  // Create the map_old_to_new_position, or get a reference to it if it exists
  // on _REORDER_MAPS.
  string name = _reordering_to_string(map_old_to_new_idxpos, old_dimensions);
  if (_REORDER_MAPS.find(name)==_REORDER_MAPS.end())
  {
    _REORDER_MAPS[name] = vector<int>(total_dim);
    _generate_binary_reordering_map(map_old_to_new_idxpos,
                                    _REORDER_MAPS.at(name));
  }
  const vector<int> & map_old_to_new_position = _REORDER_MAPS.at(name);

  // With the map_old_to_new_position, we are ready to move small chunks.
  int dim_left = total_dim;
  int tensor_dim = size();
  int dim_right = tensor_dim / dim_left; // Remember, it's all powers of 2,
                                         // so OK.
  // Copy.
  #pragma omp parallel for schedule(static, MAX_RIGHT_DIM)
  for (int p=0; p<tensor_dim; ++p)
  {
      *(scratch_copy+p) = *(_data+p);
  }
  // Move back.
  #pragma omp parallel
  {
    #pragma omp for schedule(static)
    for (int pl=0; pl<dim_left; ++pl)
    {
      int old_offset = pl * dim_right;
      int new_offset = map_old_to_new_position[pl] * dim_right;
      for (int pr=0; pr<dim_right; ++pr)
      {
        *(_data+new_offset+pr) = *(scratch_copy+old_offset+pr);
      }
    }
  }
  scratch_copy = NULL;
}

void MKLTensor::reorder(vector<string> new_ordering, s_type * scratch_copy)
{
  bool fast = true;
  for (int i=0; i<_dimensions.size(); ++i) 
  {
    if (_LOG_2.find(_dimensions[i])==_LOG_2.end())
    {
      fast = false;
      break;
    }
  }
  if (fast)
  {
    _fast_reorder(new_ordering, scratch_copy);
  }
  else
  {
    _naive_reorder(new_ordering, scratch_copy);
  }
}

void MKLTensor::scalar_multiply(s_type scalar)
{
  int this_size = size();
  #pragma omp parallel for schedule(static, MAX_RIGHT_DIM)
  for (int p=0; p<this_size; ++p)
    *(_data+p) = (*(_data+p)) * scalar;
}

double MKLTensor::tensor_norm() const
{
  double total_norm(0.0);
  int this_size = size();
  for (int p=0; p<this_size; ++p)
  {
    total_norm += norm(*(_data+p));
  }
  return total_norm;
}

int MKLTensor::num_zeros() const
{
  int count(0);
  int this_size = size();
  s_type complex_0(0.);
  for (int p=0; p<this_size; ++p)
  {
    if (*(_data+p)==complex_0)
      ++count;
  }
  return count;
}

void MKLTensor::print() const
{
  string print_str("MKLTensor of rank ");
  print_str += to_string(_indices.size());
  print_str += ": ";
  for (int i=0; i<_indices.size(); ++i)
  {
    print_str += _indices[i];
    print_str += " -> ";
    print_str += to_string(_dimensions[i]);
    if (i<_indices.size()-1) print_str += ", ";
  }
  cout << print_str << endl;
}

void MKLTensor::print_data() const
{
  for (int p=0; p<size(); ++p)
    cout << *(_data+p) << " ";
  cout << endl;
}

/////////////////////////// EXTERNAL FUNCTIONS ////////////////////////////////

// use mkl if complexity < some value.
void _multiply_MM(const s_type * A_data, const s_type * B_data,
                 s_type * C_data, int m, int n, int k)
{
  s_type alpha = 1.0;
  s_type beta = 0.0;
  cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, &alpha,
              A_data, max(1,k), B_data, max(1,n), &beta, C_data, max(1,n));
}

void _multiply_Mv(const s_type * A_data, const s_type * B_data,
                 s_type * C_data, int m, int k)
{
  s_type alpha = 1.0;
  s_type beta = 0.0;
  cblas_cgemv(CblasRowMajor, CblasNoTrans, m, k, &alpha,
              A_data, max(1,k), B_data, 1, &beta, C_data, 1);
}

void _multiply_vM(const s_type * A_data, const s_type * B_data,
                 s_type * C_data, int n, int k)
{
  s_type alpha = 1.0;
  s_type beta = 0.0;
  cblas_cgemv(CblasRowMajor, CblasTrans, k, n, &alpha,
              A_data, max(1,n), B_data, 1, &beta, C_data, 1);
}

void _multiply_vv(const s_type * A_data, const s_type * B_data,
                 s_type * C_data, int k)
{
  cblas_cdotu_sub(k, A_data, 1, B_data, 1, C_data);
}

void multiply(MKLTensor & A, MKLTensor & B,
                         MKLTensor & C, s_type * scratch_copy)
{
  high_resolution_clock::time_point t0, t1;
  duration<double> time_span;
  t0 = high_resolution_clock::now();
  // Define left_indices, left_dim, right_indices, right_dim, and
  // common_indices, common_dim. Also C_size.
  vector<string> left_indices = _vector_subtraction(A.get_indices(),
                                                    B.get_indices());
  vector<string> right_indices = _vector_subtraction(B.get_indices(),
                                                     A.get_indices());
  vector<string> common_indices = _vector_intersection(A.get_indices(),
                                                       B.get_indices());
  int left_dim = 1, right_dim = 1, common_dim = 1;
  for (int i=0; i<left_indices.size(); ++i)
  {
    left_dim *= A.get_index_to_dimension().at(left_indices[i]);
  }
  for (int i=0; i<right_indices.size(); ++i)
  {
    right_dim *= B.get_index_to_dimension().at(right_indices[i]);
  }
  for (int i=0; i<common_indices.size(); ++i)
  {
    common_dim *= A.get_index_to_dimension().at(common_indices[i]);
  }
  t1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t1 - t0);
  //cout << "Time preparing variables: " << time_span.count() << "s\n";

  // Assert.
  assert ((left_dim*right_dim<=C.size())
          && "C doesn't have enough space for the product of A*B.");

  // Reorder.
  t0 = high_resolution_clock::now();
  vector<string> A_new_ordering = _vector_union(left_indices, common_indices);
  A.reorder(A_new_ordering, scratch_copy);
  t1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t1 - t0);
  //cout << "R " << time_span.count() << " s\n";
  //cout << "Time reordering A: " << time_span.count() << "s\n";
  t0 = high_resolution_clock::now();
  vector<string> B_new_ordering = _vector_union(common_indices, right_indices);
  B.reorder(B_new_ordering, scratch_copy);
  t1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t1 - t0);
  //cout << "R " << time_span.count() << " s\n";
  //cout << "Time reordering B: " << time_span.count() << "s\n";

  // Multiply. Four cases: MxM, Mxv, vxM, vxv.
  t0 = high_resolution_clock::now();
  if (left_indices.size()>0 && right_indices.size()>0)
  {
    _multiply_MM(A.data(), B.data(), C.data(), left_dim, right_dim,
                 common_dim);
  }
  else if (left_indices.size()>0 && right_indices.size()==0)
  {
    _multiply_Mv(A.data(), B.data(), C.data(), left_dim, common_dim);
  }
  else if (left_indices.size()==0 && right_indices.size()>0)
  {
    // Very import to switch A and B to use cgemv with transpose for this case.
    _multiply_vM(B.data(), A.data(), C.data(), right_dim, common_dim);
  }
  else if (left_indices.size()==0 && right_indices.size()==0)
  {
    _multiply_vv(A.data(), B.data(), C.data(), common_dim);
  }
  t1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t1 - t0);
  //cout << "M " << time_span.count() << " s\n";
  //cout << "Time multiplying A*B: " << time_span.count() << "s\n";

  // Set indices and dimensions of C.
  t0 = high_resolution_clock::now();
  vector<string> C_indices = _vector_union(left_indices, right_indices);
  vector<int> C_dimensions(C_indices.size());
  for (int i=0; i<left_indices.size(); ++i)
    C_dimensions[i] = A.get_index_to_dimension().at(left_indices[i]);
  for (int i=0; i<right_indices.size(); ++i)
    C_dimensions[i+left_indices.size()] =
                             B.get_index_to_dimension().at(right_indices[i]);
  C.set_indices(C_indices);
  C.set_dimensions(C_dimensions);
  C.generate_index_to_dimension();
  t1 = high_resolution_clock::now();
  time_span = duration_cast<duration<double>>(t1 - t0);
  //cout << "Time updating C's variables: " << time_span.count() << "s\n";
}

// Split it in parts, as before? Nah, it was all about generating small maps.
// Here I am generating THE map. This is done only once per map anyway.
void _generate_binary_reordering_map(const vector<int>& map_old_to_new_idxpos,
                                     vector<int> & map_old_to_new_position)
{
  int dim = 2; // Hard coded!
  int num_indices = map_old_to_new_idxpos.size();
  // Assert
  assert ((int)pow(dim,num_indices)==map_old_to_new_position.size()
          && "Size of map has to be equal to 2^num_indices");

  // Define super dimensions. See _naive_reorder().
  vector<int> old_dimensions(num_indices, dim);
  vector<int> new_dimensions(num_indices, dim);
  vector<int> old_super_dimensions(num_indices);
  vector<int> new_super_dimensions(num_indices);
  old_super_dimensions[num_indices-1] = 1;
  new_super_dimensions[num_indices-1] = 1;
  for (int i=num_indices-2; i>=0; --i)
  {
    old_super_dimensions[i] = old_super_dimensions[i+1] * dim;
    new_super_dimensions[i] = new_super_dimensions[i+1] * dim;
  }

  // Iterate and generate map.
  vector<int> old_counter(num_indices, 0);
  int po, pn; // Position of the data, old and new.
  int i, j;
  while(true)
  {
      po = 0;
      pn = 0;
      for (i=0; i<num_indices; ++i)
      {
        po += old_super_dimensions[i] * old_counter[i];
        pn += new_super_dimensions[map_old_to_new_idxpos[i]] * old_counter[i];
      }
      map_old_to_new_position[po] = pn;
      for (j=num_indices-1; j>=0 ; --j)
      {
        if (++old_counter[j]<old_dimensions[j])
          break;
        else
          old_counter[j]=0;
      }
      if (j<0)
        break;
  }
}

string _reordering_to_string(const vector<int> & map_old_to_new_idxpos,
                             const vector<int> & old_dimensions)
{
  int num_indices = map_old_to_new_idxpos.size();
  string name(""); 
  for (int i=0; i<num_indices; ++i)
    name += _ALPHABET[i];
  name += "->";
  for (int i=0; i<num_indices; ++i)
    name += _ALPHABET[map_old_to_new_idxpos[i]];
  for (int i=0; i<num_indices; ++i)
    name += "," + to_string(old_dimensions[i]);
  return name;
}

bool _string_in_vector(const string & s, const vector<string> & v)
{
  if (find(v.cbegin(), v.cend(), s) != v.cend())
    return true;
  else
    return false;
}

bool _vector_s_in_vector_s(const vector<string> & v, const vector<string> & w)
{
  for (int i=0; i<v.size(); ++i)
  {
    if (!_string_in_vector(v[i], w))
      return false;
  }
  return true;
}

vector<string> _vector_intersection(const vector<string> & v,
                                    const vector<string> & w)
{
  vector<string> temp;
  for (auto it=v.cbegin(); it!=v.cend(); ++it)
  {
    if (_string_in_vector(*it, w))
      temp.push_back(*it);
  }
  return temp;
}

vector<string> _vector_union(const vector<string> & v,
                             const vector<string> & w)
{
  vector<string> temp(v);
  for (auto it=w.cbegin(); it!=w.cend(); ++it)
  {
    if (!_string_in_vector(*it, v))
      temp.push_back(*it);
  }
  return temp;
}

vector<string> _vector_subtraction(const vector<string> & v,
                                   const vector<string> & w)
{
  vector<string> temp;
  for (auto it=v.cbegin(); it!=v.cend(); ++it)
  {
    if (!_string_in_vector(*it, w))
      temp.push_back(*it);
  }
  return temp;
}

vector<string> _vector_concatenation(const vector<string> & v,
                                     const vector<string> & w)
{
  vector<string> temp(v.size()+w.size());
  for (int i=0; i<v.size(); ++i)
    temp[i] = v[i];
  for (int i=0; i<w.size(); ++i)
    temp[i+v.size()] = w[i];
  return temp;
}
