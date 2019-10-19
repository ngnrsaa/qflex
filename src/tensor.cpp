/**
 * @file tensor.cpp
 * Implementation of the Tensor class.
 * @see https://github.com/benjaminvillalonga/optimized_parallel_QC_with_TN
 *
 * @author Benjamin Villalonga (main contributor), Bron Nelson, Sergio Boixo and
 * Salvatore Mandra
 * @date Created: August 2018
 * @date Modified: August 2018
 *
 * @copyright: Copyright © 2019, United States Government, as represented
 * by the Administrator of the National Aeronautics and Space Administration.
 * All rights reserved.
 * @licence: Apache License, Version 2.0
 */

#include "tensor.h"

#include <omp.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iterator>

// Time
#include <chrono>
#include <ctime>

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

namespace qflex {

// clang-format off
/**
 * Global vector<string> with the alphabet.
 */
const std::vector<std::string> _ALPHABET(
    {"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m",
     "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z",
     "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M",
     "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"});

/**
 * unordered_map<int,int> with the log2 of powers of 2 up to 2^30, in order to
 * quickly switch to smart reordering and look up the logs.
 */
const std::unordered_map<int, int> _LOG_2(
    {{2, 1},          {4, 2},           {8, 3},          {16, 4},
     {32, 5},         {64, 6},          {128, 7},        {256, 8},
     {512, 9},        {1024, 10},       {2048, 11},      {4096, 12},
     {8192, 13},      {16384, 14},      {32768, 15},     {65536, 16},
     {131072, 17},    {262144, 18},     {524288, 19},    {1048576, 20},
     {2097152, 21},   {4194304, 22},    {8388608, 23},   {16777216, 24},
     {33554432, 25},  {67108864, 26},   {134217728, 27}, {268435456, 28},
     {536870912, 29}, {1073741824, 30}, {1073741824, 30}});
// clang-format on

/**
 * Global unordered_map<string,vector<int>> of reordering maps.
 */
std::unordered_map<std::string, std::vector<int>> _REORDER_MAPS;

///////////////////////////// CLASS FUNCTIONS /////////////////////////////////

void Tensor::_init(const std::vector<std::string>& indices,
                   const std::vector<size_t>& dimensions) {
  if (indices.size() != dimensions.size()) {
    std::cout << "The number of indices: " << indices.size()
              << ", and number of dimensions: " << dimensions.size()
              << ", should be equal." << std::endl;
    assert(indices.size() == dimensions.size());
  }
  _indices = indices;
  _dimensions = dimensions;
  for (int i = 0; i < _dimensions.size(); ++i)
    _index_to_dimension[indices[i]] = dimensions[i];
}

void Tensor::_clear() {
  delete[] _data;
  _data = NULL;
}

void Tensor::_copy(const Tensor& other) {
  if (_indices.empty()) {
    _capacity = other.size();
    _data = new s_type[_capacity];
  } else {
    // The following line takes care of the total size of the dimensions.
    set_dimensions(other.get_dimensions());
  }
  _init(other.get_indices(), other.get_dimensions());
#pragma omp parallel for schedule(static, MAX_RIGHT_DIM)
  for (size_t p = 0; p < other.size(); ++p) *(_data + p) = *(other.data() + p);
}

Tensor::Tensor() { _data = NULL; }

Tensor::Tensor(std::vector<std::string> indices,
               std::vector<size_t> dimensions) {
  _init(indices, dimensions);
  _capacity = size();
  _data = new s_type[_capacity];
}

Tensor::Tensor(std::vector<std::string> indices, std::vector<size_t> dimensions,
               const std::vector<s_type>& data)
    : Tensor(indices, dimensions) {
  // Check that the data has the same length as this Tensor's size().
  size_t this_size = size();
  if (this_size != data.size()) {
    std::cout << "The vector data size: " << data.size()
              << ", has to match the size of the Tensor: " << this_size << "."
              << std::endl;
    assert(this_size == data.size());
  }
  _capacity = this_size;
  // Fill in the _data.
  for (size_t i = 0; i < this_size; ++i) *(_data + i) = data[i];
}

Tensor::Tensor(std::vector<std::string> indices, std::vector<size_t> dimensions,
               s_type* data) {
  if (data == nullptr) {
    std::cout << "Data must be non-null." << std::endl;
    assert(data != nullptr);
  }

  _init(indices, dimensions);
  _capacity = size();
  _data = data;
}

Tensor::Tensor(const Tensor& other) { _copy(other); }

Tensor::~Tensor() { _clear(); }

const Tensor& Tensor::operator=(const Tensor& other) {
  if (this != &other) {
    _copy(other);
  }
  return *this;
}

const std::vector<std::string>& Tensor::get_indices() const { return _indices; }

void Tensor::set_indices(const std::vector<std::string>& indices) {
  _indices = indices;
}

const std::vector<size_t>& Tensor::get_dimensions() const {
  return _dimensions;
}

void Tensor::set_dimensions(const std::vector<size_t>& dimensions) {
  // Assert.
  if (_data) {
    size_t total_dim = 1;
    for (int i = 0; i < dimensions.size(); ++i) total_dim *= dimensions[i];
    if (capacity() < total_dim) {
      std::cout << "The total allocated space: " << capacity()
                << ", is insufficient for the requested tensor dimensions: "
                << total_dim << "." << std::endl;
      assert(capacity() >= total_dim);
    }
  }
  _dimensions = dimensions;
}

void Tensor::set_indices_and_dimensions(const std::vector<std::string>& indices,
                                        const std::vector<size_t>& dimensions) {
  // The following line takes care of the total size of the dimensions.
  set_dimensions(dimensions);
  _init(indices, dimensions);
}

const std::unordered_map<std::string, size_t>& Tensor::get_index_to_dimension()
    const {
  return _index_to_dimension;
}

void Tensor::generate_index_to_dimension() {
  for (int i = 0; i < _indices.size(); ++i)
    _index_to_dimension[_indices[i]] = _dimensions[i];
}

size_t Tensor::size() const {
  size_t total_dim = 1;
  for (int i = 0; i < _dimensions.size(); ++i) total_dim *= _dimensions[i];
  return total_dim;
}

size_t Tensor::capacity() const { return _capacity; }

s_type* Tensor::data() { return _data; }

const s_type* Tensor::data() const { return _data; }

void Tensor::project(std::string index, size_t index_value,
                     Tensor& projection_tensor) const {
  if (index != _indices[0]) {
    std::cout << "Index: '" << index << "' has to be equal to indices[0]: '"
              << _indices[0] << "'." << std::endl;
    assert(index == _indices[0]);
  }
  if (index_value < 0 || index_value >= _dimensions[0]) {
    std::cout << "index_value: " << index_value << " must be contained in [0, "
              << _dimensions[0] << ")." << std::endl;
    assert((index_value >= 0 && index_value < _dimensions[0]));
  }
  // Resize projection_tensor first.
  std::vector<std::string> projection_indices(_indices.begin() + 1,
                                              _indices.end());
  std::vector<size_t> projection_dimensions(_dimensions.begin() + 1,
                                            _dimensions.end());
  projection_tensor.set_indices(projection_indices);
  projection_tensor.set_dimensions(projection_dimensions);
  projection_tensor.generate_index_to_dimension();

  // Fill projection_tensor with result of projection.
  s_type* projection_data = projection_tensor.data();
  int projection_size = projection_tensor.size();
  int projection_begin = projection_size * index_value;
#pragma omp parallel for schedule(static, MAX_RIGHT_DIM)
  for (size_t p = 0; p < projection_size; ++p)
    *(projection_data + p) = *(_data + projection_begin + p);
}

void Tensor::rename_index(std::string old_name, std::string new_name) {
  auto it = find(_indices.begin(), _indices.end(), old_name);
  if (it == _indices.end()) {
    std::cout << "old_name: " << old_name << ", has to be a valid index."
              << std::endl;
    assert(it != _indices.end());
  }
  bool new_name_is_existing_index =
      (find(_indices.begin(), _indices.end(), new_name) != _indices.end());
  if (new_name_is_existing_index) {
    std::cout << "new_name: " << new_name << ", cannot be an existing index."
              << std::endl;
    assert(!new_name_is_existing_index);
  }
  *it = new_name;
  _index_to_dimension[new_name] = _index_to_dimension[old_name];
  _index_to_dimension.erase(old_name);
}

void Tensor::bundle(std::vector<std::string> indices_to_bundle,
                    std::string bundled_index) {
  // Asserts.
  bool indices_to_bundle_in_indices =
      _vector_s_in_vector_s(indices_to_bundle, _indices);
  if (!indices_to_bundle_in_indices) {
    std::cout << "indices_to_bundle: "
              << _string_vector_to_string(indices_to_bundle)
              << " has to be contained in indices: "
              << _string_vector_to_string(_indices) << "." << std::endl;
    assert(indices_to_bundle_in_indices);
  }
  std::vector<std::string> subtracted_indices(
      _vector_subtraction(_indices, indices_to_bundle));
  std::vector<std::string> indices_to_bundled_original_order(
      _vector_subtraction(_indices, subtracted_indices));
  if (indices_to_bundled_original_order != indices_to_bundle) {
    std::cout << "indices_to_bundle: "
              << _string_vector_to_string(indices_to_bundle)
              << " must be in its original order: "
              << _string_vector_to_string(indices_to_bundled_original_order)
              << "." << std::endl;
    assert(indices_to_bundled_original_order == indices_to_bundle);
  }

  int bundled_dim = 1;
  for (int i = 0; i < indices_to_bundle.size(); ++i) {
    bundled_dim *= _index_to_dimension[indices_to_bundle[i]];
    _index_to_dimension.erase(indices_to_bundle[i]);
  }
  _index_to_dimension[bundled_index] = bundled_dim;
  int bundled_idxpos = 0;
  for (int i = 0; i < _indices.size(); ++i) {
    if (_string_in_vector(_indices[i], indices_to_bundle)) {
      bundled_idxpos = i;
      break;
    }
  }
  std::vector<std::string> new_indices(subtracted_indices);
  new_indices.insert(new_indices.begin() + bundled_idxpos, bundled_index);
  std::vector<size_t> new_dimensions(new_indices.size());
  for (int i = 0; i < new_dimensions.size(); ++i)
    new_dimensions[i] = _index_to_dimension[new_indices[i]];
  _indices = new_indices;
  _dimensions = new_dimensions;
}

void Tensor::_naive_reorder(std::vector<std::string> new_ordering,
                            s_type* scratch_copy) {
  if (scratch_copy == nullptr) {
    std::cout << "Scratch copy must be non-null." << std::endl;
    assert(scratch_copy != nullptr);
  }

  // Don't do anything if there is nothing to reorder.
  if (new_ordering == _indices) return;

  std::vector<std::string> old_ordering(_indices);
  std::vector<size_t> old_dimensions(_dimensions);
  int num_indices = old_ordering.size();
  size_t total_dim = size();

  // Create map_old_to_new_idxpos from old to new indices, and new_dimensions.
  std::vector<int> map_old_to_new_idxpos(num_indices);
  std::vector<size_t> new_dimensions(num_indices);
  for (int i = 0; i < num_indices; ++i) {
    for (int j = 0; j < num_indices; ++j) {
      if (old_ordering[i] == new_ordering[j]) {
        map_old_to_new_idxpos[i] = j;
        new_dimensions[j] = old_dimensions[i];
        break;
      }
    }
  }

  // Create super dimensions (combined dimension of all to the right of i).
  std::vector<size_t> old_super_dimensions(num_indices);
  std::vector<size_t> new_super_dimensions(num_indices);
  old_super_dimensions[num_indices - 1] = 1;
  new_super_dimensions[num_indices - 1] = 1;
  for (int i = old_dimensions.size() - 2; i >= 0; --i) {
    old_super_dimensions[i] =
        old_super_dimensions[i + 1] * old_dimensions[i + 1];
    new_super_dimensions[i] =
        new_super_dimensions[i + 1] * new_dimensions[i + 1];
  }

  // Allocating small_map_old_to_new_position.
  std::vector<unsigned short int> small_map_old_to_new_position(MAX_RIGHT_DIM);

// Start moving data around.
// First copy all data into scratch.
#pragma omp parallel for schedule(static, MAX_RIGHT_DIM)
  for (size_t p = 0; p < total_dim; ++p) *(scratch_copy + p) = *(_data + p);

  // No combined efficient mapping from old to new positions with actual
  // copies in memory, all in small cache friendly (for old data, not new,
  // which could be very scattered) blocks.
  // Define i and j once for the whole iteration.
  int i, j;
  // Position old and new.
  int po = 0, pn;
  // Counter of the values of each indices in the iteration (old ordering).
  std::vector<size_t> old_counter(num_indices, 0);
  // offset is important when doing this in blocks, as it's indeed implemented.
  int offset = 0;
  // internal_po keeps track of interations within a block.
  // Blocks have size MAX_RIGHT_DIM.
  int internal_po = 0;
  // External loop loops over blocks.
  while (true) {
    // If end of entire opration, break.
    if (po == total_dim - 1) break;

    internal_po = 0;
    // Each iteration of the while block goes through a new position.
    // Inside the while, j takes care of increasing indices properly.
    while (true) {
      po = 0;
      pn = 0;
      for (i = 0; i < num_indices; ++i) {
        po += old_super_dimensions[i] * old_counter[i];
        pn += new_super_dimensions[map_old_to_new_idxpos[i]] * old_counter[i];
      }
      small_map_old_to_new_position[po - offset] = pn;
      for (j = num_indices - 1; j >= 0; --j) {
        if (++old_counter[j] < old_dimensions[j])
          break;
        else
          old_counter[j] = 0;
      }
      // If end of block or end of entire operation, break.
      if ((++internal_po == MAX_RIGHT_DIM) || (po == total_dim - 1)) break;
      // If last index (0) was increased, then go back to fastest index.
      if (j < 0) break;
    }
    // Copy data for this block, taking into account offset of small_map...
    // The following line is to avoid casting MAX_RIGHT_DIM to size_t
    // every iteration. Note that it has to be size_t for min to work,
    // since total_dim is size_t.
    size_t effective_max = std::min((size_t)MAX_RIGHT_DIM, total_dim);
    for (size_t p = 0; p < effective_max; ++p)
      *(_data + small_map_old_to_new_position[p]) =
          *(scratch_copy + offset + p);

    offset += MAX_RIGHT_DIM;
  }

  _init(new_ordering, new_dimensions);

  scratch_copy = NULL;
}

void Tensor::_fast_reorder(std::vector<std::string> new_ordering,
                           s_type* scratch_copy) {
  if (scratch_copy == nullptr) {
    std::cout << "Scratch copy must be non-null." << std::endl;
    assert(scratch_copy != nullptr);
  }

  // Create binary orderings.
  std::vector<std::string> old_ordering(_indices);
  std::vector<size_t> old_dimensions(_dimensions);
  int num_indices = old_ordering.size();
  size_t total_dim = 1;
  for (int i = 0; i < num_indices; ++i) total_dim *= old_dimensions[i];
  // Create map_old_to_new_idxpos from old to new indices, and new_dimensions.
  std::vector<int> map_old_to_new_idxpos(num_indices);
  std::vector<size_t> new_dimensions(num_indices);
  for (int i = 0; i < num_indices; ++i) {
    for (int j = 0; j < num_indices; ++j) {
      if (old_ordering[i] == new_ordering[j]) {
        map_old_to_new_idxpos[i] = j;
        new_dimensions[j] = old_dimensions[i];
        break;
      }
    }
  }
  // Create binary orderings:
  std::vector<int> old_logs(num_indices);
  for (int i = 0; i < num_indices; ++i) {
    old_logs[i] = _LOG_2.at(old_dimensions[i]);
  }
  int num_binary_indices = _LOG_2.at(total_dim);
  // Create map from old letter to new group of letters.
  std::unordered_map<std::string, std::vector<std::string>> binary_groups;
  int alphabet_position = 0;
  for (int i = 0; i < num_indices; ++i) {
    std::vector<std::string> group(old_logs[i]);
    for (int j = 0; j < old_logs[i]; ++j) {
      group[j] = _ALPHABET[alphabet_position];
      ++alphabet_position;
    }
    binary_groups[old_ordering[i]] = group;
  }
  // Create old and new binary ordering in letters.
  std::vector<std::string> old_binary_ordering(num_binary_indices);
  std::vector<std::string> new_binary_ordering(num_binary_indices);
  int binary_position = 0;
  for (int i = 0; i < num_indices; ++i) {
    std::string old_index = old_ordering[i];
    for (int j = 0; j < binary_groups[old_index].size(); ++j) {
      old_binary_ordering[binary_position] = binary_groups[old_index][j];
      ++binary_position;
    }
  }
  binary_position = 0;
  for (int i = 0; i < num_indices; ++i) {
    std::string new_index = new_ordering[i];
    for (int j = 0; j < binary_groups[new_index].size(); ++j) {
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
  if (num_binary_indices <= _LOG_2.at(MAX_RIGHT_DIM)) {
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
    std::vector<std::string> Ll_old_indices(old_binary_ordering.begin(),
                                            old_binary_ordering.begin() + Ll);
    std::vector<std::string> Ll_new_indices(new_binary_ordering.begin(),
                                            new_binary_ordering.begin() + Ll);
    // Only one R10.
    if (Ll_old_indices == Ll_new_indices) {
      std::vector<std::string> Lr_old_indices(old_binary_ordering.begin() + Ll,
                                              old_binary_ordering.end());
      std::vector<std::string> Lr_new_indices(new_binary_ordering.begin() + Ll,
                                              new_binary_ordering.end());
      _right_reorder(Lr_old_indices, Lr_new_indices, Lr);
      scratch_copy = NULL;
      return;
    }
    // Only one L\nu move.
    for (int i = 5; i >= -1; --i) {
      int extended_Rr = Rr + i;
      std::vector<std::string> Rr_old_indices(
          old_binary_ordering.end() - extended_Rr, old_binary_ordering.end());
      std::vector<std::string> Rr_new_indices(
          new_binary_ordering.end() - extended_Rr, new_binary_ordering.end());
      if (Rr_old_indices == Rr_new_indices) {
        std::vector<std::string> Rl_old_indices(
            old_binary_ordering.begin(),
            old_binary_ordering.end() - extended_Rr);
        std::vector<std::string> Rl_new_indices(
            new_binary_ordering.begin(),
            new_binary_ordering.end() - extended_Rr);
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
    std::vector<std::string> Lr_indices(Lr), Ll_indices(Ll), Rr_indices(Rr),
        Rl_indices(Rl);
    for (int i = 0; i < Rr; ++i) Rr_indices[i] = new_binary_ordering[i + Rl];
    for (int i = 0; i < Rl; ++i) Rl_indices[i] = old_binary_ordering[i];
    std::vector<std::string> Rr_new_in_Rl_old =
        _vector_intersection(Rl_indices, Rr_indices);
    std::vector<std::string> Rl_old_not_in_Rr_new =
        _vector_subtraction(Rl_indices, Rr_new_in_Rl_old);
    std::vector<std::string> Rl_first_step =
        _vector_concatenation(Rl_old_not_in_Rr_new, Rr_new_in_Rl_old);
    std::vector<std::string> Rl_zeroth_step(Rl);
    for (int i = 0; i < Rl; ++i) Rl_zeroth_step[i] = old_binary_ordering[i];
    _left_reorder(Rl_zeroth_step, Rl_first_step, Rr, scratch_copy);
    // Done with 1).
    // Let's go with 2).
    std::vector<std::string> Lr_first_step = _vector_concatenation(
        std::vector<std::string>(Rl_first_step.begin() + Ll,
                                 Rl_first_step.end()),
        std::vector<std::string>(old_binary_ordering.begin() + Rl,
                                 old_binary_ordering.end()));
    Rr_indices = std::vector<std::string>(new_binary_ordering.begin() + Rl,
                                          new_binary_ordering.end());
    std::vector<std::string> Lr_second_step =
        _vector_concatenation(_vector_subtraction(Lr_first_step, Rr_indices),
                              std::vector<std::string>(Rr_indices));
    _right_reorder(Lr_first_step, Lr_second_step, Lr);
    // Done with 2).
    // Let's go with 3).
    std::vector<std::string> Rl_second_step = _vector_concatenation(
        std::vector<std::string>(Rl_first_step.begin(),
                                 Rl_first_step.begin() + Ll),
        std::vector<std::string>(Lr_second_step.begin(),
                                 Lr_second_step.begin() + Lr - Rr));
    std::vector<std::string> Rl_thrid_step(new_binary_ordering.begin(),
                                           new_binary_ordering.begin() + Rl);
    _left_reorder(Rl_second_step, Rl_thrid_step, Rr, scratch_copy);
    // done with 3).

    scratch_copy = NULL;
  }
}

// Assuming all indices are binary for old_ordering and new_ordering.
// old_ordering and new_ordering refer to the right.
void Tensor::_right_reorder(const std::vector<std::string>& old_ordering,
                            const std::vector<std::string>& new_ordering,
                            int num_indices_right) {
  // Don't do anything if there is nothing to reorder.
  if (new_ordering == old_ordering) return;

  // Create dim, num_indices, map_old_to_new_idxpos from old to new indices,
  // old_dimensions, new_dimensions, and total_dim.
  int dim = 2;
  int num_indices = old_ordering.size();
  std::vector<int> map_old_to_new_idxpos(num_indices);
  std::vector<size_t> old_dimensions(num_indices, dim);
  std::vector<size_t> new_dimensions(num_indices, dim);
  size_t total_dim = 1;
  for (int i = 0; i < num_indices; ++i) total_dim *= old_dimensions[i];
  for (int i = 0; i < num_indices; ++i) {
    for (int j = 0; j < num_indices; ++j) {
      if (old_ordering[i] == new_ordering[j]) {
        map_old_to_new_idxpos[i] = j;
        new_dimensions[j] = old_dimensions[i];
        break;
      }
    }
  }

  // Create the map_old_to_new_position, or get a reference to it if it exists
  // on _REORDER_MAPS.
  std::string name =
      _reordering_to_string(map_old_to_new_idxpos, old_dimensions);
  if (_REORDER_MAPS.find(name) == _REORDER_MAPS.end()) {
    _REORDER_MAPS[name] = std::vector<int>(total_dim);
    _generate_binary_reordering_map(map_old_to_new_idxpos,
                                    _REORDER_MAPS.at(name));
  }
  const std::vector<int>& map_old_to_new_position = _REORDER_MAPS.at(name);

  // With the map_old_to_new_position, we are ready to reorder within
  // small chuncks.
  int dim_right = total_dim;
  int dim_left = size() / dim_right;  // Remember, it's all powers of 2, so OK.
#pragma omp parallel
  {
    // For some reason, allocating these spaces and using them is about 2
    // times faster than bringing a pointer to a scratch space and using
    // different chunks of it.
    s_type* temp_data = new s_type[dim_right];
#pragma omp for schedule(static)
    for (int pl = 0; pl < dim_left; ++pl) {
      int current_thread = omp_get_thread_num();
      int offset = pl * dim_right;
      for (int pr = 0; pr < dim_right; ++pr)
        *(temp_data + pr) = *(_data + offset + pr);
      for (int pr = 0; pr < dim_right; ++pr)
        *(_data + offset + map_old_to_new_position[pr]) = *(temp_data + pr);
    }
    delete[] temp_data;
    temp_data = NULL;
  }
}

// Assuming all indices are binary for old_ordering and new_ordering.
// old_ordering and new_ordering refer to the left.
void Tensor::_left_reorder(const std::vector<std::string>& old_ordering,
                           const std::vector<std::string>& new_ordering,
                           int num_indices_right, s_type* scratch_copy) {
  if (scratch_copy == nullptr) {
    std::cout << "Scratch copy must be non-null." << std::endl;
    assert(scratch_copy != nullptr);
  }

  // Don't do anything if there is nothing to reorder.
  if (new_ordering == old_ordering) return;

  // Create dim, num_indices, map_old_to_new_idxpos from old to new indices,
  // old_dimensions, new_dimensions, and total_dim.
  int dim = 2;
  int num_indices = old_ordering.size();
  std::vector<int> map_old_to_new_idxpos(num_indices);
  std::vector<size_t> old_dimensions(num_indices, dim);
  std::vector<size_t> new_dimensions(num_indices, dim);
  size_t total_dim = 1;
  for (int i = 0; i < num_indices; ++i) total_dim *= old_dimensions[i];
  for (int i = 0; i < num_indices; ++i) {
    for (int j = 0; j < num_indices; ++j) {
      if (old_ordering[i] == new_ordering[j]) {
        map_old_to_new_idxpos[i] = j;
        new_dimensions[j] = old_dimensions[i];
        break;
      }
    }
  }
  // total_dim is the total dimension of the left!

  // Create the map_old_to_new_position, or get a reference to it if it exists
  // on _REORDER_MAPS.
  std::string name =
      _reordering_to_string(map_old_to_new_idxpos, old_dimensions);
  if (_REORDER_MAPS.find(name) == _REORDER_MAPS.end()) {
    _REORDER_MAPS[name] = std::vector<int>(total_dim);
    _generate_binary_reordering_map(map_old_to_new_idxpos,
                                    _REORDER_MAPS.at(name));
  }
  const std::vector<int>& map_old_to_new_position = _REORDER_MAPS.at(name);

  // With the map_old_to_new_position, we are ready to move small chunks.
  int dim_left = total_dim;
  size_t tensor_dim = size();
  int dim_right = tensor_dim / dim_left;  // Remember, it's all powers of 2,
                                          // so OK.
// Copy.
#pragma omp parallel for schedule(static, MAX_RIGHT_DIM)
  for (size_t p = 0; p < tensor_dim; ++p) {
    *(scratch_copy + p) = *(_data + p);
  }
// Move back.
#pragma omp parallel
  {
#pragma omp for schedule(static)
    for (int pl = 0; pl < dim_left; ++pl) {
      int old_offset = pl * dim_right;
      int new_offset = map_old_to_new_position[pl] * dim_right;
      for (int pr = 0; pr < dim_right; ++pr) {
        *(_data + new_offset + pr) = *(scratch_copy + old_offset + pr);
      }
    }
  }
  scratch_copy = NULL;
}

void Tensor::reorder(std::vector<std::string> new_ordering,
                     s_type* scratch_copy) {
  if (scratch_copy == nullptr) {
    std::cout << "Scratch copy must be non-null." << std::endl;
    assert(scratch_copy != nullptr);
  }

  // Asserts.
  bool new_ordering_in_indices = _vector_s_in_vector_s(new_ordering, _indices);
  bool indices_in_new_ordering = _vector_s_in_vector_s(_indices, new_ordering);
  if (!new_ordering_in_indices || !indices_in_new_ordering) {
    std::cout << "new_ordering: " << _string_vector_to_string(new_ordering)
              << " must be a reordering of current indices: "
              << _string_vector_to_string(_indices) << "." << std::endl;
    assert(new_ordering_in_indices && indices_in_new_ordering);
  }
  bool fast = true;
  for (int i = 0; i < _dimensions.size(); ++i) {
    if (_LOG_2.find(_dimensions[i]) == _LOG_2.end()) {
      fast = false;
      break;
    }
  }
  if (fast) {
    _fast_reorder(new_ordering, scratch_copy);
  } else {
    _naive_reorder(new_ordering, scratch_copy);
  }
}

void Tensor::scalar_multiply(s_type scalar) {
  size_t this_size = size();
#pragma omp parallel for schedule(static, MAX_RIGHT_DIM)
  for (size_t p = 0; p < this_size; ++p) *(_data + p) = (*(_data + p)) * scalar;
}

double Tensor::tensor_norm() const {
  double total_norm(0.0);
  size_t this_size = size();
  for (size_t p = 0; p < this_size; ++p) {
    total_norm += norm(*(_data + p));
  }
  return total_norm;
}

size_t Tensor::num_zeros() const {
  size_t count(0);
  size_t this_size = size();
  s_type complex_0(0.);
  for (size_t p = 0; p < this_size; ++p) {
    if (*(_data + p) == complex_0) ++count;
  }
  return count;
}

void Tensor::print() const {
  std::string print_str("Tensor of rank ");
  print_str += std::to_string(_indices.size());
  print_str += ": ";
  for (int i = 0; i < _indices.size(); ++i) {
    print_str += _indices[i];
    print_str += " -> ";
    print_str += std::to_string(_dimensions[i]);
    if (i < _indices.size() - 1) print_str += ", ";
  }
  std::cout << print_str << std::endl;
}

void Tensor::print_data() const {
  for (size_t p = 0; p < size(); ++p) std::cout << *(_data + p) << " ";
  std::cout << std::endl;
}

/////////////////////////// EXTERNAL FUNCTIONS ////////////////////////////////

// use  if complexity < some value.
void _multiply_MM(const s_type* A_data, const s_type* B_data, s_type* C_data,
                  int m, int n, int k) {
  if (A_data == nullptr) {
    std::cout << "Data from Tensor A must be non-null." << std::endl;
    assert(A_data != nullptr);
  }
  if (B_data == nullptr) {
    std::cout << "Data from Tensor B must be non-null." << std::endl;
    assert(B_data != nullptr);
  }
  if (C_data == nullptr) {
    std::cout << "Data from Tensor C must be non-null." << std::endl;
    assert(C_data != nullptr);
  }
  s_type alpha = 1.0;
  s_type beta = 0.0;
  cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, &alpha,
              A_data, std::max(1, k), B_data, std::max(1, n), &beta, C_data,
              std::max(1, n));
}

void _multiply_Mv(const s_type* A_data, const s_type* B_data, s_type* C_data,
                  int m, int k) {
  if (A_data == nullptr) {
    std::cout << "Data from Tensor A must be non-null." << std::endl;
    assert(A_data != nullptr);
  }
  if (B_data == nullptr) {
    std::cout << "Data from Tensor B must be non-null." << std::endl;
    assert(B_data != nullptr);
  }
  if (C_data == nullptr) {
    std::cout << "Data from Tensor C must be non-null." << std::endl;
    assert(C_data != nullptr);
  }
  s_type alpha = 1.0;
  s_type beta = 0.0;
  cblas_cgemv(CblasRowMajor, CblasNoTrans, m, k, &alpha, A_data, std::max(1, k),
              B_data, 1, &beta, C_data, 1);
}

void _multiply_vM(const s_type* A_data, const s_type* B_data, s_type* C_data,
                  int n, int k) {
  if (A_data == nullptr) {
    std::cout << "Data from Tensor A must be non-null." << std::endl;
    assert(A_data != nullptr);
  }
  if (B_data == nullptr) {
    std::cout << "Data from Tensor B must be non-null." << std::endl;
    assert(B_data != nullptr);
  }
  if (C_data == nullptr) {
    std::cout << "Data from Tensor C must be non-null." << std::endl;
    assert(C_data != nullptr);
  }
  s_type alpha = 1.0;
  s_type beta = 0.0;
  cblas_cgemv(CblasRowMajor, CblasTrans, k, n, &alpha, A_data, std::max(1, n),
              B_data, 1, &beta, C_data, 1);
}

void _multiply_vv(const s_type* A_data, const s_type* B_data, s_type* C_data,
                  int k) {
  if (A_data == nullptr) {
    std::cout << "Data from Tensor A must be non-null." << std::endl;
    assert(A_data != nullptr);
  }
  if (B_data == nullptr) {
    std::cout << "Data from Tensor B must be non-null." << std::endl;
    assert(B_data != nullptr);
  }
  if (C_data == nullptr) {
    std::cout << "Data from Tensor C must be non-null." << std::endl;
    assert(C_data != nullptr);
  }
  cblas_cdotu_sub(k, A_data, 1, B_data, 1, C_data);
}

void multiply(Tensor& A, Tensor& B, Tensor& C, s_type* scratch_copy) {
  if (scratch_copy == nullptr) {
    std::cout << "Scratch copy must be non-null." << std::endl;
    assert(scratch_copy != nullptr);
  }

  if (A.data() == C.data()) {
    std::cout << "A and C cannot be the same tensor: ";
    C.print();
    assert(A.data() != C.data());
  }
  if (B.data() == C.data()) {
    std::cout << "B and C cannot be the same tensor: ";
    C.print();
    assert(B.data() != C.data());
  }

  std::chrono::high_resolution_clock::time_point t0, t1;
  std::chrono::duration<double> time_span;
  t0 = std::chrono::high_resolution_clock::now();
  // Define left_indices, left_dim, right_indices, right_dim, and
  // common_indices, common_dim. Also C_size.
  std::vector<std::string> left_indices =
      _vector_subtraction(A.get_indices(), B.get_indices());
  std::vector<std::string> right_indices =
      _vector_subtraction(B.get_indices(), A.get_indices());
  std::vector<std::string> common_indices =
      _vector_intersection(A.get_indices(), B.get_indices());
  int left_dim = 1, right_dim = 1, common_dim = 1;
  for (int i = 0; i < left_indices.size(); ++i) {
    left_dim *= A.get_index_to_dimension().at(left_indices[i]);
  }
  for (int i = 0; i < right_indices.size(); ++i) {
    right_dim *= B.get_index_to_dimension().at(right_indices[i]);
  }
  for (int i = 0; i < common_indices.size(); ++i) {
    int a_dim = A.get_index_to_dimension().at(common_indices[i]);
    if (a_dim != B.get_index_to_dimension().at(common_indices[i])) {
      std::cout
          << "Common indices must have matching dimensions, but at index: " << i
          << ", the dimensions are A: " << a_dim
          << ", B: " << B.get_index_to_dimension().at(common_indices[i]) << "."
          << std::endl;
      assert(a_dim == B.get_index_to_dimension().at(common_indices[i]));
    }
    common_dim *= a_dim;
  }
  t1 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
  // std::cout << "Time preparing variables: " << time_span.count() << "s\n";

  // Assert.
  if (left_dim * right_dim > C.capacity()) {
    std::cout << "C: " << C.capacity()
              << " doesn't have enough space for the product of A*B: "
              << left_dim * right_dim << "." << std::endl;
    assert((left_dim * right_dim <= C.capacity()));
  }
  // Reorder.
  t0 = std::chrono::high_resolution_clock::now();
  std::vector<std::string> A_new_ordering =
      _vector_union(left_indices, common_indices);
  A.reorder(A_new_ordering, scratch_copy);
  t1 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
  // std::cout << "R " << time_span.count() << " s\n";
  // std::cout << "Time reordering A: " << time_span.count() << "s\n";
  t0 = std::chrono::high_resolution_clock::now();
  std::vector<std::string> B_new_ordering =
      _vector_union(common_indices, right_indices);
  B.reorder(B_new_ordering, scratch_copy);
  t1 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
  // std::cout << "R " << time_span.count() << " s\n";
  // std::cout << "Time reordering B: " << time_span.count() << "s\n";

  // Multiply. Four cases: MxM, Mxv, vxM, vxv.
  t0 = std::chrono::high_resolution_clock::now();
  if (left_indices.size() > 0 && right_indices.size() > 0) {
    _multiply_MM(A.data(), B.data(), C.data(), left_dim, right_dim, common_dim);
  } else if (left_indices.size() > 0 && right_indices.size() == 0) {
    _multiply_Mv(A.data(), B.data(), C.data(), left_dim, common_dim);
  } else if (left_indices.size() == 0 && right_indices.size() > 0) {
    // Very import to switch A and B to use cgemv with transpose for this case.
    _multiply_vM(B.data(), A.data(), C.data(), right_dim, common_dim);
  } else if (left_indices.size() == 0 && right_indices.size() == 0) {
    _multiply_vv(A.data(), B.data(), C.data(), common_dim);
  }
  t1 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
  // std::cout << "M " << time_span.count() << " s\n";
  // std::cout << "Time multiplying A*B: " << time_span.count() << "s\n";

  // Set indices and dimensions of C.
  t0 = std::chrono::high_resolution_clock::now();
  std::vector<std::string> C_indices =
      _vector_union(left_indices, right_indices);
  std::vector<size_t> C_dimensions(C_indices.size());
  for (int i = 0; i < left_indices.size(); ++i)
    C_dimensions[i] = A.get_index_to_dimension().at(left_indices[i]);
  for (int i = 0; i < right_indices.size(); ++i)
    C_dimensions[i + left_indices.size()] =
        B.get_index_to_dimension().at(right_indices[i]);
  C.set_indices(C_indices);
  C.set_dimensions(C_dimensions);
  C.generate_index_to_dimension();
  t1 = std::chrono::high_resolution_clock::now();
  time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
  // std::cout << "Time updating C's variables: " << time_span.count() << "s\n";
}

size_t result_size(Tensor& A, Tensor& B) {
  std::vector<std::string> left_indices =
      _vector_subtraction(A.get_indices(), B.get_indices());
  std::vector<std::string> right_indices =
      _vector_subtraction(B.get_indices(), A.get_indices());
  size_t left_dim = 1, right_dim = 1, result_dim;
  for (int i = 0; i < left_indices.size(); ++i) {
    left_dim *= A.get_index_to_dimension().at(left_indices[i]);
  }
  for (int i = 0; i < right_indices.size(); ++i) {
    right_dim *= B.get_index_to_dimension().at(right_indices[i]);
  }
  result_dim = left_dim * right_dim;
  return result_dim;
}

// Split it in parts, as before? Nah, it was all about generating small maps.
// Here I am generating THE map. This is done only once per map anyway.
void _generate_binary_reordering_map(
    const std::vector<int>& map_old_to_new_idxpos,
    std::vector<int>& map_old_to_new_position) {
  int dim = 2;  // Hard coded!
  int num_indices = map_old_to_new_idxpos.size();
  // Assert
  if ((size_t)pow(dim, num_indices) != map_old_to_new_position.size()) {
    std::cout << "Size of map: " << map_old_to_new_position.size()
              << " must be equal to 2^num_indices: " << pow(dim, num_indices)
              << "." << std::endl;
    assert((size_t)pow(dim, num_indices) == map_old_to_new_position.size());
  }

  // Define super dimensions. See _naive_reorder().
  std::vector<int> old_dimensions(num_indices, dim);
  std::vector<int> new_dimensions(num_indices, dim);
  std::vector<int> old_super_dimensions(num_indices);
  std::vector<int> new_super_dimensions(num_indices);
  old_super_dimensions[num_indices - 1] = 1;
  new_super_dimensions[num_indices - 1] = 1;
  for (int i = num_indices - 2; i >= 0; --i) {
    old_super_dimensions[i] = old_super_dimensions[i + 1] * dim;
    new_super_dimensions[i] = new_super_dimensions[i + 1] * dim;
  }

  // Iterate and generate map.
  std::vector<int> old_counter(num_indices, 0);
  size_t po, pn;  // Position of the data, old and new.
  int i, j;
  while (true) {
    po = 0;
    pn = 0;
    for (i = 0; i < num_indices; ++i) {
      po += old_super_dimensions[i] * old_counter[i];
      pn += new_super_dimensions[map_old_to_new_idxpos[i]] * old_counter[i];
    }
    map_old_to_new_position[po] = pn;
    for (j = num_indices - 1; j >= 0; --j) {
      if (++old_counter[j] < old_dimensions[j])
        break;
      else
        old_counter[j] = 0;
    }
    if (j < 0) break;
  }
}

// convert int vector to string
std::string _int_vector_to_string(std::vector<int> input) {
  std::ostringstream temp;
  std::string output;
  if (!input.empty()) {
    std::copy(input.begin(), input.end() - 1,
              std::ostream_iterator<int>(temp, ", "));
    temp << input.back();
  }
  output = "{" + temp.str() + "}";
  return output;
}

// convert string vector to string
std::string _string_vector_to_string(std::vector<std::string> input) {
  std::string output;
  output += "{";
  if (!input.empty()) {
    for (std::vector<std::string>::const_iterator i = input.begin();
         i < input.end(); ++i) {
      output += *i;
      if (i != input.end() - 1) {
        output += ", ";
      }
    }
  }
  output += "}";
  return output;
}

std::string _reordering_to_string(const std::vector<int>& map_old_to_new_idxpos,
                                  const std::vector<size_t>& old_dimensions) {
  int num_indices = map_old_to_new_idxpos.size();
  std::string name("");
  for (int i = 0; i < num_indices; ++i) name += _ALPHABET[i];
  name += "->";
  for (int i = 0; i < num_indices; ++i)
    name += _ALPHABET[map_old_to_new_idxpos[i]];
  for (int i = 0; i < num_indices; ++i)
    name += "," + std::to_string(old_dimensions[i]);
  return name;
}

bool _string_in_vector(const std::string& s,
                       const std::vector<std::string>& v) {
  if (std::find(v.cbegin(), v.cend(), s) != v.cend())
    return true;
  else
    return false;
}

bool _vector_s_in_vector_s(const std::vector<std::string>& v,
                           const std::vector<std::string>& w) {
  for (int i = 0; i < v.size(); ++i) {
    if (!_string_in_vector(v[i], w)) return false;
  }
  return true;
}

std::vector<std::string> _vector_intersection(
    const std::vector<std::string>& v, const std::vector<std::string>& w) {
  std::vector<std::string> temp;
  for (auto it = v.cbegin(); it != v.cend(); ++it) {
    if (_string_in_vector(*it, w)) temp.push_back(*it);
  }
  return temp;
}

std::vector<std::string> _vector_union(const std::vector<std::string>& v,
                                       const std::vector<std::string>& w) {
  std::vector<std::string> temp(v);
  for (auto it = w.cbegin(); it != w.cend(); ++it) {
    if (!_string_in_vector(*it, v)) temp.push_back(*it);
  }
  return temp;
}

std::vector<std::string> _vector_subtraction(
    const std::vector<std::string>& v, const std::vector<std::string>& w) {
  std::vector<std::string> temp;
  for (auto it = v.cbegin(); it != v.cend(); ++it) {
    if (!_string_in_vector(*it, w)) temp.push_back(*it);
  }
  return temp;
}

std::vector<std::string> _vector_concatenation(
    const std::vector<std::string>& v, const std::vector<std::string>& w) {
  std::vector<std::string> temp(v.size() + w.size());
  for (size_t i = 0; i < v.size(); ++i) temp[i] = v[i];
  for (size_t i = 0; i < w.size(); ++i) temp[i + v.size()] = w[i];
  return temp;
}

}  // namespace qflex
