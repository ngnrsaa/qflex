/**
 * @file tensor.cpp
 * Implementation of the Tensor class.
 *
 * @author Benjamin Villalonga (main contributor), Bron Nelson, Sergio Boixo and
 * Salvatore Mandra
 * @contributors: The qFlex Developers (see CONTRIBUTORS.md)
 * @date Created: August 2018
 *
 * @copyright: Copyright © 2019, United States Government, as represented
 * by the Administrator of the National Aeronautics and Space Administration.
 * All rights reserved.
 * @licence: Apache License, Version 2.0
 */

#include "tensor.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
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
 * unordered_map<std::size_t,std::size_t> with the log2 of powers of 2 up to 2^30, in order to
 * quickly switch to smart reordering and look up the logs.
 */
const std::unordered_map<std::size_t, std::size_t> _LOG_2(
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
 * Global unordered_map<string,vector<std::size_t>> of reordering maps.
 */
std::unordered_map<std::string, std::vector<std::size_t>> _REORDER_MAPS;

///////////////////////////// CLASS FUNCTIONS /////////////////////////////////

void Tensor::_init(const std::vector<std::string>& indices,
                   const std::vector<std::size_t>& dimensions) {
  if (indices.size() != dimensions.size()) {
    throw ERROR_MSG("The number of indices: ", indices.size(),
                    ", and number of dimensions: ", dimensions.size(),
                    ", should be equal.");
  }
  _indices = indices;
  _dimensions = dimensions;
  for (std::size_t i = 0; i < _dimensions.size(); ++i)
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
    // The line "set_dimensions(other.get_dimensions());" takes care of the
    // total size of the dimensions.
    try {
      set_dimensions(other.get_dimensions());
    } catch (const std::string& err_msg) {
      throw ERROR_MSG("Failed to call set_dimensions(). Error:\n\t[", err_msg,
                      "]");
    }
  }
  try {
    _init(other.get_indices(), other.get_dimensions());
  } catch (const std::string& err_msg) {
    throw ERROR_MSG("Failed to call _init(). Error:\n\t[", err_msg, "]");
  }

#pragma omp parallel for schedule(static, MAX_RIGHT_DIM)
  for (std::size_t p = 0; p < other.size(); ++p)
    *(_data + p) = *(other.data() + p);
}

Tensor::Tensor() { _data = NULL; }

Tensor::Tensor(std::vector<std::string> indices,
               std::vector<std::size_t> dimensions) {
  try {
    _init(indices, dimensions);
  } catch (const std::string& err_msg) {
    throw ERROR_MSG("Failed to call _init(). Error:\n\t[", err_msg, "]");
  }
  _capacity = size();
  _data = new s_type[_capacity];
}

Tensor::Tensor(std::vector<std::string> indices,
               std::vector<std::size_t> dimensions,
               const std::vector<s_type>& data)
    : Tensor(indices, dimensions) {
  // Check that the data has the same length as this Tensor's size().
  std::size_t this_size = size();
  if (this_size != data.size()) {
    throw ERROR_MSG("The vector data size: ", data.size(),
                    ", has to match the size of the Tensor: ", this_size);
  }
  _capacity = this_size;
  // Fill in the _data.
  for (std::size_t i = 0; i < this_size; ++i) *(_data + i) = data[i];
}

Tensor::Tensor(std::vector<std::string> indices,
               std::vector<std::size_t> dimensions, s_type* data) {
  if (data == nullptr) {
    throw ERROR_MSG("Data must be non-null.");
  }
  try {
    _init(indices, dimensions);
  } catch (const std::string& err_msg) {
    throw ERROR_MSG("Failed to call _init(). Error:\n\t[", err_msg, "]");
  }
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

const std::vector<std::size_t>& Tensor::get_dimensions() const {
  return _dimensions;
}

void Tensor::set_dimensions(const std::vector<std::size_t>& dimensions) {
  // Exception.
  if (_data) {
    std::size_t total_dim = 1;
    for (std::size_t i = 0; i < dimensions.size(); ++i)
      total_dim *= dimensions[i];
    if (capacity() < total_dim) {
      throw ERROR_MSG("The total allocated space: ", capacity(),
                      ", is insufficient for the requested tensor dimensions: ",
                      total_dim, ".");
    }
  }
  _dimensions = dimensions;
}

void Tensor::set_indices_and_dimensions(
    const std::vector<std::string>& indices,
    const std::vector<std::size_t>& dimensions) {
  // The following line takes care of the total size of the dimensions.
  try {
    set_dimensions(dimensions);
  } catch (const std::string& err_msg) {
    throw ERROR_MSG("Failed to call set_dimensions(). Error:\n\t[", err_msg,
                    "]");
  }
  try {
    _init(indices, dimensions);
  } catch (const std::string& err_msg) {
    throw ERROR_MSG("Failed to call _init(). Error:\n\t[", err_msg, "]");
  }
}

const std::unordered_map<std::string, std::size_t>&
Tensor::get_index_to_dimension() const {
  return _index_to_dimension;
}

void Tensor::generate_index_to_dimension() {
  for (std::size_t i = 0; i < _indices.size(); ++i)
    _index_to_dimension[_indices[i]] = _dimensions[i];
}

std::size_t Tensor::size() const {
  std::size_t total_dim = 1;
  for (std::size_t i = 0; i < _dimensions.size(); ++i)
    total_dim *= _dimensions[i];
  return total_dim;
}

std::size_t Tensor::capacity() const { return _capacity; }

s_type* Tensor::data() { return _data; }

const s_type* Tensor::data() const { return _data; }

void Tensor::project(std::string index, std::size_t index_value,
                     Tensor& projection_tensor) const {
  if (index != _indices[0]) {
    throw ERROR_MSG("Index: '", index, "' has to be equal to indices[0]: '",
                    _indices[0], "'.");
  }
  if (index_value >= _dimensions[0]) {
    throw ERROR_MSG("index_value: ", index_value, " must be contained in [0, ",
                    _dimensions[0], ").");
  }
  // Resize projection_tensor first.
  std::vector<std::string> projection_indices(_indices.begin() + 1,
                                              _indices.end());
  std::vector<std::size_t> projection_dimensions(_dimensions.begin() + 1,
                                                 _dimensions.end());
  projection_tensor.set_indices(projection_indices);
  try {
    projection_tensor.set_dimensions(projection_dimensions);
  } catch (const std::string& err_msg) {
    throw ERROR_MSG("Failed to call set_dimensions(). Error:\n\t[", err_msg,
                    "]");
  }
  projection_tensor.generate_index_to_dimension();

  // Fill projection_tensor with result of projection.
  s_type* projection_data = projection_tensor.data();
  std::size_t projection_size = projection_tensor.size();
  std::size_t projection_begin = projection_size * index_value;
#pragma omp parallel for schedule(static, MAX_RIGHT_DIM)
  for (std::size_t p = 0; p < projection_size; ++p)
    *(projection_data + p) = *(_data + projection_begin + p);
}

void Tensor::rename_index(std::string old_name, std::string new_name) {
  auto it = find(_indices.begin(), _indices.end(), old_name);
  if (it == _indices.end()) {
    throw ERROR_MSG("old_name: ", old_name, ", has to be a valid index.");
  }
  bool new_name_is_existing_index =
      (find(_indices.begin(), _indices.end(), new_name) != _indices.end());
  if (new_name_is_existing_index) {
    throw ERROR_MSG("new_name: ", new_name, ", cannot be an existing index.");
  }
  *it = new_name;
  _index_to_dimension[new_name] = _index_to_dimension[old_name];
  _index_to_dimension.erase(old_name);
}

void Tensor::bundle(std::vector<std::string> indices_to_bundle,
                    std::string bundled_index) {
  // Checks.
  bool indices_to_bundle_in_indices =
      _vector_s_in_vector_s(indices_to_bundle, _indices);
  if (!indices_to_bundle_in_indices) {
    throw ERROR_MSG(
        "indices_to_bundle: ", _string_vector_to_string(indices_to_bundle),
        " has to be contained in indices: ", _string_vector_to_string(_indices),
        ".");
  }
  std::vector<std::string> subtracted_indices(
      _vector_subtraction(_indices, indices_to_bundle));
  std::vector<std::string> indices_to_bundled_original_order(
      _vector_subtraction(_indices, subtracted_indices));
  if (indices_to_bundled_original_order != indices_to_bundle) {
    throw ERROR_MSG(
        "indices_to_bundle: ", _string_vector_to_string(indices_to_bundle),
        " must be in its original order: ",
        _string_vector_to_string(indices_to_bundled_original_order), ".");
  }

  std::size_t bundled_dim = 1;
  for (std::size_t i = 0; i < indices_to_bundle.size(); ++i) {
    bundled_dim *= _index_to_dimension[indices_to_bundle[i]];
    _index_to_dimension.erase(indices_to_bundle[i]);
  }
  _index_to_dimension[bundled_index] = bundled_dim;
  std::size_t bundled_idxpos = 0;
  for (std::size_t i = 0; i < _indices.size(); ++i) {
    if (_string_in_vector(_indices[i], indices_to_bundle)) {
      bundled_idxpos = i;
      break;
    }
  }
  std::vector<std::string> new_indices(subtracted_indices);
  new_indices.insert(new_indices.begin() + bundled_idxpos, bundled_index);
  std::vector<std::size_t> new_dimensions(new_indices.size());
  for (std::size_t i = 0; i < new_dimensions.size(); ++i)
    new_dimensions[i] = _index_to_dimension[new_indices[i]];
  _indices = new_indices;
  _dimensions = new_dimensions;
}

void Tensor::_naive_reorder(std::vector<std::string> new_ordering,
                            s_type* scratch_copy) {
  if (scratch_copy == nullptr) {
    throw ERROR_MSG("Scratch copy must be non-null.");
  }

  // Don't do anything if there is nothing to reorder.
  if (new_ordering == _indices) return;

  std::vector<std::string> old_ordering(_indices);
  std::vector<std::size_t> old_dimensions(_dimensions);
  std::size_t num_indices = old_ordering.size();
  std::size_t total_dim = size();

  if (num_indices == 0) throw ERROR_MSG("Number of indices cannot be zero.");

  // Create map_old_to_new_idxpos from old to new indices, and new_dimensions.
  std::vector<std::size_t> map_old_to_new_idxpos(num_indices);
  std::vector<std::size_t> new_dimensions(num_indices);
  for (std::size_t i = 0; i < num_indices; ++i) {
    for (std::size_t j = 0; j < num_indices; ++j) {
      if (old_ordering[i] == new_ordering[j]) {
        map_old_to_new_idxpos[i] = j;
        new_dimensions[j] = old_dimensions[i];
        break;
      }
    }
  }

  // Create super dimensions (combined dimension of all to the right of i).
  std::vector<std::size_t> old_super_dimensions(num_indices);
  std::vector<std::size_t> new_super_dimensions(num_indices);
  old_super_dimensions[num_indices - 1] = 1;
  new_super_dimensions[num_indices - 1] = 1;
  if (std::size_t size = old_dimensions.size(); size >= 2)
    for (std::size_t i = size; --i;) {
      old_super_dimensions[i - 1] = old_super_dimensions[i] * old_dimensions[i];
      new_super_dimensions[i - 1] = new_super_dimensions[i] * new_dimensions[i];
    }

  // Allocating small_map_old_to_new_position.
  std::vector<unsigned short int> small_map_old_to_new_position(MAX_RIGHT_DIM);

// Start moving data around.
// First copy all data into scratch.
#pragma omp parallel for schedule(static, MAX_RIGHT_DIM)
  for (std::size_t p = 0; p < total_dim; ++p)
    *(scratch_copy + p) = *(_data + p);

  // No combined efficient mapping from old to new positions with actual
  // copies in memory, all in small cache friendly (for old data, not new,
  // which could be very scattered) blocks.

  // Position old and new.
  std::size_t po = 0, pn;
  // Counter of the values of each indices in the iteration (old ordering).
  std::vector<std::size_t> old_counter(num_indices, 0);
  // offset is important when doing this in blocks, as it's indeed implemented.
  std::size_t offset = 0;
  // internal_po keeps track of interations within a block.
  // Blocks have size MAX_RIGHT_DIM.
  std::size_t internal_po = 0;

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
      for (std::size_t i = 0; i < num_indices; ++i) {
        po += old_super_dimensions[i] * old_counter[i];
        pn += new_super_dimensions[map_old_to_new_idxpos[i]] * old_counter[i];
      }
      small_map_old_to_new_position[po - offset] = pn;

      bool complete{true};
      for (std::size_t j = num_indices; j--;) {
        if (++old_counter[j] < old_dimensions[j]) {
          complete = false;
          break;
        } else
          old_counter[j] = 0;
      }
      // If end of block or end of entire operation, break.
      if ((++internal_po == MAX_RIGHT_DIM) || (po == total_dim - 1)) break;
      // If last index (0) was increased, then go back to fastest index.
      if (complete) break;
    }
    // Copy data for this block, taking into account offset of small_map...
    // The following line is to avoid casting MAX_RIGHT_DIM to std::size_t
    // every iteration. Note that it has to be std::size_t for min to work,
    // since total_dim is std::size_t.
    std::size_t effective_max = std::min((std::size_t)MAX_RIGHT_DIM, total_dim);
    for (std::size_t p = 0; p < effective_max; ++p)
      *(_data + small_map_old_to_new_position[p]) =
          *(scratch_copy + offset + p);

    offset += MAX_RIGHT_DIM;
  }
  try {
    _init(new_ordering, new_dimensions);
  } catch (const std::string& err_msg) {
    throw ERROR_MSG("Failed to call _init(). Error:\n\t[", err_msg, "]");
  }

  scratch_copy = NULL;
}

void Tensor::_fast_reorder(std::vector<std::string> new_ordering,
                           s_type* scratch_copy) {
  if (scratch_copy == nullptr) {
    throw ERROR_MSG("Scratch copy must be non-null.");
  }

  // Create binary orderings.
  std::vector<std::string> old_ordering(_indices);
  std::vector<std::size_t> old_dimensions(_dimensions);
  std::size_t num_indices = old_ordering.size();
  std::size_t total_dim = 1;
  for (std::size_t i = 0; i < num_indices; ++i) total_dim *= old_dimensions[i];
  // Create map_old_to_new_idxpos from old to new indices, and new_dimensions.
  std::vector<std::size_t> map_old_to_new_idxpos(num_indices);
  std::vector<std::size_t> new_dimensions(num_indices);
  for (std::size_t i = 0; i < num_indices; ++i) {
    for (std::size_t j = 0; j < num_indices; ++j) {
      if (old_ordering[i] == new_ordering[j]) {
        map_old_to_new_idxpos[i] = j;
        new_dimensions[j] = old_dimensions[i];
        break;
      }
    }
  }
  // Create binary orderings:
  std::vector<std::size_t> old_logs(num_indices);
  for (std::size_t i = 0; i < num_indices; ++i) {
    old_logs[i] = _LOG_2.at(old_dimensions[i]);
  }
  std::size_t num_binary_indices = _LOG_2.at(total_dim);
  // Create map from old letter to new group of letters.
  std::unordered_map<std::string, std::vector<std::string>> binary_groups;
  std::size_t alphabet_position = 0;
  for (std::size_t i = 0; i < num_indices; ++i) {
    std::vector<std::string> group(old_logs[i]);
    for (std::size_t j = 0; j < old_logs[i]; ++j) {
      group[j] = _ALPHABET[alphabet_position];
      ++alphabet_position;
    }
    binary_groups[old_ordering[i]] = group;
  }
  // Create old and new binary ordering in letters.
  std::vector<std::string> old_binary_ordering(num_binary_indices);
  std::vector<std::string> new_binary_ordering(num_binary_indices);
  std::size_t binary_position = 0;
  for (std::size_t i = 0; i < num_indices; ++i) {
    std::string old_index = old_ordering[i];
    for (std::size_t j = 0; j < binary_groups[old_index].size(); ++j) {
      old_binary_ordering[binary_position] = binary_groups[old_index][j];
      ++binary_position;
    }
  }
  binary_position = 0;
  for (std::size_t i = 0; i < num_indices; ++i) {
    std::string new_index = new_ordering[i];
    for (std::size_t j = 0; j < binary_groups[new_index].size(); ++j) {
      new_binary_ordering[binary_position] = binary_groups[new_index][j];
      ++binary_position;
    }
  }
  // Up to here, I have created old_binary_ordering and new_binary_ordering.

  // Change _indices and _dimensions, as well as _index_to_dimension.
  // This is common to all cases, special or default (worst case).
  try {
    _init(new_ordering, new_dimensions);
  } catch (const std::string& err_msg) {
    throw ERROR_MSG("Failed to call _init(). Error:\n\t[", err_msg, "]");
  }

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
    if (new_binary_ordering.size() < _LOG_2.at(MAX_RIGHT_DIM))
      throw ERROR_MSG("New ordering is too small to be used at this point.");

    std::size_t Lr = _LOG_2.at(MAX_RIGHT_DIM);
    std::size_t Ll = new_binary_ordering.size() - Lr;
    std::size_t Rr = _LOG_2.at(MIN_RIGHT_DIM);
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

    if (Rr == 0) throw ERROR_MSG("Rr move cannot be zero.");

    // TODO: This loop has been tested to make sure that extended_Rr is
    // consistent with its previous implementation. However, no simulations so
    // far seem to use this loop, so I cannot check it!
    //
    // Only one L\nu move.
    // for (long int i = 5; i >= -1; --i) {
    //  long int extended_Rr = Rr + i;
    for (std::size_t i = 7; i--;) {
      std::size_t extended_Rr = Rr + i - 1;
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
        try {
          _left_reorder(Rl_old_indices, Rl_new_indices, extended_Rr,
                        scratch_copy);
        } catch (const std::string& err_msg) {
          throw ERROR_MSG("Failed to call _left_reorder(). Error:\n\t[",
                          err_msg, "]");
        }
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

    if (new_binary_ordering.size() < _LOG_2.at(MAX_RIGHT_DIM))
      throw ERROR_MSG("New ordering is too small to be used at this point.");
    if (new_binary_ordering.size() < _LOG_2.at(MIN_RIGHT_DIM))
      throw ERROR_MSG("New ordering is too small to be used at this point.");

    std::size_t Lr = _LOG_2.at(MAX_RIGHT_DIM);
    std::size_t Ll = new_binary_ordering.size() - Lr;
    std::size_t Rr = _LOG_2.at(MIN_RIGHT_DIM);
    std::size_t Rl = new_binary_ordering.size() - Rr;
    // Helper vectors that can be reused.
    std::vector<std::string> Lr_indices(Lr), Ll_indices(Ll), Rr_indices(Rr),
        Rl_indices(Rl);
    for (std::size_t i = 0; i < Rr; ++i)
      Rr_indices[i] = new_binary_ordering[i + Rl];
    for (std::size_t i = 0; i < Rl; ++i) Rl_indices[i] = old_binary_ordering[i];
    std::vector<std::string> Rr_new_in_Rl_old =
        _vector_intersection(Rl_indices, Rr_indices);
    std::vector<std::string> Rl_old_not_in_Rr_new =
        _vector_subtraction(Rl_indices, Rr_new_in_Rl_old);
    std::vector<std::string> Rl_first_step =
        _vector_concatenation(Rl_old_not_in_Rr_new, Rr_new_in_Rl_old);
    std::vector<std::string> Rl_zeroth_step(Rl);
    for (std::size_t i = 0; i < Rl; ++i)
      Rl_zeroth_step[i] = old_binary_ordering[i];
    try {
      _left_reorder(Rl_zeroth_step, Rl_first_step, Rr, scratch_copy);
    } catch (const std::string& err_msg) {
      throw ERROR_MSG("Failed to call _left_reorder(). Error:\n\t[", err_msg,
                      "]");
    }
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
    try {
      _left_reorder(Rl_second_step, Rl_thrid_step, Rr, scratch_copy);
    } catch (const std::string& err_msg) {
      throw ERROR_MSG("Failed to call _left_reorder(). Error:\n\t[", err_msg,
                      "]");
    }
    // done with 3).

    scratch_copy = NULL;
  }
}

// Assuming all indices are binary for old_ordering and new_ordering.
// old_ordering and new_ordering refer to the right.
void Tensor::_right_reorder(const std::vector<std::string>& old_ordering,
                            const std::vector<std::string>& new_ordering,
                            std::size_t num_indices_right) {
  // Don't do anything if there is nothing to reorder.
  if (new_ordering == old_ordering) return;

  // Create dim, num_indices, map_old_to_new_idxpos from old to new indices,
  // old_dimensions, new_dimensions, and total_dim.
  std::size_t dim = 2;
  std::size_t num_indices = old_ordering.size();
  std::vector<std::size_t> map_old_to_new_idxpos(num_indices);
  std::vector<std::size_t> old_dimensions(num_indices, dim);
  std::vector<std::size_t> new_dimensions(num_indices, dim);
  std::size_t total_dim = 1;
  for (std::size_t i = 0; i < num_indices; ++i) total_dim *= old_dimensions[i];
  for (std::size_t i = 0; i < num_indices; ++i) {
    for (std::size_t j = 0; j < num_indices; ++j) {
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
    _REORDER_MAPS[name] = std::vector<std::size_t>(total_dim);
    try {
      _generate_binary_reordering_map(map_old_to_new_idxpos,
                                      _REORDER_MAPS.at(name));
    } catch (const std::string& err_msg) {
      throw ERROR_MSG(
          "Failed to call _generate_binary_reordering_map(). Error:\n\t[",
          err_msg, "]");
    }
  }
  const std::vector<std::size_t>& map_old_to_new_position =
      _REORDER_MAPS.at(name);

  // With the map_old_to_new_position, we are ready to reorder within
  // small chuncks.
  std::size_t dim_right = total_dim;
  std::size_t dim_left =
      size() / dim_right;  // Remember, it's all powers of 2, so OK.
#pragma omp parallel
  {
    // For some reason, allocating these spaces and using them is about 2
    // times faster than bringing a pointer to a scratch space and using
    // different chunks of it.
    s_type* temp_data = new s_type[dim_right];
#pragma omp for schedule(static)
    for (std::size_t pl = 0; pl < dim_left; ++pl) {
#ifdef _OPENMP
      std::size_t current_thread = omp_get_thread_num();
#endif
      std::size_t offset = pl * dim_right;
      for (std::size_t pr = 0; pr < dim_right; ++pr)
        *(temp_data + pr) = *(_data + offset + pr);
      for (std::size_t pr = 0; pr < dim_right; ++pr)
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
                           std::size_t num_indices_right,
                           s_type* scratch_copy) {
  if (scratch_copy == nullptr) {
    throw ERROR_MSG("Scratch copy must be non-null.");
  }

  // Don't do anything if there is nothing to reorder.
  if (new_ordering == old_ordering) return;

  // Create dim, num_indices, map_old_to_new_idxpos from old to new indices,
  // old_dimensions, new_dimensions, and total_dim.
  std::size_t dim = 2;
  std::size_t num_indices = old_ordering.size();
  std::vector<std::size_t> map_old_to_new_idxpos(num_indices);
  std::vector<std::size_t> old_dimensions(num_indices, dim);
  std::vector<std::size_t> new_dimensions(num_indices, dim);
  std::size_t total_dim = 1;
  for (std::size_t i = 0; i < num_indices; ++i) total_dim *= old_dimensions[i];
  for (std::size_t i = 0; i < num_indices; ++i) {
    for (std::size_t j = 0; j < num_indices; ++j) {
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
    _REORDER_MAPS[name] = std::vector<std::size_t>(total_dim);
    try {
      _generate_binary_reordering_map(map_old_to_new_idxpos,
                                      _REORDER_MAPS.at(name));
    } catch (const std::string& err_msg) {
      throw ERROR_MSG(
          "Failed to call _generate_binary_reordering_map(). Error:\n\t[",
          err_msg, "]");
    }
  }
  const std::vector<std::size_t>& map_old_to_new_position =
      _REORDER_MAPS.at(name);

  // With the map_old_to_new_position, we are ready to move small chunks.
  std::size_t dim_left = total_dim;
  std::size_t tensor_dim = size();
  std::size_t dim_right = tensor_dim / dim_left;  // Remember, it's all powers
                                                  // of 2, so OK.
// Copy.
#pragma omp parallel for schedule(static, MAX_RIGHT_DIM)
  for (std::size_t p = 0; p < tensor_dim; ++p) {
    *(scratch_copy + p) = *(_data + p);
  }
// Move back.
#pragma omp parallel
  {
#pragma omp for schedule(static)
    for (std::size_t pl = 0; pl < dim_left; ++pl) {
      std::size_t old_offset = pl * dim_right;
      std::size_t new_offset = map_old_to_new_position[pl] * dim_right;
      for (std::size_t pr = 0; pr < dim_right; ++pr) {
        *(_data + new_offset + pr) = *(scratch_copy + old_offset + pr);
      }
    }
  }
  scratch_copy = NULL;
}

void Tensor::reorder(std::vector<std::string> new_ordering,
                     s_type* scratch_copy) {
  if (scratch_copy == nullptr) {
    throw ERROR_MSG("Scratch copy must be non-null.");
  }

  // Checks.
  bool new_ordering_in_indices = _vector_s_in_vector_s(new_ordering, _indices);
  bool indices_in_new_ordering = _vector_s_in_vector_s(_indices, new_ordering);
  if (!new_ordering_in_indices || !indices_in_new_ordering) {
    throw ERROR_MSG("new_ordering: ", _string_vector_to_string(new_ordering),
                    " must be a reordering of current indices: ",
                    _string_vector_to_string(_indices), ".");
  }
  bool fast = true;
  for (std::size_t i = 0; i < _dimensions.size(); ++i) {
    if (_LOG_2.find(_dimensions[i]) == _LOG_2.end()) {
      fast = false;
      break;
    }
  }
  if (fast) {
    try {
      _fast_reorder(new_ordering, scratch_copy);
    } catch (const std::string& err_msg) {
      throw ERROR_MSG("Failed to call _fast_reorder(). Error:\n\t[", err_msg,
                      "]");
    }
  } else {
    try {
      _naive_reorder(new_ordering, scratch_copy);
    } catch (const std::string& err_msg) {
      throw ERROR_MSG("Failed to call _naive_reorder(). Error:\n\t[", err_msg,
                      "]");
    }
  }
}

void Tensor::scalar_multiply(s_type scalar) {
  std::size_t this_size = size();
#pragma omp parallel for schedule(static, MAX_RIGHT_DIM)
  for (std::size_t p = 0; p < this_size; ++p)
    *(_data + p) = (*(_data + p)) * scalar;
}

double Tensor::tensor_norm() const {
  double total_norm(0.0);
  std::size_t this_size = size();
  for (std::size_t p = 0; p < this_size; ++p) {
    total_norm += norm(*(_data + p));
  }
  return total_norm;
}

std::size_t Tensor::num_zeros() const {
  std::size_t count(0);
  std::size_t this_size = size();
  s_type complex_0(0.);
  for (std::size_t p = 0; p < this_size; ++p) {
    if (*(_data + p) == complex_0) ++count;
  }
  return count;
}

std::string Tensor::tensor_to_string() const {
  std::string output_str("Tensor of rank ");
  output_str += std::to_string(_indices.size());
  output_str += ": ";
  for (std::size_t i = 0; i < _indices.size(); ++i) {
    output_str += _indices[i];
    output_str += " -> ";
    output_str += std::to_string(_dimensions[i]);
    if (i < _indices.size() - 1) output_str += ", ";
  }
  return output_str;
}

void Tensor::print_data() const {
  for (std::size_t p = 0; p < size(); ++p) std::cout << *(_data + p) << " ";
  std::cout << std::endl;
}

/////////////////////////// EXTERNAL FUNCTIONS ////////////////////////////////

// use  if complexity < some value.
void _multiply_MM(const s_type* A_data, const s_type* B_data, s_type* C_data,
                  std::size_t m, std::size_t n, std::size_t k) {
  if (A_data == nullptr) {
    throw ERROR_MSG("Data from Tensor A must be non-null.");
  }
  if (B_data == nullptr) {
    throw ERROR_MSG("Data from Tensor B must be non-null.");
  }
  if (C_data == nullptr) {
    throw ERROR_MSG("Data from Tensor C must be non-null.");
  }
  s_type alpha = 1.0;
  s_type beta = 0.0;
  cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, &alpha,
              A_data, std::max(1ul, k), B_data, std::max(1ul, n), &beta, C_data,
              std::max(1ul, n));
}

void _multiply_Mv(const s_type* A_data, const s_type* B_data, s_type* C_data,
                  std::size_t m, std::size_t k) {
  if (A_data == nullptr) {
    throw ERROR_MSG("Data from Tensor A must be non-null.");
  }
  if (B_data == nullptr) {
    throw ERROR_MSG("Data from Tensor B must be non-null.");
  }
  if (C_data == nullptr) {
    throw ERROR_MSG("Data from Tensor C must be non-null.");
  }
  s_type alpha = 1.0;
  s_type beta = 0.0;
  cblas_cgemv(CblasRowMajor, CblasNoTrans, m, k, &alpha, A_data,
              std::max(1ul, k), B_data, 1, &beta, C_data, 1);
}

void _multiply_vM(const s_type* A_data, const s_type* B_data, s_type* C_data,
                  std::size_t n, std::size_t k) {
  if (A_data == nullptr) {
    throw ERROR_MSG("Data from Tensor A must be non-null.");
  }
  if (B_data == nullptr) {
    throw ERROR_MSG("Data from Tensor B must be non-null.");
  }
  if (C_data == nullptr) {
    throw ERROR_MSG("Data from Tensor C must be non-null.");
  }

  s_type alpha = 1.0;
  s_type beta = 0.0;
  cblas_cgemv(CblasRowMajor, CblasTrans, k, n, &alpha, A_data, std::max(1ul, n),
              B_data, 1, &beta, C_data, 1);
}

void _multiply_vv(const s_type* A_data, const s_type* B_data, s_type* C_data,
                  std::size_t k) {
  if (A_data == nullptr) {
    throw ERROR_MSG("Data from Tensor A must be non-null.");
  }
  if (B_data == nullptr) {
    throw ERROR_MSG("Data from Tensor B must be non-null.");
  }
  if (C_data == nullptr) {
    throw ERROR_MSG("Data from Tensor C must be non-null.");
  }

  cblas_cdotu_sub(k, A_data, 1, B_data, 1, C_data);
}

void multiply(Tensor& A, Tensor& B, Tensor& C, s_type* scratch_copy) {
  if (scratch_copy == nullptr) {
    throw ERROR_MSG("Scratch copy must be non-null.");
  }

  if (A.data() == C.data()) {
    throw ERROR_MSG("A and C cannot be the same tensor: ",
                    C.tensor_to_string());
  }
  if (B.data() == C.data()) {
    throw ERROR_MSG("B and C cannot be the same tensor: ",
                    C.tensor_to_string());
  }

  std::chrono::high_resolution_clock::time_point t0, t1;
  std::chrono::duration<double> time_span;
  if (global::verbose > 1) t0 = std::chrono::high_resolution_clock::now();

  // Define left_indices, left_dim, right_indices, right_dim, and
  // common_indices, common_dim. Also C_size.
  std::vector<std::string> left_indices =
      _vector_subtraction(A.get_indices(), B.get_indices());
  std::vector<std::string> right_indices =
      _vector_subtraction(B.get_indices(), A.get_indices());
  std::vector<std::string> common_indices =
      _vector_intersection(A.get_indices(), B.get_indices());
  std::size_t left_dim = 1, right_dim = 1, common_dim = 1;
  for (std::size_t i = 0; i < left_indices.size(); ++i) {
    left_dim *= A.get_index_to_dimension().at(left_indices[i]);
  }
  for (std::size_t i = 0; i < right_indices.size(); ++i) {
    right_dim *= B.get_index_to_dimension().at(right_indices[i]);
  }
  for (std::size_t i = 0; i < common_indices.size(); ++i) {
    std::size_t a_dim = A.get_index_to_dimension().at(common_indices[i]);
    if (a_dim != B.get_index_to_dimension().at(common_indices[i])) {
      throw ERROR_MSG(
          "Common indices must have matching dimensions, but at index: ", i,
          ", the dimensions are A: ", a_dim,
          ", B: ", B.get_index_to_dimension().at(common_indices[i]), ".");
    }
    common_dim *= a_dim;
  }

  if (global::verbose > 1) {
    t1 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    std::cerr << "Time preparing variables: " << time_span.count() << "s\n";
  }

  // Check.
  if (left_dim * right_dim > C.capacity()) {
    throw ERROR_MSG("C: ", C.capacity(),
                    " doesn't have enough space for the product of A*B: ",
                    (left_dim * right_dim), ".");
  }

  if (global::verbose > 1) t0 = std::chrono::high_resolution_clock::now();

  // Reorder.
  std::vector<std::string> A_new_ordering =
      _vector_union(left_indices, common_indices);
  try {
    A.reorder(A_new_ordering, scratch_copy);
  } catch (const std::string& err_msg) {
    throw ERROR_MSG("Failed to call reorder(). Error:\n\t[", err_msg, "]");
  }

  if (global::verbose > 1) {
    t1 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    std::cerr << "R " << time_span.count() << " s\n";
    std::cerr << "Time reordering A: " << time_span.count() << "s\n";
    t0 = std::chrono::high_resolution_clock::now();
  }

  std::vector<std::string> B_new_ordering =
      _vector_union(common_indices, right_indices);
  try {
    B.reorder(B_new_ordering, scratch_copy);
  } catch (const std::string& err_msg) {
    throw ERROR_MSG("Failed to call reorder(). Error:\n\t[", err_msg, "]");
  }

  if (global::verbose > 1) {
    t1 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    std::cerr << "R " << time_span.count() << " s\n";
    std::cerr << "Time reordering B: " << time_span.count() << "s\n";
    t0 = std::chrono::high_resolution_clock::now();
  }

  // Multiply. Four cases: MxM, Mxv, vxM, vxv.
  if (left_indices.size() > 0 && right_indices.size() > 0) {
    try {
      _multiply_MM(A.data(), B.data(), C.data(), left_dim, right_dim,
                   common_dim);
    } catch (const std::string& err_msg) {
      throw ERROR_MSG("Failed to call _multiply_MM(). Error:\n\t[", err_msg,
                      "]");
    }
  } else if (left_indices.size() > 0 && right_indices.size() == 0) {
    try {
      _multiply_Mv(A.data(), B.data(), C.data(), left_dim, common_dim);
    } catch (const std::string& err_msg) {
      throw ERROR_MSG("Failed to call _multiply_Mv(). Error:\n\t[", err_msg,
                      "]");
    }
  } else if (left_indices.size() == 0 && right_indices.size() > 0) {
    // Very import to switch A and B to use cgemv with transpose for this case.
    try {
      _multiply_vM(B.data(), A.data(), C.data(), right_dim, common_dim);
    } catch (const std::string& err_msg) {
      throw ERROR_MSG("Failed to call _multiply_vM(). Error:\n\t[", err_msg,
                      "]");
    }
  } else if (left_indices.size() == 0 && right_indices.size() == 0) {
    try {
      _multiply_vv(A.data(), B.data(), C.data(), common_dim);
    } catch (const std::string& err_msg) {
      throw ERROR_MSG("Failed to call _multiply_vv(). Error:\n\t[", err_msg,
                      "]");
    }
  }

  if (global::verbose > 1) {
    t1 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    std::cerr << "M " << time_span.count() << " s\n";
    std::cerr << "Time multiplying A*B: " << time_span.count() << "s\n";
    t0 = std::chrono::high_resolution_clock::now();
  }

  // Set indices and dimensions of C.
  std::vector<std::string> C_indices =
      _vector_union(left_indices, right_indices);
  std::vector<std::size_t> C_dimensions(C_indices.size());
  for (std::size_t i = 0; i < left_indices.size(); ++i)
    C_dimensions[i] = A.get_index_to_dimension().at(left_indices[i]);
  for (std::size_t i = 0; i < right_indices.size(); ++i)
    C_dimensions[i + left_indices.size()] =
        B.get_index_to_dimension().at(right_indices[i]);
  C.set_indices(C_indices);
  try {
    C.set_dimensions(C_dimensions);
  } catch (const std::string& err_msg) {
    throw ERROR_MSG("Failed to call set_dimensions(). Error:\n\t[", err_msg,
                    "]");
  }
  C.generate_index_to_dimension();

  if (global::verbose > 1) {
    t1 = std::chrono::high_resolution_clock::now();
    time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
    std::cerr << "Time updating C's variables: " << time_span.count() << "s\n";
  }
}

// TODO: write tests for this function.
std::size_t result_size(Tensor& A, Tensor& B) {
  std::vector<std::string> left_indices =
      _vector_subtraction(A.get_indices(), B.get_indices());
  std::vector<std::string> right_indices =
      _vector_subtraction(B.get_indices(), A.get_indices());
  std::size_t left_dim = 1, right_dim = 1, result_dim;
  for (std::size_t i = 0; i < left_indices.size(); ++i) {
    left_dim *= A.get_index_to_dimension().at(left_indices[i]);
  }
  for (std::size_t i = 0; i < right_indices.size(); ++i) {
    right_dim *= B.get_index_to_dimension().at(right_indices[i]);
  }
  result_dim = left_dim * right_dim;
  return result_dim;
}

// TODO: write tests for this function.
void bundle_between(Tensor& A, Tensor& B, std::string bundled_index,
                    s_type* scratch_copy) {
  std::vector<std::string> left_indices =
      _vector_subtraction(A.get_indices(), B.get_indices());
  std::vector<std::string> right_indices =
      _vector_subtraction(B.get_indices(), A.get_indices());
  std::vector<std::string> common_indices =
      _vector_intersection(A.get_indices(), B.get_indices());
  std::vector<std::string> A_new_ordering =
      _vector_union(left_indices, common_indices);
  std::vector<std::string> B_new_ordering =
      _vector_union(common_indices, right_indices);
  try {
    A.reorder(A_new_ordering, scratch_copy);
  } catch (const std::string& err_msg) {
    throw ERROR_MSG("Failed to call reorder(). Error:\n\t[", err_msg, "]");
  }
  try {
    B.reorder(B_new_ordering, scratch_copy);
  } catch (const std::string& err_msg) {
    throw ERROR_MSG("Failed to call reorder(). Error:\n\t[", err_msg, "]");
  }
  try {
    A.bundle(common_indices, bundled_index);
  } catch (const std::string& err_msg) {
    throw ERROR_MSG("Failed to call bundle(). Error:\n\t[", err_msg, "]");
  }
  try {
    B.bundle(common_indices, bundled_index);
  } catch (const std::string& err_msg) {
    throw ERROR_MSG("Failed to call bundle(). Error:\n\t[", err_msg, "]");
  }
}

// Split it in parts, as before? Nah, it was all about generating small maps.
// Here I am generating THE map. This is done only once per map anyway.
void _generate_binary_reordering_map(
    const std::vector<std::size_t>& map_old_to_new_idxpos,
    std::vector<std::size_t>& map_old_to_new_position) {
  std::size_t dim = 2;  // Hard coded!
  std::size_t num_indices = map_old_to_new_idxpos.size();

  // Check
  if (num_indices == 0) throw ERROR_MSG("Number of indices cannot be zero.");

  // Check
  if ((std::size_t)std::pow(dim, num_indices) !=
      map_old_to_new_position.size()) {
    throw ERROR_MSG(
        "Size of map: ", map_old_to_new_position.size(),
        " must be equal to 2^num_indices: ", std::pow(dim, num_indices), ".");
  }

  // Define super dimensions. See _naive_reorder().
  std::vector<std::size_t> old_dimensions(num_indices, dim);
  std::vector<std::size_t> new_dimensions(num_indices, dim);
  std::vector<std::size_t> old_super_dimensions(num_indices);
  std::vector<std::size_t> new_super_dimensions(num_indices);
  old_super_dimensions[num_indices - 1] = 1;
  new_super_dimensions[num_indices - 1] = 1;

  if (num_indices >= 2)
    for (std::size_t i = num_indices; --i;) {
      old_super_dimensions[i - 1] = old_super_dimensions[i] * dim;
      new_super_dimensions[i - 1] = new_super_dimensions[i] * dim;
    }

  // Iterate and generate map.
  std::vector<std::size_t> old_counter(num_indices, 0);

  while (true) {
    std::size_t po{0}, pn{0};  // Position of the data, old and new.

    for (std::size_t i = 0; i < num_indices; ++i) {
      po += old_super_dimensions[i] * old_counter[i];
      pn += new_super_dimensions[map_old_to_new_idxpos[i]] * old_counter[i];
    }
    map_old_to_new_position[po] = pn;

    bool complete{true};
    for (std::size_t j = num_indices; j--;) {
      if (++old_counter[j] < old_dimensions[j]) {
        complete = false;
        break;
      } else
        old_counter[j] = 0;
    }
    if (complete) break;
  }
}

// convert int vector to string
std::string _int_vector_to_string(std::vector<std::size_t> input) {
  std::ostringstream temp;
  std::string output;
  if (!input.empty()) {
    std::copy(input.begin(), input.end() - 1,
              std::ostream_iterator<std::size_t>(temp, ", "));
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

std::string _reordering_to_string(
    const std::vector<std::size_t>& map_old_to_new_idxpos,
    const std::vector<std::size_t>& old_dimensions) {
  std::size_t num_indices = map_old_to_new_idxpos.size();
  std::string name("");
  for (std::size_t i = 0; i < num_indices; ++i) name += _ALPHABET[i];
  name += "->";
  for (std::size_t i = 0; i < num_indices; ++i)
    name += _ALPHABET[map_old_to_new_idxpos[i]];
  for (std::size_t i = 0; i < num_indices; ++i)
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
  for (std::size_t i = 0; i < v.size(); ++i) {
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
  for (std::size_t i = 0; i < v.size(); ++i) temp[i] = v[i];
  for (std::size_t i = 0; i < w.size(); ++i) temp[i + v.size()] = w[i];
  return temp;
}

}  // namespace qflex
