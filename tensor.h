/**
 * @file tensor.h
 * Definition of the Tensor class, which implements tensors with
 * matrix multiplication and a self-written entry reordering algorithm.
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

#ifndef TENSOR_H
#define TENSOR_H

#ifdef MKL_TENSOR
#include <mkl.h>
#else
#include <gsl/gsl_cblas.h>
#endif

#include <complex>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace qflex {

/**
 * Scalar type.
 */
typedef std::complex<float> s_type;

/**
 * Represents an Tensor.
 */
class Tensor {
 public:
  /**
   * Creates an uninitialized Tensor.
   */
  Tensor();

  /**
   * Creates an Tensor. New space is allocated.
   * @param indices std::vector<std::string> with the names of the indices in
   * order.
   * @param dimensions std::vector<size_t> with the ordered dimensions of the
   * indices.
   */
  Tensor(std::vector<std::string> indices, std::vector<size_t> dimensions);

  /**
   * Creates an Tensor. New space is allocated and filled with a copy of
   * the vector's data. Useful for small tensors where the copying time is
   * negligible.
   * @param indices std::vector<std::string> with the names of the indices in
   * order.
   * @param dimensions std::vector<size_t> with the ordered dimensions of the
   * indices.
   * @param data std::vector<s_type> with the data to be copied. It has to
   * match in length the dimension of the Tensor, as given by the
   * dimensions.
   */
  Tensor(std::vector<std::string> indices, std::vector<size_t> dimensions,
         const std::vector<s_type>& data);

  /**
   * Creates an Tensor. A pointer to the data is passed.
   * @param indices std::vector<std::string> with the names of the indices in
   * order.
   * @param dimensions std::vector<int> with the ordered dimensions of the
   * indices.
   * @param data pointer to the data of the tensor. It is responsibility of
   * the user to provide enough allocated memory to store the Tensor.
   */
  Tensor(std::vector<std::string> indices, std::vector<size_t> dimensions,
         s_type* data);

  /**
   * Copy constructor: creates a new Tensor that is a copy of another.
   * @param other Tensor to be copied.
   */
  Tensor(const Tensor& other);

  /**
   * Destructor: frees all memory associated with a given Tensor object.
   * Invoked by the syste.
   */
  ~Tensor();

  /**
   * Assignment operator for setting two Tensor equal to one another.
   * It is responsibility of the user to copy onto an Tensor with the same
   * total dimension as other. If there is space allocated, no new space will
   * be allocated. Changing the size of an Tensor is not allowed, although
   * if this Tensor has at least as much space allocated as other, then
   * everything will run smoothly, with a non-optimal usage of memory.
   * @param other Tensor to copy into the current Tensor.
   * @return The current Tensor for assignment chaining.
   */
  const Tensor& operator=(const Tensor& other);

  /**
   * Get inidices.
   * @return const reference to std::vector<std::string> of indices.
   */
  const std::vector<std::string>& get_indices() const;

  /**
   * Set inidices. This function is deprecated. Use rename_index() or
   * set_dimensions_and_indices().
   * @param const reference to std::vector<std::string> of indices.
   */
  void set_indices(const std::vector<std::string>& indices);

  /**
   * Get dimensions.
   * @return const reference to std::vector<int> of dimensions.
   */
  const std::vector<size_t>& get_dimensions() const;

  /**
   * Set dimensions. This function is deprecated. Use rename_index() or
   * set_dimensions_and_indices().
   * @param const reference to std::vector<int> of dimensions.
   */
  void set_dimensions(const std::vector<size_t>& dimensions);

  /**
   * Set dimensions and indices.
   * @param const reference to std::vector<std::string> of indices.
   * @param const reference to std::vector<int> of dimensions.
   */
  void set_indices_and_dimensions(const std::vector<std::string>& indices,
                                  const std::vector<size_t>& dimensions);

  /**
   * Get index_to_dimension dictionary (or unordered_map).
   * @return const reference to unordered_map of index to dimensions.
   */
  const std::unordered_map<std::string, size_t>& get_index_to_dimension() const;

  /**
   * Generate index_to_dimension map from the object's indices and
   * dimensions.
   */
  void generate_index_to_dimension();

  /**
   * Get size (total dimension) of the Tensor.
   * @return int with the number of elements in the Tensor.
   */
  size_t size() const;

  /**
   * Get the allocated space of the Tensor. This tensor can be resized to
   * any shape with total dimension less than this value.
   * @return int with the capacity of the Tensor.
   */
  size_t capacity() const;

  /**
   * Get data.
   * @return s_type * to the data of the Tensor.
   */
  s_type* data();

  /**
   * Get data.
   * @return const s_type * to the data of the Tensor.
   */
  const s_type* data() const;

  /**
   * Project tensor to a value on a particular index. The current Tensor
   * has to be ordered with the index to project in first place. This forces
   * the user to follow efficient practices, and be aware of the costs
   * involved in operations.
   * @param index std::string with the name of the index to fix.
   * @param index_value int with the value to which the index is fixed.
   * @param projection_tensor Reference to the Tensor to which the current
   * Tensor will be projected. It is responsibility of the user to provide
   * enough allocated space to store the projected Tensor. A trivial tensor
   * with the right amount of allocated space can be passed. The indices and
   * dimensions will be initialized correctly from the project() function.
   */
  void project(std::string index, size_t index_value,
               Tensor& projection_tensor) const;

  /**
   * Rename an index of the Tensor.
   * @param old_name std::string with the old name of an index.
   * @param new_name std::string with the new name of an index.
   */
  void rename_index(std::string old_name, std::string new_name);

  /**
   * Bundle indices onto a single big index. The ordering before the bundling
   * operation has to be the right one, with the indices to bundle properly
   * ordered and contiguously placed. This forces the user to be aware of the
   * cost of each operation.
   * @param indices_to_bundle std::vector<std::string> with the names of the
   * indices to bundle. They have to be contiguous and on the right order.
   * @param bundled_index std::string with the name of the big new index.
   */
  void bundle(std::vector<std::string> indices_to_bundle,
              std::string bundled_index);

  /**
   * Reorder the indices of the Tensor. This is an intensive operation,
   * And the second leading bottleneck in the contraction of a tensor network
   * after the multiplication of properly ordered tensors.
   * @param new_ordering std::vector<std::string> with the new ordering of the
   * indices.
   * @param scratch_copy Pointer to s_type with space allocated for scratch
   * work. Allocate at least as much space as the size of the Tensor.
   */
  void reorder(std::vector<std::string> new_ordering, s_type* scratch_copy);

  /**
   * Multiply Tensor by scalar.
   * @param s_type scalar that multiplies the Tensor.
   */
  void scalar_multiply(s_type scalar);

  /**
   * Compute the L2 norm (squared) of this Tensor.
   */
  double tensor_norm() const;

  /**
   * Count the number of zeroed out entries in the tensor. This is useful
   * to debug contractions that zero out entries due to the range of values
   * s_type admits.
   */
  size_t num_zeros() const;

  /**
   * Prints information about the Tensor.
   */
  void print() const;

  /**
   * Prints the data of the Tensor.
   */
  void print_data() const;

 private:
  // Storage.
  std::vector<std::string> _indices;
  std::vector<size_t> _dimensions;
  std::unordered_map<std::string, size_t> _index_to_dimension;
  s_type* _data;

  // Allocated data space. This value does not change after initialization.
  size_t _capacity;

  // Private helper functions.
  /**
   * Helper function for the constructor, copy constructor and assignment
   * operator. It only initializes indices and dimensions; it does nothing
   * to the pointer to the data.
   * @param indices std::vector<std::string> with the names of the indices in
   * order.
   * @param dimensions std::vector<int> with the ordered dimensions of the
   * indices.
   */
  void _init(const std::vector<std::string>& indices,
             const std::vector<size_t>& dimensions);

  /**
   * Helper function for the destructor. Clear memory.
   */
  void _clear();

  /**
   * Helper function for the copy constructor and the assignment operator.
   * It is responsibility of the user to copy onto an Tensor with the same
   * total dimension as other. If there is space allocated, no new space will
   * be allocated. Changing the size of an Tensor is not allowed.
   * @param other Tensor to copy into the current Tensor.
   */
  void _copy(const Tensor& other);

  /**
   * Helper function for reorder(). It is called when smart reordering doesn't
   * apply.
   * @param new_ordering std::vector<std::string> with the new ordering.
   * @param scratch_copy Pointer to an s_type array for scracth copying work.
   */
  void _naive_reorder(std::vector<std::string> new_ordering,
                      s_type* scratch_copy);

  /**
   * Helper function for reorder(). It is called when smart reordering
   * applies.
   * @param new_ordering std::vector<std::string> with the new ordering.
   * @param scratch_copy Pointer to an s_type array for scratch copying work.
   */
  void _fast_reorder(std::vector<std::string> new_ordering,
                     s_type* scratch_copy);

  /**
   * Helper function for reorder(). Only right moves are taken. For some
   * reason, allocating the small chunks in the function is faster than
   * passing a big space and using different chunks of it.
   * @param old_ordering const reference to a std::vector<std::string> with the
   * old ordering of the right indices.
   * @param new_ordering const reference to a std::vector<std::string> with the
   * new ordering of the right indices.
   * @param num_indices_right number of indices being reordered on the right.
   */
  void _right_reorder(const std::vector<std::string>& old_ordering,
                      const std::vector<std::string>& new_ordering,
                      int num_indices_right);

  /**
   * Helper function for reorder(). Only left moves are taken. As oposed
   * to right moves, for left moves passing a scratch space is faster than
   * allocating memory dynamically, at least for the current approach of
   * copying the entire array and then reordering back onto data.
   * @param old_ordering const reference to a std::vector<std::string> with the
   * old ordering of the left indices.
   * @param new_ordering const reference to a std::vector<std::string> with the
   * new ordering of the left indices.
   * @param num_indices_right number of indices being left on the right.
   * @param scratch_copy Pointer to an s_type array where scratch copy work
   * is done.
   */
  void _left_reorder(const std::vector<std::string>& old_ordering,
                     const std::vector<std::string>& new_ordering,
                     int num_indices_right, s_type* scratch_copy);
};

/**
 * Call  matrix x matrix multiplication C = A * B, or self-written if the
 * matrices are small.
 * @param A_data const pointer to s_type array with the data of matrix A.
 * @param B_data const pointer to s_type array with the data of matrix B.
 * @param C_data pointer to s_type array with the data of matrix C.
 * @param m int with the left dimension of A.
 * @param n int with the right dimension of B.
 * @param k int with the left dimension of C.
 */
void _multiply_MM(const s_type* A_data, const s_type* B_data, s_type* C_data,
                  int m, int n, int k);

/**
 * Call matrix x vector multiplication C = A * B, or self-written if the
 * objects are small.
 * @param A_data const pointer to s_type array with the data of matrix A.
 * @param B_data const pointer to s_type array with the data of vector B.
 * @param C_data pointer to s_type array with the data of vector C.
 * @param m int with the left dimension of A.
 * @param k int with the left dimension of C.
 */
void _multiply_Mv(const s_type* A_data, const s_type* B_data, s_type* C_data,
                  int m, int k);

/**
 * Call vector x matrix multiplication C = A * B, or self-written if the
 * objects are small.
 * @param A_data const pointer to s_type array with the data of vector A.
 * @param B_data const pointer to s_type array with the data of matrix B.
 * @param C_data pointer to s_type array with the data of vector C.
 * @param n int with the right dimension of B.
 * @param k int with the left dimension of C.
 */
void _multiply_vM(const s_type* A_data, const s_type* B_data, s_type* C_data,
                  int n, int k);

/**
 * Call vector x vector multiplication C = A * B, or self-written if the
 * matrices are small.
 * @param A_data const pointer to s_type array with the data of vector A.
 * @param B_data const pointer to s_type array with the data of vector B.
 * @param C_data pointer to s_type with the data of scalar C.
 * @param k int with the left dimension of C.
 */
void _multiply_vv(const s_type* A_data, const s_type* B_data, s_type* C_data,
                  int k);

/**
 * Implement the tensor multiplication C = A * B. The ordering of the incides of
 * A and B can have a big impact on performance. All three tensors have to be
 * distinct, and have their own allocated spaces. It is erroneous to call:
 * multiply(A, B, A) or multiply(A, B, B). In the case that indices are not
 * ordered correctly, reordering will be applied, which can be costly; relative
 * ordering will be preserved for the final indices, and A's ordering for the
 * contracted indices.
 * @param A Reference to Tensor A. A can be reordered, and therfore modified.
 * @param B Reference to Tensor B. B can be reordered, and therfore modified.
 * @param C Reference to Tensor C, where the product of A and B will be
 * stored. It is responsibility of the user to allocate sufficient memory for C.
 * C can allocate more memory than needed, which allows for flexible reuse of
 * Tensors; this should be used only occasionally in cases where critical
 * optimization is needed. C can be initialized trivially, and the function
 * multiply will initialize indices and dimensions correctly.
 * @param scratch_copy Pointer to s_type array for scratch work while
 * reordering. It has to allocate at least as much max(A.size(), B.size())
 * memory.
 */
void multiply(Tensor& A, Tensor& B, Tensor& C, s_type* scratch_copy);

/**
 * Creates a reordering map for the data of a tensor with binary indices
 * between old_ordering and new_ordering.
 * @param map_old_to_new_idxpos const reference to a std::vector<std::string> of
 * new indices.
 * @param map_old_to_new_position Reference to a std::vector<int> where the map
 * will be stored. The size of the vector has to be equal to
 * 2^(old_indices.size()). \todo Think about passing some scratch space to store
 * small maps. Probably not worth it, since they are small, and the big ones
 * will be memoized anyway. If the NO_MEMO_MAPS flag is passed, then it might
 * make sense. Not a priority. \todo When dimensions relax to something
 * different than 2 in the future, then old_dimensions have to be passed as an
 * argument.
 */
void _generate_binary_reordering_map(
    const std::vector<int>& map_old_to_new_idxpos,
    std::vector<int>& map_old_to_new_position);

/**
 * Converts an int vector into a string
 * @param input vector<int> int vector to convert.
 * @return string containing input vector contents
 **/
std::string _int_vector_to_string(std::vector<int> input);

/**
 * Converts a string vector into a string
 * @param input vector<string> string vector to convert.
 * @return string containing input vector contents
 **/
std::string _string_vector_to_string(std::vector<std::string> input);

/**
 * Generates the standard name of the reordering as a std::string:
 * "abc...->fbe...,dim_a,dim_b,...".
 * @param map_old_to_new_idxpos const reference to a std::vector<int> index
 * mapping (old to new).
 * @param old_dimensions const reference to a std::vector<int> with the old
 * dimensions.
 * @return std::string with the standard name of the ordering.
 */
std::string _reordering_to_string(const std::vector<int>& map_old_to_new_idxpos,
                                  const std::vector<size_t>& old_dimensions);

/**
 * Checks whether a particular std::string is in a std::vector<std::string>.
 * @param s const reference to std::string element.
 * @param v const reference to std::vector<std::string> where the element might
 * be in.
 */
bool _string_in_vector(const std::string& s, const std::vector<std::string>& v);

/**
 * Checks whether a particular std::string is in a std::vector<std::string>.
 * @param s const reference to std::string element.
 * @param v const reference to std::vector<std::string> where the element might
 * be in.
 */
bool _vector_s_in_vector_s(const std::vector<std::string>& v,
                           const std::vector<std::string>& w);

/**
 * Gets the intersection between two std::vector<std::string>.
 * @param v const reference to first std::vector<std::string>.
 * @param w const reference to second std::vector<std::string>.
 * @return std::vector<std::string> intersection of v and w. Preserves the
 * ordering of v.
 */
std::vector<std::string> _vector_intersection(
    const std::vector<std::string>& v, const std::vector<std::string>& w);

/**
 * Gets the union between two std::vector<std::string>.
 * @param v const reference to first std::vector<std::string>.
 * @param w const reference to second std::vector<std::string>.
 * @return std::vector<std::string> intersection of v and w. The ordering is
 * that one of. v followed by the elements that only appear in w, in w's order.
 */
std::vector<std::string> _vector_union(const std::vector<std::string>& v,
                                       const std::vector<std::string>& w);

/**
 * Gets the subtraction of two std::vector<std::string>.
 * @param v const reference to first std::vector<std::string>.
 * @param w const reference to second std::vector<std::string>.
 * @return std::vector<std::string> with v-w in v's order.
 */
std::vector<std::string> _vector_subtraction(const std::vector<std::string>& v,
                                             const std::vector<std::string>& w);

/**
 * Gets concatenation of two vectors.
 * @param v const reference to first std::vector<std::string>.
 * @param w const reference to second std::vector<std::string>.
 * @return std::vector<std::string> with concatenation of v and w.
 */
std::vector<std::string> _vector_concatenation(
    const std::vector<std::string>& v, const std::vector<std::string>& w);

}  // namespace qflex

#endif  // TENSOR_H
