/**
* @file mkl_tensor.h
* Definition of the MKLTensor class, which implements tensors with MKL's
* matrix multiplication and a self-written entry reordering algorithm.
* @see https://github.com/benjaminvillalonga/optimized_parallel_QC_with_TN
*
* @author Benjamin Villalonga
* @date Created: August 2018
* @date Modified: August 2018
*/

#ifndef MKL_TENSOR_H
#define MKL_TENSOR_H

#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>
#include <complex>
#include <mkl.h>

using namespace std;

/**
* Scalar type.
*/
typedef complex<float> s_type;

/**
* Represents an MKLTensor.
*/
class MKLTensor
{
  public:
    /**
    * Creates an uninitialized MKLTensor.
    */
    MKLTensor();

    /**
    * Creates an MKLTensor. New space is allocated.
    * @param indices vector<string> with the names of the indices in order.
    * @param dimensions vector<size_t> with the ordered dimensions of the indices.
    */
    MKLTensor(vector<string> indices, vector<size_t> dimensions);

    /**
    * Creates an MKLTensor. New space is allocated and filled with a copy of
    * the vector's data. Useful for small tensors where the copying time is
    * negligible.
    * @param indices vector<string> with the names of the indices in order.
    * @param dimensions vector<size_t> with the ordered dimensions of the indices.
    * @param data vector<s_type> with the data to be copied. It has to match
    * in length the dimension of the MKLTensor, as given by the dimensions.
    */
    MKLTensor(vector<string> indices, vector<size_t> dimensions,
              const vector<s_type> & data);

    /**
    * Creates an MKLTensor. A pointer to the data is passed.
    * @param indices vector<string> with the names of the indices in order.
    * @param dimensions vector<int> with the ordered dimensions of the indices.
    * @param data pointer to the data of the tensor. It is responsibility of
    * the user to provide enough allocated memory to store the MKLTensor.
    */
    MKLTensor(vector<string> indices, vector<size_t> dimensions, s_type * data);

    /**
    * Copy constructor: creates a new MKLTensor that is a copy of another.
    * @param other MKLTensor to be copied.
    */
    MKLTensor(const MKLTensor & other);

    /**
    * Destructor: frees all memory associated with a given MKLTensor object.
    * Invoked by the syste.
    */
    ~MKLTensor();

    /**
    * Assignment operator for setting two MKLTensor equal to one another.
    * It is responsibility of the user to copy onto an MKLTensor with the same
    * total dimension as other. If there is space allocated, no new space will
    * be allocated. Changing the size of an MKLTensor is not allowed, although
    * if this MKLTensor has at least as much space allocated as other, then
    * everything will run smoothly, with a non-optimal usage of memory.
    * @param other MKLTensor to copy into the current MKLTensor.
    * @return The current MKLTensor for assignment chaining.
    */
    const MKLTensor & operator=(const MKLTensor & other);

    /**
    * Get inidices.
    * @return const reference to vector<string> of indices.
    */
    const vector<string> & get_indices() const;

    /**
    * Set inidices. This function is deprecated. Use rename_index() or
    * set_dimensions_and_indices().
    * @param const reference to vector<string> of indices.
    */
    void set_indices(const vector<string> & indices);

    /**
    * Get dimensions.
    * @return const reference to vector<int> of dimensions.
    */
    const vector<size_t> & get_dimensions() const;

    /**
    * Set dimensions. This function is deprecated. Use rename_index() or
    * set_dimensions_and_indices().
    * @param const reference to vector<int> of dimensions.
    */
    void set_dimensions(const vector<size_t> & dimensions);

    /**
    * Set dimensions and indices.
    * @param const reference to vector<string> of indices.
    * @param const reference to vector<int> of dimensions.
    */
    void set_indices_and_dimensions(const vector<string> & indices,
                                    const vector<size_t> & dimensions);

    /**
    * Get index_to_dimension dictionary (or unordered_map).
    * @return const reference to unordered_map of index to dimensions.
    */
    const unordered_map<string,size_t> & get_index_to_dimension() const;

    /**
    * Generate index_to_dimension map from the object's indices and
    * dimensions.
    */
    void generate_index_to_dimension();

    /**
    * Get size (total dimension) of the MKLTensor.
    * @return int with the number of elements in the MKLTensor.
    */
    size_t size() const;

    /**
    * Get data.
    * @return s_type * to the data of the MKLTensor.
    */
    s_type * data();

    /**
    * Get data.
    * @return const s_type * to the data of the MKLTensor.
    */
    const s_type * data() const;

    /**
    * Project tensor to a value on a particular index. The current MKLTensor
    * has to be ordered with the index to project in first place. This forces
    * the user to follow efficient practices, and be aware of the costs
    * involved in operations.
    * @param index string with the name of the index to fix.
    * @param index_value int with the value to which the index is fixed.
    * @param projection_tensor Reference to the MKLTensor to which the current
    * MKLTensor will be projected. It is responsibility of the user to provide
    * enough allocated space to store the projected MKLTensor. A trivial tensor
    * with the right amount of allocated space can be passed. The indices and
    * dimensions will be initialized correctly from the project() function.
    */
    void project(string index, size_t index_value,
                 MKLTensor & projection_tensor) const;

    /**
    * Rename an index of the MKLTensor.
    * @param old_name string with the old name of an index.
    * @param new_name string with the new name of an index.
    */
    void rename_index(string old_name, string new_name);

    /**
    * Bundle indices onto a single big index. The ordering before the bundling
    * operation has to be the right one, with the indices to bundle properly
    * ordered and contiguously placed. This forces the user to be aware of the
    * cost of each operation.
    * @param indices_to_bundle vector<string> with the names of the indices
    * to bundle. They have to be contiguous and on the right order.
    * @param bundled_index string with the name of the big new index.
    */
    void bundle(vector<string> indices_to_bundle, string bundled_index);

    /**
    * Reorder the indices of the MKLTensor. This is an intensive operation,
    * And the second leading bottleneck in the contraction of a tensor network
    * after the multiplication of properly ordered tensors.
    * @param new_ordering vector<string> with the new ordering of the indices.
    * @param scratch_copy Pointer to s_type with space allocated for scratch
    * work. Allocate at least as much space as the size of the MKLTensor.
    */
    void reorder(vector<string> new_ordering, s_type * scratch_copy);

    /**
    * Multiply MKLTensor by scalar.
    * @param s_type scalar that multiplies the MKLTensor.
    */
    void scalar_multiply(s_type scalar);

    /**
    * Compute the L2 norm (squared) of this MKLTensor.
    */
    double tensor_norm() const;

    /**
    * Count the number of zeroed out entries in the tensor. This is useful
    * to debug contractions that zero out entries due to the range of values
    * s_type admits.
    */
    size_t num_zeros() const;

    /**
    * Prints information about the MKLTensor.
    */
    void print() const;

    /**
    * Prints the data of the MKLTensor.
    */
    void print_data() const;

  private:
    // Storage.
    vector<string> _indices;
    vector<size_t> _dimensions;
    unordered_map<string,size_t> _index_to_dimension;
    s_type * _data;

    // Private helper functions.
    /**
    * Helper function for the constructor, copy constructor and assignment
    * operator. It only initializes indices and dimensions; it does nothing
    * to the pointer to the data.
    * @param indices vector<string> with the names of the indices in order.
    * @param dimensions vector<int> with the ordered dimensions of the indices.
    */
    void _init(const vector<string> & indices, const vector<size_t> & dimensions);

    /**
    * Helper function for the destructor. Clear memory.
    */
    void _clear();

    /**
    * Helper function for the copy constructor and the assignment operator.
    * It is responsibility of the user to copy onto an MKLTensor with the same
    * total dimension as other. If there is space allocated, no new space will
    * be allocated. Changing the size of an MKLTensor is not allowed.
    * @param other MKLTensor to copy into the current MKLTensor.
    */
    void _copy(const MKLTensor & other);

    /**
    * Helper function for reorder(). It is called when smart reordering doesn't
    * apply.
    * @param new_ordering vector<string> with the new ordering.
    * @param scratch_copy Pointer to an s_type array for scracth copying work.
    */
    void _naive_reorder(vector<string> new_ordering, s_type * scratch_copy);

    /**
    * Helper function for reorder(). It is called when smart reordering
    * applies.
    * @param new_ordering vector<string> with the new ordering.
    * @param scratch_copy Pointer to an s_type array for scratch copying work.
    */
    void _fast_reorder(vector<string> new_ordering, s_type * scratch_copy);

    /**
    * Helper function for reorder(). Only right moves are taken. For some
    * reason, allocating the small chunks in the function is faster than
    * passing a big space and using different chunks of it.
    * @param old_ordering const reference to a vector<string> with the old
    * ordering of the right indices.
    * @param new_ordering const reference to a vector<string> with the new
    * ordering of the right indices.
    * @param num_indices_right number of indices being reordered on the right.
    */
    void _right_reorder(const vector<string> & old_ordering,
                        const vector<string> & new_ordering,
                        int num_indices_right);

    /**
    * Helper function for reorder(). Only left moves are taken. As oposed
    * to right moves, for left moves passing a scratch space is faster than
    * allocating memory dynamically, at least for the current approach of
    * copying the entire array and then reordering back onto data.
    * @param old_ordering const reference to a vector<string> with the old
    * ordering of the left indices.
    * @param new_ordering const reference to a vector<string> with the new
    * ordering of the left indices.
    * @param num_indices_right number of indices being left on the right.
    * @param scratch_copy Pointer to an s_type array where scratch copy work
    * is done.
    */
    void _left_reorder(const vector<string> & old_ordering,
                       const vector<string> & new_ordering,
                       int num_indices_right, s_type * scratch_copy);

};

/**
* Call MKL matrix x matrix multiplication C = A * B, or self-written if the
* matrices are small.
* @param A_data const pointer to s_type array with the data of matrix A.
* @param B_data const pointer to s_type array with the data of matrix B.
* @param C_data pointer to s_type array with the data of matrix C.
* @param m int with the left dimension of A.
* @param n int with the right dimension of B.
* @param k int with the left dimension of C.
*/
void _multiply_MM(const s_type * A_data, const s_type * B_data,
                 s_type * C_data, int m, int n, int k);

/**
* Call MKL matrix x vector multiplication C = A * B, or self-written if the
* objects are small.
* @param A_data const pointer to s_type array with the data of matrix A.
* @param B_data const pointer to s_type array with the data of vector B.
* @param C_data pointer to s_type array with the data of vector C.
* @param m int with the left dimension of A.
* @param k int with the left dimension of C.
*/
void _multiply_Mv(const s_type * A_data, const s_type * B_data,
                 s_type * C_data, int m, int k);

/**
* Call MKL vector x matrix multiplication C = A * B, or self-written if the
* objects are small.
* @param A_data const pointer to s_type array with the data of vector A.
* @param B_data const pointer to s_type array with the data of matrix B.
* @param C_data pointer to s_type array with the data of vector C.
* @param n int with the right dimension of B.
* @param k int with the left dimension of C.
*/
void _multiply_vM(const s_type * A_data, const s_type * B_data,
                 s_type * C_data, int n, int k);

/**
* Call MKL vector x vector multiplication C = A * B, or self-written if the
* matrices are small.
* @param A_data const pointer to s_type array with the data of vector A.
* @param B_data const pointer to s_type array with the data of vector B.
* @param C_data pointer to s_type with the data of scalar C.
* @param k int with the left dimension of C.
*/
void _multiply_vv(const s_type * A_data, const s_type * B_data,
                 s_type * C_data, int k);

/**
* Implement the tensor multiplication C = A * B. The ordering of the incides of
* A and B can have a big impact on performance. All three tensors have to be
* distinct, and have their own allocated spaces. It is erroneous to call:
* multiply(A, B, A) or multiply(A, B, B). In the case that indices are not
* ordered correctly, reordering will be applied, which can be costly; relative
* ordering will be preserved for the final indices, and A's ordering for the
* contracted indices.
* @param A Reference to MKLTensor A. A can be reordered, and therfore modified.
* @param B Reference to MKLTensor B. B can be reordered, and therfore modified.
* @param C Reference to MKLTensor C, where the product of A and B will be
* stored. It is responsibility of the user to allocate sufficient memory for C.
* C can allocate more memory than needed, which allows for flexible reuse of
* MKLTensors; this should be used only occasionally in cases where critical
* optimization is needed. C can be initialized trivially, and the function
* multiply will initialize indices and dimensions correctly.
* @param scratch_copy Pointer to s_type array for scratch work while
* reordering. It has to allocate at least as much max(A.size(), B.size())
* memory.
*/
void multiply(MKLTensor & A, MKLTensor & B, MKLTensor & C,
              s_type * scratch_copy);

/**
* Creates a reordering map for the data of a tensor with binary indices
* between old_ordering and new_ordering.
* @param map_old_to_new_idxpos const reference to a vector<string> of new
* indices.
* @param map_old_to_new_position Reference to a vector<int> where the map will
* be stored. The size
* of the vector has to be equal to 2^(old_indices.size()).
* \todo Think about passing some scratch space to store small maps. Probably
* not worth it, since they are small, and the big ones will be memoized anyway.
* If the NO_MEMO_MAPS flag is passed, then it might make sense. Not a priority.
* \todo When dimensions relax to something different than 2 in the future, then
* old_dimensions have to be passed as an argument.
*/
void _generate_binary_reordering_map(
                       const vector<int> & map_old_to_new_idxpos,
                       vector<int> & map_old_to_new_position);

/**
* Generates the standard name of the reordering as a string:
* "abc...->fbe...,dim_a,dim_b,...".
* @param map_old_to_new_idxpos const reference to a vector<int> index mapping
* (old to new).
* @param old_dimensions const reference to a vector<int> with the old
* dimensions.
* @return string with the standard name of the ordering.
*/
string _reordering_to_string(const vector<int> & map_old_to_new_idxpos,
                             const vector<size_t> & old_dimensions);

/**
* Checks whether a particular string is in a vector<string>.
* @param s const reference to string element.
* @param v const reference to vector<string> where the element might be in.
*/
bool _string_in_vector(const string & s, const vector<string> & v);

/**
* Checks whether a particular string is in a vector<string>.
* @param s const reference to string element.
* @param v const reference to vector<string> where the element might be in.
*/
bool _vector_s_in_vector_s(const vector<string> & v, const vector<string> & w);

/**
* Gets the intersection between two vector<string>.
* @param v const reference to first vector<string>.
* @param w const reference to second vector<string>.
* @return vector<string> intersection of v and w. Preserves the ordering of v.
*/
vector<string> _vector_intersection(const vector<string> & v,
                                   const vector<string> & w);

/**
* Gets the union between two vector<string>.
* @param v const reference to first vector<string>.
* @param w const reference to second vector<string>.
* @return vector<string> intersection of v and w. The ordering is that one of.
* v followed by the elements that only appear in w, in w's order.
*/
vector<string> _vector_union(const vector<string> & v,
                            const vector<string> & w);

/**
* Gets the subtraction of two vector<string>.
* @param v const reference to first vector<string>.
* @param w const reference to second vector<string>.
* @return vector<string> with v-w in v's order.
*/
vector<string> _vector_subtraction(const vector<string> & v,
                                  const vector<string> & w);

/**
* Gets concatenation of two vectors.
* @param v const reference to first vector<string>.
* @param w const reference to second vector<string>.
* @return vector<string> with concatenation of v and w.
*/
vector<string> _vector_concatenation(const vector<string> & v,
                                     const vector<string> & w);

#endif
