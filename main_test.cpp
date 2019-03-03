#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <complex>
#include <ctime>
#include <chrono>

#include <iostream>
#include <fstream>
#include <sstream>

#include "talshxx.hpp"
#include "talsh_wrapper.h"

#include <omp.h>

using namespace std;
using namespace chrono;


// Input: I J K fidelity filename initial_conf (optional) final_conf (optional)
int main(int argc, char **argv) {

  // Initialize TALSH
  talsh::initialize();
  {
    int errc;
    size_t super_dim = (size_t)pow(2,3);

    // Allocate talsh::Tensors involved in the contraction
    vector<int> dims_2(2, super_dim);
    vector<int> dims_3(3, super_dim);
    vector<int> dims_4(4, super_dim);
    vector<int> dims_5(5, super_dim);
    vector<int> dims_6(6, super_dim);
    vector<int> dims_7(7, super_dim);
    size_t vol_2 = (size_t)pow(super_dim,2);
    size_t vol_3 = (size_t)pow(super_dim,3);
    size_t vol_4 = (size_t)pow(super_dim,4);
    size_t vol_5 = (size_t)pow(super_dim,5);
    size_t vol_6 = (size_t)pow(super_dim,6);
    size_t vol_7 = (size_t)pow(super_dim,7);
    s_type * data_L = (s_type *) malloc(vol_3*sizeof(s_type));
    errc = talsh::pinHostMemory(data_L,vol_3*sizeof(s_type));
    assert(errc==0);
    s_type * data_R = (s_type *) malloc(vol_3*sizeof(s_type));
    errc = talsh::pinHostMemory(data_R,vol_3*sizeof(s_type));
    assert(errc==0);
    s_type * data_D = (s_type *) malloc(vol_2*sizeof(s_type));
    errc = talsh::pinHostMemory(data_D,vol_2*sizeof(s_type));
    assert(errc==0);

    // Set values
    for (int p=0; p<vol_3; ++p)
    {
      *(data_L+p) = 1.0;
      *(data_R+p) = 2.0;
    }

    talsh::Tensor L(dims_3, data_L);
    talsh::Tensor R(dims_3, data_R);
    talsh::Tensor D(dims_2, data_D);
    
    TensContraction tc("D(b,d)+=L(a,b,c)*R(c,a,d)",
                        &D, &L, &R);
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_NVIDIA_GPU,0));

    for (int p=0; p<vol_2; ++p)
    {
      cout << *(data_D+p) << " ";
    }

    errc = talsh::unpinHostMemory(data_L); assert(errc == 0);
    delete [] data_L; data_L = nullptr;
    errc = talsh::unpinHostMemory(data_R); assert(errc == 0);
    delete [] data_R; data_R = nullptr;
    errc = talsh::unpinHostMemory(data_D); assert(errc == 0);
    delete [] data_D; data_D = nullptr;
  }
  // Shut down TALSH
  talsh::shutdown();

  return 0;
} 
