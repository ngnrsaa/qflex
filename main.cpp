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
    size_t vol_2 = (size_t)pow(super_dim,2);
    size_t vol_3 = (size_t)pow(super_dim,3);
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
    
    TensContraction tcg("D(b,d)+=L(a,b,c)*R(c,a,d)",
                        &D, &L, &R);
    errc = tcg.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tcg.sync(DEV_NVIDIA_GPU,0));

    for (size_t p=0; p<vol_2; ++p)
    {
      cout << *(data_D+p) << " ";
    }
    cout << "\n\n";

    TensContraction tc("D(b,d)+=L(a,b,c)*R(c,a,d)",
                        &D, &L, &R);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));

    for (size_t p=0; p<vol_2; ++p)
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
