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


int main(int argc, char **argv) {

  // Initialize TALSH
  talsh::initialize();
  {
    int errc;
    size_t super_dim = (size_t)pow(2,4);

    // Allocate talsh::Tensors involved in the contraction
    vector<int> dims_2(2, super_dim);
    vector<int> dims_4(4, super_dim);
    vector<int> dims_7(7, super_dim);
    size_t vol_2 = (size_t)pow(super_dim,2);
    size_t vol_4 = (size_t)pow(super_dim,4);
    size_t vol_7 = (size_t)pow(super_dim,7);
    
    talsh::Tensor L(dims_7, s_type(0.0));
    talsh::Tensor M(dims_7, s_type(0.0));
    talsh::Tensor R(dims_4, s_type(0.0));
    talsh::Tensor D(dims_7, s_type(0.0));
    
    // First product
    TensContraction tc("D(a,b,c,d,e,h,i)+=L(a,b,c,d,e,f,g)*R(f,g,h,i)", &D, &L, &R);
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    s_type const * ptr_D;
    D.getDataAccessHostConst(&ptr_D);

    // Second product
    tc = TensContraction("D(a,b,c,d,e,h,i)+=L(a,b,c,d,e,f,g)*R(f,g,h,i)", &L, &D, &R);
    errc = tc.execute(DEV_NVIDIA_GPU,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    D.getDataAccessHostConst(&ptr_D);

  }
  // Shut down TALSH
  talsh::shutdown();

  return 0;
} 
