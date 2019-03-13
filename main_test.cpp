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
    size_t super_dim = (size_t)pow(2,3);

    // Allocate talsh::Tensors involved in the contraction
    vector<int> dims_2(2, super_dim);
    vector<int> dims_3(3, super_dim);
    size_t vol_2 = (size_t)pow(super_dim,2);
    size_t vol_3 = (size_t)pow(super_dim,3);
    vector<s_type> data_L(vol_3, 1.0);
    vector<s_type> data_R(vol_3, 1.0);
    vector<s_type> data_D(vol_2, 1.0);

    
    talsh::Tensor L(dims_3, data_L);
    talsh::Tensor R(dims_3, data_R);
    talsh::Tensor D(dims_2, s_type(1.0));
    
    // First product
    TensContraction tc("D(b,d)+=L(a,b,c)*R(c,a,d)", &D, &L, &R);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    s_type const * ptr_D;
    D.getDataAccessHostConst(&ptr_D);
    for (int p=0; p<D.getVolume(); ++p)
    {
      cout << ptr_D[p] << endl;
    }

    // Second product
    tc = TensContraction("D(b,d)+=L(a,b,c)*R(c,a,d)", &D, &L, &R);
    errc = tc.execute(DEV_HOST,0); assert(errc==TALSH_SUCCESS);
    assert(tc.sync(DEV_HOST,0));
    D.getDataAccessHostConst(&ptr_D);
    for (int p=0; p<D.getVolume(); ++p)
    {
      cout << ptr_D[p] << endl;
    }
    

  }
  // Shut down TALSH
  talsh::shutdown();

  return 0;
} 
