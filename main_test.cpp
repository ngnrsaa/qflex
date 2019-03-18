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

#include "mkl_tensor.h"

#include "talshxx.hpp"
#include "talsh_wrapper.h"

#include <omp.h>

using namespace std;
using namespace chrono;

// Input: I J K fidelity filename initial_conf (optional) final_conf (optional)
int main(int argc, char **argv) {


  s_type * scratch = new s_type[100];

  talsh::initialize();
  {

    vector<s_type> data_A({{0.8,0.1},{-0.3,0.5},{0.3,-0.9},{0.1,-0.1},
                           {-0.9,0.9},{-0.1,0.8},{0.7,0.7},{0.1,0.3}});
    vector<s_type> data_B({{0.1,0.2},{-0.2,-0.5},{-0.3,-0.7},{-0.7,-0.2},
                           {-0.2,-0.3},{-0.9,0.9},{-0.4,0.3},{0.2,0.2}});
    MKLTensor A({"a","b","c"}, {2,2,2}, data_A);
    MKLTensor B({"c","a","d"}, {2,2,2}, data_B);
    MKLTensor C({""}, {8});
    // MKLTensor
    //multiply(A, B, C, scratch);
    // MKLTensor up to here
    // TALSH
    C.set_indices_and_dimensions({"b","d"}, {2,2});
    talsh::Tensor Rsh(A.get_dimensions(), A.data());
    talsh::Tensor Ssh(B.get_dimensions(), B.data());
    talsh::Tensor Tsh(C.get_dimensions(), C.data());
    //TensContraction contraction("D(b,d)+=L(a,b,c)*R(c,a,d)", &Tsh, &Rsh, &Ssh);
    TensContraction contraction("D(d,b)+=L(c,b,a)*R(d,a,c)", &Tsh, &Rsh, &Ssh);
    int errc = contraction.execute(DEV_NVIDIA_GPU,0);
    assert(errc==TALSH_SUCCESS);
    assert(contraction.sync(DEV_NVIDIA_GPU,0));
    // TALSH up to here
    C.print();
    C.print_data();
  
    

  }
  talsh::shutdown();


  // Freeing scratch data: delete and NULL.
  delete[] scratch;
  scratch = NULL;

  return 0;
} 
