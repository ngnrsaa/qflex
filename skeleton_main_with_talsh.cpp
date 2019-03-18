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

#include "talshxx.hpp"

#include <omp.h>

using namespace std;
using namespace chrono;

using s_type = std::complex<float>;

int main(int argc, char **argv) {

  cout << "Started execution...\n";
  talsh::initialize();
  {
    cout << "Initialized TALSH.\n";
    int dim = (int)pow(2,5);
    vector<int> signature({dim,dim,dim,dim,dim}); 
    s_type * data = new s_type[(int)pow(dim,6)];
    cout << "Allocated memory.\n";
    int errc = talsh::pinHostMemory(data,(int)pow(dim,6)*sizeof(s_type)); assert(errc==0);
    cout << "Pinned memory.\n";
    talsh::Tensor T(signature, data);
    cout << "Created Tensor.\n";
    T.print();
    errc = talsh::unpinHostMemory(data); assert(errc == 0);
    cout << "Unpinned memory.\n";
  }
  talsh::shutdown();
  cout << "...shutdown TALSH.\n";

  return 0;
} 
