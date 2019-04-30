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


#ifndef CONTRACTION_H
#define CONTRACTION_H

#include <vector>
#include <string>
#include <cmath>
#include <complex>
#include <ctime>
#include <chrono>
#include <cassert>

#include <iostream>
#include <fstream>
#include <sstream>

#include <omp.h>

#include "talshxx.hpp"
#include "talsh_wrapper.h"


using namespace std;

/**
* Scalar type.
*/
typedef complex<float> s_type;

/**
* Represents an MKLTensor.
*/
class Contraction
{
  public:

    /**
    * Creates a Contraction, defining parameters and allocating all the helper
    * tensors.
    */
    Contraction(string input_string, int _num_args, int _num_amps);

    Contraction(const Contraction & another) = default;
    Contraction & operator=(const Contraction & another) = default;
    Contraction(Contraction && another) = default;
    Contraction & operator=(Contraction && another) = default;
    ~Contraction() = default;

    /**
    * Load circuit from file place it in open_tensor_grid. Only output open.
    **/
    void load_circuit(string input_string);

    /**
    * Fills amplitudes with the appropriate contractions.
    * Read from file, apply initial and final strings, contract for some paths.
    **/
    void contract(string input_string);

    /**
    * Get const reference to amplitudes.
    **/
    vector<s_type> const & get_amplitudes() const;

    /**
    * Get double with time of largest contraction.
    **/
    double get_time_largest_contraction() const;

  private:
    // Storage.

    // Parameters
    int I, J, K;
    double fidelity;
    string filename;
    size_t super_dim;
    int num_Cs;
    int num_qubits;
    vector<vector<int>> qubits_off, qubits_A;
    vector<s_type> amplitudes;
    int num_args, num_amps;

    // open_tensor_grid, universal to this object. Will be closed over and over
    // again for each amplitude computation.
    vector<vector<shared_ptr<talsh::Tensor>>> open_tensor_grid;

    // First, tensors for region C. Done by hand right now. Change in future.
    vector<shared_ptr<talsh::Tensor>> Cs; 

    // Second, helper tensors.
    shared_ptr<talsh::Tensor> H_2_legs_a, H_2_legs_b,
                              H_3_legs_a, H_3_legs_b,
                              H_4_legs_a, H_4_legs_b,
                              H_5_legs_a, H_5_legs_b, H_5_legs_c,
                              H_6_legs_a, H_6_legs_b, H_6_legs_c,
                              H_7_legs_a, H_7_legs_b,
                              H_8_legs_a, H_8_legs_b, H_8_legs_c, H_8_legs_d,
                              H_9_legs_a, H_9_legs_b,
                              H_10_legs_a, H_10_legs_b, H_10_legs_c;

    // Third, region (reused) tensors.
    shared_ptr<talsh::Tensor> AB; // It looks like we can get rid of this
                                    //one for low fidelity runs :)
    shared_ptr<talsh::Tensor> pE;

    // Tensors to hold the slices.
    shared_ptr<talsh::Tensor> S_4_10, S_5_10, S_10_4, S_10_5;
    // These ones loose the trivial leg.

    // Finally, scalar S
    shared_ptr<talsh::Tensor> S;

    // Slice configurations
    vector<vector<int>> slice_confs;

    // Time largest contraction
    double time_largest_contraction;


    // Normalization factor
    s_type norm_factor, global_norm_factor;

};


#endif
