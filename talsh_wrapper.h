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


/**
* @file talsh_wrapper.h
* Helper classes and functions for the use of TAL_SH in C++.
* @see https://github.com/ngnrsaa/qflex
*
* @author Dmitry Lyakh
* @author Benjamin Villalonga
* @date Created: February 2019
* @date Modified: February 2019
*/


#ifndef TALSH_WRAPPER_
#define TALSH_WRAPPER_

#include <iostream>
#include <memory>
#include <string>
#include <complex>
#include "talshxx.hpp"

using s_type = std::complex<float>;

//QC application tensor class:
class QCTensor{
  public:

    unsigned int getRank()
    {
      return static_cast<unsigned int>(shape_.size());
    }

    std::size_t getVolume()
    {
      std::size_t tvol = 1;
      for(const auto & dim: shape_) tvol*=static_cast<std::size_t>(dim);
      return tvol;
    }

    const std::vector<int> & getShape()
    {
      return shape_;
    }

    s_type * getDataPtr()
    {
      return tdata_;
    }

    QCTensor(const std::vector<int> & dims):
      shape_(dims)
  {
    std::size_t tvol = this->getVolume();
    tdata_ = new s_type[tvol];
    int errc = talsh::pinHostMemory(tdata_,tvol*sizeof(s_type)); assert(errc == 0);
  }

    QCTensor(const QCTensor & another) = delete;
    QCTensor & operator=(const QCTensor & another) = delete;

    QCTensor(QCTensor && another)
    {
      if(this != &another){
        this->shape_ = another.shape_;
        this->tdata_ = another.tdata_;
        another.shape_.clear();
        another.tdata_ = nullptr;
      }
    }

    QCTensor & operator=(QCTensor && another){
      if(this != &another){
        this->shape_ = another.shape_;
        this->tdata_ = another.tdata_;
        another.shape_.clear();
        another.tdata_ = nullptr;
      }
      return *this;
    }

    ~QCTensor()
    {
      if(tdata_ != nullptr){
        //std::cout << "Deleting tensor data " << (void*)tdata_ << std::endl; //debug
        int errc = talsh::unpinHostMemory(tdata_); assert(errc == 0);
        delete [] tdata_;
        tdata_ = nullptr;
        shape_.clear();
      }
    };

  private:

    std::vector<int> shape_;
    s_type * tdata_;
};


//TAL-SH tensor contraction specification class:
class TensContraction{
  public:

    TensContraction(const std::string & pattern,
                    talsh::Tensor * tens0,
                    talsh::Tensor * tens1,
                    talsh::Tensor * tens2,
                    s_type alpha = s_type{1.0f,0.0f}):
      index_pattern_(pattern),tensor0_(tens0),tensor1_(tens1),tensor2_(tens2),alpha_(alpha),task_hl_(new talsh::TensorTask())
    {
    }

    TensContraction(const TensContraction & another) = default;
    TensContraction & operator=(const TensContraction & another) = default;
    TensContraction(TensContraction && another) = default;
    TensContraction & operator=(TensContraction && another) = default;
    ~TensContraction() = default;


    //Schedules a tensor contraction for execution on the given device:
    int execute(int device_kind, int device_id)
    {
      int ierr = tensor0_->contractAccumulate(task_hl_.get(),index_pattern_,*tensor1_,*tensor2_,device_kind,device_id,alpha_);
      return ierr;
    }

    //Synchronizes tensor contraction execution and places the result to the given device:
    bool sync(const int dev_kind = DEV_HOST, //device kind on which the tensor-result should become available
              const int dev_id = 0,          //device id of a given kind on which the tensor-result should become available
              void * host_ptr = nullptr)     //external host memory pointer where to place the data from the tensor-result
    {
      bool done = tensor0_->sync(dev_kind,dev_id,host_ptr);
      return done;
    }

    //Checks for completion of a tensor contraction and places the result on the given device:
    bool ready(int * status,                  //status of the tensor contraction execution
               const int dev_kind = DEV_HOST, //device kind on which the tensor-result should become available
               const int dev_id = 0,          //device id of a given kind on which the tensor-result should become available
               void * host_ptr = nullptr)     //external host memory pointer where to place the data from the tensor-result
    {
      bool done = tensor0_->ready(status,dev_kind,dev_id,host_ptr);
      return done;
    }

  private:

    std::string index_pattern_;   //symbolic tensor contraction pattern (Einstein summation)
    talsh::Tensor * tensor0_;     //non-owning pointer to the tensor-result
    talsh::Tensor * tensor1_;     //non-owning pointer to the 1st input tensor
    talsh::Tensor * tensor2_;     //non-owning pointer to the 2nd input tensor
    s_type alpha_;           //optional scalar factor
    std::shared_ptr<talsh::TensorTask> task_hl_; //owning pointer to the TAL-SH task handle for this tensor contraction
};

#endif
