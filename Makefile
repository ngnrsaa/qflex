TARGET1 = qflex

#Set these:
CXX = mpicxx
BLAS = ESSL
FORT_LIB = -lgfortran
TALSH_ROOT = /ccs/home/villalonga/TAL_SH
MKL_ROOT = /home/div/intel
CUDA_ROOT = /sw/summit/cuda/9.2.148
#Done.


OLCF_ESSL_ROOT = /sw/summit/essl/6.1.0-2/essl/6.1
OLCF_XLF_ROOT = /sw/summit/xl/16.1.1-beta4/xlf/16.1.1
BLAS_LIB_ESSL = -L$(OLCF_ESSL_ROOT)/lib64 -lessl -L$(OLCF_XLF_ROOT)/lib -lxlf90_r -lxlfmath
BLAS_INC_ESSL = -DBLAS_ESSL

BLAS_LIB_MKL = -L$(MKL_ROOT)/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lpthread -lm -ldl
BLAS_INC_MKL = -I.

BLAS_LIB_ATLAS = -L/usr/lib/x86_64-linux-gnu -lblas
BLAS_INC_ATLAS = -I.

BLAS_LIB=$(BLAS_LIB_$(BLAS))
BLAS_INC=$(BLAS_INC_$(BLAS))

CUDA_LIB = -L$(CUDA_ROOT)/lib64 -lcublas -lcudart -lnvToolsExt
CUDA_INC = -I$(CUDA_ROOT)/include

TALSH_LIB = -L$(TALSH_ROOT) -ltalsh
#TALSH_LIB = $(TALSH_ROOT)/libtalsh.a
TALSH_INC = -I$(TALSH_ROOT)


OBJS1 = main.o

$(TARGET1): $(OBJS1)
	$(CXX) -o $(TARGET1).x -fopenmp -fPIC $(OBJS1) $(TALSH_LIB) $(BLAS_LIB) $(CUDA_LIB) $(FORT_LIB)

main.o: main.cpp
	$(CXX) -c main.cpp $(TALSH_INC) $(BLAS_INC) $(CUDA_INC) -fopenmp -O3 -std=c++11 -fPIC

.PHONY: clean
clean:
	rm -f ./*.x ./*.a ./*.so ./*.o ./*.mod
