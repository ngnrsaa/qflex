TARGET1 = qflex

CXX = mpicxx

RT_LIB = -lgfortran

OLCF_XLF_ROOT = /sw/summit/xl/16.1.1-beta4/xlf/16.1.1
BLAS_LIB = -L$(OLCF_ESSL_ROOT)/lib64 -lessl -L$(OLCF_XLF_ROOT)/lib -lxlf90_r -lxlfmath
BLAS_INC = -D BLAS_ESSL

CUDA_LIB = -L$(OLCF_CUDA_ROOT)/lib64 -lcublas -lcudart -lnvToolsExt
CUDA_INC = -I$(OLCF_CUDA_ROOT)/include

TALSH_ROOT = /ccs/home/div/src/TAL_SH
TALSH_LIB = -L$(TALSH_ROOT) -ltalsh
TALSH_INC = -I$(TALSH_ROOT)


OBJS1 = main.o

$(TARGET1): $(OBJS1)
	$(CXX) -o $(TARGET1).x -fopenmp -fPIC $(OBJS1) $(TALSH_LIB) $(BLAS_LIB) $(CUDA_LIB) $(RT_LIB)

main.o: main.cpp
	$(CXX) -c main.cpp $(TALSH_INC) $(BLAS_INC) $(CUDA_INC) -fopenmp -O3 -std=c++11 -fPIC

.PHONY: clean
clean:
	rm -f ./*.x ./*.a ./*.so ./*.o ./*.mod
