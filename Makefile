TARGET1 = qflex

#Set these:
CXX = icpc
BLAS = MKL

# Contraction can be bristlecone_70x4, 7x7x5, ...
CONTRACTION = bristlecone_70x4
#Done.


# Contraction
CONTRACTION_FILENAME = main_$(CONTRACTION)

FLAGS = -mkl  -qopenmp  -O3  -std=c++17  -march=native

OBJS1 = $(CONTRACTION_FILENAME).o mkl_tensor.o

$(TARGET1): $(OBJS1)
	$(CXX) -o $(TARGET1).x $(FLAGS) $(OBJS1) $(TALSH_LIB) $(BLAS_LIB) $(CUDA_LIB) $(FORT_LIB) -O3

$(CONTRACTION_FILENAME).o: $(CONTRACTION_FILENAME).cpp
	$(CXX) -c scheduler.cpp  $(TALSH_INC) $(BLAS_INC) $(CUDA_INC) -fopenmp -O3 -std=c++11 -fPIC -D$(CONTRACTION)

mkl_tensor.o: mkl_tensor.cpp
	$(CXX) -c mkl_tensor.cpp $(FLAGS)
.PHONY: clean
clean:
	rm -f ./*.x ./*.a ./*.so ./*.o ./*.mod
