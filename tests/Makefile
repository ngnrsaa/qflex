TARGET1 = qflex

#Set these:
CXX = icpc

# Contraction can be bristlecone_70x4, bristlecone_48x4, ...
CONTRACTION = bristlecone_48x4
#Done.

# Contraction
CONTRACTION_FILENAME = main_$(CONTRACTION)

FLAGS =  -mkl  -qopenmp  -O3  -std=c++17  -march=native

ROOT_DIR = ..

GTEST_DIR = $(ROOT_DIR)/googletest/googletest
GMOCK_DIR = $(ROOT_DIR)/googletest/googlemock

TESTFLAGS = -I$(GTEST_DIR)/include -L$(GTEST_DIR) -I$(GMOCK_DIR)/include -L$(GMOCK_DIR)

all: read_circuit_test mkl_tensor_test

mkl_tensor:
	$(MAKE) -C $(ROOT_DIR) mkl_tensor.o

read_circuit:
	$(MAKE) -C $(ROOT_DIR) read_circuit.o

gtest-all.o:
	$(MAKE) -C $(GTEST_DIR)/make gtest-all.o

gmock-all.o:
	$(MAKE) -C $(GMOCK_DIR)/make gmock-all.o

read_circuit_test: read_circuit_test.cpp mkl_tensor read_circuit gtest-all.o gmock-all.o
	$(CXX) -o read_circuit_test.x read_circuit_test.cpp $(ROOT_DIR)/mkl_tensor.o \
	$(ROOT_DIR)/read_circuit.o $(GTEST_DIR)/make/gtest-all.o \
	$(GMOCK_DIR)/make/gmock-all.o $(FLAGS) $(TESTFLAGS)

mkl_tensor_test: mkl_tensor_test.cpp mkl_tensor gtest-all.o gmock-all.o
	$(CXX) -o mkl_tensor_test.x mkl_tensor_test.cpp $(ROOT_DIR)/mkl_tensor.o \
	$(GTEST_DIR)/make/gtest-all.o $(GMOCK_DIR)/make/gmock-all.o $(FLAGS) \
	$(TESTFLAGS)
.PHONY: clean
clean:
	rm -f ./*.x ./*.a ./*.so ./*.o ./*.mod
	$(MAKE) -C $(GTEST_DIR)/make clean
