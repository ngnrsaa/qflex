TARGET1 = qflex

#Set these:
CXX = icpc

FLAGS =  -mkl  -qopenmp  -O3  -std=c++17  -march=native

ROOT_DIR = ..

GTEST_DIR = $(ROOT_DIR)/googletest/googletest
GMOCK_DIR = $(ROOT_DIR)/googletest/googlemock

TESTFLAGS = -I$(GTEST_DIR)/include -L$(GTEST_DIR) -I$(GMOCK_DIR)/include -L$(GMOCK_DIR)

all: evaluate_circuit_test read_circuit_test contraction_utils_test mkl_tensor_test

evaluate_circuit:
	$(MAKE) -C $(ROOT_DIR) evaluate_circuit.o

mkl_tensor:
	$(MAKE) -C $(ROOT_DIR) mkl_tensor.o

contraction_utils:
	$(MAKE) -C $(ROOT_DIR) contraction_utils.o

read_circuit:
	$(MAKE) -C $(ROOT_DIR) read_circuit.o

gtest-all.o:
	$(MAKE) -C $(GTEST_DIR)/make gtest-all.o

gmock-all.o:
	$(MAKE) -C $(GMOCK_DIR)/make gmock-all.o

evaluate_circuit_test: evaluate_circuit_test.cpp evaluate_circuit mkl_tensor contraction_utils read_circuit gtest-all.o gmock-all.o
	$(CXX) -o evaluate_circuit_test.x evaluate_circuit_test.cpp \
	$(ROOT_DIR)/evaluate_circuit.o $(ROOT_DIR)/mkl_tensor.o \
	$(ROOT_DIR)/contraction_utils.o $(ROOT_DIR)/read_circuit.o \
	$(GTEST_DIR)/make/gtest-all.o $(GMOCK_DIR)/make/gmock-all.o \
	$(FLAGS) $(TESTFLAGS)

read_circuit_test: read_circuit_test.cpp mkl_tensor contraction_utils read_circuit gtest-all.o gmock-all.o
	$(CXX) -o read_circuit_test.x read_circuit_test.cpp $(ROOT_DIR)/mkl_tensor.o \
	$(ROOT_DIR)/read_circuit.o $(ROOT_DIR)/contraction_utils.o \
	$(GTEST_DIR)/make/gtest-all.o $(GMOCK_DIR)/make/gmock-all.o \
	$(FLAGS) $(TESTFLAGS)

contraction_utils_test: contraction_utils_test.cpp mkl_tensor contraction_utils gtest-all.o gmock-all.o
	$(CXX) -o contraction_utils_test.x contraction_utils_test.cpp \
	$(ROOT_DIR)/mkl_tensor.o $(ROOT_DIR)/contraction_utils.o \
	$(GTEST_DIR)/make/gtest-all.o $(GMOCK_DIR)/make/gmock-all.o \
	$(FLAGS) $(TESTFLAGS)

mkl_tensor_test: mkl_tensor_test.cpp mkl_tensor gtest-all.o gmock-all.o
	$(CXX) -o mkl_tensor_test.x mkl_tensor_test.cpp $(ROOT_DIR)/mkl_tensor.o \
	$(GTEST_DIR)/make/gtest-all.o $(GMOCK_DIR)/make/gmock-all.o $(FLAGS) \
	$(TESTFLAGS)
.PHONY: clean
clean:
	rm -f ./*.x ./*.a ./*.so ./*.o ./*.mod
	$(MAKE) -C $(GTEST_DIR)/make clean
