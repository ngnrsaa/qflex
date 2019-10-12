TARGET1 = qflex

#Set these:
CXX = g++

FLAGS =  -O3  -std=c++17  -march=native

ifeq ($(CXX), icpc)
	FLAGS += -mkl -qopenmp -DMKL_TENSOR
else
  FLAGS += -fopenmp -lgsl -lgslcblas
endif

ROOT_DIR = ..
OBJS1 = evaluate_circuit.o tensor.o contraction_utils.o read_circuit.o docopt.cpp/docopt.o

GTEST_DIR = $(ROOT_DIR)/googletest/googletest
GMOCK_DIR = $(ROOT_DIR)/googletest/googlemock

TESTFLAGS = -I$(GTEST_DIR)/include -L$(GTEST_DIR) -I$(GMOCK_DIR)/include -L$(GMOCK_DIR) -I$(ROOT_DIR)/src

%.o:
	$(MAKE) -C $(ROOT_DIR) $@

gtest-all.o:
	$(MAKE) -C $(GTEST_DIR)/make gtest-all.o

gmock-all.o:
	$(MAKE) -C $(GMOCK_DIR)/make gmock-all.o

%.x: %.cpp $(OBJS1) gtest-all.o gmock-all.o 
	$(CXX) -o $@ $< $(OBJS1:%=../%) $(TESTFLAGS) $(FLAGS) $(GMOCK_DIR)/make/gmock-all.o $(GTEST_DIR)/make/gtest-all.o

.PHONY: clean
clean:
	rm -f ./*.x ./*.a ./*.so ./*.o ./*.mod
	$(MAKE) -C $(ROOT_DIR) clean
	$(MAKE) -C $(GMOCK_DIR)/make clean
	$(MAKE) -C $(GTEST_DIR)/make clean
