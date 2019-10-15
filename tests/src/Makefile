TARGETS = $(shell find . -name '*_test.cpp')
TARGETS := $(TARGETS:%.cpp=%.x)

OBJS1 = evaluate_circuit.o tensor.o contraction_utils.o read_circuit.o

SRC_DIR = ../../src
GTEST_DIR = $(CURDIR)/googletest/googletest
GMOCK_DIR = $(CURDIR)/googletest/googlemock

TESTFLAGS = -I$(GTEST_DIR)/include -L$(GTEST_DIR) -I$(GMOCK_DIR)/include -L$(GMOCK_DIR) -I$(SRC_DIR)

.PHONY: all
all: $(TARGETS)

.PHONE: run-all
run-all: all
	for exe in $(TARGETS); do ../$$exe; done

%.o:
	$(MAKE) -C $(SRC_DIR) $@

gtest-all.o:
	$(MAKE) -C $(GTEST_DIR)/make gtest-all.o

gmock-all.o:
	$(MAKE) -C $(GMOCK_DIR)/make gmock-all.o

%.x: %.cpp $(OBJS1) gtest-all.o gmock-all.o 
	$(CXX) -o ../$@ $< $(OBJS1:%=$(SRC_DIR)/%) $(TESTFLAGS) $(FLAGS) $(GMOCK_DIR)/make/gmock-all.o $(GTEST_DIR)/make/gtest-all.o

.PHONY: clean
clean:
	rm -f ./*.x ./*.a ./*.so ./*.mod $(OBJS1)
	$(MAKE) -C $(GMOCK_DIR)/make clean
	$(MAKE) -C $(GTEST_DIR)/make clean
