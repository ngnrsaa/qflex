TARGET1 = qflex

CXX = g++
#CXX = icpc

FLAGS =  -O3  -std=c++17  -march=native -Idocopt.cpp/

ifeq ($(CXX), icpc)
	FLAGS += -mkl -qopenmp -DMKL_TENSOR
else
  FLAGS += -fopenmp -lgsl -lgslcblas
endif

TEST_DIR = tests

OBJS1 = evaluate_circuit.o tensor.o contraction_utils.o read_circuit.o docopt.cpp/docopt.o

$(TARGET1): src/main.cpp $(OBJS1)
	$(CXX) -o $(@).x $< $(OBJS1) $(FLAGS)

%.o: src/%.cpp src/%.h
	$(CXX) -c $< $(FLAGS) -o $@

.PHONY: clean
clean:
	rm -f ./*.x ./*.a ./*.so ./*.o ./*.mod
