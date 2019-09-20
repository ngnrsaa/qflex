TARGET1 = qflex

#Set these:
CXX = g++

FLAGS =  -O3  -std=c++17  -march=native

ifeq ($(CXX), icpc)
	FLAGS += -mkl -qopenmp -DMKL_TENSOR
else
  FLAGS += -fopenmp -lgsl -lgslcblas
endif

TEST_DIR = tests

OBJS1 = main.o evaluate_circuit.o tensor.o contraction_utils.o read_circuit.o

$(TARGET1): $(OBJS1)
	$(CXX) -o $(TARGET1).x $(OBJS1) $(FLAGS)

main.o: main.cpp
	$(CXX) -c main.cpp $(FLAGS)

evaluate_circuit.o: evaluate_circuit.cpp
	$(CXX) -c evaluate_circuit.cpp $(FLAGS)

read_circuit.o: read_circuit.cpp
	$(CXX) -c read_circuit.cpp $(FLAGS)

contraction_utils.o: contraction_utils.cpp
	$(CXX) -c contraction_utils.cpp $(FLAGS)

tensor.o: tensor.cpp
	$(CXX) -c tensor.cpp $(FLAGS)


.PHONY: clean
clean:
	rm -f ./*.x ./*.a ./*.so ./*.o ./*.mod
	$(MAKE) -C $(TEST_DIR) clean
