TARGET1 = qflex

#Set these:
CXX = icpc

# Contraction can be bristlecone_70x4, bristlecone_48x4, ...
CONTRACTION = bristlecone_48x4
#Done.

# Contraction
CONTRACTION_FILENAME = main_$(CONTRACTION)

FLAGS =  -mkl  -qopenmp  -O3  -std=c++17  -march=native

TEST_DIR = tests

OBJS1 = $(CONTRACTION_FILENAME).o mkl_tensor.o contraction_utils.o read_circuit.o

$(TARGET1): $(OBJS1)
	$(CXX) -o $(TARGET1).x $(FLAGS) $(OBJS1)

$(CONTRACTION_FILENAME).o: $(CONTRACTION_FILENAME).cpp
	$(CXX) -c $(CONTRACTION_FILENAME).cpp $(FLAGS)

read_circuit.o: read_circuit.cpp
	$(CXX) -c read_circuit.cpp $(FLAGS)

contraction_utils.o: contraction_utils.cpp
	$(CXX) -c contraction_utils.cpp $(FLAGS)

mkl_tensor.o: mkl_tensor.cpp
	$(CXX) -c mkl_tensor.cpp $(FLAGS)


.PHONY: clean
clean:
	rm -f ./*.x ./*.a ./*.so ./*.o ./*.mod
	$(MAKE) -C $(TEST_DIR) clean
