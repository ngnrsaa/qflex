TARGET1 = qflex

#Set these:
CXX = g++

MKL_DIR = /home/alan

FLAGS = -fopenmp -O3 -std=c++17 -march=native \
-I${MKL_DIR}/intel/compilers_and_libraries_2019.3.199/linux/mkl/include \
-L${MKL_DIR}/intel/compilers_and_libraries_2019.3.199/linux/mkl/lib/intel64 \
-lmkl_intel_ilp64 -lmkl_core -lmkl_gnu_thread -Wno-narrowing

TEST_DIR = tests

OBJS1 = main.o evaluate_circuit.o mkl_tensor.o contraction_utils.o read_circuit.o

$(TARGET1): $(OBJS1)
	$(CXX) -o $(TARGET1).x $(FLAGS) $(OBJS1)

main.o: main.cpp
	$(CXX) -c main.cpp $(FLAGS)

evaluate_circuit.o: evaluate_circuit.cpp
	$(CXX) -c evaluate_circuit.cpp $(FLAGS)

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
