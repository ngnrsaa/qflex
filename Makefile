TARGET1 = qflex

CXX = g++
#CXX = icpc

FLAGS =  -O3  -std=c++17  -march=native -Idocopt.cpp/

ifeq ($(CXX), icpc)
	FLAGS += -mkl -qopenmp -DMKL_TENSOR
else
	# see instructions at https://pybind11.readthedocs.io/en/stable/basics.html
	# python3-dev is required for pybind11 to work
  FLAGS += -Ipybind11/include -fPIC `python3 -m pybind11 --includes` -fopenmp -lgsl -lgslcblas
endif

TEST_DIR = tests

OBJS1 = main.o evaluate_circuit.o tensor.o contraction_utils.o read_circuit.o docopt.o

$(TARGET1): $(OBJS1)
	$(CXX) -o $(TARGET1).x $(OBJS1) $(FLAGS)

pybind: pybind_main.o $(OBJS1)
	$(CXX) -shared -o ./cirqinterface/$(TARGET1)`python3-config --extension-suffix` pybind_main.o $(OBJS1) $(FLAGS)

pybind_main.o: pybind_main.cpp
	$(CXX) -c pybind_main.cpp $(FLAGS)

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

docopt.o:
	$(CXX) -c docopt.cpp/docopt.cpp $(FLAGS)


.PHONY: clean
clean:
	rm -f ./*.x ./*.a ./*.so ./*.o ./*.mod
	$(MAKE) -C $(TEST_DIR) clean


# c++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` example.cpp -o example`python3-config --extension-suffix`
