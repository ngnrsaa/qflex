TARGET1 = qflex

CXX = g++
#CXX = icpc

FLAGS =  -O3  -std=c++17  -march=native -I$(CURDIR)/docopt.cpp/

ifeq ($(CXX), icpc)
	FLAGS += -mkl -qopenmp -DMKL_TENSOR
else
  FLAGS += -fopenmp -lgsl -lgslcblas
endif

export CXX
export FLAGS

$(TARGET1):
	$(MAKE) -C src/
	ln -s src/$(TARGET1).x .

.PHONY: tests
tests:
	$(MAKE) -C tests/

.PHONY: run-tests
run-tests:
	$(MAKE) -C tests/ run-all

.PHONY: clean
clean:
	rm -f $(TARGET1).x
	$(MAKE) -C src/ clean
	$(MAKE) -C tests/ clean
