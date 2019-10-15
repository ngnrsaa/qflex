TARGET1 = qflex

CXX = g++
#CXX = icpc

FLAGS =  -O3  -std=c++17  -march=native

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
	$(MAKE) -C tests/src/

.PHONY: run-tests
run-tests: tests
	$(MAKE) -C tests/src/ run-all

.PHONY: clean
clean:
	rm -f $(TARGET1).x
	$(MAKE) -C src/ clean
	$(MAKE) -C tests/src/ clean
