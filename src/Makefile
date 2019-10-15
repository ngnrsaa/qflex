TARGET1 = qflex

OBJS1 = evaluate_circuit.o tensor.o contraction_utils.o read_circuit.o docopt/docopt.o

FLAGS += -Idocopt/

$(TARGET1): main.cpp $(OBJS1)
	$(CXX) -o $(@).x $< $(OBJS1) $(FLAGS)

%.o: %.cpp %.h
	$(CXX) -c $< $(FLAGS) -o $@

.PHONY: clean
clean:
	-rm -f ./*.x ./*.a ./*.so ./*.mod $(OBJS1)
