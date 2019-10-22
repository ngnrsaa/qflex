# The Cirq Interface


HOWTO use QFlex from Google Cirq.
* Note: For the moment, this the Cirq-QFlex is experimental code.

The currently supported version of Cirq is 0.5.0, and newer stable 
versions will be supported, once they are available. Therefore, the 
prerequisites are:
- Cirq 0.5.0
- QFlex C++ code (the code from this repo)

This file is an example of how to use the QFlex Python interface with Cirq.

The goal is to simulate circuits, and there are two major options:
* loading QFlex text files into Cirq objects
* native Cirq objects

The above options can be hybridized. This means that, for example, one can load
a circuit from a file, but map the circuit to a grid of qubits which is
generated in the code (and not read from a file).




A simulator and a device are necessary:

* The Device is a virtual device, that essentially makes sure that circuit gates
are decomposed (TODO) to supported gate types (for the moment,
 a restricted gate set)

* The Simulator calls the Python pybind11-based interface and checks that the
circuit to simulate was compiled (validated) against the device.

## Setting up

1. Compile QFlex with Pybind support and in the `./python` folder a shared library 
will be available as module to Python

2. If Cirq is not installed globally on the machine, create a virtual environment
and install Cirq by calling `./scripts/construct_cirq_environment.sh`

3. Run the example (Note: set PYTHONPATH before)

`export PYTHONPATH=. && python python/qflex_cirq_example.py`