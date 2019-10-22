# The Cirq Interface


QFlex can be used from Google Cirq. This is the HOWTO.
* Note: For the moment, this the Cirq-QFlex is experimental code.

The currently supported version of Cirq is 0.5.0, and newer stable 
versions will be supported, once they are available. Therefore, the 
prerequisites are:
- Cirq 0.5.0
- QFlex C++ code (the code from this repo)



This file is an example of how to use the QFlex Python interface with Cirq.

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