# The Cirq Python interface

Google Cirq can be used with QFlex. The currently supported version is 
Cirq 0.5.0 and new stable versions will be supported, once they are available.

In order to use the Cirq interface to QFlex, the following are needed:
- Cirq 0.5.0
- QFlex

## Setting up

1. Compile QFlex with Pybind support and in the `./python` folder a shared library 
will be available as module to Python

2. If Cirq is not installed globally on the machine, create a virtual environment
and install Cirq by calling `./scripts/construct_cirq_environment.sh`

3. Run the example from `./python/cirq_if_example.py`