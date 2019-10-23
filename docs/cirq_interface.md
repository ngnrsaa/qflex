# The Cirq Python interface


[Cirq](https://github.com/quantumlib/cirq) can be used with QFlex. The currently
supported version is Cirq 0.5.0, and new stable versions will be supported once
they are available.

Prerequisites are:
- Cirq 0.5.0
- QFlex

This file is an example of how to use the QFlex Python interface with Cirq.

A simulator and a device are necessary:

* The Device is a virtual device which ensures that circuit gates are decomposed
into qFlex-supported gate types. For the moment this is a restricted gate set;
for a full listing of valid gates see [QflexVirtualDevice](/python/cirq_interface/qflex_virtual_device.py).

* The Simulator calls the Python pybind11-based interface and checks that the
circuit to simulate was compiled (validated) against the device.

## Setting up

1. Compile QFlex with Pybind support and in the `./python` folder a shared library
will be available as module to Python.

2. Cirq can either be [installed on your machine](https://cirq.readthedocs.io/en/stable/install.html),
or you can install Cirq in a virtual environment by calling
`./scripts/construct_cirq_environment.sh`

3. Ensure that Python modules load correctly by modifying PYTHONPATH:
`export PYTHONPATH=$PYTHONPATH:.`

4. Run this example to verify that the interface is correctly enabled:
`python python/qflex_cirq_example.py`
