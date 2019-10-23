# The Cirq Interface


HOWTO use QFlex from Google Cirq.
* Note: For the moment, this the Cirq-QFlex is experimental code.

The currently supported version of Cirq is 0.5.0, and newer stable 
versions will be supported, once they are available. Therefore, the 
prerequisites are:
- Cirq 0.5.0
- QFlex C++ code (the code from this repo)

This file is an example of how to use the QFlex Python interface with Cirq.

## Setting up

1. Compile QFlex with Pybind support (see install.md). The `./python` folder 
should contain a shared library that is available as module to Python.

2. If Cirq is not installed globally on the machine, create a virtual environment
and install Cirq by calling `./scripts/construct_cirq_environment.sh`

3. If the environment existed, make sure to `source .venv/bin/activate`
* Note `deactivate` exits the environment

3. Run the example about how to simulate a small and a large quantum circuit.

`export PYTHONPATH=. && python python/qflex_cirq_example.py`

* Note: set PYTHONPATH before -- it will start the search of Python modules 
from the current working directory, aka root.
* Note: All the paths in this document, assume the root is the directory where the
repository was cloned. For example, `/home/user/github/qflex`.


## Interface design and operations

The goal is to simulate circuits, and there are two major options:
* loading QFlex text files into Cirq objects
* using native Cirq objects from the beginning

The above options can be hybridized. This means that, for example, one can load
a circuit from a file, but map the circuit to a grid of qubits which is
generated in the code (and not read from a file).

### Classes

The interface includes QFlexSimulator which is communicates through a Pybind11 
interface with QFlex. The simulator accepts only QFlexCircuits, which are 
effectively a specialized kind of `cirq.Circuit`. 

Architectural constraints such as permitting only grid qubits, QFlex supported 
gate sets, and different circuit validations are performed by the 
QFlexVirtualDevice.

A QFlexCircuit can work only with a QFlexVirtualDevice.

Finally, because the QFlex uses tensor contraction operations, the class
QFlexOrdering is used to load or to generate by a heuristic (see 
`circuit_to_ordering()` in 
`python/ordering/ordering/order_circuit_simulation.py`)


### Usage procedure

* Create a grid of qubits starting from a string like (see input_formats.md)
```
111
010
111
```
The string can be loaded from a file, or created (e.g. `create_rectangular()`).
For example, 
```
my_grid = qgrid.QFlexGrid.from_existing_file(grid_file_path)
```

The QflexGrid is useful for defining the QFlexVirtualDevice.
```
my_device = qdevice.QFlexVirtualDevice(qflex_grid = my_grid)
```

And the device can output a list of indexed qubits, which are useful later
during circuit compilation 
```
my_qubits = my_device.get_indexed_grid_qubits()
```

A QFlexCircuit can be created from a Cirq circuit. However, if that file is not
 available, and an example from the Qflex project has to be used First, from a Qflex file
```
my_circuit = qflexutils.GetCircuitOfMoments(config["circuit_filename"], my_qubits)
```

At this point, `my_circuit` is a cirq.Circuit, and can be used to construct a
QFlexCircuit
```
my_qflex_circuit = qcirc.QFlexCircuit(cirq_circuit = my_circuit, device = my_device)
```

The tensor contraction order was not specified, and the internal heuristic
will be used. It is possible to specify the ordering, and to also specify if
gate decompositions should be performed or not.

Finally, the QFlexCircuit can be simulated, and the `bitstrings` is a list of
strings representing the state vector for which the amplitudes are simulated.
In the following example, the simulation is performed for a single state
of 70 qubits being |1>.

```
my_sim = qsim.QFlexSimulator()
myres = my_sim.compute_amplitudes(program = my_qflex_circuit,
                                  bitstrings=["1" * 70])
```


## Notes

This version includes preliminary support for gate decompositions and 
parametrized operations. Users relying on them should know how Cirq works behind 
the scenes, in order to write the code that is still required.

### Parametrized circuits

The QFlexSimulator accepts parametrized circuits and these are resolved, but
for the moment this are not really supported. The QFlexVirtualDevice does not validate
parametrized gates, and circuit composition will fail when such gates are added.

If needed, file an issue.

### Gate decompositions

The QFlexCircuit is capable of decomposing arbitrary Cirq gates to the
elementary gate set of QFlex. However, if the user does not specify valid
decompositions, the QFlexCircuit composition will raise exceptions.

### Use QFlex from Python without Cirq

This is possible by using the Pybind11 interface, which, currently, can be used
exclusively with files. Sample files are in the `config` directory. An 
example of how to use this interface is the `run_pybind_interface()` method from
`./python/qflex_cirq_example.py`

A simulator and a device are necessary:

* The Device is a virtual device, that essentially makes sure that circuit gates
are decomposed (TODO) to supported gate types (for the moment,
 a restricted gate set)

* The Simulator calls the Python pybind11-based interface and checks that the
circuit to simulate was compiled (validated) against the device.
