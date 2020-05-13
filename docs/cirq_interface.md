# The Cirq Interface

HOWTO use qFlex from Google [Cirq](https://github.com/quantumlib/cirq).
* Note: For the moment, the Cirq-qFlex integration is experimental code.

The currently supported version of Cirq is 0.8.0, and newer stable
versions will be supported, once they are available. Therefore, the
prerequisites for running qFlex via Cirq are:
- Cirq 0.8.0 or later
- qFlex C++ code (the code from this repo)

This file is an example of how to use the qFlex Python interface with Cirq.


## Setting up

1. Compile qFlex with Pybind support (see [install.md](/docs/install.md)).
The `./python` folder
should contain a shared library that is available as module to Python.

2. Cirq can either be [installed on your machine](https://cirq.readthedocs.io/en/stable/install.html),
or you can install Cirq in a virtual environment by calling
`./scripts/create_cirq_environment.sh`

3. If the environment existed, make sure to `source .venv/bin/activate`
* Note `deactivate` exits the environment

4. Ensure that Python modules load correctly by modifying PYTHONPATH:
`export PYTHONPATH=$PYTHONPATH:.`
* Note: All the paths in this document, assume the root is the directory where the
repository was cloned. For example, `/home/user/github/qflex`.

5. Run the example about how to simulate a small and a large quantum circuit.
`python qflexcirq/qflex_cirq_example.py`


## Interface design and operations

The goal is to simulate circuits, and there are two major options:
* loading qFlex text files into Cirq objects
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
QFlexOrdering is used to manually load, or to automatically generate orderings
by a heuristic (see Notes).

### Usage procedure

* Create a grid of qubits starting from a string like (see [input_formats.md](/docs/input_formats.md))
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

A QFlexCircuit can be created from a Cirq circuit. However, if a cirq.Circuit
is not available, an example file from the qFlex (e.g.
[a circuit for Bristlecone](/config/circuits/bristlecone_48_1-16-1_0.txt))
can be loaded:
```
my_circuit = qflexutils.GetCircuitOfMoments(config["circuit_filename"], my_qubits)
```

At this point, `my_circuit` is a cirq.Circuit, and can be used to construct a
QFlexCircuit
```
my_qflex_circuit = qcirc.QFlexCircuit(cirq_circuit = my_circuit, device = my_device)
```

The tensor contraction order was not specified, and the internal heuristic
will be used. It is possible to specify the ordering (see Notes), and to also specify if
gate decompositions should be performed or not (see Notes).

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

### Gate decompositions
* TODO: Enable the full functionality.

The QFlexCircuit is capable of decomposing arbitrary Cirq gates to the
elementary gate set of qFlex. However, if the user does not specify valid
decompositions, the QFlexCircuit composition will raise exceptions.

The constructor takes the argument `allow_decomposition` which is `False` by
default. If set `True`, the gates from the original circuit will decomposed,
using their ` _decompose_()` contract. For example, the [cirq.CNOT is decomposed
into](https://github.com/quantumlib/Cirq/blob/49b2f193ad99ce6770831330c19963bfa5c66f19/cirq/ops/common_gates.py#L829):
```
yield YPowGate(exponent=-0.5).on(t)
yield CZ(c, t)**self._exponent
yield YPowGate(exponent=0.5).on(t)
```

qFlex accepts the `CZ` and the `YPowGate(exponent=0.5)`, but
`YPowGate(exponent=-0.5)` is not supported for the moment,
and an exception will be raised.

### Tensor contractions: automatic and manual

The QFlexCircuit constructor accepts a QFlexOrder object specified as
`qflex_order`. If the param is None, the circuit is used to automatically compute
an ordering.

Behind the scenes, the QFlexOrdering can be constructed either:
* automatically by using `circuit_to_ordering()` from
[order_circuit_simulation.py](python/ordering/ordering/order_circuit_simulation.py)

* manual loading from a file
```
my_order = qorder.QFlexOrder.from_existing_file(config["ordering_filename"])
```

### Parametrized circuits
* TODO: Enable the full functionality.

The QFlexSimulator accepts parametrized circuits and these are resolved, but
for the moment support is available for the Rz and FSim gates. 

he QFlexVirtualDevice does not validate parametrized gates except fSim and Rz, 
and circuit composition will fail when such gates are added.

If needed, file an issue.

### Use qFlex from Python without Cirq

This is possible by using the Pybind11 interface, which, currently, can be used
exclusively with files. Sample files are in the `config` directory. An
example of how to use this interface is the `run_pybind_interface()` method from
`./python/qflex_cirq_example.py`

A simulator and a device are necessary:

* The Device is a virtual device which ensures that circuit gates are decomposed
into qFlex-supported gate types. For the moment this is a restricted gate set;
for a full listing of valid gates see [QflexVirtualDevice](/qflexcirq/interface/qflex_virtual_device.py).

* The Simulator calls the Python pybind11-based interface and checks that the
circuit to simulate was compiled (validated) against the device.
