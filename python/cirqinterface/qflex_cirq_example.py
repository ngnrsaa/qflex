"""
This file is an example of how to use the QFlex Python interface with Cirq
For this, a simulator and a device are necessary

The Device is a virtual device, that essentially makes sure that circuit gates
are decomposed to supported gate types (for the moment a restricted gate set
that was used for supremacy circuit simulations)

The Simulator calls the Python pybind11-based interface and checks that the
circuit to simulate was compiled (validated) against the device.

TODO: Improve this documentation

In order to call this example (on Linux machines)
* Create a venv in which Cirq is installed.
    ./construct_cirq_envinronment.sh
* Start the venv
    source source .venv/bin/activate
* Run the example
    python qflex_cirq_example.py
"""

# The Cirq framework
import cirq

# The interface between Cirq and the Python interface to the C++ QFlex
from cirqinterface.qflex_virtual_device import QFlexVirtualDevice
from cirqinterface.qflex_simulator import QFlexSimulator
from cirqinterface.qflex_grids import QFlexGrid

mysim = QFlexSimulator()

mydevice = QFlexVirtualDevice(arrangement=QFlexGrid.BRISTLECONE70)
# mydevice = QFlexVirtualDevice(arrangement=QFlexGrid.create_rectangular(11, 12))

# This is a simple circuit consisting of a single GridQubit
# and a Hadamard gate
# a = cirq.GridQubit(0, 5)
# moment = cirq.Moment([cirq.H(a)])

# Take a QFlex circuit and generate a Cirq circuit from it
# The Cirq circuit will be afterwards transformed into a Qflex circuit
from cirqinterface.transform_to_cirq_circuit import convert_qflex_circuit_file
mycirc = convert_qflex_circuit_file("circuits/bristlecone_70_1-40-1_0.txt")

circ_on_device = cirq.Circuit(mycirc._moments, device=mydevice)

print("\n\nRunning QFlex simulation")
# Run the simulation

q_off = QFlexGrid.get_qubits_off(QFlexGrid.BRISTLECONE70)
nr_q = len(mydevice._qubits) - len(q_off)
nr_q = 62

myres = mysim.compute_amplitudes(program=circ_on_device,
                                 bitstrings=["1" * nr_q, "0" * nr_q])

print("\n\nAfter QFlexCirq")
print(myres)

