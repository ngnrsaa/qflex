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

# The interface between Cirq and the Python interface to the C++ QFlex
import python.cirq_interface.qflex_simulator as qsim
import python.cirq_interface.qflex_virtual_device as qdevice
import python.cirq_interface.qflex_grid as qgrid
import python.cirq_interface.qflex_circuit as qcirc
import python.cirq_interface.qflex_order as qorder

import python.utils as qflexutils

my_sim = qsim.QFlexSimulator()

# mydevice = QFlexVirtualDevice(qflex_grid_string = QFlexGrid.create_rectangular(11, 12))
my_device = qdevice.QFlexVirtualDevice(qflex_grid_string = qgrid.QFlexGrid.BRISTLECONE70)

# The qubits are collected and indexed from the underlying grid_string
# that was passed as constructor to the Device
my_qubits = my_device.get_indexed_grid_qubits()

# Take a QFlex circuit and generate a Cirq circuit from it
# The Cirq circuit will be afterwards transformed into a Qflex circuit

# You can construct a Cirq circuit from an existing QFlex circuit
# Note that circuits provided in files were designed for a specific arrangement
qflex_file_path = "config/circuits/bristlecone_70_1-40-1_0.txt"
my_circuit = qflexutils.GetCircuitOfMoments(qflex_file_path, my_qubits)

# my_order = None
my_order = qorder.QFlexOrder.from_existing_file("config/ordering/bristlecone_70.txt")

circuit_on_device = qcirc.QFlexCircuit(cirq_circuit = my_circuit,
                                 device = my_device,
                                 qflex_order = my_order)

print("\n\nRunning QFlex simulation")

nr_q = 62
myres = my_sim.compute_amplitudes(program = circuit_on_device,
                                  bitstrings=["1" * nr_q, "0" * nr_q])

print("\n\nAfter QFlexCirq")
print(myres)