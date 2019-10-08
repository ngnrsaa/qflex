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
from qflex_virtual_device import QFlexVirtualDevice, _BRISTLECONE70
from qflex_simulator import QFlexSimulator

mysim = QFlexSimulator()

mydevice = QFlexVirtualDevice(arrangement=_BRISTLECONE70)

# This is a simple circuit consisting of a single GridQubit
# and a Hadamard gate
a = cirq.GridQubit(0, 5)
moment = cirq.Moment([cirq.H(a)])

# Run the simulation
mysim.simulate(cirq.Circuit(moments=[moment], device=mydevice),
               initial_state="YYYY")

