# The Cirq framework
import cirq

# The interface between Cirq and the Python interface to the C++ QFlex
from qflex_virtual_device import QFlexVirtualDevice, _BRISTLECONE70
from qflex_simulator import QFlexSimulator

mysim = QFlexSimulator()

mydevice = QFlexVirtualDevice(arrangement=_BRISTLECONE70)

a = cirq.GridQubit(0, 5)
moment = cirq.Moment([cirq.H(a)])

mysim.simulate(cirq.Circuit(moments=[moment], device=mydevice),
               initial_state="YYYY")

