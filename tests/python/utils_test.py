import pytest
import cirq

# Due to the directory structure of QFlex
# go up in hierarchy twice
import sys
sys.path.insert(1, '../../')

from qflexcirq import utils


def test_GetNrQubits():

    a = cirq.GridQubit(0, 0)
    b = cirq.GridQubit(0, 1)
    c = cirq.GridQubit(0, 2)

    circuit = cirq.Circuit()
    circuit.append(cirq.H.on(a))
    circuit.append(cirq.H.on(b))
    circuit.append(cirq.CNOT.on(b, c))
    circuit.append(cirq.H.on(a))

    nrq = utils.GetNumberOfQubits(circuit)

    # Cirq supported allready this functionality
    nrq2 = len(circuit.all_qubits())
    assert (nrq2 == 3)

    assert (nrq == 3)
