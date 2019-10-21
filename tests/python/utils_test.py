import pytest
import cirq

# Due to the directory structure of QFlex
# go up in hierarchy twice
import sys
sys.path.insert(1, '../../')

from python import utils


def test_GetNrQubits():

    a = cirq.GridQubit(0, 0)
    b = cirq.GridQubit(0, 1)
    c = cirq.GridQubit(0, 2)

    circ = cirq.Circuit()
    circ.append(cirq.H.on(a))
    circ.append(cirq.H.on(b))
    circ.append(cirq.CNOT.on(b, c))
    circ.append(cirq.H.on(a))

    nrq = utils.GetNrQubits(circ)

    assert(nrq == 3)