import random
import cirq
import pytest

# Due to the directory structure of QFlex
# go up in hierarchy twice
import sys
sys.path.insert(1, '../../')

from qflexcirq import QFlexOrder
from qflexcirq.ordering import order_circuit_simulation as order_lib


def test_qorder_from_file():

    test_file_contents = "expand A 0\nmerge A B"

    import tempfile
    fd = tempfile.mkstemp()
    with open(fd[1], "w") as file:
        print(test_file_contents, file=file)

    qorder = QFlexOrder.from_existing_file(fd[1])

    read_contents = None
    with open(qorder.order_data, "r") as file:
        read_contents = "".join(file.readlines())

    # delete temp file
    import os
    os.unlink(fd[1])

    assert (read_contents.strip() == test_file_contents.strip())


def test_only_for_cirq_circuits():
    """
    If no strings are specified, the only other option is a Cirq circuit
    :return:
    """
    with pytest.raises(ValueError):
        assert QFlexOrder(qflex_order_strings=None, cirq_circuit="Should fail")


def test_cirquit_without_qubits():
    """
    If the Cirq circuit is specified a list of qubits is also necessary
    :return:
    """
    with pytest.raises(ValueError):
        assert QFlexOrder(qflex_order_strings=None,
                          cirq_circuit=cirq.Circuit.from_ops([]))


def test_automatic_ordering():

    test_size = 2
    test_circuit = generate_circuit(size=test_size)

    # This is a known order
    random.seed(0)
    order = order_lib.circuit_to_ordering(test_circuit)

    # Now I am creating a QFlexOrder
    qorder = QFlexOrder(qflex_order_strings=None,
                        cirq_circuit=test_circuit,
                        qubits=list(range(test_size**2)))

    read_contents = None
    with open(qorder.order_data, "r") as file:
        read_contents = "".join(file.readlines())

    assert (read_contents.strip() == "\n".join(order).strip())


def generate_circuit(size=4):
    qubits = []
    for x in range(size):
        qubits.append([cirq.GridQubit(x, y) for y in range(size)])

    moments = ()
    for x in range(size - 1):
        for y in range(size):
            for _ in range(random.randint(2, 5)):
                moments += (cirq.Moment(
                    [cirq.FSimGate(1, 2)(qubits[x][y], qubits[x + 1][y])]),)
    for y in range(size - 1):
        for x in range(size):
            for _ in range(random.randint(2, 5)):
                moments += (cirq.Moment(
                    [cirq.FSimGate(1, 2)(qubits[x][y], qubits[x][y + 1])]),)

    circuit = cirq.Circuit(moments)

    return circuit
