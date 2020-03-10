# Lint as: python3
"""Tests for google3.experimental.users.martinop.git.qflex.python.order_circuit_simulation."""

import random
import cirq
import pytest

import sys
sys.path.insert(1, '../../')

from qflexcirq.ordering import order_circuit_simulation as order_lib


def test_circuit_to_ordering():
    """Tests the circuit-to-ordering conversion.

  This is mostly a change-detecting test, as the heuristic is known to be
  non-optimal. If this test fails, verify that the new output is an improvement
  over the old output (either in time or memory cost) and update the test.
  """
    qubits = []
    size = 4
    for x in range(size):
        qubits.append([cirq.GridQubit(x, y) for y in range(size)])

    random.seed(0)

    # Build moments with ISWAPs and CZs on each adjacent pair of qubits.
    # Each bond has a pseudorandom dimension between 2 and 5:
    # (Dots represent qubits, numbers are bond dimension, and Xs are cuts)
    #
    #  .3.4.5.         Canonical qubit numbering goes in row-major order:
    #  5 5 2 4           - The top-left qubit, (0, 0), is index 0.
    #  .2.2.2.           - The top-right qubit, (0, 3), is index 3.
    #  X 5 4 5           - The bottom-left qubit, (3, 0), is index 12.
    #  .4.2.4.
    #  4 3 3 4
    #  .3.X.5.
    #
    # ISWAPs add 2 to the bond dimension, and CZs add 1.
    moments = ()
    for x in range(size - 1):
        for y in range(size):
            r = random.randint(2, 5)
            while r > 1:
                moments += (cirq.Moment(
                    [cirq.ISWAP(qubits[x][y], qubits[x + 1][y])]),)
                r -= 2
            while r > 0:
                moments += (cirq.Moment(
                    [cirq.CZ(qubits[x][y], qubits[x + 1][y])]),)
                r -= 1
    for y in range(size - 1):
        for x in range(size):
            r = random.randint(2, 5)
            while r > 1:
                moments += (cirq.Moment(
                    [cirq.ISWAP(qubits[x][y], qubits[x][y + 1])]),)
                r -= 2
            while r > 0:
                moments += (cirq.Moment(
                    [cirq.CZ(qubits[x][y], qubits[x][y + 1])]),)
                r -= 1

    circuit = cirq.Circuit(moments)
    order = order_lib.circuit_to_ordering(circuit=circuit, max_cuts=2)
    # Log the ordering file for debugging purposes.
    print("\n".join(order))
    assert len(order) == 39
    # The order of these two operations doesn't matter, and cuts may have their
    # indices Iswapped.
    # cut bond (cirq.GridQubit(3, 1), cirq.GridQubit(3, 2)) of dim 16
    cut1 = set(["cut () 13 14", "cut () 14 13"])
    # cut bond (cirq.GridQubit(1, 0), cirq.GridQubit(2, 0)) of dim 4
    cut2 = set(["cut () 9 10", "cut () 10 9"])
    # approximate fidelity: 1.0
    assert cut1.intersection(set([order[1], order[3]]))
    assert cut2.intersection(set([order[1], order[3]]))
    # These are the lines with commands.
    command_indices = [
        6, 7, 9, 10, 12, 13, 15, 17, 19, 21, 23, 25, 27, 29, 30, 32
    ]
    order5_33 = [order[i] for i in command_indices]
    assert order5_33 == [
        # 1: ['ISWAP((3, 2), (3, 3))', 'ISWAP((3, 2), (3, 3))', 'CZ((3, 2), (3, 3))']
        "expand A 14",
        "expand A 15",
        # 2: ['ISWAP((0, 0), (1, 0))', 'ISWAP((0, 0), (1, 0))', 'CZ((0, 0), (1, 0))']
        "expand B 0",
        "expand B 4",
        # 3: ['ISWAP((0, 2), (0, 3))', 'ISWAP((0, 2), (0, 3))', 'CZ((0, 2), (0, 3))']
        "expand C 2",
        "expand C 3",
        # 4: ['ISWAP((2, 3), (3, 3))', 'ISWAP((2, 3), (3, 3))']
        "expand A 11",
        # 5: ['ISWAP((2, 2), (2, 3))', 'ISWAP((2, 2), (2, 3))', 'ISWAP((2, 2), (3, 2))', 'CZ((2, 2), (3, 2))']
        "expand A 10",
        # 6: ['ISWAP((1, 3), (2, 3))', 'ISWAP((1, 3), (2, 3))', 'CZ((1, 3), (2, 3))']
        "expand A 7",
        # 7: ['ISWAP((1, 2), (2, 2))', 'ISWAP((1, 2), (2, 2))', 'ISWAP((1, 2), (1, 3))']
        "expand A 6",
        # 8: ['ISWAP((0, 2), (1, 2))', 'ISWAP((0, 3), (1, 3))', 'ISWAP((0, 3), (1, 3))']
        "merge A C",
        # 9: ['ISWAP((0, 1), (0, 2))', 'ISWAP((0, 1), (0, 2))']
        "expand C 1",
        # 10: ['ISWAP((0, 1), (1, 1))', 'ISWAP((0, 1), (1, 1))', 'CZ((0, 1), (1, 1))', 'ISWAP((1, 1), (1, 2))']
        "expand C 5",
        # 11: ['ISWAP((2, 0), (3, 0))', 'ISWAP((2, 0), (3, 0))']
        "expand D 8",
        "expand D 12",
        # 12: ['ISWAP((0, 0), (0, 1))', 'CZ((0, 0), (0, 1))', 'ISWAP((1, 0), (1, 1))']
        "merge B C"
    ]
    # The order of these two operations doesn't matter.
    # 13: ['ISWAP((1, 0), (2, 0))', 'ISWAP((1, 0), (2, 0))', 'CZ((1, 0), (2, 0))']
    # 14: ['ISWAP((1, 1), (2, 1))', 'ISWAP((1, 1), (2, 1))', 'CZ((1, 1), (2, 1))', 'ISWAP((2, 0), (2, 1))', 'ISWAP((2, 0), (2, 1))']
    assert (["expand C 9", "merge C D"] == [order[34], order[36]] or
            ["merge C D", "expand D 9"] == [order[34], order[36]])
    # 15: ['ISWAP((2, 1), (3, 1))', 'CZ((2, 1), (3, 1))', 'ISWAP((3, 0), (3, 1))', 'CZ((3, 0), (3, 1))']
    assert order[38] == "expand D 13"


def test_max_cuts_negative_fails():
    """Tests the circuit-to-ordering conversion."""
    qubits = [cirq.GridQubit(0, 0), cirq.GridQubit(0, 1)]
    size = 4
    for x in range(size):
        qubits.append([cirq.GridQubit(x, y) for y in range(size)])

    moments = (cirq.Moment([cirq.CZ(qubits[0], qubits[1])]),
               cirq.Moment([cirq.CZ(qubits[0], qubits[1])]))
    circuit = cirq.Circuit(moments)

    # max_cuts cannot be less than zero.
    with pytest.raises(ValueError):
        order_lib.circuit_to_ordering(circuit=circuit, max_cuts=-1)


def test_max_cuts_zero_succeeds():
    """Tests the circuit-to-ordering conversion."""
    qubits = [cirq.GridQubit(0, 0), cirq.GridQubit(0, 1)]
    size = 4
    for x in range(size):
        qubits.append([cirq.GridQubit(x, y) for y in range(size)])

    moments = (cirq.Moment([cirq.CZ(qubits[0], qubits[1])]),
               cirq.Moment([cirq.CZ(qubits[0], qubits[1])]))
    circuit = cirq.Circuit(moments)

    # max_cuts of zero will generate an ordering.
    order = order_lib.circuit_to_ordering(circuit=circuit, max_cuts=0)
    assert len(order) > 1


def test_match_fidelity():
    """Tests the fidelity-matching method."""
    qubits = [
        cirq.GridQubit(0, 0),
        cirq.GridQubit(0, 1),
        cirq.GridQubit(1, 1),
        cirq.GridQubit(1, 0)
    ]
    cuts = [(frozenset([qubits[0], qubits[1]]), 4),
            (frozenset([qubits[1], qubits[2]]), 8),
            (frozenset([qubits[2], qubits[3]]), 16),
            (frozenset([qubits[3], qubits[0]]), 32)]

    # Check that cut dimension affects resulting fidelity.
    target_fidelity = 17.0 / 32.0
    for cut in cuts:
        dim = cut[1]
        fidelity, unused_data = order_lib.match_fidelity(
            target_fidelity, dict([cut]))
        assert fidelity == float(dim / 2 + 1) / dim

    # Check that different sets of cuts can have different result fidelities,
    # even if the total cut dimension is the same.
    target_fidelity = (3.0 / 4.0) * (19.0 / 32.0)
    fidelity, unused_data = order_lib.match_fidelity(target_fidelity,
                                                     dict([cuts[0], cuts[3]]))
    assert fidelity == target_fidelity
    fidelity, unused_data = order_lib.match_fidelity(target_fidelity,
                                                     dict([cuts[1], cuts[2]]))
    # This is the closest possible fidelity with dimensions (8, 16).
    assert fidelity == (6.0 / 8.0) * (10.0 / 16.0)

    # Check that cut order does not affect fidelity.
    target_fidelity = (3.0 / 4.0) * (5.0 / 8.0) * (11.0 / 16.0) * (23.0 / 32.0)
    fidelity, unused_data = order_lib.match_fidelity(target_fidelity,
                                                     dict(cuts))
    assert fidelity == target_fidelity
    fidelity, unused_data = order_lib.match_fidelity(target_fidelity,
                                                     dict(cuts[::-1]))
    assert fidelity == target_fidelity
