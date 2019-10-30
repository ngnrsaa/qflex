"""Tests for RQC generator.

Invocation:
    $ pytest generator_test.py
"""

from typing import Dict, Set, Tuple

import generator


def read_patterns(filename):
    with open(filename, 'r') as f:
        return eval(f.read())


TEST_PATTERNS_FILENAME = '../../config/patterns/test.txt'
RIGETTI_ASPEN_PATTERNS_FILENAME = '../../config/patterns/rigetti_aspen.txt'
IBM_ROCHESTER_PATTERNS_FILENAME = '../../config/patterns/ibm_rochester.txt'

TEST_DEVICE = generator.Device(read_patterns(TEST_PATTERNS_FILENAME))
RIGETTI_ASPEN = generator.Device(read_patterns(RIGETTI_ASPEN_PATTERNS_FILENAME))
IBM_ROCHESTER = generator.Device(read_patterns(IBM_ROCHESTER_PATTERNS_FILENAME))


def test_qubit_numbers():
    """Verify all devices have the correct number of qubits."""
    assert TEST_DEVICE.n_qubits == 4
    assert RIGETTI_ASPEN.n_qubits == 16
    assert IBM_ROCHESTER.n_qubits == 53


def test_coupler_numbers():
    """Verify all devices have the correct number of couplers."""
    assert TEST_DEVICE.n_couplers == 4
    assert RIGETTI_ASPEN.n_couplers == 18
    assert IBM_ROCHESTER.n_couplers == 58


def compute_edges(device):
    """Returns mapping from vertices to edges."""
    activation_patterns = {n: device.pattern(n) for n in device.pattern_names()}
    edges: Dict[int, Set[int]] = {}
    for couplers in activation_patterns.values():
        for q0, q1 in couplers:
            if q0 not in edges:
                edges[q0] = set()
            if q1 not in edges:
                edges[q1] = set()
            edges[q0].add(q1)
            edges[q1].add(q0)
    return edges


def degree(edges: Dict[int, Set[int]], vertex: int) -> int:
    """Returns degree of vertex."""
    return len(edges[vertex])


def test_graphs():
    """Verify all devices have the correct number of qubits with each degree."""
    test_edges = compute_edges(TEST_DEVICE)
    assert degree(test_edges, 0) == 2
    assert degree(test_edges, 1) == 2
    assert degree(test_edges, 2) == 2
    assert degree(test_edges, 3) == 2

    aspen_edges = compute_edges(RIGETTI_ASPEN)
    for qubit in range(16):
        if qubit in {1, 2, 13, 14}:
            assert degree(aspen_edges, qubit) == 3
        else:
            assert degree(aspen_edges, qubit) == 2

    rochester_edges = compute_edges(IBM_ROCHESTER)
    for qubit in range(53):
        if qubit in {51, 52}:
            assert degree(rochester_edges, qubit) == 1
        elif qubit in {9, 11, 13, 21, 23, 25, 32, 34, 36, 44, 46, 48}:
            assert degree(rochester_edges, qubit) == 3
        else:
            assert degree(rochester_edges, qubit) == 2


def test_qsim_output():
    """Verify qsim output."""
    circuit = generator.Circuit()
    circuit.append(
        generator.Cycle().append(generator.Gate(name='x_1_2', qubits=(
            0,))).append(generator.Gate(name='y_1_2', qubits=(1,))))
    circuit.append(
        generator.Cycle().append(generator.Gate(name='fs', qubits=(0, 1))))
    assert circuit.n_qubits == 2

    qsim = circuit.to_qsim_lines()
    assert len(qsim) == 3
    assert qsim[0] == '0 x_1_2 0'
    assert qsim[1] == '0 y_1_2 1'
    assert qsim[2] == '1 fs 0 1'
