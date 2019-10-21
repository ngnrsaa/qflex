#!/usr/bin/env python3

from subprocess import Popen, PIPE
from tempfile import mkstemp
from io import StringIO
from warnings import warn
import numpy as np
import pytest
import cirq

import sys, os
sys.path.insert(
    1,
    os.path.realpath(os.path.dirname(os.path.realpath(__file__)) + '/../../'))
from python.ordering import order_circuit_simulation as auto_order
from python import utils
from python import qflex

num_runs = 50

circuit_test = """16
0 h 0
0 h 1
0 h 2
0 h 3
0 h 4
0 h 5
0 h 6
0 h 7
0 h 8
0 h 9
0 h 10
0 h 11
0 h 12
0 h 13
0 h 14
0 h 15
1 cz 0 1
1 cz 6 7
1 cz 8 9
1 cz 14 15
1 t 2
1 t 3
1 t 4
1 t 5
1 t 10
1 t 11
1 t 12
1 t 13
2 cz 4 8
2 cz 6 10
2 y_1_2 0
2 y_1_2 1
2 x_1_2 7
2 x_1_2 9
2 y_1_2 14
2 x_1_2 15
3 cz 1 2
3 cz 9 10
3 t 0
3 y_1_2 4
3 x_1_2 6
3 t 7
3 x_1_2 8
3 t 14
3 t 15
4 cz 0 4
4 cz 9 13
4 cz 2 6
4 cz 11 15
4 y_1_2 1
4 t 8
4 y_1_2 10
5 cz 2 3
5 cz 4 5
5 cz 10 11
5 cz 12 13
5 y_1_2 0
5 t 1
5 x_1_2 6
5 y_1_2 9
5 y_1_2 15
6 cz 5 9
6 cz 7 11
6 t 0
6 x_1_2 2
6 y_1_2 3
6 y_1_2 4
6 t 6
6 y_1_2 10
6 y_1_2 12
6 x_1_2 13
6 t 15
7 cz 5 6
7 cz 13 14
7 t 2
7 t 3
7 t 4
7 y_1_2 7
7 y_1_2 9
7 t 10
7 y_1_2 11
7 t 12
8 cz 8 12
8 cz 1 5
8 cz 10 14
8 cz 3 7
8 x_1_2 6
8 t 9
8 t 11
8 x_1_2 13
9 cz 0 1
9 cz 6 7
9 cz 8 9
9 cz 14 15
9 x_1_2 3
9 y_1_2 5
9 y_1_2 10
9 y_1_2 12
9 t 13
10 cz 4 8
10 cz 6 10
10 x_1_2 0
10 y_1_2 1
10 t 3
10 t 5
10 x_1_2 7
10 y_1_2 9
10 t 12
10 y_1_2 14
10 x_1_2 15
11 cz 1 2
11 cz 9 10
11 t 0
11 y_1_2 4
11 x_1_2 6
11 t 7
11 y_1_2 8
11 t 14
11 t 15
12 cz 0 4
12 cz 9 13
12 cz 2 6
12 cz 11 15
12 y_1_2 1
12 t 8
12 x_1_2 10
13 cz 2 3
13 cz 4 5
13 cz 10 11
13 cz 12 13
13 x_1_2 0
13 t 1
13 y_1_2 6
13 x_1_2 9
13 x_1_2 15
14 cz 5 9
14 cz 7 11
14 t 0
14 y_1_2 2
14 x_1_2 3
14 y_1_2 4
14 t 6
14 x_1_2 10
14 y_1_2 12
14 y_1_2 13
14 t 15
15 cz 5 6
15 cz 13 14
15 t 2
15 t 3
15 t 4
15 x_1_2 7
15 x_1_2 9
15 t 10
15 x_1_2 11
15 t 12
16 cz 8 12
16 cz 1 5
16 cz 10 14
16 cz 3 7
16 y_1_2 6
16 t 9
16 t 11
16 y_1_2 13
17 h 0
17 h 1
17 h 2
17 h 3
17 h 4
17 h 5
17 h 6
17 h 7
17 h 8
17 h 9
17 h 10
17 h 11
17 h 12
17 h 13
17 h 14
17 h 15
"""

grid_test = """1111
1111
1111
1111
"""

ordering_test = """
cut () 4 8
cut () 5 9
cut () 6 10

expand A 0
expand A 1
expand A 2
expand A 3
expand A 4
expand A 5
expand A 6
expand A 7

expand B 8
expand B 9
expand B 10
expand B 11
expand B 12
expand B 13
expand B 14
expand B 15

merge A B
"""

ordering_with_cuts_test = """
cut () 4 8
cut () 5 9
cut () 6 10

expand A 0
expand A 1
expand A 2
expand A 3
expand A 4
expand A 5
expand A 6
expand A 7

expand B 8
expand B 9
expand B 10
expand B 11
expand B 12
expand B 13
expand B 14
expand B 15

merge A B
"""

# Simulate circuit
qubits = utils.GetGridQubits(StringIO(grid_test))
circuit = utils.GetCircuit(StringIO(circuit_test), qubits)
auto_ordering = auto_order.circuit_to_ordering(circuit,
                                               qubit_names=sorted(qubits))
results = cirq.Simulator().simulate(circuit)


@pytest.mark.parametrize(
    'x', [np.random.randint(0, 2**len(qubits)) for _ in range(num_runs)])
def test_simulation(x):

    # Get configuration as a string
    final_conf = bin(x)[2:].zfill(len(qubits))

    options = {
        'circuit': circuit_test.split('\n'),
        'ordering': ordering_test.split('\n'),
        'grid': grid_test.split('\n'),
        'final_state': final_conf
    }

    # Get output from qFlex
    qflex_amplitude = qflex.simulate(options)[0][1]

    # Compare the amplitudes
    assert (np.abs(results.final_state[x] - qflex_amplitude) < 1.e-6)


@pytest.mark.parametrize(
    'x', [np.random.randint(0, 2**len(qubits)) for _ in range(num_runs)])
def test_simulation_with_cuts(x):

    # Get configuration as a string
    final_conf = bin(x)[2:].zfill(len(qubits))

    options = {
        'circuit': circuit_test.split('\n'),
        'ordering': ordering_with_cuts_test.split('\n'),
        'grid': grid_test.split('\n'),
        'final_state': final_conf
    }

    # Get output from qFlex
    qflex_amplitude = qflex.simulate(options)[0][1]

    # Compare the amplitudes
    assert (np.abs(results.final_state[x] - qflex_amplitude) < 1.e-6)


@pytest.mark.parametrize(
    'x', [np.random.randint(0, 2**len(qubits)) for _ in range(num_runs)])
def test_simulation_with_auto_order(x):

    # Get configuration as a string
    final_conf = bin(x)[2:].zfill(len(qubits))

    options = {
        'circuit': circuit_test.split('\n'),
        'ordering': auto_ordering,
        'grid': grid_test.split('\n'),
        'final_state': final_conf
    }

    # Get output from qFlex
    qflex_amplitude = qflex.simulate(options)[0][1]

    # Compare the amplitudes
    assert (np.abs(results.final_state[x] - qflex_amplitude) < 1.e-6)
