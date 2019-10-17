"""Tests for simulator.py

Invocation:
  $ pytest simulator_test.py
"""

from python.circuits import simulator
from io import StringIO
import numpy as np
import cirq

grid_test = """11
11
"""

circuit_test = """4
0 h 0
0 h 1
0 h 2
0 h 3
1 t 0
1 t 1
1 t 2
1 t 3
2 cx 0 1
2 cz 2 3
3 t 0
3 t 1
3 t 2
3 t 3
4 cz 0 2
4 cx 1 3
8 t 0
8 t 1
8 t 2
8 t 3
9 cz 0 1
9 cx 2 3
10 t 0
10 t 1
10 t 2
10 t 3
11 cx 0 2
11 cz 1 3
17 h 0
17 h 1
17 h 2
17 h 3
"""

final_state_test = [(-0.21338831-0.036611654j),
                    (-0.088388324-0.08838833j),
                    (-0.036611643-0.21338831j),
                    (0.08838834+0.088388324j),
                    (-0.08838833-0.088388324j),
                    (0.036611654+0.21338831j),
                    (0.088388324+0.08838834j),
                    (0.21338831+0.036611646j),
                    (0.03661166+0.21338834j),
                    (0.515165+0.015165047j),
                    (0.21338831+0.036611617j),
                    (-0.015165037-0.515165j),
                    (-0.08838833-0.088388324j),
                    (0.036611654+0.21338831j),
                    (0.088388324+0.08838834j),
                    (0.21338831+0.036611646j)]

def test_simulation():

  qubits = simulator.GetGridQubit(StringIO(grid_test))
  circuit = simulator.GetCircuit(StringIO(circuit_test), qubits)

  results = cirq.Simulator().simulate(circuit, qubit_order=qubits)
  assert np.all(np.abs(results.final_state - final_state_test) < 1e-6)

