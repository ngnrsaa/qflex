"""Tests for simulator.py

Invocation:
  $ pytest simulator_test.py
"""

from python.circuits import simulator
import numpy as np
import cirq

grid_filename = 'config/grid/rectangular_2x2.txt'
circuit_filename = 'config/circuits/rectangular_2x2_1-2-1_0.txt'
final_state = [(-0.21338831-0.036611654j),
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
  qubits = simulator.GetGridQubit(grid_filename)
  circuit = simulator.GetCircuit(circuit_filename, qubits)
  results = cirq.Simulator().simulate(circuit, qubit_order=qubits)
  assert np.all(np.abs(results.final_state - final_state) < 1e-6)

