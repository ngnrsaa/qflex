"""Tests for simulator.py

Invocation:
  $ pytest simulator_test.py
"""

from python.circuits import simulator
import numpy as np
import cirq

grid_filename = 'config/grid/rectangular_2x2.txt'
circuit_filename = 'config/circuits/rectangular_2x2.txt'
final_state = np.array([(-0.051776696+0.17677666j),
                       (-0.051776696-0.17677666j),
                       (0.124999985+0.24999996j),
                       (0.124999985-0.24999996j),
                       0.30177665j,
                       0.30177665j,
                       (0.17677668-0.12499996j),
                       (-0.17677668-0.12499996j),
                       (0.12499996+0.24999993j),
                       (0.12499996-0.24999993j),
                       (0.30177665-0.17677666j),
                       (0.30177665+0.17677666j),
                       (-0.17677668-0.124999985j),
                       (0.17677668-0.124999985j),
                       -0.051776696j,
                       -0.051776696j])

def test_simulation():
  qubits = simulator.GetGridQubit(grid_filename)
  circuit = simulator.GetCircuit(circuit_filename, qubits)
  results = cirq.Simulator().simulate(circuit, qubit_order=qubits)
  assert np.all(np.abs(results.final_state - final_state) < 1e-8)

