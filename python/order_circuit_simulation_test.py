# Lint as: python3
"""Tests for google3.experimental.users.martinop.git.qflex.python.order_circuit_simulation."""

import random
import cirq
import pytest

from python import order_circuit_simulation as order_lib


def test_circuit_to_ordering():
  """Tests the circuit-to-ordering conversion."""
  qubits = []
  size = 4
  for x in range(size):
    qubits.append([cirq.GridQubit(x, y) for y in range(size)])

  random.seed(0)

  # Build moments with CZs on each adjacent pair of qubits.
  # Each bond has a pseudorandom rank between 2 and 5:
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
  moments = ()
  for x in range(size - 1):
    for y in range(size):
      for _ in range(random.randint(2, 5)):
        moments += (cirq.Moment([cirq.CZ(qubits[x][y], qubits[x + 1][y])]),)
  for y in range(size - 1):
    for x in range(size):
      for _ in range(random.randint(2, 5)):
        moments += (cirq.Moment([cirq.CZ(qubits[x][y], qubits[x][y + 1])]),)

  circuit = cirq.Circuit(moments)
  order = order_lib.circuit_to_ordering(circuit)
  # Log the ordering file for debugging purposes.
  print("\n".join(order))
  assert len(order) == 36
  # The order of these two operations doesn't matter, and cuts may have their
  # indices swapped.
  cut1 = set(["cut () 13 14", "cut () 14 13"])
  cut2 = set(["cut () 9 10", "cut () 10 9"])
  assert cut1.intersection(set(order[0:2]))
  assert cut2.intersection(set(order[0:2]))
  # These are the lines with commands.
  command_indices = [3, 4, 6, 7, 9, 10, 12, 14, 16, 18, 20, 22, 24, 26, 27, 29]
  order2_30 = [order[i] for i in command_indices]
  assert order2_30 == [
      # 1: ['CZ((3, 2), (3, 3))', 'CZ((3, 2), (3, 3))', 'CZ((3, 2), (3, 3))',
      #     'CZ((3, 2), (3, 3))', 'CZ((3, 2), (3, 3))']
      "expand A 14",
      "expand A 15",
      # 2: ['CZ((0, 0), (1, 0))', 'CZ((0, 0), (1, 0))', 'CZ((0, 0), (1, 0))',
      #     'CZ((0, 0), (1, 0))', 'CZ((0, 0), (1, 0))']
      "expand B 0",
      "expand B 4",
      # 3: ['CZ((0, 2), (0, 3))', 'CZ((0, 2), (0, 3))', 'CZ((0, 2), (0, 3))',
      #     'CZ((0, 2), (0, 3))', 'CZ((0, 2), (0, 3))']
      "expand C 2",
      "expand C 3",
      # 4: ['CZ((2, 3), (3, 3))', 'CZ((2, 3), (3, 3))', 'CZ((2, 3), (3, 3))',
      #     'CZ((2, 3), (3, 3))']
      "expand A 11",
      # 5: ['CZ((2, 2), (2, 3))', 'CZ((2, 2), (2, 3))', 'CZ((2, 2), (2, 3))',
      #     'CZ((2, 2), (2, 3))', 'CZ((2, 2), (3, 2))', 'CZ((2, 2), (3, 2))',
      #     'CZ((2, 2), (3, 2))']
      "expand A 10",
      # 6: ['CZ((1, 3), (2, 3))', 'CZ((1, 3), (2, 3))', 'CZ((1, 3), (2, 3))',
      #     'CZ((1, 3), (2, 3))', 'CZ((1, 3), (2, 3))']
      "expand A 7",
      # 7: ['CZ((1, 2), (2, 2))', 'CZ((1, 2), (2, 2))', 'CZ((1, 2), (2, 2))',
      #     'CZ((1, 2), (2, 2))', 'CZ((1, 2), (1, 3))', 'CZ((1, 2), (1, 3))']
      "expand A 6",
      # 8: ['CZ((0, 2), (1, 2))', 'CZ((0, 2), (1, 2))', 'CZ((0, 3), (1, 3))',
      #     'CZ((0, 3), (1, 3))', 'CZ((0, 3), (1, 3))', 'CZ((0, 3), (1, 3))']
      "merge A C",
      # 9: ['CZ((0, 1), (0, 2))', 'CZ((0, 1), (0, 2))', 'CZ((0, 1), (0, 2))',
      #     'CZ((0, 1), (0, 2))']
      "expand C 1",
      # 10: ['CZ((0, 1), (1, 1))', 'CZ((0, 1), (1, 1))', 'CZ((0, 1), (1, 1))',
      #      'CZ((0, 1), (1, 1))', 'CZ((0, 1), (1, 1))', 'CZ((1, 1), (1, 2))',
      #      'CZ((1, 1), (1, 2))']
      "expand C 5",
      # 11: ['CZ((2, 0), (3, 0))', 'CZ((2, 0), (3, 0))', 'CZ((2, 0), (3, 0))',
      #      'CZ((2, 0), (3, 0))']
      "expand D 8",
      "expand D 12",
      # 12: ['CZ((1, 0), (1, 1))', 'CZ((1, 0), (1, 1))', 'CZ((0, 0), (0, 1))',
      #      'CZ((0, 0), (0, 1))', 'CZ((0, 0), (0, 1))']
      "merge B C"
  ]
  # The order of these two operations doesn't matter.
  assert (["expand C 9", "merge C D"] == [order[31], order[33]] or
          ["merge C D", "expand D 9"] == [order[31], order[33]])
  assert order[35] == "expand D 13"


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
