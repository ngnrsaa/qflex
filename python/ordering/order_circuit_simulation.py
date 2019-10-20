# Lint as: python3
"""
Provides a method for converting Cirq circuits to qFlex simulation order.

Usage:
  order_circuit_simulation.py (-h | --help)
  order_circuit_simulation.py <grid_filename> <circuit_filename>
  order_circuit_simulation.py -g <grid_filename> -c <circuit_filename> [-v -o <output_file>]

Options:
  -h, --help                           Show this help
  -g, --grid=<grid_filename>           Grid filename
  -c, --circuit=<circuit_filename>     Circuit filename
  -o, --output-file=<output_file>      Print on <output_file> instead of stdout
  -v, --verbose                        Verbose output
"""

import itertools
import docopt
import math
import re

from typing import Iterable, Tuple

import cirq


def powerset(iterable):
    """Calculates the powerset of a set (e.g. [1,2] --> [[],[1],[2],[1,2]])."""
    s = list(iterable)
    return itertools.chain.from_iterable(
        itertools.combinations(s, r) for r in range(len(s) + 1))


class Bond:
    """Represents a wire in the graph."""

    def __init__(self, nodes):
        self.nodes = frozenset(nodes)
        self.ops = []

    def dim(self):
        """Returns the bond dimension for all operations this represents."""
        return len(self.ops)

    def add_op(self, op):
        self.ops.append(op)

    def swap_nodes(self, old_nodes, new_nodes):
        """Removes old_nodes from, and adds new_nodes to, this bond."""
        self.nodes = frozenset([n for n in self.nodes if n not in old_nodes] +
                               new_nodes)

    def get_remote_bonds(self):
        """Returns all bonds with nodes connected to this bond."""
        remotes = set()
        for n in self.nodes:
            for b in n.bonds:
                remotes.add(b)
        remotes.remove(self)
        return remotes

    def get_remote_nodes(self):
        """Returns all nodes adjacent to nodes connected to this bond."""
        remotes = []
        for b in self.get_remote_bonds():
            for n in b.nodes:
                if n not in self.nodes:
                    remotes.append(n)
        return set(remotes)

    def get_cost(self):
        """Returns (cost heuristic, future estimate, true cost)."""
        # Heuristic is size of resulting tensor divided by size of components.
        result_size = 2**sum([rb.dim() for rb in self.get_remote_bonds()])
        component_sizes = [2**n.rank() for n in self.nodes]
        heuristic = result_size / sum(component_sizes)

        # Estimate is 1 / average bond dimension of resulting tensor.
        # Maximizing average bond dimension leads to cheaper future contractions.
        estimate = 0
        remote_nodes = len(self.get_remote_nodes())
        if remote_nodes > 0:
            estimate = remote_nodes / result_size

        # Time cost is O(nmp) for nxm and mxp matrices. All tensors involved are
        # effectively matrices of dimension 2^n x 2^m; the shared dimension between
        # two tensors is the dimension of the shared bond.
        op_size = sum([n.rank() - self.dim() for n in self.nodes]) + self.dim()
        time_cost = 2**op_size
        return (heuristic, estimate, time_cost)

    def __hash__(self):
        return hash(id(self))

    def __eq__(self, other):
        return id(self) == id(other)

    def __repr__(self):
        return '{}'.format([str(op) for op in self.ops])

    def __str__(self):
        return '{}'.format([str(op) for op in self.ops])


class Node:
    """Represents a node in the graph."""

    def __init__(self, qubits):
        self.qubits = frozenset(qubits)
        self.bonds = set()

    def rank(self):
        """Returns the rank of the tensor this node represents."""
        return sum([b.dim() for b in self.bonds])

    def add_bond(self, bond):
        self.bonds.add(bond)

    def remove_bond(self, bond):
        self.bonds.remove(bond)

    def __hash__(self):
        return hash(id(self))

    def __eq__(self, other):
        return id(self) == id(other)


class Graph:
    """Represents a circuit graph."""

    def __init__(self, qubits):
        self.bonds = {}
        self.nodes = {frozenset([q]): Node([q]) for q in qubits}

    def create_node(self, qubits):
        """Adds a node to the graph and returns it."""
        self.nodes[frozenset(qubits)] = Node(frozenset(qubits))
        return self.nodes[frozenset(qubits)]

    def add_bond(self, qubits, op):
        """Adds a bond to the graph, or adds an op to an existing bond."""
        nodes = []
        for qs in powerset(qubits):
            if frozenset(qs) in self.nodes:
                nodes.append(self.nodes[frozenset(qs)])
        nodes = frozenset(nodes)
        if nodes in self.bonds:
            self.bonds[nodes].add_op(op)
            return
        new_bond = self.bonds[nodes] = Bond(nodes)
        new_bond.add_op(op)
        for n in nodes:
            n.add_bond(new_bond)

    def update_bonds(self, bonds):
        """Merge and re-insert modified bonds."""
        for b in bonds:
            if b.nodes in self.bonds.keys():
                for op in b.ops:
                    self.bonds[b.nodes].add_op(op)
            else:
                self.bonds[b.nodes] = b
                for n in b.nodes:
                    n.add_bond(b)

    def contract(self, bond):
        """Removes the given bond from the graph and merges adjacent nodes."""
        # Make new node.
        patch_qubits = []
        for n in bond.nodes:
            for q in n.qubits:
                patch_qubits.append(q)
        new_node = self.create_node(patch_qubits)
        # Update remote bonds.
        remote_bonds = bond.get_remote_bonds()
        for rb in remote_bonds:
            self.bonds.pop(rb.nodes)
            for n in rb.nodes:
                n.remove_bond(rb)
            rb.swap_nodes(bond.nodes, [new_node])
        self.update_bonds(remote_bonds)
        # Remove old nodes.
        for n in bond.nodes:
            self.nodes.pop(n.qubits)
        # Remove contracted bond.
        self.bonds.pop(bond.nodes)


def create_ordering_data(
        contraction_steps: Iterable[cirq.ops.raw_types.Operation],
        qubit_names: Iterable[int],
        qubit_order: Tuple[cirq.ops.raw_types.Qid, ...]):
    """Converts a sequence of Cirq operations to a qFlex contraction ordering.

  Args:
    contraction_steps: a list of cirq.raw_type.Operation(s).
    qubit_names: a list of integer IDs for each Qid in qubit_order.
    qubit_order: a tuple of qubits (cirq.ops.raw_types.Qid) in canonical order.

  Returns:
    A list of string-formatted qFlex contraction commands (e.g. 'expand 1 2')
    for the provided operations, with comments indicating the operation-to-
    command mapping. This does not include cut commands.
  """
    new_patch_name = ord('A')
    patches = {}
    output = []
    k = 0
    for step in contraction_steps:
        k += 1
        current_patch_1 = ''
        current_patch_2 = ''
        expand_qubits = []
        for q in step.ops[0].qubits:
            if q not in patches:
                expand_qubits.append(q)
            elif (not current_patch_1) or patches[q] == current_patch_1:
                current_patch_1 = patches[q]
            elif (not current_patch_2) or patches[q] == current_patch_2:
                current_patch_2 = patches[q]
            else:
                print('failed in create-ordering')
                return
        output.append('# {}: {}'.format(k, step))
        if not current_patch_1:
            current_patch_1 = chr(new_patch_name)
            new_patch_name += 1
        if not current_patch_2:
            for q in expand_qubits:
                patches[q] = current_patch_1
                output.append('expand {} {}'.format(
                    current_patch_1, qubit_names[qubit_order.index(q)]))
        else:
            if current_patch_1 > current_patch_2:
                current_patch_1, current_patch_2 = current_patch_2, current_patch_1
            output.append('merge {} {}'.format(current_patch_1,
                                               current_patch_2))
            for qubit in patches:
                if patches[qubit] == current_patch_1:
                    patches[qubit] = current_patch_2
    return output


def get_steps_for_graph(g: Graph):
    """Generates a qFlex circuit ordering from a Graph.

  Args:
    g: a Graph representing a Cirq circuit.

  Returns:
    A tuple of (total time cost in multiplications for this contraction,
                list of Bonds representing the contraction).
  """
    # Loop over the graph selecting bonds to contract.
    contraction_steps = []
    k = 0
    time_cost = 0
    while g.bonds:
        k += 1
        # (cost heuristic, future estimate, true cost)
        min_cost = (math.inf, math.inf, math.inf)
        min_bond = None
        for _, bond in g.bonds.items():
            cost = bond.get_cost()
            if cost < min_cost:
                min_cost = cost
                min_bond = bond
        g.contract(min_bond)
        time_cost += min_cost[2]
        contraction_steps.append(min_bond)
    return (time_cost, contraction_steps)


def circuit_to_ordering(circuit: cirq.circuits.Circuit,
                        qubit_names: Iterable[int] = None,
                        qubit_order_method: cirq.ops.QubitOrderOrList = cirq.
                        ops.QubitOrder.DEFAULT,
                        max_cuts: int = 2):
    """Generates a qFlex circuit ordering (with cuts) from a Cirq circuit.

  Args:
    circuit: a cirq.circuits.Circuit to be simulated.
    qubit_names: a list of integer IDs for each Qid in qubit_order. If left
      empty, qubits are assumed to have IDs {0, 1, ..., N}.
    qubit_order_method: an ops.QubitOrderOrList, which can call its order_for()
      method to get a canonical ordering of the qubits in the circuit.
    max_cuts: Maximum number of cuts attempted. Cuts are made using a greedy
      algorithm; making one or more cuts multiplies the time cost of this method
      by O(# of edges in circuit).

  Returns:
    A list of string-formatted qFlex contraction commands (e.g. 'expand 1 2')
    for the provided circuit, with comments indicating the operation-to-
    command mapping. This includes cut commands.

  Throws:
      ValueError: The gate can't be applied to the qubits.
  """
    if max_cuts < 0:
        raise ValueError('max_cuts must be positive!')

    if qubit_names is None:
        qubit_names = range(len(circuit.all_qubits()))

    qubit_order = qubit_order_method.order_for(circuit.all_qubits())

    cut_indices = set()
    min_cost = math.inf
    min_steps = []
    for _ in range(max_cuts):
        # Loop over possible cuts.
        min_cut = None
        tested_pairs = set()
        for cut_op in circuit.all_operations():
            cut_index = frozenset(cut_op.qubits)
            if cut_index in tested_pairs or cut_index in cut_indices:
                continue
            tested_pairs.add(cut_index)
            # Construct the graph from uncut operations in the circuit.
            g = Graph(circuit.all_qubits())
            for op in circuit.all_operations():
                index = frozenset(op.qubits)
                if index != cut_index and index not in cut_indices:
                    g.add_bond(op.qubits, op)
            cost, steps = get_steps_for_graph(g)
            if cost < min_cost:
                min_cost = cost
                min_cut = cut_op.qubits
                min_steps = steps
        if min_cut:
            cut_indices.add(frozenset(min_cut))
        else:  # Cut failed to reduce contraction cost; stop early.
            break
    order_data = []
    for cut in cut_indices:
        order_data.append('cut () %d %d' %
                          tuple(qubit_names[qubit_order.index(c)] for c in cut))
    order_data += create_ordering_data(min_steps, qubit_names, qubit_order)
    return order_data


def GetGridQubits(grid_stream):

    grid = [[y
             for y in x.strip()
             if y == '0' or y == '1']
            for x in grid_stream.readlines()
            if len(x) and x[0] != '#']

    # Get the number of rows
    grid_I = len(grid)

    # Check that the number of columns is consistent
    if len(set(len(x) for x in grid)) != 1:
        raise AssertionError("Number of columns in grid are not consistent.")

    # Get the number of columns
    grid_J = len(grid[0])

    # Return cirq.GridQubit
    return {
        I * grid_J + J: cirq.GridQubit(I, J)
        for I in range(grid_I) for J in range(grid_J) if grid[I][J] == '1'
    }


def GetGate(line, qubits):

    # Get map from gate name to cirq
    # Ordering NOT SUPPORT two-qubit gates with Schmidt decomposition != 2
    gates_map = {}
    gates_map['h'] = cirq.H
    gates_map['x'] = cirq.X
    gates_map['z'] = cirq.Z
    gates_map['t'] = cirq.Z**(0.25)
    gates_map['x_1_2'] = cirq.X**(0.5)
    gates_map['y_1_2'] = cirq.Y**(0.5)
    gates_map['h_1_2'] = cirq.H**(0.5)
    gates_map['cz'] = cirq.CZ
    gates_map['cx'] = cirq.CNOT
    gates_map['rz'] = cirq.Rz

    # Remove last cr
    line = line.strip()

    # Remove everything after '#'
    line = re.sub(r"#.*", r"", line)

    # Remove any special character
    line = re.sub(r"[^)(\s\ta-zA-Z0-9_.,-]", r"", line)

    # Convert tabs to spaces
    line = re.sub(r"[\t]", " ", line)

    # Remove multiple spaces
    line = re.sub(r"[\s]{2,}", r" ", line)

    # Remove last space
    line = re.sub(r"\s+$", r"", line)

    # Remove any space between a non-space char and '('
    line = re.sub(r"[\s]+[(]", r"(", line)

    # Remove spaces between parentheses
    line = re.sub(r"\s+(?=[^()]*\))", r"", line)

    # After stripping, line should follow the format
    # 0 gate(p1,p2,...) q1 q2 ...

    # Only one open parenthesis is allowed, followed by one closed
    if line.count('(') == 1:
        if line.count(')') != 1:
            raise AssertionError('ERROR: Open parenthesis is not matched.')
    elif line.count('(') > 1:
        raise AssertionError('ERROR: Too many open parentheses.')
    elif line.count(')') != 0:
        raise AssertionError('ERROR: Too many close parentheses.')

    line = line.split()
    cycle = int(line[0])
    gate_qubits = [int(x) for x in line[2:]]
    if line[1].count('('):
        gate_name, params = line[1].split('(')
        params = [float(x) for x in params.replace(')', '').split(',')]
    else:
        gate_name = line[1]
        params = None

    if not gate_name in gates_map:
        raise AssertionError(
            "ERROR: Gate {} not supported yet.".format(gate_name))

    if params == None:
        return gates_map[gate_name](*[qubits[q] for q in gate_qubits])
    else:
        return gates_map[gate_name](*params)(*[qubits[q] for q in gate_qubits])


def GetCircuit(circuit_stream, qubits):

    circuit = cirq.Circuit()
    circuit.append(
        gate for gate in (GetGate(line, qubits)
                          for line in circuit_stream
                          if len(line) and len(line.strip().split()) > 1)
        if len(gate.qubits) == 2)
    return circuit


if __name__ == "__main__":

    from docopt import docopt
    from sys import stderr

    args = docopt(__doc__)

    grid_filename = args['--grid'] if args['--grid'] != None else args[
        '<grid_filename>']
    circuit_filename = args['--circuit'] if args['--circuit'] != None else args[
        '<circuit_filename>']
    output_filename = args['--output-file']
    verbose = args['--verbose']

    if (verbose):
        print('Get grid.', file=stderr)
    with open(grid_filename) as f:
        qubits = GetGridQubits(f)

    if (verbose):
        print('Get circuit.', file=stderr)
    with open(circuit_filename) as f:
        circuit = GetCircuit(f, qubits)

    if (verbose):
        print('Compute ordering.', file=stderr)

    if output_filename != None:
        with open(output_filename, 'w') as f:
            for line in circuit_to_ordering(circuit):
                print(line, file=f)
    else:
        for line in circuit_to_ordering(circuit):
            print(line)
