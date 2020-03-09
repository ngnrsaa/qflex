# Lint as: python3
"""
Provides a method for converting Cirq circuits to qFlex simulation order.

Usage:
  order_circuit_simulation.py (-h | --help)
  order_circuit_simulation.py <grid_filename> <circuit_filename> [-v -x <cut_count> -f <fidelity_ratio> -o <output_file>]
  order_circuit_simulation.py -g <grid_filename> -c <circuit_filename> [-v -x <cut_count> -f <fidelity_ratio> -o <output_file>]

Options:
  -h, --help                           Show this help
  -g, --grid=<grid_filename>           Grid filename
  -c, --circuit=<circuit_filename>     Circuit filename
  -x, --max_cuts=<cut_count>           Number of cuts to attempt [default: 2]
  -f, --fidelity=<fidelity_ratio>      Circuit fidelity (enforced with cuts)
  -o, --output-file=<output_file>      Print on <output_file> instead of stdout
  -v, --verbose                        Verbose output
"""

import itertools
import math

from typing import Dict, Iterable, Tuple

import cirq
import os, sys

# Get real path of qflex modules
sys.path.insert(
    1, os.path.realpath(os.path.dirname(os.path.realpath(__file__)) + '/../'))
import utils


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
        self._dim = 0

    def dim(self):
        """Returns the bond dimension for all operations this represents.
        Bond dimension is cached to reduce cost from computing operation ranks.
        """
        if not self._dim:
            self._dim = sum(
                [math.log2(utils.ComputeSchmidtRank(op)) for op in self.ops])
        return self._dim

    def add_op(self, op):
        """Add 'op' to this bond and cache the new bond dimension."""
        self.ops.append(op)
        self._dim = sum(
            [math.log2(utils.ComputeSchmidtRank(op)) for op in self.ops])

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
    qubit_names: Iterable[int], qubit_order: Tuple[cirq.ops.raw_types.Qid,
                                                   ...]):
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


def match_fidelity(target_fidelity: float, cuts: Dict[frozenset, int]):
    """Determines the "width" of each cut (what values to evaluate over) to
    match the target fidelity as closely as possible.

  Args:
    target_fidelity: the requested fidelity. The final fidelity must be
      greater than this, but should otherwise be as close as possible.
    cuts: a dict of (cut_qubits, bond_dim) for all cuts.

  Returns:
    The final fidelity approximation and a list of tuples
      (cut qubits, cut width, cut_dimension).
  """
    if (target_fidelity == 1.0 or not cuts):
        return (1.0,
                zip(list(cuts.keys()), list(cuts.values()),
                    list(cuts.values())))

    value_sets = []
    denominator = 1.0
    for cut, dim in cuts.items():
        # Using min() prevents float math from producing dim+1.
        min_width = min(dim, math.ceil(dim * target_fidelity))
        value_sets.append(list(range(min_width, dim + 1)))
        denominator *= dim

    combinations = itertools.product(*value_sets)
    closest_set = []
    closest_fidelity = 1.0
    for width_set in combinations:
        fidelity = 1.0
        for width in width_set:
            fidelity *= width
        fidelity /= denominator
        if fidelity >= target_fidelity and fidelity < closest_fidelity:
            closest_set = width_set
            closest_fidelity = fidelity
            if fidelity == target_fidelity:
                break

    return (closest_fidelity,
            zip(list(cuts.keys()), closest_set, list(cuts.values())))


def circuit_to_ordering(
    circuit: cirq.circuits.Circuit,
    qubit_names: Iterable[int] = None,
    qubit_order_method: cirq.ops.QubitOrderOrList = cirq.ops.QubitOrder.DEFAULT,
    max_cuts: int = 2,
    fidelity: float = 1.0):
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
    fidelity: If this is less than 1.0, cuts generated will only be evaluated
      on a subset of the possible values. Subset size is controlled to provide
      closest fidelity greater than the provided number.

  Returns:
    A list of string-formatted qFlex contraction commands (e.g. 'expand 1 2')
    for the provided circuit, with comments indicating the operation-to-
    command mapping. This includes cut commands.

  Throws:
      ValueError: The gate can't be applied to the qubits.
  """

    # Check that there are no k-qubit gates with k > 2
    if sum(len(g.qubits) > 2 for g in circuit.all_operations()):
        raise AssertionError(
            "Auto-ordering now working for k-qubit gates with k > 2.")

    # Strip single qubit gates
    if sum(len(g.qubits) == 1 for g in circuit.all_operations()):
        circuit = cirq.Circuit([
            cirq.Moment([g])
            for g in circuit.all_operations()
            if len(g.qubits) > 1
        ])

    if max_cuts < 0:
        raise ValueError('max_cuts must be positive!')

    if qubit_names is None:
        qubit_names = range(len(circuit.all_qubits()))

    qubit_order = qubit_order_method.order_for(circuit.all_qubits())

    cut_indices = {}
    # Make initial pass prior to cuts.
    g = Graph(circuit.all_qubits())
    for op in circuit.all_operations():
        g.add_bond(op.qubits, op)
    min_cost, min_steps = get_steps_for_graph(g)

    for cut_num in range(max_cuts):
        # Loop over possible cuts.
        min_cut = None
        tested_pairs = set()
        min_cut_dim = 1
        for cut_op in circuit.all_operations():
            cut_dim = 1
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
                elif index == cut_index:
                    cut_dim *= utils.ComputeSchmidtRank(op)
            cost, steps = get_steps_for_graph(g)
            if cost < min_cost:
                min_cost = cost
                min_cut = cut_op.qubits
                min_cut_dim = cut_dim
                min_steps = steps
        if min_cut:
            cut_indices[frozenset(min_cut)] = min_cut_dim
        else:  # Cut failed to reduce contraction cost; stop early.
            break
    order_data = []
    final_fidelity, cut_ratios = match_fidelity(fidelity, cut_indices)
    for cut, width, dim in cut_ratios:
        cut_qubits = tuple(cut)
        cut_ids = tuple(qubit_names[qubit_order.index(c)] for c in cut)
        order_data.append('# cut bond {} of dim={}'.format(cut_qubits, dim))
        # If cut evaluates all values, value list can be removed.
        cut_values = tuple(range(width)) if width < dim else tuple()
        order_data.append('cut {} {} {}'.format(cut_values, cut_ids[0],
                                                cut_ids[1]))
    order_data.append('# approximate fidelity: {}'.format(final_fidelity))
    order_data += create_ordering_data(min_steps, qubit_names, qubit_order)

    return order_data


if __name__ == "__main__":

    from docopt import docopt
    from sys import stdout, stderr

    args = docopt(__doc__)

    grid_filename = args['--grid'] if args['--grid'] != None else args[
        '<grid_filename>']
    circuit_filename = args['--circuit'] if args['--circuit'] != None else args[
        '<circuit_filename>']
    max_cuts = int(args['--max_cuts']) if args['--max_cuts'] != None else 2
    fidelity = float(args['--fidelity']) if args['--fidelity'] != None else 1.0
    output_filename = args['--output-file']
    verbose = args['--verbose']

    if verbose:
        print('Get grid.', file=stderr)
    with open(grid_filename) as f:
        qubits = utils.GetGridQubits(f)

    if verbose:
        print('Get circuit.', file=stderr)
    with open(circuit_filename) as f:
        circuit = utils.GetCircuit(f, qubits)

    if verbose:
        print('Compute ordering.', file=stderr)

    ordering = circuit_to_ordering(circuit=circuit,
                                   qubit_names=sorted(qubits),
                                   max_cuts=max_cuts,
                                   fidelity=fidelity)
    with (stdout
          if output_filename is None else open(output_filename, 'w')) as f:
        print('\n'.join(ordering), file=f)
