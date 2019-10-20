#!/usr/bin/env python3
"""
Usage:
  simulator.py (-h | --help)
  simulator.py <grid_filename> <circuit_filename>
  simulator.py -g <grid_filename> -c <circuit_filename> [-v]

Options:
  -h, --help                           Show this help
  -g, --grid=<grid_filename>           Grid filename
  -c, --circuit=<circuit_filename>     Circuit filename
  -v, --verbose                        Verbose output
"""

from sys import stderr
import cirq

# Avaiable single/two qubit gates
single_qubit_gates = {'h', 'x_1_2', 'y_1_2', 'h_1_2', 't'}
two_qubit_gates = {'cx', 'cz'}

# Get map from gate name to cirq
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


def GetGridQubit(grid_stream):

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
    return [
        cirq.GridQubit(I, J)
        for I in range(grid_I)
        for J in range(grid_J)
        if grid[I][J] == '1'
    ]


def GetCircuit(circuit_stream, qubits):

    # Create cirq.Circuit
    circuit = cirq.Circuit()

    # The first line must be the number of active qubits
    if len(qubits) != int(circuit_stream.readline()):
        raise AssertionError("Number of active qubits is incosistent.")

    # Read circuit
    for line in circuit_stream.readlines():
        if len(line) and line[0] != '#':

            line = line.strip().split()[1:]
            gate = line[0]

            if gate in single_qubit_gates:
                q1 = int(line[1])
                circuit.append([gates_map[gate](qubits[q1])])
            elif gate in two_qubit_gates:
                q1 = int(line[1])
                q2 = int(line[2])
                circuit.append([gates_map[gate](qubits[q1], qubits[q2])])
            else:
                raise AssertionError('Gate \'{}\' not recognized.'.format(gate))

    return circuit


if __name__ == "__main__":

    from docopt import docopt

    args = docopt(__doc__)

    grid_filename = args['--grid'] if args['--grid'] != None else args[
        '<grid_filename>']
    circuit_filename = args['--circuit'] if args['--circuit'] != None else args[
        '<circuit_filename>']
    verbose = args['--verbose']

    # Get grid from file
    with open(grid_filename) as grid_stream:
        qubits = GetGridQubit(grid_stream)

    # Get circuit from file
    with open(circuit_filename) as circuit_stream:
        circuit = GetCircuit(circuit_stream, qubits)

    # Print quantum circuit
    if verbose: print(circuit, file=stderr)

    # Print amplitudes
    simulator = cirq.Simulator()
    results = simulator.simulate(circuit, qubit_order=qubits)

    for x, v in enumerate(results.final_state):
        print(('{{: {}d}}: |{{}}>'.format(len(str(len(qubits))) + 2)).format(
            x,
            bin(x)[2:].zfill(len(qubits))[::-1]), v)
