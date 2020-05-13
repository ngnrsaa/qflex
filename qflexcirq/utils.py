# Lint as: python3
"""
Provides utils for qFlex.
"""

import numpy as np
import cirq
import re


def ComputeSchmidtRank(gate):

    if len(gate.qubits) == 1:
        return 1
    if len(gate.qubits) > 2:
        raise AssertionError("Not yet implemented.")

    V, S, W = np.linalg.svd(
        np.reshape(
            np.einsum('abcd->acbd', np.reshape(cirq.unitary(gate),
                                               [2, 2, 2, 2])), [4, 4]))

    return sum(S != 0)


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
    return GetMomentAndGate(line, qubits)[1]


def GetMomentAndGate(line, qubits):

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
    gates_map['rz'] = cirq.ZPowGate
    gates_map['hz_1_2'] = cirq.PhasedXPowGate(phase_exponent=0.25, exponent=0.5)
    gates_map['fsim'] = cirq.FSimGate

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

    if params is None:
        return cycle, gates_map[gate_name](*[qubits[q] for q in gate_qubits])
    else:
        if gate_name == "fsim":
            # the Cirq Fsim gate takes angles and not exponents
            # Transform the params
            params = [p * np.pi for p in params]
            return cycle, gates_map[gate_name](
                theta=params[0],
                phi=params[1])(*[qubits[q] for q in gate_qubits])

        if gate_name == "rz":
            # the Cirq Fsim gate takes angles and not exponents
            # Transform the params
            # params = [p * np.pi for p in params]
            return cycle, gates_map[gate_name](exponent=params[0])(
                qubits[gate_qubits[0]])

        return cycle, gates_map[gate_name](*params)(
            *[qubits[q] for q in gate_qubits])


def GetCircuit(circuit_stream, qubits):

    circuit = cirq.Circuit()
    circuit.append(
        gate for gate in (GetGate(line, qubits)
                          for line in circuit_stream
                          if len(line) and len(line.strip().split()) > 1))
    return circuit


def GetCircuitOfMoments(file_name, qubits):

    with open(file_name, "r") as circuit_stream:

        moment_index = -1
        current_moment = []
        moments = []

        for line in circuit_stream:
            if not (len(line) and len(line.strip().split()) > 1):
                continue

            parts = GetMomentAndGate(line, qubits)

            moment_idx_dif = int(parts[0]) - moment_index

            if moment_idx_dif != 0:

                moment_index = int(parts[0])
                if len(current_moment) > 0:
                    moments.append(cirq.Moment(current_moment))

                for mi in range(moment_idx_dif - 1):
                    # add empty moments
                    moments.append(cirq.Moment([]))

                current_moment = []

            current_moment.append(parts[1])

        if len(current_moment) > 0:
            moments.append(cirq.Moment(current_moment))

        return cirq.Circuit(moments)


def GetNumberOfQubits(cirq_circuit):
    """
    Determine the number of qubits from an unknown Cirq circuit
    :param cirq_circuit:
    :return:
    """
    known_qubits = {}
    size = 0
    for operation in cirq_circuit:
        for qub in operation.qubits:
            if not qub in known_qubits:
                size += 1
                known_qubits[qub] = size
    return size


def GetGridQubitFromIndex(index, rows=11, cols=12):

    row = index // cols
    col = index % cols

    if row >= rows:
        raise ValueError("Wrong maximum of rows?")

    qub = cirq.GridQubit(row, col)

    return qub


def GetIndexFromGridQubit(grid_qubit, rows=11, cols=12):

    if grid_qubit.row >= rows:
        raise ValueError("This GridQubit seems to have wrong row coordinate...")

    return grid_qubit.row * cols + grid_qubit.col
