import cirq

qubit_dict = {}

def from_index_to_grid_qubit(index, rows=11, cols=12):

    if index in qubit_dict.keys():
        return qubit_dict[index]

    row = index // cols
    col = index % cols

    if row >= rows:
        raise ValueError("Wrong maximum of rows?")

    qub = cirq.GridQubit(row, col)
    qubit_dict[index] = qub

    return qub

def convert_qflex_circuit_file(fname):

    SQRTX = cirq.XPowGate(exponent=0.5)
    SQRTY = cirq.YPowGate(exponent=0.5)

    moment_index = -1
    current_moment = []
    moments = []

    with open(fname, "r") as file:
        for line in file:

            parts = line.strip().split(" ")

            if len(parts) <= 1:
                # skip number of qubits line
                # skip empty lines
                continue

            if moment_index != int(parts[0]):
                moment_index = int(parts[0])
                if len(current_moment) > 0:
                    moments.append(cirq.Moment(current_moment))
                current_moment = []

            qubit1 = from_index_to_grid_qubit(int(parts[2]), rows=11, cols=12)

            if parts[1] == "h":
                current_moment.append(cirq.H.on(qubit1))
            elif parts[1] =="x_1_2":
                current_moment.append(SQRTX.on(qubit1))
            elif parts[1] == "y_1_2":
                current_moment.append(SQRTY.on(qubit1))
            elif parts[1] =="t":
                current_moment.append(cirq.T.on(qubit1))
            elif parts[1] == "cz":
                qubit2 = from_index_to_grid_qubit(int(parts[3]), rows=11, cols=12)
                current_moment.append(cirq.CZ.on(qubit1, qubit2))
            else:
                print("not parsed:" + line)

    print("Constructed {} moments".format(moment_index))

    circuit = cirq.Circuit(moments=moments)
    return circuit
