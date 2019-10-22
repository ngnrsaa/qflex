import cirq

# Due to the directory structure of QFlex
# go up in hierarchy twice
import sys
sys.path.insert(1, '../../')

from python.cirq_interface.qflex_circuit import QFlexCircuit

def test_resolve_parameters():
    # After a circuit with symbols was resolved
    # Does the resulting circuit still have the same device and grid?


def test_translate_cirq_to_qflex():

    qubit_1 = cirq.GridQubit(0, 0)
    qubit_2 = cirq.GridQubit(0, 1)

    circuit = cirq.Circuit()
    circuit.append(cirq.Moment([cirq.ops.H.on(qubit_1)]))
    circuit.append(cirq.Moment([cirq.ops.T.on(qubit_1)]))
    circuit.append(cirq.Moment([cirq.ops.XPowGate(exponent=0.5).on(qubit_1)]))
    circuit.append(cirq.Moment([cirq.ops.YPowGate(exponent=0.5).on(qubit_1)]))
    circuit.append(cirq.Moment([cirq.ops.CNOT.on(qubit_1, qubit_2)]))
    circuit.append(cirq.Moment([cirq.ops.CZ.on(qubit_1, qubit_2)]))

    file_lines = [
        "2",
        "0 h 1",
        "1 t 1",
        "2 x_1_2 1",
        "3 y_1_2 1",
        "4 cx 1 2",
        "5 cz 1 2"
    ]

    qubits_to_index_dict = { qubit_1: 1, qubit_2: 2}
    translation = QFlexCircuit.translate_cirq_to_qflex(cirq_circuit=circuit,
                                                       qubit_to_index_dict = qubits_to_index_dict)

    assert(translation.strip() == "\n".join(file_lines))



