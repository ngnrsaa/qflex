import pytest
import cirq

import numpy as np

# Due to the directory structure of QFlex
# go up in hierarchy twice
import sys
sys.path.insert(1, '../../')

from python.cirq_interface.qflex_circuit import QFlexCircuit

# Used to include a class which does not exist in Cirq 0.5.0
import python.cirq_interface.fsim_gate as cirqtmp

def test_constructor_decomposition():

    qubit_1 = cirq.GridQubit(0, 0)
    qubit_2 = cirq.GridQubit(0, 1)

    circuit = cirq.Circuit()
    circuit.append(cirq.Moment([cirq.ops.SWAP.on(qubit_1, qubit_2)]))

    from python.cirq_interface.qflex_virtual_device import QFlexVirtualDevice
    dummy_device_1 = QFlexVirtualDevice()
    # Override the function to accept anything
    dummy_device_1.get_grid_qubits_as_keys = lambda: {qubit_1: 1, qubit_2: 2}

    from python.cirq_interface.qflex_order import QFlexOrder
    dummy_order = QFlexOrder("This is dummy order")

    qcircuit1 = QFlexCircuit(circuit,
                            dummy_device_1,
                            dummy_order,
                            allow_decomposition=True)
    # There are three CNOTs
    assert(len(list(qcircuit1.all_operations())) == 3)

    # This should crash because SWAP is not supported by the device
    with pytest.raises(ValueError):
        qcircuit2 = QFlexCircuit(circuit,
                                 dummy_device_1,
                                 dummy_order,
                                 allow_decomposition=False)



def test_resolve_parameters():
    # After a circuit with symbols was resolved
    # Does the resulting circuit still have the same device and grid?
    # In general, resolving parameters will result in a strange circuit
    # But in this test I am placing a single gate and setting its exponent

    circuit = cirq.Circuit()

    import sympy
    s = sympy.Symbol("symbol")
    operation = cirq.X.on(cirq.GridQubit(0,0)) ** s
    circuit.append(operation)

    from python.cirq_interface.qflex_virtual_device import QFlexVirtualDevice
    dummy_device = QFlexVirtualDevice()

    # Override the function to accept anything
    dummy_device.is_qflex_virt_dev_op = lambda x : True
    # Return a single pair gridqubit : index
    dummy_device.get_grid_qubits_as_keys = lambda : {cirq.GridQubit(0,0) : 99}

    from python.cirq_interface.qflex_order import QFlexOrder
    dummy_order = QFlexOrder("This is a dummy order")

    # save the translator
    original_translator = QFlexCircuit.translate_cirq_to_qflex

    # place a dummy translator
    QFlexCircuit.translate_cirq_to_qflex = lambda x, y: "dummy contents"

    # Create a QFlexCircuit
    qcircuit = QFlexCircuit(circuit, dummy_device, dummy_order)

    # put back the translator
    QFlexCircuit.translate_cirq_to_qflex =original_translator

    param_solver = {"symbol" : 0.5}
    solved_circuit = cirq.protocols.resolve_parameters(qcircuit, param_solver)

    assert (solved_circuit[0].operations[0].gate.exponent == 0.5)
    assert (solved_circuit.device is dummy_device)
    assert (solved_circuit._own_order is dummy_order)


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
    circuit.append(cirq.Moment([cirq.ops.ZPowGate(exponent=0.7).on(qubit_1)]))
    circuit.append(cirq.Moment([cirq.ops.PhasedXPowGate(phase_exponent=0.25, exponent=0.5).on(qubit_2)]))
    circuit.append(cirq.Moment([cirqtmp.FSimGate(0, np.pi).on(qubit_2, qubit_1)]))

    file_lines = [
        "2",
        "0 h 1",
        "1 t 1",
        "2 x_1_2 1",
        "3 y_1_2 1",
        "4 cx 1 2",
        "5 cz 1 2",
        "6 rz(0.7) 1",
        "7 hz_1_2 2",
        "8 fsim(0.0, 1.0) 2 1"
    ]

    qubits_to_index_dict = { qubit_1: 1, qubit_2: 2}
    translation = QFlexCircuit.translate_cirq_to_qflex(cirq_circuit=circuit,
                                                       qubit_to_index_dict = qubits_to_index_dict)

    assert(translation.strip() == "\n".join(file_lines))