import pytest
import cirq

import numpy as np

# Due to the directory structure of QFlex
# go up in hierarchy twice
import sys
sys.path.insert(1, '../../')

from qflexcirq import QFlexCircuit


def test_qflexcircuit_equality():
    qubit_1 = cirq.GridQubit(0, 0)
    qubit_2 = cirq.GridQubit(0, 1)

    from qflexcirq.interface.qflex_virtual_device import QFlexVirtualDevice
    dummy_device_1 = QFlexVirtualDevice()
    dummy_device_1.get_grid_qubits_as_keys = lambda: {qubit_1: 1, qubit_2: 2}

    from qflexcirq.interface.qflex_order import QFlexOrder
    dummy_order = QFlexOrder("This is dummy order")

    circuit = cirq.Circuit()
    circuit.append(cirq.Moment([cirq.ops.SWAP.on(qubit_1, qubit_2)]))

    qflexcirc1 = QFlexCircuit(circuit,
                              dummy_device_1,
                              dummy_order,
                              allow_decomposition=True)

    qflexcirc2 = QFlexCircuit(circuit,
                              dummy_device_1,
                              dummy_order,
                              allow_decomposition=True)

    assert (qflexcirc1 == qflexcirc2)

    circuit = cirq.Circuit()
    circuit.append(cirq.Moment([cirq.ops.CNOT.on(qubit_1, qubit_2)]))
    qflexcirc3 = QFlexCircuit(circuit,
                              dummy_device_1,
                              dummy_order,
                              allow_decomposition=True)
    assert (qflexcirc1 != qflexcirc3)


def test_constructor_decomposition():

    qubit_1 = cirq.GridQubit(0, 0)
    qubit_2 = cirq.GridQubit(0, 1)

    circuit = cirq.Circuit()
    circuit.append(cirq.Moment([cirq.ops.SWAP.on(qubit_1, qubit_2)]))

    from qflexcirq.interface.qflex_virtual_device import QFlexVirtualDevice
    dummy_device_1 = QFlexVirtualDevice()
    dummy_device_1.get_grid_qubits_as_keys = lambda: {qubit_1: 1, qubit_2: 2}

    from qflexcirq.interface.qflex_order import QFlexOrder
    dummy_order = QFlexOrder("This is dummy order")

    qcircuit1 = QFlexCircuit(circuit,
                             dummy_device_1,
                             dummy_order,
                             allow_decomposition=True)
    # There are three CNOTs
    assert (len(list(qcircuit1.all_operations())) == 3)

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
    s1 = sympy.Symbol("theta")
    s2 = sympy.Symbol("phi")

    operation1 = cirq.Z.on(cirq.GridQubit(0, 0))**s1
    circuit.append(operation1)

    operation2 = cirq.FSimGate(s1, s2).on(cirq.GridQubit(0, 0),
                                          cirq.GridQubit(0, 1))
    circuit.append(operation2)

    from qflexcirq.interface.qflex_virtual_device import QFlexVirtualDevice
    dummy_device = QFlexVirtualDevice()

    # Override the function to accept anything
    # dummy_device.is_qflex_virt_dev_op = lambda x : True
    # Return a single pair gridqubit : index
    dummy_device.get_grid_qubits_as_keys = lambda: {
        cirq.GridQubit(0, 0): 1,
        cirq.GridQubit(0, 1): 2
    }

    from qflexcirq.interface.qflex_order import QFlexOrder
    dummy_order = QFlexOrder("This is a dummy order")

    # save the translator
    original_translator = QFlexCircuit.translate_cirq_to_qflex

    # place a dummy translator
    QFlexCircuit.translate_cirq_to_qflex = lambda x, y: "dummy contents"

    # Create a QFlexCircuit
    qcircuit = QFlexCircuit(circuit, dummy_device, dummy_order)

    # put back the translator
    QFlexCircuit.translate_cirq_to_qflex = original_translator

    param_solver = {"theta": 0.3, "phi": 0.7}
    solved_circuit = cirq.protocols.resolve_parameters(qcircuit, param_solver)

    assert (solved_circuit[0].operations[0].gate.exponent == 0.3)
    assert (solved_circuit[1].operations[0].gate.theta == 0.3)
    assert (solved_circuit[1].operations[0].gate.phi == 0.7)

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
    circuit.append(
        cirq.Moment([
            cirq.ops.PhasedXPowGate(phase_exponent=0.25,
                                    exponent=0.5).on(qubit_2)
        ]))
    circuit.append(
        cirq.Moment(
            [cirq.ops.FSimGate(np.pi / 2, np.pi / 16).on(qubit_2, qubit_1)]))

    file_lines = [
        "2", "0 h 1", "1 t 1", "2 x_1_2 1", "3 y_1_2 1", "4 cx 1 2", "5 cz 1 2",
        "6 rz(0.7) 1", "7 hz_1_2 2", "8 fsim(0.5, 0.0625) 2 1"
    ]

    qubits_to_index_dict = {qubit_1: 1, qubit_2: 2}
    translation = QFlexCircuit.translate_cirq_to_qflex(
        cirq_circuit=circuit, qubit_to_index_dict=qubits_to_index_dict)

    assert (translation.strip() == "\n".join(file_lines))


def test_supported_gate_set():
    """
    Is the circuit validated correctly on append?
    """

    import qflexcirq.interface.qflex_simulator as qsim
    import qflexcirq.interface.qflex_virtual_device as qdevice
    import qflexcirq.interface.qflex_grid as qgrid
    import qflexcirq.interface.qflex_circuit as qcirc
    import qflexcirq.interface.qflex_order as qorder

    qdev = qdevice.QFlexVirtualDevice(
        qflex_grid=qgrid.QFlexGrid.create_rectangular(4, 4))
    qord = qorder.QFlexOrder("Dummy Ordering")

    qubits = qdev.get_indexed_grid_qubits()

    cirq_circuit = cirq.Circuit()
    # cirq_circuit.append(cirq.ops.TOFFOLI.on(qubits[0], qubits[1], qubits[2]))

    my_circuit = QFlexCircuit(cirq_circuit=cirq_circuit,
                              device=qdev,
                              qflex_order=qord,
                              allow_decomposition=False)

    # Unparameterized gates
    my_circuit.append(cirq.Moment([cirq.ops.CZ.on(qubits[0], qubits[1])]))
    my_circuit.append(cirq.Moment([cirq.ops.CNOT.on(qubits[0], qubits[1])]))
    my_circuit.append(cirq.Moment([cirq.ops.H.on(qubits[0])]))
    my_circuit.append(cirq.Moment([cirq.ops.T.on(qubits[0])]))
    my_circuit.append(
        cirq.Moment([cirq.ops.ZPowGate(exponent=0.7).on(qubits[0])]))
    my_circuit.append(
        cirq.Moment([cirq.ops.XPowGate(exponent=0.5).on(qubits[0])]))
    my_circuit.append(
        cirq.Moment([cirq.ops.YPowGate(exponent=0.5).on(qubits[0])]))
    my_circuit.append(
        cirq.Moment([
            cirq.ops.PhasedXPowGate(phase_exponent=0.25,
                                    exponent=0.5).on(qubits[0])
        ]))
    my_circuit.append(
        cirq.Moment(
            [cirq.ops.FSimGate(np.pi / 2, np.pi / 16).on(qubits[0],
                                                         qubits[1])]))

    # Parameterized gates
    import sympy
    s1 = sympy.Symbol("theta")
    s2 = sympy.Symbol("phi")
    my_circuit.append(
        cirq.Moment([cirq.ops.ZPowGate(exponent=s1).on(qubits[0])]))
    my_circuit.append(
        cirq.Moment([cirq.ops.FSimGate(s1, s2).on(qubits[0], qubits[1])]))

    # Should fail
    with pytest.raises(ValueError):
        my_circuit.append(
            cirq.Moment([cirq.ops.TOFFOLI.on(qubits[0], qubits[1], qubits[2])]))

    with pytest.raises(ValueError):
        s3 = sympy.Symbol("phail")
        my_circuit.append(
            cirq.Moment([cirq.ops.XPowGate(exponent=s3).on(qubits[0])]))
