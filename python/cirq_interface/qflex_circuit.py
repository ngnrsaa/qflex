import tempfile
import os

import numpy as np

import cirq

import python.cirq_interface.qflex_virtual_device as qdevice
import python.cirq_interface.qflex_order as qorder

import python.utils as qflexutils

# Used to include a class which does not exist in Cirq 0.5.0
import python.cirq_interface.fsim_gate as cirqtmp

class QFlexCircuit(cirq.Circuit):

    def __init__(self,
                 cirq_circuit,
                 device,
                 qflex_order = None,
                 allow_decomposition = False):

        if (device is None) \
                or (not isinstance(device, qdevice.QFlexVirtualDevice)):
            raise ValueError("QFlexVirtualDevice necessary for constructor!")

        if qflex_order is None:
            # No order was specified.
            # Construct and order from a circuit
            qubits = device.get_indexed_grid_qubits()
            self._own_order = qorder.QFlexOrder(cirq_circuit = cirq_circuit,
                                                qubits = qubits)
        elif isinstance(qflex_order, qorder.QFlexOrder):
            self._own_order = qflex_order
        else:
            raise ValueError("{!r} is not of a QFlexOrder!")


        if allow_decomposition:
            super().__init__([], device)
            for moment in cirq_circuit:
                for op in moment:
                    # This should call decompose on the gates
                    self.append(op)
        else:
            # This super constructor does not call decompose?
            super().__init__(cirq_circuit, device)

        # Behind the scene, this class creates a temporary file for each object
        self._file_handle = tempfile.mkstemp()

        with open(self._file_handle[1], "w") as f:
            # The cirq_circuit has QFlexVirtualDevice
            qubit_to_index_dict = self.device.get_grid_qubits_as_keys()
            print(QFlexCircuit.translate_cirq_to_qflex(self,
                                                       qubit_to_index_dict),
                  file = f)


    @property
    def circuit_data(self):
        return self._file_handle[1]

    @property
    def ordering_data(self):
        return self._own_order._file_handle[1]


    def _resolve_parameters_(self, param_resolver: cirq.study.ParamResolver):

        qflex_circuit = super()._resolve_parameters_(param_resolver)

        qflex_circuit.device = self.device
        qflex_circuit._own_order = self._own_order

        return qflex_circuit


    def __del__(self):
        # The destructor removes the temporary file

        try:
            os.close(self._file_handle[0])
        except OSError as e:
            if e.errno == 9:
                # if it was closed before
                pass
            else:
                raise e

        # remove the temporary file from disk
        os.remove(self._file_handle[1])


    @staticmethod
    def translate_cirq_to_qflex(cirq_circuit, qubit_to_index_dict):
        """
        Translates a Cirq/QFlex circuit to the QFlex representation
        :param cirq_circuit: the initial circuit
        :param qubit_to_index_dict: a dictionary mapping Cirq Qubits to integers
        :return: the string representing line-by-line the QFlex circuit
        """

        number_qubits = qflexutils.GetNumberOfQubits(cirq_circuit)

        first_line = str(number_qubits) + "\n"
        circuit_data = [first_line]

        for mi, moment in enumerate(cirq_circuit):
            for op in moment:

                qub_str = ""
                for qub in op.qubits:
                    qub_str += "{} ".format(qubit_to_index_dict[qub])

                qflex_gate = ""
                if isinstance(op.gate, cirq.ops.CZPowGate)\
                        and op.gate.exponent == 1.0:
                    qflex_gate = "cz"
                elif isinstance(op.gate, cirq.ops.CNotPowGate) \
                        and op.gate.exponent == 1.0:
                    qflex_gate = "cx"
                elif isinstance(op.gate, cirq.ops.HPowGate) \
                        and op.gate.exponent == 1.0:
                    qflex_gate = "h"
                elif isinstance(op.gate, cirq.ops.XPowGate) \
                        and op.gate.exponent == 0.5:
                    qflex_gate = "x_1_2"
                elif isinstance(op.gate, cirq.ops.YPowGate) \
                        and op.gate.exponent == 0.5:
                    qflex_gate = "y_1_2"
                elif isinstance(op.gate, cirq.ops.ZPowGate) \
                        and op.gate.exponent == 0.25:
                    qflex_gate = "t"
                elif isinstance(op.gate, cirq.ops.PhasedXPowGate) \
                        and op.gate.phase_exponent == 0.25 \
                        and op.gate.exponent == 0.5:
                    qflex_gate = "hz_1_2"
                elif isinstance(op.gate, cirq.ops.ZPowGate):
                    qflex_gate = "rz({})".format(op.gate.exponent)
                elif isinstance(op.gate, cirqtmp.FSimGate):
                    # qFlex uses fractions of pi instead of radians
                    exponent1 = op.gate.theta / np.pi
                    exponent2 = op.gate.phi / np.pi
                    qflex_gate = "fsim({}, {})".format(exponent1, exponent2)
                else:
                    raise ValueError("{!r} No translation for ".format(op))


                # The moment is missing
                qflex_gate = "{} {} {}\n".format(mi, qflex_gate,
                                                 qub_str.strip())
                circuit_data.append(qflex_gate)

        return "".join(circuit_data)