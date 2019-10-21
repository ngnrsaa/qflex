import tempfile
import os

import cirq

import python.cirq_interface.qflex_virtual_device as qdevice
import python.cirq_interface.qflex_order as qorder

import python.utils as qflexutils

class QFlexCircuit(cirq.Circuit):

    def __init__(self,
                 cirq_circuit,
                 device,
                 qflex_order = None):

        if device is None or not isinstance(device, qdevice.QFlexVirtualDevice):
            raise ValueError("QFlexVirtualDevice necessary for constructor!")

        if qflex_order is None:
            # No order was specified. Construct and order from a circuit
            qubits = device.get_indexed_grid_qubits()
            self._own_order = qorder.QFlexOrder(cirq_circuit = cirq_circuit,
                                                qubits = qubits)
        elif isinstance(qflex_order, qorder.QFlexOrder):
            self._own_order = qflex_order
        else:
            raise ValueError("{!r} is not of a QFlexOrder!")

        # The super constructor
        super().__init__(cirq_circuit, device)

        # Behind the scene, this class creates a temporary file for each object
        self._file_handle = tempfile.mkstemp()

        with open(self._file_handle[1], "w") as f:
            # I do have the file handle anyway...
            print(self.translate_cirq_to_qflex(self), file = f)


    @property
    def circuit_data(self):
        return self._file_handle[1]


    @property
    def ordering_data(self):
        return self._own_order._file_handle[1]


    def _resolve_parameters_(self,
                             param_resolver: cirq.study.ParamResolver):

        qflex_circuit = super()._resolve_parameters_(param_resolver)

        qflex_circuit.device = self.device
        qflex_circuit._own_order = self._own_order


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


    @property
    def circuit_data(self):
        return self._file_handle[1]


    def translate_cirq_to_qflex(self, cirq_circuit):

        number_qubits = qflexutils.GetNumberOfQubits(cirq_circuit)

        first_line = str(number_qubits) + "\n"
        circuit_data = [first_line]


        # Assume that the cirq_circuit has QFlexVirtualDevice
        assert(isinstance(cirq_circuit.device, qdevice.QFlexVirtualDevice))

        grid_qubits = cirq_circuit.device.get_grid_qubits_as_keys()

        for mi, moment in enumerate(cirq_circuit):
            for op in moment:

                qub_str = ""
                for qub in op.qubits:
                    qub_str += "{} ".format(grid_qubits[qub])

                qflex_gate = ""
                if isinstance(op.gate, cirq.ops.CZPowGate)\
                        and op.gate.exponent == 1.0:
                    qflex_gate = "cz"
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

                # The moment is missing
                qflex_gate = "{} {} {}\n".format(mi, qflex_gate,
                                                 qub_str.strip())
                circuit_data.append(qflex_gate)

        return "".join(circuit_data)