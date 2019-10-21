import tempfile

import cirq

from python.cirq_interface.qflex_virtual_device import  QFlexVirtualDevice

import python.utils as qflexutils

class QFlexCircuit():

    def __init__(self, cirq_circuit):
        # Behind the scene, this class creates a temporary file for each object
        self._file_handle = tempfile.mkstemp()

        with open(self._file_handle[0], "w") as f:
            # I do have the file handle anyway...
            print(self.translate_cirq_to_qflex(cirq_circuit), file = f)

    def __del__(self):
        # The destructor removes the temporary file

        import os

        # if open, close the file handle
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


    def is_qflex_compatible(self):

        if not isinstance(self.device, QFlexVirtualDevice):
            # The circuit was not validated against the device
            # TODO: Make it compatible? Validate, but for which grid?
            raise ValueError('{!r} is not a QFlexVirtualDevice'.format(self.device))

        return True

    def translate_cirq_to_qflex(self, cirq_circuit):

        number_qubits = qflexutils.GetNumberOfQubits(cirq_circuit)

        first_line = str(number_qubits) + "\n"
        circuit_data = [first_line]


        # Assume that the cirq_circuit has QFlexVirtualDevice
        assert(isinstance(cirq_circuit.device, QFlexVirtualDevice))

        grid_qubits = cirq_circuit.device.get_grid_qubits_as_keys()

        # Access moments which are private
        for mi, moment in enumerate(cirq_circuit._moments):
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