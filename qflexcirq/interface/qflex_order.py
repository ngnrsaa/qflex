import cirq

import qflexcirq.interface.data_storage_interface as tmpi

from qflexcirq.ordering import order_circuit_simulation as auto_order


class QFlexOrder():

    def __init__(self,
                 qflex_order_strings=None,
                 cirq_circuit=None,
                 qubits=None):

        if (qflex_order_strings is None) and (cirq_circuit is None):
            # TODO: More serious checking
            raise ValueError("No order specified to constructor!")

        if (cirq_circuit is not None) and (not isinstance(
                cirq_circuit, cirq.Circuit)):
            raise ValueError("Order accepts only QFlexCircuits")

        ord_list = qflex_order_strings
        if cirq_circuit is not None:

            if qubits is None:
                raise ValueError("Qubits have to be specified!")

            # The device has to be QFlex
            # qubits = qflex_circuit.device.get_indexed_grid_qubits()
            # List of ordering commands
            print("Ordering is being computed from the provided circuit ...")
            ord_list = auto_order.circuit_to_ordering(
                cirq_circuit, qubit_names=sorted(qubits))
            print("... Done!")

        # Long string of ordering commands
        _local_order_string = '\n'.join([x.strip() for x in ord_list])

        # Behind the scene, this class creates a temporary file for each object
        self.temp_file_if = tmpi.DataStorageInterface()

        with open(self.temp_file_if._file_handle[1], "w") as f:
            # I do have the file handle anyway...
            print(_local_order_string, file=f)

    @property
    def order_data(self):
        return self.temp_file_if._file_handle[1]

    @staticmethod
    def from_existing_file(file_path):
        with open(file_path, "r") as f:
            lines = f.readlines()
            return QFlexOrder(qflex_order_strings=lines)
