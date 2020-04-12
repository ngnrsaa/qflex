import cirq
import cirq.ops as ops

from qflexcirq.interface.qflex_grid import QFlexGrid


class QFlexVirtualDevice(cirq.Device):

    def __init__(self, qflex_grid=QFlexGrid.BRISTLECONE70):

        if isinstance(qflex_grid, str):
            self._qflex_grid = QFlexGrid(qflex_grid_strings=qflex_grid)
        elif isinstance(qflex_grid, QFlexGrid):
            self._qflex_grid = qflex_grid
        else:
            raise TypeError("{!r} is not a valid QFlexGrid".format(qflex_grid))

    @property
    def grid_data(self):
        return self._qflex_grid.temp_file_if.fullpath

    def duration_of(self, operation: 'cirq.Operation'):
        # No duration
        return 0

    def get_indexed_grid_qubits(self):
        return self._qflex_grid.get_grid_qubits()

    def get_grid_qubits_as_keys(self):
        original_dict = self.get_indexed_grid_qubits()
        return {original_dict[y]: y for y in original_dict}

    def decompose_operation(self, operation):

        # Known gate name
        if not isinstance(operation, ops.GateOperation):
            raise TypeError("{!r} is not a gate operation.".format(operation))

        for qub in operation.qubits:
            if not isinstance(qub, cirq.GridQubit):
                raise TypeError(
                    "This device accepts only GridQubit, but {}".format(qub))

        # default value
        decomposition = [operation]
        """
            Try to decompose the operation into elementary device operations
            TODO: Test how this works for different circuits
        """
        if not self.is_qflex_virt_dev_op(operation):
            decomposition = cirq.decompose(operation,
                                           keep=self.is_qflex_virt_dev_op)

        for dec in decomposition:
            if not self.is_qflex_virt_dev_op(dec):
                raise TypeError("Don't know how to work with {!r}.".format(
                    operation.gate))

        return decomposition

    def is_qflex_virt_dev_op(self, op):
        """
            This checks for the currently supported gate set
        """
        if not isinstance(op, ops.GateOperation):
            return False

        keep = False

        keep = keep or (isinstance(op.gate, ops.CZPowGate) and
                        (op.gate.exponent == 1))

        keep = keep or (isinstance(op.gate, ops.CNotPowGate) and
                        (op.gate.exponent == 1))

        keep = keep or (isinstance(op.gate, ops.HPowGate) and
                        (op.gate.exponent == 1))

        keep = keep or (isinstance(op.gate, ops.XPowGate) and
                        (op.gate.exponent == 0.5))

        keep = keep or (isinstance(op.gate, ops.YPowGate) and
                        (op.gate.exponent == 0.5))

        keep = keep or (isinstance(op.gate, ops.ZPowGate) and
                        (op.gate.exponent == 0.25))

        keep = keep or (isinstance(op.gate, ops.ZPowGate))

        keep = keep or (isinstance(op.gate, ops.PhasedXPowGate) and
                        (op.gate.exponent == 0.5) and
                        (op.gate.phase_exponent == 0.25))

        keep = keep or (isinstance(op.gate, cirq.ops.FSimGate))

        return keep

    def validate_operation(self, operation):
        if not isinstance(operation, cirq.GateOperation):
            raise ValueError(
                '{!r} is not a supported operation'.format(operation))

        if not self.is_qflex_virt_dev_op(operation):
            raise ValueError('{!r} is not a supported gate'.format(
                operation.gate))

        # Are the qubits GridQubits?
        # If no -> Error
        # TODO: Translate to GridQubits automatically?
        for qub in operation.qubits:
            if not isinstance(qub, cirq.GridQubit):
                raise ValueError('{} is not a grid qubit for gate {!r}'.format(
                    qub, operation.gate))

        # here we check connectivity?
        if len(operation.qubits) == 2:
            p, q = operation.qubits
            if not p.is_adjacent(q):
                # we could introduce automatic swap network
                raise ValueError('Non-local interaction: {}'.format(
                    repr(operation)))

    def validate_scheduled_operation(self, schedule, scheduled_operation):
        self.validate_operation(scheduled_operation.operation)

    def validate_circuit(self, circuit):
        #
        # Circuit and grid should have same number of qubits?
        # Otherwise -> Problem?
        #
        for moment in circuit:
            for operation in moment.operations:
                self.validate_operation(operation)

    def validate_schedule(self, schedule):
        for scheduled_operation in schedule.scheduled_operations:
            self.validate_scheduled_operation(schedule, scheduled_operation)
