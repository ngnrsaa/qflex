import cirq
import cirq.ops as ops

class QFlexVirtualDevice(cirq.Device):

    def __init__(self, arrangement = None):

        # default value
        self._arrangement = _BRISTLECONE48

        if arrangement is not None:
            self._arrangement = arrangement

        lines = self._arrangement.split("\n")

        self._sizex = len(lines)
        self._sizey = len(lines[0].strip())

        self._qubits = []

        for row, line in enumerate(lines):
            for col, c in enumerate(line.strip()):
                if c == "1":
                    self._qubits.append(cirq.GridQubit(row, col))

        self._circuit_data = []

    def duration_of(self, operation: 'cirq.Operation'):
        # No duration
        return 0

    @property
    def sizex(self):
        return self._sizex

    @property
    def sizey(self):
        return self._sizey

    def compute_index(self, row, col):
        return row * self.sizex + col


    def compute_circuit_data(self, program):
        """
        Decomposition in the device worked fine and no errors were raised
        Perform translation into QFlex instructions

        :param program: Decomposed circuit
        :return: Lists of strings to be sent to the C++ QFlex Simulator
        """
        # TODO: Hard coded number of qubits in the circuit
        circuit_data = ["1\n"]

        # Access moment in unorthodox manner?
        for mi, moment in enumerate(program._moments):
            for op in moment:

                qub_str = ""
                for qub in op.qubits:
                    qub_str += "{} ".format(self.compute_index(qub.row, qub.col))

                qflex_gate = ""
                if isinstance(op.gate, ops.CZPowGate):
                    qflex_gate = "cz"
                elif isinstance(op.gate, ops.HPowGate):
                    qflex_gate = "h"
                elif isinstance(op.gate, ops.XPowGate):
                    qflex_gate = "x_1_2"
                elif isinstance(op.gate, ops.YPowGate):
                    qflex_gate = "y_1_2"
                elif isinstance(op.gate, ops.ZPowGate):
                    qflex_gate = "t"

                # The moment is missing
                qflex_gate = "{} {} {}\n".format(mi, qflex_gate, qub_str)
                circuit_data.append(qflex_gate)

        # circuit_data = []
        # with open("circuits/ben_11_16_0.txt", "r") as file:
        #     circuit_data = file.readlines()

        return circuit_data

    @property
    def grid_data(self):
        gdata = self._arrangement.replace("0", "0 ")\
                        .replace("1", "1 ")

        grid_data = [x.strip()+"\n" for x in gdata.split("\n")]

        return grid_data

    @property
    def ordering_data(self):
        """
            For the moment we have hardcoded ordering if the device is a grid
            with 48 or 70 qubits
        """
        file_name = "no_file"

        if self._arrangement == _BRISTLECONE48:
            file_name = "ordering/bristlecone_48.txt"
        elif self._arrangement == _BRISTLECONE70:
            file_name = "ordering/bristlecone_70.txt"

        lines = []
        if file_name != "no_file":
            with open(file_name, "r") as file:
                lines = file.readlines()
        else:
            # TODO: Use Python lib for ordering heuristic
            lines = ["TODO: heuristic based list"]

        return lines

    def decompose_operation(self, operation):

        # Known gate name
        if not isinstance(operation, ops.GateOperation):
            raise TypeError("{!r} is not a gate operation.".format(operation))

        for qub in operation.qubits:
            if not isinstance(qub, cirq.GridQubit):
                raise TypeError("This device accepts only GridQubit, but {}"
                                .format(qub))

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
                raise TypeError("Don't know how to work with {!r}."
                                .format(operation.gate))

        return decomposition

    def is_qflex_virt_dev_op(self, op):
        """
            For the moment allow only CZ, sqrt(X), sqrt(Y), T, H
        """
        if not isinstance(op, ops.GateOperation):
            return False

        keep = True

        keep = keep or (isinstance(op.gate, ops.CZPowGate)
                        and
                        (op.gate.exponent == 1))

        keep = keep or (isinstance(op.gate, ops.HPowGate)
                        and
                        (op.gate.exponent == 1))

        keep = keep or (isinstance(op.gate, ops.XPowGate)
                        and
                        (op.gate.exponent == 0.5))

        keep = keep or (isinstance(op.gate, ops.YPowGate)
                        and
                        (op.gate.exponent == 0.5))

        keep = keep or (isinstance(op.gate, ops.ZPowGate)
                        and
                        (op.gate.exponent == 0.25))

        return keep


    def validate_operation(self, operation):
        if not isinstance(operation, cirq.GateOperation):
            raise ValueError('{!r} is not a supported operation'.format(operation))

        if not self.is_qflex_virt_dev_op(operation):
            raise ValueError('{!r} is not a supported gate'.format(operation.gate))

        # Are the qubits GridQubits?
        # If no -> Error
        # TODO: Translate to GridQubits automatically?
        for qub in operation.qubits:
            if not isinstance(qub, cirq.GridQubit):
                raise ValueError('{} is not a grid qubit for gate {!r}'.format(qub, operation.gate))

        # here we check connectivity?
        if len(operation.qubits) == 2:
            p, q = operation.qubits
            if not p.is_adjacent(q):
                # we could introduce automatic swap network
                raise ValueError('Non-local interaction: {}'.format(repr(operation)))

    def validate_scheduled_operation(self, schedule, scheduled_operation):
        self.validate_operation(scheduled_operation.operation)

    def validate_circuit(self, circuit):
        for moment in circuit:
            for operation in moment.operations:
                self.validate_operation(operation)

    def validate_schedule(self, schedule):
        for scheduled_operation in schedule.scheduled_operations:
            self.validate_scheduled_operation(schedule, scheduled_operation)


_BRISTLECONE48 = """000001100000
                    000011110000
                    000111111000
                    001111111100
                    001111111100
                    001111111100
                    000111111000
                    000011110000
                    000001100000
                    000000000000
                    000000000000"""#11 lines of 12 cols

_BRISTLECONE70 = """000001100000
                    000011110000
                    000111111000
                    001111111100
                    011111111110
                    011111111110
                    011111111110
                    001111111100
                    000111111000
                    000011110000
                    000001100000"""#11 lines of 12 cols
