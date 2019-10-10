"""Generator for Random Quantum Circuits.

Example usage:
    $ python generator.py --device=rochester \
                          --single_qubit_gates=x_1_2,y_1_2,hz_1_2 \
                          --two_qubit_gate=cx \
                          --sequence=ABAC \
                          --depth=20 \
                          --seed=0 \
                          --output=rochester_abac_m20.qsim
"""

from typing import Dict, Iterator, List, NamedTuple, Optional, Sequence, Set, Tuple

import argparse
import itertools
import random


class Device:
    """Stores device-dependent configuration parameters for RQC generation.

    Device-dependent parameters of interest include:
        * number of qubits,
        * number of couples,
        * coupler activation patterns.
    """

    def __init__(self,
                 interaction_patterns: Dict[str, Set[Tuple[int, int]]]) -> None:
        qubits: Set[int] = set()
        couplers: Set[Tuple[int, int]] = set()
        for pattern in interaction_patterns.values():
            for coupler in pattern:
                qubits.add(coupler[0])
                qubits.add(coupler[1])
                couplers.add(coupler)
        self.n_qubits = len(qubits)
        self.n_couplers = len(couplers)
        self.n_patterns = len(interaction_patterns)
        self._interaction_patterns = interaction_patterns
        assert list(qubits) == list(range(self.n_qubits))
        assert len(couplers) <= self.n_qubits * (self.n_qubits - 1) / 2

    def pattern_names(self) -> Set[str]:
        """Returns set of names of coupler activation patterns."""
        return set(self._interaction_patterns.keys())

    def pattern(self, name: str) -> Set[Tuple[int, int]]:
        """Returns coupler activation pattern of a given name."""
        return self._interaction_patterns[name]


class TestDevice(Device):
    """Device for testing with three qubits and two couplers."""
    ACTIVATION_PATTERNS = {'A': {(0, 1)}, 'B': {(1, 2)}}

    def __init__(self) -> None:
        super().__init__(TestDevice.ACTIVATION_PATTERNS)


class Aspen(Device):
    """Rigetti's Aspen QPU with two octagonal rings of qubits."""
    ACTIVATION_PATTERNS = {
        'A': {(0, 7), (1, 14), (15, 8), (9, 10), (11, 12), (13, 2), (3, 4),
              (5, 6)},
        'B': {(0, 1), (14, 15), (8, 9), (10, 11), (12, 13), (2, 3), (4, 5),
              (6, 7)},
        'C': {(0, 7), (1, 2), (14, 13), (15, 8), (9, 10), (11, 12), (3, 4),
              (5, 6)},
    }

    def __init__(self) -> None:
        super().__init__(Aspen.ACTIVATION_PATTERNS)


class Rochester(Device):
    """IBM's 53-qubit Rochester device."""
    ACTIVATION_PATTERNS = {
        'A': {(0, 1), (2, 3), (4, 6), (5, 9), (7, 8), (10, 11), (12, 13),
              (14, 15), (16, 19), (18, 27), (20, 21), (23, 24), (25, 26),
              (28, 32), (30, 31), (33, 34), (36, 37), (38, 41), (39, 42),
              (43, 44), (45, 46), (47, 48), (49, 50)},
        'B': {(0, 5), (1, 2), (3, 4), (6, 13), (7, 16), (9, 10), (11, 12),
              (15, 18), (17, 23), (19, 20), (21, 22), (24, 25), (26, 27),
              (29, 36), (30, 39), (32, 33), (34, 35), (37, 38), (40, 46),
              (41, 50), (42, 43), (44, 45), (48, 49)},
        'C': {(0, 5), (1, 2), (3, 4), (7, 16), (8, 9), (11, 17), (13, 14),
              (15, 18), (19, 20), (21, 28), (22, 23), (25, 29), (26, 27),
              (30, 39), (31, 32), (34, 40), (35, 36), (37, 38), (41, 50),
              (42, 43), (44, 51), (46, 47), (48, 52)},
    }

    def __init__(self) -> None:
        super().__init__(Rochester.ACTIVATION_PATTERNS)


Gate = NamedTuple('Gate', [('name', str), ('qubits', Tuple[int, ...])])


class Cycle:
    """Circuit cycle consisting of gates applied to disjoint sets of qubits."""
    def __init__(self) -> None:
        self.gates: List[Gate] = []

    @property
    def active_qubits(self) -> Set[int]:
        """Return qubits to which gates are applied in this cycle."""
        active_qubits: Set[int] = set()
        for gate in self.gates:
            assert len(gate.qubits) in {1, 2}
            for qubit in gate.qubits:
                active_qubits.add(qubit)
        return active_qubits

    def append(self, gate: Gate) -> 'Cycle':
        """Adds gate to this cycle."""
        assert not self.active_qubits.intersection(set(gate.qubits))
        self.gates.append(gate)
        return self

    def to_qsim_lines(self) -> Sequence[str]:
        """Returns list of lines describing this cycle in qsim format."""
        def gate_to_qsim_line(gate: Gate) -> str:
            return gate.name + ' ' + ' '.join(str(q) for q in gate.qubits)

        return [gate_to_qsim_line(gate) for gate in self.gates]


class Circuit:
    """Quantum circuit."""
    def __init__(self) -> None:
        self.cycles: List[Cycle] = []
        self.n_qubits = 0

    def append(self, cycle: Cycle) -> None:
        """Adds cycle to this circuit."""
        self.cycles.append(cycle)
        self.n_qubits = max(self.n_qubits, max(cycle.active_qubits) + 1)

    def to_qsim_lines(self) -> Sequence[str]:
        """Returns list of lines describing this circuit in qsim format."""
        return [
            ' '.join((str(i), line))
            for i, cycle in enumerate(self.cycles)
            for line in cycle.to_qsim_lines()
        ]

    def save_as_qsim(self, filename: str):
        """Saves this circuit to file in qsim format."""
        with open(filename, 'w') as output_file:
            print(str(self.n_qubits), file=output_file)
            for line in self.to_qsim_lines():
                print(line, file=output_file)


class PseudoRandomGateGenerator(Iterator[str]):
    """Generator of pseudo-random single-qubit gates."""
    def __init__(self, valid_gates: Sequence[str]) -> None:
        self._valid_gates = valid_gates
        self._last_gate: Optional[str] = None

    def __next__(self) -> str:
        gates = [gate for gate in self._valid_gates if gate != self._last_gate]
        selected_gate = gates[random.randrange(len(gates))]
        self._last_gate = selected_gate
        return selected_gate


class PseudoRandomCircuitGenerator:
    """Generator of pseudo-random quantum circuits."""
    def __init__(self,
                 device: Device,
                 single_qubit_gates: Sequence[str],
                 two_qubit_gate: str) -> None:
        self._device = device
        self._single_qubit_gates = single_qubit_gates
        self._two_qubit_gate = two_qubit_gate

    def _generate_local_operations_cycle(
            self, prggs: Sequence[PseudoRandomGateGenerator]) -> Cycle:
        """Makes a cycle of pseudo-random single-qubit gates."""
        cycle = Cycle()
        for qubit in range(self._device.n_qubits):
            cycle.append(Gate(name=next(prggs[qubit]), qubits=(qubit,)))
        return cycle

    def _generate_interaction_cycle(self, pattern_name: str) -> Cycle:
        """Makes a cycle of two-qubit gates."""
        cycle = Cycle()
        for qubits in self._device.pattern(pattern_name):
            cycle.append(Gate(name=self._two_qubit_gate, qubits=qubits))
        return cycle

    def generate(self, sequence: str, depth: int, seed: int) -> Circuit:
        """Makes a pseudo-random quantum circuit."""
        prng_state = random.getstate()
        random.seed(seed)
        prggs = [
            PseudoRandomGateGenerator(self._single_qubit_gates)
            for _ in range(self._device.n_qubits)
        ]
        pattern_names = itertools.cycle(sequence)

        circuit = Circuit()
        for _ in range(depth):
            circuit.append(self._generate_local_operations_cycle(prggs))
            circuit.append(
                self._generate_interaction_cycle(next(pattern_names)))
        circuit.append(self._generate_local_operations_cycle(prggs))

        random.setstate(prng_state)
        return circuit


DEVICES = {'test': TestDevice(), 'aspen': Aspen(), 'rochester': Rochester()}


def main(args):
    """Generates an RQC and saves it to a file."""
    device = DEVICES[args.device]
    single_qubit_gates = args.single_qubit_gates.split(',')
    prcg = PseudoRandomCircuitGenerator(device, single_qubit_gates,
                                        args.two_qubit_gate)
    rqc = prcg.generate(args.sequence, args.depth, args.seed)
    rqc.save_as_qsim(args.output)


if __name__ == '__main__':
    AVAILABLE_DEVICES = ', '.join(DEVICES.keys())
    AVAILABLE_PATTERNS = '; '.join(
        str(d) + ': ' + ''.join(sorted(d.pattern_names()))
        for d in DEVICES.values())

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument(
        '--device',
        type=str,
        required=True,
        choices=DEVICES.keys(),
        help=f'device name; one of: {AVAILABLE_DEVICES}')
    arg_parser.add_argument(
        '--single_qubit_gates',
        type=str,
        required=True,
        help=
        'comma-separated list of single-qubit gates known to qsim, e.g. h, x_1_2'
    )
    arg_parser.add_argument(
        '--two_qubit_gate',
        type=str,
        required=True,
        help='name of the two-qubit gate known to qsim, e.g. cz, fs')
    arg_parser.add_argument(
        '--sequence',
        type=str,
        required=True,
        help=
        f'sequence of device-specific coupler activation patterns; available '
        f'patterns by device: {AVAILABLE_PATTERNS}'
    )
    arg_parser.add_argument(
        '--depth',
        type=int,
        required=True,
        help='number of layers of two-qubit gates')
    arg_parser.add_argument(
        '--seed',
        type=int,
        required=True,
        help=
        'random seed to initialize PRNG used in selection of single-qubit gates'
    )
    arg_parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='name of the file to write the generated RQC to in qsim format')

    main(arg_parser.parse_args())
