
from typing import (
    Any,
    Dict,
    Hashable,
    Iterator,
    List,
    Sequence,
    Tuple,
    Union,
    Optional,
)

import abc

import collections
import numpy as np

from cirq import circuits, ops, protocols, schedules, study, value


class SimulatesAmplitudes(metaclass=abc.ABCMeta):
    """Simulator that computes final amplitudes of given bitstrings.
    Given a circuit and a list of bitstrings, computes the amplitudes
    of the given bitstrings in the state obtained by applying the circuit
    to the all zeros state. Implementors of this interface should implement
    the compute_amplitudes_sweep method.
    """

    def compute_amplitudes(
            self,
            program: Union[circuits.Circuit, schedules.Schedule],
            bitstrings: Sequence[int],
            param_resolver: 'study.ParamResolverOrSimilarType' = None,
            qubit_order: ops.QubitOrderOrList = ops.QubitOrder.DEFAULT,
    ) -> Sequence[complex]:
        """Computes the desired amplitudes.
        The initial state is assumed to be the all zeros state.
        Args:
            program: The circuit or schedule to simulate.
            bitstrings: The bitstrings whose amplitudes are desired, input
                as an integer array where each integer is formed from measured
                qubit values according to `qubit_order` from most to least
                significant qubit, i.e. in big-endian ordering.
            param_resolver: Parameters to run with the program.
            qubit_order: Determines the canonical ordering of the qubits. This
                is often used in specifying the initial state, i.e. the
                ordering of the computational basis states.
        Returns:
            List of amplitudes.
        """
        return self.compute_amplitudes_sweep(
            program, bitstrings, study.ParamResolver(param_resolver),
            qubit_order)[0]

    @abc.abstractmethod
    def compute_amplitudes_sweep(
            self,
            program: Union[circuits.Circuit, schedules.Schedule],
            bitstrings: Sequence[int],
            params: study.Sweepable,
            qubit_order: ops.QubitOrderOrList = ops.QubitOrder.DEFAULT,
    ) -> Sequence[Sequence[complex]]:
        """Computes the desired amplitudes.
        The initial state is assumed to be the all zeros state.
        Args:
            program: The circuit or schedule to simulate.
            bitstrings: The bitstrings whose amplitudes are desired, input
                as an integer array where each integer is formed from measured
                qubit values according to `qubit_order` from most to least
                significant qubit, i.e. in big-endian ordering.
            params: Parameters to run with the program.
            qubit_order: Determines the canonical ordering of the qubits. This
                is often used in specifying the initial state, i.e. the
                ordering of the computational basis states.
        Returns:
            List of lists of amplitudes. The outer dimension indexes the
            circuit parameters and the inner dimension indexes the bitstrings.
        """
        raise NotImplementedError()