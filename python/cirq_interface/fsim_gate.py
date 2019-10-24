# Copyright 2019 The Cirq Developers
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Defines the fermionic simulation gate family.

This is the family of two-qubit gates that preserve excitations (number of ON
qubits), ignoring single-qubit gates and global phase. For example, when using
the second quantized representation of electrons to simulate chemistry, this is
a natural gateset because each ON qubit corresponds to an electron and in the
context of chemistry the electron count is conserved over time. This property
applies more generally to fermions, thus the name of the gate.
"""

import cmath
import math
from typing import Optional

import numpy as np

import cirq
from cirq import protocols, value
from cirq._compat import proper_repr
from cirq.ops import gate_features


@value.value_equality(approximate=True)
class FSimGate(gate_features.TwoQubitGate,
               gate_features.InterchangeableQubitsGate):
    """Fermionic simulation gate family.

    Contains all two qubit interactions that preserve excitations, up to
    single-qubit rotations and global phase.

    The unitary matrix of this gate is:

        [[1, 0, 0, 0],
         [0, a, b, 0],
         [0, b, a, 0],
         [0, 0, 0, c]]

    where:

        a = cos(theta)
        b = -i·sin(theta)
        c = exp(-i·phi)

    Note the difference in sign conventions between FSimGate and the
    ISWAP and CZPowGate:

        FSimGate(θ, φ) = ISWAP**(-2θ/π) CZPowGate(exponent=-φ/π)
    """

    def __init__(self, theta: float, phi: float):
        """
        Args:
            theta: Swap angle on the span(|01⟩, |10⟩) subspace, in radians.
                Determined by the strength and duration of the XX+YY
                interaction. Note: uses opposite sign convention to the
                iSWAP gate.
            phi: Controlled phase angle, in radians. Determines how much the
                |11⟩ state is phased. Note: uses opposite sign convention to
                the CZPowGate.
        """
        self.theta = theta
        self.phi = phi

    def _value_equality_values_(self):
        return self.theta, self.phi

    def _is_parameterized_(self):
        return cirq.is_parameterized(self.theta) or cirq.is_parameterized(
            self.phi)

    def _unitary_(self) -> Optional[np.ndarray]:
        if cirq.is_parameterized(self):
            return None
        a = math.cos(self.theta)
        b = -1j * math.sin(self.theta)
        c = cmath.exp(-1j * self.phi)
        return np.array([
            [1, 0, 0, 0],
            [0, a, b, 0],
            [0, b, a, 0],
            [0, 0, 0, c],
        ])

    def _pauli_expansion_(self) -> value.LinearDict[str]:
        if protocols.is_parameterized(self):
            return NotImplemented
        a = math.cos(self.theta)
        b = -1j * math.sin(self.theta)
        c = cmath.exp(-1j * self.phi)
        return value.LinearDict({
            'II': (1 + c) / 4 + a / 2,
            'IZ': (1 - c) / 4,
            'ZI': (1 - c) / 4,
            'ZZ': (1 + c) / 4 - a / 2,
            'XX': b / 2,
            'YY': b / 2,
        })

    def _resolve_parameters_(self, param_resolver: 'cirq.ParamResolver'
                            ) -> 'cirq.FSimGate':
        return FSimGate(
            protocols.resolve_parameters(self.theta, param_resolver),
            protocols.resolve_parameters(self.phi, param_resolver))

    def _apply_unitary_(self,
                        args: 'cirq.ApplyUnitaryArgs') -> Optional[np.ndarray]:
        if cirq.is_parameterized(self):
            return None
        if self.theta != 0:
            inner_matrix = protocols.unitary(cirq.Rx(2 * self.theta))
            oi = args.subspace_index(0b01)
            io = args.subspace_index(0b10)
            out = cirq.apply_matrix_to_slices(args.target_tensor,
                                              inner_matrix,
                                              slices=[oi, io],
                                              out=args.available_buffer)
        else:
            out = args.target_tensor
        if self.phi != 0:
            ii = args.subspace_index(0b11)
            out[ii] *= cmath.exp(-1j * self.phi)
        return out

    def _decompose_(self, qubits) -> 'cirq.OP_TREE':
        a, b = qubits
        xx = cirq.XXPowGate(exponent=self.theta / np.pi, global_shift=-0.5)
        yy = cirq.YYPowGate(exponent=self.theta / np.pi, global_shift=-0.5)
        yield xx(a, b)
        yield yy(a, b)
        yield cirq.CZ(a, b)**(-self.phi / np.pi)

    def _circuit_diagram_info_(self, args: 'cirq.CircuitDiagramInfoArgs'):
        t = _format_rads(args, self.theta)
        p = _format_rads(args, self.phi)
        return 'fsim({}, {})'.format(t, p), '#2'

    def __pow__(self, power):
        return FSimGate(cirq.mul(self.theta, power), cirq.mul(self.phi, power))

    def __repr__(self):
        return 'cirq.FSimGate(theta={}, phi={})'.format(proper_repr(self.theta),
                                                        proper_repr(self.phi))

    def _json_dict_(self):
        return protocols.obj_to_dict_helper(self, ['theta', 'phi'])


def _format_rads(args: 'cirq.CircuitDiagramInfoArgs', radians: float) -> str:
    if cirq.is_parameterized(radians):
        return str(radians)
    unit = 'π' if args.use_unicode_characters else 'pi'
    if radians == np.pi:
        return unit
    if radians == 0:
        return '0'
    if radians == -np.pi:
        return '-' + unit
    if args.precision is not None:
        quantity = '{{:.{}}}'.format(args.precision).format(radians / np.pi)
        return quantity + unit
    return repr(radians)
