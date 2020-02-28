from typing import Sequence

from cirq import study, ops, circuits, SimulatesAmplitudes

from qflexcirq import qflex

import qflexcirq.interface.qflex_virtual_device as qdevice
import qflexcirq.interface.qflex_circuit as qcirc


class QFlexSimulator(SimulatesAmplitudes):

    def __init__(self):
        return

    def compute_amplitudes_sweep(
        self,
        program: circuits.Circuit,
        bitstrings: Sequence[int],
        params: study.Sweepable,
        qubit_order: ops.QubitOrderOrList = ops.QubitOrder.DEFAULT,
    ) -> Sequence[Sequence[complex]]:

        if not isinstance(program, qcirc.QFlexCircuit):
            raise ValueError('{!r} is not a QFlexCircuit'.format(program))

        if not isinstance(program.device, qdevice.QFlexVirtualDevice):
            # The circuit was not validated against the device
            # TODO: Make it compatible? Validate, but for which grid?
            raise ValueError('{!r} is not a QFlexVirtualDevice'.format(
                program.device))

        param_resolvers = study.to_resolvers(params)

        trials_results = []
        for prs in param_resolvers:

            from cirq import protocols
            # The QFlexVirtualDevice device is "inherited" from the original program
            solved_circuit = protocols.resolve_parameters(program, prs)

            sweep_return = []
            # Not sure what these params could look like...for the moment
            for bitstring in bitstrings:

                options = {
                    'circuit_filename': solved_circuit.circuit_data,
                    'ordering_filename': solved_circuit.ordering_data,
                    'grid_filename': program.device.grid_data,
                    'final_state': bitstring
                }

                amplitudes = qflex.simulate(options)

                for amp in amplitudes:
                    amplitude = complex(amp[1])

                    # For debugging purposes commented the following
                    # state = amp[0]
                    # print(input_initial_state + " --> " + state + ": " + \
                    #       str(amplitude.real) + " " + str(amplitude.imag))

                    # the amplitude for bitstring is stored
                    sweep_return.append(amplitude)

            trials_results.append(sweep_return)

        return trials_results
