from typing import Union, Sequence

from cirq import study, schedules, ops, circuits

import qflex

from cirqinterface.cirq_amplitudes_sim import SimulatesAmplitudes
from cirqinterface.qflex_virtual_device import QFlexVirtualDevice

class QFlexSimulator(SimulatesAmplitudes):

    def __init__(self):
        return

    def compute_amplitudes_sweep(
            self,
            program: Union[circuits.Circuit, schedules.Schedule],
            bitstrings: Sequence[int],
            params: study.Sweepable,
            qubit_order: ops.QubitOrderOrList = ops.QubitOrder.DEFAULT,
    ) -> Sequence[Sequence[complex]]:

        if not isinstance(program, circuits.Circuit):
            raise ValueError('{!r} is not a Circuit'.format(program))

        # circuit = (program if isinstance(program, circuits.Circuit)
        #            else program.to_circuit())

        if not isinstance(program.device, QFlexVirtualDevice):
            raise ValueError('{!r} is not a QFlexVirtualDevice'.format(program.device))
        else:
            print("The circuits's device is a QFlexVirtualDevice...OK")

        param_resolvers = study.to_resolvers(params)

        trials_results = []
        for prs in param_resolvers:

            from cirq import protocols
            solved_circuit = protocols.resolve_parameters(program, prs)

            sweep_return = []
            # Not sure what these params could look like...for the moment
            for bitstring in bitstrings:

                amplitudes = qflex.simulate(
                    program.device.compute_circuit_data(solved_circuit),
                    program.device.ordering_data,
                    program.device.grid_data,
                    program.device.sizex,
                    program.device.sizey,
                    "0" * 70,
                    "1" * 62,
                    2)

                # hard coded
                input_initial_state = "XXXX" #000000000

                for amp in amplitudes:
                    state = amp[0]
                    amplitude = complex(amp[1])

                    print(input_initial_state + " --> " + state + ": " + \
                          str(amplitude.real) + " " + str(amplitude.imag))

                    # the amplitude for bitstring
                    sweep_return.append(amplitude)

            trials_results.append(sweep_return)

        return trials_results
