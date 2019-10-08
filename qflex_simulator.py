from typing import Union, List, Any

import cirq
from cirq import study, schedules, ops, circuits
from qflex_virtual_device import QFlexVirtualDevice

import qflex


class QFlexSimulator(cirq.SimulatesFinalState):

    def __init__(self):
        return

    def simulate_sweep(
        self,
        program: Union[circuits.Circuit, schedules.Schedule],
        params: study.Sweepable,
        qubit_order: ops.QubitOrderOrList = ops.QubitOrder.DEFAULT,
        initial_state: Any = None,
    ) -> List['SimulationTrialResult']:

        if not isinstance(program, circuits.Circuit):
            raise ValueError('{!r} is not a Circuit'.format(program))

        if not isinstance(program.device, QFlexVirtualDevice):
            raise ValueError('{!r} is not a QFlexVirtualDevice'.format(program.device))
        else:
            print("OK device")

        # A strange way...
        amplitudes = qflex.simulate(program.device.compute_circuit_data(program),
                                    program.device.ordering_data,
                                    program.device.grid_data,
                                    program.device.sizex,
                                    program.device.sizey,
                                    2)

        # hard coded
        input_initial_state = "XXXX" + initial_state
        measurements = {}
        for amp in amplitudes:
            state = amp[0]
            amplitude = complex(amp[1])

            print(input_initial_state + " --> " + state + ": " + \
                  str(amplitude.real) + " " + str(amplitude.imag))

            # measurements[]

        params = {}
        measurements = {1: 1}

        trial_res = cirq.TrialResult(params=params, measurements=measurements)

        return [trial_res]