#!/usr/bin/env python3

import qflex

simulation = {
    'circuit_filename': '../config/circuits/rectangular_2x2_1-2-1_0.txt',
    'ordering_filename': '../config/ordering/rectangular_2x2.txt',
    'grid_filename': '../config/grid/rectangular_2x2.txt',
    'final_state': "0110"
}

print(qflex.simulate(simulation))

with open(simulation['circuit_filename']) as f:
    simulation['circuit'] = f.readlines()

with open(simulation['ordering_filename']) as f:
    simulation['ordering'] = f.readlines()

with open(simulation['grid_filename']) as f:
    simulation['grid'] = f.readlines()

del (simulation['circuit_filename'])
del (simulation['grid_filename'])
del (simulation['ordering_filename'])

print(qflex.simulate(simulation))
