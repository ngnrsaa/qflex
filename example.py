#!/usr/bin/env python3

from python import qflex

options = {
    'circuit_filename': 'config/circuits/rectangular_2x2_1-2-1_0.txt',
    'ordering_filename': 'config/ordering/rectangular_2x2.txt',
    'grid_filename': 'config/grid/rectangular_2x2.txt',
    'initial_states': ['0000', '1111', '0011'],
    'final_states': ['0000', '0101', '1111']
}

print(qflex.simulate(options))

options = {
    'circuit_filename': 'config/circuits/rectangular_2x2_1-2-1_0.txt',
    'ordering_filename': 'config/ordering/rectangular_2x2.txt',
    'grid_filename': 'config/grid/rectangular_2x2.txt',
    'initial_state': '0011',
    'final_state': '0000'
}

print(qflex.simulate(options))
