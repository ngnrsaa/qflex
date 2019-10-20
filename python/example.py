#!/usr/bin/env python3

import qflex

simulation = {
    'circuit_filename': '../config/circuits/rectangular_2x2_1-2-1_0.txt',
    'ordering_filename': '../config/ordering/rectangular_2x2.txt',
    'grid_filename': '../config/grid/rectangular_2x2.txt',
    'final_state': "0110"
}

print(qflex.simulate(simulation))
