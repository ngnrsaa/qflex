#!/usr/bin/env python3

import qflex

# simulation = {
#     'circuit_filename': '../config/circuits/rectangular_2x2_1-2-1_0.txt',
#     'ordering_filename': '../config/ordering/rectangular_2x2.txt',
#     'grid_filename': '../config/grid/rectangular_2x2.txt',
#     'final_state': "0110"
# }
#
# print(qflex.simulate(simulation))

options = {
    'circuit_filename': "../config/circuits/bristlecone_70_1-40-1_0.txt",
    'ordering_filename': "../config/ordering/bristlecone_70.txt",
    'grid_filename': "../config/grid/bristlecone_70.txt",
    'final_state': "1" * 62
}

# Get output from qFlex
print(qflex.simulate(options))
