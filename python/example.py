#!/usr/bin/env python3

import qflex

options = {
    'circuit_filename': '../config/circuits/rectangular_2x2_1-2-1_0.txt',
    'ordering_filename': '../config/ordering/rectangular_2x2.txt',
    'grid_filename': '../config/grid/rectangular_2x2.txt',
    'final_state': "0110"
}

#options = {
#    'circuit_filename': "../config/circuits/bristlecone_70_1-24-1_0.txt",
#    'ordering_filename': "../config/ordering/bristlecone_70.txt",
#    'grid_filename': "../config/grid/bristlecone_70.txt",
#    'final_state': "1" * 70
#}

# Get output from qFlex
print(qflex.simulate(options))

with open(options['circuit_filename']) as f:
    options['circuit'] = f.readlines()

with open(options['ordering_filename']) as f:
    options['ordering'] = f.readlines()

with open(options['grid_filename']) as f:
    options['grid'] = f.readlines()

del (options['circuit_filename'])
del (options['grid_filename'])
del (options['ordering_filename'])

print(qflex.simulate(options))
