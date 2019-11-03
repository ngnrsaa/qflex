#!/usr/bin/env python3

from python import qflex

options = {
    'circuit-filename' : 'config/circuits/rectangular_2x2_1-2-1_0.txt',
    'ordering-filename' : 'config/ordering/rectangular_2x2.txt',
    'grid-filename' : 'config/grid/rectangular_2x2.txt',
    'initial-states' : ['0000', '1111', '0011'],
    'final-states' : ['0000', '0101', '1111']
}

print(qflex.simulate(options))
