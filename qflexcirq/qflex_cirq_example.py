# The interface between Cirq and the Python interface to the C++ QFlex
import sys, os
sys.path.insert(
    1, os.path.realpath(os.path.dirname(os.path.realpath(__file__)) + '/../'))
import qflexcirq.interface.qflex_simulator as qsim
import qflexcirq.interface.qflex_virtual_device as qdevice
import qflexcirq.interface.qflex_grid as qgrid
import qflexcirq.interface.qflex_circuit as qcirc
import qflexcirq.interface.qflex_order as qorder

import qflexcirq.utils as qflexutils

from qflexcirq import qflex
"""
Example HOWTO

Running the example requires a compiled version of QFlex.

Two possibilities are available:
1. compilation from a clone of the github repository
2. pip installing the qflexcirq package


In the examples below different simulation configurations are illustrated. 
A configuration includes three types of information: grid, order of tensor 
contractions, quantum circuit. For example, config_small or config_large use 
different grid  arrangements, different circuits and different input states. 

When using the pip install version of qflexcirq, the 
   !!! config files are not copied onto the  local machine !!!

Config files have to be explicitly downloaded from the
github repository. Correspondingly, the paths in the configurations below need
to be adapted to where the config files are stored.
"""

config_small = {
    'circuit_filename': 'config/circuits/rectangular_2x2_1-2-1_0.txt',
    'ordering_filename': 'config/ordering/rectangular_2x2.txt',
    'grid_filename': 'config/grid/rectangular_2x2.txt',
    'final_state': "0110"
}

config_mem_crash = {
    'circuit_filename': "config/circuits/bristlecone_70_1-40-1_0.txt",
    'ordering_filename': "config/ordering/bristlecone_70.txt",
    'grid_filename': 'config/grid/bristlecone_70.txt',
    'final_state': "1" * 70
}

config_large = {
    'circuit_filename': "config/circuits/bristlecone_70_1-16-1_0.txt",
    'ordering_filename': "config/ordering/bristlecone_70.txt",
    'grid_filename': 'config/grid/bristlecone_70.txt',
    'final_state': "1" * 70
}

config_sycamore = {
    'circuit_filename': "config/circuits/sycamore_53_4_0.txt",
    'ordering_filename': "config/ordering/sycamore_53.txt",
    'grid_filename': 'config/grid/sycamore_53.txt',
    'final_state': "1" * 53
}


def run_qflex_simulator(config):

    my_grid = qgrid.QFlexGrid.from_existing_file(config['grid_filename'])
    my_device = qdevice.QFlexVirtualDevice(qflex_grid=my_grid)

    # The qubits are collected and indexed from the underlying grid_string
    # that was passed as constructor to the Device
    my_qubits = my_device.get_indexed_grid_qubits()

    # Take a QFlex circuit and generate a Cirq circuit from it
    # The Cirq circuit will be afterwards transformed into a Qflex circuit
    # You can construct a Cirq circuit from an existing QFlex circuit
    # Note that circuits provided in files were designed for a specific arrangement
    my_circuit = qflexutils.GetCircuitOfMoments(config["circuit_filename"],
                                                my_qubits)

    my_order = qorder.QFlexOrder.from_existing_file(config["ordering_filename"])

    circuit_on_device = qcirc.QFlexCircuit(cirq_circuit=my_circuit,
                                           device=my_device,
                                           qflex_order=my_order)

    print("\nRunning QFlex simulation\n")

    my_sim = qsim.QFlexSimulator()
    myres = my_sim.compute_amplitudes(program=circuit_on_device,
                                      bitstrings=[config['final_state']])
    print(myres)


def run_pybind_interface(config):
    print("\nRunning Pybind Interface\n")
    print(qflex.simulate(config))


def main():
    #
    print("\n\n  === Simulation 1" + str(config_small))
    run_qflex_simulator(config_small)
    run_pybind_interface(config_small)

    print("\n\n  === Simulation 2" + str(config_large))
    run_qflex_simulator(config_large)
    run_pybind_interface(config_large)

    print("\n\n  === Simulation 3" + str(config_sycamore))
    run_qflex_simulator(config_sycamore)
    run_pybind_interface(config_sycamore)

    #
    # TODO: This simulation fails due to insufficient memory
    #
    # print("  === Simulation 3" + str(config_mem_crash))
    # run_qflex_simulator(config_mem_crash)
    # run_pybind_interface(config_mem_crash)


if __name__ == "__main__":
    main()
