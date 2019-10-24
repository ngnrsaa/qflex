# The interface between Cirq and the Python interface to the C++ QFlex
import python.cirq_interface.qflex_simulator as qsim
import python.cirq_interface.qflex_virtual_device as qdevice
import python.cirq_interface.qflex_grid as qgrid
import python.cirq_interface.qflex_circuit as qcirc
import python.cirq_interface.qflex_order as qorder

import python.utils as qflexutils

from python import qflex

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
    my_device = qdevice.QFlexVirtualDevice(qflex_grid = my_grid)

    # The qubits are collected and indexed from the underlying grid_string
    # that was passed as constructor to the Device
    my_qubits = my_device.get_indexed_grid_qubits()

    # Take a QFlex circuit and generate a Cirq circuit from it
    # The Cirq circuit will be afterwards transformed into a Qflex circuit
    # You can construct a Cirq circuit from an existing QFlex circuit
    # Note that circuits provided in files were designed for a specific arrangement
    my_circuit = qflexutils.GetCircuitOfMoments(config["circuit_filename"], my_qubits)

    my_order = qorder.QFlexOrder.from_existing_file(config["ordering_filename"])

    circuit_on_device = qcirc.QFlexCircuit(cirq_circuit = my_circuit,
                                     device = my_device,
                                     qflex_order = my_order)

    print("\nRunning QFlex simulation\n")

    my_sim = qsim.QFlexSimulator()
    myres = my_sim.compute_amplitudes(program = circuit_on_device,
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

