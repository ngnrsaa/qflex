#include "pybind_main.h"

std::vector<std::pair<std::string, std::complex<double>>> simulate(
    std::vector<std::string> circuit_content,
    std::vector<std::string> ordering_content,
    std::vector<std::string> grid_content,
    int grid_height,
    int grid_width,
    std::string initial_state,
    std::string final_state,
    int super_cycles) {

    qflex::QflexInput input;

    input.grid.I = grid_height;
    input.grid.J = grid_width;

    input.K = super_cycles;

    // Creating streams for input files.
    std::stringstream circuit_stream;
    for (std::vector<std::string>::size_type i = 0; i < circuit_content.size();
        i++) {
        circuit_stream << circuit_content[i];
    }
    input.circuit_data = &circuit_stream;

    std::stringstream ordering_stream;
        for (std::vector<std::string>::size_type i = 0; i < ordering_content.size();
        i++) {
        ordering_stream << ordering_content[i];
    }
    input.ordering_data = &ordering_stream;

    std::stringstream grid_stream;
    for (std::vector<std::string>::size_type i = 0; i < grid_content.size();
        i++) {
        grid_stream << grid_content[i];
    }

    input.grid.load(grid_stream);

    // Setting initial and final circuit states.
    input.initial_state = initial_state;
    input.final_state_A = final_state;

    // Evaluating circuit.
    std::vector<std::pair<std::string, std::complex<double>>> amplitudes =
        qflex::EvaluateCircuit(&input);

    // Printing output.
    // For debugging purposes kept here commented
    // std::cout << "This is from C++" << std::endl;
    // for (int c = 0; c < amplitudes.size(); ++c) {
    //     const auto &state = amplitudes[c].first;
    //     const auto &amplitude = amplitudes[c].second;
    //     std::cout << "XXXX"
    //             << " --> " << state << ": " << std::real(amplitude) << " "
    //             << std::imag(amplitude) << std::endl;
    // }
    // std::cout << std::endl;

    return amplitudes;
}