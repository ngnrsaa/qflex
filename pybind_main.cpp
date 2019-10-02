#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <omp.h>

#include <vector>

#include "evaluate_circuit.h"

// Input: ./qflex.x I J K fidelity circuit_filename ordering_filename \
//            grid_filename [initial_conf] [final_conf]
//
// Example:
// $ ./qflex.x 11 12 2 0.005 ./circuits/ben_11_16_0.txt \
//       ./ordering/bristlecone_48.txt ./grid/bristlecone_48.txt

std::vector<std::pair<std::string, std::complex<double>>> simulate(
    std::vector<std::string> circuit_content,
    std::vector<std::string> ordering_content,
    std::vector<std::string> grid_content) {
  // Reading input.
  // if (argc < 8) throw std::logic_error("ERROR: Not enough arguments.");
  int current_arg = 1;
  qflex::QflexInput input;

  /*
      TODO: for the moment, hard coded
  */
  input.grid_height = 11;
  input.grid_width = 12;
  input.super_cycles = 2;
  input.fidelity = 0.005;

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
  input.grid_data = &grid_stream;

  /*
      TODO: for the moment, the following two are not sent as parameters
  */
  // Setting initial and final circuit states.
  // if (argc > current_arg) {
  //     input.initial_state = std::string(argv[current_arg++]);
  // }
  // if (argc > current_arg) {
  //     input.final_state_A = std::string(argv[current_arg++]);
  // }

  // Evaluating circuit.
  std::vector<std::pair<std::string, std::complex<double>>> amplitudes =
      qflex::EvaluateCircuit(&input);

  // Printing output.
  std::cout << "This is from C++" << std::endl;
  for (int c = 0; c < amplitudes.size(); ++c) {
    const auto &state = amplitudes[c].first;
    const auto &amplitude = amplitudes[c].second;
    std::cout << "XXXX"
              << " --> " << state << ": " << std::real(amplitude) << " "
              << std::imag(amplitude) << std::endl;
  }
  std::cout << std::endl;

  return amplitudes;
}

PYBIND11_MODULE(qflex, m) {
  m.doc() = "pybind11 example plugin";  // optional module docstring

  m.def("simulate", &simulate, "A function which adds two numbers");
}
