#include <omp.h>

#include "evaluate_circuit.h"

// Input: ./qflex.x I J K fidelity circuit_filename ordering_filename \
//            grid_filename [initial_conf] [final_conf]
//
// Example:
// $ ./qflex.x 11 12 2 0.005 ./circuits/ben_11_16_0.txt \
//       ./ordering/bristlecone_48.txt ./grid/bristlecone_48.txt
int main(int argc, char** argv) {
  // Reading input.
  if (argc < 8) throw std::logic_error("ERROR: Not enough arguments.");
  int current_arg = 1;
  qflex::QflexInput input;
  input.I = atoi(argv[current_arg++]);
  input.J = atoi(argv[current_arg++]);
  input.K = atoi(argv[current_arg++]);
  input.fidelity = atof(argv[current_arg++]);

  // Creating streams for input files.
  std::string circuit_filename = std::string(argv[current_arg++]);
  auto circuit_data = std::ifstream(circuit_filename);
  if (!circuit_data.good()) {
    std::cout << "Cannot open circuit data file: " << circuit_filename
              << std::endl;
    assert(circuit_data.good());
  }
  input.circuit_data = &circuit_data;
  std::string ordering_filename = std::string(argv[current_arg++]);
  auto ordering_data = std::ifstream(ordering_filename);
  if (!ordering_data.good()) {
    std::cout << "Cannot open ordering data file: " << ordering_filename
              << std::endl;
    assert(ordering_data.good());
  }
  input.ordering_data = &ordering_data;
  std::string grid_filename = std::string(argv[current_arg++]);
  auto grid_data = std::ifstream(grid_filename);
  if (!grid_data.good()) {
    std::cout << "Cannot open grid data file: " << grid_filename << std::endl;
    assert(grid_data.good());
  }
  input.grid_data = &grid_data;

  // Setting initial and final circuit states.
  if (argc > current_arg) {
    input.initial_state = std::string(argv[current_arg++]);
  }
  if (argc > current_arg) {
    input.final_state_A = std::string(argv[current_arg++]);
  }

  // Evaluating circuit.
  std::vector<std::pair<std::string, std::complex<double>>> amplitudes =
      qflex::EvaluateCircuit(&input);

  // Printing output.
  for (int c = 0; c < amplitudes.size(); ++c) {
    const auto& state = amplitudes[c].first;
    const auto& amplitude = amplitudes[c].second;
    std::cout << input.initial_state << " --> " << state << ": "
              << std::real(amplitude) << " " << std::imag(amplitude)
              << std::endl;
  }
  return 0;
}
