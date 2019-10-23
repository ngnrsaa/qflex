#include "pybind_main.h"

std::vector<std::pair<std::string, std::complex<double>>> simulate(
    const py::dict &options) {
  qflex::QflexInput input;

  // Temporary streams
  std::ifstream fs_ordering_data;

  std::stringstream ss_circuit_data;
  std::stringstream ss_ordering_data;

  auto get_stream_from_vector = [](const auto &vector) {
    std::stringstream stream;
    for (const auto &v : vector) stream << v << std::endl;
    return stream;
  };

  // Get circuit
  if (options.contains("circuit_filename")) {
    input.circuit.load(options["circuit_filename"].cast<std::string>());

  } else if (options.contains("circuit")) {
    input.circuit.load(get_stream_from_vector(
        options["circuit"].cast<std::vector<std::string>>()));

  } else {
    std::cerr << "ERROR: Either --circuit_filename or --circuit must be specified."
              << std::endl;
    return {};
  }

  // Get ordering
  if (options.contains("ordering_filename")) {
    fs_ordering_data.open(options["ordering_filename"].cast<std::string>());
    if (not fs_ordering_data.good()) {
      std::cerr << "ERROR: cannot open file: "
                << options["ordering_filename"].cast<std::string>() << "."
                << std::endl;
      return {};
    }
    input.ordering_data = &fs_ordering_data;

  } else if (options.contains("ordering")) {
    input.ordering_data =
        &(ss_ordering_data = get_stream_from_vector(
              options["ordering"].cast<std::vector<std::string>>()));

  } else {
    std::cerr
        << "ERROR: Either --ordering_filename or --ordering must be specified."
        << std::endl;
    return {};
  }

  // Get grid
  if (options.contains("grid_filename")) {
    input.grid.load(options["grid_filename"].cast<std::string>());

  } else if (options.contains("grid")) {
    input.grid.load(get_stream_from_vector(
        options["grid"].cast<std::vector<std::string>>()));

  } else {
    std::cerr << "ERROR: Either --grid_filename or --grid must be specified."
              << std::endl;
    return {};
  }

  // Setting initial and final circuit states.
  if (options.contains("initial_state"))
    input.initial_state = options["initial_state"].cast<std::string>();

  if (options.contains("final_state"))
    input.final_state = options["final_state"].cast<std::string>();

  // Evaluating circuit.
  std::vector<std::pair<std::string, std::complex<double>>> amplitudes =
      qflex::EvaluateCircuit(&input);

  return amplitudes;
}
