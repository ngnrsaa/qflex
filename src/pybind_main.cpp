#include "pybind_main.h"

std::vector<std::pair<std::string, std::complex<double>>> simulate(
    const py::dict &options) {
  qflex::QflexInput input;

  // Check options for circuit
  switch (options.contains("circuit_filename") + options.contains("circuit")) {
    case 0:
    case 2:
      std::cerr
          << "ERROR: either 'circuit_filename' or 'circuit' must be specified"
          << std::endl;
      return {};
  }

  // Check options for ordering
  switch (options.contains("ordering_filename") +
          options.contains("ordering")) {
    case 0:
    case 2:
      std::cerr
          << "ERROR: either 'ordering_filename' or 'ordering' must be specified"
          << std::endl;
      return {};
  }

  // Check options for grid
  switch (options.contains("grid_filename") + options.contains("grid")) {
    case 0:
    case 2:
      std::cerr << "ERROR: either 'grid_filename' or 'grid' must be specified"
                << std::endl;
      return {};
  }

  // Temporary streams
  std::ifstream fs_ordering_data;

  // Get circuit
  if (options.contains("circuit_filename")) {
    input.circuit.load(options["circuit_filename"].cast<std::string>());
  } else {
    std::cerr << "ERROR: not yet implemented." << std::endl;
    return {};
  }

  // Get ordering
  if (options.contains("ordering_filename")) {
    input.ordering.load(options["ordering_filename"].cast<std::string>());
  } else {
    std::cerr << "ERROR: not yet implemented." << std::endl;
    return {};
  }

  // Get grid
  if (options.contains("grid_filename")) {
    input.grid.load(options["grid_filename"].cast<std::string>());
  } else {
    std::cerr << "ERROR: not yet implemented." << std::endl;
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
