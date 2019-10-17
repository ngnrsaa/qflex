#ifndef __PYBIND_MAIN
#define __PYBIND_MAIN

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <omp.h>

#include <vector>

#include "evaluate_circuit.h"

std::vector<std::pair<std::string, std::complex<double>>> simulate(
    std::vector<std::string> circuit_content,
    std::vector<std::string> ordering_content,
    std::vector<std::string> grid_content,
    int grid_height,
    int grid_width,
    std::string initial_state,
    std::string final_state,
    int super_cycles);


PYBIND11_MODULE(qflex, m) {
  m.doc() = "pybind11 plugin";  // optional module docstring

  m.def("simulate", &simulate, "Call the simulator with the normal parameters");
}

#endif