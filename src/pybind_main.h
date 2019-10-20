#ifndef __PYBIND_MAIN
#define __PYBIND_MAIN

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include <vector>

#include "evaluate_circuit.h"
#include "read_circuit.h"

std::vector<std::pair<std::string, std::complex<double>>> simulate(const py::dict &options);

PYBIND11_MODULE(qflex, m) {
  m.doc() = "pybind11 plugin";  // optional module docstring

  m.def("simulate", &simulate, "Call the simulator with the normal parameters");
}

#endif
