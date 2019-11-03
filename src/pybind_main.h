#ifndef __PYBIND_MAIN
#define __PYBIND_MAIN

#include <Python.h>
#include "PyCPP/py.h"
#include "PyCPP/py_parse.h"
#include "PyCPP/py_build.h"
#include "input.h"

namespace qflex {
  extern std::vector<std::pair<std::string, std::complex<double>>> EvaluateCircuit(QflexInput*);
}

inline PyObject* simulate(PyObject *, PyObject *q);

static PyMethodDef methods[] = {
    { "simulate", (PyCFunction)simulate, METH_O, "Simulate circuit using qFlex." },
    { nullptr, nullptr, 0, nullptr }
};

static PyModuleDef module = {
    PyModuleDef_HEAD_INIT,
    "qFlex",
    "qFlex porting to Python.",
    0,
    methods
};

PyMODINIT_FUNC PyInit_qflex() {
    return PyModule_Create(&module);
}

#endif
