# qFlex Directory Structure

Libraries in the qFlex repository are roughly organized based on language and
the type of data that each is responsible for handling.

## src

Contains the core behavior of qFlex, written in C++.

- **evaluate_circuit** is the top-level library which provides APIs for circuit
  simulation. It also handles parsing of "grid" files.
- **contraction_utils** is responsible for the step-by-step contraction of the
  tensor grid, as well as parsing the "ordering" files which specify the order
  of these contractions.
- **read_circuit** handles both parsing the "circuit" file and collapsing the
  3D tensor network down to a 2D grid.
- **tensor** contains the qFlex "Tensor" class, which is the basis for all
  operations performed in qFlex.
- **pybind_main** uses [pybind](https://github.com/pybind/pybind11) to
  provide a Python wrapper for qFlex. The resulting module is stored under
  [qflex/qflexcirq](/qflexcirq).

## config

Contains an assortment of ready-made input files for qFlex.

- **circuits** defines the sequence of gates to be simulated.
- **ordering** specifies the ordering of the simulation steps.
- **grid** maps the qubits of a device to a square lattice.
- **patterns** contains sets of 2-qubit pairs for various devices. These
  files are used for generating random quantum circuits.

## docs

Contains documentation for building and running qFlex.

- **nasa-cla** contains the Contributing License Agreement for qFlex.

## python

Contains useful Python scripts and integration with other Python-based tools.

- **cirq_interface** contains the qFlex-[Cirq](https://github.com/quantumlib/cirq)
  integration layer.
- **circuits** contains utility scripts for generating qFlex circuits.
- **ordering** provides a tool for generating qFlex "ordering" files from a Cirq
  "Circuit" object.

## scripts

Contains various shell scripts, mostly oriented towards qFlex contributors.

## tests

Contains test code for other qFlex libraries. Subdirectories are organized by
language and the library under test.
