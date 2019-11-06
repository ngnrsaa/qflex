#include "pybind_main.h"

std::vector<std::pair<std::string, std::complex<double>>> simulate(
    const py::dict &options) {
  qflex::QflexInput input;

  auto GetStream = [](const py::iterable &data, const std::string &caller) {
    std::stringstream out;
    for (const auto &line : data)
      if (py::isinstance<py::str>(line))
        out << line.cast<std::string>() << std::endl;
      else
        throw caller + " must be a list of strings.";
    return out;
  };

  auto LoadData = [&GetStream](const py::dict &options,
                               const std::string &argument, auto &data) {
    if (options.contains(argument.c_str())) {
      const auto &data_iterable = options[argument.c_str()];
      if (py::isinstance<py::iterable>(data_iterable))
        data.load(GetStream(data_iterable.cast<py::iterable>(), argument));
      else
        throw "'" + argument + "' must be a list of strings.";
    } else if (options.contains((argument + "_filename").c_str())) {
      const auto &filename = options[(argument + "_filename").c_str()];
      if (py::isinstance<py::str>(filename))
        data.load(filename.cast<std::string>());
      else
        throw "'" + argument + "_filename' must be a string.";
    }
  };

  auto LoadStates = [](const py::dict &options, const std::string &argument,
                       const std::size_t num_active_qubits = 0) {
    std::vector<std::string> states;
    if (options.contains((argument + "_states").c_str())) {
      const auto &in_states = options[(argument + "_states").c_str()];
      if (py::isinstance<py::str>(in_states)) {
        states.push_back(in_states.cast<std::string>());
      } else if (py::isinstance<py::iterable>(in_states)) {
        for (const auto &x : in_states.cast<py::iterable>())
          if (py::isinstance<py::str>(x))
            states.push_back(x.cast<std::string>());
          else
            throw "states must be strings.";
      } else
        throw argument + "_states must be a list of strings.";
    } else if (options.contains((argument + "_state").c_str())) {
      const auto &in_state = options[(argument + "_state").c_str()];
      if (py::isinstance<py::str>(in_state)) {
        states.push_back(in_state.cast<std::string>());
      } else
        throw "states must be strings.";
    } else {
      if (argument == "initial")
        states.push_back(std::string(num_active_qubits, '0'));
      else
        states.push_back("");
    }
    return states;
  };

  try {
    // Get qFlex input
    qflex::QflexInput input;

    // Load grid, circuit and ordering
    LoadData(options, "grid", input.grid);
    LoadData(options, "circuit", input.circuit);
    LoadData(options, "ordering", input.ordering);
    const auto initial_states =
        LoadStates(options, "initial", input.circuit.num_active_qubits);
    const auto final_states = LoadStates(options, "final");

    // Define container for amplitudes
    std::vector<std::pair<
        std::string, std::vector<std::pair<std::string, std::complex<double>>>>>
        amplitudes;

    auto run_simulation = [&](const std::string &initial_state = "",
                              const std::string &final_state = "") {
      input.initial_state = initial_state;
      input.final_state = final_state;
      return EvaluateCircuit(&input);
    };

    // TODO:
    if (std::size(initial_states) != 1 or std::size(final_states) != 1)
      throw std::string("Not yet supported");

    for (const auto &is : initial_states) {
      for (const auto &fs : final_states) {
        amplitudes.push_back({is, run_simulation(is, fs)});
      }
    }

    if (std::size(initial_states) == 1 and std::size(final_states) == 1)
      return std::get<1>(amplitudes[0]);
    else {
      return {};
      // return amplitudes;
    }

    return {};

  } catch (std::string err_msg) {
    std::cerr << err_msg << std::endl;
    return {};
  } catch (const char *err_msg) {
    std::cerr << err_msg << std::endl;
    return {};
  } catch (...) {
    std::cerr << "Something went wrong" << std::endl;
    return {};
  }
}
