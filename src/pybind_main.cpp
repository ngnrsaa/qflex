#include "pybind_main.h"

/**
 * Given a pybind11::iterable, return a stream.
 * @param data to create the stream from.
 * @param caller string representing the name of the iterable.
 */
std::stringstream GetStream(const py::iterable &data,
                            const std::string &caller = "") {
  std::stringstream out;
  for (const auto &line : data)
    if (py::isinstance<py::str>(line))
      out << line.cast<std::string>() << std::endl;
    else
      throw ERROR_MSG(caller, " must be a list of strings.");
  return out;
}

/**
 * Given options in pybind11::dict format, load the corresponding data to
 * data_type.
 * @param options in pybind11::dict format.
 * @param argument key corresponding to the desired option.
 * @param data object where to upload the option to.
 */
template <typename data_type>
void LoadData(const py::dict &options, const std::string &argument,
              data_type &data) {
  if (options.contains(argument.c_str())) {
    const auto &data_iterable = options[argument.c_str()];
    if (py::isinstance<py::iterable>(data_iterable))
      data.load(GetStream(data_iterable.cast<py::iterable>(), argument));
    else
      throw ERROR_MSG("'", argument, "' must be a list of strings.");
  } else if (options.contains((argument + "_filename").c_str())) {
    const auto &filename = options[(argument + "_filename").c_str()];
    if (py::isinstance<py::str>(filename))
      data.load(filename.cast<std::string>());
    else
      throw ERROR_MSG("'", argument, "_filename' must be a string.");
  }
}

/**
 * Given options in pybind11::dict format, load states.
 * @param options in pybind11::dict format.
 * @param argument key corresponding to the desired option.
 */
std::vector<std::string> LoadStates(const py::dict &options,
                                    const std::string &argument) {
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
          throw ERROR_MSG("states must be strings.");
    } else
      throw ERROR_MSG(argument, "_states must be a list of strings.");
  } else if (options.contains((argument + "_state").c_str())) {
    const auto &in_state = options[(argument + "_state").c_str()];
    if (py::isinstance<py::str>(in_state)) {
      states.push_back(in_state.cast<std::string>());
    } else
      throw ERROR_MSG("states must be strings.");
  }
  return states;
}

std::vector<std::pair<std::string, std::complex<double>>> simulate(
    const py::dict &options) {
  qflex::QflexInput input;

  try {
    // Get qFlex input
    qflex::QflexInput input;

    // Load grid, circuit and ordering
    LoadData(options, "grid", input.grid);
    LoadData(options, "circuit", input.circuit);
    LoadData(options, "ordering", input.ordering);

    // Load initial/final states
    input.initial_states = [&options, &input]() {
      std::vector<std::string> states = LoadStates(options, "initial");
      if (std::empty(states))
        states.push_back(std::string(input.circuit.num_active_qubits, '0'));
      return states;
    }();
    input.final_states = [&options, &input]() {
      std::vector<std::string> states = LoadStates(options, "final");
      if (std::empty(states))
        states.push_back(std::string(input.circuit.num_active_qubits, '0'));
      return states;
    }();

    // TODO: Not yet supported in Cirq
    if (std::size(input.initial_states) != 1 ||
        std::size(input.final_states) != 1)
      throw ERROR_MSG("Not yet supported");

    // Remove duplicate initial/final states
    input.initial_states.erase(std::unique(std::begin(input.initial_states),
                                           std::end(input.initial_states)),
                               std::end(input.initial_states));
    input.final_states.erase(std::unique(std::begin(input.final_states),
                                         std::end(input.final_states)),
                             std::end(input.final_states));

    // Define container for amplitudes
    std::vector<std::tuple<std::string, std::string, std::complex<double>>>
        amplitudes = EvaluateCircuit(input);

    return {{std::get<1>(amplitudes[0]), std::get<2>(amplitudes[0])}};

    // Gently return an error msg if exception is known. Otherwise, rethrow
    // exception.
  } catch (std::string err_msg) {
    std::cerr << err_msg << std::endl;
  } catch (const char *err_msg) {
    std::cerr << err_msg << std::endl;
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
  } catch (...) {
    std::rethrow_exception(std::current_exception());
  }

  return {};
}
