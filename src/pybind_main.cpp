#include "pybind_main.h"

/**
 * Given a PyCPP::iterable, return a stream.
 * @param data to create the stream from.
 * @param caller string representing the name of the iterable.
 */
template <typename data_type>
std::stringstream GetStream(const data_type &data,
                            const std::string &caller = "") {
  std::stringstream out;
  for (const auto &line : data)
    if (PyCPP::Is<PyCPP::String>(line))
      out << line << std::endl;
    else
      throw ERROR_MSG(caller, " must be a list of strings.");
  return out;
}

/**
 * Given options in PyCPP::Map format, load the corresponding data to
 * data_type.
 * @param options in PyCPP::Map format.
 * @param argument key corresponding to the desired option.
 * @param data object where to upload the option to.
 */
template <typename options_type, typename data_type>
void LoadData(const options_type &options, const std::string &argument,
              data_type &data) {
  if (auto w = options.find(PyCPP::String(argument)); w != std::end(options)) {
    const auto &data_iterable = std::get<1>(*w);
    if (PyCPP::Is<PyCPP::Sequence>(data_iterable))
      data.load(GetStream(PyCPP::As<PyCPP::Sequence>(data_iterable), argument));
    else
      throw ERROR_MSG("'", argument, "' must be a list of strings.");
  } else if (auto w = options.find(PyCPP::String(argument + "_filename"));
             w != std::end(options)) {
    const auto &filename = std::get<1>(*w);
    if (PyCPP::Is<PyCPP::String>(filename))
      data.load(PyCPP::As<PyCPP::String>(filename));
    else
      throw ERROR_MSG("'", argument, "_filename' must be a string.");
  }
}

/**
 * Given options in pybind11::dict format, load states.
 * @param options in pybind11::dict format.
 * @param argument key corresponding to the desired option.
 */
template <typename options_type>
std::vector<std::string> LoadStates(const options_type &options,
                                    const std::string &argument) {
  std::vector<std::string> states;
  if (auto w = options.find(PyCPP::String(argument + "_states"));
      w != std::end(options)) {
    const auto &in_states = std::get<1>(*w);
    if (PyCPP::Is<PyCPP::String>(in_states)) {
      states.push_back(PyCPP::As<PyCPP::String>(in_states));
    } else if (PyCPP::Is<PyCPP::Sequence>(in_states)) {
      for (const auto &x : PyCPP::As<PyCPP::Sequence>(in_states))
        if (PyCPP::Is<PyCPP::String>(x))
          states.push_back(PyCPP::As<PyCPP::String>(x));
        else
          throw ERROR_MSG("states must be strings.");
    } else
      throw ERROR_MSG(argument, "_states must be a list of strings.");
  } else if (auto w = options.find(PyCPP::String(argument + "_state"));
             w != std::end(options)) {
    const auto &in_state = std::get<1>(*w);
    if (PyCPP::Is<PyCPP::String>(in_state)) {
      states.push_back(PyCPP::As<PyCPP::String>(in_state));
    } else
      throw ERROR_MSG("states must be strings.");
  } else
    states.push_back("");
  return states;
}

inline PyObject *simulate(PyObject *, PyObject *arg) {
  try {
    // The passed argument must be a dictionary
    const auto v_options = PyCPP::Parse(arg);
    if (not PyCPP::Is<PyCPP::Map>(v_options[0])) {
      throw ERROR_MSG(
          "Simulation takes only one argument and it must be a dictionary");
    }
    const PyCPP::Map &options = PyCPP::As<PyCPP::Map>(v_options[0]);

    // Get qFlex input
    qflex::QflexInput input;

    // Load grid, circuit and ordering
    LoadData(options, "grid", input.grid);
    LoadData(options, "circuit", input.circuit);
    LoadData(options, "ordering", input.ordering);
    const auto initial_states = LoadStates(options, "initial");
    const auto final_states = LoadStates(options, "final");

    // Define container for amplitudes
    std::vector<std::pair<
        std::string, std::vector<std::pair<std::string, std::complex<double>>>>>
        amplitudes;

    for (const auto &is : initial_states) {
      for (const auto &fs : final_states) {
        input.initial_state = is;
        input.final_state = fs;
        amplitudes.push_back({is, EvaluateCircuit(&input)});
      }
    }

    if (std::size(initial_states) == 1 and std::size(final_states) == 1)
      return PyCPP::Build(std::get<1>(amplitudes[0]));
    else
      return PyCPP::Build(amplitudes);

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

  return nullptr;
}
