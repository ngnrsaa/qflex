#include "pybind_main.h"

inline PyObject* simulate(PyObject *, PyObject *arg) {

  auto GetStream = [](const auto &data, const std::string &caller) {
    std::stringstream out;
    for(const auto &line: data)
      if(PyCPP::Is<PyCPP::String>(line))
        out << line << std::endl;
      else
        throw PyCPP::SetError(caller + " must be a list of strings.");
    return out;
  };

  auto LoadData = [&GetStream](const auto &options, const std::string &argument, auto &data) {
    if(auto w = options.find(PyCPP::String(argument)); w != std::end(options)) {
      const auto &q = std::get<1>(*w);
      if(PyCPP::Is<PyCPP::Sequence>(q))
        data.load(GetStream(PyCPP::As<PyCPP::Sequence>(q), argument));
      else
        throw PyCPP::SetError("'" + argument + "' must be a list of strings.");
    } else if(auto w = options.find(PyCPP::String(argument + "-filename")); w != std::end(options)) {
      const auto &q = std::get<1>(*w);
      if(PyCPP::Is<PyCPP::String>(q))
        data.load(PyCPP::As<PyCPP::String>(q));
      else
        throw PyCPP::SetError("'" + argument + "-filename' must be a string.");
    }
  };

  auto LoadStates = [](const auto &options, const std::string &argument, const std::size_t num_active_qubits = 0) {
    std::vector<std::string> states;
    if(auto w = options.find(PyCPP::String(argument + "-states")); w != std::end(options)) {
      const auto &q = std::get<1>(*w);
      if(PyCPP::Is<PyCPP::String>(q)) { 
        states.push_back(PyCPP::As<PyCPP::String>(q));
      } else if(PyCPP::Is<PyCPP::Sequence>(q)) {
        for(const auto &x: PyCPP::As<PyCPP::Sequence>(q))
          if(PyCPP::Is<PyCPP::String>(x)) states.push_back(PyCPP::As<PyCPP::String>(x));
          else throw PyCPP::SetError("states must be strings.");
      } else throw PyCPP::SetError(argument + "-states must be a list of strings.");
    } else {
      if(argument == "initial") states.push_back(std::string(num_active_qubits, '0'));
      else states.push_back("");
    }
    return states;
  };

  try {

    // The passed argument must be a dictionary
    const auto v_options = PyCPP::Parse(arg);
    if(not PyCPP::Is<PyCPP::Map>(v_options[0])) {
      throw PyCPP::SetError("Simulation takes only one argument and it must be a dictionary");
    }
    const PyCPP::Map &options = PyCPP::As<PyCPP::Map>(v_options[0]);

    // Get qFlex input
    qflex::QflexInput input;

    // Load grid, circuit and ordering
    LoadData(options, "grid", input.grid);
    LoadData(options, "circuit", input.circuit);
    LoadData(options, "ordering", input.ordering);
    const auto initial_states = LoadStates(options, "initial", input.circuit.num_active_qubits);
    const auto final_states = LoadStates(options, "final");

    // Define container for amplitudes
    std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::complex<double>>>>> amplitudes;

    auto run_simulation = [&](const std::string &initial_state = "", const std::string &final_state = "") {
      input.initial_state = initial_state;
      input.final_state = final_state;
      return EvaluateCircuit(&input);
    };

    for(const auto &is: initial_states) {
      for(const auto &fs: final_states) {
        amplitudes.push_back({is, run_simulation(is, fs)});
      }
    }

    return PyCPP::Build(amplitudes);

  } catch(...) {

    return nullptr;

  }

}
