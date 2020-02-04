#include <stdlib.h>

#include "circuit.h"
#include "docopt.h"
#include "evaluate_circuit.h"
#include "global.h"
#include "grid.h"
#include "input.h"
#include "memory.h"
#include "ordering.h"
#include "utils.h"

static const char VERSION[] = "qFlex v0.1";
const std::string USAGE = qflex::utils::concat(
    R"(Flexible Quantum Circuit Simulator (qFlex) implements an efficient
tensor network, CPU-based simulator of large quantum circuits.

  Usage:
    qflex <circuit_filename> <ordering_filename> <grid_filename> [<initial_conf> <final_conf> <verbosity_level> <memory_limit> <track_memory_seconds>]
    qflex -c <circuit_filename> -o <ordering_filename> -g <grid_filename> [--initial-conf=<initial_conf> --final-conf=<final_conf> --verbosity=<verbosity_level> --memory=<memory_limit> --track-memory=<seconds>]
    qflex (-h | --help)
    qflex --version

  Options:
    -h,--help                              Show this help.
    -c,--circuit=<circuit_filename>        Circuit filename.
    -o,--ordering=<ordering_filename>      Ordering filename.
    -g,--grid=<grid_filename>              Grid filename.
    -v,--verbosity=<verbosity_level>       Verbosity level [default: )",
    qflex::global::verbose, R"(].
    -m,--memory=<memory_limit>             Memory limit [default: )",
    qflex::utils::readable_memory_string(qflex::global::memory_limit), R"(].
    -t,--track-memory=<seconds>            If <verbosity_level> > 0, track memory usage [default: )",
    qflex::global::track_memory_seconds, R"(].
    --initial-conf=<initial_conf>          Initial configuration [default: 00...00].
    --final-conf=<final_conf>              Final configuration [default: 00...00].
    --version                              Show version.
)");

/*
 * Example:
 * $ src/qflex.x config/circuits/bristlecone_48_1-24-1_0.txt \
 *               config/ordering/bristlecone_48.txt \
 *               config/grid/bristlecone_48.txt
 *
 */
int main(int argc, char** argv) {
  try {
    std::map<std::string, docopt::value> args =
        docopt::docopt(USAGE.c_str(), {argv + 1, argv + argc}, true, VERSION);

    // Reading input
    qflex::QflexInput input;

    // Update global qflex::global::verbose
    if (static_cast<bool>(args["--verbosity"]))
      qflex::global::verbose = args["--verbosity"].asLong();
    else if (static_cast<bool>(args["<verbosity_level>"]))
      qflex::global::verbose = args["<verbosity_level>"].asLong();

    // Update global qflex::global::memory_limit
    if (static_cast<bool>(args["--memory"]))
      qflex::global::memory_limit = qflex::utils::from_readable_memory_string(
          args["--memory"].asString());
    else if (static_cast<bool>(args["<memory_limit>"]))
      qflex::global::memory_limit = qflex::utils::from_readable_memory_string(
          args["<memory_limit>"].asString());

    if (static_cast<bool>(args["--track-memory"]))
      qflex::global::track_memory_seconds = args["--track-memory"].asLong();
    else if (static_cast<bool>(args["<track_memory_seconds>"]))
      qflex::global::track_memory_seconds =
          args["<track_memory_seconds>"].asLong();

    // Getting filenames
    std::string circuit_filename = static_cast<bool>(args["--circuit"])
                                       ? args["--circuit"].asString()
                                       : args["<circuit_filename>"].asString();
    std::string ordering_filename =
        static_cast<bool>(args["--ordering"])
            ? args["--ordering"].asString()
            : args["<ordering_filename>"].asString();
    std::string grid_filename = static_cast<bool>(args["--grid"])
                                    ? args["--grid"].asString()
                                    : args["<grid_filename>"].asString();

    // Load circuit
    input.circuit.load(std::ifstream(circuit_filename));

    // Load ordering
    input.ordering.load(std::ifstream(ordering_filename));

    // Load grid
    input.grid.load(grid_filename);

    // Get initial/final configurations
    for (const auto& arg : {"--initial-conf", "<initial_conf>"})
      if (static_cast<bool>(args[arg])) {
        std::stringstream ss(args[arg].asString());
        std::string val;
        while (std::getline(ss, val, ',')) {
          if (val == "00...00")
            input.initial_states.push_back(
                std::string(input.circuit.num_active_qubits, '0'));
          else if (val.find_first_not_of("01") != std::string::npos)
            throw ERROR_MSG(
                "Initial configurations can only have 0 or 1 characters.");
          else if (std::size(val) != input.circuit.num_active_qubits)
            throw ERROR_MSG(
                "Initial configurations must have the number of active "
                "qubits.");
          else
            input.initial_states.push_back(val);
        }
      }

    for (const auto& arg : {"--final-conf", "<final_conf>"})
      if (static_cast<bool>(args[arg])) {
        std::stringstream ss(args[arg].asString());
        std::string val;
        while (std::getline(ss, val, ',')) {
          if (val == "00...00")
            input.final_states.push_back(
                std::string(input.circuit.num_active_qubits, '0'));
          else if (val.find_first_not_of("01") != std::string::npos)
            throw ERROR_MSG(
                "Final configurations can only have 0 or 1 characters.");
          else if (std::size(val) != input.circuit.num_active_qubits)
            throw ERROR_MSG(
                "Final configurations must have the number of active qubits.");
          else
            input.final_states.push_back(val);
        }
      }

    // Delete duplicate initial/final configuration
    input.initial_states.erase(std::unique(std::begin(input.initial_states),
                                           std::end(input.initial_states)),
                               std::end(input.initial_states));
    input.final_states.erase(std::unique(std::begin(input.final_states),
                                         std::end(input.final_states)),
                             std::end(input.final_states));

    // Print OMP_NUM_THREADS and MKL_NUM_THREADS
    if (qflex::global::verbose > 0)
      for (const char* var : {"OMP_NUM_THREADS", "MKL_NUM_THREADS"})
        if (const char* value = getenv(var); value != nullptr)
          std::cerr << WARN_MSG(var, " = ", value) << std::endl;

    // Print info on maximum memory
    if (qflex::global::verbose > 0)
      std::cerr << WARN_MSG("Maximum allowed memory: ",
                            qflex::utils::readable_memory_string(
                                qflex::global::memory_limit))
                << std::endl;

    // set alarms to get memory usage in real time
    if (qflex::global::verbose > 0 && qflex::global::track_memory_seconds > 0) {
      signal(SIGALRM, qflex::memory::print_memory_usage);
      alarm(qflex::global::track_memory_seconds);
    }

    // Evaluating circuit.
    std::vector<std::tuple<std::string, std::string, std::complex<double>>>
        amplitudes;
    try {
      amplitudes = qflex::EvaluateCircuit(input);
    } catch (const std::string& err_msg) {
      throw ERROR_MSG("Failed to call EvaluateCircuit(). Error:\n\t[", err_msg,
                      "]");
    }
    // If no error is caught, amplitudes will be initialized.

    if (qflex::global::verbose > 0 && qflex::global::track_memory_seconds > 0)
      qflex::memory::print_memory_usage();

    // Printing output.
    for (const auto& [initial_state, final_state, amplitude] : amplitudes)
      std::cout << initial_state << " --> " << final_state << ": " << amplitude
                << std::endl;

  } catch (const std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    return 1;

  } catch (const std::string& msg) {
    std::cerr << msg << std::endl;
    return 2;
  }

  return 0;
}
