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
    --initial-conf=<initial_conf>          Initial configuration.
    --final-conf=<final_conf>              Final configuration.
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

    // Get initial/final configurations
    if (static_cast<bool>(args["--initial-conf"]))
      input.initial_state = args["--initial-conf"].asString();
    else if (static_cast<bool>(args["<initial_conf>"]))
      input.initial_state = args["<initial_conf>"].asString();

    if (static_cast<bool>(args["--final-conf"]))
      input.final_state = args["--final-conf"].asString();
    else if (static_cast<bool>(args["<final_conf>"]))
      input.final_state = args["<final_conf>"].asString();

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

    // Load circuit
    input.circuit.load(std::ifstream(circuit_filename));

    // Load ordering
    input.ordering.load(std::ifstream(ordering_filename));

    // Load grid
    input.grid.load(grid_filename);

    // Evaluating circuit.
    std::vector<std::pair<std::string, std::complex<double>>> amplitudes;
    try {
      amplitudes = qflex::EvaluateCircuit(&input);
    } catch (const std::string& err_msg) {
      throw ERROR_MSG("Failed to call EvaluateCircuit(). Error:\n\t[", err_msg,
                      "]");
    }
    // If no error is caught, amplitudes will be initialized.

    if (qflex::global::verbose > 0 && qflex::global::track_memory_seconds > 0)
      qflex::memory::print_memory_usage();

    // Printing output.
    for (std::size_t c = 0; c < amplitudes.size(); ++c) {
      const auto& state = amplitudes[c].first;
      const auto& amplitude = amplitudes[c].second;
      std::cout << input.initial_state << " --> " << state << ": "
                << std::real(amplitude) << " " << std::imag(amplitude)
                << std::endl;
    }

  } catch (const std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    return 1;

  } catch (const std::string& msg) {
    std::cerr << msg << std::endl;
    return 2;
  }

  return 0;
}
