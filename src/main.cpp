#include "memory.h"
#include "circuit.h"
#include "docopt.h"
#include "evaluate_circuit.h"
#include "global.h"
#include "grid.h"
#include "input.h"
#include "ordering.h"

static const char VERSION[] = "qFlex v0.1";
static const char USAGE[] =
    R"(Flexible Quantum Circuit Simulator (qFlex) implements an efficient
tensor network, CPU-based simulator of large quantum circuits.

  Usage:
    qflex <circuit_filename> <ordering_filename> <grid_filename> [<initial_conf> <final_conf> --verbosity <verbosity_level> --memory <memory_limit>]
    qflex -c <circuit_filename> -o <ordering_filename> -g <grid_filename> [--initial-conf <initial_conf> --final-conf <final_conf> --verbosity <verbosity_level> --memory <memory_limit>]
    qflex (-h | --help)
    qflex --version

  Options:
    -h,--help                              Show this help.
    -c,--circuit=<circuit_filename>        Circuit filename.
    -o,--ordering=<ordering_filename>      Ordering filename.
    -g,--grid=<grid_filename>              Grid filename.
    -v,--verbosity=<verbosity_level>       Verbosity level.
    -m,--memory=<memory_limit>             Memory limit (default 1GB).
    --initial-conf=<initial_conf>          Initial configuration.
    --final-conf=<final_conf>              Final configuration.
    --version                              Show version.

)";

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
        docopt::docopt(USAGE, {argv + 1, argv + argc}, true, VERSION);

    // Reading input
    qflex::QflexInput input;

    // Update global qflex::global::verbose
    if (static_cast<bool>(args["--verbosity"]))
      qflex::global::verbose = args["--verbosity"].asLong();
    else if (static_cast<bool>(args["<verbosity_level>"]))
      qflex::global::verbose = args["<verbosity_level>"].asLong();
    else
      qflex::global::verbose = 0;

    // Update global qflex::global::memory_limit
    if (static_cast<bool>(args["--memory"]))
      qflex::global::memory_limit = args["--memory"].asLong();
    else if (static_cast<bool>(args["<memory_limit>"]))
      qflex::global::memory_limit = args["<memory_limit>"].asLong();
    else
      // Default limit of one gigabyte.
      qflex::global::memory_limit = 1L << 30;

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

    // set alarms to get memory usage in real time
    if(qflex::global::verbose > 0) {
      signal(SIGALRM, qflex::memory::print_peak_memory_usage);
      ualarm(1e5,1e5);
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
    
    if(qflex::global::verbose > 0)
      std::cerr << "Peak memory usage: " << qflex::memory::get_peak_memory_usage() << std::endl;

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
