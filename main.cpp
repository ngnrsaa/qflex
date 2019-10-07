#include <omp.h>

#include "evaluate_circuit.h"
#include "docopt.h"

static const char VERSION[] = "qFlex v1.0";
static const char USAGE[] =
R"(Flexible Quantum Circuit Simulator (qFlex) implements an efficient
tensor network, CPU-based simulator of large quantum circuits.

  Usage:
    qflex <grid_x> <grid_y> <depth> <fidelity> <circuit_filename> <ordering_filename> <grid_filename> [<initial_conf> <final_conf>]
    qflex -x <grid_x> -y <grid_y> -d <depth> -f <fidelity> -c <circuit_filename> -o <ordering_filename> -g <grid_filename> [--initial-conf <initial_conf> --final-conf <final_conf>]
    qflex (-h | --help)
    qflex --version

  Options:
    -h,--help                              Show this help.
    -x <grid_x>                            Size of x-axis of grid.
    -y <grid_y>                            Size of y-axis of grid.
    -d,--depth=<depth>                     Target circuit depth.
    -f,--fidelity=<fidelity>               Target circuit fidelity.
    -c,--circuit=<circuit_filename>        Circuit filename.
    -o,--ordering=<ordering_filename>      Ordering filename.
    -g,--grid=<grid_filename>              Grid filename.
    --initial-conf=<initial_conf>          Initial configuration.
    --final-conf=<final_conf>              Final configuration.
    --version                              Show version.

)";

// Example:
// $ ./qflex.x 11 12 2 0.005 ./circuits/bristlecone_48_1-40-1_0.txt \
//       ./ordering/bristlecone_48.txt ./grid/bristlecone_48.txt
//
int main(int argc, char** argv) {

  std::map<std::string, docopt::value> args = docopt::docopt(USAGE, {argv+1, argv+argc}, true, VERSION);

  // Reading input
  qflex::QflexInput input;
  input.I = bool(args["-x"])?args["-x"].asLong():args["<grid_x>"].asLong();
  input.J = bool(args["-y"])?args["-y"].asLong():args["<grid_y>"].asLong();
  input.K = bool(args["--depth"])?args["--depth"].asLong():args["<depth>"].asLong();
  input.fidelity = std::stof(bool(args["--fidelity"])?args["--fidelity"].asString():args["<fidelity>"].asString());

  // Get initial/final configurations
  if(bool(args["--initial-conf"])) input.initial_state = args["--initial-conf"].asString();
  else if(bool(args["<initial_conf>"])) input.initial_state = args["<initial_conf>"].asString();

  if(bool(args["--final-conf"])) input.final_state_A = args["--final-conf"].asString();
  else if(bool(args["<final_conf>"])) input.final_state_A = args["<final_conf>"].asString();

  // Getting filenames
  std::string circuit_filename = bool(args["--circuit"])?args["--circuit"].asString():args["<circuit_filename>"].asString();
  std::string ordering_filename = bool(args["--ordering"])?args["--ordering"].asString():args["<ordering_filename>"].asString();
  std::string grid_filename = bool(args["--grid"])?args["--grid"].asString():args["<grid_filename>"].asString();

  // Creating streams for input files.
  auto circuit_data = std::ifstream(circuit_filename);
  if (!circuit_data.good()) {
    std::cout << "Cannot open circuit data file: " << circuit_filename
              << std::endl;
    assert(circuit_data.good());
  }
  input.circuit_data = &circuit_data;
  auto ordering_data = std::ifstream(ordering_filename);
  if (!ordering_data.good()) {
    std::cout << "Cannot open ordering data file: " << ordering_filename
              << std::endl;
    assert(ordering_data.good());
  }
  input.ordering_data = &ordering_data;
  auto grid_data = std::ifstream(grid_filename);
  if (!grid_data.good()) {
    std::cout << "Cannot open grid data file: " << grid_filename << std::endl;
    assert(grid_data.good());
  }
  input.grid_data = &grid_data;

  // Evaluating circuit.
  std::vector<std::pair<std::string, std::complex<double>>> amplitudes =
      qflex::EvaluateCircuit(&input);

  // Printing output.
  for (int c = 0; c < amplitudes.size(); ++c) {
    const auto& state = amplitudes[c].first;
    const auto& amplitude = amplitudes[c].second;
    std::cout << input.initial_state << " --> " << state << ": "
              << std::real(amplitude) << " " << std::imag(amplitude)
              << std::endl;
  }
  return 0;
}
