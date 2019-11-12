#include "circuit.h"

namespace qflex {

std::ostream &QflexGate::operator<<(std::ostream &out) const {
  out << "gate_name: " << this->name << std::endl;
  out << "qubits: ";
  for (const auto &q : this->qubits) out << q << " ";
  out << std::endl;
  if (std::size(this->params)) {
    out << "params: ";
    for (const auto &p : this->params) out << p << " ";
    out << std::endl;
  }
  return out;
}
std::ostream &operator<<(std::ostream &out, const QflexGate &gate) {
  return gate.operator<<(out);
}

void QflexCircuit::clear() {
  this->num_active_qubits = 0;
  this->depth = 0;
  gates.clear();
}

void QflexCircuit::load(std::istream &istream) {
  this->load(std::move(istream));
}
void QflexCircuit::load(std::istream &&istream) {
  // Check if token is a number
  auto is_number = [](const std::string &token) {
    try {
      std::stol(token);
    } catch (...) {
      return false;
    }
    return true;
  };

  // Check if token is an integer
  auto is_integer = [&is_number](const std::string &token) {
    return is_number(token) and std::stol(token) == std::stod(token);
  };

  // Given a line, strip everything that doesn't follow the right format
  auto strip_line = [](std::string line) {
    // Remove everything after '#'
    line = std::regex_replace(line, std::regex("#.*"), "");

    // Remove any special character
    line = std::regex_replace(line, std::regex("[^)(\\s\\ta-zA-Z0-9_.,-]"), "");

    // Convert tabs to spaces
    line = std::regex_replace(line, std::regex("[\\t]"), " ");

    // Remove multiple spaces
    line = std::regex_replace(line, std::regex("[\\s]{2,}"), " ");

    // Remove last space
    line = std::regex_replace(line, std::regex("\\s+$"), "");

    // Remove any space between a non-space char and '('
    line = std::regex_replace(line, std::regex("[\\s]+[(]"), "(");

    // Remove spaces between parentheses
    line = std::regex_replace(line, std::regex("\\s+(?=[^()]*\\))"), "");

    // After stripping, line should follow the format
    // 0 gate(p1,p2,...) q1 q2 ...

    return line;
  };

  // Given a line, tokenize it
  auto tokenize = [](const std::string &line,
                     const std::string &regex_expr = "[^\\s]+") {
    std::vector<std::string> tokens;
    auto word_regex = std::regex(regex_expr);
    for (auto w =
             std::sregex_iterator(std::begin(line), std::end(line), word_regex);
         w != std::sregex_iterator(); ++w)
      tokens.push_back(w->str());
    return tokens;
  };

  // Clear this circuit
  this->clear();

  std::size_t line_counter{0}, last_cycle_number{0};
  std::unordered_set<std::size_t> used_qubits;
  std::string line;

  auto error_msg = [&line, &line_counter](const std::string &msg) {
    return ERROR_MSG("[" + std::to_string(line_counter + 1) + ": " + line +
                     "] " + msg);
  };

  while (std::getline(istream, line)) {
    if (std::size(line = strip_line(line))) {
      // Check number of parentheses
      if (std::size_t n_open =
              std::count(std::begin(line), std::end(line), '('),
          n_close = std::count(std::begin(line), std::end(line), ')');
          n_open != n_close or n_open > 1)
        throw ERROR_MSG("Not-matching parentheses.");

      // Tokenize the line
      auto tokens = tokenize(line);

      // Enforce first line to be a single number which correspond to the number
      // of qubits
      if (line_counter == 0) {
        if (std::size(tokens) != 1 or std::stol(tokens[0]) <= 0)
          throw ERROR_MSG(
              "First line in circuit must be the number of active qubits.");
        this->num_active_qubits = std::stol(tokens[0]);
      } else {
        // Check the correct number of tokens
        if (std::size(tokens) < 3)
          throw ERROR_MSG(
              "Gate must be specified as: cycle gate_name[(p1[,p2,...])] q1 "
              "[q2, ...]");

        // Get gate
        this->gates.emplace_back();
        auto &gate = gates.back();

        // Add raw line to gate
        gate.raw = line;

        // Check the first token is actually a number
        if (not is_integer(tokens[0]))
          throw ERROR_MSG("First token must be a valid cycle number.");
        gate.cycle = std::stol(tokens[0]);

        // Check that cycle number is monotonically increasing
        if (gate.cycle < last_cycle_number)
          throw ERROR_MSG("Cycle number can only increase.");

        // If cycle number change, reset used qubits
        if (gate.cycle != last_cycle_number) {
          used_qubits.clear();
          last_cycle_number = gate.cycle;
        }

        // Get gate name and params
        if (std::size_t beg = tokens[1].find_first_of('(');
            beg != std::string::npos) {
          gate.name = tokens[1].substr(0, beg);
          std::string params =
              tokens[1].substr(beg + 1, tokens[1].find_first_of(')') - beg - 1);
          for (const auto &p : tokenize(params, "[^,]+")) {
            if (not is_number(p))
              throw ERROR_MSG("Params must be valid numbers.");
            gate.params.push_back(std::stod(p));
          }
        } else {
          gate.name = tokens[1];
        }

        // Add all the qubits
        for (std::size_t i = 2; i < std::size(tokens); ++i) {
          if (not is_integer(tokens[i]))
            throw ERROR_MSG("Qubit must be a valid number.");
          gate.qubits.push_back(std::stol(tokens[i]));
        }

        // Check that qubits are not already used in the same cycle
        for (const auto &qubit : gate.qubits)
          if (used_qubits.find(qubit) != std::end(used_qubits))
            throw ERROR_MSG("Qubits can only used one for each cycle");
          else
            used_qubits.insert(qubit);
      }
    }

    // Increment line counter
    ++line_counter;
  }

  // Compute circuit depth
  // TODO: have scratch space be allocated not relying on depth.
  {
    std::unordered_map<std::size_t,
                       std::unordered_map<std::size_t, std::size_t>>
        layers;
    for (const auto &gate : this->gates) {
      if (std::size(gate.qubits) > 2)
        throw ERROR_MSG(
            "Depth calculation does not handle k-qubit gates with k > 2.");

      if (std::size(gate.qubits) == 2) {
        auto q1 = gate.qubits[0];
        auto q2 = gate.qubits[1];
        if (q2 > q1) std::swap(q1, q2);

        size_t pair_depth_increment = 1;
        if (gate.name == "fsim") pair_depth_increment = 2;
        layers[q1][q2] += pair_depth_increment;
        size_t new_depth = layers[q1][q2];
        if (new_depth > depth) this->depth = new_depth;
      }
    }
  }
}

void QflexCircuit::load(const std::string &filename) {
  if (auto in = std::ifstream(filename); in.good())
    this->load(in);
  else
    throw ERROR_MSG("Cannot open circuit file " + filename + ".");
};

}  // namespace qflex
