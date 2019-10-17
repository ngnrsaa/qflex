#include "circuit.h"

namespace qflex {

std::ostream &QflexGate::operator<<(std::ostream &out) const { 
  out << "gate_name: " << this->name << std::endl;
  out << "qubits: ";
  for(const auto &q : this->qubits) out << q << " ";
  out << std::endl;
  if(std::size(this->params)) {
    out << "params: ";
    for(const auto &p : this->params) out << p << " ";
    out << std::endl;
  }
  return out;
}
std::ostream &operator<<(std::ostream &out, const QflexGate &gate) { return gate.operator<<(out); }

void QflexCircuit::load(std::istream& istream) {

  auto is_number = [](const std::string &token) {
    try { std::stol(token); } catch(...) { return false; }
    return true;
  };

  auto is_integer = [&is_number](const std::string &token) {
    return is_number(token) and std::stol(token) == std::stod(token);
  };

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

    // Remove any space before '('
    line = std::regex_replace(line, std::regex("[\\s]+[(]"), "(");

    //line = std::regex_replace(line, std::regex("\((?=\))"), "[$&]");
    line = std::regex_replace(line, std::regex("\\s+(?=[^()]*\\))"), "");

    return line;
  };

  auto tokenize = [](const std::string &line, const std::string &regex_expr = "[^\\s]+") {
    std::vector<std::string> tokens;
    auto word_regex = std::regex(regex_expr);
    for(auto w = std::sregex_iterator(std::begin(line), std::end(line), word_regex); w != std::sregex_iterator(); ++w)
      tokens.push_back(w->str());
    return tokens;
  };

  std::size_t line_counter{0}, last_cycle_number{0};
  std::string line;

  auto error_msg = [&line, &line_counter](const std::string &msg) {
    return ERROR_MSG("[" + std::to_string(line_counter+1) + ": " + line + "] " + msg);
  };

  while(std::getline(istream, line)) {
    
    if(std::size(line = strip_line(line))) {

      // Tokenize the line
      auto tokens = tokenize(line);

      // Enforce first line to be a single number which correspond to the number of qubits
      if(line_counter == 0) {
        if(std::size(tokens) != 1 or std::stol(tokens[0]) <= 0)
          throw error_msg("First line in circuit must be the number of active qubits.");
        this->num_active_qubits = std::stol(tokens[0]);
      } else {
        // Check the correct number of tokens
        if(std::size(tokens) < 3)
          throw error_msg("Gate must be specified as: cycle gate_name[(p1[,p2,...])] q1 [q2, ...]");

        // Get gate
        this->gates.emplace_back();
        auto &gate = gates.back();

        // Check the first token is actually a number
        if(not is_integer(tokens[0])) throw error_msg("First token must be a valid cycle number.");
        gate.cycle = std::stol(tokens[0]);

        // Get gate name and params
        if(std::size_t beg = tokens[1].find_first_of('('); beg != std::string::npos) {
          gate.name = tokens[1].substr(0, beg);
          std::string params = tokens[1].substr(beg+1, tokens[1].find_first_of(')')-beg-1);
          for(const auto &p: tokenize(params, "[^,]+")) {
            if(not is_number(p)) throw ERROR_MSG("Params must be valid numbers.");
            gate.params.push_back(std::stod(p));
          }
        } else {
          gate.name = tokens[1];
        }

        // Add all the qubits
        for(std::size_t i = 2; i < std::size(tokens); ++i) {
          if(not is_integer(tokens[i])) throw error_msg("Qubit must be a valid number.");
          gate.qubits.push_back(std::stol(tokens[i]));  
        }
      }
    }

    // Increment line counter
    ++line_counter;

  }

}

void QflexCircuit::load(const std::string& filename) {
  if (auto in = std::ifstream(filename); in.good())
    this->load(in);
  else
    throw ERROR_MSG("Cannot open circuit file " + filename + ".");
};

}
