#include "utils.h"

namespace qflex {

std::size_t compute_depth(std::istream &&istream) {
  auto is_number = [](const std::string &token) {
    try {
      std::stol(token);
    } catch (...) {
      return false;
    }
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

    // Remove spaces between parentheses
    line = std::regex_replace(line, std::regex("\\s+(?=[^()]*\\))"), "");

    return line;
  };

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

  std::size_t line_counter{0}, last_cycle_number{0};
  std::string line;

  auto error_msg = [&line, &line_counter](const std::string &msg) {
    std::string err_msg =
        "[" + std::to_string(line_counter + 1) + ": " + line + "] " + msg;
    std::cerr << err_msg << std::endl;
    return err_msg;
  };

  std::size_t depth{0};
  std::unordered_map<std::pair<std::size_t, std::size_t>, std::size_t> layers;

  while (std::getline(istream, line)) {
    if (std::size(line = strip_line(line))) {
      // Check that there are only one '(' and one ')'
      if (std::count(std::begin(line), std::end(line), '(') > 1 or
          std::count(std::begin(line), std::end(line), ')') > 1)
        throw error_msg("Wrong format.");

      // Tokenize the line
      auto tokens = tokenize(line);

      // Enforce first line to be a single number which correspond to the number
      // of qubits
      if (line_counter == 0) {
        if (std::size(tokens) != 1 or std::stol(tokens[0]) <= 0)
          throw error_msg(
              "First line in circuit must be the number of active qubits.");
      } else {
        // Check the correct number of tokens
        if (std::size(tokens) < 3)
          throw error_msg(
              "Gate must be specified as: cycle gate_name[(p1[,p2,...])] q1 "
              "[q2, ...]");

        // Check the first token is actually a number
        if (not is_integer(tokens[0]))
          throw error_msg("First token must be a valid cycle number.");

        std::size_t cycle = std::stol(tokens[0]);

        // Check that cycle number is monotonically increasing
        if (cycle < last_cycle_number)
          throw error_msg("Cycle number can only increase.");

        // Add all the qubits
        if (std::size_t num_qubits = std::size(tokens) - 2; num_qubits == 2) {
          std::size_t q1 = std::stol(tokens[2]);
          std::size_t q2 = std::stol(tokens[3]);
          if (q1 > q2) std::swap(q1, q2);

          if (auto new_depth = ++layers[{q1, q2}]; new_depth > depth)
            depth = new_depth;
          else if (num_qubits > 2)
            throw error_msg("k-qubit gates are not supported.");
        }
      }
    }

    // Increment line counter
    ++line_counter;
  }

  return depth;
}

}  // namespace qflex
