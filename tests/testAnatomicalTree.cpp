#include "AnatomicalTree.h"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}

const std::string valid =
    "OpenCCOAnatomy 1\narterial metres\n4 3\n"
    "0 0 0 0.001 0 0\n"
    "0 1 0 0.002 1 0\n"
    "-1 -1 0 0.0005 0 1\n"
    "1 -1 0 0.0006 0 0\n"
    "1 0\n0 2\n0 3\n";

OpenCCO::AnatomicalTree read(const std::string &text) {
  std::istringstream input(text);
  return OpenCCO::readAnatomicalTree(input);
}

void rejects(const std::string &text) {
  bool rejected = false;
  try { read(text); } catch (const std::runtime_error &) { rejected = true; }
  require(rejected, "Malformed anatomy was accepted");
}

std::string replace(std::string text, const std::string &from,
                    const std::string &to) {
  const auto pos = text.find(from);
  require(pos != std::string::npos, "Invalid regression fixture");
  return text.replace(pos, from.size(), to);
}
}

int main() {
  try {
    const auto tree = read(valid);
    require(tree.root() == 1, "Root identity changed");
    require(tree.nodes.size() == 4 && tree.edges.size() == 3, "Wrong counts");
    require(tree.nodes[2].grow && !tree.nodes[3].grow, "Growth flags changed");
    require(tree.edges[0][0] == 1 && tree.edges[0][1] == 0, "Orientation changed");
    require(tree.nodes[3].radius == 0.0006, "Radius changed");
    std::ostringstream output;
    OpenCCO::writeAnatomicalTree(output, tree);
    const auto copy = read(output.str());
    require(copy.role == tree.role && copy.edges == tree.edges, "Connectivity roundtrip");
    for (std::size_t i = 0; i < tree.nodes.size(); ++i) {
      require(copy.nodes[i].position == tree.nodes[i].position &&
              copy.nodes[i].radius == tree.nodes[i].radius &&
              copy.nodes[i].root == tree.nodes[i].root &&
              copy.nodes[i].grow == tree.nodes[i].grow, "Node roundtrip");
    }
    rejects(replace(valid, "arterial metres", "arterial mm"));
    rejects(replace(valid, "OpenCCOAnatomy 1", "OpenCCOAnatomy 2"));
    rejects(replace(valid, "0.002 1 0", "0.002 0 0"));
    rejects(replace(valid, "0.001 0 0", "0.001 1 0"));
    rejects(replace(valid, "0.002 1 0", "0.002 1 1"));
    rejects(replace(valid, "0.001 0 0", "0.001 0 1"));
    rejects(replace(valid, "0.0005 0 1", "-1 0 1"));
    rejects(replace(valid, "0.0005 0 1", "nan 0 1"));
    rejects(replace(valid, "0.0005 0 1", "inf 0 1"));
    rejects(replace(valid, "0.0005 0 1", "0.0005 0 2"));
    rejects(replace(valid, "0 3\n", "0 4\n"));
    rejects(replace(valid, "0 3\n", "0 -1\n"));
    rejects(replace(valid, "0 3\n", "0 2\n"));
    rejects(replace(valid, "0 3\n", "2 0\n"));
    rejects(replace(valid, "0 3\n", "3 3\n"));
    rejects(replace(valid, "-1 -1 0", "0 0 0"));
    rejects(replace(valid, "1 0\n0 2\n0 3", "0 2\n2 3\n3 0"));
    rejects(valid + "unexpected\n");
    rejects(valid.substr(0, valid.size() - 3));
    rejects("OpenCCOAnatomy 1\narterial metres\n-1 0\n");
    std::cout << "Anatomical tree validation and exact roundtrip passed\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
