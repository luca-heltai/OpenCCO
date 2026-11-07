// SPDX-License-Identifier: GPL-3.0-or-later
#include "AnatomicalTree.h"

#include <filesystem>
#include <fstream>
#include <iostream>

int main(int argc, char **argv) {
  if (argc != 2 && argc != 3) {
    std::cerr << "Usage: inspectAnatomy INPUT [ROUNDTRIP_OUTPUT]\n";
    return 2;
  }
  try {
    std::ifstream input(argv[1]);
    if (!input) throw std::runtime_error("Cannot open anatomical input");
    const auto tree = OpenCCO::readAnatomicalTree(input);
    if (argc == 3) {
      if (std::filesystem::exists(argv[2]))
        throw std::runtime_error("Refusing to overwrite roundtrip output");
      std::ofstream output(argv[2]);
      OpenCCO::writeAnatomicalTree(output, tree);
    }
    std::size_t grow = 0;
    for (const auto &node : tree.nodes) grow += node.grow;
    std::cout << "{\"role\":\"" << tree.role << "\",\"units\":\"metres\",\"nodes\":"
              << tree.nodes.size() << ",\"edges\":" << tree.edges.size()
              << ",\"root\":" << tree.root() << ",\"growth_endpoints\":" << grow
              << ",\"connected_tree\":true,\"growth_performed\":false}\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
