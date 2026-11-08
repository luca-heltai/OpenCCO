// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <istream>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace OpenCCO {

// Imported physical data, deliberately separate from mutable CCO solver state.
// Node and edge order and supplied edge orientation are retained verbatim.
struct AnatomicalTree {
  struct Node {
    std::array<double, 3> position;
    double radius;
    bool root;
    bool grow;
  };
  std::string role;
  std::vector<Node> nodes;
  std::vector<std::array<std::size_t, 2>> edges;

  std::size_t root() const {
    for (std::size_t i = 0; i < nodes.size(); ++i)
      if (nodes[i].root) return i;
    throw std::runtime_error("Anatomy has no root");
  }
};

inline void validateAnatomicalTree(const AnatomicalTree &tree) {
  if (tree.role != "arterial" && tree.role != "venous" && tree.role != "portal")
    throw std::runtime_error("Unknown anatomical vascular role");
  if (tree.nodes.size() < 2 || tree.edges.size() != tree.nodes.size() - 1)
    throw std::runtime_error("Anatomy must have N >= 2 nodes and N-1 edges");
  std::size_t roots = 0;
  for (const auto &node : tree.nodes) {
    for (const double x : node.position)
      if (!std::isfinite(x)) throw std::runtime_error("Nonfinite anatomical position");
    if (!std::isfinite(node.radius) || node.radius <= 0.)
      throw std::runtime_error("Anatomical radius must be finite and positive");
    roots += node.root;
  }
  if (roots != 1) throw std::runtime_error("Anatomy must identify exactly one root");
  std::vector<std::vector<std::size_t>> adjacency(tree.nodes.size());
  for (const auto &edge : tree.edges) {
    const auto a = edge[0], b = edge[1];
    if (a >= tree.nodes.size() || b >= tree.nodes.size() || a == b)
      throw std::runtime_error("Invalid anatomical edge endpoint");
    const auto &x = tree.nodes[a].position, &y = tree.nodes[b].position;
    const double length = std::hypot(x[0] - y[0], x[1] - y[1], x[2] - y[2]);
    if (!std::isfinite(length) || length <= 0.)
      throw std::runtime_error("Anatomical edge must have finite positive length");
    adjacency[a].push_back(b);
    adjacency[b].push_back(a);
  }
  for (std::size_t i = 0; i < tree.nodes.size(); ++i) {
    if (tree.nodes[i].root && adjacency[i].size() != 1)
      throw std::runtime_error("Anatomical root must be a degree-one node");
    if (tree.nodes[i].grow && (tree.nodes[i].root || adjacency[i].size() != 1))
      throw std::runtime_error("Grow must identify nonroot degree-one endpoints");
  }
  // Connectedness plus E=N-1 excludes cycles, parallel edges and isolated nodes.
  std::vector<bool> seen(tree.nodes.size(), false);
  std::vector<std::size_t> pending{tree.root()};
  seen[pending[0]] = true;
  for (std::size_t i = 0; i < pending.size(); ++i)
    for (const auto other : adjacency[pending[i]])
      if (!seen[other]) {
        seen[other] = true;
        pending.push_back(other);
      }
  if (pending.size() != tree.nodes.size())
    throw std::runtime_error("Anatomy is disconnected or contains a cycle");
}

inline AnatomicalTree readAnatomicalTree(std::istream &input) {
  AnatomicalTree tree;
  std::string magic, units;
  unsigned int version = 0;
  std::size_t nodes = 0, edges = 0;
  if (!(input >> magic >> version >> tree.role >> units >> nodes >> edges) ||
      magic != "OpenCCOAnatomy" || version != 1 || units != "metres")
    throw std::runtime_error("Invalid OpenCCOAnatomy v1 header (units must be metres)");
  // Bound allocations before parsing untrusted counts, including negative tokens.
  if (nodes < 2 || nodes > 10000000 || edges != nodes - 1)
    throw std::runtime_error("Invalid or oversized anatomical counts");
  tree.nodes.reserve(nodes);
  for (std::size_t i = 0; i < nodes; ++i) {
    AnatomicalTree::Node node;
    int root = -1, grow = -1;
    if (!(input >> node.position[0] >> node.position[1] >> node.position[2] >>
          node.radius >> root >> grow) ||
        (root != 0 && root != 1) || (grow != 0 && grow != 1))
      throw std::runtime_error("Invalid anatomical node record");
    node.root = root == 1;
    node.grow = grow == 1;
    tree.nodes.push_back(node);
  }
  tree.edges.resize(edges);
  for (auto &edge : tree.edges)
    if (!(input >> edge[0] >> edge[1]))
      throw std::runtime_error("Invalid anatomical edge record");
  input >> std::ws;
  if (!input.eof()) throw std::runtime_error("Unexpected trailing anatomical data");
  validateAnatomicalTree(tree);
  return tree;
}

inline void writeAnatomicalTree(std::ostream &output, const AnatomicalTree &tree) {
  validateAnatomicalTree(tree);
  const auto precision = output.precision();
  const auto flags = output.flags();
  output << std::defaultfloat << std::setprecision(std::numeric_limits<double>::max_digits10);
  output << "OpenCCOAnatomy 1\n" << tree.role << " metres\n"
         << tree.nodes.size() << ' ' << tree.edges.size() << '\n';
  for (const auto &node : tree.nodes)
    output << node.position[0] << ' ' << node.position[1] << ' ' << node.position[2]
           << ' ' << node.radius << ' ' << int(node.root) << ' ' << int(node.grow) << '\n';
  for (const auto &edge : tree.edges) output << edge[0] << ' ' << edge[1] << '\n';
  output.precision(precision);
  output.flags(flags);
  if (!output) throw std::runtime_error("Cannot write anatomical tree");
}

} // namespace OpenCCO
