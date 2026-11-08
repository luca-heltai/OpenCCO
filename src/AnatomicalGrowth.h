// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once

#include "AnatomicalTree.h"
#include "CoronaryArteryTree.h"
#include <set>

namespace OpenCCO {

struct PhysicalBranch {
  std::size_t a, b;
  double radiusA, radiusB;
  bool generated;
  std::size_t sourceEdge;
};

// Source records are copied verbatim and only exposed through const accessors.
// New CCO branches are appended; no imported coordinate or radius is optimized.
class GrownAnatomy {
public:
  explicit GrownAnatomy(const AnatomicalTree &source) : source_(source) {
    validateAnatomicalTree(source_);
    for (const auto &node : source.nodes) nodes_.push_back(node.position);
    for (std::size_t i = 0; i < source.edges.size(); ++i) {
      const auto &edge = source.edges[i];
      branches_.push_back({edge[0], edge[1], source.nodes[edge[0]].radius,
                           source.nodes[edge[1]].radius, false, i});
    }
  }
  const AnatomicalTree &source() const { return source_; }
  const std::vector<std::array<double, 3>> &nodes() const { return nodes_; }
  const std::vector<PhysicalBranch> &branches() const { return branches_; }
  const std::vector<std::size_t> &tissueTerminals() const { return terminals_; }

  template <typename Tree>
  void append(std::size_t attachment, const Tree &tree,
              const typename Tree::TPointD &origin, double pitch) {
    if (attachment >= source_.nodes.size() || !source_.nodes[attachment].grow ||
        attachments_.count(attachment) || !std::isfinite(pitch) || pitch <= 0.)
      throw std::runtime_error("Invalid or repeated anatomical attachment");
    const auto count = tree.myVectSegments.size();
    if (count < 2 || tree.myVectParent.size() != count || tree.myVectParent[0] != 0)
      throw std::runtime_error("Malformed generated subtree");
    const auto &originalRoot = source_.nodes[attachment].position;
    for (unsigned int j = 0; j < 3; ++j) {
      const double actual = origin[j] + pitch*tree.myVectSegments[0].myCoordinate[j];
      if (!std::isfinite(actual) || std::abs(actual-originalRoot[j]) >
          64.*std::numeric_limits<double>::epsilon()*std::max(1., std::abs(originalRoot[j])))
        throw std::runtime_error("Generated root moved away from anatomical attachment");
    }
    std::vector<std::vector<std::size_t>> children(count);
    for (std::size_t i = 1; i < count; ++i) {
      const auto parent = tree.myVectParent[i];
      if (parent >= count || parent == i) throw std::runtime_error("Invalid generated parent");
      children[parent].push_back(i);
    }
    std::vector<bool> seen(count, false);
    std::vector<std::size_t> pending{0};
    seen[0] = true;
    for (std::size_t i = 0; i < pending.size(); ++i)
      for (const auto child : children[pending[i]]) {
        if (seen[child]) throw std::runtime_error("Generated subtree contains a cycle");
        seen[child] = true;
        pending.push_back(child);
      }
    if (pending.size() != count || children[0].size() != 1)
      throw std::runtime_error("Generated subtree is not rooted and connected");
    std::set<std::array<double, 3>> distinctTips;
    std::vector<std::array<double, 3>> generated;
    std::size_t leaves = 0;
    for (std::size_t i = 1; i < count; ++i) {
      const auto &segment = tree.myVectSegments[i];
      if (segment.myIndex != i || !std::isfinite(segment.myRadius*pitch) || segment.myRadius*pitch <= 0.)
        throw std::runtime_error("Invalid generated index/radius");
      std::array<double, 3> point;
      for (unsigned int j = 0; j < 3; ++j) {
        point[j] = origin[j] + pitch*segment.myCoordinate[j];
        if (!std::isfinite(point[j])) throw std::runtime_error("Nonfinite generated point");
      }
      if ((segment.myCoordinate - tree.myVectSegments[tree.myVectParent[i]].myCoordinate).norm() <= 0.)
        throw std::runtime_error("Zero-length generated branch");
      if (children[i].empty()) {
        ++leaves;
        if (!distinctTips.insert(point).second) throw std::runtime_error("Coincident generated terminals");
      } else if (children[i].size() != 2) throw std::runtime_error("Nonbinary generated branch");
      generated.push_back(point);
    }
    if (leaves != tree.myKTerm) throw std::runtime_error("Generated terminal count mismatch");
    for (std::size_t i = 1; i < count; ++i) {
      const auto parent = tree.myVectParent[i];
      const auto &a = parent == 0 ? originalRoot : generated[parent-1];
      const auto &b = generated[i-1];
      const double length = std::hypot(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
      if (!std::isfinite(length) || length <= 0.)
        throw std::runtime_error("Generated branch lost physical length in world coordinates");
    }
    // Commit only after validation: failed append cannot mutate the original or
    // leave an incomplete prefix in the physical export.
    const auto first = nodes_.size();
    nodes_.insert(nodes_.end(), generated.begin(), generated.end());
    const auto mapped = [&](std::size_t i) { return i == 0 ? attachment : first+i-1; };
    for (std::size_t i = 1; i < count; ++i) {
      const double radius = pitch*tree.myVectSegments[i].myRadius;
      branches_.push_back({mapped(tree.myVectParent[i]), mapped(i), radius, radius, true,
                           std::numeric_limits<std::size_t>::max()});
      if (children[i].empty()) terminals_.push_back(mapped(i));
    }
    attachments_.insert(attachment);
  }

private:
  AnatomicalTree source_;
  std::vector<std::array<double, 3>> nodes_;
  std::vector<PhysicalBranch> branches_;
  std::vector<std::size_t> terminals_;
  std::set<std::size_t> attachments_;
};

struct SubtreeParameters {
  unsigned int terminals;
  double volumeVoxels, flowVoxelCubedPerSecond, pressureDropPa, viscosityPaSecond;
};

template <typename Domain>
CoronaryArteryTree<Domain, 3>
seedAnatomicalSubtree(Domain &domain,
                      const DGtal::PointVector<3, double> &root,
                      const DGtal::PointVector<3, double> &terminal,
                      const SubtreeParameters &parameters) {
  using Tree = CoronaryArteryTree<Domain, 3>;
  if (parameters.terminals < 1 || !std::isfinite(parameters.volumeVoxels) || !(parameters.volumeVoxels > 0.) ||
      !std::isfinite(parameters.flowVoxelCubedPerSecond) || !(parameters.flowVoxelCubedPerSecond > 0.) ||
      !std::isfinite(parameters.pressureDropPa) || !(parameters.pressureDropPa > 0.) ||
      !std::isfinite(parameters.viscosityPaSecond) || !(parameters.viscosityPaSecond > 0.) ||
      !std::isfinite(root.norm()) || !std::isfinite(terminal.norm()) || root == terminal)
    throw std::runtime_error("Invalid anatomical subtree seed parameters");
  Tree tree(domain);
  tree.my_NTerm = parameters.terminals;
  tree.myKTerm = 1;
  tree.my_aPerf = parameters.volumeVoxels;
  tree.my_rPerf = std::cbrt(3.*tree.my_aPerf/(4.*M_PI));
  tree.myRsupp = std::cbrt(3.*tree.my_aPerf/(4.*M_PI*tree.my_NTerm));
  tree.my_qPerf = parameters.flowVoxelCubedPerSecond;
  tree.my_qTerm = tree.my_qPerf/tree.my_NTerm;
  tree.my_pDrop = parameters.pressureDropPa;
  tree.my_pTerm = 0.;
  tree.my_pPerf = parameters.pressureDropPa;
  tree.my_mu = parameters.viscosityPaSecond;
  tree.myTreeCenter = terminal;
  tree.myLengthFactor = 1.;
  tree.myVectSegments.resize(2);
  tree.myVectSegments[0].myCoordinate = root;
  tree.myVectSegments[0].myIndex = 0;
  tree.myVectSegments[0].myKTerm = 0;
  tree.myVectSegments[1].myCoordinate = terminal;
  tree.myVectSegments[1].myIndex = 1;
  tree.myVectSegments[1].myFlow = tree.my_qTerm;
  tree.myVectParent = {0, 0};
  tree.myVectChildren = {{0, 0}, {0, 0}};
  tree.myVectTerminals = {1};
  tree.updateResistanceFromRoot();
  tree.updateRootRadius();
  return tree;
}

} // namespace OpenCCO
