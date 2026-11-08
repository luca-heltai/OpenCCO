// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "CoronaryArteryTree.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace OpenCCO {
struct GeneratedRelaxation {
  bool valid=false;
  unsigned int sweeps=0;
  double maximumDisplacement=0.;
};

/** Repair a disposable generated candidate, never its anatomical root0.
 * Coordinates and radii must share units (length factor1). A failed result is
 * not committed by the caller. This is only a length repair: the caller must
 * still validate its complete mask, wall and collision predicates afterwards.
 * Motion is bounded by one initial inlet diameter; no physical input changes. */
template<class Tree>
GeneratedRelaxation relaxGeneratedLengths(Tree &tree,unsigned int maximumSweeps) {
  GeneratedRelaxation result;
  const auto count=tree.myVectSegments.size();
  if (count<2 || tree.myVectParent.size()!=count || tree.myVectChildren.size()!=count ||
      tree.myVectParent[0]!=0 || tree.myVectParent[1]!=0 || tree.myLengthFactor!=1.)
    throw std::runtime_error("Generated relaxation requires a valid physical-coordinate tree");
  std::vector<unsigned int> order{1};
  std::vector<bool> seen(count,false);seen[0]=seen[1]=true;
  for (std::size_t j=0; j<order.size(); ++j) {
    const auto i=order[j];
    const auto children=tree.myVectChildren[i];
    if ((children.first==0)!=(children.second==0))
      throw std::runtime_error("Generated relaxation requires binary junctions");
    if (children.first!=0)
      for (const auto child : {children.first,children.second}) {
        if (child>=count || seen[child] || tree.myVectParent[child]!=i)
          throw std::runtime_error("Invalid generated relaxation topology");
        seen[child]=true;order.push_back(child);
      }
  }
  if (order.size()+1!=count) throw std::runtime_error("Disconnected generated relaxation tree");
  std::vector<typename Tree::TPointD> original;
  for (const auto &segment : tree.myVectSegments) {
    if (!std::isfinite(segment.myCoordinate.norm())) return result;
    original.push_back(segment.myCoordinate);
  }
  tree.updateResistanceFromRoot();tree.updateRootRadius();
  const double motionBudget=2.*tree.myVectSegments[1].myRadius;
  if (!std::isfinite(motionBudget) || motionBudget<=0.) return result;
  for (unsigned int sweep=0;;++sweep) {
    tree.updateResistanceFromRoot();tree.updateRootRadius();
    bool valid=true;
    for (const auto i : order) {
      const auto &segment=tree.myVectSegments[i];
      const auto &parent=tree.myVectSegments[tree.myVectParent[i]].myCoordinate;
      const double length=(segment.myCoordinate-parent).norm();
      if (!std::isfinite(length) || length<=0. || !std::isfinite(segment.myRadius) || segment.myRadius<=0.)
        return result;
      if (length<2.*segment.myRadius) valid=false;
    }
    result.sweeps=sweep;
    if (valid) { result.valid=true;return result; }
    if (sweep==maximumSweeps) return result;
    // Parent-first projection may shorten a child edge. Visiting the entire
    // generated tree and recomputing radii handles that coupling explicitly.
    for (const auto i : order) {
      auto &point=tree.myVectSegments[i].myCoordinate;
      const auto &parent=tree.myVectSegments[tree.myVectParent[i]].myCoordinate;
      const auto direction=point-parent;
      const double length=direction.norm(),required=2.*tree.myVectSegments[i].myRadius;
      if (!(length>0.) || !std::isfinite(length)) return result;
      if (length>=required) continue;
      const double guard=128.*std::numeric_limits<double>::epsilon()*
                         std::max({1.,point.norm(),parent.norm(),required});
      point=parent+direction*((required+guard)/length);
      const double moved=(point-original[i]).norm();
      result.maximumDisplacement=std::max(result.maximumDisplacement,moved);
      if (!std::isfinite(moved) || moved>motionBudget+guard) return result;
    }
  }
}
}
