// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "AnatomicalGrowth.h"
#include <iomanip>
#include <istream>
#include <ostream>

namespace OpenCCO {
/** Geometry continuation, not an opaque optimizer/RNG restart. Source binding
 * belongs to the caller; the anatomical CLI also compares source.anatomy.
 * Current requested flow/quota/viscosity determine restored generated physics. */
template<class Tree>
void writeGeneratedSubtree(std::ostream &out, const AnatomicalTree &source,
                          std::size_t attachment, const Tree &tree,
                          const typename Tree::TPointD &origin, double pitch) {
  GrownAnatomy checked(source);
  checked.append(attachment,tree,origin,pitch);
  out << std::setprecision(17) << "OpenCCOGeneratedSubtree 1\n" << source.role
      << " metres " << attachment << ' ' << tree.myVectSegments.size() << ' '
      << tree.my_gamma << ' ' << (tree.myFixedRootRadius>0.) << '\n';
  for (std::size_t i=0; i<tree.myVectSegments.size(); ++i) {
    out << i << ' ' << tree.myVectParent[i];
    for (unsigned int j=0; j<3; ++j)
      out << ' ' << (i==0 ? source.nodes[attachment].position[j] :
                            origin[j]+pitch*tree.myVectSegments[i].myCoordinate[j]);
    out << ' ' << (i==0 ? source.nodes[attachment].radius :
                          pitch*tree.myVectSegments[i].myRadius) << '\n';
  }
  if (!out) throw std::runtime_error("Cannot write generated subtree geometry");
}

template<class Domain>
CoronaryArteryTree<Domain,3>
readGeneratedSubtree(std::istream &in, Domain &domain, const AnatomicalTree &source,
                     std::size_t attachment,
                     const DGtal::PointVector<3,double> &origin, double pitch,
                     const SubtreeParameters &parameters, bool matchRadius) {
  using Tree=CoronaryArteryTree<Domain,3>;
  using Point=typename Tree::TPointD;
  std::string magic,role,units,extra;
  unsigned int version;
  std::size_t savedAttachment,count;
  double gamma;
  int matched;
  if (!(in>>magic>>version>>role>>units>>savedAttachment>>count>>gamma>>matched) ||
      magic!="OpenCCOGeneratedSubtree" || version!=1 || role!=source.role || units!="metres" ||
      attachment>=source.nodes.size() || !source.nodes[attachment].grow || savedAttachment!=attachment ||
      count<2 || count>2ull*parameters.terminals || count>std::numeric_limits<unsigned int>::max() || (matched!=0 && matched!=1) ||
      bool(matched)!=matchRadius || !std::isfinite(gamma) || gamma!=3. ||
      !std::isfinite(pitch) || pitch<=0. || !std::isfinite(origin.norm()))
    throw std::runtime_error("Incompatible generated subtree geometry header");
  std::vector<Point> points(count);
  std::vector<unsigned int> parents(count);
  std::vector<std::vector<unsigned int>> children(count);
  for (std::size_t i=0; i<count; ++i) {
    std::size_t id,parent;
    Point world;
    double savedRadius;
    if (!(in>>id>>parent>>world[0]>>world[1]>>world[2]>>savedRadius) || id!=i || parent>=count ||
        (i!=0 && parent==i) || !std::isfinite(world.norm()) ||
        !std::isfinite(savedRadius) || savedRadius<=0.)
      throw std::runtime_error("Invalid generated subtree geometry row");
    if (i==0) {
      const auto &root=source.nodes[attachment];
      if (parent!=0 || world!=Point(root.position[0],root.position[1],root.position[2]) ||
          savedRadius!=root.radius) throw std::runtime_error("Saved anatomical attachment changed");
    } else children[parent].push_back(static_cast<unsigned int>(i));
    parents[i]=static_cast<unsigned int>(parent);
    points[i]=(world-origin)/pitch;
    if (i!=0) // Restore integer proposals within arithmetic roundoff only.
      for (unsigned int j=0; j<3; ++j) {
        const double nearest=std::round(points[i][j]);
        if (std::abs(points[i][j]-nearest)<=64.*std::numeric_limits<double>::epsilon()*
                                              std::max(1.,std::abs(points[i][j])))
          points[i][j]=nearest;
      }
  }
  if (in>>extra) throw std::runtime_error("Trailing generated subtree geometry records");
  if (children[0]!=std::vector<unsigned int>{1}) throw std::runtime_error("Saved root must own segment1");
  std::vector<unsigned int> order{0};
  std::vector<bool> seen(count,false);seen[0]=true;
  for (std::size_t j=0; j<order.size(); ++j)
    for (const auto child : children[order[j]]) {
      if (seen[child]) throw std::runtime_error("Saved subtree contains a cycle");
      seen[child]=true;order.push_back(child);
    }
  if (order.size()!=count) throw std::runtime_error("Saved subtree is disconnected");
  std::set<std::array<double,3>> distinct;
  for (std::size_t i=0; i<count; ++i) {
    if (!distinct.insert(std::array<double,3>{{points[i][0],points[i][1],points[i][2]}}).second)
      throw std::runtime_error("Saved generated coordinates coincide");
    if (i!=0 && (points[i]-points[parents[i]]).norm()<=0.)
      throw std::runtime_error("Saved generated edge has zero length");
  }
  auto tree=seedAnatomicalSubtree(domain,points[0],points[1],parameters);
  tree.myVectSegments.resize(count);
  tree.myVectParent=parents;
  tree.myVectChildren.assign(count,{0,0});
  tree.myVectTerminals.clear();
  for (auto it=order.rbegin(); it!=order.rend(); ++it) {
    const auto i=*it;
    auto &segment=tree.myVectSegments[i];
    segment.myCoordinate=points[i];segment.myIndex=i;
    if (i==0) { segment.myKTerm=0;segment.myFlow=0.;continue; }
    if (children[i].empty()) { segment.myKTerm=1;tree.myVectTerminals.push_back(i); }
    else if (children[i].size()==2) {
      tree.myVectChildren[i]={children[i][0],children[i][1]};
      segment.myKTerm=tree.myVectSegments[children[i][0]].myKTerm+
                      tree.myVectSegments[children[i][1]].myKTerm;
    } else throw std::runtime_error("Saved generated junction is not binary");
    segment.myFlow=segment.myKTerm*tree.my_qTerm;
  }
  std::sort(tree.myVectTerminals.begin(),tree.myVectTerminals.end());
  tree.myKTerm=static_cast<unsigned int>(tree.myVectTerminals.size());
  if (count!=2ull*tree.myKTerm) throw std::runtime_error("Saved terminal count is inconsistent");
  tree.myFixedRootRadius=matchRadius ? source.nodes[attachment].radius/pitch : 0.;
  tree.updateResistanceFromRoot();tree.updateRootRadius();
  GrownAnatomy checked(source);
  checked.append(attachment,tree,origin,pitch);
  return tree;
}
}
