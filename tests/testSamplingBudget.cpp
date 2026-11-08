#include "CoronaryArteryTree.h"
#include "AnatomicalGrowth.h"
#include <set>

#include <iostream>
#include <stdexcept>

struct BudgetDomain {
  using TPoint = DGtal::PointVector<3, double>;
  unsigned int calls = 0;
  unsigned int reject = 100;
  unsigned int guard = 8;
  TPoint accepted{10., 10., 10.};
  TPoint randomPoint() {
    if (++calls > guard) throw std::runtime_error("Guard: sampling exceeded total budget");
    return calls <= reject ? TPoint(0., 0., 0.) : accepted;
  }
  bool isInside(const TPoint &) { return true; }
  bool checkNoIntersectDomain(const TPoint &, const TPoint &) { return true; }
};

void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}

template<int dim> struct FixedCandidateDomain {
  using TPoint = DGtal::PointVector<dim,double>;
  unsigned int calls=0;
  TPoint randomPoint() {
    ++calls;
    TPoint p; p[0]=2.5; p[1]=std::sqrt(18.75); return p;
  }
  bool isInside(const TPoint &) { return true; }
  bool checkNoIntersectDomain(const TPoint &,const TPoint &) { return true; }
};

template<int dim> void testLastPlacementSearch() {
  FixedCandidateDomain<dim> domain;
  CoronaryArteryTree<FixedCandidateDomain<dim>,dim> tree(domain);
  tree.myVectSegments.resize(2);
  tree.myVectSegments[0].myKTerm=0;
  tree.myVectSegments[1].myIndex=1;
  tree.myVectSegments[1].myCoordinate[0]=5.;
  tree.myVectSegments[1].myFlow=.001;
  tree.myVectParent={0,0}; tree.myVectChildren={{0,0},{0,0}};
  tree.myVectTerminals={1}; tree.my_NTerm=2; tree.my_qTerm=.001;
  tree.my_qPerf=.002; tree.my_pDrop=400.; tree.my_mu=.00371;
  tree.my_pTerm=0.; tree.my_pPerf=400.; tree.myRsupp=1.;
  tree.updateResistanceFromRoot(); tree.updateRootRadius();
  require(ExpandTreeHelpers::expandTree(tree,false,1,1) &&
          tree.myKTerm==2 && domain.calls==1,
          "Success on the sole/final placement search was not reported");
}

void testTerminalBookkeeping() {
  BudgetDomain domain;
  auto tree = OpenCCO::seedAnatomicalSubtree(domain, BudgetDomain::TPoint(),
      BudgetDomain::TPoint(5.,0.,0.), {2,1000.,.002,400.,.00371});
  require(tree.isAddable(BudgetDomain::TPoint(0.,5.,0.),1,100,.01,20,false),
          "Controlled leaf bifurcation did not succeed");
  std::set<unsigned int> leaves;
  for (unsigned int i=1; i<tree.myVectSegments.size(); ++i)
    if (tree.myVectChildren[i] == std::make_pair(0u,0u)) leaves.insert(i);
  require(leaves == std::set<unsigned int>(tree.myVectTerminals.begin(),tree.myVectTerminals.end()),
          "Bifurcation retained an internal node and lost the copied tissue leaf");
  require(tree.myVectTerminals.size() == tree.myKTerm, "Terminal index count mismatch");
}

void testEveryTerminalAndSegment() {
  BudgetDomain domain;
  domain.reject = 0;
  CoronaryArteryTree<BudgetDomain,3> tree(domain);
  const std::vector<BudgetDomain::TPoint> points{{-2.,0.,0.},{-1.,0.,0.},{0.,2.,0.},
                                               {0.,-2.,0.},{2.,-3.,0.},{-2.,-3.,0.}};
  tree.myVectSegments.resize(points.size());
  for (unsigned int i=0; i<points.size(); ++i) {
    tree.myVectSegments[i].myCoordinate = points[i];
    tree.myVectSegments[i].myIndex = i;
    tree.myVectSegments[i].myRadius = .1;
  }
  tree.myVectParent = {0,0,1,1,3,3};
  tree.myVectChildren = {{0,0},{2,3},{0,0},{4,5},{0,0},{0,0}};
  tree.myVectTerminals = {2,4,5};
  tree.myKTerm = 3;
  domain.accepted = points[2];
  require(!tree.generateALocation(.5).second, "Repeated first terminal was sampled");
  domain.accepted = {-1.5,0.,0.};
  require(!tree.generateALocation(.5).second, "Proximal segment was skipped after terminal loop");
  domain.accepted = {0.,2.3,0.};
  require(!tree.generateALocation(.5).second, "First-terminal spacing threshold was ignored");
}

int main() {
  try {
    testLastPlacementSearch<2>();
    testLastPlacementSearch<3>();
    testTerminalBookkeeping();
    testEveryTerminalAndSegment();
    BudgetDomain domain;
    CoronaryArteryTree<BudgetDomain, 3> tree(domain);
    tree.myVectSegments.resize(3);
    for (unsigned int i = 0; i < 3; ++i) {
      tree.myVectSegments[i].myIndex = i;
      tree.myVectSegments[i].myCoordinate = {0., 0., double(i == 0)};
    }
    tree.myVectParent = {0, 0, 1};
    tree.myVectChildren = {{0, 0}, {0, 0}, {0, 0}};
    tree.myVectTerminals = {1, 2};
    bool exhausted = false;
    try { tree.generateNewLocation(3); }
    catch (const std::runtime_error &) { exhausted = true; }
    require(exhausted && domain.calls == 3, "Sampler exceeded the requested total budget");
    domain.calls = 0;
    domain.reject = 2;
    const auto point = tree.generateNewLocation(3);
    require(domain.calls == 3 && point == BudgetDomain::TPoint(10., 10., 10.),
            "Success on the final allowed sample was lost");
    domain.calls = 0;
    exhausted = false;
    try { tree.generateNewLocation(0); }
    catch (const std::runtime_error &) { exhausted = true; }
    require(exhausted && domain.calls == 0, "Zero budget must not sample");
    domain.calls = 0;
    const auto empty = tree.tryGenerateNewLocation(0);
    require(!empty.second && domain.calls == 0, "Nonthrowing zero-budget result failed");
    domain.reject = 0;
    domain.guard = 120;
    domain.accepted = {1.55, 0., 0.};
    require(tree.tryGenerateNewLocation(101).second && domain.calls == 101,
            "Finite adaptive-spacing schedule changed");
    domain.calls = 0;
    tree.my_NTerm = 2;
    require(!ExpandTreeHelpers::expandTree(tree, false, 0, 1) && domain.calls == 0,
            "Zero outer-search budget did not stop expansion");
    std::cout << "Finite sampling budgets and final-attempt success passed\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
