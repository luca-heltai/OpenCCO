#include "CoronaryArteryTree.h"

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

int main() {
  try {
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
