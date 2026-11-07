#include "CoronaryArteryTree.h"
#include "ExpandTreeHelpers.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

// Reject immediately after observing the query; this tests coordinate transfer,
// not the optimizer's convergence or the anatomical acceptance predicates.
template <typename Component>
struct RecordingDomain {
  using TPoint = DGtal::PointVector<3, Component>;
  TPoint myCenter;
  TPoint observed;
  TPoint randomPoint() { return myCenter; }
  TPoint firtCandidatePoint() const { return myCenter; }
  bool isInside(const TPoint &point) { observed = point; return false; }
  bool checkNoIntersectDomain(const TPoint &, const TPoint &) { return false; }
};

void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}

int main() {
  static_assert(std::is_same<
      CoronaryArteryTree<ImageMaskDomainCtrl<3>, 3>::TDomainPoint,
      ImageMaskDomainCtrl<3>::TPointI>::value, "Mask queries must remain integer voxels");
  static_assert(std::is_same<
      CoronaryArteryTree<CircularDomainCtrl<3>, 3>::TDomainPoint,
      CircularDomainCtrl<3>::TPoint>::value, "Continuous queries must remain real points");
  try {
    RecordingDomain<double> domain;
    domain.myCenter = {0.123, 0.157, 0.219};
    using Tree = CoronaryArteryTree<RecordingDomain<double>, 3>;
    Tree tree(domain);
    tree.myVectSegments.resize(2);
    tree.myVectSegments[0].myCoordinate = {0.101, 0.131, 0.173};
    tree.myVectSegments[1].myCoordinate = {0.113, 0.149, 0.191};
    tree.myVectParent = {0, 0};
    tree.myVectChildren = {{0, 0}, {0, 0}};
    const Tree::TPointD candidate{0.137, 0.163, 0.227};
    require(!tree.isAddable(candidate, 1, 1, 0.01, 1, false), "Recording domain must reject");
    const auto expected = (candidate + tree.myVectSegments[0].myCoordinate +
                           tree.myVectSegments[1].myCoordinate) / 3.;
    require((domain.observed - expected).norm() < 1e-15,
            "Continuous candidate coordinates were truncated");
    const auto center = tree.getDomainCenter();
    for (unsigned int i = 0; i < 3; ++i)
      require(center[i] == domain.myCenter[i], "Continuous domain center was truncated");
    ExpandTreeHelpers::initFirtElemTree(tree);
    require(tree.myTreeCenter == domain.myCenter, "Initialization center was truncated");
    require(tree.myVectSegments[0].myCoordinate == domain.myCenter,
            "Root candidate coordinates changed");

    RecordingDomain<int> voxel;
    voxel.myCenter = {10, 20, 30};
    CoronaryArteryTree<RecordingDomain<int>, 3> integerTree(voxel);
    integerTree.myVectSegments.resize(2);
    integerTree.myVectSegments[0].myCoordinate = {10.1, 20.2, 30.3};
    integerTree.myVectSegments[1].myCoordinate = {11.2, 21.3, 31.4};
    integerTree.myVectParent = {0, 0};
    integerTree.myVectChildren = {{0, 0}, {0, 0}};
    require(!integerTree.isAddable({12.3, 22.4, 32.5}, 1, 1, .01, 1, false),
            "Recording voxel domain must reject");
    require(voxel.observed == RecordingDomain<int>::TPoint(11, 21, 31),
            "Integer domain query conversion changed");
    std::cout << "Continuous and integer domain coordinate contracts passed\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
