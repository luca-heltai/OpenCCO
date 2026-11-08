#include "AnatomicalGrowth.h"
#include "ExpandTreeHelpers.h"
#include <iostream>

void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}

void testMeasuredRadius() {
  using Point = DGtal::PointVector<3, double>;
  CircularDomainCtrl<3> domain(10., Point());
  auto tree = OpenCCO::seedAnatomicalSubtree(domain, Point(), Point(1.,0.,0.),
                                            {2,1000.,1.,400.,.00371});
  tree.myKTerm = 2;
  tree.myVectSegments.resize(4);
  tree.myVectParent = {0,0,1,1};
  tree.myVectChildren = {{0,0},{2,3},{0,0},{0,0}};
  tree.myVectTerminals = {2,3};
  tree.myVectSegments[1].myKTerm = 2;
  tree.myVectSegments[1].myFlow = 1.;
  for (unsigned int i=2; i<4; ++i) {
    tree.myVectSegments[i].myIndex = i;
    tree.myVectSegments[i].myCoordinate = Point(2.,i==2 ? 1. : -1.,0.);
    tree.myVectSegments[i].myFlow = .5;
  }
  tree.myFixedRootRadius = .2;
  tree.updateResistanceFromRoot();
  tree.updateRootRadius();
  const double child = .2/std::cbrt(2.);
  const double kappa = 8.*.00371/M_PI;
  const double expected = kappa/std::pow(.2,4) + .5*kappa*std::sqrt(2.)/std::pow(child,4);
  require(tree.myVectSegments[1].myRadius == .2 &&
          std::abs(tree.myVectSegments[2].myRadius-child)<1e-14,
          "Measured inlet or analytic Murray child radius changed");
  require(std::abs(tree.my_pDrop-expected)<1e-12 &&
          std::abs(tree.my_pPerf-tree.my_pTerm-expected)<1e-12,
          "Derived binary pressure disagrees with independent Poiseuille path");
  tree.my_qTerm *= 2.;
  tree.myVectSegments[1].myFlow = 2.;
  tree.myVectSegments[2].myFlow = tree.myVectSegments[3].myFlow = 1.;
  tree.updateResistanceFromRoot();
  tree.updateRootRadius();
  require(tree.myVectSegments[1].myRadius == .2 && std::abs(tree.my_pDrop-2.*expected)<1e-12,
          "Measured lumen was resized instead of deriving flow-dependent pressure");
  tree.myFixedRootRadius = -1.;
  bool rejected=false;
  try { tree.updateRootRadius(); } catch(const std::runtime_error &) { rejected=true; }
  require(rejected,"Invalid radius constraint accepted");
  tree.myFixedRootRadius = 0.;
  tree.my_pDrop = 400.;
  tree.updateRootRadius();
  const double resistance = kappa*(1.+std::cbrt(2.)*std::sqrt(2.));
  require(std::abs(tree.myVectSegments[1].myRadius-std::pow(2.*resistance/400.,.25))<1e-14,
          "Default fixed-pressure behavior changed");
}

int main() {
  try {
    testMeasuredRadius();
    OpenCCO::AnatomicalTree source;
    source.role = "portal";
    source.nodes = {{{{.1, .2, .3}}, .001, true, false},
                    {{{.102, .2, .3}}, .0008, false, true}};
    source.edges = {{{0, 1}}};
    const auto original = source;
    using Point = DGtal::PointVector<3, double>;
    CircularDomainCtrl<3> domain(100., Point());
    const Point origin(.05, .1, .2);
    const double pitch = .002;
    const Point root = (Point(.102, .2, .3)-origin)/pitch;
    const OpenCCO::SubtreeParameters parameters{2, 1000., 1., 400., .00371};
    auto tree = OpenCCO::seedAnatomicalSubtree(domain, root, root+Point(3., 0., 0.), parameters);
    const double expectedRadius = std::pow(8.*parameters.viscosityPaSecond*3./M_PI*
                                           (parameters.flowVoxelCubedPerSecond/2.)/400., .25);
    require(std::abs(tree.myVectSegments[1].myRadius-expectedRadius) < 1e-14,
            "Seed radius disagrees with independent Poiseuille calculation");
    require(tree.getLengthSegment(1) == 3. && tree.myLengthFactor == 1., "Physical length seed changed");
    OpenCCO::GrownAnatomy grown(source);
    grown.append(1, tree, origin, pitch);
    require(grown.nodes()[0] == original.nodes[0].position && grown.nodes()[1] == original.nodes[1].position,
            "Imported coordinates changed");
    const auto &edge = grown.branches()[0];
    require(edge.a == 0 && edge.b == 1 && !edge.generated &&
            edge.radiusA == .001 && edge.radiusB == .0008, "Imported endpoint traces changed");
    require(grown.nodes().size() == 3 && grown.branches().size() == 2 &&
            grown.tissueTerminals() == std::vector<std::size_t>{2}, "Append topology incorrect");
    require(grown.source().nodes[1].grow && grown.source().root() == 0, "Source metadata changed");
    bool rejected = false;
    try { grown.append(1, tree, origin, pitch); } catch (const std::runtime_error &) { rejected = true; }
    require(rejected, "Repeated attachment accepted");
    OpenCCO::GrownAnatomy unchanged(source);
    tree.myVectSegments[0].myCoordinate[0] += 1.;
    rejected = false;
    try { unchanged.append(1, tree, origin, pitch); } catch (const std::runtime_error &) { rejected = true; }
    require(rejected && unchanged.nodes().size() == 2 && unchanged.branches().size() == 1,
            "Invalid attachment changed materialized source");

    CircularDomainCtrl<3> sampling(50., Point());
    auto candidate = OpenCCO::seedAnatomicalSubtree(sampling, Point(0., 25., 0.), Point(),
                              {2, 20000., 1., 400., .00371});
    candidate.myCandidateDiagnostics = std::make_shared<std::map<std::string, std::size_t>>();
    unsigned int calls = 0;
    bool sentinelStable = true;
    srand(23);
    const bool complete = ExpandTreeHelpers::expandTree(candidate, false, 10, 1000,
                                                       [&](auto &state) {
                                                         ++calls;
                                                         sentinelStable = sentinelStable && state.myVectParent[0] == 0;
                                                         return false;
                                                       });
    require(!complete && calls > 0 && candidate.myVectSegments.size() == 2 && candidate.myKTerm == 1,
            "Rejected candidate was committed or acceptance hook was not tested");
    require(sentinelStable, "Bifurcation overwrote the special root parent sentinel");
    require((*candidate.myCandidateDiagnostics)["accepted"] > 0,
            "Shared native telemetry lost decisions from copied candidates");
    std::cout << "Immutable anatomical append, physical seed and candidate gate passed\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
