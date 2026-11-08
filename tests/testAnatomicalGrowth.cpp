#include "AnatomicalGrowth.h"
#include "GeneratedSubtreeState.h"
#include "GeneratedGeometry.h"
#include <sstream>
#include "ExpandTreeHelpers.h"
#include "FrustumGeometry.h"
#include <iostream>

void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}

void testGeneratedLengthRelaxation() {
  using Point=DGtal::PointVector<3,double>;
  CircularDomainCtrl<3> domain(100.,Point());
  auto tree=OpenCCO::seedAnatomicalSubtree(domain,Point(),Point(.3,0.,0.),{2,1000.,1.,400.,.00371});
  tree.myVectSegments.resize(4);tree.myVectParent={0,0,1,1};
  tree.myVectChildren={{0,0},{2,3},{0,0},{0,0}};tree.myVectTerminals={2,3};tree.myKTerm=2;
  tree.myVectSegments[1].myKTerm=2;tree.myVectSegments[1].myFlow=1.;
  for (unsigned int i=2;i<4;++i) {
    tree.myVectSegments[i].myIndex=i;
    tree.myVectSegments[i].myCoordinate=Point(2.,i==2 ? 1. : -1.,0.);
    tree.myVectSegments[i].myFlow=.5;
  }
  tree.myFixedRootRadius=.2;
  const auto original=tree;
  const auto exhausted=OpenCCO::relaxGeneratedLengths(tree,0);
  require(!exhausted.valid && tree.myVectSegments[1].myCoordinate==Point(.3,0.,0.),
          "Zero repair budget moved generated geometry");
  const auto fixed=OpenCCO::relaxGeneratedLengths(tree,32);
  require(fixed.valid && fixed.sweeps>0 && fixed.maximumDisplacement<=.4+1e-12 &&
          tree.myVectSegments[0].myCoordinate==Point(),
          "Bounded repair moved its anatomical root or exceeded its motion budget");
  const double kappa=8.*.00371/M_PI;
  const auto check=[&](const auto &state) {
    for (unsigned int i=1;i<4;++i)
      require((state.myVectSegments[i].myCoordinate-state.myVectSegments[state.myVectParent[i]].myCoordinate).norm()
                >=2.*state.myVectSegments[i].myRadius,"Repair left an invalid actual diameter");
    const double rootLength=state.myVectSegments[1].myCoordinate.norm();
    const double childLength=(state.myVectSegments[2].myCoordinate-state.myVectSegments[1].myCoordinate).norm();
    const double pressure=kappa*rootLength/std::pow(state.myVectSegments[1].myRadius,4)+
                          .5*kappa*childLength/std::pow(state.myVectSegments[2].myRadius,4);
    require(std::abs(pressure-state.my_pDrop)<1e-12,"Repaired binary pressure violates independent Poiseuille oracle");
  };
  check(tree);
  auto pressureSized=original;
  pressureSized.myFixedRootRadius=0.;pressureSized.my_pDrop=.5;pressureSized.my_pPerf=.5;
  const auto canonical=OpenCCO::relaxGeneratedLengths(pressureSized,32);
  require(canonical.valid && pressureSized.my_pDrop==.5 &&
          pressureSized.myVectSegments[0].myCoordinate==Point(),
          "Canonical repair changed pressure input or anatomical root");
  check(pressureSized);
}

void testAcceptedCandidateScoring() {
  using Point=DGtal::PointVector<3,double>;
  CircularDomainCtrl<3> domain(50.,Point());
  auto tree=OpenCCO::seedAnatomicalSubtree(domain,Point(),Point(5.,0.,0.),{3,1000.,.003,400.,.00371});
  require(ExpandTreeHelpers::tryAddPoint(tree,Point(2.5,std::sqrt(18.75),0.),false,
                                       [](const auto &) { return true; }),"Scoring fixture initialization");
  unsigned int calls=0;
  double maximum=0., minimum=std::numeric_limits<double>::infinity();
  require(ExpandTreeHelpers::tryAddPoint(tree,Point(2.5,1.4,5.),false,[&](auto &candidate) {
    ++calls;
    const double before=candidate.computeTotalVolume(1);
    maximum=std::max(maximum,before); minimum=std::min(minimum,before);
    // Invert the original ordering: stale-volume scoring must pick wrongly.
    const double scale=std::sqrt((1e-8/before)/before);
    for (unsigned int i=1; i<candidate.myVectSegments.size(); ++i)
      candidate.myVectSegments[i].myRadius*=scale;
    candidate.myVectSegments[0].myRadius=before; // sentinel-only selection witness
    return true;
  }),"Scoring fixture has no admissible third site");
  require(calls>1 && maximum-minimum>1e-12 && tree.myVectSegments[0].myRadius==maximum,
          "Selection used stale pre-callback volume instead of accepted candidate score");
  require(std::abs(tree.computeTotalVolume(1)-1e-8/maximum)<1e-18,
          "Selected physical volume disagrees with callback scoring oracle");
}

void testNativeSeparationCertificate() {
  using Point = DGtal::PointVector<3,double>;
  CircularDomainCtrl<3> domain(10.,Point());
  auto tree = OpenCCO::seedAnatomicalSubtree(domain,Point(-1.,0.,0.),Point(),{1,1000.,1.,400.,.00371});
  tree.myVectSegments[1].myRadius = 1.;
  const Point a(.05,0.,0.), b(1.,0.,0.);
  require(tree.isIntersectingTree(a,b,1.,std::make_tuple(-1,-1,-1)),
          "Legacy rounded-capsule behavior changed");
  tree.mySegmentSeparationCertificate = [](const Point &p, const Point &q, double r,
                                           const Point &s, const Point &t, double u) {
    return FrustumGeometry::separated(p,q,{{r,r}},s,t,{{u,u}});
  };
  require(!tree.isIntersectingTree(a,b,1.,std::make_tuple(-1,-1,-1)),
          "Native precheck discarded an analytic flat-cap separating plane");
  require(tree.isIntersectingTree(Point(-.5,-2.,0.),Point(-.5,2.,0.),1.,std::make_tuple(-1,-1,-1)),
          "Native separation hook admitted crossing cylinders");
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
  tree.my_qPerf *= 2.;
  tree.myVectSegments[1].myFlow = 2.;
  tree.myVectSegments[2].myFlow = tree.myVectSegments[3].myFlow = 1.;
  tree.updateResistanceFromRoot();
  tree.updateRootRadius();
  require(tree.myVectSegments[1].myRadius == .2 && std::abs(tree.my_pDrop-2.*expected)<1e-12,
          "Measured lumen was resized instead of deriving flow-dependent pressure");
  OpenCCO::AnatomicalTree source;
  source.role="portal";
  source.nodes={{{{-.001,0.,0.}},.0002,true,false},
                {{{0.,0.,0.}},.0002,false,true}};
  source.edges={{{0,1}}};
  std::stringstream snapshot;
  OpenCCO::writeGeneratedSubtree(snapshot,source,1,tree,Point(),.001);
  const auto encoding=snapshot.str();
  auto restored=OpenCCO::readGeneratedSubtree(snapshot,domain,source,1,Point(),.001,
                                             {4,1000.,4.,400.,.00371},true);
  require(restored.myKTerm==2 && restored.my_NTerm==4 && restored.my_qTerm==1. &&
          restored.myVectTerminals==std::vector<unsigned int>({2,3}),
          "Continuation lost topology, true leaves or explicit current sizing inputs");
  for (unsigned int i=1; i<tree.myVectSegments.size(); ++i)
    require(restored.myVectSegments[i].myCoordinate==tree.myVectSegments[i].myCoordinate &&
            restored.myVectSegments[i].myRadius==tree.myVectSegments[i].myRadius,
            "Measured-radius geometry continuation changed its hydraulic shape");
  require(restored.my_pDrop==tree.my_pDrop,"Continuation changed the analytic binary pressure");
  const auto mustReject=[&](std::string text,bool matched) {
    bool rejected=false;
    try {
      std::istringstream input(text);
      (void)OpenCCO::readGeneratedSubtree(input,domain,source,1,Point(),.001,
                                         {4,1000.,4.,400.,.00371},matched);
    } catch (const std::runtime_error &) { rejected=true; }
    require(rejected,"Malformed or model-incompatible continuation was accepted");
  };
  mustReject(encoding,false);
  std::istringstream explicitResize(encoding);
  auto resized=OpenCCO::readGeneratedSubtree(explicitResize,domain,source,1,Point(),.001,
                                            {4,1000.,4.,400.,.00371},false,true);
  require(resized.myFixedRootRadius==0. && resized.my_pDrop==400.,
          "Explicit resize did not restore canonical prescribed-drop normalization");
  const double oracleResistance=kappa*(1.+std::cbrt(2.)*std::sqrt(2.));
  require(std::abs(resized.myVectSegments[1].myRadius-std::pow(2.*oracleResistance/400.,.25))<1e-14,
          "Explicit resize disagrees with independent Poiseuille radius");
  OpenCCO::GrownAnatomy resizedNetwork(source);
  resizedNetwork.append(1,resized,Point(),.001);
  require(resizedNetwork.branches()[0].radiusB==source.nodes[1].radius &&
          resizedNetwork.branches()[1].radiusA!=source.nodes[1].radius,
          "Explicit generated resizing overwrote or averaged the imported radius trace");
  mustReject(encoding+"unexpected",true);
  auto moved=encoding;
  moved.replace(moved.find("\n0 0 ")+5,1,".001");
  mustReject(moved,true);
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
    testNativeSeparationCertificate();
    testAcceptedCandidateScoring();
    testGeneratedLengthRelaxation();
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
