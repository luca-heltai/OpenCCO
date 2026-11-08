#include "AnatomicalGrowth.h"
#include "ExpandTreeHelpers.h"
#include <iostream>

void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}

int main() {
  try {
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
    std::cout << "Immutable anatomical append, physical seed and candidate gate passed\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
