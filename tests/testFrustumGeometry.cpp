#include "FrustumGeometry.h"
#include <iostream>
#include <stdexcept>

void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}

int main() {
  try {
    using namespace FrustumGeometry;
    Point closest;
    require(std::abs(distanceSquare(Point(-1., 0., 0.), Point(1., 0., 0.),
                                    Point(0., -1., 1.), Point(0., 1., 1.), closest)-1.) < 1e-14,
            "Skew-segment interior distance disagrees with analytic oracle");
    // The plane x=.05 separates these flat caps although radius-one rounded
    // capsules overlap. No sampled intersection oracle is used here.
    require(separated(Point(-1., 0., 0.), Point(), {{1., 1.}},
                       Point(.1, 0., 0.), Point(1.1, 0., 0.), {{1., 1.}}),
            "Flat-cap separation certificate lost");
    require(!separated(Point(-1., 0., 0.), Point(), {{1., 1.}},
                        Point(.1, 0., 0.), Point(1.1, 0., 0.), {{1., 1.}}, .2),
            "Clearance larger than actual cap gap was waived");
    // Cone side: .45*x+y <= .55 for A; B's minimum .45*x+y is .81.
    require(separated(Point(-1., 0., 0.), Point(1., 0., 0.), {{1., .1}},
                       Point(.8, .5, 0.), Point(1.2, .5, 0.), {{.05, .05}}),
            "Analytic tapered-side separator lost");
    require(!separated(Point(-1., 0., 0.), Point(1., 0., 0.), {{1., .1}},
                        Point(.8, .15, 0.), Point(1.2, .15, 0.), {{.1, .1}}),
            "Intersecting tapered vessels accepted");
    require(!separated(Point(-1., 0., 0.), Point(1., 0., 0.), {{.1, .1}},
                        Point(0., -1., 0.), Point(0., 1., 0.), {{.1, .1}}),
            "Crossing cylinders accepted");
    require(!separated(Point(), Point(), {{.1, .1}},
                        Point(2., 0., 0.), Point(3., 0., 0.), {{.1, .1}}),
            "Degenerate vessel certified");
    const Point shift(1000., -2000., 3000.);
    require(separated(shift+Point(-.001, 0., 0.), shift, {{.001, .001}},
                       shift+Point(.0001, 0., 0.), shift+Point(.0011, 0., 0.), {{.001, .001}}),
            "Translated/scaled cap certificate lost");
    std::cout << "Analytic frustum separation and overlap oracles passed\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
