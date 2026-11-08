#include "CoronaryArteryTree.h"

#include <iostream>
#include <stdexcept>

void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}

struct GuardedCircle : CircularDomainCtrl<3> {
  using CircularDomainCtrl<3>::CircularDomainCtrl;
  unsigned int queries = 0;
  bool isInside(const TPoint &point) override {
    if (++queries > 10000) throw std::runtime_error("Translated sampler did not terminate");
    return CircularDomainCtrl<3>::isInside(point);
  }
};

int main() {
  try {
    using Point = CircularDomainCtrl<3>::TPoint;
    const Point center(1., 2., 3.);
    GuardedCircle sphere(.5, center);
    require(sphere.isInside(Point(1., 2.5, 3.)), "Closed sphere root was rejected");
    require(!sphere.isInside(Point(1., 2.6, 3.)), "Outside sphere point was accepted");
    require(sphere.firtCandidatePoint() == Point(1., 2.5, 3.),
            "First candidate ignored translated center");
    srand(23);
    for (unsigned int i = 0; i < 100; ++i) {
      const auto point = sphere.randomPoint();
      require((point - center).norm() <= .5, "Sample outside translated sphere");
    }
    SquareDomainCtrl<3> cube(.5, center);
    require(cube.isInside(Point(.5, 1.5, 2.5)), "Closed cube lower corner rejected");
    require(cube.isInside(Point(1.5, 2.5, 3.5)), "Closed cube upper corner rejected");
    require(!cube.isInside(Point(.49, 2., 3.)), "Cube lower bound ignored");
    require(!cube.isInside(Point(1.51, 2., 3.)), "Cube upper bound ignored");
    for (unsigned int i = 0; i < 100; ++i) {
      const auto point = cube.randomPoint();
      for (unsigned int j = 0; j < 3; ++j)
        require(std::abs(point[j] - center[j]) <= .5, "Sample outside translated cube");
    }
    std::cout << "Closed and translated implicit domains passed\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
