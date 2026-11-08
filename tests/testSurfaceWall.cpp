#include "SurfaceWall.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using Wall = OpenCCO::SurfaceWall;
using Point = Wall::Point;
void require(bool value, const char *message) { if (!value) throw std::runtime_error(message); }

int main() {
  try {
    const std::string fixtures = std::string(OPENCCO_TEST_SOURCE_DIR)+"/data/";
    Wall wall(fixtures+"unit_wall.stl");
    require(std::abs(wall.clearance(Point(.5,.5,.5))-.5)<1e-12,"Cube centre distance");
    require(std::abs(wall.clearance(Point(-.001,.5,.5))+.001)<1e-12,"Outside cube sign");
    require(std::abs(wall.clearance(Point(.0003,.005,.005))-.0003)<1e-12,"Near-wall distance");
    require(wall.tubeInside(Point(.005,.005,.005),Point(.02,.02,.02),.001)==Wall::Result::inside,"Interior capsule");
    // A0.2mm margin needs no submillimetre full-organ mask.
    const Point a(.0003,.005,.005), b(.0003,.025,.005);
    require(wall.tubeInside(a,b,.0001)==Wall::Result::inside,"Adaptive close-wall certificate");
    require(wall.tubeInside(a,b,.0001,2)==Wall::Result::unresolved,"Budget must not certify unresolved capsule");
    require(wall.tubeInside(a,b,.0004)==Wall::Result::outside,"Outside radius");
    bool invalid=false;
    try { wall.tubeInside(a,b,-1.); } catch(const std::runtime_error &) { invalid=true; }
    require(invalid,"Nonpositive radius accepted");
    Wall disconnected(fixtures+"two_wall_cubes.stl");
    const Point c(.005,.015,.015), d(.025,.015,.015);
    require(disconnected.clearance(c)>0. && disconnected.clearance(d)>0.,"Disconnected inside endpoints");
    require(disconnected.tubeInside(c,d,.001)==Wall::Result::outside,"Segment crosses exterior between components");
    bool open=false;
    try { Wall bad(fixtures+"open_wall.stl"); } catch(const std::runtime_error &) { open=true; }
    require(open,"Open surface accepted");
    std::cout << "All original-surface wall tests passed\n";
    return 0;
  } catch(const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
