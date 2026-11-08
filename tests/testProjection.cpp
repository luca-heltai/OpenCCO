#include "GeomHelpers.h"
#include <iostream>
#include <limits>
#include <stdexcept>

void require(bool value, const char *message) { if (!value) throw std::runtime_error(message); }

template<int dim> void testAxis() {
  using Point = DGtal::PointVector<dim,double>;
  Point a, b, c, expected, projected;
  a[0]=2.; c[0]=expected[0]=1.; c[1]=3.;
  require(GeomHelpers::projectOnStraightLine(a,b,c,projected),"Reverse-axis interior projection rejected");
  require((projected-expected).norm()<1e-14,"Reverse-axis projection is wrong");
  require(GeomHelpers::projectOnStraightLine(b,a,c,projected),"Forward-axis interior projection rejected");
  require((projected-expected).norm()<1e-14,"Forward-axis projection is wrong");
  c[0]=3.; expected[0]=3.;
  require(!GeomHelpers::projectOnStraightLine(a,b,c,projected),"Exterior projection accepted");
  require((projected-expected).norm()<1e-14,"Infinite-line projection output was clamped");
  require(!GeomHelpers::projectOnStraightLine(b,a,c,projected),"Reversed exterior projection accepted");
  require(GeomHelpers::projectOnStraightLine(a,a,c,projected) && projected==a,
          "Degenerate point segment must have a finite endpoint projection");
}

int main() {
  try {
    testAxis<2>(); testAxis<3>();
    using Point = DGtal::PointVector<3,double>;
    const Point a(0.,2.,3.), b(2.,0.,1.), c(4.,4.,2.), expected(1.,1.,2.);
    Point projected;
    require(GeomHelpers::projectOnStraightLine(a,b,c,projected) && (projected-expected).norm()<1e-14,
            "Mixed-direction transverse projection disagrees with analytic oracle");
    require(GeomHelpers::projectOnStraightLine(b,a,c,projected) && (projected-expected).norm()<1e-14,
            "Endpoint order changed mixed-direction projection");
    using Point2 = DGtal::PointVector<2,double>;
    Point2 decimalProjection;
    require(GeomHelpers::projectOnStraightLine(Point2(),Point2(10.,0.),Point2(5.3,15.3),decimalProjection) &&
            (decimalProjection-Point2(5.3,0.)).norm() <= 16.*std::numeric_limits<double>::epsilon()*10.,
            "Non-dyadic projection exceeds arithmetic roundoff bound");
    std::cout << "Analytic2D/3D segment projections passed\n";
    return 0;
  } catch(const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
