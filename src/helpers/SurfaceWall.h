#pragma once

#include "DGtal/helpers/StdDefs.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef OPENCCO_WITH_VTK
#include <vtkFeatureEdges.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkSTLReader.h>
#include <vtkSmartPointer.h>
#endif

namespace OpenCCO {
/** Optional original-surface wall certificate; the existing mask still samples
 * tissue and runs the native optimizer. No new tissue-domain controller. */
class SurfaceWall {
public:
  using Point = DGtal::PointVector<3, double>;
  enum class Result { inside, outside, unresolved };
  explicit SurfaceWall(const std::string &path) {
#ifdef OPENCCO_WITH_VTK
    auto reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(path.c_str());
    reader->MergingOn();
    reader->Update();
    if (reader->GetErrorCode() || reader->GetOutput()->GetNumberOfCells() == 0)
      throw std::runtime_error("Cannot read wall surface");
    auto edges = vtkSmartPointer<vtkFeatureEdges>::New();
    edges->SetInputConnection(reader->GetOutputPort());
    edges->BoundaryEdgesOn(); edges->NonManifoldEdgesOn();
    edges->FeatureEdgesOff(); edges->ManifoldEdgesOff();
    edges->Update();
    if (edges->GetOutput()->GetNumberOfCells() != 0)
      throw std::runtime_error("Wall surface must be closed and manifold");
    auto normals = vtkSmartPointer<vtkPolyDataNormals>::New();
    normals->SetInputConnection(reader->GetOutputPort());
    normals->ConsistencyOn(); normals->AutoOrientNormalsOn(); normals->SplittingOff();
    normals->Update();
    distance = vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
    distance->SetInput(normals->GetOutput());
#else
    (void)path;
    throw std::runtime_error("Original-surface walls require a build with VTK");
#endif
  }

  double clearance(const Point &p) const {
#ifdef OPENCCO_WITH_VTK
    double x[3] = {p[0], p[1], p[2]};
    const double value = -distance->EvaluateFunction(x);
    if (!std::isfinite(value)) throw std::runtime_error("Nonfinite surface distance");
    return value;
#else
    (void)p;
    throw std::runtime_error("Surface distance unavailable without VTK");
#endif
  }

  /** Certify the entire constant-radius capsule, not just sampled centres.
   * Signed distance is 1-Lipschitz: each subinterval is bounded by its nearest
   * endpoint distance minus half its length. Unresolved budgets stay rejected. */
  Result tubeInside(const Point &a, const Point &b, const double radius,
                    const unsigned int budget = 4096) const {
    if (!std::isfinite(radius) || radius <= 0.) throw std::runtime_error("Invalid wall radius");
    if (budget < 2) throw std::runtime_error("Wall certificate needs at least two evaluations");
    for (unsigned int j=0; j<3; ++j)
      if (!std::isfinite(a[j]) || !std::isfinite(b[j])) throw std::runtime_error("Invalid wall point");
    const double guard = 64.*std::numeric_limits<double>::epsilon()*(1.+std::max(a.norm(), b.norm()));
    struct Interval { Point a, b; double da, db; };
    const double da = clearance(a), db = clearance(b);
    std::vector<Interval> pending{{a, b, da, db}};
    unsigned int evaluations = 2;
    while (!pending.empty()) {
      const auto part = pending.back(); pending.pop_back();
      if (std::min(part.da, part.db) <= radius+guard) return Result::outside;
      if (std::min(part.da, part.db)-.5*(part.b-part.a).norm() > radius+guard) continue;
      if (evaluations >= budget) return Result::unresolved;
      const auto middle = (part.a+part.b)*.5;
      const double dm = clearance(middle); ++evaluations;
      pending.push_back({middle, part.b, dm, part.db});
      pending.push_back({part.a, middle, part.da, dm});
    }
    return Result::inside;
  }
private:
#ifdef OPENCCO_WITH_VTK
  vtkSmartPointer<vtkImplicitPolyDataDistance> distance;
#endif
};
}
