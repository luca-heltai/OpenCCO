// Copyright (C) 2026 by Luca Heltai
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// Adapted to DGtal points from ImmersX growth_geometry.h at
// c873cd6081a0aa1185564face956cbb7cbebd7fd, segment_distance_square/compatible.
// https://github.com/luca-heltai/immersx/blob/c873cd6081a0aa1185564face956cbb7cbebd7fd/source/vascularization/growth_geometry.h
#pragma once

#include "DGtal/helpers/StdDefs.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace FrustumGeometry {
using Point = DGtal::PointVector<3, double>;
inline Point cross(const Point &a, const Point &b) {
  return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}

inline double distanceSquare(const Point &a0, const Point &a1,
                              const Point &b0, const Point &b1, Point &closest) {
  const auto u = a1-a0, v = b1-b0, w = a0-b0;
  const double vv = v.dot(v), uv = u.dot(v), wv = w.dot(v);
  std::vector<double> cuts{0., 1.};
  if (uv != 0.)
    for (const double t : {-wv/uv, (vv-wv)/uv})
      if (t > 0. && t < 1.) cuts.push_back(t);
  std::sort(cuts.begin(), cuts.end());
  double best = std::numeric_limits<double>::infinity();
  for (std::size_t i = 1; i < cuts.size(); ++i) {
    const double middle = .5*(cuts[i-1]+cuts[i]);
    const double projection = (wv+middle*uv)/vv;
    auto offset = w, slope = u;
    if (projection >= 1.) offset -= v;
    else if (projection > 0. && projection < 1.) {
      offset -= v*(wv/vv);
      slope -= v*(uv/vv);
    }
    const double norm = slope.dot(slope);
    const double t = norm > 0. ? std::clamp(-offset.dot(slope)/norm, cuts[i-1], cuts[i]) : middle;
    const auto delta = offset+slope*t;
    const double square = delta.dot(delta);
    if (square < best) { best = square; closest = delta; }
  }
  return best;
}

// True means a conservative separation witness was found for the two closed,
// flat-ended linear-taper frusta. False is not proof of intersection: preserve
// conservative rejection when this finite axis set cannot certify separation.
inline bool separated(const Point &a0, const Point &a1, std::array<double, 2> ra,
                       const Point &b0, const Point &b1, std::array<double, 2> rb,
                       double clearance = 0.) {
  const auto da = a1-a0, db = b1-b0;
  const double la = da.norm(), lb = db.norm();
  if (!(la > 0.) || !(lb > 0.) || !std::isfinite(la) || !std::isfinite(lb) ||
      !std::isfinite(clearance) || clearance < 0.) return false;
  for (const double radius : {ra[0], ra[1], rb[0], rb[1]})
    if (!std::isfinite(radius) || radius <= 0.) return false;
  const double radiusA = std::max(ra[0], ra[1]), radiusB = std::max(rb[0], rb[1]);
  const double guard = 64.*std::numeric_limits<double>::epsilon()*
                       std::max({1., a0.norm(), a1.norm(), b0.norm(), b1.norm()});
  Point closest;
  const double square = distanceSquare(a0, a1, b0, b1, closest);
  const double capsule = radiusA+radiusB+clearance+guard;
  if (std::isfinite(square) && square > capsule*capsule) return true;
  const Point ua = da/la, ub = db/lb;
  const double condition = 1.+radiusA/la+radiusB/lb;
  const double required = clearance + std::max(guard, 64.*std::numeric_limits<double>::epsilon()*
                                              std::max({1., radiusA, radiusB, clearance}))*condition;
  if (!std::isfinite(required)) return false;
  const auto witness = [&](Point axis) {
    const double norm = axis.norm();
    if (!(norm > 0.) || !std::isfinite(norm)) return false;
    axis /= norm;
    const double radialA = cross(axis, ua).norm(), radialB = cross(axis, ub).norm();
    const double a = axis.dot(da), b = axis.dot(b0-a0), c = axis.dot(b1-a0);
    const double lowA = std::min(-ra[0]*radialA, a-ra[1]*radialA);
    const double highA = std::max(ra[0]*radialA, a+ra[1]*radialA);
    const double lowB = std::min(b-rb[0]*radialB, c-rb[1]*radialB);
    const double highB = std::max(b+rb[0]*radialB, c+rb[1]*radialB);
    return std::isfinite(lowA) && std::isfinite(highA) && std::isfinite(lowB) && std::isfinite(highB) &&
           std::max(lowA-highB, lowB-highA) > required;
  };
  if (witness(da) || witness(db) || witness(cross(ua, ub)) || witness(closest)) return true;
  for (const auto &a : {a0, a1})
    for (const auto &b : {b0, b1}) if (witness(a-b)) return true;
  const auto lateral = [](const Point &first, const Point &direction, double length,
                          const std::array<double, 2> &radii, const Point &p) {
    const auto delta = p-first;
    const auto radial = delta-direction*delta.dot(direction);
    const double norm = radial.norm();
    if (!(norm > 0.) || !std::isfinite(norm)) return Point();
    return radial/norm-direction*((radii[1]-radii[0])/length);
  };
  for (const auto &p : {b0, b1}) if (witness(lateral(a0, ua, la, ra, p))) return true;
  for (const auto &p : {a0, a1}) if (witness(lateral(b0, ub, lb, rb, p))) return true;
  return false;
}
} // namespace FrustumGeometry
