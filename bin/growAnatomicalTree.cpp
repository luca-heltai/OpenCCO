// SPDX-License-Identifier: GPL-3.0-or-later
#include "AnatomicalGrowth.h"
#include "ExpandTreeHelpers.h"
#include "CLI11.hpp"
#include "FrustumGeometry.h"
#include "SurfaceWall.h"
#include <functional>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <set>
#include <sstream>

using Point = DGtal::PointVector<3, double>;
using Mask = ImageMaskDomainCtrl<3>;
using Tree = CoronaryArteryTree<Mask, 3>;

namespace {
int systemId(const std::string &role) {
  return role == "arterial" ? 0 : role == "venous" ? 1 : 2;
}
Point point(const std::array<double, 3> &p) { return {p[0], p[1], p[2]}; }
OpenCCO::AnatomicalTree readAnatomy(const std::string &path) {
  std::ifstream stream(path);
  if (!stream) throw std::runtime_error("Cannot open anatomy: " + path);
  return OpenCCO::readAnatomicalTree(stream);
}
struct Obstacle {
  Point a, b;
  double radius;
  int system;
  long long sourceEdge;
  double radiusA, radiusB;
};

bool boxesOverlap(const Point &a, const Point &b, double radius,
                  const Point &c, const Point &d, double otherRadius) {
  for (unsigned int j = 0; j < 3; ++j)
    if (std::max(a[j], b[j])+radius < std::min(c[j], d[j])-otherRadius ||
        std::max(c[j], d[j])+otherRadius < std::min(a[j], b[j])-radius) return false;
  return true;
}

void writeNetwork(const std::filesystem::path &directory, const OpenCCO::GrownAnatomy &network) {
  std::filesystem::create_directories(directory);
  std::ofstream nodes(directory / "nodes.csv"), edges(directory / "segments.csv");
  if (!nodes || !edges) throw std::runtime_error("Cannot write materialized anatomy");
  nodes << std::setprecision(17) << "id,x,y,z,source_node,source_root,source_grow,active_tissue\n";
  const std::set<std::size_t> terminals(network.tissueTerminals().begin(), network.tissueTerminals().end());
  for (std::size_t i = 0; i < network.nodes().size(); ++i) {
    const auto &p = network.nodes()[i];
    const bool original = i < network.source().nodes.size();
    nodes << i << ',' << p[0] << ',' << p[1] << ',' << p[2] << ','
          << (original ? static_cast<long long>(i) : -1) << ','
          << (original && network.source().nodes[i].root) << ','
          << (original && network.source().nodes[i].grow) << ',' << terminals.count(i) << '\n';
  }
  edges << std::setprecision(17)
        << "system,id,a,b,x0,y0,z0,x1,y1,z1,r0,r1,generated,source_edge\n";
  for (std::size_t i = 0; i < network.branches().size(); ++i) {
    const auto &branch = network.branches()[i];
    const auto &a = network.nodes()[branch.a], &b = network.nodes()[branch.b];
    edges << systemId(network.source().role) << ',' << i << ',' << branch.a << ',' << branch.b
          << ',' << a[0] << ',' << a[1] << ',' << a[2] << ',' << b[0] << ',' << b[1] << ',' << b[2]
          << ',' << branch.radiusA << ',' << branch.radiusB << ',' << branch.generated << ','
          << (branch.generated ? -1 : static_cast<long long>(branch.sourceEdge)) << '\n';
  }
  if (!nodes || !edges) throw std::runtime_error("Materialized anatomy write failed");
}

struct Filter {
  Mask &mask;
  const std::vector<Obstacle> &fixed;
  int system;
  long long incident;
  Point outward;
  std::map<std::string, std::size_t> counts;
  std::vector<std::string> collisionTraces;
  std::function<OpenCCO::SurfaceWall::Result(const Point &, const Point &, double)> surfaceCertificate;
  bool reject(const std::string &reason) { ++counts[reason]; return false; }
  bool wall(const Point &a, const Point &b, double radius) {
    if (surfaceCertificate) {
      const auto result = surfaceCertificate(a, b, radius);
      if (result == OpenCCO::SurfaceWall::Result::unresolved) ++counts["wall_unresolved"];
      return result == OpenCCO::SurfaceWall::Result::inside;
    }
    const double length = (b-a).norm();
    const auto steps = std::max(1u, static_cast<unsigned int>(std::ceil(length/.5)));
    const double required = radius + std::sqrt(3.)/2. + length/(2.*steps);
    for (unsigned int i = 0; i <= steps; ++i) {
      const Point p = a + (b-a)*(double(i)/steps);
      Mask::TPointI voxel;
      for (unsigned int j = 0; j < 3; ++j) voxel[j] = static_cast<int>(std::lround(p[j]));
      if (!mask.myImage.domain().isInside(voxel) || mask.myImage(voxel) <= 128) return false;
      const Point center(voxel[0], voxel[1], voxel[2]);
      // The native 3D distance uses an L3 metric, a lower bound on Euclidean
      // distance. Subtract the voxel half-diagonal and sample-gap bound.
      if (mask.myDistanceImage(voxel)-(p-center).norm() < required) return false;
    }
    return true;
  }
  bool operator()(Tree &tree) {
    tree.myLengthFactor = 1.;
    tree.updateResistanceFromRoot();
    tree.updateRootRadius();
    std::set<std::array<double, 3>> positions;
    for (const auto &segment : tree.myVectSegments) {
      const auto &p = segment.myCoordinate;
      if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2]) ||
          !std::isfinite(segment.myRadius) || segment.myRadius <= 0.) return reject("nonfinite");
      if (!positions.insert(std::array<double, 3>{{p[0], p[1], p[2]}}).second) return reject("duplicate");
    }
    for (std::size_t i = 1; i < tree.myVectSegments.size(); ++i) {
      const auto parent = tree.myVectParent[i];
      const auto &a = tree.myVectSegments[parent].myCoordinate;
      const auto &b = tree.myVectSegments[i].myCoordinate;
      const double radius = tree.myVectSegments[i].myRadius;
      if ((b-a).norm() <= 0.) return reject("degenerate");
      if (tree.myFixedRootRadius > 0. && 2*radius > (b-a).norm()) return reject("short_diameter");
      if (parent == 0 && (b-a).dot(outward) <= 0.) return reject("foldback");
      if (!wall(a, b, radius)) return reject("wall");
      for (const auto &obstacle : fixed) {
        if (parent == 0 && obstacle.system == system && obstacle.sourceEdge == incident) continue;
        if (boxesOverlap(a, b, radius, obstacle.a, obstacle.b, obstacle.radius) &&
            GeomHelpers::segment2segmentDistance<3>(a, b, obstacle.a, obstacle.b) <= radius+obstacle.radius) {
          const bool separated = FrustumGeometry::separated(a, b, {{radius, radius}},
                                  obstacle.a, obstacle.b, {{obstacle.radiusA, obstacle.radiusB}});
          if (collisionTraces.size() < 128) {
            std::ostringstream row;
            row << std::setprecision(17) << i << ',' << parent << ',' << obstacle.system << ','
                << obstacle.sourceEdge << ',' << separated;
            for (const auto &p : {a, b, obstacle.a, obstacle.b})
              for (unsigned int j = 0; j < 3; ++j) row << ',' << p[j];
            row << ',' << radius << ',' << obstacle.radiusA << ',' << obstacle.radiusB;
            collisionTraces.push_back(row.str());
          }
          if (!separated) return reject("fixed_collision");
          ++counts["frustum_certificates"];
        }
      }
      for (std::size_t j = 1; j < i; ++j) {
        const auto otherParent = tree.myVectParent[j];
        if (parent == otherParent || parent == j || otherParent == i) continue;
        const auto &c = tree.myVectSegments[otherParent].myCoordinate;
        const auto &d = tree.myVectSegments[j].myCoordinate;
        const double otherRadius = tree.myVectSegments[j].myRadius;
        if (boxesOverlap(a, b, radius, c, d, otherRadius) &&
            GeomHelpers::segment2segmentDistance<3>(a, b, c, d) <= radius+otherRadius &&
            !FrustumGeometry::separated(a, b, {{radius, radius}}, c, d, {{otherRadius, otherRadius}}))
          return reject("self_collision");
      }
    }
    ++counts["valid_candidates"];
    return true;
  }
};
}

int main(int argc, char **argv) {
  CLI::App app;
  std::string input, maskPath, output, allocationPath, surfacePath;
  std::vector<std::string> obstacleFiles, grownFiles;
  std::vector<double> originValues;
  double pitch = .002, flow = 0., pressureDrop = 0., viscosity = .00371;
  unsigned int target = 64, seed = 23, searches = 100, samples = 1000;
  bool matchRadius = false;
  app.add_option("--input", input)->required();
  app.add_option("--mask", maskPath)->required();
  app.add_option("--output", output)->required();
  app.add_option("--origin", originValues)->expected(3)->required();
  app.add_option("--pitch", pitch);
  app.add_option("--flow", flow)->required();
  app.add_option("--pressure-drop", pressureDrop)->required();
  app.add_option("--viscosity", viscosity);
  app.add_option("--allocation", allocationPath, "CSV node,terminals,flow_fraction for every Grow endpoint");
  app.add_option("--surface", surfacePath, "Original closed STL in metres (requires VTK build)");
  app.add_flag("--match-attachment-radius", matchRadius, "Preserve measured inlet lumen and derive required pressure");
  app.add_option("--terminals", target);
  app.add_option("--seed", seed);
  app.add_option("--searches", searches);
  app.add_option("--samples", samples);
  app.add_option("--obstacles", obstacleFiles)->expected(-1);
  app.add_option("--grown-obstacles", grownFiles)->expected(-1);
  CLI11_PARSE(app, argc, argv);
  try {
    if (!std::isfinite(pitch) || pitch <= 0. || !std::isfinite(flow) || flow <= 0. ||
        !std::isfinite(pressureDrop) || pressureDrop <= 0. || !std::isfinite(viscosity) || viscosity <= 0. ||
        target < 1 || target > 1000000 || searches < 1 || searches > 1000000 ||
        samples < 1 || samples > 1000000)
      throw std::runtime_error("Require finite positive physical parameters and budgets");
    const Point origin(originValues[0], originValues[1], originValues[2]);
    if (!std::isfinite(origin.norm())) throw std::runtime_error("Invalid physical origin");
    const auto source = readAnatomy(input);
    OpenCCO::GrownAnatomy network(source);
    if (std::filesystem::exists(output)) throw std::runtime_error("Output directory must be new");
    std::filesystem::create_directories(output);
    std::ofstream preserved(std::filesystem::path(output)/"source.anatomy");
    OpenCCO::writeAnatomicalTree(preserved, source);
    preserved.close();
    Mask mask(maskPath, 128, 100);
    mask.myMinDistanceToBorder = 1.;
    std::size_t foreground = 0;
    for (const auto &p : mask.myImage.domain()) foreground += mask.myImage(p) > 128;
    std::unique_ptr<OpenCCO::SurfaceWall> surface;
    decltype(Filter::surfaceCertificate) surfaceCertificate;
    if (!surfacePath.empty()) {
      surface = std::make_unique<OpenCCO::SurfaceWall>(surfacePath);
      surfaceCertificate = [&](const Point &a, const Point &b, double radius) {
        return surface->tubeInside(origin+a*pitch, origin+b*pitch, radius*pitch);
      };
    }
    std::vector<Obstacle> obstacles;
    if (std::find(obstacleFiles.begin(), obstacleFiles.end(), input) == obstacleFiles.end())
      obstacleFiles.push_back(input);
    std::map<std::string, std::string> roleFingerprints;
    std::ostringstream sourceFingerprint;
    OpenCCO::writeAnatomicalTree(sourceFingerprint, source);
    for (const auto &file : obstacleFiles) {
      const auto anatomy = readAnatomy(file);
      std::ostringstream fingerprint;
      OpenCCO::writeAnatomicalTree(fingerprint, anatomy);
      if (anatomy.role == source.role && fingerprint.str() != sourceFingerprint.str())
        throw std::runtime_error("Growth source and same-role obstacle anatomy disagree");
      const auto found = roleFingerprints.find(anatomy.role);
      if (found != roleFingerprints.end()) {
        if (found->second != fingerprint.str()) throw std::runtime_error("Conflicting same-role obstacles");
        continue;
      }
      roleFingerprints[anatomy.role] = fingerprint.str();
      for (std::size_t i = 0; i < anatomy.edges.size(); ++i) {
        const auto &edge = anatomy.edges[i];
        obstacles.push_back({(point(anatomy.nodes[edge[0]].position)-origin)/pitch,
                             (point(anatomy.nodes[edge[1]].position)-origin)/pitch,
                             std::max(anatomy.nodes[edge[0]].radius, anatomy.nodes[edge[1]].radius)/pitch,
                             systemId(anatomy.role), static_cast<long long>(i),
                             anatomy.nodes[edge[0]].radius/pitch, anatomy.nodes[edge[1]].radius/pitch});
      }
    }
    for (const auto &file : grownFiles) {
      std::ifstream stream(file);
      if (!stream) throw std::runtime_error("Cannot read grown obstacles");
      std::string line;
      std::getline(stream, line);
      while (std::getline(stream, line)) {
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream row(line);
        int system, generated;
        std::size_t id, a, b;
        long long sourceEdge;
        Point x, y;
        double r0, r1;
        if (!(row >> system >> id >> a >> b >> x[0] >> x[1] >> x[2] >> y[0] >> y[1] >> y[2]
                  >> r0 >> r1 >> generated >> sourceEdge)) throw std::runtime_error("Invalid grown obstacle record");
        if (system < 0 || system > 2 || (generated != 0 && generated != 1) ||
            !std::isfinite(x.norm()) || !std::isfinite(y.norm()) || (x-y).norm() <= 0. ||
            !std::isfinite(r0) || !std::isfinite(r1) || r0 <= 0. || r1 <= 0.)
          throw std::runtime_error("Nonfinite or invalid grown obstacle geometry");
        if (generated) obstacles.push_back({(x-origin)/pitch, (y-origin)/pitch,
                                            std::max(r0, r1)/pitch, system, -1, r0/pitch, r1/pitch});
      }
    }
    std::vector<std::size_t> attachments;
    for (std::size_t i = 0; i < source.nodes.size(); ++i) if (source.nodes[i].grow) attachments.push_back(i);
    if (attachments.empty() || target < attachments.size()) throw std::runtime_error("Insufficient target/Grow endpoints");
    std::map<std::size_t, std::pair<unsigned int, double>> allocation;
    if (!allocationPath.empty()) {
      std::ifstream stream(allocationPath);
      if (!stream) throw std::runtime_error("Cannot read source allocation");
      std::string line;
      std::getline(stream, line);
      std::size_t total = 0;
      double totalFraction = 0.;
      while (std::getline(stream, line)) {
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream row(line);
        std::size_t node;
        unsigned int quota;
        double fraction;
        std::string extra;
        if (!(row >> node >> quota >> fraction) || (row >> extra) || quota < 1 || quota > target ||
            !std::isfinite(fraction) || fraction <= 0. || fraction > 1. ||
            std::find(attachments.begin(), attachments.end(), node) == attachments.end() ||
            !allocation.emplace(node, std::make_pair(quota, fraction)).second)
          throw std::runtime_error("Invalid source allocation row");
        total += quota;
        totalFraction += fraction;
      }
      if (allocation.size() != attachments.size() || total != target ||
          std::abs(totalFraction-1.) > 128.*std::numeric_limits<double>::epsilon()*attachments.size())
        throw std::runtime_error("Source allocation must cover every Grow node, total target and unit flow fraction");
    }
    srand(seed);
    bool complete = true;
    std::size_t processed = 0;
    for (std::size_t part = 0; part < attachments.size(); ++part) {
      const auto attachment = attachments[part];
      std::size_t parent = source.nodes.size(), incident = source.edges.size();
      for (std::size_t i = 0; i < source.edges.size(); ++i) {
        if (source.edges[i][0] == attachment) { parent = source.edges[i][1]; incident = i; }
        if (source.edges[i][1] == attachment) { parent = source.edges[i][0]; incident = i; }
      }
      const Point root = (point(source.nodes[attachment].position)-origin)/pitch;
      Point direction = root-(point(source.nodes[parent].position)-origin)/pitch;
      direction /= direction.norm();
      const unsigned int quota = allocation.empty() ?
          target/attachments.size() + (part < target%attachments.size()) : allocation.at(attachment).first;
      const double fraction = allocation.empty() ? double(quota)/target : allocation.at(attachment).second;
      OpenCCO::SubtreeParameters parameters{quota, double(foreground)*quota/target,
          flow/(pitch*pitch*pitch)*fraction, pressureDrop, viscosity};
      Filter filter{mask, obstacles, systemId(source.role), static_cast<long long>(incident), direction, {}, {}, surfaceCertificate};
      std::unique_ptr<Tree> tree;
      const auto nativeDiagnostics = std::make_shared<std::map<std::string, std::size_t>>();
      const auto makeSeed = [&](const Point &terminal) {
        auto candidate = OpenCCO::seedAnatomicalSubtree(mask, root, terminal, parameters);
        candidate.myCandidateDiagnostics = nativeDiagnostics;
        if (matchRadius) candidate.myFixedRootRadius = source.nodes[attachment].radius/pitch;
        return candidate;
      };
      // Preserve a valid tangent seed before searching other outward directions.
      for (const double length : {5., 2.5, 1.25, .625}) {
        auto candidate = makeSeed(root+direction*length);
        if (filter(candidate)) { tree = std::make_unique<Tree>(candidate); break; }
      }
      Mask::TPointI center;
      for (unsigned int j = 0; j < 3; ++j) center[j] = static_cast<int>(std::lround(root[j]));
      for (const double length : {5., 2.5, 1.25, .625}) {
        if (tree) break;
        // Reuse the native mask initializer's digital-sphere search. Only the
        // generated endpoint is sampled; the anatomical root remains exact.
        const auto endpoints = GeomHelpers::pointsOnSphere<Mask::TPointI, Mask::TDGset>(center, length);
        double bestClearance = -std::numeric_limits<double>::infinity();
        for (const auto &endpoint : endpoints) {
          const Point terminal(endpoint[0], endpoint[1], endpoint[2]);
          if (terminal == root || (terminal-root).dot(direction) <= 0.) continue;
          auto candidate = makeSeed(terminal);
          if (!filter(candidate)) continue;
          const double clearance = mask.myDistanceImage(endpoint);
          if (clearance > bestClearance) {
            bestClearance = clearance;
            tree = std::make_unique<Tree>(candidate);
          }
        }
        if (tree) break;
      }
      if (!tree) {
        complete = false;
        std::cerr << "Cannot initialize anatomical attachment " << attachment << '\n';
      } else {
        for (unsigned int chunk = std::min(16u, quota);; chunk = std::min(chunk+16, quota)) {
          tree->my_NTerm = chunk;
          const bool grown = ExpandTreeHelpers::expandTree(*tree, false, searches, samples,
                                                          [&](Tree &candidate) { return filter(candidate); });
          auto snapshot = network;
          snapshot.append(attachment, *tree, origin, pitch);
          writeNetwork(std::filesystem::path(output)/("checkpoint-"+std::to_string(part)+"-"+
                                                     std::to_string(tree->myKTerm)), snapshot);
          if (!grown || chunk == quota) { complete = complete && grown; break; }
        }
        network.append(attachment, *tree, origin, pitch);
        for (std::size_t i = 1; i < tree->myVectSegments.size(); ++i) {
          const auto &segment = tree->myVectSegments[i];
          obstacles.push_back({tree->myVectSegments[tree->myVectParent[i]].myCoordinate,
                               segment.myCoordinate, segment.myRadius, systemId(source.role), -1,
                               segment.myRadius, segment.myRadius});
        }
        ++processed;
      }
      for (const auto &entry : *nativeDiagnostics) filter.counts["native_"+entry.first] += entry.second;
      std::ofstream diagnostic(std::filesystem::path(output)/("attachment-"+std::to_string(attachment)+".json"));
      diagnostic << "{\"attachment\":" << attachment << ",\"requested\":" << quota
                 << ",\"actual\":" << (tree ? tree->myKTerm : 0) << ",\"design_flow_fraction\":"
                 << std::setprecision(17) << fraction
                 << ",\"required_pressure_drop_pa\":" << (tree ? tree->my_pDrop : -1.)
                 << ",\"inlet_radius_m\":" << (tree ? tree->myVectSegments[1].myRadius*pitch : -1.)
                 << ",\"candidate_counts\":{";
      bool comma = false;
      for (const auto &entry : filter.counts) {
        if (comma) diagnostic << ',';
        diagnostic << '"' << entry.first << "\":" << entry.second;
        comma = true;
      }
      diagnostic << "}}\n";
      std::ofstream traces(std::filesystem::path(output)/("collision-"+std::to_string(attachment)+".csv"));
      traces << "candidate,parent,obstacle_system,source_edge,separated,a0,a1,a2,b0,b1,b2,c0,c1,c2,d0,d1,d2,radius,source_r0,source_r1\n";
      for (const auto &row : filter.collisionTraces) traces << row << '\n';
      if (!complete) break;
    }
    writeNetwork(output, network);
    complete = complete && network.tissueTerminals().size() == target;
    std::ofstream summary(std::filesystem::path(output)/"summary.json");
    summary << std::setprecision(17) << "{\"role\":\"" << source.role << "\",\"requested\":" << target
            << ",\"actual\":" << network.tissueTerminals().size() << ",\"role_growth_complete\":"
            << (complete ? "true" : "false") << ",\"processed_attachments\":" << processed
            << ",\"imported_edges\":" << source.edges.size() << ",\"flow_m3_s\":" << flow
            << ",\"local_pressure_drop_pa\":" << pressureDrop << ",\"viscosity_pa_s\":" << viscosity
            << ",\"pressure_drop_prescribed\":" << (matchRadius ? "false" : "true")
            << ",\"wall_model\":\"" << (surface ? "original-stl-lipschitz" : "mask-l3-lower-bound")
            << "\",\"physical_length_factor\":1,\"coupled_ready\":false}\n";
    std::cout << source.role << ": " << network.tissueTerminals().size() << '/' << target << " tissue tips\n";
    return complete ? 0 : 1;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
