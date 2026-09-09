#include "CageGenerator.hh"
#include "CageSimplifier/SimplifyStages/FlipStage.hh"
#include "CageSimplifier/SimplifyStages/RelocateStage.hh"
#include "CageSimplifier/SimplifyStages/CollapseStage.hh"
#include "CageSimplifier/Topo/EdgeCollapser.h"
#include "CageSimplifier/Topo/VertexRelocater.h"
#include "CageSimplifier/CageSimplifier.hh"
#include "Geometry/Exact/TriTriIntersect.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <map>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>
#include "spdlog/sinks/ostream_sink.h"
#include "omp.h"

namespace
{
using namespace Cage;
using namespace Cage::CageSimp;
using Cage::Geometry::Vec3d;
using Cage::Geometry::ExactPoint;
using Cage::SurfaceMesh::SMeshT;
using OpenMesh::VertexHandle;
using OpenMesh::EdgeHandle;
using OpenMesh::FaceHandle;

void require(bool condition, const char* message)
{
  if (!condition)
    throw std::runtime_error(message);
}

bool close(const Vec3d& a, const Vec3d& b, double tolerance = 1e-11)
{
  return (a - b).length() <= tolerance;
}

void add_face(SMeshT& mesh, VertexHandle a, VertexHandle b, VertexHandle c)
{
  require(mesh.add_face(a, b, c).is_valid(), "fixture must form a manifold triangle mesh");
}

void add_triangle(SMeshT& mesh, const Vec3d& a, const Vec3d& b, const Vec3d& c)
{
  const auto va = mesh.add_vertex(a);
  const auto vb = mesh.add_vertex(b);
  const auto vc = mesh.add_vertex(c);
  add_face(mesh, va, vb, vc);
}

double triangle_quality(const Vec3d& a, const Vec3d& b, const Vec3d& c)
{
  const double denominator = (b - a).sqrnorm() + (c - b).sqrnorm() + (a - c).sqrnorm();
  return denominator > 0.0 ?
    2.0 * std::sqrt(3.0) * (b - a).cross(c - a).length() / denominator : 0.0;
}

double min_quality(const SMeshT& mesh)
{
  double result = 1.0;
  for (FaceHandle face : mesh.faces())
  {
    std::array<Vec3d, 3> points;
    size_t index = 0;
    for (VertexHandle vertex : mesh.fv_range(face))
      points[index++] = mesh.point(vertex);
    result = std::min(result, triangle_quality(points[0], points[1], points[2]));
  }
  return result;
}

struct Fixture
{
  SMeshT original;
  SMeshT cage;
  ParamCageGenerator parameters;
  std::unique_ptr<VertexTree> vertices;
  std::unique_ptr<DFaceTree> source_tree;
  std::unique_ptr<LightDFaceTree> cage_tree;
  std::unique_ptr<FaceGrid> source_grid;

  Fixture()
  {
    // A real, nonempty source with extent in every axis keeps all production
    // spatial structures active. It is remote from the small test cages.
    const auto a = original.add_vertex(Vec3d(10.0, 10.0, 10.0));
    const auto b = original.add_vertex(Vec3d(11.0, 10.0, 10.0));
    const auto c = original.add_vertex(Vec3d(10.0, 11.0, 10.0));
    const auto d = original.add_vertex(Vec3d(10.0, 10.0, 11.0));
    add_face(original, a, c, b);
    add_face(original, a, b, d);
    add_face(original, a, d, c);
    add_face(original, b, c, d);
  }

  void initialize()
  {
    pre_calculate_edge_length(&original);
    pre_calculate_face_area(&original);
    pre_calculate_edge_length(&cage);
    pre_calculate_face_area(&cage);
    init_one_ring_faces(&cage);
    cage.update_normals();
    vertices = std::make_unique<VertexTree>(original);
    source_tree = std::make_unique<DFaceTree>(original);
    source_grid = std::make_unique<FaceGrid>(original);
    cage_tree = std::make_unique<LightDFaceTree>(cage);
  }

  FlipStage flips()
  {
    return FlipStage(&original, &cage, &parameters.paramCageSimplifier.paramFlip,
      source_tree.get(), cage_tree.get(), source_grid.get(), 1.0);
  }

  RelocateStage relocations()
  {
    return RelocateStage(&original, &cage, &parameters.paramCageSimplifier.paramRelocate,
      vertices.get(), source_tree.get(), cage_tree.get(), source_grid.get(), 1.0);
  }

  VertexRelocater relocater(VertexHandle center)
  {
    VertexRelocater result(&original, &cage, source_tree.get(), cage_tree.get(),
      vertices.get(), source_grid.get());
    require(result.init(center), "relocation fixture center must initialize");
    return result;
  }

  void require_source_clear()
  {
    for (FaceHandle face : cage.faces())
    {
      std::array<VertexHandle, 3> corners;
      size_t index = 0;
      for (VertexHandle vertex : cage.fv_range(face))
        corners[index++] = vertex;
      require(!source_tree->do_intersect(cage.point(corners[0]), cage.point(corners[1]),
        cage.point(corners[2]), cage.data(corners[0]).ep.get(), cage.data(corners[1]).ep.get(),
        cage.data(corners[2]).ep.get()), "fixture or result unexpectedly intersects source");
    }
  }
};

std::array<VertexHandle, 4> add_flip_quad(SMeshT& mesh)
{
  // BD is the poorer diagonal (minimum quality about .632); AC raises it
  // to about .797. The nonplanarity makes collision rejection observable.
  std::array<VertexHandle, 4> v = {
    mesh.add_vertex(Vec3d(0.0, 0.0, 0.0)),
    mesh.add_vertex(Vec3d(2.0, 0.0, 0.0)),
    mesh.add_vertex(Vec3d(1.5, 1.0, 0.2)),
    mesh.add_vertex(Vec3d(0.0, 1.0, 0.0))
  };
  add_face(mesh, v[0], v[1], v[3]);
  add_face(mesh, v[1], v[2], v[3]);
  return v;
}

EdgeHandle diagonal(SMeshT& mesh, const std::array<VertexHandle, 4>& v)
{
  return mesh.edge_handle(mesh.find_halfedge(v[1], v[3]));
}

void improving_flip()
{
  Fixture fixture;
  const auto v = add_flip_quad(fixture.cage);
  fixture.initialize();
  fixture.require_source_clear();
  const double before = min_quality(fixture.cage);
  std::vector<Vec3d> before_points;
  for (VertexHandle vertex : fixture.cage.vertices())
    before_points.push_back(fixture.cage.point(vertex));
  // The new pass must choose quality independently of legacy valence mode.
  fixture.parameters.paramCageSimplifier.paramFlip.requireRegularValence = true;
  require(fixture.flips().do_quality_flip() == 1, "quality-improving diagonal must flip once");
  require(min_quality(fixture.cage) > before + 0.1, "flip must improve worst triangle quality");
  require(fixture.cage.n_vertices() == before_points.size(), "flip must preserve vertex count");
  require(fixture.cage.find_halfedge(v[0], v[2]).is_valid(), "new diagonal AC must exist");
  require(!fixture.cage.find_halfedge(v[1], v[3]).is_valid(), "old diagonal BD must disappear");
  for (VertexHandle vertex : fixture.cage.vertices())
    require(close(fixture.cage.point(vertex), before_points[vertex.idx()]), "flip must preserve positions");
  require(fixture.flips().do_quality_flip() == 0, "quality-optimal diagonal must remain stable");
  fixture.require_source_clear();
}

void rail_flip_rejected()
{
  Fixture fixture;
  const auto v = add_flip_quad(fixture.cage);
  const auto edge = diagonal(fixture.cage, v);
  fixture.cage.data(edge).boundary_rail_id = 7;
  fixture.cage.data(v[1]).boundary_rail_id = 7;
  fixture.cage.data(v[3]).boundary_rail_id = 7;
  fixture.initialize();
  require(fixture.flips().do_quality_flip() == 0, "rail diagonal must never flip");
  require(fixture.cage.find_halfedge(v[1], v[3]).is_valid(), "rejected rail must retain connectivity");
  require(fixture.cage.data(edge).boundary_rail_id == 7, "rail label must survive rejection");
}

void overlapping_flip_rejected()
{
  Fixture fixture;
  const auto a = fixture.cage.add_vertex(Vec3d(0.0, 0.0, 0.0));
  const auto b = fixture.cage.add_vertex(Vec3d(2.0, 0.0, 0.0));
  const auto c = fixture.cage.add_vertex(Vec3d(0.1, 0.1, 0.0));
  const auto d = fixture.cage.add_vertex(Vec3d(0.0, 2.0, 0.0));
  add_face(fixture.cage, a, b, c);
  add_face(fixture.cage, a, c, d);
  fixture.initialize();
  const auto edge = fixture.cage.edge_handle(fixture.cage.find_halfedge(a, c));
  require(fixture.cage.is_flip_ok(edge), "concave quad must be topologically flippable");
  const double proposed_quality = std::min(
    triangle_quality(fixture.cage.point(a), fixture.cage.point(b), fixture.cage.point(d)),
    triangle_quality(fixture.cage.point(b), fixture.cage.point(c), fixture.cage.point(d)));
  require(proposed_quality > min_quality(fixture.cage) + 0.5,
    "overlapping replacement must look attractive to the quality metric");
  // C lies strictly inside ABD. The replacement faces ABD and BCD would
  // overlap even though they improve both triangle shapes individually.
  require(fixture.flips().do_quality_flip() == 0, "concave quad flip must reject overlapping new faces");
  require(fixture.cage.find_halfedge(a, c).is_valid(), "overlap rejection must keep original diagonal");
}

void intersecting_flip_rejected(bool source_obstacle)
{
  Fixture fixture;
  const auto v = add_flip_quad(fixture.cage);
  const Vec3d p(3.5 / 3.0, 1.0 / 3.0, 0.2 / 3.0);
  SMeshT& obstacle_mesh = source_obstacle ? fixture.original : fixture.cage;
  add_triangle(obstacle_mesh, p + Vec3d(0.0, 0.0, -0.005),
    p + Vec3d(0.0, 0.0, 0.005), p + Vec3d(0.005, 0.0, 0.0));
  fixture.initialize();
  fixture.require_source_clear();
  EdgeFlipper flipper(&fixture.original, &fixture.cage, fixture.source_tree.get(),
    fixture.cage_tree.get(), fixture.source_grid.get());
  require(flipper.init(diagonal(fixture.cage, v)), "colliding diagonal must be topologically flippable");
  require(flipper.local_Hausdorff_after_flipping() == DBL_MAX,
    "fixture must intersect only after proposed flip");
  const size_t before_vertices = fixture.cage.n_vertices();
  require(fixture.flips().do_quality_flip() == 0,
    "quality-improving flip must reject source or self-intersection");
  require(fixture.cage.find_halfedge(v[1], v[3]).is_valid(), "rejected flip must retain old diagonal");
  require(fixture.cage.n_vertices() == before_vertices, "rejected flip must preserve vertex count");
}

std::array<VertexHandle, 5> add_fan(SMeshT& mesh, const Vec3d& center)
{
  std::array<VertexHandle, 5> v = {
    mesh.add_vertex(center),
    mesh.add_vertex(Vec3d(1.0, 0.0, 0.0)),
    mesh.add_vertex(Vec3d(0.0, 1.0, 0.0)),
    mesh.add_vertex(Vec3d(-1.0, 0.0, 0.0)),
    mesh.add_vertex(Vec3d(0.0, -1.0, 0.0))
  };
  for (size_t i = 1; i <= 4; ++i)
    add_face(mesh, v[0], v[i], v[i == 4 ? 1 : i + 1]);
  // Keep a labeled, closed rail around the movable center.
  for (size_t i = 1; i <= 4; ++i)
  {
    mesh.data(v[i]).boundary_rail_id = 3;
    mesh.data(mesh.edge_handle(mesh.find_halfedge(v[i], v[i == 4 ? 1 : i + 1]))).boundary_rail_id = 3;
  }
  return v;
}

void voronoi_target_uses_neighbor_faces()
{
  Fixture fixture;
  const auto v = add_fan(fixture.cage, Vec3d(0.0, 0.0, 0.0));
  const auto outer = fixture.cage.add_vertex(Vec3d(1.0, 1.0, 0.0));
  add_face(fixture.cage, v[2], v[1], outer);
  fixture.initialize();
  auto relocater = fixture.relocater(v[0]);
  Vec3d target;
  require(relocater.find_area_equalizing_tangential_smooth_target(target), "Voronoi target must exist");
  // Each diamond neighbor initially owns 1/4 of area. The exterior right
  // triangle adds 1/8 to east and north. Their total area is now 5/4,
  // giving exactly (1/10, 1/10, 0). The center fan alone is symmetric.
  require(close(target, Vec3d(0.1, 0.1, 0.0)), "target must include neighbor faces outside center fan");
  require(close(relocater.find_tangential_smooth_target(), Vec3d(0.0, 0.0, 0.0)),
    "fixture must distinguish Voronoi weighting from uniform averaging");
  require(close(relocater.find_weighted_tangential_smooth_target(), Vec3d(0.0, 0.0, 0.0)),
    "fixture must distinguish Voronoi weighting from face-centroid averaging");
}

void improving_relocation()
{
  Fixture fixture;
  fixture.parameters.paramCageSimplifier.paramRelocate.surfaceWeight = 0.0;
  const auto v = add_fan(fixture.cage, Vec3d(0.6, 0.0, 0.4));
  fixture.cage.data(v[0]).ep = std::make_unique<ExactPoint>(fixture.cage.point(v[0]));
  fixture.initialize();
  const double before = min_quality(fixture.cage);
  const Vec3d old_center = fixture.cage.point(v[0]);
  std::array<Vec3d, 4> rail_points;
  for (size_t i = 1; i <= 4; ++i)
    rail_points[i - 1] = fixture.cage.point(v[i]);
  require(fixture.relocations().do_quality_relocate(8) == 1, "only the free center should relocate");
  require(min_quality(fixture.cage) > before + 1e-6, "relocation must improve minimum fan quality");
  require(!close(fixture.cage.point(v[0]), old_center), "successful relocation must change center");
  require(std::abs(fixture.cage.point(v[0])[2] - old_center[2]) < 1e-12,
    "relocation direction must remain in center tangent plane");
  require(!fixture.cage.data(v[0]).ep, "relocation must clear stale exact coordinates");
  require(fixture.cage.n_vertices() == 5 && fixture.cage.n_faces() == 4,
    "relocation must preserve mesh topology and vertex count");
  for (size_t i = 1; i <= 4; ++i)
  {
    require(close(fixture.cage.point(v[i]), rail_points[i - 1]), "rail vertices must remain fixed");
    require(fixture.cage.data(v[i]).boundary_rail_id == 3, "rail vertex labels must survive relocation");
  }
  for (EdgeHandle edge : fixture.cage.edges())
    require(std::abs(fixture.cage.data(edge).edge_length - fixture.cage.calc_edge_length(edge)) < 1e-12,
      "relocation must refresh edge-length cache");
  fixture.require_source_clear();
}

void symmetric_relocation_is_noop()
{
  Fixture fixture;
  fixture.parameters.paramCageSimplifier.paramRelocate.surfaceWeight = 0.0;
  const auto v = add_fan(fixture.cage, Vec3d(0.0, 0.0, 0.0));
  fixture.initialize();
  require(fixture.relocations().do_quality_relocate(8) == 0, "symmetric fan has no improving smoothing step");
  require(close(fixture.cage.point(v[0]), Vec3d(0.0, 0.0, 0.0)), "failed relocation must retain center");
}

void tiny_surface_weight_preserves_tangential_descent()
{
  const auto relocated_center = [](double surface_weight)
  {
    Fixture fixture;
    fixture.parameters.paramCageSimplifier.paramRelocate.surfaceWeight = surface_weight;
    for (VertexHandle vertex : fixture.original.vertices())
      fixture.original.set_point(vertex,
        fixture.original.point(vertex) + Vec3d(1e6, 1e6, 1e6));
    const Vec3d old_point(0.2, 0.0, 0.4);
    const auto v = add_fan(fixture.cage, old_point);
    fixture.initialize();
    require(fixture.relocations().do_quality_relocate(1) == 1,
      "a distant surface with negligible weight must not block meaningful tangential descent");
    require((fixture.cage.point(v[0]) - old_point).length() > 0.01,
      "tiny-weight regression must perform a substantial tangential move");
    fixture.require_source_clear();
    return fixture.cage.point(v[0]);
  };
  const Vec3d pure_tangent = relocated_center(0.0);
  const Vec3d tiny_attraction = relocated_center(1e-20);
  require(close(pure_tangent, tiny_attraction, 1e-10),
    "vanishing surface weight must approach the pure tangential solution continuously");
}

void invalid_relocation_is_rejected()
{
  Fixture fixture;
  const auto v = add_fan(fixture.cage, Vec3d(0.2, 0.0, 0.0));
  fixture.initialize();
  auto relocater = fixture.relocater(v[0]);
  const Vec3d before = fixture.cage.point(v[0]);
  require(!relocater.try_relocate_vertex(Vec3d(std::numeric_limits<double>::infinity(), 0.0, 0.0)),
    "infinite relocation target must reject before exact predicates");
  require(!relocater.try_relocate_vertex(Vec3d(0.0, std::numeric_limits<double>::quiet_NaN(), 0.0)),
    "NaN relocation target must reject before exact predicates");
  require(!relocater.try_relocate_vertex(fixture.cage.point(v[1])),
    "relocation that collapses fan triangles must reject");
  require(close(fixture.cage.point(v[0]), before), "invalid relocation must preserve original position");
}

void collision_backtracking()
{
  Fixture fixture;
  fixture.parameters.paramCageSimplifier.paramRelocate.surfaceWeight = 0.0;
  const auto v = add_fan(fixture.cage, Vec3d(0.6, 0.0, 0.4));
  fixture.initialize();
  const Vec3d old_center = fixture.cage.point(v[0]);
  Vec3d target;
  require(fixture.relocater(v[0]).find_area_equalizing_tangential_smooth_target(target),
    "backtracking fixture requires a smoothing target");
  const Vec3d half_target = (old_center + target) * 0.5;
  const Vec3d west = fixture.cage.point(v[3]);
  const Vec3d south = fixture.cage.point(v[4]);
  // Put a tiny vertical obstacle through the full-target face interior.
  // Its distance from the half-target plane determines a safe fixture size.
  const Vec3d point = (target + west + south) / 3.0;
  const Vec3d half_normal = (west - half_target).cross(south - half_target);
  const double separation = std::abs((half_normal | (point - half_target)) / half_normal.length());
  require(separation > 1e-4, "full and half-target surfaces must be distinct");
  const double radius = separation * 0.1;
  add_triangle(fixture.original, point + Vec3d(0.0, 0.0, -radius),
    point + Vec3d(0.0, 0.0, radius), point + Vec3d(radius, 0.0, 0.0));
  fixture.initialize();
  fixture.require_source_clear();
  auto relocater = fixture.relocater(v[0]);
  require(relocater.relocate_would_cause_intersection(target), "full step must encounter source collision");
  require(!relocater.relocate_would_cause_intersection(half_target), "half step must clear source collision");
  const double before = min_quality(fixture.cage);
  require(fixture.relocations().do_quality_relocate(1) == 0,
    "one-candidate search must reject colliding full step");
  require(close(fixture.cage.point(v[0]), old_center), "exhausted search must not mutate mesh");
  require(fixture.relocations().do_quality_relocate(8) == 1,
    "backtracking must recover with a collision-free smaller step");
  require(close(fixture.cage.point(v[0]), half_target), "first valid half step should be accepted");
  require(min_quality(fixture.cage) > before + 1e-6, "recovered step must still improve triangle quality");
  fixture.require_source_clear();
}

std::array<VertexHandle, 5> add_attraction_fan(
  Fixture& fixture, double center_height = 1.0, double source_height = 0.4)
{
  const auto v = add_fan(fixture.cage, Vec3d(0.0, 0.0, center_height));
  // A small source face directly below (or above) the center has an
  // analytical nearest point, without intersecting the sloped cage fan.
  add_triangle(fixture.original, Vec3d(-0.1, -0.1, source_height),
    Vec3d(0.1, -0.1, source_height), Vec3d(0.0, 0.1, source_height));
  fixture.initialize();
  fixture.require_source_clear();
  return v;
}

void surface_attraction_decreases_frozen_energy()
{
  Fixture fixture;
  const auto v = add_attraction_fan(fixture);
  const Vec3d old_point = fixture.cage.point(v[0]);
  Vec3d tangent_target;
  require(fixture.relocater(v[0]).find_area_equalizing_tangential_smooth_target(tangent_target),
    "surface-attraction fixture requires a tangential target");
  require(close(tangent_target, old_point), "symmetric fan must have no tangential displacement");
  const Vec3d surface_target = fixture.source_tree->closest_point(old_point).first;
  require(close(surface_target, Vec3d(0.0, 0.0, 0.4)),
    "closest point must lie inside the small source face, not at a source vertex");
  const double old_quality = min_quality(fixture.cage);
  const auto energy = [&](const Vec3d& point)
  {
    return (point - tangent_target).sqrnorm() + (point - surface_target).sqrnorm();
  };
  require(fixture.relocations().do_quality_relocate(1) == 1,
    "equal weights must accept the blended full step without any backtracking");
  const Vec3d new_point = fixture.cage.point(v[0]);
  require(close(new_point, Vec3d(0.0, 0.0, 0.7)),
    "equal-weight unconstrained solution must be the midpoint of the two fixed targets");
  require(energy(new_point) < energy(old_point) - 1e-12,
    "accepted relocation must decrease the quadratic energy with frozen targets");
  require((new_point - surface_target).sqrnorm() < (old_point - surface_target).sqrnorm(),
    "source attraction must reduce actual source distance through normal movement");
  require(min_quality(fixture.cage) < old_quality - 1e-6 &&
    min_quality(fixture.cage) >= fixture.parameters.paramCageSimplifier.paramRelocate.minTriangleQuality,
    "energy descent must permit a controlled triangle-quality decrease above the floor");
  for (size_t i = 1; i <= 4; ++i)
    require(fixture.cage.point(v[i])[2] == 0.0 && fixture.cage.data(v[i]).boundary_rail_id == 3,
      "surface attraction must leave rail vertices fixed");
  fixture.require_source_clear();
}

void relocation_weights_and_source_collision()
{
  Fixture fixture;
  const auto v = add_attraction_fan(fixture);
  auto& parameters = fixture.parameters.paramCageSimplifier.paramRelocate;
  const Vec3d old_point = fixture.cage.point(v[0]);
  parameters.surfaceWeight = 0.0;
  require(fixture.relocations().do_quality_relocate(12) == 0 &&
    close(fixture.cage.point(v[0]), old_point),
    "zero surface weight must preserve a fan with no tangential displacement");
  parameters.tangentialWeight = 0.0;
  require(fixture.relocations().do_quality_relocate(12) == 0 &&
    close(fixture.cage.point(v[0]), old_point),
    "both zero weights must disable relocation");
  parameters.surfaceWeight = 1.0;
  const Vec3d source_target(0.0, 0.0, 0.4);
  require(fixture.relocater(v[0]).relocate_would_cause_intersection(source_target),
    "distance-only full step must touch the source and fail the collision guard");
  require(fixture.relocations().do_quality_relocate(1) == 0 &&
    close(fixture.cage.point(v[0]), old_point),
    "a colliding distance-only step must not mutate the mesh when backtracking is exhausted");
  require(fixture.relocations().do_quality_relocate(12) == 1 &&
    close(fixture.cage.point(v[0]), Vec3d(0.0, 0.0, 0.7)),
    "zero tangential weight must allow distance-only attraction with a safe half step");
  fixture.require_source_clear();
}

void relocation_quality_floor_backtracks()
{
  Fixture fixture;
  const auto v = add_attraction_fan(fixture);
  auto& parameters = fixture.parameters.paramCageSimplifier.paramRelocate;
  parameters.minTriangleQuality = 0.99;
  const Vec3d old_point = fixture.cage.point(v[0]);
  const Vec3d full_target(0.0, 0.0, 0.7);
  const Vec3d half_target(0.0, 0.0, 0.85);
  require(triangle_quality(full_target, fixture.cage.point(v[1]), fixture.cage.point(v[2])) <
    parameters.minTriangleQuality &&
    triangle_quality(half_target, fixture.cage.point(v[1]), fixture.cage.point(v[2])) >
    parameters.minTriangleQuality,
    "quality-floor fixture must reject only the full step before the half step");
  require(!fixture.relocater(v[0]).relocate_would_cause_intersection(full_target),
    "quality-floor full step must be collision-free so this test isolates quality rejection");
  require(fixture.relocations().do_quality_relocate(1) == 0 &&
    close(fixture.cage.point(v[0]), old_point),
    "full-step quality-floor rejection must retain the original vertex");
  require(fixture.relocations().do_quality_relocate(12) == 1 &&
    close(fixture.cage.point(v[0]), half_target),
    "quality-floor rejection must backtrack to the first sufficiently good candidate");
  require(min_quality(fixture.cage) >= parameters.minTriangleQuality,
    "accepted backtracked relocation must satisfy the quality floor");
  fixture.require_source_clear();
}

void relocation_preserves_quality_below_floor()
{
  Fixture fixture;
  const auto v = add_attraction_fan(fixture, 20.0, 30.0);
  const Vec3d old_point = fixture.cage.point(v[0]);
  const double old_quality = min_quality(fixture.cage);
  require(old_quality < fixture.parameters.paramCageSimplifier.paramRelocate.minTriangleQuality,
    "tall fan must already have quality below the default floor");
  require(close(fixture.source_tree->closest_point(old_point).first, Vec3d(0.0, 0.0, 30.0)),
    "low-quality fixture must attract the center farther from the boundary plane");
  require(fixture.relocations().do_quality_relocate(12) == 0 &&
    close(fixture.cage.point(v[0]), old_point) && min_quality(fixture.cage) == old_quality,
    "energy descent must not worsen a fan that already lies below the quality floor");
}

void quality_configuration_roundtrip()
{
  ParamCageGenerator defaults;
  require(defaults.paramCageSimplifier.phase2QualityPolishIterations == 30,
    "linear simplification must default to at most thirty outer cycles");
  require(defaults.paramCageSimplifier.paramRelocate.qualitySweeps == 20,
    "relocation must default to at most twenty sweeps per outer cycle");
  require(defaults.paramCageSimplifier.paramRelocate.lineSearchMaxIter == 12,
    "relocation must default to twelve backtracking candidates");
  require(defaults.paramCageSimplifier.paramRelocate.tangentialWeight == 1.0 &&
    defaults.paramCageSimplifier.paramRelocate.surfaceWeight == 1.0 &&
    defaults.paramCageSimplifier.paramRelocate.minTriangleQuality == 0.2,
    "relocation must default to equal energy weights and a 0.2 quality floor");
  auto& configured = defaults.paramCageSimplifier;
  configured.phase2QualityPolishIterations = 5;
  configured.paramRelocate.qualitySweeps = 7;
  configured.paramRelocate.lineSearchMaxIter = 9;
  configured.paramRelocate.tangentialWeight = 2.5;
  configured.paramRelocate.surfaceWeight = 0.75;
  configured.paramRelocate.minTriangleQuality = 0.35;
  // Exercise the actual textual JSON boundary, including JSON integer types.
  auto serialized = boost::json::parse(boost::json::serialize(configured.serialize())).as_object();
  ParamCageGenerator restored;
  restored.paramCageSimplifier.deserialize(serialized);
  require(restored.paramCageSimplifier.phase2QualityPolishIterations == 5 &&
    restored.paramCageSimplifier.paramRelocate.qualitySweeps == 7 &&
    restored.paramCageSimplifier.paramRelocate.lineSearchMaxIter == 9 &&
    restored.paramCageSimplifier.paramRelocate.tangentialWeight == 2.5 &&
    restored.paramCageSimplifier.paramRelocate.surfaceWeight == 0.75 &&
    restored.paramCageSimplifier.paramRelocate.minTriangleQuality == 0.35,
    "quality parameters must round-trip through JSON");
  serialized.erase("phase2QualityPolishIterations");
  serialized.at("paramRelocate").as_object().erase("qualitySweeps");
  serialized.at("paramRelocate").as_object().erase("lineSearchMaxIter");
  serialized.at("paramRelocate").as_object().erase("tangentialWeight");
  serialized.at("paramRelocate").as_object().erase("surfaceWeight");
  serialized.at("paramRelocate").as_object().erase("minTriangleQuality");
  ParamCageGenerator legacy;
  legacy.paramCageSimplifier.deserialize(serialized);
  require(legacy.paramCageSimplifier.phase2QualityPolishIterations == 30 &&
    legacy.paramCageSimplifier.paramRelocate.qualitySweeps == 20 &&
    legacy.paramCageSimplifier.paramRelocate.lineSearchMaxIter == 12 &&
    legacy.paramCageSimplifier.paramRelocate.tangentialWeight == 1.0 &&
    legacy.paramCageSimplifier.paramRelocate.surfaceWeight == 1.0 &&
    legacy.paramCageSimplifier.paramRelocate.minTriangleQuality == 0.2,
    "legacy JSON must retain quality defaults when new fields are absent");
  serialized["phase2QualityPolishIterations"] = 0;
  serialized.at("paramRelocate").as_object()["qualitySweeps"] = 0;
  serialized.at("paramRelocate").as_object()["lineSearchMaxIter"] = 0;
  legacy.paramCageSimplifier.deserialize(serialized);
  require(legacy.paramCageSimplifier.phase2QualityPolishIterations == 0 &&
    legacy.paramCageSimplifier.paramRelocate.qualitySweeps == 0 &&
    legacy.paramCageSimplifier.paramRelocate.lineSearchMaxIter == 0,
    "zero must explicitly disable quality cycles, relocation sweeps, or line search");
  serialized["phase2QualityPolishIterations"] = -2;
  serialized.at("paramRelocate").as_object()["qualitySweeps"] = -3;
  serialized.at("paramRelocate").as_object()["lineSearchMaxIter"] = -4;
  legacy.paramCageSimplifier.deserialize(serialized);
  require(legacy.paramCageSimplifier.phase2QualityPolishIterations == 0 &&
    legacy.paramCageSimplifier.paramRelocate.qualitySweeps == 0 &&
    legacy.paramCageSimplifier.paramRelocate.lineSearchMaxIter == 0,
    "negative iteration limits must clamp to zero rather than wrap around");
  serialized.at("paramRelocate").as_object()["tangentialWeight"] = 0;
  serialized.at("paramRelocate").as_object()["surfaceWeight"] = 1;
  serialized.at("paramRelocate").as_object()["minTriangleQuality"] = 0;
  legacy.paramCageSimplifier.deserialize(serialized);
  require(legacy.paramCageSimplifier.paramRelocate.tangentialWeight == 0.0 &&
    legacy.paramCageSimplifier.paramRelocate.surfaceWeight == 1.0 &&
    legacy.paramCageSimplifier.paramRelocate.minTriangleQuality == 0.0,
    "relocation scalar settings must accept numeric JSON integers");
}

void invalid_relocation_energy_settings_rejected()
{
  Fixture fixture;
  const auto v = add_attraction_fan(fixture);
  const Vec3d old_point = fixture.cage.point(v[0]);
  auto& parameters = fixture.parameters.paramCageSimplifier.paramRelocate;
  const auto defaults = parameters;
  const auto require_invalid = [&]
  {
    bool rejected_direct = false;
    try { fixture.relocations().do_quality_relocate(12); }
    catch (const std::invalid_argument&) { rejected_direct = true; }
    require(rejected_direct && close(fixture.cage.point(v[0]), old_point),
      "invalid direct relocation settings must fail before mutating the mesh");
    bool rejected_json = false;
    try
    {
      ParamRelocateStage restored = defaults;
      restored.deserialize(parameters.serialize());
    }
    catch (const std::invalid_argument&) { rejected_json = true; }
    require(rejected_json, "invalid relocation settings must also fail JSON deserialization");
    parameters = defaults;
  };
  for (double invalid : {-1.0, std::numeric_limits<double>::infinity(),
    std::numeric_limits<double>::quiet_NaN()})
  {
    parameters.tangentialWeight = invalid;
    require_invalid();
    parameters.surfaceWeight = invalid;
    require_invalid();
  }
  for (double invalid : {-0.1, 1.1, std::numeric_limits<double>::infinity(),
    std::numeric_limits<double>::quiet_NaN()})
  {
    parameters.minTriangleQuality = invalid;
    require_invalid();
  }
}

double linear_solve_integration(bool enable_polish, size_t sweeps = 20, size_t line_attempts = 12)
{
  Fixture fixture;
  // Place the source tetrahedron strictly inside the closed octahedron.
  for (VertexHandle vertex : fixture.original.vertices())
    fixture.original.set_point(vertex,
      (fixture.original.point(vertex) - Vec3d(10.25, 10.25, 10.25)) * 0.2);
  const auto fan = add_fan(fixture.cage, Vec3d(0.6, 0.0, 1.0));
  const auto bottom = fixture.cage.add_vertex(Vec3d(0.0, 0.0, -1.0));
  for (size_t i = 1; i <= 4; ++i)
    add_face(fixture.cage, bottom, fan[i == 4 ? 1 : i + 1], fan[i]);
  for (VertexHandle vertex : fixture.cage.vertices())
    require(!fixture.cage.is_boundary(vertex), "integration cage must be closed");

  std::vector<Vec3d> before_points;
  std::set<std::pair<int, int>> before_edges;
  for (VertexHandle vertex : fixture.cage.vertices())
    before_points.push_back(fixture.cage.point(vertex));
  for (EdgeHandle edge : fixture.cage.edges())
  {
    const auto halfedge = fixture.cage.halfedge_handle(edge, 0);
    const int a = fixture.cage.from_vertex_handle(halfedge).idx();
    const int b = fixture.cage.to_vertex_handle(halfedge).idx();
    before_edges.emplace(std::min(a, b), std::max(a, b));
  }
  const double before_quality = min_quality(fixture.cage);
  const std::filesystem::path output_directory = std::filesystem::path("quality_polish_test_output") /
    (enable_polish ? "enabled_" + std::to_string(sweeps) + "_" + std::to_string(line_attempts) : "disabled");
  std::filesystem::create_directories(output_directory);
  fixture.parameters.setOutputPath(output_directory.generic_string() + "/", "octahedron");
  fixture.parameters.setCageLabel(0);
  auto& parameters = fixture.parameters.paramCageSimplifier;
  parameters.phase2Mode = "linear_solve";
  parameters.targetVerticesNum = 6;
  parameters.enableBoundaryRails = true;
  parameters.phase2QualityPolishIterations = enable_polish ? 1 : 0;
  parameters.paramRelocate.qualitySweeps = sweeps;
  parameters.paramRelocate.lineSearchMaxIter = line_attempts;
  // This integration isolates repeated tangential smoothing. Separate
  // operator and full-pipeline tests exercise the default source attraction.
  parameters.paramRelocate.surfaceWeight = 0.0;
  CageSimplifier simplifier(&fixture.original, &fixture.cage, &parameters);
  simplifier.simplify();

  require(fixture.cage.n_vertices() == 6 && fixture.cage.n_faces() == 8,
    "already-at-target integration must preserve vertex and face counts");
  bool position_changed = false;
  for (VertexHandle vertex : fixture.cage.vertices())
    position_changed = position_changed || !close(fixture.cage.point(vertex), before_points[vertex.idx()]);
  if (enable_polish && sweeps > 0 && line_attempts > 0)
  {
    require(position_changed, "linear mode must run relocation even when collapse target was already met");
    require(min_quality(fixture.cage) > before_quality + 1e-6,
      "integrated linear quality polish must improve minimum triangle quality");
  }
  else
  {
    require(!position_changed, "disabled relocation must preserve at-target vertex positions");
  }
  if (!enable_polish)
  {
    for (EdgeHandle edge : fixture.cage.edges())
    {
      const auto halfedge = fixture.cage.halfedge_handle(edge, 0);
      const int a = fixture.cage.from_vertex_handle(halfedge).idx();
      const int b = fixture.cage.to_vertex_handle(halfedge).idx();
      require(before_edges.count({std::min(a, b), std::max(a, b)}) == 1,
        "zero polish rounds must preserve connectivity");
    }
  }
  for (size_t i = 1; i <= 4; ++i)
  {
    require(close(fixture.cage.point(fan[i]), before_points[fan[i].idx()]),
      "integrated polish must keep rail positions fixed");
    const auto halfedge = fixture.cage.find_halfedge(fan[i], fan[i == 4 ? 1 : i + 1]);
    require(halfedge.is_valid() && fixture.cage.data(fixture.cage.edge_handle(halfedge)).boundary_rail_id == 3,
      "integrated polish must preserve the labeled rail cycle");
  }
  require(std::filesystem::exists(output_directory / "octahedron_debug_phase2_linear_solve.obj"),
    "integration must complete the production Phase 2 output path");
  fixture.initialize();
  fixture.require_source_clear();
  return min_quality(fixture.cage);
}

void repeated_relocation_sweeps_improve_quality()
{
  // One outer cycle isolates the number of inner smoothing sweeps. This
  // closed cage is already at target, and its four rail vertices stay fixed.
  const double disabled = linear_solve_integration(true, 0);
  const double no_line_search = linear_solve_integration(true, 20, 0);
  const double single = linear_solve_integration(true, 1);
  const double repeated = linear_solve_integration(true, 20);
  require(std::abs(disabled - no_line_search) < 1e-12,
    "zero sweeps and zero line-search attempts must both disable relocation");
  require(single > disabled + 1e-6, "the first relocation sweep must improve quality");
  require(repeated > single + 1e-6,
    "additional relocation sweeps must improve beyond the one-sweep result");
}

size_t require_closed_rail_cycles(SMeshT& mesh)
{
  std::map<int, std::set<int>> rail_vertices;
  std::map<int, std::map<int, std::set<int>>> adjacency;
  for (VertexHandle vertex : mesh.vertices())
    if (mesh.data(vertex).boundary_rail_id >= 0)
      rail_vertices[mesh.data(vertex).boundary_rail_id].insert(vertex.idx());
  for (EdgeHandle edge : mesh.edges())
  {
    const int rail = mesh.data(edge).boundary_rail_id;
    if (rail < 0)
      continue;
    const auto halfedge = mesh.halfedge_handle(edge, 0);
    const auto a = mesh.from_vertex_handle(halfedge);
    const auto b = mesh.to_vertex_handle(halfedge);
    require(mesh.data(a).boundary_rail_id == rail && mesh.data(b).boundary_rail_id == rail,
      "generated rail edge must join vertices with its label");
    adjacency[rail][a.idx()].insert(b.idx());
    adjacency[rail][b.idx()].insert(a.idx());
  }
  require(adjacency.size() == rail_vertices.size(), "every rail vertex must belong to a labeled edge cycle");
  for (const auto& rail : rail_vertices)
  {
    require(rail.second.size() >= 3, "rail cycles require at least three vertices");
    for (int vertex : rail.second)
      require(adjacency[rail.first][vertex].size() == 2, "rail vertices must have exactly two rail neighbors");
    std::set<int> visited;
    std::vector<int> pending = {*rail.second.begin()};
    while (!pending.empty())
    {
      const int vertex = pending.back();
      pending.pop_back();
      if (!visited.insert(vertex).second)
        continue;
      for (int next : adjacency[rail.first][vertex])
        pending.push_back(next);
    }
    require(visited == rail.second, "each rail label must form one connected closed cycle");
  }
  return rail_vertices.size();
}

std::array<VertexHandle, 6> add_collapse_octahedron(SMeshT& mesh, bool rail_first = false,
  bool subdivide_top_face = false)
{
  std::array<VertexHandle, 6> v = {
    mesh.add_vertex(Vec3d(0.0, 0.0, 1.0)),
    mesh.add_vertex(Vec3d(1.0, 0.0, 0.0)),
    mesh.add_vertex(Vec3d(0.0, 1.0, 0.0)),
    mesh.add_vertex(Vec3d(-1.0, 0.0, 0.0)),
    mesh.add_vertex(Vec3d(0.0, -1.0, 0.0)),
    mesh.add_vertex(Vec3d(0.0, 0.0, -1.0))
  };
  // Rotating face insertion makes halfedge 0 of top/south point either
  // rail -> ordinary or ordinary -> rail without changing the surface.
  for (size_t j = 0; j < 4; ++j)
  {
    const size_t i = (j + (rail_first ? 3 : 0)) % 4 + 1;
    const size_t next = i == 4 ? 1 : i + 1;
    if (subdivide_top_face && i == 1)
    {
      const auto center = mesh.add_vertex((mesh.point(v[0]) + mesh.point(v[1]) +
        mesh.point(v[2])) / 3.0);
      add_face(mesh, v[0], v[1], center);
      add_face(mesh, v[1], v[2], center);
      add_face(mesh, v[2], v[0], center);
    }
    else
      add_face(mesh, v[0], v[i], v[next]);
  }
  for (size_t i = 1; i <= 4; ++i)
    add_face(mesh, v[5], v[i == 4 ? 1 : i + 1], v[i]);
  return v;
}

void label_rail_cycle(SMeshT& mesh, const std::vector<VertexHandle>& vertices, int rail)
{
  for (size_t i = 0; i < vertices.size(); ++i)
  {
    mesh.data(vertices[i]).boundary_rail_id = rail;
    const auto halfedge = mesh.find_halfedge(vertices[i], vertices[(i + 1) % vertices.size()]);
    require(halfedge.is_valid(), "fixture rail cycle must follow mesh edges");
    mesh.data(mesh.edge_handle(halfedge)).boundary_rail_id = rail;
  }
}

EdgeCollapser collapse_operator(Fixture& fixture)
{
  EdgeCollapser collapser(&fixture.original, &fixture.cage, fixture.source_tree.get(),
    fixture.cage_tree.get(), fixture.source_grid.get());
  collapser.set_flags(false, false, true, true, true, true);
  return collapser;
}

void mixed_collapse_keeps_rail_vertex(bool rail_first, bool with_chord = false)
{
  Fixture fixture;
  const auto v = add_collapse_octahedron(fixture.cage, rail_first);
  const std::vector<VertexHandle> rail = with_chord ?
    std::vector<VertexHandle>{v[0], v[1], v[2], v[3]} :
    std::vector<VertexHandle>{v[0], v[1], v[2]};
  label_rail_cycle(fixture.cage, rail, 3);

  // An actual rail-edge collapse would be projected onto these remote
  // source support patches at z = 20. A mixed collapse must ignore them.
  add_triangle(fixture.original, Vec3d(20.0, 20.0, 20.0),
    Vec3d(22.0, 20.0, 20.0), Vec3d(20.0, 22.0, 20.0));
  for (EdgeHandle edge : fixture.original.edges())
    if (fixture.original.is_boundary(edge))
    {
      fixture.original.data(edge).boundary_rail_id = 3;
      fixture.original.data(edge).boundary_rail_outer_direction = Vec3d(0.0, 1.0, 0.0);
    }
  std::vector<Vec3d> before_points;
  for (VertexHandle vertex : v)
    before_points.push_back(fixture.cage.point(vertex));
  const ExactPoint original_exact(before_points[0]);
  fixture.cage.data(v[0]).ep = std::make_unique<ExactPoint>(original_exact);
  fixture.initialize();
  fixture.require_source_clear();
  require(require_closed_rail_cycles(fixture.cage) == 1, "fixture must have one closed rail");

  const auto edge = fixture.cage.edge_handle(fixture.cage.find_halfedge(v[4], v[0]));
  const auto halfedge0 = fixture.cage.halfedge_handle(edge, 0);
  require((fixture.cage.from_vertex_handle(halfedge0) == v[0]) == rail_first,
    "fixture must exercise the requested halfedge-0 direction");
  auto collapser = collapse_operator(fixture);
  require(collapser.init(edge), "mixed edge must initialize even on a three-vertex rail");
  require(collapser.has_fixed_rail_target(), "mixed edge must select a fixed rail target");
  const Vec3d proposed(0.3, -0.2, 0.4);
  const ExactPoint unrelated_exact(proposed);
  require(collapser.constrained_target_point(proposed) == before_points[0],
    "mixed target must be the original rail point, independent of proposal and source support");
  require(collapser.target_point_is_valid(proposed, &unrelated_exact),
    "validity checks must evaluate the fixed rail position");
  require(collapser.try_collapse_edge(proposed, &unrelated_exact),
    "valid mixed edge must collapse into its existing rail vertex");
  require(collapser.get_collapsed_center() == v[0] && !fixture.cage.status(v[0]).deleted() &&
    fixture.cage.status(v[4]).deleted(), "collapse must retain the rail handle and delete the ordinary handle");
  require(fixture.cage.data(v[0]).ep && fixture.cage.data(v[0]).ep->exact() == original_exact.exact(),
    "mixed collapse must preserve the rail's exact coordinate, ignoring the proposal's exact point");
  for (VertexHandle vertex : fixture.cage.vertices())
    require(fixture.cage.point(vertex) == before_points[vertex.idx()],
      "mixed collapse must leave every surviving vertex at its original position");
  for (size_t i = 0; i < rail.size(); ++i)
  {
    const auto halfedge = fixture.cage.find_halfedge(rail[i], rail[(i + 1) % rail.size()]);
    require(fixture.cage.data(rail[i]).boundary_rail_id == 3 && halfedge.is_valid() &&
      fixture.cage.data(fixture.cage.edge_handle(halfedge)).boundary_rail_id == 3,
      "mixed collapse must preserve the rail ID and each original cycle connection");
  }
  if (with_chord)
  {
    const auto chord = fixture.cage.edge_handle(fixture.cage.find_halfedge(v[0], v[2]));
    require(fixture.cage.data(chord).boundary_rail_id < 0,
      "a same-rail chord must not become a rail edge after mixed collapse");
  }
  require(require_closed_rail_cycles(fixture.cage) == 1, "mixed collapse must preserve a closed rail cycle");
  fixture.require_source_clear();
  fixture.cage.garbage_collection();
  require(fixture.cage.n_vertices() == 5 && fixture.cage.n_faces() == 6,
    "mixed collapse must remove exactly one vertex and two faces");
}

void prohibited_rail_collapses_stay_rejected()
{
  for (int kind = 0; kind < 3; ++kind)
  {
    Fixture fixture;
    const auto v = add_collapse_octahedron(fixture.cage);
    EdgeHandle edge;
    if (kind == 0)
    {
      label_rail_cycle(fixture.cage, {v[0], v[1], v[2]}, 3);
      edge = fixture.cage.edge_handle(fixture.cage.find_halfedge(v[0], v[1]));
    }
    else if (kind == 1)
    {
      label_rail_cycle(fixture.cage, {v[0], v[1], v[2], v[3]}, 3);
      edge = fixture.cage.edge_handle(fixture.cage.find_halfedge(v[0], v[2]));
    }
    else
    {
      label_rail_cycle(fixture.cage, {v[0], v[1], v[2]}, 3);
      label_rail_cycle(fixture.cage, {v[5], v[3], v[4]}, 4);
      edge = fixture.cage.edge_handle(fixture.cage.find_halfedge(v[0], v[4]));
    }
    fixture.initialize();
    require(fixture.cage.is_collapse_ok(fixture.cage.halfedge_handle(edge, 0)),
      "forbidden fixture edge must be otherwise topologically collapsible");
    auto collapser = collapse_operator(fixture);
    require(!collapser.init(edge),
      "three-vertex rail edges, same-rail chords, and different-rail edges must remain forbidden");
    require(require_closed_rail_cycles(fixture.cage) == (kind == 2 ? 2 : 1),
      "rejected rail collapse must retain its original cycles");
  }
}

void mixed_collapse_collision_is_rejected()
{
  Fixture fixture;
  const auto v = add_collapse_octahedron(fixture.cage);
  label_rail_cycle(fixture.cage, {v[0], v[1], v[2]}, 3);
  // The fixed collapse creates triangle top/east/bottom in y = 0. This
  // small source triangle is strictly inside the original octahedron.
  const Vec3d obstacle(1.0 / 3.0, 0.0, 0.0);
  add_triangle(fixture.original, obstacle + Vec3d(0.0, -0.01, 0.0),
    obstacle + Vec3d(0.0, 0.01, 0.0), obstacle + Vec3d(0.01, 0.0, 0.0));
  fixture.initialize();
  fixture.require_source_clear();
  const auto edge = fixture.cage.edge_handle(fixture.cage.find_halfedge(v[4], v[0]));
  auto collapser = collapse_operator(fixture);
  require(collapser.init(edge) && collapser.has_fixed_rail_target(),
    "colliding mixed edge must first pass topology and rail checks");
  require(!collapser.target_point_is_valid(fixture.cage.point(v[4]), nullptr) &&
    !collapser.try_collapse_edge(fixture.cage.point(v[4]), nullptr),
    "fixed rail placement must still reject source intersection");
  require(!fixture.cage.status(v[0]).deleted() && !fixture.cage.status(v[4]).deleted() &&
    require_closed_rail_cycles(fixture.cage) == 1, "collision rejection must preserve both endpoints and the rail");
}

void linear_solve_collapses_only_mixed_edges()
{
  Fixture fixture;
  const auto v = add_collapse_octahedron(fixture.cage, false, true);
  label_rail_cycle(fixture.cage, {v[0], v[1], v[2]}, 3);
  label_rail_cycle(fixture.cage, {v[5], v[3], v[4]}, 4);
  std::array<Vec3d, 6> before;
  for (size_t i = 0; i < v.size(); ++i)
    before[i] = fixture.cage.point(v[i]);
  fixture.initialize();
  // Every old vertex belongs to one of two three-vertex rails. Their edges
  // and cross-rail connections are forbidden; only the face subdivision's
  // ordinary vertex can collapse. The zero proxy has no unique solve.
  auto& parameters = fixture.parameters.paramCageSimplifier.paramCollapse;
  parameters.phase2PlacementStrategy = "linear_solve";
  parameters.qemWeight = 0.0;
  parameters.triangleQualityWeight = 0.0;
  parameters.uniformityWeight = 0.0;
  parameters.phase2LinearSolveCollisionReject = true;
  parameters.lineSearchMaxIter = 0;
  CollapseStage stage(&fixture.original, &fixture.cage, &parameters, fixture.source_tree.get(),
    fixture.cage_tree.get(), fixture.source_grid.get(), 1.0);
  require(stage.do_phase2_energy_simplification(6) == 1,
    "linear Phase 2 must queue and commit a mixed edge even with a singular zero-weight proxy");
  require(fixture.cage.n_vertices() == 6 && fixture.cage.n_faces() == 8 &&
    require_closed_rail_cycles(fixture.cage) == 2,
    "mixed-only Phase 2 must remove the free vertex while retaining both minimal rail cycles");
  for (size_t i = 0; i < v.size(); ++i)
    require(fixture.cage.point(v[i]) == before[i], "linear mixed collapse must keep all six rail positions exactly");
  fixture.require_source_clear();
}

void require_clear_closed_mesh(SMeshT& mesh, SMeshT& original)
{
  DFaceTree source_tree(original);
  std::vector<std::array<VertexHandle, 3>> faces;
  for (VertexHandle vertex : mesh.vertices())
  {
    const auto& p = mesh.point(vertex);
    require(std::isfinite(p[0]) && std::isfinite(p[1]) && std::isfinite(p[2]),
      "generated cage coordinates must remain finite");
    require(mesh.is_manifold(vertex) && !mesh.is_boundary(vertex),
      "generated cage must remain a closed manifold");
  }
  for (FaceHandle face : mesh.faces())
  {
    std::array<VertexHandle, 3> v;
    size_t index = 0;
    for (VertexHandle vertex : mesh.fv_range(face))
      v[index++] = vertex;
    require(!Geometry::are_points_colinear(mesh.point(v[0]), mesh.point(v[1]), mesh.point(v[2]),
      mesh.data(v[0]).ep.get(), mesh.data(v[1]).ep.get(), mesh.data(v[2]).ep.get()),
      "generated cage faces must not be exactly degenerate");
    require(!source_tree.do_intersect(mesh.point(v[0]), mesh.point(v[1]), mesh.point(v[2]),
      mesh.data(v[0]).ep.get(), mesh.data(v[1]).ep.get(), mesh.data(v[2]).ep.get()),
      "generated cage must not intersect the original tube");
    faces.push_back(v);
  }
  // Check all triangle pairs, including adjacent pairs: shared edges and
  // vertices are allowed, overlap beyond the shared simplex is forbidden.
  for (size_t i = 0; i < faces.size(); ++i)
    for (size_t j = i + 1; j < faces.size(); ++j)
    {
      std::vector<VertexHandle> common, a, b;
      for (VertexHandle vertex : faces[i])
      {
        if (std::find(faces[j].begin(), faces[j].end(), vertex) != faces[j].end())
          common.push_back(vertex);
        else
          a.push_back(vertex);
      }
      for (VertexHandle vertex : faces[j])
        if (std::find(common.begin(), common.end(), vertex) == common.end())
          b.push_back(vertex);
      const auto p = [&](VertexHandle v) -> const Vec3d& { return mesh.point(v); };
      const auto ep = [&](VertexHandle v) { return mesh.data(v).ep.get(); };
      bool intersects = true;
      if (common.empty())
        intersects = Geometry::triangle_do_intersect(p(a[0]), p(a[1]), p(a[2]), p(b[0]), p(b[1]), p(b[2]),
          ep(a[0]), ep(a[1]), ep(a[2]), ep(b[0]), ep(b[1]), ep(b[2]));
      else if (common.size() == 1)
        intersects = Geometry::triangle_do_intersect(p(a[0]), p(a[1]), p(common[0]), p(b[0]), p(b[1]),
          ep(a[0]), ep(a[1]), ep(common[0]), ep(b[0]), ep(b[1]));
      else if (common.size() == 2)
        intersects = Geometry::triangle_do_overlap(p(a[0]), p(common[0]), p(common[1]), p(b[0]),
          ep(a[0]), ep(common[0]), ep(common[1]), ep(b[0]));
      require(!intersects, "generated cage must not contain adjacent or nonadjacent self-intersections");
    }
}

struct PipelineSnapshot
{
  std::vector<Vec3d> points;
  std::vector<std::array<int, 3>> faces;
  std::vector<std::array<int, 3>> rail_edges;
  std::string log;
};

class ScopedUserLogCapture
{
  std::shared_ptr<spdlog::logger> previous_logger;
  std::ostringstream output;
public:
  ScopedUserLogCapture() : previous_logger(Logger::user_logger)
  {
    auto sink = std::make_shared<spdlog::sinks::ostream_sink_mt>(output);
    Logger::user_logger = std::make_shared<spdlog::logger>("quality-test-capture", sink);
    Logger::user_logger->set_pattern("%v");
    Logger::user_logger->set_level(spdlog::level::info);
  }
  ~ScopedUserLogCapture() { Logger::user_logger = previous_logger; }
  std::string text() const { return output.str(); }
};

size_t full_topological_offset_rail_pipeline(size_t outer_cycles = 3,
  size_t target_vertices = 40, double surface_weight = 1.0, size_t relocation_sweeps = 20,
  PipelineSnapshot* snapshot = nullptr)
{
  CageGenerator generator;
  generator.originalMesh = std::make_unique<SMeshT>();
  std::array<std::array<VertexHandle, 6>, 2> rings;
  for (size_t ring = 0; ring < 2; ++ring)
    for (size_t side = 0; side < 6; ++side)
    {
      const double angle = 2.0 * M_PI * static_cast<double>(side) / 6.0;
      rings[ring][side] = generator.originalMesh->add_vertex(
        Vec3d(std::cos(angle), std::sin(angle), ring == 0 ? -0.5 : 0.5));
    }
  for (size_t side = 0; side < 6; ++side)
  {
    const size_t next = (side + 1) % 6;
    add_face(*generator.originalMesh, rings[0][side], rings[0][next], rings[1][next]);
    add_face(*generator.originalMesh, rings[0][side], rings[1][next], rings[1][side]);
  }
  const std::filesystem::path output_directory = std::filesystem::path("quality_polish_test_output") /
    ("full_pipeline_" + std::to_string(outer_cycles) + "_" + std::to_string(target_vertices) +
      "_surface_" + std::to_string(surface_weight) + "_sweeps_" + std::to_string(relocation_sweeps));
  std::filesystem::create_directories(output_directory);
  generator.param.setOutputPath(output_directory.generic_string() + "/", "open_tube");
  generator.param.setCageLabel(0);
  generator.param.paramCageInitializer.phase1Mode = "topological_offset";
  auto& lattice = generator.param.paramCageInitializer.paramTetrahedralizer.paramLatticePointsGenerator;
  lattice.minDistanceFactor = 0.1;
  lattice.pointNumAlongAxis = 3;
  lattice.areaThresholdRate = 16.0;
  auto& simplification = generator.param.paramCageSimplifier;
  simplification.phase2Mode = "linear_solve";
  simplification.enableBoundaryRails = true;
  simplification.targetVerticesNum = target_vertices;
  simplification.phase2QualityPolishIterations = outer_cycles;
  simplification.paramRelocate.surfaceWeight = surface_weight;
  simplification.paramRelocate.qualitySweeps = relocation_sweeps;
  // Match the CLI linear-solve preset; only Phase 1 sampling is coarser.
  simplification.paramCollapse.uniformityMode = "global";
  simplification.paramCollapse.qemWeight = 1.0;
  simplification.paramCollapse.triangleQualityWeight = 1.0;
  simplification.paramCollapse.uniformityWeight = 1.0;
  omp_set_num_threads(2);
  const auto start = std::chrono::steady_clock::now();
  generator.stageInitialize();
  generator.stageBuildBoundaryRails();
  const auto phase1_end = std::chrono::steady_clock::now();
  require(require_closed_rail_cycles(*generator.cage) == 2,
    "both open-tube source boundaries must produce rail cycles in Phase 1");
  const size_t phase1_vertices = generator.cage->n_vertices();
  require(phase1_vertices > simplification.targetVerticesNum,
    "full pipeline fixture must exercise actual Phase 2 collapses");
  generator.cageInitializer.reset();
  generator.VMesh.reset();
  ScopedUserLogCapture captured_log;
  generator.stageSimplify();
  const auto phase2_end = std::chrono::steady_clock::now();
  require(generator.cage->n_vertices() < phase1_vertices,
    "full pipeline must perform collapses before quality polish");
  require(require_closed_rail_cycles(*generator.cage) == 2,
    "both generated rail cycles must survive collapse and quality polish");
  require_clear_closed_mesh(*generator.cage, *generator.originalMesh);
  require(std::filesystem::exists(output_directory / "open_tube_debug_phase2_linear_solve.obj"),
    "full pipeline must write the final linear-solve cage OBJ");
  if (snapshot)
  {
    *snapshot = PipelineSnapshot();
    snapshot->log = captured_log.text();
    for (VertexHandle vertex : generator.cage->vertices())
      snapshot->points.push_back(generator.cage->point(vertex));
    for (FaceHandle face : generator.cage->faces())
    {
      std::array<int, 3> corners;
      size_t index = 0;
      for (VertexHandle vertex : generator.cage->fv_range(face))
        corners[index++] = vertex.idx();
      snapshot->faces.push_back(corners);
    }
    for (EdgeHandle edge : generator.cage->edges())
    {
      const int rail = generator.cage->data(edge).boundary_rail_id;
      if (rail < 0)
        continue;
      const auto halfedge = generator.cage->halfedge_handle(edge, 0);
      const int a = generator.cage->from_vertex_handle(halfedge).idx();
      const int b = generator.cage->to_vertex_handle(halfedge).idx();
      snapshot->rail_edges.push_back({std::min(a, b), std::max(a, b), rail});
    }
    std::sort(snapshot->rail_edges.begin(), snapshot->rail_edges.end());
  }
  std::cout << "Full pipeline (cycles " << outer_cycles << ", target " << target_vertices << "): "
    << phase1_vertices << " -> " << generator.cage->n_vertices()
    << " vertices, min quality " << min_quality(*generator.cage)
    << ", Phase 1 " << std::chrono::duration<double>(phase1_end - start).count()
    << " s, Phase 2 " << std::chrono::duration<double>(phase2_end - phase1_end).count() << " s\n";
  return generator.cage->n_vertices();
}

void require_final_relocation_order(const PipelineSnapshot& snapshot, const char* stop_reason)
{
  const std::string& log = snapshot.log;
  const size_t first_relocation = log.find("energy relocation:");
  const size_t last_cycle = log.rfind("phase 2 linear-solve cycle ");
  const size_t stop = log.find(stop_reason);
  const std::string summary_prefix = "phase 2 linear-solve final relocation: moves ";
  const size_t final_summary = log.find(summary_prefix);
  require(first_relocation != std::string::npos && last_cycle < first_relocation &&
    stop < first_relocation && final_summary > first_relocation && final_summary != std::string::npos,
    "all collapse/flip cycles must finish before the single final relocation stage");
  require(log.find("phase 2 linear-solve cycle ", first_relocation) == std::string::npos &&
    log.find("phase 2 energy simplification", first_relocation) == std::string::npos &&
    log.find("quality flip:", first_relocation) == std::string::npos &&
    log.find(summary_prefix, final_summary + summary_prefix.size()) == std::string::npos,
    "collapse and flip must not resume after relocation, and the final stage must run only once");
  size_t moves = 0;
  std::istringstream(log.substr(final_summary + summary_prefix.size())) >> moves;
  require(moves > 0, "final-only ordering fixture must accept actual vertex moves");
}

void relocation_runs_only_after_topology_selection()
{
  PipelineSnapshot disabled, attracted, tangent_only, capped, collapse_only;
  full_topological_offset_rail_pipeline(6, 8, 1.0, 0, &disabled);
  full_topological_offset_rail_pipeline(6, 8, 1.0, 20, &attracted);
  full_topological_offset_rail_pipeline(6, 8, 0.0, 20, &tangent_only);
  require(attracted.points.size() > 8,
    "ordering fixture must exercise a collapse/flip loop that stalls above its target");
  require(attracted.log.find("phase 2 linear-solve cycle 2:") != std::string::npos,
    "ordering fixture must exercise more than one collapse/flip cycle");
  require(disabled.points.size() == attracted.points.size() &&
    disabled.points.size() == tangent_only.points.size() &&
    disabled.faces == attracted.faces && disabled.faces == tangent_only.faces &&
    disabled.rail_edges == attracted.rail_edges && disabled.rail_edges == tangent_only.rail_edges,
    "final mesh topology and rail cycles must be independent of relocation settings");
  bool moved = false;
  for (size_t i = 0; i < disabled.points.size(); ++i)
    moved = moved || !close(disabled.points[i], attracted.points[i]);
  require(moved, "default final relocation must change positions after topology selection");
  for (const auto& rail : disabled.rail_edges)
    for (size_t endpoint = 0; endpoint < 2; ++endpoint)
      require(close(disabled.points[rail[endpoint]], attracted.points[rail[endpoint]]) &&
        close(disabled.points[rail[endpoint]], tangent_only.points[rail[endpoint]]),
        "final relocation must preserve every rail vertex position");
  require_final_relocation_order(attracted, "no accepted collapses or flips.");
  require_final_relocation_order(tangent_only, "no accepted collapses or flips.");
  require(disabled.log.find("energy relocation:") == std::string::npos,
    "zero sweeps must disable final relocation while retaining collapse/flip cycles");

  full_topological_offset_rail_pipeline(1, 8, 1.0, 20, &capped);
  require_final_relocation_order(capped, "cycle limit reached.");
  full_topological_offset_rail_pipeline(0, 8, 1.0, 20, &collapse_only);
  require(collapse_only.log.find("phase 2 energy simplification") != std::string::npos &&
    collapse_only.log.find("quality flip:") == std::string::npos &&
    collapse_only.log.find("energy relocation:") == std::string::npos,
    "zero outer cycles must preserve the collapse-only compatibility path");
}
}

int main()
{
  Logger::InitLogger(spdlog::level::err, false);
  const std::vector<std::pair<const char*, void(*)()>> tests = {
    {"quality flip improves triangles without relocating vertices", improving_flip},
    {"rail edge cannot flip", rail_flip_rejected},
    {"concave flip cannot overlap replacement faces", overlapping_flip_rejected},
    {"source-intersecting flip is rejected", [] { intersecting_flip_rejected(true); }},
    {"self-intersecting flip is rejected", [] { intersecting_flip_rejected(false); }},
    {"mixed collapse keeps the rail endpoint when halfedge 0 removes the ordinary vertex", [] { mixed_collapse_keeps_rail_vertex(false); }},
    {"mixed collapse reverses halfedge 0 to retain the rail endpoint", [] { mixed_collapse_keeps_rail_vertex(true); }},
    {"mixed collapse preserves rail cycles without promoting chords", [] { mixed_collapse_keeps_rail_vertex(true, true); }},
    {"minimal rail edges, chords, and cross-rail edges remain forbidden", prohibited_rail_collapses_stay_rejected},
    {"mixed collapse still rejects collisions at the fixed rail point", mixed_collapse_collision_is_rejected},
    {"linear solve queues and commits mixed edges with fixed rail positions", linear_solve_collapses_only_mixed_edges},
    {"Voronoi smoothing includes neighbors' exterior faces", voronoi_target_uses_neighbor_faces},
    {"tangential relocation improves quality and fixes rail vertices", improving_relocation},
    {"symmetric smoothing is a no-op", symmetric_relocation_is_noop},
    {"tiny distant-surface weight preserves tangential energy descent", tiny_surface_weight_preserves_tangential_descent},
    {"nonfinite and degenerate relocation targets are rejected", invalid_relocation_is_rejected},
    {"relocation backtracks after a source collision", collision_backtracking},
    {"surface attraction moves normally and decreases frozen quadratic energy", surface_attraction_decreases_frozen_energy},
    {"relocation weights disable terms and distance-only moves backtrack on contact", relocation_weights_and_source_collision},
    {"relocation backtracks to preserve its configured quality floor", relocation_quality_floor_backtracks},
    {"relocation cannot worsen fans already below the quality floor", relocation_preserves_quality_below_floor},
    {"quality configuration round-trips and supports legacy JSON", quality_configuration_roundtrip},
    {"invalid relocation weights and quality floors are rejected", invalid_relocation_energy_settings_rejected},
    {"linear mode polishes an already-at-target closed rail cage", [] { linear_solve_integration(true); }},
    {"zero polish rounds disable integrated quality operations", [] { linear_solve_integration(false); }},
    {"relocation sweep limits control repeated improvement at target", repeated_relocation_sweeps_improve_quality},
    {"topological offset, linear solve, and generated rails complete together", [] { full_topological_offset_rail_pipeline(); }},
    {"relocation runs once after topology selection and cannot affect its result", relocation_runs_only_after_topology_selection}
  };
  size_t failures = 0;
  for (const auto& test : tests)
  {
    try
    {
      test.second();
      std::cout << "PASS: " << test.first << '\n';
    }
    catch (const std::exception& exception)
    {
      ++failures;
      std::cerr << "FAIL: " << test.first << ": " << exception.what() << '\n';
    }
  }
  std::cout << (tests.size() - failures) << '/' << tests.size() << " quality-polish tests passed\n";
  return failures == 0 ? 0 : 1;
}
