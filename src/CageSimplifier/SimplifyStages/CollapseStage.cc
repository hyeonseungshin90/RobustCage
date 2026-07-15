#include "CollapseStage.hh"
#include "Dense"
#include <autodiff/reverse/var.hpp>
#include <autodiff/reverse/var/eigen.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <random>

#define Np 50
#define MUTATION_SCALE 0.9
#define CROSSOVER_RATE 0.9
#define CONVERGENCE_RATE 1e-4
#define MAX_CONSECUTIVE_ITER 5

namespace Cage
{
namespace CageSimp
{
namespace
{
double sqr(double v)
{
  return v * v;
}

double clamp_value(double v, double lo, double hi)
{
  return std::max(lo, std::min(hi, v));
}

bool finite_vec(const Vec3d& v)
{
  return std::isfinite(v.x()) && std::isfinite(v.y()) && std::isfinite(v.z());
}

Vec3d normalized_or_fallback(const Vec3d& v, const Vec3d& fallback)
{
  const double n = v.norm();
  if (n <= 0.0 || !std::isfinite(n))
    return fallback;
  return v / n;
}

Vec3d triangle_normal(const Vec3d& p0, const Vec3d& p1, const Vec3d& p2)
{
  return normalized_or_fallback((p1 - p0).cross(p2 - p0), Vec3d(0.0, 0.0, 1.0));
}

double triangle_area_from_points(const Vec3d& p0, const Vec3d& p1, const Vec3d& p2)
{
  return 0.5 * (p1 - p0).cross(p2 - p0).norm();
}

double triangle_quality_from_points(const Vec3d& p0, const Vec3d& p1, const Vec3d& p2)
{
  const double a = (p1 - p0).length();
  const double b = (p2 - p1).length();
  const double c = (p0 - p2).length();
  const double denom = a * a + b * b + c * c;
  if (denom <= 0.0)
    return 0.0;
  return 4.0 * std::sqrt(3.0) * triangle_area_from_points(p0, p1, p2) / denom;
}

Eigen::Vector3d to_eigen(const Vec3d& v)
{
  return Eigen::Vector3d(v.x(), v.y(), v.z());
}

Vec3d from_eigen(const Eigen::Vector3d& v)
{
  return Vec3d(v.x(), v.y(), v.z());
}

double face_average_edge_length(SMeshT* mesh, FaceHandle fh)
{
  double length = 0.0;
  size_t count = 0;
  for (EdgeHandle eh : mesh->fe_range(fh))
  {
    length += mesh->data(eh).edge_length;
    count++;
  }
  return count == 0 ? 0.0 : length / static_cast<double>(count);
}

void collect_face_points(SMeshT* mesh, FaceHandle fh, Vec3d pts[3])
{
  size_t idx = 0;
  for (VertexHandle vh : mesh->fv_range(fh))
  {
    if (idx < 3)
      pts[idx++] = mesh->point(vh);
  }
}

template<typename Scalar>
struct DiffVec3
{
  Scalar x;
  Scalar y;
  Scalar z;

  DiffVec3() : x(0.0), y(0.0), z(0.0) {}
  DiffVec3(const Scalar& _x, const Scalar& _y, const Scalar& _z) : x(_x), y(_y), z(_z) {}
};

template<typename Scalar>
DiffVec3<Scalar> diff_vec(const Vec3d& v)
{
  return DiffVec3<Scalar>(v.x(), v.y(), v.z());
}

template<typename Scalar>
DiffVec3<Scalar> operator+(const DiffVec3<Scalar>& a, const DiffVec3<Scalar>& b)
{
  return DiffVec3<Scalar>(a.x + b.x, a.y + b.y, a.z + b.z);
}

template<typename Scalar>
DiffVec3<Scalar> operator-(const DiffVec3<Scalar>& a, const DiffVec3<Scalar>& b)
{
  return DiffVec3<Scalar>(a.x - b.x, a.y - b.y, a.z - b.z);
}

template<typename Scalar>
DiffVec3<Scalar> operator*(const DiffVec3<Scalar>& v, double s)
{
  return DiffVec3<Scalar>(v.x * s, v.y * s, v.z * s);
}

template<typename Scalar>
DiffVec3<Scalar> operator*(const DiffVec3<Scalar>& v, const Scalar& s)
{
  return DiffVec3<Scalar>(v.x * s, v.y * s, v.z * s);
}

template<typename Scalar>
DiffVec3<Scalar> operator*(double s, const DiffVec3<Scalar>& v)
{
  return v * s;
}

template<typename Scalar>
DiffVec3<Scalar> operator*(const Scalar& s, const DiffVec3<Scalar>& v)
{
  return v * s;
}

template<typename Scalar>
DiffVec3<Scalar> operator/(const DiffVec3<Scalar>& v, double s)
{
  return DiffVec3<Scalar>(v.x / s, v.y / s, v.z / s);
}

template<typename Scalar>
Scalar diff_sqr(const Scalar& v)
{
  return v * v;
}

template<typename Scalar>
Scalar diff_dot(const DiffVec3<Scalar>& a, const DiffVec3<Scalar>& b)
{
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

template<typename Scalar>
DiffVec3<Scalar> diff_cross(const DiffVec3<Scalar>& a, const DiffVec3<Scalar>& b)
{
  return DiffVec3<Scalar>(
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x);
}

template<typename Scalar>
Scalar diff_norm_sqr(const DiffVec3<Scalar>& v)
{
  return diff_dot(v, v);
}

template<typename Scalar>
Scalar diff_safe_sqrt(const Scalar& v)
{
  using std::sqrt;
  const double eps = 1e-30;
  if (v <= eps)
    return Scalar(std::sqrt(eps));
  return sqrt(v);
}

template<typename Scalar>
Scalar diff_norm(const DiffVec3<Scalar>& v)
{
  return diff_safe_sqrt(diff_norm_sqr(v));
}

template<typename Scalar>
Scalar diff_clamp(const Scalar& v, double lo, double hi)
{
  if (v < lo)
    return Scalar(lo);
  if (v > hi)
    return Scalar(hi);
  return v;
}

template<typename Scalar>
Scalar diff_point_segment_distance_sqr(
  const DiffVec3<Scalar>& p,
  const DiffVec3<Scalar>& a,
  const DiffVec3<Scalar>& b)
{
  const DiffVec3<Scalar> ab = b - a;
  const Scalar denom = diff_norm_sqr(ab);
  if (denom <= 1e-30)
    return diff_norm_sqr(p - a);

  Scalar t = diff_dot(p - a, ab) / denom;
  t = diff_clamp(t, 0.0, 1.0);
  const DiffVec3<Scalar> closest = a + ab * t;
  return diff_norm_sqr(p - closest);
}

template<typename Scalar>
Scalar diff_point_triangle_distance_sqr(
  const DiffVec3<Scalar>& p,
  const DiffVec3<Scalar>& a,
  const DiffVec3<Scalar>& b,
  const DiffVec3<Scalar>& c)
{
  const DiffVec3<Scalar> ab = b - a;
  const DiffVec3<Scalar> ac = c - a;
  const DiffVec3<Scalar> tri_normal = diff_cross(ab, ac);
  if (diff_norm_sqr(tri_normal) <= 1e-30)
  {
    const Scalar d01 = diff_point_segment_distance_sqr(p, a, b);
    const Scalar d12 = diff_point_segment_distance_sqr(p, b, c);
    const Scalar d20 = diff_point_segment_distance_sqr(p, c, a);
    if (d01 <= d12 && d01 <= d20)
      return d01;
    if (d12 <= d20)
      return d12;
    return d20;
  }

  const DiffVec3<Scalar> ap = p - a;
  const Scalar d1 = diff_dot(ab, ap);
  const Scalar d2 = diff_dot(ac, ap);
  if (d1 <= 0.0 && d2 <= 0.0)
    return diff_norm_sqr(ap);

  const DiffVec3<Scalar> bp = p - b;
  const Scalar d3 = diff_dot(ab, bp);
  const Scalar d4 = diff_dot(ac, bp);
  if (d3 >= 0.0 && d4 <= d3)
    return diff_norm_sqr(bp);

  const Scalar vc = d1 * d4 - d3 * d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
  {
    const Scalar v = d1 / (d1 - d3);
    const DiffVec3<Scalar> closest = a + ab * v;
    return diff_norm_sqr(p - closest);
  }

  const DiffVec3<Scalar> cp = p - c;
  const Scalar d5 = diff_dot(ab, cp);
  const Scalar d6 = diff_dot(ac, cp);
  if (d6 >= 0.0 && d5 <= d6)
    return diff_norm_sqr(cp);

  const Scalar vb = d5 * d2 - d1 * d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
  {
    const Scalar w = d2 / (d2 - d6);
    const DiffVec3<Scalar> closest = a + ac * w;
    return diff_norm_sqr(p - closest);
  }

  const Scalar va = d3 * d6 - d5 * d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
  {
    const DiffVec3<Scalar> bc = c - b;
    const Scalar w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    const DiffVec3<Scalar> closest = b + bc * w;
    return diff_norm_sqr(p - closest);
  }

  const Scalar denom = Scalar(1.0) / (va + vb + vc);
  const Scalar v = vb * denom;
  const Scalar w = vc * denom;
  const DiffVec3<Scalar> closest = a + ab * v + ac * w;
  return diff_norm_sqr(p - closest);
}

template<typename Scalar>
Scalar diff_ipc_barrier(const Scalar& distance, double dhat)
{
  if (dhat <= 0.0 || distance >= dhat)
    return Scalar(0.0);

  using std::log;
  const double dmin = dhat * 1e-9;
  const Scalar d = distance <= dmin ? Scalar(dmin) : distance;
  return -diff_sqr(d - dhat) * log(d / dhat) / diff_sqr(dhat);
}

template<typename Scalar>
Scalar diff_triangle_quality(
  const DiffVec3<Scalar>& p0,
  const DiffVec3<Scalar>& p1,
  const DiffVec3<Scalar>& p2)
{
  const Scalar a = diff_norm(p1 - p0);
  const Scalar b = diff_norm(p2 - p1);
  const Scalar c = diff_norm(p0 - p2);
  const Scalar denom = a * a + b * b + c * c;
  if (denom <= 1e-30)
    return Scalar(0.0);

  const Scalar area = Scalar(0.5) * diff_norm(diff_cross(p1 - p0, p2 - p0));
  return Scalar(4.0 * std::sqrt(3.0)) * area / denom;
}

template<typename Scalar>
DiffVec3<Scalar> diff_triangle_normal(
  const DiffVec3<Scalar>& p0,
  const DiffVec3<Scalar>& p1,
  const DiffVec3<Scalar>& p2)
{
  const DiffVec3<Scalar> n = diff_cross(p1 - p0, p2 - p0);
  const Scalar inv_norm = Scalar(1.0) / diff_safe_sqrt(diff_norm_sqr(n));
  return n * inv_norm;
}

struct DiffBarrierTerm
{
  double x_weight = 0.0;
  Vec3d bias;
  Vec3d tri[3];
};

struct DiffNormalTerm
{
  Vec3d from;
  Vec3d to;
  Vec3d source_normal;
};

struct DiffUniformityTerm
{
  Vec3d neighbor;
  double target_length = 1.0;
};

struct DiffPositionFidelityTerm
{
  double x_weight = 0.0;
  Vec3d bias;
  Vec3d target;
};
}

CollapseStage::CollapseStage(
  SMeshT* original, SMeshT* cage, ParamCollapseStage* p,
  DFaceTree* original_tree, LightDFaceTree* remeshing_tree,
  FaceGrid* original_grid,
  double _original_diagonal_length)
  :om(original), rm(cage), param(p),
  ot(original_tree), lrt(remeshing_tree), og(original_grid),
  original_diagonal_length(_original_diagonal_length),
  avg_edge_length(0.0),
  avg_original_edge_length(0.0)
{}

void CollapseStage::update(size_t _candidate_points_size, bool _allow_negtive, double _max_distance_error)
{
  candidate_points_size = _candidate_points_size;
  allow_negtive = _allow_negtive;
  max_distance_error = _max_distance_error;
}

/// @brief generate random points around smooth target of edge.
std::vector<Vec3d> CollapseStage::generate_candidate_points_for_collapse(EdgeHandle e, EdgeCollapser& edge_collapser)
{
  // calculate radius, which is average length of adjacent edges.
  double radius = 0.0;
  HalfedgeHandle he = rm->halfedge_handle(e, 0);
  for (EdgeHandle ve : rm->ve_range(rm->to_vertex_handle(he)))
    radius += rm->data(ve).edge_length;
  for (EdgeHandle ve : rm->ve_range(rm->from_vertex_handle(he)))
    radius += rm->data(ve).edge_length;
  radius /= rm->valence(rm->to_vertex_handle(he)) + rm->valence(rm->from_vertex_handle(he));
  radius *= 0.2;
  radius = std::min(radius, original_diagonal_length * 0.01);

  // tangential relaxation target.
  auto& halfedges = edge_collapser.get_halfedges();

  Vec3d edge_midpoint = rm->calc_edge_midpoint(he);
  Vec3d new_vertex_normal, vertex_target;
  edge_collapser.predict_tangential_weighted_smooth_target(edge_midpoint, new_vertex_normal, vertex_target);

  // local coordinate system on target point.
  Vec3d local_axis_x, local_axis_y;
  make_coordinate_system(new_vertex_normal, local_axis_x, local_axis_y);

  std::vector<Vec3d> points; points.reserve(candidate_points_size);
  // generate regular points
  points.push_back(vertex_target);
  points.push_back(rm->calc_edge_midpoint(he));
  points.push_back(rm->point(rm->to_vertex_handle(he)));
  points.push_back(rm->point(rm->from_vertex_handle(he)));

  // generate random points
  while (points.size() < candidate_points_size)
  {
    double height = (double)rand() / (double)RAND_MAX;
    double beta = (double)rand() / (double)RAND_MAX;
    double len = (double)rand() / (double)RAND_MAX;
    height = (height - 0.5) * original_diagonal_length * 0.005;
    beta = beta * 2 * M_PI;

    Vec3d new_point = vertex_target + (local_axis_x * cos(beta) + local_axis_y * sin(beta)) * len * radius + height * new_vertex_normal;
    points.push_back(new_point);
  }
  points.resize(candidate_points_size);

  return points;
}

bool CollapseStage::is_optimization_collapse_placement_method() const
{
  return param->collapsePlacementMethod == "optimization" ||
    param->collapsePlacementMethod == "energy" ||
    param->collapsePlacementMethod == "newton";
}

bool CollapseStage::is_phase2_adaptive_strategy() const
{
  return param->phase2PlacementStrategy == "adaptive" || param->phase2PlacementStrategy.empty();
}

bool CollapseStage::is_phase2_linear_only_strategy() const
{
  return param->phase2PlacementStrategy == "linear_only" || param->phase2PlacementStrategy == "linear-only" ||
    param->phase2PlacementStrategy == "qem_only" || param->phase2PlacementStrategy == "qem-only";
}

bool CollapseStage::is_phase2_final_newton_strategy() const
{
  return param->phase2PlacementStrategy == "final_newton" || param->phase2PlacementStrategy == "final-newton";
}

bool CollapseStage::is_phase2_newton_only_strategy() const
{
  return param->phase2PlacementStrategy == "newton_only" || param->phase2PlacementStrategy == "newton-only";
}

bool CollapseStage::is_trust_region_solver_mode() const
{
  return param->newtonSolverMode == "trust_region" || param->newtonSolverMode == "trust-region";
}

bool CollapseStage::is_exact_reject_robustness_mode() const
{
  return param->robustnessMode == "exact_reject" || param->robustnessMode == "exact-reject";
}

bool CollapseStage::is_exact_backtracking_robustness_mode() const
{
  return param->robustnessMode == "exact_backtracking" || param->robustnessMode == "exact-backtracking";
}

bool CollapseStage::is_ipc_line_search_robustness_mode() const
{
  return param->robustnessMode == "ipc_line_search" || param->robustnessMode == "ipc-line-search" ||
    param->robustnessMode == "ccd" || param->robustnessMode == "ipc";
}

// Build the local data needed to evaluate candidate positions for this collapse.
// No topology is changed here; the context describes the mesh patch that would
// surround the new vertex after collapsing edge eh.
CollapseStage::Phase2PlacementContext CollapseStage::make_phase2_placement_context(
  EdgeHandle eh, EdgeCollapser& edge_collapser) const
{
  Phase2PlacementContext ctx;
  ctx.edge = eh;
  ctx.halfedges = edge_collapser.get_halfedges();
  ctx.global_target_length = avg_edge_length > 0.0 ? avg_edge_length : original_diagonal_length * 0.01;

  HalfedgeHandle heh = rm->halfedge_handle(eh, 0);
  ctx.midpoint = rm->calc_edge_midpoint(heh);
  Vec3d smooth_normal;
  edge_collapser.predict_tangential_weighted_smooth_target(ctx.midpoint, smooth_normal, ctx.start_point);
  if (!finite_vec(ctx.start_point))
    ctx.start_point = ctx.midpoint;

  double local_scale = rm->data(eh).edge_length;
  std::set<VertexHandle> neighbor_vertices;
  for (HalfedgeHandle h : ctx.halfedges)
  {
    const VertexHandle from_v = rm->from_vertex_handle(h);
    const VertexHandle to_v = rm->to_vertex_handle(h);
    const Vec3d& from = rm->point(from_v);
    const Vec3d& to = rm->point(to_v);
    ctx.fan_edges.push_back({ from, to });
    neighbor_vertices.insert(from_v);
    neighbor_vertices.insert(to_v);
    local_scale += (from - to).length();
    if (rm->face_handle(h).is_valid())
      ctx.ignored_faces.insert(rm->face_handle(h).idx());
  }
  if (!rm->is_boundary(heh))
    ctx.ignored_faces.insert(rm->face_handle(heh).idx());
  const HalfedgeHandle heh_opp = rm->opposite_halfedge_handle(heh);
  if (!rm->is_boundary(heh_opp))
    ctx.ignored_faces.insert(rm->face_handle(heh_opp).idx());

  for (VertexHandle vh : neighbor_vertices)
    ctx.neighbor_points.push_back(rm->point(vh));

  // Average local edge length used only for scale normalization and step sizing.
  const size_t scale_count = ctx.fan_edges.empty() ? 1 : ctx.fan_edges.size() + 1;
  ctx.local_scale = std::max(local_scale / static_cast<double>(scale_count), original_diagonal_length * 1e-6);

  std::set<FaceHandle> qem_faces;
  one_ring_faces_around_edge(rm, heh, qem_faces);
  for (FaceHandle fh : qem_faces)
  {
    Vec3d pts[3];
    collect_face_points(rm, fh, pts);
    const Vec3d normal = triangle_normal(pts[0], pts[1], pts[2]);
    const double area = std::max(triangle_area_from_points(pts[0], pts[1], pts[2]), 1e-16);
    ctx.qem_planes.push_back({ normal, -(normal | pts[0]), area });
  }

  std::vector<Vec3d> adjacent_normals;
  for (FaceHandle fh : faces_adjacent_to_edge(rm, eh))
  {
    Vec3d pts[3];
    collect_face_points(rm, fh, pts);
    adjacent_normals.push_back(triangle_normal(pts[0], pts[1], pts[2]));
  }
  if (adjacent_normals.size() == 2)
  {
    const double dot_value = clamp_value(adjacent_normals[0] | adjacent_normals[1], -1.0, 1.0);
    ctx.edge_curvature = std::acos(dot_value) / std::max(rm->data(eh).edge_length, ctx.local_scale * 1e-3);
  }

  return ctx;
}

// Build the same x-dependent local fan used by collapse placement, but for
// moving an existing vertex after all collapses are finished.
CollapseStage::Phase2PlacementContext CollapseStage::make_phase2_vertex_relocation_context(
  VertexHandle vh, VertexRelocater& vertex_relocater) const
{
  Phase2PlacementContext ctx;
  ctx.halfedges = vertex_relocater.get_halfedges();
  ctx.global_target_length = avg_edge_length > 0.0 ? avg_edge_length : original_diagonal_length * 0.01;
  ctx.midpoint = rm->point(vh);
  ctx.start_point = vertex_relocater.find_weighted_tangential_smooth_target();
  if (!finite_vec(ctx.start_point))
    ctx.start_point = ctx.midpoint;

  double local_scale = 0.0;
  size_t local_scale_count = 0;
  for (EdgeHandle eh : rm->ve_range(vh))
  {
    local_scale += rm->data(eh).edge_length;
    local_scale_count++;
  }

  for (VertexHandle vv : rm->vv_range(vh))
    ctx.neighbor_points.push_back(rm->point(vv));

  for (HalfedgeHandle h : ctx.halfedges)
  {
    const FaceHandle fh = rm->face_handle(h);
    if (!fh.is_valid())
      continue;

    const VertexHandle from_v = rm->from_vertex_handle(h);
    const VertexHandle to_v = rm->to_vertex_handle(h);
    const Vec3d& from = rm->point(from_v);
    const Vec3d& to = rm->point(to_v);
    ctx.fan_edges.push_back({ from, to });
    ctx.ignored_faces.insert(fh.idx());

    Vec3d pts[3];
    collect_face_points(rm, fh, pts);
    const Vec3d normal = triangle_normal(pts[0], pts[1], pts[2]);
    const double area = std::max(triangle_area_from_points(pts[0], pts[1], pts[2]), 1e-16);
    ctx.qem_planes.push_back({ normal, -(normal | pts[0]), area });
  }

  if (local_scale_count == 0)
    local_scale = original_diagonal_length * 0.01;
  else
    local_scale /= static_cast<double>(local_scale_count);
  ctx.local_scale = std::max(local_scale, original_diagonal_length * 1e-6);

  return ctx;
}

double CollapseStage::ipc_barrier(double distance, double dhat) const
{
  if (dhat <= 0.0 || distance >= dhat)
    return 0.0;

  const double d = std::max(distance, dhat * 1e-9);
  return -sqr(d - dhat) * std::log(d / dhat) / sqr(dhat);
}

double CollapseStage::evaluate_qem_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const
{
  double energy = 0.0;
  for (const QEMPlane& plane : ctx.qem_planes)
  {
    const double value = (plane.normal | x) + plane.offset;
    energy += plane.weight * sqr(value / ctx.local_scale);
  }

  if (param->curvatureMode == "weighted_qem" || param->curvatureMode == "weighted-qem")
  {
    const double scaled_curvature = ctx.edge_curvature * ctx.local_scale;
    energy *= 1.0 + param->curvatureWeight * sqr(scaled_curvature);
  }

  return energy;
}

double CollapseStage::evaluate_original_barrier_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const
{
  if (param->originalBarrierWeight <= 0.0)
    return 0.0;

  const double dhat = std::max(
    original_diagonal_length * param->barrierActivationDistanceFactor,
    ctx.local_scale * 1e-3);
  double energy = 0.0;

  for (const CollapseFanEdge& fan_edge : ctx.fan_edges)
  {
    const Vec3d samples[6] = {
      x,
      fan_edge.from,
      fan_edge.to,
      (fan_edge.from + fan_edge.to) * 0.5,
      (fan_edge.from + x) * 0.5,
      (fan_edge.to + x) * 0.5
    };

    for (const Vec3d& sample : samples)
    {
#ifdef USE_TREE_SEARCH
      const double distance = ot ? ot->closest_distance(sample).first : dhat;
#else
      const double distance = og ? og->closest_point(sample).first.second : dhat;
#endif
      energy += ipc_barrier(distance, dhat);
    }
  }
  return energy;
}

double CollapseStage::evaluate_self_barrier_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const
{
  if (param->selfBarrierWeight <= 0.0 || !lrt)
    return 0.0;

  const double dhat = std::max(
    original_diagonal_length * param->barrierActivationDistanceFactor,
    ctx.local_scale * 1e-3);
  double energy = 0.0;

  for (const CollapseFanEdge& fan_edge : ctx.fan_edges)
  {
    const Vec3d samples[4] = {
      (fan_edge.from + fan_edge.to + x) / 3.0,
      (fan_edge.from + x) * 0.5,
      (fan_edge.to + x) * 0.5,
      x
    };

    for (const Vec3d& sample : samples)
    {
      double min_distance_sqr = DBL_MAX;
      for (size_t prim_idx = 0; prim_idx < lrt->m_primitives.size(); prim_idx++)
      {
        if (prim_idx < lrt->m_primitive_deleted.size() && lrt->m_primitive_deleted[prim_idx])
          continue;
        const auto& prim = lrt->m_primitives[prim_idx];
        if (ctx.ignored_faces.count(static_cast<int>(prim.tri.index)))
          continue;
        HeavyTriangle tri(prim.tri.ver0, prim.tri.ver1, prim.tri.ver2);
        Vec3d closest;
        const double distance_sqr = tri.closest_point(sample, closest);
        if (distance_sqr < min_distance_sqr)
          min_distance_sqr = distance_sqr;
      }
      if (min_distance_sqr != DBL_MAX)
        energy += ipc_barrier(std::sqrt(min_distance_sqr), dhat);
    }
  }
  return energy;
}

bool CollapseStage::closest_original_point(const Vec3d& p, Vec3d& closest) const
{
#ifdef USE_TREE_SEARCH
  if (!ot)
    return false;
  const auto result = ot->closest_point(p);
  if (result.second < 0)
    return false;
  closest = result.first;
#else
  if (!og)
    return false;
  const auto result = og->closest_point(p);
  if (result.second < 0)
    return false;
  closest = result.first.first;
#endif
  return finite_vec(closest);
}

double CollapseStage::evaluate_position_fidelity_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const
{
  if (param->positionFidelityWeight <= 0.0 || !finite_vec(x))
    return 0.0;

  const double inv_scale_sqr = 1.0 / std::max(sqr(ctx.local_scale), 1e-24);
  double energy = 0.0;
  size_t count = 0;

  const auto append_sample = [&](double x_weight, const Vec3d& bias)
  {
    const Vec3d sample = bias + x * x_weight;
    Vec3d target;
    if (!closest_original_point(sample, target))
      return;
    energy += (sample - target).squaredNorm() * inv_scale_sqr;
    count++;
  };

  append_sample(1.0, Vec3d(0.0, 0.0, 0.0));
  for (const CollapseFanEdge& fan_edge : ctx.fan_edges)
  {
    append_sample(1.0 / 3.0, (fan_edge.from + fan_edge.to) / 3.0);
    append_sample(0.5, fan_edge.from * 0.5);
    append_sample(0.5, fan_edge.to * 0.5);
  }

  return count == 0 ? 0.0 : energy / static_cast<double>(count);
}

double CollapseStage::evaluate_curvature_normal_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const
{
  if (param->curvatureMode != "normal_matching" && param->curvatureMode != "normal-matching")
    return 0.0;

  double energy = 0.0;
  for (const CollapseFanEdge& fan_edge : ctx.fan_edges)
  {
    const Vec3d normal = triangle_normal(fan_edge.from, fan_edge.to, x);
    const Vec3d centroid = (fan_edge.from + fan_edge.to + x) / 3.0;

    FaceHandle closest_face;
#ifdef USE_TREE_SEARCH
    if (!ot)
      continue;
    closest_face = FaceHandle(ot->closest_point(centroid).second);
#else
    if (!og)
      continue;
    closest_face = FaceHandle(og->closest_point(centroid).second);
#endif
    if (!closest_face.is_valid() || closest_face.idx() < 0 ||
      closest_face.idx() >= static_cast<int>(om->n_faces()))
      continue;

    Vec3d pts[3];
    collect_face_points(om, closest_face, pts);
    const Vec3d source_normal = triangle_normal(pts[0], pts[1], pts[2]);
    const double alignment = std::abs(clamp_value(normal | source_normal, -1.0, 1.0));
    energy += sqr(1.0 - alignment);
  }
  return energy;
}

double CollapseStage::evaluate_triangle_quality_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const
{
  if (param->triangleQualityWeight <= 0.0)
    return 0.0;

  double energy = 0.0;
  for (const CollapseFanEdge& fan_edge : ctx.fan_edges)
  {
    const double quality = clamp_value(triangle_quality_from_points(fan_edge.from, fan_edge.to, x), 0.0, 1.0);
    energy += sqr(1.0 - quality);
  }
  return energy;
}

double CollapseStage::source_size_at(const Vec3d& p) const
{
  FaceHandle closest_face;
#ifdef USE_TREE_SEARCH
  if (!ot)
    return avg_edge_length;
  closest_face = FaceHandle(ot->closest_point(p).second);
#else
  if (!og)
    return avg_edge_length;
  closest_face = FaceHandle(og->closest_point(p).second);
#endif

  if (!closest_face.is_valid() || closest_face.idx() < 0 ||
    closest_face.idx() >= static_cast<int>(om->n_faces()))
    return avg_edge_length;

  const double source_size = face_average_edge_length(om, closest_face);
  return source_size > 0.0 ? source_size : avg_edge_length;
}

double CollapseStage::evaluate_uniformity_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const
{
  if (param->uniformityWeight <= 0.0 || param->uniformityMode == "none")
    return 0.0;

  double energy = 0.0;
  for (const Vec3d& neighbor : ctx.neighbor_points)
  {
    double target_length = source_uniformity_target_length(ctx, (x + neighbor) * 0.5);
    target_length = std::max(target_length, ctx.local_scale * 1e-3);
    energy += sqr(((x - neighbor).length() - target_length) / target_length);
  }
  return energy;
}

double CollapseStage::evaluate_newton_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const
{
  if (!finite_vec(x))
    return DBL_MAX;

  double energy = 0.0;
  energy += param->qemWeight * evaluate_qem_energy(ctx, x);
  energy += param->originalBarrierWeight * evaluate_original_barrier_energy(ctx, x);
  energy += param->selfBarrierWeight * evaluate_self_barrier_energy(ctx, x);
  energy += param->positionFidelityWeight * evaluate_position_fidelity_energy(ctx, x);
  energy += param->curvatureWeight * evaluate_curvature_normal_energy(ctx, x);
  energy += param->triangleQualityWeight * evaluate_triangle_quality_energy(ctx, x);
  energy += param->uniformityWeight * evaluate_uniformity_energy(ctx, x);
  return std::isfinite(energy) ? energy : DBL_MAX;
}

CollapseStage::NewtonDerivatives CollapseStage::finite_difference_newton_derivatives(
  const Phase2PlacementContext& ctx, const Vec3d& x) const
{
  NewtonDerivatives deriv;
  const double h = std::max(ctx.local_scale * param->newtonFiniteDiffScale, original_diagonal_length * 1e-9);
  deriv.energy = evaluate_newton_energy(ctx, x);

  Vec3d basis[3] = { Vec3d(h, 0.0, 0.0), Vec3d(0.0, h, 0.0), Vec3d(0.0, 0.0, h) };
  double e_plus[3];
  double e_minus[3];

  for (int i = 0; i < 3; i++)
  {
    e_plus[i] = evaluate_newton_energy(ctx, x + basis[i]);
    e_minus[i] = evaluate_newton_energy(ctx, x - basis[i]);
    deriv.gradient[i] = (e_plus[i] - e_minus[i]) / (2.0 * h);
    deriv.hessian[i][i] = (e_plus[i] - 2.0 * deriv.energy + e_minus[i]) / (h * h);
  }

  for (int i = 0; i < 3; i++)
  {
    for (int j = i + 1; j < 3; j++)
    {
      const double e_pp = evaluate_newton_energy(ctx, x + basis[i] + basis[j]);
      const double e_pm = evaluate_newton_energy(ctx, x + basis[i] - basis[j]);
      const double e_mp = evaluate_newton_energy(ctx, x - basis[i] + basis[j]);
      const double e_mm = evaluate_newton_energy(ctx, x - basis[i] - basis[j]);
      const double hij = (e_pp - e_pm - e_mp + e_mm) / (4.0 * h * h);
      deriv.hessian[i][j] = hij;
      deriv.hessian[j][i] = hij;
    }
  }

  return deriv;
}

CollapseStage::NewtonDerivatives CollapseStage::approximate_newton_derivatives(
  const Phase2PlacementContext& ctx, const Vec3d& x) const
{
  return autodiff_newton_derivatives(ctx, x);
}

CollapseStage::NewtonDerivatives CollapseStage::autodiff_newton_derivatives(
  const Phase2PlacementContext& ctx, const Vec3d& x) const
{
  NewtonDerivatives deriv;
  deriv.energy = evaluate_newton_energy(ctx, x);
  if (!std::isfinite(deriv.energy))
    return deriv;

  const double dhat = std::max(
    original_diagonal_length * param->barrierActivationDistanceFactor,
    ctx.local_scale * 1e-3);

  std::vector<DiffBarrierTerm> original_barrier_terms;
  std::vector<DiffBarrierTerm> self_barrier_terms;
  std::vector<DiffPositionFidelityTerm> position_fidelity_terms;
  std::vector<DiffNormalTerm> normal_terms;
  std::vector<DiffUniformityTerm> uniformity_terms;

  const auto append_original_barrier_term =
    [&](double x_weight, const Vec3d& bias)
  {
    if (param->originalBarrierWeight <= 0.0)
      return;

    const Vec3d sample = bias + x_weight * x;
    FaceHandle closest_face;
#ifdef USE_TREE_SEARCH
    if (!ot)
      return;
    closest_face = FaceHandle(ot->closest_point(sample).second);
#else
    if (!og)
      return;
    closest_face = FaceHandle(og->closest_point(sample).second);
#endif
    if (!closest_face.is_valid() || closest_face.idx() < 0 ||
      closest_face.idx() >= static_cast<int>(om->n_faces()))
      return;

    DiffBarrierTerm term;
    term.x_weight = x_weight;
    term.bias = bias;
    collect_face_points(om, closest_face, term.tri);
    original_barrier_terms.push_back(term);
  };

  const auto append_self_barrier_term =
    [&](double x_weight, const Vec3d& bias)
  {
    if (param->selfBarrierWeight <= 0.0 || !lrt)
      return;

    const Vec3d sample = bias + x_weight * x;
    double min_distance_sqr = DBL_MAX;
    DiffBarrierTerm best_term;
    bool found = false;
    for (size_t prim_idx = 0; prim_idx < lrt->m_primitives.size(); prim_idx++)
    {
      if (prim_idx < lrt->m_primitive_deleted.size() && lrt->m_primitive_deleted[prim_idx])
        continue;
      const auto& prim = lrt->m_primitives[prim_idx];
      if (ctx.ignored_faces.count(static_cast<int>(prim.tri.index)))
        continue;

      HeavyTriangle tri(prim.tri.ver0, prim.tri.ver1, prim.tri.ver2);
      Vec3d closest;
      const double distance_sqr = tri.closest_point(sample, closest);
      if (distance_sqr < min_distance_sqr)
      {
        min_distance_sqr = distance_sqr;
        best_term.x_weight = x_weight;
        best_term.bias = bias;
        best_term.tri[0] = prim.tri.ver0;
        best_term.tri[1] = prim.tri.ver1;
        best_term.tri[2] = prim.tri.ver2;
        found = true;
      }
    }

    if (found)
      self_barrier_terms.push_back(best_term);
  };

  const auto append_position_fidelity_term =
    [&](double x_weight, const Vec3d& bias)
  {
    if (param->positionFidelityWeight <= 0.0)
      return;

    const Vec3d sample = bias + x_weight * x;
    Vec3d target;
    if (!closest_original_point(sample, target))
      return;

    DiffPositionFidelityTerm term;
    term.x_weight = x_weight;
    term.bias = bias;
    term.target = target;
    position_fidelity_terms.push_back(term);
  };

  append_position_fidelity_term(1.0, Vec3d(0.0, 0.0, 0.0));
  for (const CollapseFanEdge& fan_edge : ctx.fan_edges)
  {
    append_original_barrier_term(1.0, Vec3d(0.0, 0.0, 0.0));
    append_original_barrier_term(0.0, fan_edge.from);
    append_original_barrier_term(0.0, fan_edge.to);
    append_original_barrier_term(0.0, (fan_edge.from + fan_edge.to) * 0.5);
    append_original_barrier_term(0.5, fan_edge.from * 0.5);
    append_original_barrier_term(0.5, fan_edge.to * 0.5);

    append_self_barrier_term(1.0 / 3.0, (fan_edge.from + fan_edge.to) / 3.0);
    append_self_barrier_term(0.5, fan_edge.from * 0.5);
    append_self_barrier_term(0.5, fan_edge.to * 0.5);
    append_self_barrier_term(1.0, Vec3d(0.0, 0.0, 0.0));

    append_position_fidelity_term(1.0 / 3.0, (fan_edge.from + fan_edge.to) / 3.0);
    append_position_fidelity_term(0.5, fan_edge.from * 0.5);
    append_position_fidelity_term(0.5, fan_edge.to * 0.5);
  }

  if (param->curvatureWeight > 0.0 &&
    (param->curvatureMode == "normal_matching" || param->curvatureMode == "normal-matching"))
  {
    for (const CollapseFanEdge& fan_edge : ctx.fan_edges)
    {
      const Vec3d centroid = (fan_edge.from + fan_edge.to + x) / 3.0;
      FaceHandle closest_face;
#ifdef USE_TREE_SEARCH
      if (!ot)
        continue;
      closest_face = FaceHandle(ot->closest_point(centroid).second);
#else
      if (!og)
        continue;
      closest_face = FaceHandle(og->closest_point(centroid).second);
#endif
      if (!closest_face.is_valid() || closest_face.idx() < 0 ||
        closest_face.idx() >= static_cast<int>(om->n_faces()))
        continue;

      Vec3d pts[3];
      collect_face_points(om, closest_face, pts);
      DiffNormalTerm term;
      term.from = fan_edge.from;
      term.to = fan_edge.to;
      term.source_normal = triangle_normal(pts[0], pts[1], pts[2]);
      normal_terms.push_back(term);
    }
  }

  if (param->uniformityWeight > 0.0 && param->uniformityMode != "none")
  {
    for (const Vec3d& neighbor : ctx.neighbor_points)
    {
      DiffUniformityTerm term;
      term.neighbor = neighbor;
      term.target_length = source_uniformity_target_length(ctx, (x + neighbor) * 0.5);
      uniformity_terms.push_back(term);
    }
  }

  using ADScalar = autodiff::var;
  Eigen::Matrix<ADScalar, 3, 1> ad_x;
  ad_x[0] = x.x();
  ad_x[1] = x.y();
  ad_x[2] = x.z();
  const DiffVec3<ADScalar> ad_point(ad_x[0], ad_x[1], ad_x[2]);
  ADScalar ad_energy = 0.0;

  const double inv_scale_sqr = 1.0 / std::max(sqr(ctx.local_scale), 1e-24);
  double curvature_multiplier = 1.0;
  if (param->curvatureMode == "weighted_qem" || param->curvatureMode == "weighted-qem")
  {
    const double scaled_curvature = ctx.edge_curvature * ctx.local_scale;
    curvature_multiplier += param->curvatureWeight * sqr(scaled_curvature);
  }

  for (const QEMPlane& plane : ctx.qem_planes)
  {
    const ADScalar value =
      plane.normal.x() * ad_point.x +
      plane.normal.y() * ad_point.y +
      plane.normal.z() * ad_point.z +
      plane.offset;
    const double weight =
      param->qemWeight * curvature_multiplier * plane.weight * inv_scale_sqr;
    ad_energy += weight * diff_sqr(value);
  }

  for (const DiffBarrierTerm& term : original_barrier_terms)
  {
    const DiffVec3<ADScalar> sample = diff_vec<ADScalar>(term.bias) + ad_point * term.x_weight;
    const ADScalar distance = diff_safe_sqrt(diff_point_triangle_distance_sqr(
      sample,
      diff_vec<ADScalar>(term.tri[0]),
      diff_vec<ADScalar>(term.tri[1]),
      diff_vec<ADScalar>(term.tri[2])));
    ad_energy += param->originalBarrierWeight * diff_ipc_barrier(distance, dhat);
  }

  for (const DiffBarrierTerm& term : self_barrier_terms)
  {
    const DiffVec3<ADScalar> sample = diff_vec<ADScalar>(term.bias) + ad_point * term.x_weight;
    const ADScalar distance = diff_safe_sqrt(diff_point_triangle_distance_sqr(
      sample,
      diff_vec<ADScalar>(term.tri[0]),
      diff_vec<ADScalar>(term.tri[1]),
      diff_vec<ADScalar>(term.tri[2])));
    ad_energy += param->selfBarrierWeight * diff_ipc_barrier(distance, dhat);
  }

  if (!position_fidelity_terms.empty())
  {
    const double position_weight =
      param->positionFidelityWeight /
      static_cast<double>(position_fidelity_terms.size()) *
      inv_scale_sqr;
    for (const DiffPositionFidelityTerm& term : position_fidelity_terms)
    {
      const DiffVec3<ADScalar> sample = diff_vec<ADScalar>(term.bias) + ad_point * term.x_weight;
      const DiffVec3<ADScalar> residual = sample - diff_vec<ADScalar>(term.target);
      ad_energy += position_weight * diff_norm_sqr(residual);
    }
  }

  for (const DiffNormalTerm& term : normal_terms)
  {
    const DiffVec3<ADScalar> normal = diff_triangle_normal(
      diff_vec<ADScalar>(term.from),
      diff_vec<ADScalar>(term.to),
      ad_point);
    ADScalar alignment = diff_dot(normal, diff_vec<ADScalar>(term.source_normal));
    alignment = diff_clamp(alignment, -1.0, 1.0);
    using std::abs;
    const ADScalar abs_alignment = abs(alignment);
    ad_energy += param->curvatureWeight * diff_sqr(ADScalar(1.0) - abs_alignment);
  }

  for (const CollapseFanEdge& fan_edge : ctx.fan_edges)
  {
    const ADScalar quality = diff_clamp(diff_triangle_quality(
      diff_vec<ADScalar>(fan_edge.from),
      diff_vec<ADScalar>(fan_edge.to),
      ad_point), 0.0, 1.0);
    ad_energy += param->triangleQualityWeight * diff_sqr(ADScalar(1.0) - quality);
  }

  for (const DiffUniformityTerm& term : uniformity_terms)
  {
    const double target_length = std::max(term.target_length, ctx.local_scale * 1e-3);
    const ADScalar length = diff_norm(ad_point - diff_vec<ADScalar>(term.neighbor));
    ad_energy += param->uniformityWeight *
      diff_sqr((length - target_length) / target_length);
  }

  Eigen::Vector3d gradient = Eigen::Vector3d::Zero();
  Eigen::Matrix3d hessian = autodiff::hessian(ad_energy, ad_x, gradient);
  hessian = 0.5 * (hessian + hessian.transpose());

  const double damping = std::max(
    1e-10,
    1e-6 * std::max(1.0, std::abs(deriv.energy)) / std::max(sqr(ctx.local_scale), 1e-12));
  hessian += Eigen::Matrix3d::Identity() * damping;

  if (!std::isfinite(gradient.x()) || !std::isfinite(gradient.y()) || !std::isfinite(gradient.z()) ||
    !hessian.allFinite())
    return deriv;

  deriv.gradient = from_eigen(gradient);
  for (int r = 0; r < 3; r++)
    for (int c = 0; c < 3; c++)
      deriv.hessian[r][c] = hessian(r, c);

  return deriv;
}

bool CollapseStage::collapse_target_valid(EdgeCollapser& edge_collapser, const Vec3d& x) const
{
  return finite_vec(x) && edge_collapser.target_point_is_valid(x, nullptr);
}

bool CollapseStage::sampled_path_valid(EdgeCollapser& edge_collapser, const Vec3d& from, const Vec3d& to) const
{
  const size_t samples = std::max<size_t>(param->lineSearchCcdSamples, 2);
  for (size_t i = 1; i <= samples; i++)
  {
    const double alpha = static_cast<double>(i) / static_cast<double>(samples);
    if (!collapse_target_valid(edge_collapser, from + alpha * (to - from)))
      return false;
  }
  return true;
}

// Cheap score for queue initialization and fallback candidate selection.
// It includes QEM plus inexpensive quality/uniformity/position-fit terms, but
// skips the full Newton solve unless the placement policy below decides
// refinement is needed.
double CollapseStage::evaluate_phase2_proxy_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const
{
  if (!finite_vec(x))
    return DBL_MAX;
  double energy = param->qemWeight * evaluate_qem_energy(ctx, x);

  if (param->positionFidelityWeight > 0.0)
    energy += param->positionFidelityWeight * evaluate_position_fidelity_energy(ctx, x);

  if (param->triangleQualityWeight > 0.0)
  {
    const double fan_count = static_cast<double>(std::max<size_t>(ctx.fan_edges.size(), 1));
    energy += param->triangleQualityWeight * evaluate_triangle_quality_energy(ctx, x) / fan_count;

    const double min_quality = calc_phase2_fan_min_quality(ctx, x);
    if (param->phase2NewtonQualityThreshold > 0.0 &&
      min_quality < param->phase2NewtonQualityThreshold)
    {
      const double normalized_gap =
        (param->phase2NewtonQualityThreshold - min_quality) /
        std::max(param->phase2NewtonQualityThreshold, 1e-8);
      energy += 10.0 * param->triangleQualityWeight * sqr(normalized_gap);
    }
  }

  if (param->uniformityWeight > 0.0 && param->uniformityMode != "none")
  {
    const double neighbor_count = static_cast<double>(std::max<size_t>(ctx.neighbor_points.size(), 1));
    energy += param->uniformityWeight * evaluate_uniformity_energy(ctx, x) / neighbor_count;
  }

  return std::isfinite(energy) ? energy : DBL_MAX;
}

double CollapseStage::evaluate_phase2_refinement_residual(const Phase2PlacementContext& ctx, const Vec3d& x) const
{
  if (!finite_vec(x))
    return DBL_MAX;

  double residual = 0.0;
  if (param->curvatureWeight > 0.0 &&
    (param->curvatureMode == "normal_matching" || param->curvatureMode == "normal-matching"))
  {
    const double denom = std::max<size_t>(ctx.fan_edges.size(), 1);
    residual += param->curvatureWeight * evaluate_curvature_normal_energy(ctx, x) /
      static_cast<double>(denom);
  }

  if (param->uniformityWeight > 0.0 && param->uniformityMode != "none")
  {
    const double denom = std::max<size_t>(ctx.neighbor_points.size(), 1);
    residual += param->uniformityWeight * evaluate_uniformity_energy(ctx, x) /
      static_cast<double>(denom);
  }

  if (param->positionFidelityWeight > 0.0)
    residual += param->positionFidelityWeight * evaluate_position_fidelity_energy(ctx, x);

  return std::isfinite(residual) ? residual : DBL_MAX;
}

double CollapseStage::calc_phase2_fan_min_quality(const Phase2PlacementContext& ctx, const Vec3d& x) const
{
  if (!finite_vec(x))
    return 0.0;
  double min_quality = DBL_MAX;
  for (const CollapseFanEdge& fan_edge : ctx.fan_edges)
  {
    const double quality = clamp_value(triangle_quality_from_points(fan_edge.from, fan_edge.to, x), 0.0, 1.0);
    min_quality = std::min(min_quality, quality);
  }
  return min_quality == DBL_MAX ? 1.0 : min_quality;
}

double CollapseStage::source_uniformity_target_length(const Phase2PlacementContext& ctx, const Vec3d& p) const
{
  double target_length = ctx.global_target_length;
  if (param->uniformityMode == "source")
  {
    const double source_length = source_size_at(p);
    const double source_average = avg_original_edge_length > 0.0 ? avg_original_edge_length : source_length;
    if (source_average > 0.0 && std::isfinite(source_length))
      target_length *= source_length / source_average;
  }
  return std::max(target_length, ctx.local_scale * 1e-3);
}

// Hard validity test for a proposed new vertex position. This is the exact
// rejection layer used before collapse, independent of soft energy values.
bool CollapseStage::phase2_placement_satisfies_hard_constraints(
  const Phase2PlacementContext& ctx, EdgeCollapser& edge_collapser, const Vec3d& x) const
{
  if (!collapse_target_valid(edge_collapser, x))
    return false;
  if (is_ipc_line_search_robustness_mode())
    return sampled_path_valid(edge_collapser, ctx.midpoint, x);
  return true;
}

// Decide whether the cheap QEM/proxy placement is enough, or whether this edge
// needs the slower nonlinear 3D Newton placement refinement.
bool CollapseStage::should_refine_phase2_placement_with_newton(
  const Phase2PlacementContext& ctx, const Vec3d& x,
  size_t remaining_vertices, size_t target_vertices_num,
  double min_quality, double nonlinear_residual) const
{
  if (is_phase2_newton_only_strategy())
    return finite_vec(x);
  if (is_phase2_linear_only_strategy() || is_phase2_final_newton_strategy())
    return false;

  const bool quality_bad =
    param->phase2NewtonQualityThreshold > 0.0 &&
    min_quality < param->phase2NewtonQualityThreshold;

  bool residual_bad = false;
  if (param->phase2NewtonResidualThreshold > 0.0 && std::isfinite(nonlinear_residual) &&
    nonlinear_residual > param->phase2NewtonResidualThreshold)
  {
    double reference_residual = evaluate_phase2_refinement_residual(ctx, ctx.start_point);
    if (!std::isfinite(reference_residual))
      reference_residual = 0.0;
    const double growth = std::max(1.0, param->phase2NewtonResidualGrowth);
    residual_bad = nonlinear_residual > growth * reference_residual + 1e-12;
  }

  const bool final_refine =
    param->phase2NewtonFinalRefineCollapses > 0 &&
    remaining_vertices > target_vertices_num &&
    remaining_vertices - target_vertices_num <= param->phase2NewtonFinalRefineCollapses;

  return finite_vec(x) && (quality_bad || residual_bad || final_refine);
}

// Select the new vertex position for one edge collapse before the topology
// change is committed. The policy is controlled by phase2PlacementStrategy:
// adaptive, linear_only, final_newton, or newton_only.
bool CollapseStage::choose_phase2_collapse_placement(
  EdgeHandle eh, EdgeCollapser& edge_collapser, const Vec3d& qem_point,
  size_t remaining_vertices, size_t target_vertices_num,
  Phase2PlacementDecision& decision, double& newton_seconds)
{
  decision = Phase2PlacementDecision();
  const Phase2PlacementContext ctx = make_phase2_placement_context(eh, edge_collapser);
  if (ctx.qem_planes.empty())
    return false;

  std::vector<Vec3d> candidates;
  std::vector<bool> is_fallback;
  const double duplicate_tol = std::max(ctx.local_scale * 1e-8, original_diagonal_length * 1e-12);
  const auto append_candidate = [&](const Vec3d& p, bool fallback)
  {
    if (!finite_vec(p))
      return;
    for (const Vec3d& existing : candidates)
    {
      if ((existing - p).length() <= duplicate_tol)
        return;
    }
    candidates.push_back(p);
    is_fallback.push_back(fallback);
  };

  append_candidate(qem_point, false);
  append_candidate(ctx.start_point, true);
  append_candidate(ctx.midpoint, true);
  const HalfedgeHandle heh = rm->halfedge_handle(eh, 0);
  append_candidate(rm->point(rm->to_vertex_handle(heh)), true);
  append_candidate(rm->point(rm->from_vertex_handle(heh)), true);

  bool found_valid_candidate = false;
  for (size_t i = 0; i < candidates.size(); i++)
  {
    if (!phase2_placement_satisfies_hard_constraints(ctx, edge_collapser, candidates[i]))
      continue;

    const double candidate_energy = evaluate_phase2_proxy_energy(ctx, candidates[i]);
    if (!std::isfinite(candidate_energy))
      continue;

    if (!found_valid_candidate || candidate_energy < decision.priority_energy)
    {
      decision.point = candidates[i];
      decision.priority_energy = candidate_energy;
      decision.used_fallback = is_fallback[i];
      found_valid_candidate = true;
    }
  }

  if (!found_valid_candidate)
    return false;

  decision.min_quality = calc_phase2_fan_min_quality(ctx, decision.point);
  decision.nonlinear_residual = evaluate_phase2_refinement_residual(ctx, decision.point);

  if (!should_refine_phase2_placement_with_newton(
    ctx, decision.point, remaining_vertices, target_vertices_num,
    decision.min_quality, decision.nonlinear_residual))
    return true;

  decision.attempted_newton = true;
  Vec3d refined_point;
  double refined_energy = DBL_MAX;
  const auto newton_start_time = std::chrono::steady_clock::now();
  const bool newton_ok = refine_collapse_placement_with_newton(
    eh, edge_collapser, refined_point, refined_energy, &decision.point);
  newton_seconds += std::chrono::duration<double>(
    std::chrono::steady_clock::now() - newton_start_time).count();

  if (newton_ok && phase2_placement_satisfies_hard_constraints(ctx, edge_collapser, refined_point))
  {
    const double current_energy = evaluate_newton_energy(ctx, decision.point);
    const double refined_quality = calc_phase2_fan_min_quality(ctx, refined_point);
    const double refined_residual = evaluate_phase2_refinement_residual(ctx, refined_point);
    const bool energy_improved = !std::isfinite(current_energy) || refined_energy <= current_energy + 1e-12;
    const bool quality_not_worse = refined_quality + 1e-10 >= decision.min_quality;
    const bool residual_not_worse = refined_residual <= decision.nonlinear_residual + 1e-10;
    if (energy_improved || (quality_not_worse && residual_not_worse))
    {
      decision.point = refined_point;
      decision.priority_energy = std::isfinite(refined_energy) ? refined_energy : decision.priority_energy;
      decision.min_quality = refined_quality;
      decision.nonlinear_residual = refined_residual;
      decision.accepted_newton = true;
      decision.used_fallback = false;
      return true;
    }
  }

  decision.newton_failed = true;
  return true;
}

// Optimize only the placement x for the current edge collapse. The edge has
// not been collapsed yet; topology changes happen later in do_phase2_optimization_simplification().
bool CollapseStage::refine_collapse_placement_with_newton(
  EdgeHandle eh, EdgeCollapser& edge_collapser, Vec3d& new_point, double& energy,
  const Vec3d* initial_point)
{
  const Phase2PlacementContext ctx = make_phase2_placement_context(eh, edge_collapser);
  Vec3d x = initial_point && finite_vec(*initial_point) ? *initial_point : ctx.start_point;
  if (!finite_vec(x))
    x = ctx.midpoint;

  if (!is_exact_reject_robustness_mode() && !collapse_target_valid(edge_collapser, x))
  {
    if (collapse_target_valid(edge_collapser, ctx.start_point))
      x = ctx.start_point;
    else if (collapse_target_valid(edge_collapser, ctx.midpoint))
      x = ctx.midpoint;
  }

  double trust_radius = std::max(ctx.local_scale * param->trustRegionRadiusScale, original_diagonal_length * 1e-8);
  energy = evaluate_newton_energy(ctx, x);
  if (!std::isfinite(energy))
    return false;

  for (size_t iter = 0; iter < param->newtonMaxIter; iter++)
  {
    const NewtonDerivatives deriv = autodiff_newton_derivatives(ctx, x);
    if (!std::isfinite(deriv.energy))
      break;

    const Eigen::Vector3d g = to_eigen(deriv.gradient);
    if (g.norm() <= param->newtonGradTol)
    {
      energy = deriv.energy;
      break;
    }

    Eigen::Matrix3d hessian;
    for (int r = 0; r < 3; r++)
      for (int c = 0; c < 3; c++)
        hessian(r, c) = deriv.hessian[r][c];
    hessian = 0.5 * (hessian + hessian.transpose());

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(hessian);
    Eigen::Matrix3d projected_hessian = hessian;
    if (eig.info() == Eigen::Success)
    {
      Eigen::Vector3d values = eig.eigenvalues();
      const double floor_value = std::max(1e-10, 1e-8 * std::max(1.0, values.cwiseAbs().maxCoeff()));
      for (int i = 0; i < 3; i++)
        values[i] = std::max(values[i], floor_value);
      projected_hessian = eig.eigenvectors() * values.asDiagonal() * eig.eigenvectors().transpose();
    }
    else
      projected_hessian += Eigen::Matrix3d::Identity() * 1e-8;

    Eigen::Vector3d step = -projected_hessian.ldlt().solve(g);
    if (!std::isfinite(step.x()) || !std::isfinite(step.y()) || !std::isfinite(step.z()) ||
      step.dot(g) >= 0.0)
      step = -g.normalized() * std::max(ctx.local_scale * 0.1, original_diagonal_length * 1e-8);

    const double max_step = std::max(ctx.local_scale, original_diagonal_length * 1e-8);
    if (step.norm() > max_step)
      step *= max_step / step.norm();

    if (step.norm() <= param->newtonStepTol * std::max(1.0, ctx.local_scale))
      break;

    bool accepted = false;
    if (is_trust_region_solver_mode())
    {
      if (step.norm() > trust_radius)
        step *= trust_radius / step.norm();

      const Vec3d trial = x + from_eigen(step);
      const double trial_energy = evaluate_newton_energy(ctx, trial);
      const double predicted_decrease =
        -(g.dot(step) + 0.5 * step.dot(projected_hessian * step));
      const double actual_decrease = deriv.energy - trial_energy;
      const double rho = predicted_decrease > 0.0 ? actual_decrease / predicted_decrease : -1.0;
      bool robust_ok = true;
      if (std::isfinite(trial_energy) && actual_decrease > 0.0)
      {
        if (is_exact_backtracking_robustness_mode())
          robust_ok = collapse_target_valid(edge_collapser, trial);
        else if (is_ipc_line_search_robustness_mode())
          robust_ok = sampled_path_valid(edge_collapser, x, trial);
      }

      if (robust_ok && std::isfinite(trial_energy) && actual_decrease > 0.0)
      {
        x = trial;
        energy = trial_energy;
        accepted = true;
        if (rho > 0.75)
          trust_radius *= 2.0;
        else if (rho < 0.25)
          trust_radius *= 0.5;
      }
      else
        trust_radius *= 0.5;

      if (trust_radius <= param->newtonStepTol * std::max(1.0, ctx.local_scale))
        break;
    }
    else
    {
      double alpha = 1.0;
      const double descent = std::max(1e-16, -g.dot(step));
      for (size_t ls = 0; ls < param->lineSearchMaxIter; ls++)
      {
        const Vec3d trial = x + alpha * from_eigen(step);
        const double trial_energy = evaluate_newton_energy(ctx, trial);
        if (std::isfinite(trial_energy) &&
          trial_energy <= deriv.energy - 1e-4 * alpha * descent)
        {
          bool robust_ok = true;
          if (is_exact_backtracking_robustness_mode())
            robust_ok = collapse_target_valid(edge_collapser, trial);
          else if (is_ipc_line_search_robustness_mode())
            robust_ok = sampled_path_valid(edge_collapser, x, trial);
          if (!robust_ok)
          {
            alpha *= 0.5;
            continue;
          }

          x = trial;
          energy = trial_energy;
          accepted = true;
          break;
        }
        alpha *= 0.5;
      }
    }

    if (!accepted)
      break;
  }

  if (!collapse_target_valid(edge_collapser, x))
    return false;

  new_point = x;
  energy = evaluate_newton_energy(ctx, x);
  return std::isfinite(energy);
}

// Compute the fast placement used to rank this edge in the Phase 2 queue.
// Most modes use the QEM 3x3 linear solve. newton_only deliberately skips that
// solve and queues a valid starting point for the later nonlinear placement.
bool CollapseStage::compute_phase2_queue_placement_candidate(EdgeHandle eh, Vec3d& new_point, double& energy)
{
  auto edge_collapser = new_edge_collapser();
  if (!edge_collapser.init(eh))
    return false;

  HalfedgeHandle heh = rm->halfedge_handle(eh, 0);
  size_t valence_after_collapsing =
    (rm->valence(rm->to_vertex_handle(heh)) + rm->valence(rm->from_vertex_handle(heh))) - 3;
  if (valence_after_collapsing > param->maxValence)
    return false;

  const Phase2PlacementContext ctx = make_phase2_placement_context(eh, edge_collapser);
  if (ctx.qem_planes.empty())
    return false;

  if (is_phase2_newton_only_strategy())
  {
    std::vector<Vec3d> candidates;
    candidates.push_back(ctx.start_point);
    candidates.push_back(ctx.midpoint);
    candidates.push_back(rm->point(rm->to_vertex_handle(heh)));
    candidates.push_back(rm->point(rm->from_vertex_handle(heh)));

    bool found_valid_candidate = false;
    energy = DBL_MAX;
    for (const Vec3d& candidate : candidates)
    {
      if (!phase2_placement_satisfies_hard_constraints(ctx, edge_collapser, candidate))
        continue;

      const double candidate_energy = evaluate_phase2_proxy_energy(ctx, candidate);
      if (!std::isfinite(candidate_energy))
        continue;

      if (!found_valid_candidate || candidate_energy < energy)
      {
        new_point = candidate;
        energy = candidate_energy;
        found_valid_candidate = true;
      }
    }
    return found_valid_candidate;
  }

  Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
  Eigen::Vector3d b = Eigen::Vector3d::Zero();

  double curvature_multiplier = 1.0;
  if (param->curvatureMode == "weighted_qem" || param->curvatureMode == "weighted-qem")
  {
    const double scaled_curvature = ctx.edge_curvature * ctx.local_scale;
    curvature_multiplier += param->curvatureWeight * sqr(scaled_curvature);
  }

  const double inv_scale_sqr = 1.0 / std::max(sqr(ctx.local_scale), 1e-24);
  for (const QEMPlane& plane : ctx.qem_planes)
  {
    const Eigen::Vector3d n = to_eigen(plane.normal);
    const double w = param->qemWeight * curvature_multiplier * plane.weight * inv_scale_sqr;
    A += w * (n * n.transpose());
    b += w * plane.offset * n;
  }

  const double regularization =
    std::max(1e-12, 1e-8 * std::max(1.0, std::abs(A.trace())));
  const Eigen::Matrix3d A_reg = A + regularization * Eigen::Matrix3d::Identity();
  const Eigen::Vector3d rhs = regularization * to_eigen(ctx.midpoint) - b;
  Eigen::Vector3d x = A_reg.ldlt().solve(rhs);

  if (!std::isfinite(x.x()) || !std::isfinite(x.y()) || !std::isfinite(x.z()))
    x = to_eigen(ctx.start_point);
  if (!std::isfinite(x.x()) || !std::isfinite(x.y()) || !std::isfinite(x.z()))
    x = to_eigen(ctx.midpoint);

  const double max_offset = std::max(2.0 * ctx.local_scale, original_diagonal_length * 1e-8);
  Eigen::Vector3d offset = x - to_eigen(ctx.midpoint);
  if (offset.norm() > max_offset)
    x = to_eigen(ctx.midpoint) + offset.normalized() * max_offset;

  new_point = from_eigen(x);
  energy = evaluate_phase2_proxy_energy(ctx, new_point);
  return std::isfinite(energy);
}

bool CollapseStage::enqueue_phase2_candidate(EdgeHandle eh, size_t state)
{
  if (!eh.is_valid() || rm->status(eh).deleted())
    return false;

  Vec3d new_point;
  double energy = DBL_MAX;
  if (!compute_phase2_queue_placement_candidate(eh, new_point, energy))
    return false;

  phase2_edges_to_collapse.emplace(eh, state, -energy, new_point);
  return true;
}

// Initialize the candidate queue. Each queue item stores an edge and its current
// fast placement estimate; it does not perform any collapse.
void CollapseStage::initialize_phase2_candidates()
{
  update_states.clear();
  update_states.resize(rm->n_edges(), 0);
  phase2_edges_to_collapse = Phase2EdgeRewardQueue();

  avg_edge_length = 0.0;
  for (EdgeHandle eh : rm->edges())
    avg_edge_length += rm->data(eh).edge_length;
  avg_edge_length = rm->n_edges() == 0 ? 0.0 : avg_edge_length / rm->n_edges();

  avg_original_edge_length = 0.0;
  for (EdgeHandle eh : om->edges())
    avg_original_edge_length += om->data(eh).edge_length;
  avg_original_edge_length = om->n_edges() == 0 ? 0.0 : avg_original_edge_length / om->n_edges();

  Logger::user_logger->info(
    "phase 2 optimization [{}]: initializing {} edge candidates.",
    param->phase2PlacementStrategy, rm->n_edges());

  size_t evaluated_edges = 0;
  size_t enqueued_edges = 0;
  for (EdgeHandle eh : rm->edges())
  {
    if (enqueue_phase2_candidate(eh, 0))
      enqueued_edges++;

    evaluated_edges++;
    if (evaluated_edges % 10000 == 0)
    {
      Logger::user_logger->info(
        "phase 2 optimization candidates: {}/{} edges, {} enqueued.",
        evaluated_edges, rm->n_edges(), enqueued_edges);
    }
  }

  Logger::user_logger->info(
    "phase 2 optimization candidates: {}/{} edges, {} enqueued.",
    evaluated_edges, rm->n_edges(), enqueued_edges);
}

// After one collapse, only nearby edges have stale placement scores.
// Recompute and requeue those local candidates.
void CollapseStage::update_phase2_after_collapsing(VertexHandle collapsed_center)
{
  std::vector<EdgeHandle> affected_edges = find_1rv_1re(rm, collapsed_center);
  for (EdgeHandle eh : affected_edges)
  {
    if (!eh.is_valid() || rm->status(eh).deleted())
      continue;
    update_states[eh.idx()]++;
    enqueue_phase2_candidate(eh, update_states[eh.idx()]);
  }
}

bool CollapseStage::refine_vertex_relocation_with_newton(
  VertexHandle vh, VertexRelocater& vertex_relocater, Vec3d& new_point, double& energy)
{
  const Phase2PlacementContext ctx = make_phase2_vertex_relocation_context(vh, vertex_relocater);
  if (ctx.fan_edges.empty())
    return false;

  Vec3d x = finite_vec(ctx.start_point) ? ctx.start_point : ctx.midpoint;
  if (!finite_vec(x))
    return false;

  double trust_radius = std::max(ctx.local_scale * param->trustRegionRadiusScale, original_diagonal_length * 1e-8);
  energy = evaluate_newton_energy(ctx, x);
  if (!std::isfinite(energy))
    return false;

  for (size_t iter = 0; iter < param->newtonMaxIter; iter++)
  {
    const NewtonDerivatives deriv = autodiff_newton_derivatives(ctx, x);
    if (!std::isfinite(deriv.energy))
      break;

    const Eigen::Vector3d g = to_eigen(deriv.gradient);
    if (g.norm() <= param->newtonGradTol)
    {
      energy = deriv.energy;
      break;
    }

    Eigen::Matrix3d hessian;
    for (int r = 0; r < 3; r++)
      for (int c = 0; c < 3; c++)
        hessian(r, c) = deriv.hessian[r][c];
    hessian = 0.5 * (hessian + hessian.transpose());

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(hessian);
    Eigen::Matrix3d projected_hessian = hessian;
    if (eig.info() == Eigen::Success)
    {
      Eigen::Vector3d values = eig.eigenvalues();
      const double floor_value = std::max(1e-10, 1e-8 * std::max(1.0, values.cwiseAbs().maxCoeff()));
      for (int i = 0; i < 3; i++)
        values[i] = std::max(values[i], floor_value);
      projected_hessian = eig.eigenvectors() * values.asDiagonal() * eig.eigenvectors().transpose();
    }
    else
      projected_hessian += Eigen::Matrix3d::Identity() * 1e-8;

    Eigen::Vector3d step = -projected_hessian.ldlt().solve(g);
    if (!std::isfinite(step.x()) || !std::isfinite(step.y()) || !std::isfinite(step.z()) ||
      step.dot(g) >= 0.0)
      step = -g.normalized() * std::max(ctx.local_scale * 0.1, original_diagonal_length * 1e-8);

    const double max_step = std::max(ctx.local_scale, original_diagonal_length * 1e-8);
    if (step.norm() > max_step)
      step *= max_step / step.norm();

    if (step.norm() <= param->newtonStepTol * std::max(1.0, ctx.local_scale))
      break;

    bool accepted = false;
    if (is_trust_region_solver_mode())
    {
      if (step.norm() > trust_radius)
        step *= trust_radius / step.norm();

      const Vec3d trial = x + from_eigen(step);
      const double trial_energy = evaluate_newton_energy(ctx, trial);
      const double predicted_decrease =
        -(g.dot(step) + 0.5 * step.dot(projected_hessian * step));
      const double actual_decrease = deriv.energy - trial_energy;
      const double rho = predicted_decrease > 0.0 ? actual_decrease / predicted_decrease : -1.0;

      if (std::isfinite(trial_energy) && actual_decrease > 0.0)
      {
        x = trial;
        energy = trial_energy;
        accepted = true;
        if (rho > 0.75)
          trust_radius *= 2.0;
        else if (rho < 0.25)
          trust_radius *= 0.5;
      }
      else
        trust_radius *= 0.5;

      if (trust_radius <= param->newtonStepTol * std::max(1.0, ctx.local_scale))
        break;
    }
    else
    {
      double alpha = 1.0;
      const double descent = std::max(1e-16, -g.dot(step));
      for (size_t ls = 0; ls < param->lineSearchMaxIter; ls++)
      {
        const Vec3d trial = x + alpha * from_eigen(step);
        const double trial_energy = evaluate_newton_energy(ctx, trial);
        if (std::isfinite(trial_energy) &&
          trial_energy <= deriv.energy - 1e-4 * alpha * descent)
        {
          x = trial;
          energy = trial_energy;
          accepted = true;
          break;
        }
        alpha *= 0.5;
      }
    }

    if (!accepted)
      break;
  }

  new_point = x;
  energy = evaluate_newton_energy(ctx, x);
  return finite_vec(new_point) && std::isfinite(energy);
}

void CollapseStage::do_phase2_final_newton_relocation()
{
  avg_edge_length = 0.0;
  for (EdgeHandle eh : rm->edges())
    avg_edge_length += rm->data(eh).edge_length;
  avg_edge_length = rm->n_edges() == 0 ? 0.0 : avg_edge_length / rm->n_edges();

  avg_original_edge_length = 0.0;
  for (EdgeHandle eh : om->edges())
    avg_original_edge_length += om->data(eh).edge_length;
  avg_original_edge_length = om->n_edges() == 0 ? 0.0 : avg_original_edge_length / om->n_edges();

  std::vector<VertexHandle> vertices;
  vertices.reserve(rm->n_vertices());
  for (VertexHandle vh : rm->vertices())
    vertices.push_back(vh);

  auto vertex_relocater = new_vertex_relocater();
  vertex_relocater.set_flags(
    /*update_links*/false, /*update_target_length*/false,
    /*update_normals*/true, /*check_wrinkle*/false);

  size_t attempted = 0;
  size_t relocated = 0;
  size_t failed = 0;
  double newton_seconds = 0.0;
  const auto start_time = std::chrono::steady_clock::now();

  for (VertexHandle vh : vertices)
  {
    if (!vh.is_valid() || rm->status(vh).deleted())
      continue;
    if (!vertex_relocater.init(vh))
      continue;

    const Phase2PlacementContext ctx = make_phase2_vertex_relocation_context(vh, vertex_relocater);
    const Vec3d current_point = rm->point(vh);
    const double current_energy = evaluate_newton_energy(ctx, current_point);

    Vec3d new_point;
    double new_energy = DBL_MAX;
    attempted++;
    if (!refine_vertex_relocation_with_newton(vh, vertex_relocater, new_point, new_energy) ||
      !std::isfinite(current_energy) ||
      new_energy > current_energy + 1e-12 ||
      (new_point - current_point).length() <= std::max(ctx.local_scale * 1e-8, original_diagonal_length * 1e-12))
    {
      failed++;
      continue;
    }

    if (vertex_relocater.try_relocate_vertex(new_point))
      relocated++;
    else
      failed++;
  }

  newton_seconds = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - start_time).count();
  Logger::user_logger->info(
    "phase 2 final Newton relocation: attempted {}, relocated {}, failed {}, time {:.3f}s.",
    attempted, relocated, failed, newton_seconds);

  rm->garbage_collection();
  lrt->collect_garbage();
  init_one_ring_faces(rm);
}

// Main Phase 2 replacement:
// 1. Build a priority queue from fast placement estimates.
// 2. Pop one edge, choose its actual placement according to phase2PlacementStrategy.
// 3. Commit the collapse with exact validity checks.
// 4. Requeue only the affected one-ring edges.
// 5. In final_newton mode, polish fixed-topology vertex positions once.
void CollapseStage::do_phase2_optimization_simplification(size_t target_vertices_num)
{
  if (target_vertices_num == 0 || rm->n_vertices() <= target_vertices_num)
    return;

  Logger::user_logger->info(
    "phase 2 optimization simplification [{}]: {} -> {} vertices.",
    param->phase2PlacementStrategy, rm->n_vertices(), target_vertices_num);

  size_t total_collapsed_edges = 0;
  size_t total_attempted_edges = 0;
  size_t remaining_vertices = rm->n_vertices();
  size_t no_progress_passes = 0;
  const size_t max_phase2_passes = 12;

  for (size_t pass = 1; remaining_vertices > target_vertices_num && pass <= max_phase2_passes; pass++)
  {
    Logger::user_logger->info("phase 2 optimization pass {}.", pass);
    initialize_phase2_candidates();

    auto edge_collapser = new_edge_collapser();
    edge_collapser.set_flags(
      /*update_links*/false, /*update_target_length*/false, /*update_normals*/true,
      /*check_wrinkle*/false, /*check_selfinter*/true, /*check_inter*/true);

    size_t pass_collapsed_edges = 0;
    size_t pass_attempted_edges = 0;
    size_t pass_newton_attempted_edges = 0;
    size_t pass_newton_refined_edges = 0;
    size_t pass_newton_failed_edges = 0;
    size_t pass_fallback_edges = 0;
    size_t pass_placement_failed_edges = 0;
    size_t pass_collapse_failed_edges = 0;
    double pass_newton_seconds = 0.0;
    double pass_exact_collapse_seconds = 0.0;
    auto last_progress_log_time = std::chrono::steady_clock::now();
    while (remaining_vertices > target_vertices_num && !phase2_edges_to_collapse.empty())
    {
      Phase2EdgeReward edge_reward = phase2_edges_to_collapse.top();
      phase2_edges_to_collapse.pop();

      if (edge_reward.state < update_states[edge_reward.eh.idx()])
        continue;
      if (rm->status(edge_reward.eh).deleted())
        continue;

      if (!edge_collapser.init(edge_reward.eh))
        continue;

      pass_attempted_edges++;
      total_attempted_edges++;
      Phase2PlacementDecision placement;
      double edge_newton_seconds = 0.0;
      if (!choose_phase2_collapse_placement(
        edge_reward.eh, edge_collapser, edge_reward.initial_point,
        remaining_vertices, target_vertices_num, placement, edge_newton_seconds))
      {
        pass_placement_failed_edges++;
        const auto now = std::chrono::steady_clock::now();
        if (pass_attempted_edges % 100 == 0 ||
          std::chrono::duration<double>(now - last_progress_log_time).count() >= 10.0)
        {
          Logger::user_logger->info(
            "phase 2 optimization pass {} progress: attempted {}, collapsed {}, newton_attempted {}, newton_refined {}, newton_failed {}, fallback {}, placement_failed {}, collapse_failed {}, queue {}, remain {}, newton_time {:.3f}s, exact_collapse_time {:.3f}s.",
            pass, pass_attempted_edges, pass_collapsed_edges,
            pass_newton_attempted_edges, pass_newton_refined_edges,
            pass_newton_failed_edges, pass_fallback_edges,
            pass_placement_failed_edges, pass_collapse_failed_edges,
            phase2_edges_to_collapse.size(), remaining_vertices,
            pass_newton_seconds, pass_exact_collapse_seconds);
          last_progress_log_time = now;
        }
        continue;
      }

      pass_newton_seconds += edge_newton_seconds;
      if (placement.attempted_newton)
        pass_newton_attempted_edges++;
      if (placement.accepted_newton)
        pass_newton_refined_edges++;
      if (placement.newton_failed)
        pass_newton_failed_edges++;
      if (placement.used_fallback)
        pass_fallback_edges++;

      const auto collapse_start_time = std::chrono::steady_clock::now();
      if (edge_collapser.try_collapse_edge(placement.point, nullptr))
      {
        pass_exact_collapse_seconds += std::chrono::duration<double>(
          std::chrono::steady_clock::now() - collapse_start_time).count();
        VertexHandle center_v = edge_collapser.get_collapsed_center();
        update_phase2_after_collapsing(center_v);
        pass_collapsed_edges++;
        total_collapsed_edges++;
        remaining_vertices--;

        if (total_collapsed_edges % 1000 == 0)
        {
          Logger::user_logger->info(
            "phase 2 optimization collapsed {} edges; {} vertices remain.",
            total_collapsed_edges, remaining_vertices);
        }
      }
      else
      {
        pass_exact_collapse_seconds += std::chrono::duration<double>(
          std::chrono::steady_clock::now() - collapse_start_time).count();
        pass_collapse_failed_edges++;
      }

      const auto now = std::chrono::steady_clock::now();
      if (pass_attempted_edges % 100 == 0 ||
        std::chrono::duration<double>(now - last_progress_log_time).count() >= 10.0)
      {
        Logger::user_logger->info(
          "phase 2 optimization pass {} progress: attempted {}, collapsed {}, newton_attempted {}, newton_refined {}, newton_failed {}, fallback {}, placement_failed {}, collapse_failed {}, queue {}, remain {}, newton_time {:.3f}s, exact_collapse_time {:.3f}s.",
          pass, pass_attempted_edges, pass_collapsed_edges,
          pass_newton_attempted_edges, pass_newton_refined_edges,
          pass_newton_failed_edges, pass_fallback_edges,
          pass_placement_failed_edges, pass_collapse_failed_edges,
          phase2_edges_to_collapse.size(), remaining_vertices,
          pass_newton_seconds, pass_exact_collapse_seconds);
        last_progress_log_time = now;
      }
    }

    Logger::user_logger->info(
      "phase 2 optimization pass {} collapsed {} edges after attempting {} candidates (newton_attempted {}, newton_refined {}, newton_failed {}, fallback {}, placement_failed {}, collapse_failed {}, newton_time {:.3f}s, exact_collapse_time {:.3f}s).",
      pass, pass_collapsed_edges, pass_attempted_edges,
      pass_newton_attempted_edges, pass_newton_refined_edges,
      pass_newton_failed_edges, pass_fallback_edges,
      pass_placement_failed_edges, pass_collapse_failed_edges,
      pass_newton_seconds, pass_exact_collapse_seconds);

    rm->garbage_collection();
    lrt->collect_garbage();
    init_one_ring_faces(rm);
    remaining_vertices = rm->n_vertices();
    Logger::user_logger->info("[{}] vertices and [{}] faces remained.", rm->n_vertices(), rm->n_faces());

    if (pass_collapsed_edges == 0)
    {
      no_progress_passes++;
      if (no_progress_passes >= 2)
        break;
    }
    else
      no_progress_passes = 0;
  }

  Logger::user_logger->info(
    "phase 2 optimization collapsed {} edges after attempting {} candidates.",
    total_collapsed_edges, total_attempted_edges);

  if (is_phase2_final_newton_strategy())
    do_phase2_final_newton_relocation();

  Logger::user_logger->info("[{}] vertices and [{}] faces remained.", rm->n_vertices(), rm->n_faces());
}

bool CollapseStage::find_collapse_hausdorff_deviation(
  EdgeHandle eh, double& local_hd_before, double& local_hd_after, Vec3d& new_point, double& priority_score)
{
  auto edge_collapser = new_edge_collapser();

  if (!edge_collapser.init(eh))
    return false;

  // constrain valence
  HalfedgeHandle heh = rm->halfedge_handle(eh, 0);
  size_t valence_after_collapsing =
    (rm->valence(rm->to_vertex_handle(heh)) + rm->valence(rm->from_vertex_handle(heh))) - 3;
  if (valence_after_collapsing > param->maxValence)
    return false;

  if (is_optimization_collapse_placement_method())
    edge_collapser.set_flags(
      /*update_links*/false, /*update_target_length*/false, /*update_normals*/false,
      /*check_wrinkle*/true, /*check_selfinter*/true, /*check_inter*/true);

  local_hd_before = edge_collapser.local_Hausdorff_before_collapsing();

  const bool use_post_metric = is_post_metric_priority_mode();
  const double pre_collapse_metric = use_post_metric ? calc_pre_collapse_metric(eh) : 0.0;

  if (is_optimization_collapse_placement_method())
  {
    double optimized_energy = DBL_MAX;
    if (!refine_collapse_placement_with_newton(eh, edge_collapser, new_point, optimized_energy))
      return false;

    VertexHandle local_center_v;
    auto local_mesh = construct_local_mesh(rm, edge_collapser.get_halfedges(), new_point, local_center_v);
    double hd_threshold = use_post_metric || allow_negtive ? max_distance_error : cage_infinite_fp;
    local_hd_after = edge_collapser.local_Hausdorff_after_collapsing(local_mesh.get(), new_point, hd_threshold);
    if (local_hd_after == DBL_MAX)
      return false;

    priority_score = use_post_metric ? calc_post_collapse_metric_score(pre_collapse_metric, local_mesh.get()) : -optimized_energy;
    return true;
  }

  // get candidate points
  std::vector<Vec3d> candidate_points = generate_candidate_points_for_collapse(eh, edge_collapser);

  // construct local mesh
  VertexHandle local_center_v;
  std::vector<SMeshT> local_meshes = construct_local_meshes(rm, edge_collapser.get_halfedges(), candidate_points, local_center_v);

  // find optimal point that minimize local hausdorff distance.
  size_t minimal_idx = 0;
  double minimal_local_hd = DBL_MAX;
  double maximal_priority_score = -DBL_MAX;

#if USE_TREE_SEARCH
  ot->set_hint(candidate_points[0]);
#endif
#pragma omp parallel for schedule(dynamic)
  for (int i = 0;i < (int)candidate_points_size;i++)
  {
    double hd_threshold = use_post_metric ? max_distance_error : minimal_local_hd;
    double local_hd_i = edge_collapser.local_Hausdorff_after_collapsing(&local_meshes[i], candidate_points[i], hd_threshold);
    double candidate_priority_score = 0.0;
    bool candidate_ok = local_hd_i != DBL_MAX;
    if (candidate_ok && use_post_metric && local_hd_i >= max_distance_error)
      candidate_ok = false;
    if (candidate_ok && use_post_metric)
    {
      candidate_priority_score = calc_post_collapse_metric_score(pre_collapse_metric, &local_meshes[i]);
      if (is_length_quality_priority_mode())
      {
        const double post_quality = pre_collapse_metric + candidate_priority_score;
        candidate_ok = is_length_quality_allowed(pre_collapse_metric, post_quality);
      }
      else if (is_hard_post_metric_priority_mode())
      {
        candidate_ok = candidate_priority_score >= 0.0;
      }
    }
  #pragma omp critical
    if (candidate_ok)
    {
      if (use_post_metric)
      {
        if (candidate_priority_score > maximal_priority_score ||
          (std::abs(candidate_priority_score - maximal_priority_score) <= 1e-12 && local_hd_i < minimal_local_hd))
        {
          maximal_priority_score = candidate_priority_score;
          minimal_local_hd = local_hd_i;
          minimal_idx = i;
        }
      }
      else if (local_hd_i < minimal_local_hd)
      {
        minimal_local_hd = local_hd_i;
        minimal_idx = i;
      }
    }
  }
#if USE_TREE_SEARCH
  ot->unset_hint();
#endif
  // return result
  if (minimal_local_hd != DBL_MAX)
  {
    // found a available point
    local_hd_after = minimal_local_hd;
    new_point = candidate_points[minimal_idx];
    priority_score = use_post_metric ? maximal_priority_score : 0.0;
    return true;
  }
  else return false;
}

bool CollapseStage::is_length_priority_mode() const
{
  return param->priorityMode == "length" || param->priorityMode == "edge_length";
}

bool CollapseStage::is_length_quality_priority_mode() const
{
  return param->priorityMode == "length_quality" || param->priorityMode == "length-quality";
}

bool CollapseStage::is_length_quality_weighted_submode() const
{
  return param->lengthQualitySubMode == "weighted" || param->lengthQualitySubMode == "weight";
}

bool CollapseStage::is_length_quality_relative_reject_submode() const
{
  return param->lengthQualitySubMode == "relative_reject" || param->lengthQualitySubMode == "relative-reject" ||
    param->lengthQualitySubMode == "relative";
}

bool CollapseStage::is_length_quality_absolute_reject_submode() const
{
  return param->lengthQualitySubMode == "absolute_reject" || param->lengthQualitySubMode == "absolute-reject" ||
    param->lengthQualitySubMode == "absolute";
}

bool CollapseStage::is_length_quality_lexicographic_submode() const
{
  return param->lengthQualitySubMode == "lexicographic" || param->lengthQualitySubMode == "lex";
}

bool CollapseStage::is_post_edge_length_priority_mode() const
{
  return param->priorityMode == "post_edge_length" || param->priorityMode == "post-edge-length" ||
    param->priorityMode == "post_edge_length_hard" || param->priorityMode == "post-edge-length-hard";
}

bool CollapseStage::is_post_face_area_priority_mode() const
{
  return param->priorityMode == "post_face_area" || param->priorityMode == "post-face-area" ||
    param->priorityMode == "post_face_area_hard" || param->priorityMode == "post-face-area-hard";
}

bool CollapseStage::is_triangle_quality_priority_mode() const
{
  return param->priorityMode == "triangle_quality" || param->priorityMode == "triangle-quality" ||
    param->priorityMode == "triangle_quality_hard" || param->priorityMode == "triangle-quality-hard";
}

bool CollapseStage::is_post_metric_priority_mode() const
{
  return is_post_edge_length_priority_mode() ||
    is_post_face_area_priority_mode() ||
    is_triangle_quality_priority_mode() ||
    is_length_quality_priority_mode();
}

bool CollapseStage::is_hard_post_metric_priority_mode() const
{
  return param->priorityMode == "post_edge_length_hard" || param->priorityMode == "post-edge-length-hard" ||
    param->priorityMode == "post_face_area_hard" || param->priorityMode == "post-face-area-hard" ||
    param->priorityMode == "triangle_quality_hard" || param->priorityMode == "triangle-quality-hard";
}

bool CollapseStage::is_length_quality_allowed(double pre_quality, double post_quality) const
{
  if (is_length_quality_relative_reject_submode())
    return post_quality >= pre_quality * param->lengthQualityDegradationRatio;
  if (is_length_quality_absolute_reject_submode())
    return post_quality >= param->lengthQualityMinQuality;
  return true;
}

double CollapseStage::calc_length_quality_reward(double normalized_length_score, double quality_delta) const
{
  if (is_length_quality_weighted_submode() ||
    (!is_length_quality_relative_reject_submode() &&
      !is_length_quality_absolute_reject_submode() &&
      !is_length_quality_lexicographic_submode()))
    return normalized_length_score + param->lengthQualityWeight * quality_delta;
  return normalized_length_score + quality_delta;
}

double CollapseStage::calc_pre_collapse_metric(EdgeHandle eh) const
{
  std::set<FaceHandle> faces;
  one_ring_faces_around_edge(rm, rm->halfedge_handle(eh, 0), faces);

  if (is_post_edge_length_priority_mode())
    return calc_max_edge_length(rm, faces);
  if (is_post_face_area_priority_mode())
    return calc_max_face_area(rm, faces);
  if (is_triangle_quality_priority_mode() || is_length_quality_priority_mode())
    return calc_min_triangle_quality(rm, faces);
  return 0.0;
}

double CollapseStage::calc_post_collapse_metric_score(double pre_metric, SMeshT* local_mesh) const
{
  if (is_post_edge_length_priority_mode())
    return pre_metric - calc_max_edge_length(local_mesh);
  if (is_post_face_area_priority_mode())
    return pre_metric - calc_max_face_area(local_mesh);
  if (is_triangle_quality_priority_mode() || is_length_quality_priority_mode())
    return calc_min_triangle_quality(local_mesh) - pre_metric;
  return 0.0;
}

double CollapseStage::calc_max_edge_length(SMeshT* mesh) const
{
  double max_edge_length = 0.0;
  for (EdgeHandle eh : mesh->edges())
    max_edge_length = std::max(max_edge_length, mesh->data(eh).edge_length);
  return max_edge_length;
}

double CollapseStage::calc_max_edge_length(SMeshT* mesh, const std::set<FaceHandle>& faces) const
{
  double max_edge_length = 0.0;
  for (FaceHandle fh : faces)
    for (EdgeHandle eh : mesh->fe_range(fh))
      max_edge_length = std::max(max_edge_length, mesh->data(eh).edge_length);
  return max_edge_length;
}

double CollapseStage::calc_max_face_area(SMeshT* mesh) const
{
  double max_face_area = 0.0;
  for (FaceHandle fh : mesh->faces())
    max_face_area = std::max(max_face_area, mesh->data(fh).face_area);
  return max_face_area;
}

double CollapseStage::calc_max_face_area(SMeshT* mesh, const std::set<FaceHandle>& faces) const
{
  double max_face_area = 0.0;
  for (FaceHandle fh : faces)
    max_face_area = std::max(max_face_area, mesh->data(fh).face_area);
  return max_face_area;
}

double CollapseStage::calc_triangle_quality(SMeshT* mesh, FaceHandle fh) const
{
  Vec3d pts[3];
  size_t vertex_count = 0;
  for (VertexHandle vh : mesh->fv_range(fh))
  {
    if (vertex_count >= 3)
      return 0.0;
    pts[vertex_count++] = mesh->point(vh);
  }
  if (vertex_count != 3)
    return 0.0;

  const double a = (pts[1] - pts[0]).length();
  const double b = (pts[2] - pts[1]).length();
  const double c = (pts[0] - pts[2]).length();
  const double denom = a * a + b * b + c * c;
  if (denom <= 0.0)
    return 0.0;

  const double area = 0.5 * (pts[1] - pts[0]).cross(pts[2] - pts[0]).length();
  return 4.0 * std::sqrt(3.0) * area / denom;
}

double CollapseStage::calc_min_triangle_quality(SMeshT* mesh) const
{
  double min_quality = DBL_MAX;
  for (FaceHandle fh : mesh->faces())
    min_quality = std::min(min_quality, calc_triangle_quality(mesh, fh));
  return min_quality == DBL_MAX ? 0.0 : min_quality;
}

double CollapseStage::calc_min_triangle_quality(SMeshT* mesh, const std::set<FaceHandle>& faces) const
{
  double min_quality = DBL_MAX;
  for (FaceHandle fh : faces)
    min_quality = std::min(min_quality, calc_triangle_quality(mesh, fh));
  return min_quality == DBL_MAX ? 0.0 : min_quality;
}

bool CollapseStage::try_enqueue_collapse_candidate(
  EdgeHandle eh, size_t state, double local_hd_before, double local_hd_after, const Vec3d& new_point, double priority_score)
{
  if (is_length_quality_priority_mode())
  {
    if (local_hd_after >= max_distance_error)
      return false;

    const double length_score = avg_edge_length - rm->data(eh).edge_length;
    if (length_score <= 0.0)
      return false;

    const double pre_quality = calc_pre_collapse_metric(eh);
    const double post_quality = pre_quality + priority_score;
    if (!is_length_quality_allowed(pre_quality, post_quality))
      return false;

    const double normalized_length_score = avg_edge_length > 0.0 ? length_score / avg_edge_length : length_score;
    if (is_length_quality_lexicographic_submode())
      edges_to_collapse.emplace(eh, state, normalized_length_score, new_point, priority_score);
    else
      edges_to_collapse.emplace(eh, state, calc_length_quality_reward(normalized_length_score, priority_score), new_point);
    return true;
  }

  if (is_length_priority_mode())
  {
    if (local_hd_after >= max_distance_error)
      return false;

    const double length_score = avg_edge_length - rm->data(eh).edge_length;
    if (length_score <= 0.0)
      return false;

    edges_to_collapse.emplace(eh, state, length_score, new_point);
    return true;
  }

  if (is_post_metric_priority_mode())
  {
    if (local_hd_after >= max_distance_error)
      return false;
    if (is_hard_post_metric_priority_mode() && priority_score < 0.0)
      return false;

    edges_to_collapse.emplace(eh, state, priority_score, new_point);
    return true;
  }

  if (allow_negtive)
  {
    // Original behavior: collapse smallest post-collapse local HD first.
    if (local_hd_after < max_distance_error)
    {
      double error = local_hd_after + 0.1 * (rm->data(eh).edge_length - avg_edge_length);
      edges_to_collapse.emplace(eh, state, -error, new_point);
      return true;
    }
  }
  else
  {
    // Original behavior: collapse largest local HD gain first.
    if (local_hd_before - local_hd_after >= 0.0)
    {
      double error = local_hd_before - local_hd_after + 0.1 * (avg_edge_length - rm->data(eh).edge_length);
      edges_to_collapse.emplace(eh, state, error, new_point);
      return true;
    }
  }

  return false;
}

void CollapseStage::initialize_collapse_edges_reward()
{
  // clear
  update_states.resize(rm->n_edges());
  std::fill(update_states.begin(), update_states.end(), 0);
  edges_to_collapse = CollapseEdgeRewardQueue();

  // calculate average edge length
  avg_edge_length = 0.0;
  for (EdgeHandle eh : rm->edges())
    avg_edge_length += rm->data(eh).edge_length;
  avg_edge_length /= rm->n_edges();

  Logger::user_logger->info("initializing collapse candidates over {} edges.", rm->n_edges());

  // initialize for all edges
  size_t evaluated_edges = 0;
  size_t enqueued_edges = 0;
  for (EdgeHandle eh : rm->edges())
  {
    double local_hd_after, local_hd_before;
    double priority_score;
    Vec3d new_point;
    if (find_collapse_hausdorff_deviation(eh, local_hd_before, local_hd_after, new_point, priority_score))
    {
      if (try_enqueue_collapse_candidate(eh, 0, local_hd_before, local_hd_after, new_point, priority_score))
        enqueued_edges++;
    }

    evaluated_edges++;
    if (evaluated_edges % 1000 == 0)
    {
      Logger::user_logger->info(
        "initialized collapse candidates: {}/{} edges, {} enqueued.",
        evaluated_edges, rm->n_edges(), enqueued_edges);
    }
  }
  Logger::user_logger->info(
    "initialized collapse candidates: {}/{} edges, {} enqueued.",
    evaluated_edges, rm->n_edges(), enqueued_edges);
}

void CollapseStage::update_after_collapsing(VertexHandle collapsed_center)
{
  // update in out error
  std::vector<FaceHandle> faces;
  std::vector<EdgeHandle> edges;
  for (FaceHandle vf : rm->vf_range(collapsed_center)) faces.push_back(vf);
  for (EdgeHandle ve : rm->ve_range(collapsed_center)) edges.push_back(ve);

#ifdef USE_TREE_SEARCH
  calc_face_out_error(rm, om, ot, faces, cage_infinite_fp);
  calc_edge_out_error(rm, om, ot, edges, cage_infinite_fp);
  calc_vertex_out_error(rm, om, ot, collapsed_center);
#else
  calc_face_out_error(rm, om, og, faces, cage_infinite_fp);
  calc_edge_out_error(rm, om, og, edges, cage_infinite_fp);
  calc_vertex_out_error(rm, om, og, collapsed_center);
#endif
  calc_face_in_error(rm, faces);

  calc_out_surround_error(rm, collapsed_center);
  for (VertexHandle vv : rm->vv_range(collapsed_center))
    calc_out_surround_error(rm, vv);

  // update optimal vertex position
  std::vector<EdgeHandle> affected_edges = find_1rv_1re(rm, collapsed_center);
  for (EdgeHandle eh : affected_edges)
  {
    update_states[eh.idx()]++;

    double local_hd_after, local_hd_before;
    double priority_score;
    Vec3d new_point;
    if (find_collapse_hausdorff_deviation(eh, local_hd_before, local_hd_after, new_point, priority_score))
      try_enqueue_collapse_candidate(eh, update_states[eh.idx()], local_hd_before, local_hd_after, new_point, priority_score);
  }
}

CollapseStage::EdgeSide CollapseStage::classify_edge_side(EdgeHandle eh) const
{
  if (!eh.is_valid() || rm->status(eh).deleted())
    return EdgeSide::Unknown;

  HalfedgeHandle heh = rm->halfedge_handle(eh, 0);
  const Vec3d midpoint =
    (rm->point(rm->from_vertex_handle(heh)) + rm->point(rm->to_vertex_handle(heh))) * 0.5;

  Vec3d closest_point;
  FaceHandle closest_face;
#ifdef USE_TREE_SEARCH
  if (!ot)
    return EdgeSide::Unknown;
  auto closest = ot->closest_point(midpoint);
  closest_point = closest.first;
  closest_face = FaceHandle(closest.second);
#else
  if (!og)
    return EdgeSide::Unknown;
  auto closest = og->closest_point(midpoint);
  closest_point = closest.first.first;
  closest_face = FaceHandle(closest.second);
#endif

  if (!closest_face.is_valid() || closest_face.idx() < 0 ||
    closest_face.idx() >= static_cast<int>(om->n_faces()))
    return EdgeSide::Unknown;

  const Vec3d& normal = om->normal(closest_face);
  const double signed_distance = (midpoint - closest_point) | normal;
  const double side_eps = original_diagonal_length * 1e-9;

  if (signed_distance > side_eps)
    return EdgeSide::Outside;
  if (signed_distance < -side_eps)
    return EdgeSide::Inside;
  return EdgeSide::OnSurface;
}

void CollapseStage::add_edge_side(EdgeSideStats& stats, EdgeSide side) const
{
  switch (side)
  {
  case EdgeSide::Inside:
    stats.inside++;
    break;
  case EdgeSide::Outside:
    stats.outside++;
    break;
  case EdgeSide::OnSurface:
    stats.on_surface++;
    break;
  default:
    stats.unknown++;
    break;
  }
}

CollapseStage::EdgeSideStats CollapseStage::collect_candidate_edge_side_stats() const
{
  EdgeSideStats stats;
  CollapseEdgeRewardQueue queue_copy = edges_to_collapse;
  while (!queue_copy.empty())
  {
    const auto edge_reward = queue_copy.top();
    queue_copy.pop();
    add_edge_side(stats, classify_edge_side(edge_reward.eh));
  }
  return stats;
}

void CollapseStage::log_edge_side_stats(const char* label, const EdgeSideStats& stats) const
{
  const size_t total = stats.total();
  const auto ratio = [total](size_t count)
  {
    return total == 0 ? 0.0 : 100.0 * static_cast<double>(count) / static_cast<double>(total);
  };

  Logger::user_logger->info(
    "{}: total {}, normal_positive(+N) {} ({:.2f}%), normal_negative(-N) {} ({:.2f}%), on_surface {} ({:.2f}%), unknown {} ({:.2f}%).",
    label,
    total,
    stats.outside, ratio(stats.outside),
    stats.inside, ratio(stats.inside),
    stats.on_surface, ratio(stats.on_surface),
    stats.unknown, ratio(stats.unknown));
}

void CollapseStage::do_collapse(size_t edge_num_to_collapse, size_t& total_collapsed_edge_num)
{
  auto edge_collapser = new_edge_collapser();
  edge_collapser.set_flags(
    /*update_links*/true, /*update_target_length*/false, /*update_normals*/true,
    /*check_wrinkle*/true, /*check_selfinter*/true, /*check_inter*/true);

  initialize_collapse_edges_reward();
  if (is_length_quality_priority_mode())
  {
    Logger::user_logger->info(
      "collapse priority mode: {} (submode {}, weight {}, relative ratio {}, min quality {})",
      param->priorityMode,
      param->lengthQualitySubMode,
      param->lengthQualityWeight,
      param->lengthQualityDegradationRatio,
      param->lengthQualityMinQuality);
  }
  else
    Logger::user_logger->info("collapse priority mode: {}", param->priorityMode);
  log_edge_side_stats("collapse candidates", collect_candidate_edge_side_stats());

  size_t n_vertices = rm->n_vertices();
  size_t collapsed_edge_num = 0;
  EdgeSideStats attempted_stats;
  EdgeSideStats collapsed_stats;
  while (!edges_to_collapse.empty())
  {
    auto edge_reward = edges_to_collapse.top();
    edges_to_collapse.pop();

    // out of date
    if (edge_reward.state < update_states[edge_reward.eh.idx()])
      continue;
    if (rm->status(edge_reward.eh).deleted())
      continue;

    EdgeSide edge_side = classify_edge_side(edge_reward.eh);
    add_edge_side(attempted_stats, edge_side);

    if (edge_collapser.try_collapse_edge(edge_reward.eh, edge_reward.new_point, nullptr))
    {
      add_edge_side(collapsed_stats, edge_side);
      VertexHandle center_v = edge_collapser.get_collapsed_center();
      update_after_collapsing(center_v);

      collapsed_edge_num++;
      total_collapsed_edge_num++;

      if (total_collapsed_edge_num == edge_num_to_collapse)
        break; // end collapsing

      //if (collapsed_edge_num % 1000 == 0)
      //  Logger::user_logger->info("collapsed {} edges.", collapsed_edge_num);
    }
  }
  Logger::user_logger->info("collapsed {} edges.", collapsed_edge_num);
  log_edge_side_stats("collapse attempted edges", attempted_stats);
  log_edge_side_stats("collapse succeeded edges", collapsed_stats);
  rm->garbage_collection();
  lrt->collect_garbage();
  init_one_ring_faces(rm);
}

}// namespace CageSimp
}// namespace Cage
