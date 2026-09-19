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
// Queue-only uniformity experiment. Set this to true to restore uniformity in
// placement objectives while keeping the queue formula unchanged.
constexpr bool kEnableUniformityInPhase3Solve = false;

// Square only the global target-edge ratio used by the linear-solve queue.
// Set this to false to restore the previous linear ratio without changing the
// source-adaptive uniformity path.
constexpr bool kSquareGlobalQueueUniformityScore = true;

// Include the same triangle-quality quadratic surrogate in both linear-solve
// placement and queue scoring.
constexpr bool kEnableTriangleQualityInPhase3LinearSolve = true;

// Keep uniformity out of both linear-solve placement and its line-search
// energy. The collapse queue still uses its weighted edge-length ratio.
// Retain the difference quadric for comparisons if this is re-enabled.
constexpr bool kEnableUniformityInPhase3LinearSolve = false;

double sqr(double v)
{
  return v * v;
}

// Use the same protected local length scale in energy, derivative, and solve paths.
double inverse_local_length_squared(double local_length)
{
  return 1.0 / std::max(sqr(local_length), 1e-24);
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

struct DiffUniformityTerm
{
  Vec3d neighbor;
  double target_length = 1.0;
};

struct Phase3HomogeneousQuadric
{
  Eigen::Matrix4d Q = Eigen::Matrix4d::Zero();

  void add_plane(double weight, const Vec3d& normal, double offset)
  {
    if (weight <= 0.0 || !std::isfinite(weight) || !finite_vec(normal) || !std::isfinite(offset))
      return;
    Eigen::Vector4d plane;
    plane << to_eigen(normal), offset;
    Q += weight * (plane * plane.transpose());
  }

  void add_linear_residual(double weight, const Vec3d& a, double b)
  {
    add_plane(weight, a, -b);
  }

  void add_metric_target(double weight, const Eigen::Matrix3d& metric, const Vec3d& target)
  {
    if (weight <= 0.0 || !std::isfinite(weight) || !finite_vec(target))
      return;

    const Eigen::Vector3d t = to_eigen(target);
    const Eigen::Matrix3d weighted_metric = weight * metric;
    Q.topLeftCorner<3, 3>() += weighted_metric;
    Q.topRightCorner<3, 1>() -= weighted_metric * t;
    Q.bottomLeftCorner<1, 3>() = Q.topRightCorner<3, 1>().transpose();
    Q(3, 3) += t.dot(weighted_metric * t);
  }

  void add_point_target(double weight, const Vec3d& target)
  {
    add_metric_target(weight, Eigen::Matrix3d::Identity(), target);
  }

  bool solve_unregularized(Vec3d& x) const
  {
    Eigen::Matrix3d A = Q.topLeftCorner<3, 3>();
    Eigen::Vector3d rhs = -Q.topRightCorner<3, 1>();
    A = 0.5 * (A + A.transpose());

    Eigen::FullPivLU<Eigen::Matrix3d> lu(A);
    lu.setThreshold(1e-10 * std::max(1.0, A.cwiseAbs().maxCoeff()));
    if (!lu.isInvertible())
      return false;

    const Eigen::Vector3d solved = lu.solve(rhs);
    if (!std::isfinite(solved.x()) || !std::isfinite(solved.y()) || !std::isfinite(solved.z()))
      return false;

    x = from_eigen(solved);
    return finite_vec(x);
  }

  double evaluate(const Vec3d& x) const
  {
    Eigen::Vector4d h;
    h << to_eigen(x), 1.0;
    return h.dot(Q * h);
  }

  bool solve_segment(const Vec3d& a, const Vec3d& b, Vec3d& x) const
  {
    if (!finite_vec(a) || !finite_vec(b))
      return false;

    const Eigen::Matrix3d A = 0.5 * (Q.topLeftCorner<3, 3>() + Q.topLeftCorner<3, 3>().transpose());
    const Eigen::Vector3d q = Q.topRightCorner<3, 1>();
    const Eigen::Vector3d p0 = to_eigen(a);
    const Eigen::Vector3d d = to_eigen(b - a);
    const double denom = d.dot(A * d);
    if (denom <= 1e-14 * std::max(1.0, A.cwiseAbs().maxCoeff()))
      return false;

    double t = -d.dot(A * p0 + q) / denom;
    t = clamp_value(t, 0.0, 1.0);
    x = a + (b - a) * t;
    return finite_vec(x);
  }
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
  avg_source_edge_length(0.0),
  avg_source_edge_length_initialized(false),
  phase3_target_edge_length(0.0),
  phase3_target_edge_length_initialized(false),
  phase3_target_vertices_num(0)
{}

CollapseStage::~CollapseStage()
{
  if (phase3_qem_quadric_prop_added && rm)
    rm->remove_property(phase3_qem_quadric_prop);
}

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

bool CollapseStage::is_phase3_quadric_strategy() const
{
  return is_phase3_linear_solve_strategy() || is_phase3_qem_original_strategy();
}

bool CollapseStage::is_phase3_linear_solve_strategy() const
{
  return param->phase3PlacementStrategy == "linear_solve";
}

bool CollapseStage::is_phase3_qem_original_strategy() const
{
  return param->phase3PlacementStrategy == "qem_original";
}

bool CollapseStage::is_phase3_newton_solve_strategy() const
{
  return param->phase3PlacementStrategy == "newton_solve";
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

// Build the local data needed to evaluate candidate positions for this collapse.
// No topology is changed here; the context describes the mesh patch that would
// surround the new vertex after collapsing edge eh.
CollapseStage::Phase3PlacementContext CollapseStage::make_phase3_placement_context(
  EdgeHandle eh, EdgeCollapser& edge_collapser) const
{
  Phase3PlacementContext ctx;
  ctx.edge = eh;
  ctx.halfedges = edge_collapser.get_halfedges();
  ctx.avg_cage_edge_length = avg_edge_length > 0.0 ? avg_edge_length : original_diagonal_length * 0.01;

  HalfedgeHandle heh = rm->halfedge_handle(eh, 0);
  const VertexHandle from_v = rm->from_vertex_handle(heh);
  const VertexHandle to_v = rm->to_vertex_handle(heh);
  ctx.endpoint0 = rm->point(from_v);
  ctx.endpoint1 = rm->point(to_v);
  ctx.midpoint = rm->calc_edge_midpoint(heh);
  Vec3d smooth_normal;
  edge_collapser.predict_tangential_weighted_smooth_target(
    ctx.midpoint, smooth_normal, ctx.tangential_smoothing_point);
  if (!finite_vec(ctx.tangential_smoothing_point))
    ctx.tangential_smoothing_point = ctx.midpoint;

  if (is_phase3_qem_original_strategy() && phase3_qem_quadric_prop_added)
  {
    ctx.qem_matrix = phase3_qem_quadric(from_v) + phase3_qem_quadric(to_v);
    ctx.use_qem_matrix = true;
  }

  double local_scale = rm->data(eh).edge_length;
  std::set<VertexHandle> neighbor_vertices;
  for (HalfedgeHandle h : ctx.halfedges)
  {
    const VertexHandle from_v = rm->from_vertex_handle(h);
    const VertexHandle to_v = rm->to_vertex_handle(h);
    const Vec3d& from = rm->point(from_v);
    const Vec3d& to = rm->point(to_v);
    const VertexHandle apex_v = rm->opposite_vh(h);
    const Vec3d apex = apex_v.is_valid() ? rm->point(apex_v) : ctx.midpoint;
    ctx.fan_edges.push_back({ from, to, apex });
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

  // Local mean edge length for energy normalization, tolerances, and step sizing.
  const size_t scale_count = ctx.fan_edges.empty() ? 1 : ctx.fan_edges.size() + 1;
  ctx.local_scale = std::max(local_scale / static_cast<double>(scale_count), original_diagonal_length * 1e-6);

  std::set<FaceHandle> plane_faces;
  one_ring_faces_around_edge(rm, heh, plane_faces);
  for (FaceHandle fh : plane_faces)
  {
    Vec3d pts[3];
    collect_face_points(rm, fh, pts);
    const Vec3d normal = triangle_normal(pts[0], pts[1], pts[2]);
    const double area = std::max(triangle_area_from_points(pts[0], pts[1], pts[2]), 1e-16);
    ctx.plane_constraints.push_back({ normal, -(normal | pts[0]), area });
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

double CollapseStage::evaluate_plane_energy(const Phase3PlacementContext& ctx, const Vec3d& x) const
{
  if (ctx.use_qem_matrix)
  {
    Eigen::Vector4d h;
    h << to_eigen(x), 1.0;
    const double energy = h.dot(ctx.qem_matrix * h);
    return std::isfinite(energy) ? energy : DBL_MAX;
  }

  // Normalize both face area and squared plane distance by L^2.
  const double inv_scale_fourth = sqr(inverse_local_length_squared(ctx.local_scale));
  double energy = 0.0;
  for (const PlaneConstraint& plane : ctx.plane_constraints)
  {
    const double value = (plane.normal | x) + plane.offset;
    energy += plane.weight * inv_scale_fourth * sqr(value);
  }

  if (param->curvatureMode == "weighted_plane" || param->curvatureMode == "weighted-plane")
  {
    const double scaled_curvature = ctx.edge_curvature * ctx.local_scale;
    energy *= 1.0 + sqr(scaled_curvature);
  }

  return energy;
}

// Experimental: score by plane approximation alone (ignores quality/uniformity/position
// fidelity), for testing ordering the queue purely by plane approximation
// cost while leaving the actual placement to Newton's full objective.
// Not currently wired to anything (the plane_only_score parameter on
// select_phase3_feasibility_fallback defaults to false everywhere, including
// the is_phase3_newton_solve_strategy() branch in
// compute_phase3_queue_placement_candidate) -- pass true there to re-enable.
double CollapseStage::evaluate_plane_only_energy(const Phase3PlacementContext& ctx, const Vec3d& x) const
{
  if (!finite_vec(x))
    return DBL_MAX;
  const double energy = param->planeWeight * evaluate_plane_energy(ctx, x);
  return std::isfinite(energy) ? energy : DBL_MAX;
}

void CollapseStage::initialize_phase3_qem_quadrics()
{
  if (!phase3_qem_quadric_prop_added)
  {
    rm->add_property(phase3_qem_quadric_prop, "phase3_qem_quadric");
    phase3_qem_quadric_prop_added = true;
  }

  for (VertexHandle vh : rm->vertices())
    rm->property(phase3_qem_quadric_prop, vh) = Eigen::Matrix4d::Zero();

  for (FaceHandle fh : rm->faces())
  {
    Vec3d pts[3];
    collect_face_points(rm, fh, pts);
    const double area = triangle_area_from_points(pts[0], pts[1], pts[2]);
    if (area <= 1e-16 || !std::isfinite(area))
      continue;

    const Vec3d normal = triangle_normal(pts[0], pts[1], pts[2]);
    if (!finite_vec(normal))
      continue;

    Eigen::Vector4d plane;
    plane << to_eigen(normal), -(normal | pts[0]);
    const Eigen::Matrix4d fundamental_quadric = plane * plane.transpose();

    for (VertexHandle vh : rm->fv_range(fh))
      rm->property(phase3_qem_quadric_prop, vh) += fundamental_quadric;
  }
}

Eigen::Matrix4d CollapseStage::phase3_qem_quadric(VertexHandle vh) const
{
  if (!phase3_qem_quadric_prop_added || !vh.is_valid() || rm->status(vh).deleted())
    return Eigen::Matrix4d::Zero();
  return rm->property(phase3_qem_quadric_prop, vh);
}

void CollapseStage::update_phase3_qem_quadric_after_collapse(
  VertexHandle center_vh, const Eigen::Matrix4d& quadric)
{
  if (!phase3_qem_quadric_prop_added || !center_vh.is_valid() || rm->status(center_vh).deleted())
    return;
  rm->property(phase3_qem_quadric_prop, center_vh) = quadric;
}

double CollapseStage::evaluate_triangle_quality_energy(const Phase3PlacementContext& ctx, const Vec3d& x) const
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

Eigen::Matrix4d CollapseStage::build_phase3_triangle_quality_surrogate_quadric(
  const Phase3PlacementContext& ctx) const
{
  param->validate_triangle_shape_settings();
  Phase3HomogeneousQuadric surrogate;
  if (ctx.fan_edges.empty())
    return surrogate.Q;

  // Average directional residual costs normalized by L^2, with equal face weights.
  // The caller applies triangleQualityWeight to placement and queue scoring.
  const double fan_count = static_cast<double>(ctx.fan_edges.size());
  const double w = inverse_local_length_squared(ctx.local_scale) / fan_count;

  for (const CollapseFanEdge& fan_edge : ctx.fan_edges)
  {
    const Vec3d base = fan_edge.to - fan_edge.from;
    const double base_length = base.length();
    if (base_length <= ctx.local_scale * 1e-8 ||
      !std::isfinite(base_length))
      continue;

    const Vec3d t = base / base_length;
    const Vec3d mid = (fan_edge.from + fan_edge.to) * 0.5;
    Vec3d s = fan_edge.apex - mid;
    s -= (s | t) * t;
    s = normalized_or_fallback(
      s,
      triangle_normal(
        fan_edge.from, fan_edge.to, fan_edge.apex).cross(t));
    s = normalized_or_fallback(s, Vec3d(0.0, 1.0, 0.0));

    Vec3d n = t.cross(s);
    if (!finite_vec(n) || n.length() <= 1e-12 || std::abs(t | s) > 1e-8)
    {
      // Recover an orthogonal frame only for degenerate or rounded input.
      // The least-aligned coordinate axis has a nonzero projection off t.
      const Vec3d axis = std::abs(t.x()) <= std::abs(t.y()) && std::abs(t.x()) <= std::abs(t.z())
        ? Vec3d(1.0, 0.0, 0.0)
        : (std::abs(t.y()) <= std::abs(t.z()) ? Vec3d(0.0, 1.0, 0.0) : Vec3d(0.0, 0.0, 1.0));
      s = axis - (axis | t) * t;
      s = s / s.norm();
      n = t.cross(s);
    }
    n = normalized_or_fallback(
      n, triangle_normal(fan_edge.from, fan_edge.to, fan_edge.apex));

    const double target_height =
      0.5 * std::sqrt(3.0) * base_length;

    // The fixed frame gives r_t=t.(x-mid), r_s=s.(x-mid)-height,
    // and r_n=n.(x-mid). The homogeneous builder retains their constants.
    surrogate.add_linear_residual(
      w * param->triangleShapeTangentWeight, t, t | mid);
    surrogate.add_linear_residual(
      w * param->triangleShapeHeightWeight, s, (s | mid) + target_height);
    surrogate.add_linear_residual(
      w * param->triangleShapeNormalWeight, n, n | mid);
  }

  return surrogate.Q;
}

double CollapseStage::evaluate_phase3_triangle_quality_surrogate(
  const Phase3PlacementContext& ctx, const Vec3d& x) const
{
  if (!finite_vec(x))
    return DBL_MAX;

  Eigen::Vector4d h;
  h << to_eigen(x), 1.0;
  const Eigen::Matrix4d surrogate =
    build_phase3_triangle_quality_surrogate_quadric(ctx);
  const double energy = h.dot(surrogate * h);
  if (!std::isfinite(energy))
    return DBL_MAX;
  return std::max(energy, 0.0);
}

// Uniformity as an exactly linear system, so it can join the plane and triangle
// terms in the same single linear solve.
//
// The target is |x - n_i| = L_i for every one-ring neighbor n_i. Squaring it
// gives |x|^2 - 2 n_i.x + (|n_i|^2 - L_i^2) = 0, whose only nonlinear term,
// |x|^2, is identical in every neighbor equation. Subtracting the mean
// equation therefore cancels it exactly and leaves
//
//   2 (n_mean - n_i).x = c_mean - c_i,    c_i = |n_i|^2 - L_i^2
//
// which is not a linearization but an algebraically equivalent restatement.
// Differencing against the mean, rather than against one arbitrary reference
// neighbor, keeps the rows balanced.
//
// The rows span {n_i - n_mean}, essentially the tangent plane of the one-ring,
// so this term positions x tangentially and leaves the normal direction to plane approximation
// instead of competing with it. In "global" mode every L_i is the same
// constant, the target-length term cancels as well, and the rows reduce to
// "x is equidistant from all its neighbors".
//
// Each row is normalized to a unit direction so its residual is a signed
// distance, then scaled by 1/local_scale^2. The caller applies uniformityWeight.
Eigen::Matrix4d CollapseStage::build_phase3_uniformity_difference_quadric(
  const Phase3PlacementContext& ctx) const
{
  Phase3HomogeneousQuadric uniformity;
  const size_t neighbor_count = ctx.neighbor_points.size();
  if (neighbor_count < 2)
    return uniformity.Q;

  std::vector<double> offsets;
  offsets.reserve(neighbor_count);
  Vec3d mean_neighbor(0.0, 0.0, 0.0);
  double mean_offset = 0.0;
  for (const Vec3d& neighbor : ctx.neighbor_points)
  {
    const double target_length =
      source_uniformity_target_length(ctx, (ctx.midpoint + neighbor) * 0.5);
    const double offset = neighbor.sqrnorm() - sqr(target_length);
    offsets.push_back(offset);
    mean_neighbor += neighbor;
    mean_offset += offset;
  }
  mean_neighbor /= static_cast<double>(neighbor_count);
  mean_offset /= static_cast<double>(neighbor_count);

  const double inv_scale_sqr = 1.0 / std::max(sqr(ctx.local_scale), 1e-24);
  const double direction_epsilon = ctx.local_scale * 1e-8;
  const double w = inv_scale_sqr / static_cast<double>(neighbor_count);
  for (size_t i = 0; i < neighbor_count; i++)
  {
    Vec3d direction = (mean_neighbor - ctx.neighbor_points[i]) * 2.0;
    double rhs = mean_offset - offsets[i];
    const double direction_length = direction.length();
    // A neighbor sitting on the one-ring centroid carries no direction.
    if (!(direction_length > direction_epsilon) ||
      !std::isfinite(direction_length) || !std::isfinite(rhs))
      continue;
    direction /= direction_length;
    rhs /= direction_length;

    uniformity.add_linear_residual(w, direction, rhs);
  }

  return uniformity.Q;
}

double CollapseStage::evaluate_phase3_uniformity_difference(
  const Phase3PlacementContext& ctx, const Vec3d& x) const
{
  if (!finite_vec(x))
    return DBL_MAX;

  Eigen::Vector4d h;
  h << to_eigen(x), 1.0;
  const Eigen::Matrix4d uniformity =
    build_phase3_uniformity_difference_quadric(ctx);
  const double energy = h.dot(uniformity * h);
  if (!std::isfinite(energy))
    return DBL_MAX;
  return std::max(energy, 0.0);
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

double CollapseStage::evaluate_pre_collapse_triangle_quality_penalty(
  const Phase3PlacementContext& ctx) const
{
  if (!ctx.edge.is_valid() || rm->status(ctx.edge).deleted())
    return 0.0;

  std::set<FaceHandle> faces;
  one_ring_faces_around_edge(rm, rm->halfedge_handle(ctx.edge, 0), faces);
  if (faces.empty())
    return 0.0;

  double penalty_sum = 0.0;
  for (FaceHandle fh : faces)
  {
    const double quality = clamp_value(calc_triangle_quality(rm, fh), 0.0, 1.0);
    penalty_sum += sqr(1.0 - quality);
  }
  return clamp_value(
    penalty_sum / static_cast<double>(faces.size()), 0.0, 1.0);
}

double CollapseStage::source_edge_length_at(const Vec3d& sample_point) const
{
  const double fallback_length =
    avg_source_edge_length > 0.0 ? avg_source_edge_length : avg_edge_length;
  FaceHandle closest_face;
#ifdef USE_TREE_SEARCH
  if (!ot)
    return fallback_length;
  closest_face = FaceHandle(ot->closest_point(sample_point).second);
#else
  if (!og)
    return fallback_length;
  closest_face = FaceHandle(og->closest_point(sample_point).second);
#endif

  if (!closest_face.is_valid() || closest_face.idx() < 0 ||
    closest_face.idx() >= static_cast<int>(om->n_faces()))
    return fallback_length;

  const double source_edge_length = face_average_edge_length(om, closest_face);
  return source_edge_length > 0.0 ? source_edge_length : fallback_length;
}

double CollapseStage::evaluate_uniformity_energy(const Phase3PlacementContext& ctx, const Vec3d& x) const
{
  if (param->uniformityWeight <= 0.0 || param->uniformityMode == "none")
    return 0.0;

  double energy = 0.0;
  for (const Vec3d& neighbor : ctx.neighbor_points)
  {
    const double target_length = source_uniformity_target_length(ctx, (x + neighbor) * 0.5);
    const double relative_overshoot =
      ((x - neighbor).length() - target_length) / target_length;
    energy += sqr(std::max(relative_overshoot, 0.0));
  }
  return energy;
}

double CollapseStage::evaluate_newton_energy(const Phase3PlacementContext& ctx, const Vec3d& x) const
{
  if (!finite_vec(x))
    return DBL_MAX;

  // Normalize quality/uniformity by fan/neighbor count, matching
  // evaluate_phase3_proxy_energy's convention, so a given weight means the
  // same relative strength (O(1) per vertex, not growing with valence)
  // regardless of which placement strategy (linear solve vs Newton) is in use.
  const double fan_count = static_cast<double>(std::max<size_t>(ctx.fan_edges.size(), 1));
  const double neighbor_count = static_cast<double>(ctx.neighbor_points.size());

  double energy = 0.0;
  energy += param->planeWeight * evaluate_plane_energy(ctx, x);
  if (kEnableTriangleQualityInPhase3LinearSolve &&
    is_phase3_linear_solve_strategy())
  {
    energy += param->triangleQualityWeight *
      evaluate_phase3_triangle_quality_surrogate(ctx, x);
  }
  else if (!is_phase3_quadric_strategy())
  {
    energy += param->triangleQualityWeight * evaluate_triangle_quality_energy(ctx, x) / fan_count;
  }
  if (kEnableUniformityInPhase3LinearSolve &&
    is_phase3_linear_solve_strategy())
  {
    // The same quadratic the linear solve minimizes, so Armijo backtracking
    // does not reject the solved point for optimizing a different objective.
    energy += param->uniformityWeight *
      evaluate_phase3_uniformity_difference(ctx, x);
  }
  else if (kEnableUniformityInPhase3Solve)
    energy += param->uniformityWeight * evaluate_uniformity_energy(ctx, x) / neighbor_count;
  return std::isfinite(energy) ? energy : DBL_MAX;
}

CollapseStage::NewtonDerivatives CollapseStage::finite_difference_newton_derivatives(
  const Phase3PlacementContext& ctx, const Vec3d& x) const
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
  const Phase3PlacementContext& ctx, const Vec3d& x) const
{
  return autodiff_newton_derivatives(ctx, x);
}

CollapseStage::NewtonDerivatives CollapseStage::autodiff_newton_derivatives(
  const Phase3PlacementContext& ctx, const Vec3d& x) const
{
  NewtonDerivatives deriv;
  deriv.energy = evaluate_newton_energy(ctx, x);
  if (!std::isfinite(deriv.energy))
    return deriv;

  std::vector<DiffUniformityTerm> uniformity_terms;

  if (kEnableUniformityInPhase3Solve &&
    param->uniformityWeight > 0.0 && param->uniformityMode != "none")
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

  const double inv_scale_fourth = sqr(inverse_local_length_squared(ctx.local_scale));
  // Match evaluate_newton_energy's weights, so the gradient/Hessian
  // used for the actual Newton step is consistent with the energy value
  // (deriv.energy, set above from evaluate_newton_energy) used for
  // acceptance/line-search comparisons.
  const double fan_count = static_cast<double>(std::max<size_t>(ctx.fan_edges.size(), 1));
  const double neighbor_count = static_cast<double>(ctx.neighbor_points.size());
  double curvature_multiplier = 1.0;
  if (param->curvatureMode == "weighted_plane" || param->curvatureMode == "weighted-plane")
  {
    const double scaled_curvature = ctx.edge_curvature * ctx.local_scale;
    curvature_multiplier += sqr(scaled_curvature);
  }

  for (const PlaneConstraint& plane : ctx.plane_constraints)
  {
    const ADScalar value =
      plane.normal.x() * ad_point.x +
      plane.normal.y() * ad_point.y +
      plane.normal.z() * ad_point.z +
      plane.offset;
    const double weight =
      param->planeWeight * curvature_multiplier * plane.weight * inv_scale_fourth;
    ad_energy += weight * diff_sqr(value);
  }

  if (kEnableTriangleQualityInPhase3LinearSolve || !is_phase3_quadric_strategy())
  {
    for (const CollapseFanEdge& fan_edge : ctx.fan_edges)
    {
      const ADScalar quality = diff_clamp(diff_triangle_quality(
        diff_vec<ADScalar>(fan_edge.from),
        diff_vec<ADScalar>(fan_edge.to),
        ad_point), 0.0, 1.0);
      ad_energy += param->triangleQualityWeight * diff_sqr(ADScalar(1.0) - quality) / fan_count;
    }
  }

  for (const DiffUniformityTerm& term : uniformity_terms)
  {
    const double target_length = term.target_length;
    const ADScalar length = diff_norm(ad_point - diff_vec<ADScalar>(term.neighbor));
    ADScalar relative_overshoot = (length - target_length) / target_length;
    if (relative_overshoot < 0.0)
      relative_overshoot = 0.0;
    ad_energy += param->uniformityWeight *
      diff_sqr(relative_overshoot) / neighbor_count;
  }

  Eigen::Vector3d gradient = Eigen::Vector3d::Zero();
  Eigen::Matrix3d hessian = autodiff::hessian(ad_energy, ad_x, gradient);
  hessian = 0.5 * (hessian + hessian.transpose());

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

// Objective used for candidate comparison and line search. Queue ordering uses
// evaluate_phase3_queue_score after the candidate position has been selected.
double CollapseStage::evaluate_phase3_proxy_energy(const Phase3PlacementContext& ctx, const Vec3d& x) const
{
  return evaluate_newton_energy(ctx, x);
}

// Compute per-edge queue components. The area-weighted plane approximation and triangle
// surrogate use dimensionless local costs, normalized by L^4 and L^2.
CollapseStage::Phase3QueueComponents CollapseStage::evaluate_phase3_queue_components(
  const Phase3PlacementContext& ctx, const Vec3d& x) const
{
  Phase3QueueComponents components;
  if (!finite_vec(x))
  {
    components.plane = DBL_MAX;
    return components;
  }

  if (param->planeWeight > 0.0)
    components.plane = evaluate_plane_energy(ctx, x);

  const double fan_count =
    static_cast<double>(std::max<size_t>(ctx.fan_edges.size(), 1));

  if (param->triangleQualityWeight > 0.0)
  {
    if (kEnableTriangleQualityInPhase3LinearSolve &&
      is_phase3_linear_solve_strategy())
    {
      // Evaluate exactly the same homogeneous (A, b, c) surrogate used by
    // the Phase 3 linear solve at the selected queue point.
      components.triangle_quality =
        evaluate_phase3_triangle_quality_surrogate(ctx, x);
    }
    else
    {
      const double pre_collapse_penalty =
        evaluate_pre_collapse_triangle_quality_penalty(ctx);
      const double post_collapse_penalty = clamp_value(
        evaluate_triangle_quality_energy(ctx, x) / fan_count, 0.0, 1.0);

      const double quality_improvement =
        pre_collapse_penalty - post_collapse_penalty;
      components.triangle_quality =
        clamp_value(0.5 * (1.0 - quality_improvement), 0.0, 1.0);
    }
  }

  if (param->uniformityWeight > 0.0 && param->uniformityMode != "none")
  {
    // Edge-size ratio used only for queue ordering:
    //   global mode: queue ratio = current cage edge / estimated target mean
    //   source mode: preserve the source-adaptive target that scales the cage
    //                global mean by the local/source mean-edge ratio
    // A smaller ratio gives this edge a better (smaller) queue score.
    double target_length = phase3_target_edge_length;
    if (param->uniformityMode != "global" ||
      target_length <= 0.0 || !std::isfinite(target_length))
      target_length = source_uniformity_target_length(ctx, ctx.midpoint);
    const double cage_edge_length =
      (ctx.endpoint1 - ctx.endpoint0).length();
    const double edge_length_ratio =
      cage_edge_length / std::max(target_length, 1e-24);
    components.uniformity =
      param->uniformityMode == "global" &&
      kSquareGlobalQueueUniformityScore
      ? sqr(edge_length_ratio)
      : edge_length_ratio;
  }

  return components;
}

// Combine active components as a weighted sum. Plane approximation and the linear-solve
// triangle surrogate are intentionally not pass-normalized.
double CollapseStage::evaluate_phase3_queue_score(
  const Phase3QueueComponents& components) const
{
  double weighted_score = 0.0;
  const auto normalized_by_pass_max = [](double value, double maximum)
  {
    if (!std::isfinite(value))
      return 1.0;
    if (maximum <= 1e-24)
      return 0.0;
    return clamp_value(value / maximum, 0.0, 1.0);
  };
  const auto add_component = [&](double weight, double score)
  {
    if (weight <= 0.0)
      return;
    weighted_score += weight * clamp_value(score, 0.0, 1.0);
  };
  const auto add_nonnegative_component = [&](double weight, double score)
  {
    if (weight <= 0.0)
      return;
    weighted_score += weight * std::max(score, 0.0);
  };

  add_nonnegative_component(param->planeWeight, components.plane);

  if (kEnableTriangleQualityInPhase3LinearSolve &&
    is_phase3_linear_solve_strategy())
  {
    // Use the mean weighted directional residuals, normalized by L^2.
    // triangleQualityWeight controls its contribution to the queue cost.
    add_nonnegative_component(
      param->triangleQualityWeight,
      components.triangle_quality);
  }
  else
  {
    add_component(param->triangleQualityWeight, components.triangle_quality);
  }

  if (param->uniformityMode != "none")
  {
    // The global component is already dimensionless relative to the fixed
    // target length. Do not divide it by the pass maximum: that would cancel
    // the target length and reduce the component back to edge/max-edge.
    if (param->uniformityMode == "global")
      add_nonnegative_component(param->uniformityWeight, components.uniformity);
    else
    {
      // Preserve the previous source-adaptive queue normalization.
      add_component(
        param->uniformityWeight,
        normalized_by_pass_max(
          components.uniformity,
          phase3_queue_component_maxima.uniformity));
    }
  }

  return weighted_score;
}

double CollapseStage::evaluate_phase3_refinement_residual(const Phase3PlacementContext& ctx, const Vec3d& x) const
{
  if (!finite_vec(x))
    return DBL_MAX;

  double residual = 0.0;
  if (kEnableUniformityInPhase3Solve &&
    param->uniformityWeight > 0.0 && param->uniformityMode != "none")
  {
    const double denom = static_cast<double>(ctx.neighbor_points.size());
    residual += param->uniformityWeight * evaluate_uniformity_energy(ctx, x) /
      denom;
  }

  return std::isfinite(residual) ? residual : DBL_MAX;
}

double CollapseStage::calc_phase3_fan_min_quality(const Phase3PlacementContext& ctx, const Vec3d& x) const
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

double CollapseStage::source_uniformity_target_length(const Phase3PlacementContext& ctx, const Vec3d& sample_point) const
{
  double target_length = ctx.avg_cage_edge_length;
  if (param->uniformityMode == "source")
  {
    const double source_edge_length = source_edge_length_at(sample_point);
    if (avg_source_edge_length > 0.0 &&
      source_edge_length > 0.0 && std::isfinite(source_edge_length))
    {
      const double ratio = source_edge_length / avg_source_edge_length;
      target_length *= ratio;
    }
  }

  return target_length;
}

// Hard validity test for a proposed new vertex position. Collision checks are
// disabled by the qem_original strategy; topology/degeneracy guards remain.
bool CollapseStage::phase3_placement_satisfies_hard_constraints(
  const Phase3PlacementContext& ctx, EdgeCollapser& edge_collapser, const Vec3d& x) const
{
  UNUSED(ctx);
  if (!collapse_target_valid(edge_collapser, x))
    return false;
  return true;
}

void CollapseStage::set_phase3_edge_collapser_flags(EdgeCollapser& edge_collapser) const
{
  const bool check_collision = !is_phase3_qem_original_strategy();
  edge_collapser.set_flags(
    /*update_links*/false, /*update_target_length*/false, /*update_normals*/true,
    /*check_face_orientation_violation*/check_collision, /*check_selfinter*/check_collision,
    /*check_inter*/check_collision);
}

bool CollapseStage::select_phase3_feasibility_fallback(
  EdgeHandle eh, const Phase3PlacementContext& ctx, EdgeCollapser& edge_collapser,
  Vec3d& new_point, double& energy, bool plane_only_score) const
{
  std::vector<Vec3d> candidates;
  const double duplicate_tol = std::max(ctx.local_scale * 1e-8, original_diagonal_length * 1e-12);
  const auto append_candidate = [&](const Vec3d& p)
  {
    if (!finite_vec(p))
      return;
    for (const Vec3d& existing : candidates)
    {
      if ((existing - p).length() <= duplicate_tol)
        return;
    }
    candidates.push_back(p);
  };

  if (is_phase3_qem_original_strategy())
  {
    append_candidate(ctx.endpoint0);
    append_candidate(ctx.endpoint1);
    append_candidate(ctx.midpoint);
  }
  else
  {
    append_candidate(ctx.midpoint);

    if (eh.is_valid())
    {
      const HalfedgeHandle heh = rm->halfedge_handle(eh, 0);
      if (heh.is_valid())
      {
        append_candidate(rm->point(rm->to_vertex_handle(heh)));
        append_candidate(rm->point(rm->from_vertex_handle(heh)));
      }
    }
  }

  bool found_valid_candidate = false;
  energy = DBL_MAX;
  for (const Vec3d& candidate : candidates)
  {
    if (!phase3_placement_satisfies_hard_constraints(ctx, edge_collapser, candidate))
      continue;

    const double candidate_energy = plane_only_score ?
      evaluate_plane_only_energy(ctx, candidate) : evaluate_phase3_proxy_energy(ctx, candidate);
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

bool CollapseStage::select_phase3_linear_solve_line_search_candidate(
  const Phase3PlacementContext& ctx, EdgeCollapser& edge_collapser,
  const Vec3d& quadric_point, Vec3d& selected_point, double& selected_energy) const
{
  if (!finite_vec(quadric_point))
    return false;

  // Build the Botsch-Kobbelt area-equalizing tangential smoothing endpoint
  // from the raw linear-system solution, not from the edge midpoint. Armijo
  // backtracking then moves from that endpoint toward the raw solution; its
  // alpha supplies the effective damping for the smoothing displacement.
  Vec3d smooth_normal;
  Vec3d base;
  if (!edge_collapser.predict_area_equalizing_tangential_smooth_target(
    quadric_point, smooth_normal, base, /*damping*/1.0) || !finite_vec(base))
    return false;

  const Vec3d direction = quadric_point - base;
  const double direction_length = direction.length();
  const double base_energy = evaluate_phase3_proxy_energy(ctx, base);
  if (!std::isfinite(base_energy))
    return false;

  if (direction_length <= 1e-14)
  {
    if (!phase3_placement_satisfies_hard_constraints(ctx, edge_collapser, base))
      return false;
    selected_point = base;
    selected_energy = base_energy;
    return true;
  }

  // Approximate the directional derivative at the base because the closed-form
  // Quadric solve does not provide the gradient used by the Newton line search.
  const double h = std::max(direction_length * 1e-4, ctx.local_scale * 1e-6);
  const Vec3d probe = base + (h / direction_length) * direction;
  const double probe_energy = evaluate_phase3_proxy_energy(ctx, probe);
  const double descent = std::isfinite(probe_energy) ?
    std::max(1e-16, (base_energy - probe_energy) / h * direction_length) : 1e-16;

  double alpha = 1.0;
  const size_t max_iter = std::max<size_t>(param->lineSearchMaxIter, 1);
  for (size_t ls = 0; ls < max_iter; ls++)
  {
    const Vec3d trial = base + alpha * direction;
    const double trial_energy = evaluate_phase3_proxy_energy(ctx, trial);
    if (std::isfinite(trial_energy) &&
      trial_energy <= base_energy - 1e-4 * alpha * descent &&
      phase3_placement_satisfies_hard_constraints(ctx, edge_collapser, trial))
    {
      selected_point = trial;
      selected_energy = trial_energy;
      return true;
    }
    alpha *= 0.5;
  }

  return false;
}

// Select the new vertex position for one edge collapse before the topology
// change is committed. The policy is controlled by phase3PlacementStrategy:
// linear_solve, newton_solve, or qem_original.
bool CollapseStage::choose_phase3_collapse_placement(
  EdgeHandle eh, EdgeCollapser& edge_collapser, const Vec3d& queued_point,
  Phase3PlacementDecision& decision, double& newton_seconds)
{
  decision = Phase3PlacementDecision();
  const Phase3PlacementContext ctx = make_phase3_placement_context(eh, edge_collapser);
  // A mixed edge has exactly one placement: the existing rail vertex.
  // Evaluate that point directly, without a solve, line search or projection.
  if (edge_collapser.has_fixed_rail_target())
  {
    decision.point = edge_collapser.constrained_target_point(ctx.midpoint);
    if (!phase3_placement_satisfies_hard_constraints(ctx, edge_collapser, decision.point))
      return false;
    decision.priority_energy = evaluate_phase3_proxy_energy(ctx, decision.point);
    return std::isfinite(decision.priority_energy);
  }
  if (is_phase3_qem_original_strategy() &&
    !ctx.use_qem_matrix && ctx.plane_constraints.empty())
    return false;

  std::vector<Vec3d> candidates;
  std::vector<bool> is_fallback;
  const double duplicate_tol = std::max(ctx.local_scale * 1e-8, original_diagonal_length * 1e-12);
  const auto append_candidate = [&](const Vec3d& p, bool fallback)
  {
    const Vec3d constrained = edge_collapser.constrained_target_point(p);
    if (!finite_vec(constrained))
      return;
    for (const Vec3d& existing : candidates)
    {
      if ((existing - constrained).length() <= duplicate_tol)
        return;
    }
    candidates.push_back(constrained);
    is_fallback.push_back(fallback);
  };

  const bool use_quadric_solve = is_phase3_quadric_strategy();
  append_candidate(queued_point, false);

  if (use_quadric_solve)
  {
    Vec3d quadric_point;
    double quadric_energy = DBL_MAX;
    if (solve_phase3_quadric_placement(ctx, quadric_point, quadric_energy))
    {
      if (is_phase3_qem_original_strategy())
      {
        append_candidate(quadric_point, false);
      }
      else if (param->phase3LinearSolveCollisionReject)
      {
        // Opt-in alternative to the backtracking line search below: if the
        // raw linear-solve point already fails hard validity (collision,
        // degenerate, wrinkle), reject this edge outright instead of
        // searching for a nearby fallback placement.
        if (!phase3_placement_satisfies_hard_constraints(ctx, edge_collapser, quadric_point))
          return false;
        append_candidate(quadric_point, false);
      }
      else
      {
        Vec3d line_search_point;
        double line_search_energy = DBL_MAX;
        if (select_phase3_linear_solve_line_search_candidate(
          ctx, edge_collapser, quadric_point, line_search_point, line_search_energy))
          append_candidate(line_search_point, false);
      }
    }
  }

  if (!use_quadric_solve)
  {
    append_candidate(ctx.midpoint, true);
    const HalfedgeHandle heh = rm->halfedge_handle(eh, 0);
    append_candidate(rm->point(rm->to_vertex_handle(heh)), true);
    append_candidate(rm->point(rm->from_vertex_handle(heh)), true);
  }

  bool found_valid_candidate = false;
  for (size_t i = 0; i < candidates.size(); i++)
  {
    if (!phase3_placement_satisfies_hard_constraints(ctx, edge_collapser, candidates[i]))
      continue;

    const double candidate_energy = evaluate_phase3_proxy_energy(ctx, candidates[i]);
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

  if (!found_valid_candidate && use_quadric_solve)
  {
    Vec3d fallback_point;
    double fallback_energy = DBL_MAX;
    if (select_phase3_feasibility_fallback(eh, ctx, edge_collapser, fallback_point, fallback_energy))
    {
      decision.point = fallback_point;
      decision.priority_energy = fallback_energy;
      decision.used_fallback = true;
      found_valid_candidate = true;
    }
  }

  if (!found_valid_candidate)
    return false;

  if (!is_phase3_newton_solve_strategy())
    return true;

  decision.min_quality = calc_phase3_fan_min_quality(ctx, decision.point);
  decision.nonlinear_residual = evaluate_phase3_refinement_residual(ctx, decision.point);

  decision.attempted_newton = true;
  Vec3d refined_point;
  double refined_energy = DBL_MAX;
  const auto newton_start_time = std::chrono::steady_clock::now();
  const bool newton_ok = refine_collapse_placement_with_newton(
    eh, edge_collapser, refined_point, refined_energy, &decision.point);
  newton_seconds += std::chrono::duration<double>(
    std::chrono::steady_clock::now() - newton_start_time).count();

  if (newton_ok)
    refined_point = edge_collapser.constrained_target_point(refined_point);
  if (newton_ok && phase3_placement_satisfies_hard_constraints(ctx, edge_collapser, refined_point))
  {
    const double current_energy = evaluate_newton_energy(ctx, decision.point);
    const double refined_quality = calc_phase3_fan_min_quality(ctx, refined_point);
    const double refined_residual = evaluate_phase3_refinement_residual(ctx, refined_point);
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
// not been collapsed yet; topology changes happen later in do_phase3_energy_simplification().
bool CollapseStage::refine_collapse_placement_with_newton(
  EdgeHandle eh, EdgeCollapser& edge_collapser, Vec3d& new_point, double& energy,
  const Vec3d* initial_point)
{
  const Phase3PlacementContext ctx = make_phase3_placement_context(eh, edge_collapser);
  Vec3d x = initial_point && finite_vec(*initial_point) ?
    *initial_point : ctx.tangential_smoothing_point;
  if (!finite_vec(x))
    x = ctx.midpoint;

  if (!is_exact_reject_robustness_mode() && !collapse_target_valid(edge_collapser, x))
  {
    if (collapse_target_valid(edge_collapser, ctx.tangential_smoothing_point))
      x = ctx.tangential_smoothing_point;
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

bool CollapseStage::solve_phase3_quadric_placement(
  const Phase3PlacementContext& ctx, Vec3d& new_point, double& energy) const
{
  // An edge without plane constraints is not hopeless: Garland-Heckbert already
  // prescribes falling back to the segment optimum and then to the midpoint,
  // and the triangle-quality surrogate added below can still shape the solve.
  // Only pure QEM, which has no other term to solve with, gives up here.
  if (is_phase3_qem_original_strategy() &&
    !ctx.use_qem_matrix && ctx.plane_constraints.empty())
    return false;

  double curvature_multiplier = 1.0;
  if (param->curvatureMode == "weighted_plane" || param->curvatureMode == "weighted-plane")
  {
    const double scaled_curvature = ctx.edge_curvature * ctx.local_scale;
    curvature_multiplier += sqr(scaled_curvature);
  }

  Phase3HomogeneousQuadric quadric;
  const double inv_scale_fourth = sqr(inverse_local_length_squared(ctx.local_scale));
  if (ctx.use_qem_matrix)
  {
    quadric.Q = ctx.qem_matrix;
  }
  else for (const PlaneConstraint& plane : ctx.plane_constraints)
  {
    const double w = param->planeWeight * curvature_multiplier * plane.weight * inv_scale_fourth;
    quadric.add_plane(w, plane.normal, plane.offset);
  }

  const bool pure_qem = is_phase3_qem_original_strategy();

  if ((kEnableTriangleQualityInPhase3LinearSolve || !is_phase3_quadric_strategy()) &&
    !pure_qem && param->triangleQualityWeight > 0.0 && !ctx.fan_edges.empty())
  {
    quadric.Q += param->triangleQualityWeight *
      build_phase3_triangle_quality_surrogate_quadric(ctx);
  }

  if (kEnableUniformityInPhase3LinearSolve && is_phase3_linear_solve_strategy() &&
    param->uniformityWeight > 0.0 && param->uniformityMode != "none")
  {
    quadric.Q += param->uniformityWeight *
      build_phase3_uniformity_difference_quadric(ctx);
  }
  else if (kEnableUniformityInPhase3Solve && !pure_qem &&
    param->uniformityWeight > 0.0 && param->uniformityMode != "none" &&
    !ctx.neighbor_points.empty())
  {
    const double neighbor_count = static_cast<double>(ctx.neighbor_points.size());
    for (const Vec3d& neighbor : ctx.neighbor_points)
    {
      const double target_length = source_uniformity_target_length(ctx, (ctx.midpoint + neighbor) * 0.5);

      Vec3d direction =
        normalized_or_fallback(ctx.midpoint - neighbor, ctx.tangential_smoothing_point - neighbor);
      direction = normalized_or_fallback(direction, Vec3d(1.0, 0.0, 0.0));

      const double w = param->uniformityWeight /
        neighbor_count / std::max(sqr(target_length), 1e-24);
      quadric.add_linear_residual(w, direction, (direction | neighbor) + target_length);
    }
  }

  const bool solved = quadric.solve_unregularized(new_point);
  if (!solved)
  {
    // Garland-Heckbert singular-quadric cascade (Sec. 4): first the optimum
    // constrained to the collapsed segment, then the midpoint and endpoints.
    // The midpoint leads the candidate list and seeds new_point so that a
    // quadric that is identically zero keeps the midpoint rather than an
    // arbitrary endpoint.
    if (!quadric.solve_segment(ctx.endpoint0, ctx.endpoint1, new_point))
    {
      const Vec3d candidates[3] = { ctx.midpoint, ctx.endpoint0, ctx.endpoint1 };
      double best_energy = DBL_MAX;
      new_point = ctx.midpoint;
      for (const Vec3d& candidate : candidates)
      {
        const double candidate_energy = quadric.evaluate(candidate);
        if (finite_vec(candidate) && std::isfinite(candidate_energy) && candidate_energy < best_energy)
        {
          best_energy = candidate_energy;
          new_point = candidate;
        }
      }
      if (!finite_vec(new_point))
        new_point = ctx.midpoint;
    }
  }

  // This energy compares placement candidates. Queue priority is recomputed
  // from the selected point with normalized [0, 1] components.
  energy = evaluate_phase3_proxy_energy(ctx, new_point);
  return finite_vec(new_point) && std::isfinite(energy);
}

// Compute the cheap placement used only to rank this edge in the Phase 3 queue.
// Keep queue ordering on the stable quadric placement when possible. Pure QEM uses
// only Garland-Heckbert fallback points if the quadric solve is singular or
// unusable for the mesh collapse.
bool CollapseStage::compute_phase3_queue_placement_candidate(EdgeHandle eh, Vec3d& new_point, double& energy)
{
  auto edge_collapser = new_edge_collapser();
  set_phase3_edge_collapser_flags(edge_collapser);
  if (!edge_collapser.init(eh))
    return false;

  HalfedgeHandle heh = rm->halfedge_handle(eh, 0);
  size_t valence_after_collapsing =
    (rm->valence(rm->to_vertex_handle(heh)) + rm->valence(rm->from_vertex_handle(heh))) - 3;
  // Both quadric-based strategies choose the collapsed position by solving an
  // energy instead of keeping the mesh near-regular, so the valence cap only
  // hides candidates they could still place. Newton placement keeps it.
  if (!is_phase3_quadric_strategy() &&
    valence_after_collapsing > param->maxValence)
    return false;

  const Phase3PlacementContext ctx = make_phase3_placement_context(eh, edge_collapser);
  if (edge_collapser.has_fixed_rail_target())
  {
    new_point = edge_collapser.constrained_target_point(ctx.midpoint);
    if (!phase3_placement_satisfies_hard_constraints(ctx, edge_collapser, new_point))
      return false;
    energy = evaluate_phase3_proxy_energy(ctx, new_point);
    return std::isfinite(energy);
  }
  if (is_phase3_qem_original_strategy() &&
    !ctx.use_qem_matrix && ctx.plane_constraints.empty())
    return false;

  if (is_phase3_newton_solve_strategy())
    return select_phase3_feasibility_fallback(eh, ctx, edge_collapser, new_point, energy, /*plane_only_score*/false);

  if (!solve_phase3_quadric_placement(ctx, new_point, energy))
    return select_phase3_feasibility_fallback(eh, ctx, edge_collapser, new_point, energy);

  if (is_phase3_linear_solve_strategy())
  {
    Vec3d line_search_point;
    double line_search_energy = DBL_MAX;
    if (select_phase3_linear_solve_line_search_candidate(
      ctx, edge_collapser, new_point, line_search_point, line_search_energy))
    {
      new_point = line_search_point;
      energy = line_search_energy;
      return true;
    }
    return select_phase3_feasibility_fallback(eh, ctx, edge_collapser, new_point, energy);
  }

  if (is_phase3_qem_original_strategy() &&
    !phase3_placement_satisfies_hard_constraints(ctx, edge_collapser, new_point))
    return select_phase3_feasibility_fallback(eh, ctx, edge_collapser, new_point, energy);

  return true;
}

bool CollapseStage::compute_phase3_queue_candidate_data(
  EdgeHandle eh, Vec3d& new_point, Phase3QueueComponents& components)
{
  double placement_energy = DBL_MAX;
  if (!compute_phase3_queue_placement_candidate(
    eh, new_point, placement_energy))
    return false;

  auto score_collapser = new_edge_collapser();
  set_phase3_edge_collapser_flags(score_collapser);
  if (!score_collapser.init(eh))
    return false;

  new_point = score_collapser.constrained_target_point(new_point);

  const Phase3PlacementContext score_ctx =
    make_phase3_placement_context(eh, score_collapser);
  components = evaluate_phase3_queue_components(score_ctx, new_point);
  return std::isfinite(components.plane) &&
    std::isfinite(components.triangle_quality) &&
    std::isfinite(components.uniformity);
}

void CollapseStage::refresh_phase3_average_lengths()
{
  avg_edge_length = 0.0;
  for (EdgeHandle eh : rm->edges())
    avg_edge_length += rm->data(eh).edge_length;
  avg_edge_length = rm->n_edges() == 0 ? 0.0 : avg_edge_length / rm->n_edges();

  if (avg_source_edge_length_initialized)
    return;

  avg_source_edge_length = 0.0;
  for (EdgeHandle eh : om->edges())
    avg_source_edge_length += om->data(eh).edge_length;
  avg_source_edge_length = om->n_edges() == 0 ? 0.0 : avg_source_edge_length / om->n_edges();
  avg_source_edge_length_initialized = true;
}

void CollapseStage::initialize_phase3_target_edge_length(size_t target_vertices_num)
{
  if (phase3_target_edge_length_initialized &&
    phase3_target_vertices_num == target_vertices_num)
    return;

  refresh_phase3_average_lengths();
  phase3_target_edge_length = avg_edge_length;

  // With no vertex target, keep the initial cage mean as the queue's fixed
  // length scale instead of estimating a length from a fictitious target.
  if (target_vertices_num == 0)
  {
    phase3_target_edge_length_initialized = true;
    phase3_target_vertices_num = 0;
    Logger::user_logger->info(
      "phase 3 target edge length: no vertex target; fixed initial mean {}.",
      phase3_target_edge_length);
    return;
  }

  const double initial_vertices = static_cast<double>(rm->n_vertices());
  const double initial_edges = static_cast<double>(rm->n_edges());
  const double initial_faces = static_cast<double>(rm->n_faces());
  const double euler_characteristic =
    initial_vertices - initial_edges + initial_faces;
  const double initial_topology_scale =
    initial_vertices - euler_characteristic;
  const double target_topology_scale =
    static_cast<double>(target_vertices_num) - euler_characteristic;

  if (avg_edge_length > 0.0 &&
    std::isfinite(avg_edge_length) &&
    initial_topology_scale > 0.0 &&
    target_topology_scale > 0.0)
  {
    const double length_scale =
      std::sqrt(initial_topology_scale / target_topology_scale);
    if (std::isfinite(length_scale) && length_scale > 0.0)
      phase3_target_edge_length = avg_edge_length * length_scale;
  }
  else
  {
    Logger::user_logger->warn(
      "phase 3 target edge length: invalid topology/length inputs; "
      "falling back to the initial mean edge length.");
  }

  phase3_target_edge_length_initialized = true;
  phase3_target_vertices_num = target_vertices_num;

  Logger::user_logger->info(
    "phase 3 target edge length: initial V/E/F {}/{}/{}, Euler characteristic {}, "
    "initial mean {}, target vertices {}, estimated target mean {}.",
    rm->n_vertices(), rm->n_edges(), rm->n_faces(),
    euler_characteristic, avg_edge_length, target_vertices_num,
    phase3_target_edge_length);
}

bool CollapseStage::enqueue_phase3_candidate(EdgeHandle eh, size_t state)
{
  if (!eh.is_valid() || rm->status(eh).deleted())
    return false;

  Vec3d new_point;
  Phase3QueueComponents components;
  if (!compute_phase3_queue_candidate_data(eh, new_point, components))
    return false;

  const double priority_score =
    evaluate_phase3_queue_score(components);
  if (!std::isfinite(priority_score))
    return false;

  phase3_edges_to_collapse.emplace(eh, state, -priority_score, new_point);
  return true;
}

// Initialize the candidate queue. Each queue item stores an edge and its current
// fast placement estimate; it does not perform any collapse.
void CollapseStage::initialize_phase3_candidates()
{
  update_states.clear();
  update_states.resize(rm->n_edges(), 0);
  phase3_edges_to_collapse = Phase3EdgeRewardQueue();

  refresh_phase3_average_lengths();

  Logger::user_logger->info(
    "phase 3 energy [{}]: initializing {} edge candidates.",
    param->phase3PlacementStrategy, rm->n_edges());

  struct InitialQueueCandidate
  {
    EdgeHandle edge;
    Vec3d point;
    Phase3QueueComponents components;
  };

  std::vector<InitialQueueCandidate> candidates;
  candidates.reserve(rm->n_edges());
  phase3_queue_component_maxima = Phase3QueueComponents();

  const auto include_in_maxima = [&](const Phase3QueueComponents& components)
  {
    phase3_queue_component_maxima.uniformity =
      std::max(
        phase3_queue_component_maxima.uniformity,
        components.uniformity);
  };

  size_t evaluated_edges = 0;
  for (EdgeHandle eh : rm->edges())
  {
    Vec3d new_point;
    Phase3QueueComponents components;
    if (compute_phase3_queue_candidate_data(eh, new_point, components))
    {
      candidates.push_back({ eh, new_point, components });
      include_in_maxima(components);
    }

    evaluated_edges++;
    if (evaluated_edges % 10000 == 0)
    {
      Logger::user_logger->info(
        "phase 3 energy candidates: {}/{} edges, {} enqueued.",
        evaluated_edges, rm->n_edges(), candidates.size());
    }
  }

  for (const InitialQueueCandidate& candidate : candidates)
  {
    const double priority_score =
      evaluate_phase3_queue_score(candidate.components);
    if (std::isfinite(priority_score))
      phase3_edges_to_collapse.emplace(
        candidate.edge, 0, -priority_score, candidate.point);
  }

  Logger::user_logger->info(
    "phase 3 queue pass maximum target-edge ratio: {}.",
    phase3_queue_component_maxima.uniformity);
  Logger::user_logger->info(
    "phase 3 energy candidates: {}/{} edges, {} enqueued.",
    evaluated_edges, rm->n_edges(), phase3_edges_to_collapse.size());
}

// After one collapse, only nearby edges have stale placement scores.
// Recompute and requeue those local candidates.
void CollapseStage::update_phase3_after_collapsing(VertexHandle collapsed_center)
{
  std::vector<EdgeHandle> affected_edges = find_1rv_1re(rm, collapsed_center);
  for (EdgeHandle eh : affected_edges)
  {
    if (!eh.is_valid() || rm->status(eh).deleted())
      continue;
    update_states[eh.idx()]++;
    enqueue_phase3_candidate(eh, update_states[eh.idx()]);
  }
}

// Main Phase 3 replacement:
// 1. Build a priority queue from fast placement estimates.
// 2. Pop one edge, choose its actual placement according to phase3PlacementStrategy.
// 3. Commit the collapse with exact validity checks.
// 4. Requeue only the affected one-ring edges.
size_t CollapseStage::do_phase3_energy_simplification(size_t target_vertices_num)
{
  const bool until_stalled = target_vertices_num == 0;
  if (!until_stalled && rm->n_vertices() <= target_vertices_num)
    return 0;

  if (until_stalled)
    Logger::user_logger->info(
      "phase 3 energy simplification [{}]: {} vertices, collapsing until no further progress.",
      param->phase3PlacementStrategy, rm->n_vertices());
  else
    Logger::user_logger->info(
      "phase 3 energy simplification [{}]: {} -> {} vertices.",
      param->phase3PlacementStrategy, rm->n_vertices(), target_vertices_num);
  initialize_phase3_target_edge_length(target_vertices_num);
  Logger::user_logger->info(
    "phase 3 energy terms: planeWeight {}, triangleQualityWeight {}, curvatureMode [{}], uniformityMode [{}], uniformityWeight {}, robustnessMode [{}].",
    param->planeWeight, param->triangleQualityWeight,
    param->curvatureMode,
    param->uniformityMode, param->uniformityWeight,
    param->robustnessMode);
  if (kEnableTriangleQualityInPhase3LinearSolve && is_phase3_linear_solve_strategy())
  {
    param->validate_triangle_shape_settings();
    Logger::user_logger->info(
      "phase 3 triangle shape directions: triangleShapeTangentWeight {}, triangleShapeHeightWeight {}, triangleShapeNormalWeight {}.",
      param->triangleShapeTangentWeight, param->triangleShapeHeightWeight, param->triangleShapeNormalWeight);
  }
  if (param->uniformityWeight > 0.0 && param->uniformityMode != "none")
  {
    if (kEnableUniformityInPhase3LinearSolve && is_phase3_linear_solve_strategy())
      Logger::user_logger->info(
        "phase 3 uniformity: linear-difference form solved together with plane approximation and the triangle surrogate; the queue keeps its edge-length ratio.");
    else if (!kEnableUniformityInPhase3Solve)
      Logger::user_logger->info(
        "phase 3 uniformity: queue-only experiment enabled; placement objectives exclude uniformity.");
  }
  if (!kEnableTriangleQualityInPhase3LinearSolve &&
    is_phase3_linear_solve_strategy() && param->triangleQualityWeight > 0.0)
  {
    Logger::user_logger->info(
      "phase 3 triangle quality: queue-only experiment enabled for linear_solve; placement objectives exclude triangle quality.");
  }
  if (is_phase3_qem_original_strategy())
    Logger::user_logger->info(
      "phase 3 QEM: pure Garland-Heckbert QEM is enabled; collision checks and non-QEM energies are disabled.");
  if (is_phase3_linear_solve_strategy() &&
    param->phase3LinearSolveCollisionReject)
    Logger::user_logger->info(
      "phase 3 linear solve: phase3LinearSolveCollisionReject is enabled; invalid raw placements are rejected instead of backtracked.");

  if (is_phase3_qem_original_strategy())
    initialize_phase3_qem_quadrics();

  size_t total_collapsed_edges = 0;
  size_t total_attempted_edges = 0;
  size_t remaining_vertices = rm->n_vertices();
  size_t no_progress_passes = 0;
  const size_t max_phase3_passes = 12;

  for (size_t pass = 1;
       (until_stalled || remaining_vertices > target_vertices_num) &&
         (until_stalled || pass <= max_phase3_passes);
       pass++)
  {
    Logger::user_logger->info("phase 3 energy pass {}.", pass);
    initialize_phase3_candidates();

    auto edge_collapser = new_edge_collapser();
    set_phase3_edge_collapser_flags(edge_collapser);

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
    while ((until_stalled || remaining_vertices > target_vertices_num) &&
           !phase3_edges_to_collapse.empty())
    {
      Phase3EdgeReward edge_reward = phase3_edges_to_collapse.top();
      phase3_edges_to_collapse.pop();

      if (edge_reward.state < update_states[edge_reward.eh.idx()])
        continue;
      if (rm->status(edge_reward.eh).deleted())
        continue;

      if (!edge_collapser.init(edge_reward.eh))
        continue;

      Eigen::Matrix4d merged_qem_quadric = Eigen::Matrix4d::Zero();
      bool has_merged_qem_quadric = false;
      if (is_phase3_qem_original_strategy())
      {
        const HalfedgeHandle collapse_heh = rm->halfedge_handle(edge_reward.eh, 0);
        const VertexHandle from_vh = rm->from_vertex_handle(collapse_heh);
        const VertexHandle to_vh = rm->to_vertex_handle(collapse_heh);
        merged_qem_quadric = phase3_qem_quadric(from_vh) + phase3_qem_quadric(to_vh);
        has_merged_qem_quadric = true;
      }

      pass_attempted_edges++;
      total_attempted_edges++;
      Phase3PlacementDecision placement;
      double edge_newton_seconds = 0.0;
      if (!choose_phase3_collapse_placement(
        edge_reward.eh, edge_collapser, edge_reward.initial_point,
        placement, edge_newton_seconds))
      {
        pass_placement_failed_edges++;
        const auto now = std::chrono::steady_clock::now();
        if (pass_attempted_edges % 100 == 0 ||
          std::chrono::duration<double>(now - last_progress_log_time).count() >= 10.0)
        {
          Logger::user_logger->info(
            "phase 3 energy pass {} progress: attempted {}, collapsed {}, newton_attempted {}, newton_refined {}, newton_failed {}, fallback {}, placement_failed {}, collapse_failed {}, queue {}, remain {}, newton_time {:.3f}s, exact_collapse_time {:.3f}s.",
            pass, pass_attempted_edges, pass_collapsed_edges,
            pass_newton_attempted_edges, pass_newton_refined_edges,
            pass_newton_failed_edges, pass_fallback_edges,
            pass_placement_failed_edges, pass_collapse_failed_edges,
            phase3_edges_to_collapse.size(), remaining_vertices,
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
        if (has_merged_qem_quadric)
          update_phase3_qem_quadric_after_collapse(center_v, merged_qem_quadric);
        pass_collapsed_edges++;
        total_collapsed_edges++;
        remaining_vertices--;

        update_phase3_after_collapsing(center_v);

        if (total_collapsed_edges % 1000 == 0)
        {
          Logger::user_logger->info(
            "phase 3 energy collapsed {} edges; {} vertices remain.",
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
          "phase 3 energy pass {} progress: attempted {}, collapsed {}, newton_attempted {}, newton_refined {}, newton_failed {}, fallback {}, placement_failed {}, collapse_failed {}, queue {}, remain {}, newton_time {:.3f}s, exact_collapse_time {:.3f}s.",
          pass, pass_attempted_edges, pass_collapsed_edges,
          pass_newton_attempted_edges, pass_newton_refined_edges,
          pass_newton_failed_edges, pass_fallback_edges,
          pass_placement_failed_edges, pass_collapse_failed_edges,
          phase3_edges_to_collapse.size(), remaining_vertices,
          pass_newton_seconds, pass_exact_collapse_seconds);
        last_progress_log_time = now;
      }
    }

    Logger::user_logger->info(
      "phase 3 energy pass {} collapsed {} edges after attempting {} candidates (newton_attempted {}, newton_refined {}, newton_failed {}, fallback {}, placement_failed {}, collapse_failed {}, newton_time {:.3f}s, exact_collapse_time {:.3f}s).",
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
    "phase 3 energy collapsed {} edges after attempting {} candidates.",
    total_collapsed_edges, total_attempted_edges);

  Logger::user_logger->info("[{}] vertices and [{}] faces remained.", rm->n_vertices(), rm->n_faces());
  return total_collapsed_edges;
}

bool CollapseStage::find_collapse_hausdorff_deviation(
  EdgeHandle eh, double& local_hd_before, double& local_hd_after, Vec3d& new_point)
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
      /*check_face_orientation_violation*/true, /*check_selfinter*/true, /*check_inter*/true);

  local_hd_before = edge_collapser.local_Hausdorff_before_collapsing();

  if (edge_collapser.has_fixed_rail_target() || is_optimization_collapse_placement_method())
  {
    if (edge_collapser.has_fixed_rail_target())
    {
      new_point = edge_collapser.constrained_target_point(rm->calc_edge_midpoint(heh));
      if (!edge_collapser.target_point_is_valid(new_point, nullptr))
        return false;
    }
    else
    {
      double optimized_energy = DBL_MAX;
      if (!refine_collapse_placement_with_newton(eh, edge_collapser, new_point, optimized_energy))
        return false;
    }

    VertexHandle local_center_v;
    auto local_mesh = construct_local_mesh(rm, edge_collapser.get_halfedges(), new_point, local_center_v);
    double hd_threshold = allow_negtive ? max_distance_error : cage_infinite_fp;
    local_hd_after = edge_collapser.local_Hausdorff_after_collapsing(local_mesh.get(), new_point, hd_threshold);
    if (local_hd_after == DBL_MAX)
      return false;

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

#if USE_TREE_SEARCH
  ot->set_hint(candidate_points[0]);
#endif
#pragma omp parallel for schedule(dynamic)
  for (int i = 0;i < (int)candidate_points_size;i++)
  {
    double local_hd_i = edge_collapser.local_Hausdorff_after_collapsing(
      &local_meshes[i], candidate_points[i], minimal_local_hd);
  #pragma omp critical
    if (local_hd_i < minimal_local_hd)
    {
      minimal_local_hd = local_hd_i;
      minimal_idx = i;
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
    return true;
  }
  else return false;
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
    Vec3d new_point;
    if (find_collapse_hausdorff_deviation(eh, local_hd_before, local_hd_after, new_point))
    {
      if (allow_negtive)
      {
        if (local_hd_after < max_distance_error)
        {
          double error = local_hd_after + 0.1 * (rm->data(eh).edge_length - avg_edge_length);
          edges_to_collapse.emplace(eh, 0, -error, new_point);
          enqueued_edges++;
        }
      }
      else if (local_hd_before - local_hd_after >= 0.0)
      {
        double error = local_hd_before - local_hd_after + 0.1 * (avg_edge_length - rm->data(eh).edge_length);
        edges_to_collapse.emplace(eh, 0, error, new_point);
        enqueued_edges++;
      }
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
    Vec3d new_point;
    if (find_collapse_hausdorff_deviation(eh, local_hd_before, local_hd_after, new_point))
    {
      if (allow_negtive)
      {
        if (local_hd_after < max_distance_error)
        {
          double error = local_hd_after + 0.1 * (rm->data(eh).edge_length - avg_edge_length);
          edges_to_collapse.emplace(eh, update_states[eh.idx()], -error, new_point);
        }
      }
      else if (local_hd_before - local_hd_after >= 0.0)
      {
        double error = local_hd_before - local_hd_after + 0.1 * (avg_edge_length - rm->data(eh).edge_length);
        edges_to_collapse.emplace(eh, update_states[eh.idx()], error, new_point);
      }
    }
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
    /*check_face_orientation_violation*/true, /*check_selfinter*/true, /*check_inter*/true);

  initialize_collapse_edges_reward();
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
