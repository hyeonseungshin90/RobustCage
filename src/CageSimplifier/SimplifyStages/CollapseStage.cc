#include "CollapseStage.hh"
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

CollapseStage::CollapseStage(
  SMeshT* original, SMeshT* cage, ParamCollapseStage* p,
  DFaceTree* original_tree, LightDFaceTree* remeshing_tree,
  FaceGrid* original_grid,
  double _original_diagonal_length)
  :om(original), rm(cage), param(p),
  ot(original_tree), lrt(remeshing_tree), og(original_grid),
  original_diagonal_length(_original_diagonal_length)
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

  local_hd_before = edge_collapser.local_Hausdorff_before_collapsing();

  // get candidate points
  std::vector<Vec3d> candidate_points = generate_candidate_points_for_collapse(eh, edge_collapser);

  // construct local mesh
  VertexHandle local_center_v;
  std::vector<SMeshT> local_meshes = construct_local_meshes(rm, edge_collapser.get_halfedges(), candidate_points, local_center_v);

  // find optimal point that minimize local hausdorff distance.
  size_t minimal_idx = 0;
  double minimal_local_hd = DBL_MAX;
  double maximal_priority_score = -DBL_MAX;
  const bool use_post_metric = is_post_metric_priority_mode();
  const double pre_collapse_metric = use_post_metric ? calc_pre_collapse_metric(eh) : 0.0;

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

  // initialize for all edges
  for (EdgeHandle eh : rm->edges())
  {
    double local_hd_after, local_hd_before;
    double priority_score;
    Vec3d new_point;
    if (find_collapse_hausdorff_deviation(eh, local_hd_before, local_hd_after, new_point, priority_score))
      try_enqueue_collapse_candidate(eh, 0, local_hd_before, local_hd_after, new_point, priority_score);
  }
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
