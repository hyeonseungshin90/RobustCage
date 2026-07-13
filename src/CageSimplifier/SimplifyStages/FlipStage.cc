#include "FlipStage.hh"
#include <algorithm>
#include <cmath>

namespace Cage
{
namespace CageSimp
{

FlipStage::FlipStage(
  SMeshT* original, SMeshT* cage, ParamFlipStage* p,
  DFaceTree* original_tree, LightDFaceTree* remeshing_tree,
  FaceGrid* original_grid,
  double _original_diagonal_length)
  :om(original), rm(cage), param(p),
  ot(original_tree), lrt(remeshing_tree), og(original_grid),
  original_diagonal_length(_original_diagonal_length)
{ }

void FlipStage::update(bool _allow_negtive, double _max_distance_error)
{
  allow_negtive = _allow_negtive;
  max_distance_error = _max_distance_error;
}

void FlipStage::initialize_flip_edges_reward()
{
  // clear
  update_states.resize(rm->n_edges());
  std::fill(update_states.begin(), update_states.end(), 0);
  edges_to_flip = FlipEdgeRewardQueue();

  auto edge_flipper = new_edge_flipper();
  // initialize for all edges
  for (EdgeHandle eh : rm->edges())
  {
    if (!edge_flipper.init(eh))
      continue;
    if (edge_flipper.flip_will_cause_over_valence(param->maxValence))
      continue;
    if (!is_triangle_quality_priority_mode() && !edge_flipper.flip_will_decrease_valence())
      continue;
    double local_hd_before = edge_flipper.local_Hausdorff_before_flipping();
    double local_hd_after = edge_flipper.local_Hausdorff_after_flipping();
    try_enqueue_flip_candidate(eh, 0, local_hd_before, local_hd_after);
  }
}

bool FlipStage::try_enqueue_flip_candidate(
  EdgeHandle eh, size_t state, double local_hd_before, double local_hd_after)
{
  if (local_hd_after == DBL_MAX)
    return false;

  if (allow_negtive)
  {
    if (local_hd_after >= max_distance_error)
      return false;
  }
  else if (local_hd_before - local_hd_after < 0.0)
    return false;

  if (is_triangle_quality_priority_mode())
  {
    const double quality_delta = calc_flip_quality_delta(eh);
    if (quality_delta <= 1e-12)
      return false;

    edges_to_flip.emplace(eh, state, quality_delta);
    return true;
  }

  if (allow_negtive)
    // Decrease valence and won't cause out of distance error.
    edges_to_flip.emplace(eh, state, -local_hd_after);
  else
    edges_to_flip.emplace(eh, state, local_hd_before - local_hd_after);
  return true;
}

bool FlipStage::is_triangle_quality_priority_mode() const
{
  return param->priorityMode == "triangle_quality" || param->priorityMode == "triangle-quality" ||
    is_triangle_quality_hard_priority_mode();
}

bool FlipStage::is_triangle_quality_hard_priority_mode() const
{
  return param->priorityMode == "triangle_quality_hard" || param->priorityMode == "triangle-quality-hard";
}

bool FlipStage::is_flip_quality_allowed(EdgeHandle eh) const
{
  return calc_flip_quality_delta(eh) > 1e-12;
}

double FlipStage::calc_flip_quality_delta(EdgeHandle eh) const
{
  return calc_post_flip_quality(eh) - calc_pre_flip_quality(eh);
}

double FlipStage::calc_pre_flip_quality(EdgeHandle eh) const
{
  if (!eh.is_valid() || rm->status(eh).deleted() || !rm->is_flip_ok(eh))
    return 0.0;

  HalfedgeHandle heh = rm->halfedge_handle(eh, 0);
  HalfedgeHandle heh_opp = rm->halfedge_handle(eh, 1);
  FaceHandle fh = rm->face_handle(heh);
  FaceHandle fh_opp = rm->face_handle(heh_opp);
  if (!fh.is_valid() || !fh_opp.is_valid())
    return 0.0;

  return std::min(calc_triangle_quality(fh), calc_triangle_quality(fh_opp));
}

double FlipStage::calc_post_flip_quality(EdgeHandle eh) const
{
  if (!eh.is_valid() || rm->status(eh).deleted() || !rm->is_flip_ok(eh))
    return 0.0;

  HalfedgeHandle a0 = rm->halfedge_handle(eh, 0);
  HalfedgeHandle b0 = rm->halfedge_handle(eh, 1);
  HalfedgeHandle a1 = rm->next_halfedge_handle(a0);
  HalfedgeHandle b1 = rm->next_halfedge_handle(b0);

  VertexHandle va0 = rm->to_vertex_handle(a0);
  VertexHandle va1 = rm->to_vertex_handle(a1);
  VertexHandle vb0 = rm->to_vertex_handle(b0);
  VertexHandle vb1 = rm->to_vertex_handle(b1);

  const Vec3d& pa0 = rm->point(va0);
  const Vec3d& pa1 = rm->point(va1);
  const Vec3d& pb0 = rm->point(vb0);
  const Vec3d& pb1 = rm->point(vb1);

  const double quality_a = calc_triangle_quality(pa1, pb0, pb1);
  const double quality_b = calc_triangle_quality(pa0, pa1, pb1);
  return std::min(quality_a, quality_b);
}

double FlipStage::calc_triangle_quality(FaceHandle fh) const
{
  Vec3d pts[3];
  size_t vertex_count = 0;
  for (VertexHandle vh : rm->fv_range(fh))
  {
    if (vertex_count >= 3)
      return 0.0;
    pts[vertex_count++] = rm->point(vh);
  }
  if (vertex_count != 3)
    return 0.0;

  return calc_triangle_quality(pts[0], pts[1], pts[2]);
}

double FlipStage::calc_triangle_quality(const Vec3d& p0, const Vec3d& p1, const Vec3d& p2) const
{
  const double a = (p1 - p0).length();
  const double b = (p2 - p1).length();
  const double c = (p0 - p2).length();
  const double denom = a * a + b * b + c * c;
  if (denom <= 0.0)
    return 0.0;

  const double area = 0.5 * (p1 - p0).cross(p2 - p0).length();
  return 4.0 * std::sqrt(3.0) * area / denom;
}

void FlipStage::update_after_flipping(EdgeHandle flipped_edge)
{
  // two faces, five edges and four vertices need to be updated.
  HalfedgeHandle flip_he = rm->halfedge_handle(flipped_edge, 0);
  HalfedgeHandle flip_he_opp = rm->halfedge_handle(flipped_edge, 1);
  // update in out error
  std::vector<FaceHandle> faces = { rm->face_handle(flip_he), rm->face_handle(flip_he_opp) };
  std::vector<EdgeHandle> edges = { flipped_edge };
  std::vector<VertexHandle> vertices = {
    rm->to_vertex_handle(flip_he), rm->to_vertex_handle(flip_he_opp) ,
    rm->opposite_vh(flip_he), rm->opposite_vh(flip_he_opp) };

#ifdef USE_TREE_SEARCH
  calc_face_out_error(rm, om, ot, faces, cage_infinite_fp);
  calc_edge_out_error(rm, om, ot, edges, cage_infinite_fp);
#else
  calc_face_out_error(rm, om, og, faces, cage_infinite_fp);
  calc_edge_out_error(rm, om, og, edges, cage_infinite_fp);
#endif
  calc_face_in_error(rm, faces);
  for (VertexHandle v : vertices)
    calc_out_surround_error(rm, v);

  // update affected edges 
  std::vector<EdgeHandle> affected_edges = {
    rm->edge_handle(rm->next_halfedge_handle(flip_he)),
    rm->edge_handle(rm->prev_halfedge_handle(flip_he)),
    rm->edge_handle(rm->next_halfedge_handle(flip_he_opp)),
    rm->edge_handle(rm->prev_halfedge_handle(flip_he_opp)),
  };

  auto edge_flipper = new_edge_flipper();
  for (EdgeHandle eh : affected_edges)
  {
    update_states[eh.idx()]++;

    if (!edge_flipper.init(eh))
      continue;
    if (edge_flipper.flip_will_cause_over_valence(param->maxValence))
      continue;
    if (!is_triangle_quality_priority_mode() && !edge_flipper.flip_will_decrease_valence())
      continue;
    double local_hd_before = edge_flipper.local_Hausdorff_before_flipping();
    double local_hd_after = edge_flipper.local_Hausdorff_after_flipping();
    try_enqueue_flip_candidate(eh, update_states[eh.idx()], local_hd_before, local_hd_after);
  }
}

void FlipStage::do_flip()
{
  auto edge_flipper = new_edge_flipper();
  edge_flipper.set_flags(/*update_links*/true, /*update_target_length*/false, /*update_normals*/true);

  initialize_flip_edges_reward();
  Logger::user_logger->info("flip priority mode: {}", param->priorityMode);
  size_t flipped_edge_num = 0;
  while (!edges_to_flip.empty())
  {
    auto edge_reward = edges_to_flip.top();
    edges_to_flip.pop();

    // out of date
    if (edge_reward.state < update_states[edge_reward.eh.idx()])
      continue;

    if (is_triangle_quality_priority_mode() && !is_flip_quality_allowed(edge_reward.eh))
      continue;

    if (edge_flipper.try_flip_edge(edge_reward.eh))
    {
      update_after_flipping(edge_reward.eh);
      flipped_edge_num++;
    }
  }
  Logger::user_logger->info("flipped {} edges.", flipped_edge_num);
}
}// namespace CageSimp
}// namespace Cage
