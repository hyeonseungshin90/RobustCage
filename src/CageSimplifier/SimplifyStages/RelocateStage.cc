#include "RelocateStage.hh"
#include <algorithm>
#include <cmath>
#include <random>

namespace Cage
{
namespace CageSimp
{
RelocateStage::RelocateStage(
  SMeshT* original, SMeshT* cage, ParamRelocateStage* p,
  VertexTree* vertex_tree, DFaceTree* original_tree, LightDFaceTree* remeshing_tree,
  FaceGrid* original_grid,
  double _original_diagonal_length)
  :om(original), rm(cage), param(p),
  vt(vertex_tree), ot(original_tree), lrt(remeshing_tree), og(original_grid),
  original_diagonal_length(_original_diagonal_length)
{}

void RelocateStage::update(size_t _candidate_points_size, bool _allow_negtive, double _max_distance_error)
{
  candidate_points_size = _candidate_points_size;
  allow_negtive = _allow_negtive;
  max_distance_error = _max_distance_error;
}

std::vector<Vec3d> RelocateStage::generate_candidate_points_for_relocate(VertexHandle vh, VertexRelocater& vertex_relocater)
{
  // calculate radius, which is average length of adjacent edges.
  double radius = 0.0;
  for (EdgeHandle ve : rm->ve_range(vh))
    radius += rm->data(ve).edge_length;
  radius /= rm->valence(vh);
  radius *= 0.2;
  radius = std::min(radius, original_diagonal_length * 0.01);

  // tangential relaxation target
  Vec3d vertex_target;
  vertex_target = vertex_relocater.find_weighted_tangential_smooth_target();

  // local coordinate system on target point
  const Vec3d& vertex_normal = rm->normal(vh);
  Vec3d local_axis_x, local_axis_y;
  make_coordinate_system(vertex_normal, local_axis_x, local_axis_y);

  std::vector<Vec3d> points; points.reserve(candidate_points_size);
  // generate regular points
  points.push_back(vertex_target);

  // generate random points
  while (points.size() < candidate_points_size)
  {
    double height = (double)rand() / (double)RAND_MAX;
    double beta = (double)rand() / (double)RAND_MAX;
    double len = (double)rand() / (double)RAND_MAX;
    height = (height - 0.5) * original_diagonal_length * 0.005;
    beta = beta * 2 * M_PI;

    Vec3d new_point = vertex_target + (local_axis_x * cos(beta) + local_axis_y * sin(beta)) * len * radius + height * vertex_normal;
    points.push_back(new_point);
  }

  return points;
}

bool RelocateStage::find_relocate_hausdorff_deviation(
  VertexRelocater& vertex_relocater, VertexHandle vh, double& local_hd_before, double& local_hd_after, Vec3d& new_point)
{
  if (!vertex_relocater.init(vh))
    return false;

  local_hd_before = vertex_relocater.local_Hausdorff_before_relocating();

  // get candidate points
  std::vector<Vec3d> candidate_points = generate_candidate_points_for_relocate(vh, vertex_relocater);

  // construct local mesh
  VertexHandle local_center_v;
  std::vector<SMeshT> local_meshes = construct_local_meshes(rm, vertex_relocater.get_halfedges(), candidate_points, local_center_v);

  // find optimal point that minimize local hausdorff distance.
  size_t minimal_idx = 0;
  double minimal_local_hd = DBL_MAX;
  double maximal_quality_delta = -DBL_MAX;
  const bool use_quality_hard = is_triangle_quality_hard_priority_mode();
  const double pre_relocate_quality = use_quality_hard ? calc_pre_relocate_quality(vh) : 0.0;

#if USE_TREE_SEARCH
  ot->set_hint(candidate_points[0]);
#endif
#pragma omp parallel for schedule(dynamic)
  for (int i = 0;i < (int)candidate_points_size;i++)
  {
    const double hd_threshold = use_quality_hard ? DBL_MAX : minimal_local_hd;
    double local_hd_i = vertex_relocater.local_Hausdorff_after_relocating(&local_meshes[i], candidate_points[i], hd_threshold);
    double quality_delta = 0.0;
    bool candidate_ok = local_hd_i != DBL_MAX;
    if (candidate_ok && use_quality_hard)
    {
      quality_delta = calc_min_triangle_quality(&local_meshes[i]) - pre_relocate_quality;
      if (quality_delta < 0.0)
        candidate_ok = false;
    }
  #pragma omp critical 
    if (candidate_ok)
    {
      if (use_quality_hard)
      {
        if (quality_delta > maximal_quality_delta ||
          (std::abs(quality_delta - maximal_quality_delta) <= 1e-12 && local_hd_i < minimal_local_hd))
        {
          maximal_quality_delta = quality_delta;
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
    return true;
  }
  else return false;
}

bool RelocateStage::is_triangle_quality_hard_priority_mode() const
{
  return param->priorityMode == "triangle_quality_hard" || param->priorityMode == "triangle-quality-hard";
}

double RelocateStage::calc_pre_relocate_quality(VertexHandle vh) const
{
  std::vector<FaceHandle> faces;
  for (FaceHandle vf : rm->vf_range(vh))
    faces.push_back(vf);
  return calc_min_triangle_quality(rm, faces);
}

double RelocateStage::calc_triangle_quality(SMeshT* mesh, FaceHandle fh) const
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

double RelocateStage::calc_min_triangle_quality(SMeshT* mesh) const
{
  double min_quality = DBL_MAX;
  for (FaceHandle fh : mesh->faces())
    min_quality = std::min(min_quality, calc_triangle_quality(mesh, fh));
  return min_quality == DBL_MAX ? 0.0 : min_quality;
}

double RelocateStage::calc_min_triangle_quality(SMeshT* mesh, const std::vector<FaceHandle>& faces) const
{
  double min_quality = DBL_MAX;
  for (FaceHandle fh : faces)
    min_quality = std::min(min_quality, calc_triangle_quality(mesh, fh));
  return min_quality == DBL_MAX ? 0.0 : min_quality;
}

void RelocateStage::update_after_relocating(VertexHandle relocate_center)
{
  // update in out error
  std::vector<FaceHandle> faces;
  std::vector<EdgeHandle> edges;
  for (FaceHandle vf : rm->vf_range(relocate_center)) faces.push_back(vf);
  for (EdgeHandle ve : rm->ve_range(relocate_center)) edges.push_back(ve);

  double infinite = DBL_MAX;
#ifdef USE_TREE_SEARCH
  calc_face_out_error(rm, om, ot, faces, infinite);
  calc_edge_out_error(rm, om, ot, edges, infinite);
  calc_vertex_out_error(rm, om, ot, relocate_center);
#else
  calc_face_out_error(rm, om, og, faces, infinite);
  calc_edge_out_error(rm, om, og, edges, infinite);
  calc_vertex_out_error(rm, om, og, relocate_center);
#endif
  calc_face_in_error(rm, faces);

  calc_out_surround_error(rm, relocate_center);
  for (VertexHandle vv : rm->vv_range(relocate_center))
    calc_out_surround_error(rm, vv);
}

void RelocateStage::do_relocate()
{
  auto vertex_relocater = new_vertex_relocater();
  vertex_relocater.set_flags(/*update_links*/true, /*update_target_length*/false, /*update_normals*/true, /*check_wrinkle*/false);
  Logger::user_logger->info("relocate priority mode: {}", param->priorityMode);
  size_t relocated_vertex_num = 0;
  for (size_t it = 0;it < param->smoothIter;it++)
  {
    for (VertexHandle vh : rm->vertices())
    {
      Vec3d new_point;
      double local_hd_before, local_hd_after = DBL_MAX;

      if (find_relocate_hausdorff_deviation(vertex_relocater, vh, local_hd_before, local_hd_after, new_point))
      {
        bool do_it = false;
        if (!allow_negtive && local_hd_before - local_hd_after >= 0.0)
        {
          do_it = true;
        }
        else if (allow_negtive&& local_hd_after < max_distance_error)
        {
          double local_hd_diff = (local_hd_after - local_hd_before) / local_hd_before;
          if (local_hd_diff <= 0.0)
            do_it = true;
          else if (local_hd_diff >= 1.0)
            do_it = false;
          else
            do_it = ((double)rand() / RAND_MAX) > local_hd_diff;
        }

        if (do_it)
        {
          vertex_relocater.relocate(new_point);
          relocated_vertex_num += 1;
          update_after_relocating(vh);
        }
      }
    }
  }
  Logger::user_logger->info("relocated {} vertices.", relocated_vertex_num);
}

}// namespace CageSimp
}// namespace Cage
