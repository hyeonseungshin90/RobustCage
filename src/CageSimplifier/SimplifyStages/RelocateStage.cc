#include "RelocateStage.hh"
#include "CageSimplifier/Geom/TangentialSmoothing.h"
#include <algorithm>
#include <cmath>
#include <random>

namespace Cage
{
namespace CageSimp
{
namespace
{
double max_component_magnitude(const Vec3d& vector)
{
  return std::max({ std::abs(vector.x()), std::abs(vector.y()), std::abs(vector.z()) });
}

double min_fan_quality(
  SMeshT* mesh, const std::vector<HalfedgeHandle>& halfedges, const Vec3d& center)
{
  if (halfedges.empty() || !TangentialSmoothing::finite_point(center))
    return 0.0;
  double min_quality = DBL_MAX;
  for (HalfedgeHandle h : halfedges)
  {
    Vec3d a = mesh->point(mesh->from_vertex_handle(h)) - center;
    Vec3d b = mesh->point(mesh->to_vertex_handle(h)) - center;
    if (!TangentialSmoothing::finite_point(a) || !TangentialSmoothing::finite_point(b))
      return 0.0;
    // Rescaling leaves quality unchanged and avoids overflow/underflow in
    // squared lengths and cross products on differently scaled inputs.
    const double scale = std::max({ std::abs(a.x()), std::abs(a.y()), std::abs(a.z()),
      std::abs(b.x()), std::abs(b.y()), std::abs(b.z()) });
    if (!std::isfinite(scale) || scale <= 0.0)
      return 0.0;
    a /= scale;
    b /= scale;
    const double denominator = a.sqrnorm() + b.sqrnorm() + (a - b).sqrnorm();
    const double quality = 2.0 * std::sqrt(3.0) * a.cross(b).length() / denominator;
    if (!std::isfinite(quality) || quality <= 0.0)
      return 0.0;
    min_quality = std::min(min_quality, quality);
  }
  return min_quality;
}
}

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

size_t RelocateStage::do_quality_relocate(size_t line_search_max_iter)
{
  param->validate_quality_settings();
  if (line_search_max_iter == 0 ||
    (param->tangentialWeight == 0.0 && param->surfaceWeight == 0.0))
    return 0;
  // Both trees are also required by the candidate intersection checks.
  if (!om || !rm || !ot || ot->empty() || !lrt)
    return 0;

  // A common weight scale does not change the minimizer. Normalize before
  // adding the weights so even large finite settings cannot overflow.
  const double weight_scale = std::max(param->tangentialWeight, param->surfaceWeight);
  const double scaled_tangent_weight = param->tangentialWeight / weight_scale;
  const double scaled_surface_weight = param->surfaceWeight / weight_scale;
  const double weight_sum = scaled_tangent_weight + scaled_surface_weight;
  const double tangent_weight = scaled_tangent_weight / weight_sum;
  const double surface_weight = scaled_surface_weight / weight_sum;

  auto relocater = new_vertex_relocater();
  relocater.set_flags(/*update_links*/false, /*update_target_length*/false,
    /*update_normals*/true, /*check_wrinkle*/true);

  size_t relocated = 0;
  size_t backtracked = 0;
  size_t fixed = 0;
  size_t quality_rejected = 0;
  size_t geometry_rejected = 0;
  for (VertexHandle vh : rm->vertices())
  {
    if (rm->status(vh).deleted())
      continue;
    // Rail support is generally a nonconvex union of half-strips. Keeping
    // rail vertices fixed preserves that constraint throughout smoothing.
    if (rm->data(vh).boundary_rail_id >= 0 || rm->is_boundary(vh))
    {
      ++fixed;
      continue;
    }
    if (!relocater.init(vh))
      continue;

    const Vec3d old_point = rm->point(vh);
    if (!TangentialSmoothing::finite_point(old_point))
      continue;

    Vec3d tangent_offset(0.0, 0.0, 0.0);
    Vec3d surface_offset(0.0, 0.0, 0.0);
    if (tangent_weight > 0.0)
    {
      Vec3d tangent_target;
      if (!relocater.find_area_equalizing_tangential_smooth_target(tangent_target))
        continue;
      tangent_offset = tangent_target - old_point;
    }
    if (surface_weight > 0.0)
    {
      // Query once per proposal; the source is static and its search
      // structures already exist. Keep this correspondence fixed throughout
      // backtracking. This is a point-distance surrogate, not Hausdorff.
      const Vec3d closest = ot->closest_point(old_point).first;
      surface_offset = closest - old_point;
    }
    if (!TangentialSmoothing::finite_point(tangent_offset) ||
      !TangentialSmoothing::finite_point(surface_offset))
      continue;

    // Evaluate the frozen quadratic in relative, scaled coordinates to
    // avoid cancellation from translations and squared-distance overflow.
    const double coordinate_scale = std::max(
      max_component_magnitude(tangent_offset), max_component_magnitude(surface_offset));
    if (coordinate_scale == 0.0)
      continue;
    tangent_offset /= coordinate_scale;
    surface_offset /= coordinate_scale;
    const Vec3d scaled_direction =
      tangent_weight * tangent_offset + surface_weight * surface_offset;
    const Vec3d direction = coordinate_scale * scaled_direction;
    if (!TangentialSmoothing::finite_point(direction) || max_component_magnitude(direction) == 0.0)
      continue;

    const auto energy = [&](const Vec3d& displacement)
    {
      return tangent_weight * (displacement - tangent_offset).sqrnorm() +
        surface_weight * (displacement - surface_offset).sqrnorm();
    };
    const double old_energy = energy(Vec3d(0.0, 0.0, 0.0));
    const double energy_tolerance = 1e-12 * old_energy;
    const double old_quality = min_fan_quality(rm, relocater.get_halfedges(), old_point);
    const double quality_floor = std::min(old_quality, param->minTriangleQuality);
    double alpha = 1.0;
    for (size_t step = 0; step < line_search_max_iter; ++step, alpha *= 0.5)
    {
      const Vec3d candidate = old_point + alpha * direction;
      if (alpha == 0.0 || candidate == old_point)
        break;
      if (!TangentialSmoothing::finite_point(candidate))
        continue;
      const double candidate_energy = energy((candidate - old_point) / coordinate_scale);
      if (!std::isfinite(candidate_energy) || old_energy - candidate_energy <= energy_tolerance)
        continue;
      // Permit distance-driven improvement without requiring triangle
      // quality to increase. Fans already below the floor cannot get worse.
      const double candidate_quality = min_fan_quality(rm, relocater.get_halfedges(), candidate);
      if (candidate_quality <= 0.0 || candidate_quality + 1e-12 < quality_floor)
      {
        ++quality_rejected;
        continue;
      }
      // This checks exact degeneracy, old/new orientation, source-mesh
      // intersection, and self-intersection at the candidate endpoint.
      if (relocater.try_relocate_vertex(candidate))
      {
        ++relocated;
        if (step > 0)
          ++backtracked;
        break;
      }
      ++geometry_rejected;
    }
  }
  Logger::user_logger->info(
    "energy relocation: moved {} vertices ({} backtracked), fixed {} boundary/rail vertices; "
    "rejected candidates: quality {}, geometry {}",
    relocated, backtracked, fixed, quality_rejected, geometry_rejected);
  return relocated;
}

}// namespace CageSimp
}// namespace Cage
