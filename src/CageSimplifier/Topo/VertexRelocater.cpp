#include "VertexRelocater.h"
#include "CageSimplifier/Geom/TangentialSmoothing.h"
#include "Geometry/Exact/TriTriIntersect.h"

namespace Cage
{
namespace CageSimp
{

bool VertexRelocater::init(VertexHandle _v)
{
  clear();
  if (!rm || !_v.is_valid() || _v.idx() >= static_cast<int>(rm->n_vertices()) ||
    rm->status(_v).deleted())
    return false;
  relocate_vh = _v;

  find_faces_affected();

  initialized = !halfedges.empty();

  return initialized;
}

void VertexRelocater::clear()
{
  halfedges.clear();
  one_ring_faces.clear();
  faces_in_links.clear();
  initialized = false;
}

void VertexRelocater::set_flags(
  bool _f_update_links, bool _f_update_target_length,
  bool _f_update_normals, bool _f_check_wrinkle)
{
  f_update_links = _f_update_links;
  f_update_target_length = _f_update_target_length;
  f_update_normals = _f_update_normals;
  f_check_wrinkle = _f_check_wrinkle;
}

void VertexRelocater::relocate(const Vec3d& new_point)
{
  if (f_update_links)
    backup_links();

  // real relocating
  rm->set_point(relocate_vh, new_point);
  // All guards evaluated this candidate as a double point. Keeping an old
  // rational point here would make later exact predicates see its old position.
  rm->data(relocate_vh).ep.reset();

  // update after relocating
  if (f_update_target_length)
    update_target_len();
  update_length_and_area();
  if (f_update_normals)
    update_normals();
  update_remeshing_tree();
  if (f_update_links)
    generate_links();
}

bool VertexRelocater::try_relocate_vertex(const Vec3d& new_point)
{
  ASSERT(initialized, "unintialized vertex relocater.");
  if (!TangentialSmoothing::finite_point(new_point))
    return false;
  // check and backup before relocating
  if (relocate_would_cause_degenerate(new_point))
    return false;
  if (relocate_would_cause_intersection(new_point))
    return false;
  if (f_check_wrinkle && relocate_would_cause_wrinkle(new_point))
    return false;

  relocate(new_point);

  return true;
}

/// @brief A simple wrap of try_relocate_vertex.
bool VertexRelocater::try_relocate_vertex(VertexHandle _v, const Vec3d& new_point)
{
  if (init(_v))
    return try_relocate_vertex(new_point);
  else
    return false;
}

double VertexRelocater::local_Hausdorff_before_relocating()const
{
  ASSERT(initialized, "relocater not initialized.");

  return rm->data(relocate_vh).out_surround_error;
}

/// @brief Assume relocater is initialized, calculate local Hausdorff distance.
/// @param [in] local_rm local remeshing mesh, neighbor faces of vertex.
/// @param [in] new_point new point of relocating vertex.
/// @param [in] threshold hausdorff distance threshold.
/// @return false if violate constraints.
double VertexRelocater::local_Hausdorff_after_relocating(SMeshT* local_rm, const Vec3d& new_point, double threshold)const
{
  ASSERT(initialized, "relocater not initialized.");

  if (f_check_wrinkle && relocate_would_cause_wrinkle(new_point))
    return DBL_MAX;

  if (relocate_would_cause_intersection(new_point))
    return DBL_MAX;

  double hd = 0.0;

  // 1. construct local tree, calculate "in" hausdorff distance.
  FaceTree local_rt(*local_rm);
  for (FaceHandle f : one_ring_faces)
  {
    for (Link& link : rm->data(f).face_in_links)
    {
      auto cd = local_rt.closest_distance_and_face_handle(link.first);
      if (cd.first > threshold)
        return DBL_MAX;
      if (cd.first > hd)
        hd = cd.first;
    }
  }
  // 2. calculate "out" hausdorff distance.
#ifdef USE_TREE_SEARCH
  hd = std::max(hd, calc_out_error(local_rm, om, ot, threshold));
#else
  hd = std::max(hd, calc_out_error(local_rm, om, og, threshold));
#endif
  return hd;
}

Vec3d VertexRelocater::find_smooth_target()const
{
  ASSERT(initialized, "vertex relocater not initialized.");

  Vec3d q(0.0, 0.0, 0.0);
  for (VertexHandle vv : rm->vv_range(relocate_vh))
  {
    q += rm->point(vv);
  }
  q /= rm->valence(relocate_vh);
  return q;
}

/// @brief Uniform Laplacian smoothing.
Vec3d VertexRelocater::find_tangential_smooth_target()const
{
  ASSERT(initialized, "vertex relocater not initialized.");

  Vec3d q(0.0, 0.0, 0.0);
  for (VertexHandle vv : rm->vv_range(relocate_vh))
  {
    q += rm->point(vv);
  }
  q /= rm->valence(relocate_vh);
  return q + (rm->normal(relocate_vh) | (rm->point(relocate_vh) - q)) * rm->normal(relocate_vh);
}

Vec3d VertexRelocater::find_weighted_smooth_target()const
{
  ASSERT(initialized, "vertex relocater not initialized.");

  const Vec3d& vertex_normal = rm->normal(relocate_vh);
  Vec3d q(0., 0., 0.);
  double weight = 0.0;

  for (FaceHandle f : rm->vf_range(relocate_vh))
  {
    Vec3d barycenter(0., 0., 0.);
    for (VertexHandle fv : rm->fv_range(f))
      barycenter += rm->point(fv);
    barycenter /= 3.0;
    q += barycenter * rm->data(f).face_area;
    weight += rm->data(f).face_area;
  }
  if (weight != 0.0)
  {
    q /= weight;
    return q;
  }
  else return find_smooth_target();
}

Vec3d VertexRelocater::find_weighted_tangential_smooth_target()const
{
  ASSERT(initialized, "vertex relocater not initialized.");

  const Vec3d& vertex_normal = rm->normal(relocate_vh);
  Vec3d q(0., 0., 0.);
  double weight = 0.0;

  for (FaceHandle f : rm->vf_range(relocate_vh))
  {
    Vec3d barycenter(0., 0., 0.);
    for (VertexHandle fv : rm->fv_range(f))
      barycenter += rm->point(fv);
    barycenter /= 3.0;
    q += barycenter * rm->data(f).face_area;
    weight += rm->data(f).face_area;
  }
  if (weight != 0.0)
  {
    q /= weight;
    return q + (vertex_normal | (rm->point(relocate_vh) - q)) * vertex_normal;
  }
  else return find_tangential_smooth_target();
}

bool VertexRelocater::find_area_equalizing_tangential_smooth_target(Vec3d& target)const
{
  ASSERT(initialized, "vertex relocater not initialized.");
  using TangentialSmoothing::finite_point;
  using TangentialSmoothing::mixed_voronoi_area_at_vertex;

  const Vec3d& center = rm->point(relocate_vh);
  target = center;
  if (!finite_point(center))
    return false;

  // Sum raw face cross products, exactly as in the collapse predictor. This
  // uses current geometry rather than a potentially stale cached normal.
  Vec3d normal_sum(0.0, 0.0, 0.0);
  for (HalfedgeHandle h : halfedges)
  {
    const Vec3d face_normal =
      (rm->point(rm->from_vertex_handle(h)) - center).cross(
        rm->point(rm->to_vertex_handle(h)) - center);
    if (!finite_point(face_normal))
      return false;
    normal_sum += face_normal;
  }
  const double normal_length = normal_sum.length();
  if (!std::isfinite(normal_length) || normal_length <= 0.0)
    return false;
  const Vec3d normal = normal_sum / normal_length;

  Vec3d gravity_centroid(0.0, 0.0, 0.0);
  double total_area = 0.0;
  for (VertexHandle neighbor : rm->vv_range(relocate_vh))
  {
    if (!finite_point(rm->point(neighbor)))
      return false;
    double neighbor_area = 0.0;
    // Each weight is the neighbor's FULL vertex area, including triangles
    // outside the center fan. Face-centroid weighting is a different method.
    for (FaceHandle fh : rm->vf_range(neighbor))
    {
      if (!fh.is_valid() || rm->status(fh).deleted())
        continue;
      Vec3d other_points[2];
      size_t other_count = 0;
      for (VertexHandle fv : rm->fv_range(fh))
      {
        if (fv != neighbor && other_count < 2)
          other_points[other_count++] = rm->point(fv);
      }
      if (other_count == 2)
        neighbor_area += mixed_voronoi_area_at_vertex(
          rm->point(neighbor), other_points[0], other_points[1]);
    }
    if (!std::isfinite(neighbor_area))
      return false;
    if (neighbor_area > 0.0)
    {
      gravity_centroid += neighbor_area * rm->point(neighbor);
      total_area += neighbor_area;
    }
  }
  if (!std::isfinite(total_area) || total_area <= 0.0)
    return false;
  gravity_centroid /= total_area;
  const Vec3d update = gravity_centroid - center;
  const Vec3d tangential_update = update - (normal | update) * normal;
  target = center + tangential_update;
  return finite_point(target);
}

/// @brief tangential relaxation in "Isotropic Surface Remeshing without Large
/// and Small Angles".
Vec3d VertexRelocater::find_smooth_target_with_target_length()const
{
  ASSERT(initialized, "vertex relocater not initialized.");

  Vec3d q(0.0, 0.0, 0.0);
  double w = 0.0;
  for (FaceHandle fh : rm->vf_range(relocate_vh))
  {
    double t_i = rm->data(fh).face_area;
    double L_b_i = 0.0;
    Vec3d b_i(0.0, 0.0, 0.0);
    for (VertexHandle vv : rm->fv_range(fh))
    {
      b_i += rm->point(vv);
      L_b_i += rm->data(vv).target_length;
    }
    b_i /= 3.0;
    L_b_i /= 3.0;
    w += (t_i * L_b_i);
    q += (t_i * L_b_i) * b_i;
  }
  if (w != 0.0)
  {
    q /= w;
    return q + (rm->normal(relocate_vh) | (rm->point(relocate_vh) - q)) * rm->normal(relocate_vh);
  }
  else return find_tangential_smooth_target();
}

void VertexRelocater::update_normals()
{
  for (FaceHandle vf : rm->vf_range(relocate_vh))
    rm->update_normal(vf);
  rm->update_normal(relocate_vh);
  for (VertexHandle vv : rm->vv_range(relocate_vh))
    rm->update_normal(vv);
}

void VertexRelocater::find_faces_affected()
{
  halfedges.reserve(rm->valence(relocate_vh));
  for (HalfedgeHandle voh : rm->voh_range(relocate_vh))
  {
    // The outer halfedge of an open mesh has no incident triangle.
    if (rm->is_boundary(voh))
      continue;
    halfedges.push_back(rm->next_halfedge_handle(voh));
    one_ring_faces.insert(rm->face_handle(voh));
  }
}

void VertexRelocater::backup_links()
{
  faces_in_links = backup_local_in_links(
    rm, std::vector<FaceHandle>(one_ring_faces.begin(), one_ring_faces.end()));
}

void VertexRelocater::generate_links()
{
  FaceTree local_face_tree(*rm, std::vector<FaceHandle>(one_ring_faces.begin(), one_ring_faces.end()));

  for (Link& link : faces_in_links)
  {
    auto cp = local_face_tree.closest_point_and_face_handle(link.first);
    link.second = cp.first;
    rm->data(cp.second).face_in_links.push_back(link);
  }
}

bool VertexRelocater::relocate_would_cause_wrinkle(const Vec3d& new_point)const
{
  return check_wrinkle(rm, halfedges, new_point);
}

bool VertexRelocater::relocate_would_cause_intersection(const Vec3d& new_point)const
{
  if (check_intersection(rm, halfedges, new_point, nullptr, ot, og, lrt, one_ring_faces))
    return true;

  // The shared checker excludes all old fan faces from its BVH queries.
  // Compare the proposed fan triangles too, allowing only their topological
  // common edge/vertex through the same exact predicates used by the checker.
  for (size_t i = 0; i < halfedges.size(); ++i)
  {
    const VertexHandle a = rm->from_vertex_handle(halfedges[i]);
    const VertexHandle b = rm->to_vertex_handle(halfedges[i]);
    for (size_t j = i + 1; j < halfedges.size(); ++j)
    {
      const VertexHandle c = rm->from_vertex_handle(halfedges[j]);
      const VertexHandle d = rm->to_vertex_handle(halfedges[j]);
      VertexHandle shared, other_i, other_j;
      if (a == c || a == d)
      {
        shared = a;
        other_i = b;
        other_j = a == c ? d : c;
      }
      else if (b == c || b == d)
      {
        shared = b;
        other_i = a;
        other_j = b == c ? d : c;
      }
      if (shared.is_valid())
      {
        if (Geometry::triangle_do_overlap(
          rm->point(other_i), new_point, rm->point(shared), rm->point(other_j),
          rm->data(other_i).ep.get(), nullptr, rm->data(shared).ep.get(),
          rm->data(other_j).ep.get()))
          return true;
      }
      else if (Geometry::triangle_do_intersect(
        rm->point(a), rm->point(b), new_point, rm->point(c), rm->point(d),
        rm->data(a).ep.get(), rm->data(b).ep.get(), nullptr,
        rm->data(c).ep.get(), rm->data(d).ep.get()))
        return true;
    }
  }
  return false;
}

bool VertexRelocater::relocate_would_cause_degenerate(const Vec3d& new_point)const
{
  return check_degenerate(rm, halfedges, new_point, nullptr);
}

void VertexRelocater::update_target_len()
{
  // update target length
  auto& new_point = rm->point(relocate_vh);
  VertexHandle closest_original_vertex = VertexHandle(vt->closest_vertex(new_point).idx());
  double new_target_length = om->data(closest_original_vertex).target_length;
  rm->data(relocate_vh).target_length = new_target_length;
  for (EdgeHandle ve : rm->ve_range(relocate_vh))
  {
    rm->data(ve).target_length = std::min(new_target_length, rm->data(ve).target_length);
  }
}

void VertexRelocater::update_length_and_area()
{
  // update edge length
  for (EdgeHandle ve : rm->ve_range(relocate_vh))
    rm->data(ve).edge_length = rm->calc_edge_length(ve);
  // calculate face area after updating all edges' length.
  for (FaceHandle f : one_ring_faces)
    rm->data(f).face_area = calc_face_area(rm, f);
}

void VertexRelocater::update_remeshing_tree()
{
  for (FaceHandle f : one_ring_faces)
    lrt->update(rm, f);
}

}// namespace CageSimp
}// namespace Cage
