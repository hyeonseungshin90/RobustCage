#pragma once

#include "Geometry/Basic/Types.h"
#include <algorithm>
#include <cmath>

namespace Cage
{
namespace CageSimp
{
namespace TangentialSmoothing
{
using Geometry::Vec3d;

inline bool finite_point(const Vec3d& point)
{
  return std::isfinite(point.x()) &&
    std::isfinite(point.y()) &&
    std::isfinite(point.z());
}

// Mixed Voronoi area associated with p in triangle (p, q, r), using the
// midpoint construction for obtuse triangles described by Meyer et al.
// Shared by collapse prediction and Botsch-Kobbelt tangential relocation.
inline double mixed_voronoi_area_at_vertex(
  const Vec3d& p, const Vec3d& q, const Vec3d& r)
{
  if (!finite_point(p) || !finite_point(q) || !finite_point(r))
    return 0.0;

  const Vec3d pq = q - p;
  const Vec3d pr = r - p;
  const Vec3d qr = r - q;
  const double twice_area = pq.cross(pr).length();
  const double max_edge_sqr = std::max(
    pq.sqrnorm(), std::max(pr.sqrnorm(), qr.sqrnorm()));
  if (!std::isfinite(twice_area) || !std::isfinite(max_edge_sqr) ||
    max_edge_sqr <= 0.0 || twice_area <= 1e-14 * max_edge_sqr)
    return 0.0;

  const double triangle_area = 0.5 * twice_area;
  const double dot_p = pq | pr;
  const double dot_q = (p - q) | (r - q);
  const double dot_r = (p - r) | (q - r);

  if (dot_p < 0.0)
    return 0.5 * triangle_area;
  if (dot_q < 0.0 || dot_r < 0.0)
    return 0.25 * triangle_area;

  const double cot_q = dot_q / twice_area;
  const double cot_r = dot_r / twice_area;
  const double area =
    (pr.sqrnorm() * cot_q + pq.sqrnorm() * cot_r) / 8.0;
  return std::isfinite(area) && area > 0.0 ? area : 0.0;
}
} // namespace TangentialSmoothing
} // namespace CageSimp
} // namespace Cage
