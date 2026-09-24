#include "BoundaryRailBuilder.hh"
#include "BoundaryRailNumerics.hh"

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <chrono>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <limits>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "Utils/logger.hh"

namespace Cage
{
namespace CageInit
{
namespace
{
using namespace SurfaceMesh;
using namespace SimpleUtils;
using Geometry::Vec3d;

constexpr int kNoRail = -1;

bool finite_vec(const Vec3d& p)
{
  return std::isfinite(p.x()) && std::isfinite(p.y()) && std::isfinite(p.z());
}

double clamp01(double x)
{
  return std::max(0.0, std::min(1.0, x));
}

struct Box
{
  Vec3d lo = Vec3d(DBL_MAX, DBL_MAX, DBL_MAX);
  Vec3d hi = Vec3d(-DBL_MAX, -DBL_MAX, -DBL_MAX);

  void add(const Vec3d& p)
  {
    lo.minimize(p);
    hi.maximize(p);
  }

  void add(const Box& box)
  {
    add(box.lo);
    add(box.hi);
  }
};

struct RayTriangle
{
  Vec3d p[3];
  VertexHandle vertices[3];
  FaceHandle face;
  Box box;
  Vec3d centroid;
};

struct RayHit
{
  bool hit = false;
  double t = DBL_MAX;
  double weights[3] = {};
  Vec3d point;
  size_t triangle = 0;
};

class RayTriangleBvh
{
public:
  explicit RayTriangleBvh(SMeshT& mesh)
  {
    triangles.reserve(mesh.n_faces());
    for (FaceHandle fh : mesh.faces())
    {
      RayTriangle tri;
      tri.face = fh;
      size_t i = 0;
      for (VertexHandle vh : mesh.fv_range(fh))
      {
        if (i >= 3)
          break;
        tri.vertices[i] = vh;
        tri.p[i] = mesh.point(vh);
        tri.box.add(tri.p[i]);
        i++;
      }
      if (i != 3)
        continue;
      tri.centroid = (tri.p[0] + tri.p[1] + tri.p[2]) / 3.0;
      triangles.push_back(tri);
    }

    order.resize(triangles.size());
    std::iota(order.begin(), order.end(), size_t(0));
    if (!order.empty())
      build_node(0, order.size());
  }

  const RayTriangle& triangle(size_t index) const { return triangles[index]; }

  bool first_intersection(
    const Vec3d& origin, const Vec3d& direction,
    double min_t, RayHit& result) const
  {
    result = RayHit();
    if (nodes.empty())
      return false;

    std::vector<int> stack(1, 0);
    while (!stack.empty())
    {
      const int node_index = stack.back();
      stack.pop_back();
      const Node& node = nodes[node_index];
      double box_t = 0.0;
      if (!ray_box(origin, direction, node.box, min_t, result.t, box_t))
        continue;

      if (node.left < 0)
      {
        for (size_t i = node.begin; i < node.end; i++)
        {
          const size_t tri_index = order[i];
          double t = 0.0;
          double weights[3] = {};
          if (!ray_triangle(origin, direction, triangles[tri_index], min_t, t, weights))
            continue;
          if (t < result.t)
          {
            result.hit = true;
            result.t = t;
            result.weights[0] = weights[0];
            result.weights[1] = weights[1];
            result.weights[2] = weights[2];
            result.point = origin + direction * t;
            result.triangle = tri_index;
          }
        }
      }
      else
      {
        double left_t = 0.0;
        double right_t = 0.0;
        const bool hit_left = ray_box(
          origin, direction, nodes[node.left].box, min_t, result.t, left_t);
        const bool hit_right = ray_box(
          origin, direction, nodes[node.right].box, min_t, result.t, right_t);
        // Push the farther node first so the nearer node is processed first.
        if (hit_left && hit_right)
        {
          if (left_t < right_t)
          {
            stack.push_back(node.right);
            stack.push_back(node.left);
          }
          else
          {
            stack.push_back(node.left);
            stack.push_back(node.right);
          }
        }
        else if (hit_left)
          stack.push_back(node.left);
        else if (hit_right)
          stack.push_back(node.right);
      }
    }
    return result.hit;
  }

private:
  struct Node
  {
    Box box;
    size_t begin = 0;
    size_t end = 0;
    int left = -1;
    int right = -1;
  };

  std::vector<RayTriangle> triangles;
  std::vector<size_t> order;
  std::vector<Node> nodes;

  int build_node(size_t begin, size_t end)
  {
    Node node;
    node.begin = begin;
    node.end = end;
    Box centroid_box;
    for (size_t i = begin; i < end; i++)
    {
      node.box.add(triangles[order[i]].box);
      centroid_box.add(triangles[order[i]].centroid);
    }

    const int index = static_cast<int>(nodes.size());
    nodes.push_back(node);
    if (end - begin <= 8)
      return index;

    const Vec3d extent = centroid_box.hi - centroid_box.lo;
    size_t axis = 0;
    if (extent.y() > extent.x())
      axis = 1;
    if (extent.z() > extent[axis])
      axis = 2;
    const size_t middle = begin + (end - begin) / 2;
    std::nth_element(
      order.begin() + begin, order.begin() + middle, order.begin() + end,
      [&](size_t lhs, size_t rhs)
      {
        return triangles[lhs].centroid[axis] < triangles[rhs].centroid[axis];
      });

    const int left = build_node(begin, middle);
    const int right = build_node(middle, end);
    nodes[index].left = left;
    nodes[index].right = right;
    return index;
  }

  static bool ray_box(
    const Vec3d& origin, const Vec3d& direction, const Box& box,
    double min_t, double max_t, double& near_t)
  {
    double lo = min_t;
    double hi = max_t;
    for (size_t axis = 0; axis < 3; axis++)
    {
      if (std::abs(direction[axis]) <= 1e-30)
      {
        if (origin[axis] < box.lo[axis] || origin[axis] > box.hi[axis])
          return false;
        continue;
      }
      double t0 = (box.lo[axis] - origin[axis]) / direction[axis];
      double t1 = (box.hi[axis] - origin[axis]) / direction[axis];
      if (t0 > t1)
        std::swap(t0, t1);
      lo = std::max(lo, t0);
      hi = std::min(hi, t1);
      if (lo > hi)
        return false;
    }
    near_t = lo;
    return true;
  }

  static bool ray_triangle(
    const Vec3d& origin, const Vec3d& direction, const RayTriangle& tri,
    double min_t, double& t, double weights[3])
  {
    const Vec3d e1 = tri.p[1] - tri.p[0];
    const Vec3d e2 = tri.p[2] - tri.p[0];
    const Vec3d pvec = direction.cross(e2);
    const double det = e1 | pvec;
    const double det_scale = std::max(e1.length() * e2.length(), 1e-300);
    if (std::abs(det) <= 1e-14 * det_scale)
      return false;

    const double inv_det = 1.0 / det;
    const Vec3d tvec = origin - tri.p[0];
    const double u = (tvec | pvec) * inv_det;
    const double bary_eps = 1e-10;
    if (u < -bary_eps || u > 1.0 + bary_eps)
      return false;

    const Vec3d qvec = tvec.cross(e1);
    const double v = (direction | qvec) * inv_det;
    if (v < -bary_eps || u + v > 1.0 + bary_eps)
      return false;

    t = (e2 | qvec) * inv_det;
    if (!(t > min_t) || !std::isfinite(t))
      return false;

    weights[0] = 1.0 - u - v;
    weights[1] = u;
    weights[2] = v;
    return true;
  }
};

struct AnchorRequest
{
  size_t loop_index = 0;
  size_t order_index = 0;
  Vec3d point;
  std::optional<Geometry::Point_3> exact_point;
  // This is always the co-normal of the corresponding source edge.  In the
  // vertex benchmark the ray itself follows the bisector of two co-normals,
  // so keeping the support direction separate is essential for Phase 3.
  Vec3d support_outer_direction;
  FaceHandle source_face;
  enum class Location { Vertex, Edge, Face } location = Location::Face;
  VertexHandle vertex;
  EdgeHandle edge;
  VertexHandle edge_start;
  VertexHandle edge_end;
  double edge_parameter = 0.0;
  VertexHandle inserted_vertex;
};

struct BoundaryLoop
{
  size_t detected_loop_index = 0;
  std::vector<HalfedgeHandle> halfedges;
  std::vector<size_t> request_indices;
  std::vector<VertexHandle> anchors;
};

constexpr uint64_t kInvalidWeldedEdge = std::numeric_limits<uint64_t>::max();

// Folds the source vertices that share one position into a single index.  A
// half-edge mesh cannot store a non-manifold vertex, so OpenMesh's importer
// duplicates the vertices of every face it fails to attach; welding by
// position undoes that split and recovers the connectivity of the input file.
std::vector<int> weld_vertices_by_position(SMeshT& mesh)
{
  std::vector<int> welded(mesh.n_vertices(), -1);
  std::vector<int> order;
  order.reserve(mesh.n_vertices());
  for (VertexHandle vh : mesh.vertices())
    order.push_back(vh.idx());

  // Exact comparison is what is wanted here: the importer copies the position
  // of the vertex it duplicates verbatim, so the two carry identical bits.
  const auto position_less = [&mesh](int a, int b)
  {
    const Vec3d& pa = mesh.point(mesh.vertex_handle(a));
    const Vec3d& pb = mesh.point(mesh.vertex_handle(b));
    if (pa.x() != pb.x())
      return pa.x() < pb.x();
    if (pa.y() != pb.y())
      return pa.y() < pb.y();
    return pa.z() < pb.z();
  };
  std::sort(order.begin(), order.end(), position_less);

  int next_id = 0;
  for (size_t i = 0; i < order.size(); i++)
  {
    if (i > 0 && !position_less(order[i - 1], order[i]))
      welded[order[i]] = welded[order[i - 1]];
    else
      welded[order[i]] = next_id++;
  }
  return welded;
}

uint64_t welded_edge_key(
  const std::vector<int>& welded, VertexHandle a, VertexHandle b)
{
  if (!a.is_valid() || !b.is_valid())
    return kInvalidWeldedEdge;
  const int ia = welded[a.idx()];
  const int ib = welded[b.idx()];
  if (ia < 0 || ib < 0 || ia == ib)
    return kInvalidWeldedEdge;
  const uint64_t lo = static_cast<uint64_t>(std::min(ia, ib));
  const uint64_t hi = static_cast<uint64_t>(std::max(ia, ib));
  return (lo << 32) | hi;
}

// How many faces use each welded edge.  One face means the input surface is
// genuinely open along that edge; two or more mean it is closed there and the
// half-edge boundary is only the seam left by a duplicated vertex.
std::unordered_map<uint64_t, int> count_welded_edge_faces(
  SMeshT& mesh, const std::vector<int>& welded)
{
  std::unordered_map<uint64_t, int> counts;
  counts.reserve(mesh.n_edges() * 2);
  for (FaceHandle fh : mesh.faces())
  {
    VertexHandle vertices[3];
    size_t i = 0;
    for (VertexHandle vh : mesh.fv_range(fh))
    {
      if (i >= 3)
        break;
      vertices[i++] = vh;
    }
    if (i != 3)
      continue;
    for (size_t k = 0; k < 3; k++)
    {
      const uint64_t key =
        welded_edge_key(welded, vertices[k], vertices[(k + 1) % 3]);
      if (key != kInvalidWeldedEdge)
        counts[key]++;
    }
  }
  return counts;
}

bool get_boundary_edge_outer_direction(
  SMeshT& source, HalfedgeHandle boundary_halfedge,
  double length_epsilon, Vec3d& outer, const char*& failure_reason)
{
  failure_reason = "none";
  if (!source.is_boundary(boundary_halfedge))
  {
    failure_reason = "not_boundary_halfedge";
    return false;
  }
  const HalfedgeHandle inner = source.opposite_halfedge_handle(boundary_halfedge);
  const FaceHandle face = source.face_handle(inner);
  if (!face.is_valid())
  {
    failure_reason = "missing_incident_face";
    return false;
  }

  const VertexHandle v0 = source.from_vertex_handle(boundary_halfedge);
  const VertexHandle v1 = source.to_vertex_handle(boundary_halfedge);
  VertexHandle opposite;
  for (VertexHandle vh : source.fv_range(face))
  {
    if (vh != v0 && vh != v1)
    {
      opposite = vh;
      break;
    }
  }
  if (!opposite.is_valid())
  {
    failure_reason = "missing_opposite_vertex";
    return false;
  }

  const Vec3d p0 = source.point(v0);
  const Vec3d p1 = source.point(v1);
  const Vec3d edge = p1 - p0;
  const double edge_sqr = edge.squaredNorm();
  if (edge_sqr <= length_epsilon * length_epsilon)
  {
    failure_reason = "degenerate_edge";
    return false;
  }

  const Vec3d midpoint = (p0 + p1) * 0.5;
  const Vec3d toward_opposite = source.point(opposite) - midpoint;
  Vec3d inward = toward_opposite - edge * ((toward_opposite | edge) / edge_sqr);
  if (inward.length() <= length_epsilon)
  {
    failure_reason = "degenerate_triangle";
    return false;
  }
  inward.normalize();
  outer = -inward;
  return true;
}

// These helpers are called only after an existing acceptance test fails.
// Missing rays, vertices, and hits are omitted rather than read uninitialized.
void append_diagnostic_vector(
  std::ostream& out, const std::string& name, const Vec3d& value)
{
  out << ' ' << name << "=(" << value.x() << ',' << value.y() << ','
      << value.z() << ')';
}

void append_source_edge_diagnostics(
  std::ostream& out, SMeshT& source, HalfedgeHandle halfedge,
  const std::string& role)
{
  out << ' ' << role << "_halfedge=" << halfedge.idx();
  if (!halfedge.is_valid() ||
    static_cast<size_t>(halfedge.idx()) >= source.n_halfedges())
    return;
  const VertexHandle from = source.from_vertex_handle(halfedge);
  const VertexHandle to = source.to_vertex_handle(halfedge);
  out << ' ' << role << "_from_vertex=" << from.idx()
      << ' ' << role << "_to_vertex=" << to.idx();
  if (!from.is_valid() || !to.is_valid())
    return;
  const Vec3d p0 = source.point(from);
  const Vec3d p1 = source.point(to);
  append_diagnostic_vector(out, role + "_from_point", p0);
  append_diagnostic_vector(out, role + "_to_point", p1);
  const Vec3d edge = p1 - p0;
  const double edge_sqr = edge.squaredNorm();
  out << ' ' << role << "_edge_squared_length=" << edge_sqr;
  const FaceHandle face = source.face_handle(
    source.opposite_halfedge_handle(halfedge));
  out << ' ' << role << "_incident_face=" << face.idx();
  if (!face.is_valid())
    return;
  for (VertexHandle opposite : source.fv_range(face))
  {
    if (opposite == from || opposite == to)
      continue;
    out << ' ' << role << "_opposite_vertex=" << opposite.idx();
    const Vec3d p2 = source.point(opposite);
    append_diagnostic_vector(out, role + "_opposite_point", p2);
    if (edge_sqr > 0.0 && std::isfinite(edge_sqr))
    {
      const Vec3d toward_opposite = p2 - (p0 + p1) * 0.5;
      const Vec3d inward =
        toward_opposite - edge * ((toward_opposite | edge) / edge_sqr);
      out << ' ' << role << "_triangle_height=" << inward.length();
    }
    break;
  }
}

void append_ray_hit_diagnostics(
  std::ostream& out, const RayTriangleBvh& ray_tree, const RayHit& hit,
  const std::string& role)
{
  out << ' ' << role << "_found=" << hit.hit;
  if (!hit.hit)
    return;
  out << ' ' << role << "_t=" << hit.t
      << ' ' << role << "_finite_point=" << finite_vec(hit.point)
      << ' ' << role << "_weights=(" << hit.weights[0] << ','
      << hit.weights[1] << ',' << hit.weights[2] << ')';
  append_diagnostic_vector(out, role + "_point", hit.point);
  const RayTriangle& triangle = ray_tree.triangle(hit.triangle);
  out << ' ' << role << "_triangle=" << hit.triangle
      << ' ' << role << "_cage_face=" << triangle.face.idx();
  for (size_t i = 0; i < 3; i++)
  {
    const std::string vertex_name = role + "_cage_vertex" + std::to_string(i);
    out << ' ' << vertex_name << "_id=" << triangle.vertices[i].idx();
    append_diagnostic_vector(out, vertex_name + "_point", triangle.p[i]);
  }
}

EdgeHandle edge_between(SMeshT& mesh, VertexHandle a, VertexHandle b)
{
  HalfedgeHandle halfedge = mesh.find_halfedge(a, b);
  if (!halfedge.is_valid())
    halfedge = mesh.find_halfedge(b, a);
  return halfedge.is_valid() ? mesh.edge_handle(halfedge) : EdgeHandle();
}

struct GraphPath
{
  std::vector<VertexHandle> vertices;
  double cost = DBL_MAX;
};

struct RailSearchDeadline
{
  using Clock = std::chrono::steady_clock;
  const Clock::time_point* until = nullptr;
  bool timedOut = false;

  bool expired()
  {
    if (until && Clock::now() >= *until) timedOut = true;
    return timedOut;
  }
};

bool shortest_path(
  SMeshT& mesh, VertexHandle source, VertexHandle target,
  const std::unordered_set<int>& blocked_vertices,
  const std::unordered_set<int>& blocked_edges,
  GraphPath& path, RailSearchDeadline& deadline)
{
  const size_t vertex_capacity = mesh.n_vertices();
  std::vector<double> distance(vertex_capacity, DBL_MAX);
  std::vector<int> previous(vertex_capacity, -1);
  using QueueEntry = std::pair<double, int>;
  std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;
  distance[source.idx()] = 0.0;
  queue.emplace(0.0, source.idx());

  size_t queue_pops = 0;
  while (!queue.empty())
  {
    if ((queue_pops++ & 255u) == 0u && deadline.expired())
      return false;
    const auto entry = queue.top();
    queue.pop();
    const double current_distance = entry.first;
    const VertexHandle current(entry.second);
    if (current_distance != distance[current.idx()])
      continue;
    if (current == target)
      break;

    for (HalfedgeHandle outgoing : mesh.voh_range(current))
    {
      const EdgeHandle edge = mesh.edge_handle(outgoing);
      const VertexHandle next = mesh.to_vertex_handle(outgoing);
      if (mesh.status(edge).deleted() || mesh.status(next).deleted())
        continue;
      if (blocked_edges.count(edge.idx()) > 0)
        continue;
      if (next != target && blocked_vertices.count(next.idx()) > 0)
        continue;
      if (mesh.data(edge).boundary_rail_id != kNoRail)
        continue;
      if (next != target && mesh.data(next).boundary_rail_id != kNoRail)
        continue;

      const double edge_length =
        (mesh.point(current) - mesh.point(next)).length();
      const double candidate = current_distance + edge_length;
      if (candidate < distance[next.idx()])
      {
        distance[next.idx()] = candidate;
        previous[next.idx()] = current.idx();
        queue.emplace(candidate, next.idx());
      }
    }
  }

  if (!std::isfinite(distance[target.idx()]))
    return false;

  path.vertices.clear();
  for (int current = target.idx(); current >= 0; current = previous[current])
  {
    path.vertices.emplace_back(current);
    if (current == source.idx())
      break;
  }
  if (path.vertices.empty() || path.vertices.back() != source)
    return false;
  std::reverse(path.vertices.begin(), path.vertices.end());
  path.cost = distance[target.idx()];
  return true;
}

std::vector<GraphPath> candidate_paths(
  SMeshT& mesh, VertexHandle source, VertexHandle target,
  const std::unordered_set<int>& blocked_vertices,
  const std::unordered_set<int>& blocked_edges,
  size_t max_candidates, RailSearchDeadline& deadline)
{
  std::vector<GraphPath> result;
  GraphPath shortest;
  if (!shortest_path(
    mesh, source, target, blocked_vertices, blocked_edges, shortest, deadline))
    return result;
  if (deadline.expired()) return {};
  result.push_back(shortest);

  const size_t edge_count = shortest.vertices.size() > 0 ?
    shortest.vertices.size() - 1 : 0;
  const size_t probes = std::min<size_t>(edge_count, 32);
  for (size_t probe = 0; probe < probes && result.size() < max_candidates; probe++)
  {
    if (deadline.expired()) return {};
    const size_t i = probes == edge_count ? probe :
      (probe * edge_count) / probes;
    const EdgeHandle excluded = edge_between(
      mesh, shortest.vertices[i], shortest.vertices[i + 1]);
    if (!excluded.is_valid())
      continue;
    std::unordered_set<int> trial_blocked_edges = blocked_edges;
    trial_blocked_edges.insert(excluded.idx());
    GraphPath alternative;
    const bool found = shortest_path(
      mesh, source, target, blocked_vertices,
      trial_blocked_edges, alternative, deadline);
    if (deadline.expired()) return {};
    if (!found) continue;
    bool duplicate = false;
    for (const GraphPath& existing : result)
    {
      if (existing.vertices == alternative.vertices)
      {
        duplicate = true;
        break;
      }
    }
    if (!duplicate)
      result.push_back(std::move(alternative));
  }
  std::sort(result.begin(), result.end(),
    [](const GraphPath& lhs, const GraphPath& rhs)
    {
      return lhs.cost < rhs.cost;
    });
  if (result.size() > max_candidates)
    result.resize(max_candidates);
  if (deadline.expired()) return {};
  return result;
}

bool validate_cycle(
  SMeshT& mesh, const std::vector<GraphPath>& paths,
  const std::vector<VertexHandle>& anchors,
  std::vector<VertexHandle>& cycle_vertices,
  std::vector<EdgeHandle>& cycle_edges)
{
  std::unordered_map<int, std::vector<int>> adjacency;
  std::set<int> unique_edges;
  cycle_edges.clear();
  for (const GraphPath& path : paths)
  {
    for (size_t i = 0; i + 1 < path.vertices.size(); i++)
    {
      const VertexHandle a = path.vertices[i];
      const VertexHandle b = path.vertices[i + 1];
      const EdgeHandle edge = edge_between(mesh, a, b);
      if (!edge.is_valid() || !unique_edges.insert(edge.idx()).second)
        return false;
      cycle_edges.push_back(edge);
      adjacency[a.idx()].push_back(b.idx());
      adjacency[b.idx()].push_back(a.idx());
    }
  }

  if (adjacency.size() < 3 || cycle_edges.size() != adjacency.size())
    return false;
  for (const auto& entry : adjacency)
  {
    if (entry.second.size() != 2)
      return false;
  }
  for (VertexHandle anchor : anchors)
  {
    if (adjacency.count(anchor.idx()) == 0)
      return false;
  }

  cycle_vertices.clear();
  const int start = anchors.front().idx();
  int previous = -1;
  int current = start;
  do
  {
    cycle_vertices.emplace_back(current);
    const std::vector<int>& neighbors = adjacency[current];
    const int next = neighbors[0] == previous ? neighbors[1] : neighbors[0];
    previous = current;
    current = next;
    if (cycle_vertices.size() > adjacency.size())
      return false;
  } while (current != start);

  if (cycle_vertices.size() != adjacency.size())
    return false;

  std::vector<int> encountered_anchors;
  std::unordered_set<int> anchor_set;
  for (VertexHandle anchor : anchors)
    anchor_set.insert(anchor.idx());
  for (VertexHandle vh : cycle_vertices)
  {
    if (anchor_set.count(vh.idx()) > 0)
      encountered_anchors.push_back(vh.idx());
  }
  if (encountered_anchors.size() != anchors.size())
    return false;

  const auto cyclic_order_matches = [&](bool reverse)
  {
    const size_t count = anchors.size();
    size_t offset = count;
    for (size_t i = 0; i < count; i++)
    {
      if (encountered_anchors[i] == anchors.front().idx())
      {
        offset = i;
        break;
      }
    }
    if (offset == count)
      return false;
    for (size_t i = 0; i < count; i++)
    {
      const size_t expected = reverse ? (count - i) % count : i;
      if (encountered_anchors[(offset + i) % count] != anchors[expected].idx())
        return false;
    }
    return true;
  };
  return cyclic_order_matches(false) || cyclic_order_matches(true);
}

bool construct_rail_attempt(
  SMeshT& mesh, const std::vector<VertexHandle>& input_anchors,
  const std::unordered_set<int>& all_anchor_vertices,
  int rail_id, BoundaryRailLoop& result, size_t candidate_limit,
  RailSearchDeadline& deadline, BoundaryRailSearchStats& stats)
{
  stats = BoundaryRailSearchStats{};
  std::vector<VertexHandle> anchors;
  anchors.reserve(input_anchors.size());
  for (VertexHandle anchor : input_anchors)
  {
    if (anchors.empty() || anchors.back() != anchor)
      anchors.push_back(anchor);
  }
  if (anchors.size() > 1 && anchors.front() == anchors.back())
    anchors.pop_back();
  std::unordered_set<int> unique_anchors;
  for (VertexHandle anchor : anchors)
    unique_anchors.insert(anchor.idx());
  if (anchors.size() < 3 || unique_anchors.size() != anchors.size())
  {
    stats.outcome = anchors.size() < 3 ? "too_few_anchors" : "repeated_anchor";
    return false;
  }
  // A ray hit may reuse a cage vertex that an earlier boundary component has
  // already claimed.  Endpoints are intentionally exempt from Dijkstra's
  // intermediate-vertex blocking, so reject that loop here instead of letting
  // the later rail overwrite the earlier rail id on a shared anchor.
  for (VertexHandle anchor : anchors)
  {
    if (mesh.data(anchor).boundary_rail_id != kNoRail)
    {
      stats.outcome = "anchor_already_claimed";
      return false;
    }
  }

  std::unordered_set<int> used_vertices;
  std::unordered_set<int> used_edges;
  std::vector<GraphPath> selected(anchors.size());
  size_t search_nodes = 0;
  const size_t max_search_nodes = std::max<size_t>(1024, anchors.size() * 16);

  std::function<bool(size_t)> select_paths = [&](size_t segment)
  {
    if (deadline.expired()) return false;
    if (++search_nodes > max_search_nodes)
      return false;
    if (segment == anchors.size())
      return true;

    const VertexHandle source = anchors[segment];
    const VertexHandle target = anchors[(segment + 1) % anchors.size()];
    std::unordered_set<int> blocked_vertices = used_vertices;
    for (int anchor : all_anchor_vertices)
    {
      if (anchor != source.idx() && anchor != target.idx())
        blocked_vertices.insert(anchor);
    }

    const std::vector<GraphPath> candidates = candidate_paths(
      mesh, source, target, blocked_vertices, used_edges, candidate_limit, deadline);
    if (deadline.timedOut) return false;
    for (const GraphPath& candidate : candidates)
    {
      std::vector<int> added_vertices;
      std::vector<int> added_edges;
      bool valid = true;
      for (size_t i = 1; i + 1 < candidate.vertices.size(); i++)
      {
        const int vertex = candidate.vertices[i].idx();
        if (!used_vertices.insert(vertex).second)
        {
          valid = false;
          break;
        }
        added_vertices.push_back(vertex);
      }
      if (valid)
      {
        for (size_t i = 0; i + 1 < candidate.vertices.size(); i++)
        {
          const EdgeHandle edge = edge_between(
            mesh, candidate.vertices[i], candidate.vertices[i + 1]);
          if (!edge.is_valid() || !used_edges.insert(edge.idx()).second)
          {
            valid = false;
            break;
          }
          added_edges.push_back(edge.idx());
        }
      }

      if (valid)
      {
        selected[segment] = candidate;
        if (select_paths(segment + 1))
          return true;
      }
      for (int edge : added_edges)
        used_edges.erase(edge);
      for (int vertex : added_vertices)
        used_vertices.erase(vertex);
      if (deadline.timedOut) return false;
    }
    return false;
  };

  const bool connected = select_paths(0);
  stats.recursiveCalls = search_nodes;
  stats.exploredStates = std::min(search_nodes, max_search_nodes);
  stats.stateLimitHits = search_nodes > max_search_nodes ? 1 : 0;
  stats.timedOut = deadline.timedOut;
  if (!connected)
  {
    stats.outcome = stats.timedOut ? "time_limit" :
      (stats.stateLimitHits ? "state_limit" : "candidates_exhausted");
    return false;
  }

  std::vector<VertexHandle> cycle_vertices;
  std::vector<EdgeHandle> cycle_edges;
  if (!validate_cycle(mesh, selected, anchors, cycle_vertices, cycle_edges))
  {
    stats.outcome = "cycle_validation_failed";
    return false;
  }
  if (deadline.expired())
  {
    stats.timedOut = true;
    stats.outcome = "time_limit";
    return false;
  }

  for (VertexHandle vh : cycle_vertices)
    mesh.data(vh).boundary_rail_id = rail_id;
  for (EdgeHandle eh : cycle_edges)
    mesh.data(eh).boundary_rail_id = rail_id;
  result.railId = rail_id;
  result.vertices = std::move(cycle_vertices);
  result.anchors = std::move(anchors);
  stats.outcome = "success";
  return true;
}
}// namespace

bool detail::construct_boundary_rail(
  SMeshT& mesh, const std::vector<VertexHandle>& input_anchors,
  const std::unordered_set<int>& all_anchor_vertices,
  int rail_id, BoundaryRailLoop& result,
  const BoundaryRailSearchOptions& options, BoundaryRailSearchStats& stats)
{
  if (options.candidateLimit == 0)
    throw std::invalid_argument("boundary rail candidate limit must be positive");
  if (!std::isfinite(options.maxSeconds) || options.maxSeconds < 0.0)
    throw std::invalid_argument("boundary rail search seconds must be finite and nonnegative");
  stats = BoundaryRailSearchStats{};
  const auto started = RailSearchDeadline::Clock::now();
  const auto finish = [&](bool success)
  {
    stats.seconds = std::chrono::duration<double>(RailSearchDeadline::Clock::now() - started).count();
    return success;
  };
  std::vector<VertexHandle> anchors;
  for (VertexHandle anchor : input_anchors)
    if (anchors.empty() || anchors.back() != anchor) anchors.push_back(anchor);
  if (anchors.size() > 1 && anchors.front() == anchors.back()) anchors.pop_back();
  stats.anchors = anchors.size();
  std::unordered_set<int> unique;
  for (VertexHandle anchor : anchors) unique.insert(anchor.idx());
  if (anchors.size() < 3 || unique.size() != anchors.size())
  {
    stats.outcome = anchors.size() < 3 ? "too_few_anchors" : "repeated_anchor";
    return finish(false);
  }
  for (VertexHandle anchor : anchors)
    if (mesh.data(anchor).boundary_rail_id != kNoRail)
    {
      stats.outcome = "anchor_already_claimed";
      return finish(false);
    }
  // Original start first; every later attempt is a forward cyclic shift.
  // The graph, global reservations, and previously completed rails stay fixed.
  auto end = RailSearchDeadline::Clock::time_point::max();
  if (options.maxSeconds > 0.0)
  {
    const double available = std::chrono::duration<double>(end - started).count();
    if (options.maxSeconds < available)
      end = started + std::chrono::duration_cast<RailSearchDeadline::Clock::duration>(
        std::chrono::duration<double>(options.maxSeconds));
  }
  RailSearchDeadline deadline;
  deadline.until = options.maxSeconds > 0.0 ? &end : nullptr;
  const size_t starts = options.retryCyclicStarts ? anchors.size() : 1;
  for (size_t shift = 0; shift < starts; ++shift)
  {
    if (deadline.expired())
    {
      stats.timedOut = true;
      stats.outcome = "time_limit";
      break;
    }
    auto trial_anchors = anchors;
    std::rotate(trial_anchors.begin(), trial_anchors.begin() + shift, trial_anchors.end());
    BoundaryRailLoop trial_result;
    BoundaryRailSearchStats attempt;
    ++stats.attempts;
    const bool success = construct_rail_attempt(mesh, trial_anchors, all_anchor_vertices,
      rail_id, trial_result, options.candidateLimit, deadline, attempt);
    stats.exploredStates += attempt.exploredStates;
    stats.recursiveCalls += attempt.recursiveCalls;
    stats.stateLimitHits += attempt.stateLimitHits;
    stats.timedOut = attempt.timedOut;
    stats.outcome = attempt.outcome;
    if (success)
    {
      stats.successfulShift = static_cast<int>(shift);
      result = std::move(trial_result);
      return finish(true);
    }
    if (stats.timedOut) break;
  }
  return finish(false);
}

size_t count_welded_boundary_edges(SMeshT& mesh)
{
  const std::vector<int> welded = weld_vertices_by_position(mesh);
  size_t boundary_edges = 0;
  for (const auto& edge : count_welded_edge_faces(mesh, welded))
    boundary_edges += edge.second == 1;
  return boundary_edges;
}

bool BoundaryRailBuilder::build()
{
  return build_with_stats().successful;
}

BoundaryRailBuildStats BoundaryRailBuilder::build_with_stats()
{
  BoundaryRailBuildStats stats;
  if (!source || !cage || source->n_vertices() == 0 || cage->n_faces() == 0)
    return stats;

  // The source-edge labels and directions are consumed by EdgeCollapser in
  // Phase 3.  Reset them in case a mesh instance is reused for another build.
  for (EdgeHandle eh : source->edges())
  {
    source->data(eh).boundary_rail_id = kNoRail;
    source->data(eh).boundary_rail_outer_direction = Vec3d(0.0, 0.0, 0.0);
  }

  Box source_box;
  for (VertexHandle vh : source->vertices())
    source_box.add(source->point(vh));
  Box cage_box;
  for (VertexHandle vh : cage->vertices())
    cage_box.add(cage->point(vh));
  const double source_scale = std::max((source_box.hi - source_box.lo).length(), 1.0);
  const double cage_scale = std::max((cage_box.hi - cage_box.lo).length(), 1.0);
  const double length_epsilon = source_scale * 1e-12;
  // Source rays can meet the enclosing cage arbitrarily close to their origin.
  // ray_triangle still rejects t <= 0 and non-finite intersection distances.
  const double min_ray_t = 0.0;
  const double anchor_merge_epsilon = cage_scale * 1e-9;

  // A half-edge mesh cannot represent a non-manifold vertex, so OpenMesh's
  // importer duplicates the vertices of every face it cannot attach.  That
  // tears a closed source surface open around the singular vertex, and to
  // is_boundary() the tear is indistinguishable from a real hole: a watertight
  // input with one bow-tie vertex would otherwise grow rails around nothing.
  // Welding the source vertices back together by position separates the two.
  // A half-edge borders a hole only when its welded edge carries one face;
  // when it carries two or more the surface is closed there and the boundary
  // is a seam, not a source boundary loop.
  const std::vector<int> welded_vertices = weld_vertices_by_position(*source);
  const std::unordered_map<uint64_t, int> welded_edge_faces =
    count_welded_edge_faces(*source, welded_vertices);
  const auto borders_hole = [&](HalfedgeHandle halfedge)
  {
    const uint64_t key = welded_edge_key(
      welded_vertices,
      source->from_vertex_handle(halfedge),
      source->to_vertex_handle(halfedge));
    if (key == kInvalidWeldedEdge)
      return false;
    const auto found = welded_edge_faces.find(key);
    return found != welded_edge_faces.end() && found->second == 1;
  };

  std::vector<BoundaryLoop> loops;
  size_t seam_loops = 0;
  std::vector<bool> visited(source->n_halfedges(), false);
  for (HalfedgeHandle start : source->halfedges())
  {
    if (!source->is_boundary(start) || visited[start.idx()])
      continue;
    BoundaryLoop loop;
    HalfedgeHandle current = start;
    bool valid = true;
    do
    {
      if (!current.is_valid() || !source->is_boundary(current) ||
        visited[current.idx()])
      {
        valid = current == start;
        break;
      }
      visited[current.idx()] = true;
      loop.halfedges.push_back(current);
      current = source->next_halfedge_handle(current);
      if (loop.halfedges.size() > source->n_halfedges())
      {
        valid = false;
        break;
      }
    } while (current != start);
    if (!valid || current != start || loop.halfedges.size() < 3)
      continue;

    // Every edge of the loop has to bound a hole.  A loop that mixes the two
    // runs along a seam for part of its length and does not describe a source
    // boundary that a rail could follow.
    if (!std::all_of(
        loop.halfedges.begin(), loop.halfedges.end(), borders_hole))
    {
      seam_loops++;
      continue;
    }
    loops.push_back(std::move(loop));
  }
  stats.detectedLoopCount = loops.size();

  if (seam_loops > 0)
  {
    Logger::user_logger->info(
      "boundary rail construction ignored {} boundary loops that are seams of a non-manifold vertex rather than holes in the source surface.",
      seam_loops);
  }

  if (loops.empty())
  {
    Logger::user_logger->info(
      "boundary rail construction: source mesh has no closed boundary loop.");
    stats.successful = true;
    return stats;
  }

  RayTriangleBvh ray_tree(*cage);
  std::vector<AnchorRequest> requests;
  std::vector<BoundaryLoop> valid_loops;
  valid_loops.reserve(loops.size());
  for (size_t loop_index = 0; loop_index < loops.size(); loop_index++)
  {
    BoundaryLoop candidate = loops[loop_index];
    candidate.detected_loop_index = loop_index;
    std::vector<AnchorRequest> loop_requests;
    loop_requests.reserve(candidate.halfedges.size());
    bool valid = true;
    for (size_t i = 0; i < candidate.halfedges.size(); i++)
    {
      const HalfedgeHandle boundary_halfedge = candidate.halfedges[i];
      const auto failure_context = [&](const char* reason)
      {
        std::ostringstream out;
        out << std::setprecision(17)
            << "boundary rail anchor failure: mode="
            << (mode == BoundaryRailAnchorMode::EdgeMidpoint ?
                "edge-midpoint" : "vertex-bisector")
            << " detected_loop=" << loop_index << " anchor=" << i
            << " loop_anchor_count=" << candidate.halfedges.size()
            << " reason=" << reason
            << " source_scale=" << source_scale << " cage_scale=" << cage_scale
            << " length_epsilon=" << length_epsilon
            << " min_ray_t=" << min_ray_t
            << " anchor_merge_epsilon=" << anchor_merge_epsilon;
        append_source_edge_diagnostics(out, *source, boundary_halfedge, "current");
        if (mode != BoundaryRailAnchorMode::EdgeMidpoint)
        {
          const HalfedgeHandle incoming = candidate.halfedges[
            (i + candidate.halfedges.size() - 1) % candidate.halfedges.size()];
          append_source_edge_diagnostics(out, *source, incoming, "incoming");
        }
        return out;
      };
      Vec3d support_outer_direction;
      const char* conormal_failure = "none";
      if (!get_boundary_edge_outer_direction(
          *source, boundary_halfedge, length_epsilon,
          support_outer_direction, conormal_failure))
      {
        PhaseTimer::Exclusion exclusion("diagnostic");
        auto diagnostic = failure_context("current_edge_conormal_failed");
        diagnostic << " conormal_failure=" << conormal_failure;
        Logger::user_logger->warn("{}", diagnostic.str());
        valid = false;
        break;
      }

      const VertexHandle source_from =
        source->from_vertex_handle(boundary_halfedge);
      const VertexHandle source_to =
        source->to_vertex_handle(boundary_halfedge);
      Vec3d origin;
      Vec3d ray_direction;
      if (mode == BoundaryRailAnchorMode::EdgeMidpoint)
      {
        // One anchor corresponds to one source boundary edge.
        origin =
          (source->point(source_from) + source->point(source_to)) * 0.5;
        ray_direction = support_outer_direction;
      }
      else
      {
        // Reproduce the vertex-origin alternative from commit 759563c: one
        // anchor starts at the current boundary vertex and follows the
        // normalized sum of the incoming and outgoing edge co-normals.
        const HalfedgeHandle incoming = candidate.halfedges[
          (i + candidate.halfedges.size() - 1) % candidate.halfedges.size()];
        Vec3d incoming_outer_direction;
        if (!get_boundary_edge_outer_direction(
            *source, incoming, length_epsilon,
            incoming_outer_direction, conormal_failure))
        {
          PhaseTimer::Exclusion exclusion("diagnostic");
          auto diagnostic = failure_context("incoming_edge_conormal_failed");
          diagnostic << " conormal_failure=" << conormal_failure;
          append_diagnostic_vector(diagnostic, "current_conormal", support_outer_direction);
          Logger::user_logger->warn("{}", diagnostic.str());
          valid = false;
          break;
        }
        ray_direction =
          incoming_outer_direction + support_outer_direction;
        // Both co-normals are unit vectors, so cancellation is an angular
        // condition and must not depend on the source mesh's physical scale.
        if (ray_direction.length() <= 1e-12)
        {
          PhaseTimer::Exclusion exclusion("diagnostic");
          auto diagnostic = failure_context("vertex_conormal_cancellation");
          diagnostic << " cancellation_epsilon=" << 1e-12
                     << " conormal_sum_length=" << ray_direction.length();
          append_diagnostic_vector(diagnostic, "current_conormal", support_outer_direction);
          append_diagnostic_vector(diagnostic, "incoming_conormal", incoming_outer_direction);
          append_diagnostic_vector(diagnostic, "conormal_sum", ray_direction);
          append_diagnostic_vector(diagnostic, "ray_origin", source->point(source_from));
          Logger::user_logger->warn("{}", diagnostic.str());
          valid = false;
          break;
        }
        ray_direction.normalize();
        origin = source->point(source_from);
      }

      RayHit hit;
      if (!ray_tree.first_intersection(
          origin, ray_direction, min_ray_t, hit))
      {
        PhaseTimer::Exclusion exclusion("diagnostic");
        auto diagnostic = failure_context("ray_no_hit");
        append_diagnostic_vector(diagnostic, "ray_origin", origin);
        append_diagnostic_vector(diagnostic, "ray_direction", ray_direction);
        append_diagnostic_vector(diagnostic, "cage_box_min", cage_box.lo);
        append_diagnostic_vector(diagnostic, "cage_box_max", cage_box.hi);
        Logger::user_logger->warn("{}", diagnostic.str());
        valid = false;
        break;
      }
      if (!finite_vec(hit.point))
      {
        PhaseTimer::Exclusion exclusion("diagnostic");
        auto diagnostic = failure_context("ray_nonfinite_hit");
        append_diagnostic_vector(diagnostic, "ray_origin", origin);
        append_diagnostic_vector(diagnostic, "ray_direction", ray_direction);
        append_ray_hit_diagnostics(diagnostic, ray_tree, hit, "hit");
        Logger::user_logger->warn("{}", diagnostic.str());
        valid = false;
        break;
      }

      const RayTriangle& triangle = ray_tree.triangle(hit.triangle);
      AnchorRequest request;
      request.loop_index = valid_loops.size();
      request.order_index = i;
      request.point = hit.point;
      request.support_outer_direction = support_outer_direction;
      request.source_face = triangle.face;
      const double bary_epsilon = 1e-8;
      std::vector<size_t> near_zero;
      for (size_t j = 0; j < 3; j++)
      {
        if (std::abs(hit.weights[j]) <= bary_epsilon)
          near_zero.push_back(j);
      }
      if (near_zero.size() >= 2)
      {
        size_t vertex_index = 0;
        if (hit.weights[1] > hit.weights[vertex_index])
          vertex_index = 1;
        if (hit.weights[2] > hit.weights[vertex_index])
          vertex_index = 2;
        request.location = AnchorRequest::Location::Vertex;
        request.vertex = triangle.vertices[vertex_index];
      }
      else if (near_zero.size() == 1)
      {
        const size_t zero = near_zero.front();
        const size_t a = (zero + 1) % 3;
        const size_t b = (zero + 2) % 3;
        request.location = AnchorRequest::Location::Edge;
        request.edge = edge_between(*cage, triangle.vertices[a], triangle.vertices[b]);
        if (!request.edge.is_valid())
        {
          PhaseTimer::Exclusion exclusion("diagnostic");
          auto diagnostic = failure_context("cage_edge_lookup_failed");
          diagnostic << " cage_edge_from_vertex=" << triangle.vertices[a].idx()
                     << " cage_edge_to_vertex=" << triangle.vertices[b].idx()
                     << " bary_epsilon=" << bary_epsilon;
          append_diagnostic_vector(diagnostic, "ray_origin", origin);
          append_diagnostic_vector(diagnostic, "ray_direction", ray_direction);
          append_ray_hit_diagnostics(diagnostic, ray_tree, hit, "hit");
          Logger::user_logger->warn("{}", diagnostic.str());
          valid = false;
          break;
        }
        const HalfedgeHandle edge_halfedge = cage->halfedge_handle(request.edge, 0);
        request.edge_start = cage->from_vertex_handle(edge_halfedge);
        request.edge_end = cage->to_vertex_handle(edge_halfedge);
        const Vec3d edge_vector =
          cage->point(request.edge_end) - cage->point(request.edge_start);
        request.edge_parameter = clamp01(
          ((request.point - cage->point(request.edge_start)) | edge_vector) /
          std::max(edge_vector.squaredNorm(), 1e-300));
        const auto exact_vertex = [&](VertexHandle vertex)
        {
          const auto* ep = cage->data(vertex).ep.get();
          const Vec3d& p = cage->point(vertex);
          return ep ? ep->exact() : Geometry::Point_3(p[0], p[1], p[2]);
        };
        const Geometry::Point_3 exact_start = exact_vertex(request.edge_start);
        const Geometry::Point_3 exact_end = exact_vertex(request.edge_end);
        // Preserve the original exact edge, including offsets too small to
        // retain in a double ray hit. Use the original endpoints before any
        // requests split this edge; edge_parameter is relative to those ends.
        request.exact_point = exact_start +
          (exact_end - exact_start) * Geometry::ET(request.edge_parameter);
        request.point = Geometry::ExactPoint(*request.exact_point).approx();
        if (request.edge_parameter == 0.0 || request.edge_parameter == 1.0)
        {
          request.location = AnchorRequest::Location::Vertex;
          request.vertex = request.edge_parameter == 0.0 ?
            request.edge_start : request.edge_end;
        }
      }
      else
        request.location = AnchorRequest::Location::Face;
      loop_requests.push_back(request);
    }

    if (!valid)
      continue;
    candidate.request_indices.clear();
    for (AnchorRequest& request : loop_requests)
    {
      candidate.request_indices.push_back(requests.size());
      requests.push_back(std::move(request));
    }
    valid_loops.push_back(std::move(candidate));
  }

  stats.rayValidLoopCount = valid_loops.size();
  stats.anchorRequestCount = requests.size();

  if (valid_loops.empty())
    return stats;

  const auto insertion_failure_context = [&](const char* reason, size_t request_index)
  {
    const AnchorRequest& request = requests[request_index];
    std::ostringstream out;
    out << std::setprecision(17)
        << "boundary rail insertion failure: mode="
        << (mode == BoundaryRailAnchorMode::EdgeMidpoint ?
            "edge-midpoint" : "vertex-bisector")
        << " reason=" << reason << " request=" << request_index
        << " detected_loop=" << valid_loops[request.loop_index].detected_loop_index
        << " anchor=" << request.order_index
        << " source_cage_face=" << request.source_face.idx()
        << " anchor_merge_epsilon=" << anchor_merge_epsilon;
    append_diagnostic_vector(out, "anchor_point", request.point);
    return out;
  };

  OpenMesh::FPropHandleT<int> source_face_property;
  cage->add_property(source_face_property, "boundary_rail_source_face");
  for (FaceHandle fh : cage->faces())
    cage->property(source_face_property, fh) = fh.idx();

  for (AnchorRequest& request : requests)
  {
    if (request.location == AnchorRequest::Location::Vertex)
      request.inserted_vertex = request.vertex;
  }

  std::unordered_map<int, std::vector<size_t>> edge_requests;
  for (size_t i = 0; i < requests.size(); i++)
  {
    if (requests[i].location == AnchorRequest::Location::Edge)
      edge_requests[requests[i].edge.idx()].push_back(i);
  }
  bool insertion_ok = true;
  for (auto& group : edge_requests)
  {
    std::vector<size_t>& indices = group.second;
    std::sort(indices.begin(), indices.end(), [&](size_t lhs, size_t rhs)
    {
      return requests[lhs].edge_parameter < requests[rhs].edge_parameter;
    });
    VertexHandle left = requests[indices.front()].edge_start;
    const VertexHandle right = requests[indices.front()].edge_end;
    VertexHandle previous_inserted;
    Vec3d previous_point;
    for (size_t request_index : indices)
    {
      AnchorRequest& request = requests[request_index];
      if (previous_inserted.is_valid() &&
        (request.point - previous_point).length() <= anchor_merge_epsilon)
      {
        request.inserted_vertex = previous_inserted;
        continue;
      }
      const EdgeHandle current_edge = edge_between(*cage, left, right);
      if (!current_edge.is_valid())
      {
        PhaseTimer::Exclusion exclusion("diagnostic");
        auto diagnostic = insertion_failure_context("split_edge_not_found", request_index);
        diagnostic << " original_cage_edge=" << request.edge.idx()
                   << " from_vertex=" << left.idx() << " to_vertex=" << right.idx();
        Logger::user_logger->error("{}", diagnostic.str());
        insertion_ok = false;
        break;
      }
      const VertexHandle inserted = cage->add_vertex(request.point);
      cage->data(inserted).ep =
        std::make_unique<Geometry::ExactPoint>(*request.exact_point);
      cage->split_copy(current_edge, inserted);
      request.inserted_vertex = inserted;
      previous_inserted = inserted;
      previous_point = request.point;
      left = inserted;
    }
    if (!insertion_ok)
      break;
  }

  std::vector<std::vector<FaceHandle>> descendants;
  int max_source_face = -1;
  for (FaceHandle fh : cage->faces())
    max_source_face = std::max(max_source_face, cage->property(source_face_property, fh));
  descendants.clear();
  descendants.resize(static_cast<size_t>(max_source_face + 1));
  for (FaceHandle fh : cage->faces())
  {
    const int source_face = cage->property(source_face_property, fh);
    if (source_face >= 0)
      descendants[source_face].push_back(fh);
  }

  std::vector<size_t> face_request_indices;
  for (size_t i = 0; i < requests.size(); i++)
  {
    if (requests[i].location == AnchorRequest::Location::Face)
      face_request_indices.push_back(i);
  }
  for (size_t request_index : face_request_indices)
  {
    if (!insertion_ok)
      break;
    AnchorRequest& request = requests[request_index];
    for (const AnchorRequest& existing : requests)
    {
      if (existing.inserted_vertex.is_valid() &&
        (cage->point(existing.inserted_vertex) - request.point).length() <= anchor_merge_epsilon)
      {
        request.inserted_vertex = existing.inserted_vertex;
        break;
      }
    }
    if (request.inserted_vertex.is_valid())
      continue;

    const int source_face = request.source_face.idx();
    if (source_face < 0 || static_cast<size_t>(source_face) >= descendants.size())
    {
      PhaseTimer::Exclusion exclusion("diagnostic");
      auto diagnostic = insertion_failure_context("source_face_out_of_range", request_index);
      diagnostic << " descendants_size=" << descendants.size();
      Logger::user_logger->error("{}", diagnostic.str());
      insertion_ok = false;
      break;
    }
    FaceHandle containing_face;
    double containing_weights[3] = {};
    std::array<VertexHandle, 3> containing_vertices;
    for (FaceHandle fh : descendants[source_face])
    {
      if (!fh.is_valid() || cage->status(fh).deleted())
        continue;
      size_t i = 0;
      for (VertexHandle vh : cage->fv_range(fh))
        containing_vertices[i++] = vh;
      double weights[3] = {};
      if (!detail::rail_barycentric_coordinates(
          request.point,
          cage->point(containing_vertices[0]),
          cage->point(containing_vertices[1]),
          cage->point(containing_vertices[2]), weights))
        continue;
      const double eps = 1e-8;
      if (weights[0] >= -eps && weights[1] >= -eps && weights[2] >= -eps)
      {
        containing_face = fh;
        containing_weights[0] = weights[0];
        containing_weights[1] = weights[1];
        containing_weights[2] = weights[2];
        break;
      }
    }
    if (!containing_face.is_valid())
    {
      PhaseTimer::Exclusion exclusion("diagnostic");
      auto diagnostic = insertion_failure_context("containing_face_not_found", request_index);
      diagnostic << " descendants_count=" << descendants[source_face].size()
                 << " bary_epsilon=" << 1e-8;
      for (FaceHandle fh : descendants[source_face])
      {
        diagnostic << " candidate_face=" << fh.idx();
        if (!fh.is_valid() || cage->status(fh).deleted())
        {
          diagnostic << " candidate_unavailable=1";
          continue;
        }
        std::array<Vec3d, 3> points;
        size_t count = 0;
        for (VertexHandle vh : cage->fv_range(fh))
        {
          if (count < points.size())
            points[count] = cage->point(vh);
          diagnostic << " candidate_vertex" << count << '=' << vh.idx();
          append_diagnostic_vector(diagnostic,
            "candidate_point" + std::to_string(count), cage->point(vh));
          count++;
        }
        if (count != points.size())
          continue;
        double weights[3] = {};
        const bool bary_valid = detail::rail_barycentric_coordinates(
          request.point, points[0], points[1], points[2], weights);
        diagnostic << " candidate_bary_valid=" << bary_valid;
        if (bary_valid)
          diagnostic << " candidate_weights=(" << weights[0] << ','
                     << weights[1] << ',' << weights[2] << ')';
      }
      Logger::user_logger->error("{}", diagnostic.str());
      insertion_ok = false;
      break;
    }

    size_t face_vertex_index = 0;
    for (VertexHandle vh : cage->fv_range(containing_face))
      containing_vertices[face_vertex_index++] = vh;
    const double bary_epsilon = 1e-8;
    int zero_weight = -1;
    for (size_t i = 0; i < 3; i++)
    {
      if (std::abs(containing_weights[i]) <= bary_epsilon)
      {
        zero_weight = static_cast<int>(i);
        break;
      }
    }
    const VertexHandle inserted = cage->add_vertex(request.point);
    if (zero_weight >= 0)
    {
      const VertexHandle a = containing_vertices[(zero_weight + 1) % 3];
      const VertexHandle b = containing_vertices[(zero_weight + 2) % 3];
      const EdgeHandle edge = edge_between(*cage, a, b);
      if (!edge.is_valid())
      {
        PhaseTimer::Exclusion exclusion("diagnostic");
        auto diagnostic = insertion_failure_context("containing_face_edge_not_found", request_index);
        diagnostic << " containing_face=" << containing_face.idx()
                   << " from_vertex=" << a.idx() << " to_vertex=" << b.idx();
        Logger::user_logger->error("{}", diagnostic.str());
        insertion_ok = false;
        break;
      }
      cage->split_copy(edge, inserted);
    }
    else
      cage->split_copy(containing_face, inserted);
    request.inserted_vertex = inserted;
    for (FaceHandle fh : cage->vf_range(inserted))
    {
      cage->property(source_face_property, fh) = source_face;
      if (std::find(descendants[source_face].begin(),
          descendants[source_face].end(), fh) == descendants[source_face].end())
        descendants[source_face].push_back(fh);
    }
  }
  cage->remove_property(source_face_property);

  if (!insertion_ok)
  {
    Logger::user_logger->error(
      "boundary rail construction could not insert every anchor into the cage surface.");
    return stats;
  }

  std::unordered_set<int> all_anchor_vertices;
  for (BoundaryLoop& loop : valid_loops)
  {
    for (size_t request_index : loop.request_indices)
    {
      const VertexHandle anchor = requests[request_index].inserted_vertex;
      if (!anchor.is_valid())
      {
        PhaseTimer::Exclusion exclusion("diagnostic");
        auto diagnostic = insertion_failure_context("missing_inserted_vertex", request_index);
        Logger::user_logger->error("{}", diagnostic.str());
        return stats;
      }
      loop.anchors.push_back(anchor);
      all_anchor_vertices.insert(anchor.idx());
    }
  }

  size_t built_rails = 0;
  std::vector<BoundaryRailLoop> rail_loops;
  for (BoundaryLoop& loop : valid_loops)
  {
    const int rail_id = static_cast<int>(built_rails);
    BoundaryRailLoop rail_loop;
    BoundaryRailSearchStats search_stats;
    const bool rail_built = detail::construct_boundary_rail(
        *cage, loop.anchors, all_anchor_vertices,
        rail_id, rail_loop, searchOptions, search_stats);
    Logger::user_logger->info(
      "rail_search mode={} component={} anchors={} candidates={} retry_cyclic={} time_limit_seconds={} attempts={} successful_shift={} states={} calls={} state_limit={} budget_hits={} timed_out={} outcome={} seconds={}",
      mode == BoundaryRailAnchorMode::EdgeMidpoint ? "edge-midpoint" : "vertex-bisector",
      loop.detected_loop_index, search_stats.anchors, searchOptions.candidateLimit,
      searchOptions.retryCyclicStarts, searchOptions.maxSeconds,
      search_stats.attempts, search_stats.successfulShift,
      search_stats.exploredStates, search_stats.recursiveCalls,
      std::max<size_t>(1024, search_stats.anchors * 16), search_stats.stateLimitHits,
      search_stats.timedOut, search_stats.outcome, search_stats.seconds);
    if (rail_built)
    {
      rail_loops.push_back(std::move(rail_loop));
      // Give the source boundary component the same id as its cage rail and
      // retain the per-edge outward co-normal.  Together with the boundary
      // edge tangent this defines the half-strip
      //   p(u, t) = edge(u) + t * outward, 0 <= u <= 1, t >= 0,
      // used as the Phase 3 rail support surface.
      for (size_t request_index : loop.request_indices)
      {
        const AnchorRequest& request = requests[request_index];
        const HalfedgeHandle source_halfedge =
          loop.halfedges[request.order_index];
        const EdgeHandle source_edge = source->edge_handle(source_halfedge);
        source->data(source_edge).boundary_rail_id = rail_id;
        source->data(source_edge).boundary_rail_outer_direction =
          request.support_outer_direction;
      }
      built_rails++;
    }
    else
    {
      Logger::user_logger->warn(
        "boundary rail construction ({}) could not form a simple closed cage edge loop for source boundary component {}.",
        mode == BoundaryRailAnchorMode::EdgeMidpoint ?
          "edge-midpoint" : "vertex-bisector",
        loop.detected_loop_index);
    }
  }

  // Keep all rails in one network so shortening respects the other loops.
  // The embedding routine commits only after checking the split mesh and
  // ordered closed rails; failure leaves these valid Dijkstra labels intact.
  if (!rail_loops.empty())
  {
    if (refine_boundary_rails_with_flip_geodesics(*cage, rail_loops, stats.geodesic))
    {
      Logger::user_logger->info(
        "boundary rail flip geodesics: refined {} loops, {} intrinsic flips, {} shorten iterations, length {} -> {}, inserted {} cage vertices.",
        rail_loops.size(), stats.geodesic.intrinsicFlipCount,
        stats.geodesic.shortenIterationCount, stats.geodesic.lengthBefore,
        stats.geodesic.lengthAfter, stats.geodesic.insertedVertexCount);
    }
    else
    {
      Logger::user_logger->warn(
        "boundary rail flip geodesics was not applied; keeping the Dijkstra rails and cage unchanged: {}.",
        stats.geodesic.failureReason);
    }
  }

  size_t rail_vertices = 0;
  size_t rail_edges = 0;
  for (VertexHandle vh : cage->vertices())
    rail_vertices += cage->data(vh).boundary_rail_id != kNoRail;
  for (EdgeHandle eh : cage->edges())
    rail_edges += cage->data(eh).boundary_rail_id != kNoRail;
  stats.builtLoopCount = built_rails;
  stats.railVertexCount = rail_vertices;
  stats.railEdgeCount = rail_edges;
  // The minimum triangle quality was only reported; disabled because the
  // rails do not need it.
  // double min_quality = DBL_MAX;
  // for (FaceHandle fh : cage->faces())
  // {
  //   std::array<Vec3d, 3> points;
  //   size_t i = 0;
  //   for (VertexHandle vh : cage->fv_range(fh))
  //     points[i++] = cage->point(vh);
  //   const double a = (points[1] - points[0]).length();
  //   const double b = (points[2] - points[1]).length();
  //   const double c = (points[0] - points[2]).length();
  //   const double denominator = a * a + b * b + c * c;
  //   const double quality = denominator > 0.0 ?
  //     2.0 * std::sqrt(3.0) *
  //     (points[1] - points[0]).cross(points[2] - points[0]).length() /
  //     denominator : 0.0;
  //   min_quality = std::min(min_quality, quality);
  // }
  Logger::user_logger->info(
    "boundary rail construction ({}) built {} closed rails from {} detected loops ({} ray-valid), with {} anchor requests, {} rail vertices, and {} rail edges.",
    mode == BoundaryRailAnchorMode::EdgeMidpoint ?
      "edge-midpoint" : "vertex-bisector",
    built_rails, stats.detectedLoopCount, stats.rayValidLoopCount,
    requests.size(), rail_vertices, rail_edges);
  // Match the pre-benchmark return semantics: ray-invalid source loops are
  // reported but do not make a partially successful build fail.
  stats.successful = built_rails == valid_loops.size();
  return stats;
}

}// namespace CageInit
}// namespace Cage
