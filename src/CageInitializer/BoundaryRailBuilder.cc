#include "BoundaryRailBuilder.hh"

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <queue>
#include <set>
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
  std::vector<HalfedgeHandle> halfedges;
  std::vector<size_t> request_indices;
  std::vector<VertexHandle> anchors;
};

bool get_boundary_edge_outer_direction(
  SMeshT& source, HalfedgeHandle boundary_halfedge,
  double length_epsilon, Vec3d& outer)
{
  if (!source.is_boundary(boundary_halfedge))
    return false;
  const HalfedgeHandle inner = source.opposite_halfedge_handle(boundary_halfedge);
  const FaceHandle face = source.face_handle(inner);
  if (!face.is_valid())
    return false;

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
    return false;

  const Vec3d p0 = source.point(v0);
  const Vec3d p1 = source.point(v1);
  const Vec3d edge = p1 - p0;
  const double edge_sqr = edge.squaredNorm();
  if (edge_sqr <= length_epsilon * length_epsilon)
    return false;

  const Vec3d midpoint = (p0 + p1) * 0.5;
  const Vec3d toward_opposite = source.point(opposite) - midpoint;
  Vec3d inward = toward_opposite - edge * ((toward_opposite | edge) / edge_sqr);
  if (inward.length() <= length_epsilon)
    return false;
  inward.normalize();
  outer = -inward;
  return true;
}

bool barycentric_coordinates(
  const Vec3d& p, const Vec3d& a, const Vec3d& b, const Vec3d& c,
  double weights[3])
{
  const Vec3d v0 = b - a;
  const Vec3d v1 = c - a;
  const Vec3d v2 = p - a;
  const double d00 = v0 | v0;
  const double d01 = v0 | v1;
  const double d11 = v1 | v1;
  const double d20 = v2 | v0;
  const double d21 = v2 | v1;
  const double denom = d00 * d11 - d01 * d01;
  if (std::abs(denom) <= 1e-30 * std::max(d00 * d11, 1.0))
    return false;
  weights[1] = (d11 * d20 - d01 * d21) / denom;
  weights[2] = (d00 * d21 - d01 * d20) / denom;
  weights[0] = 1.0 - weights[1] - weights[2];
  return true;
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

bool shortest_path(
  SMeshT& mesh, VertexHandle source, VertexHandle target,
  const std::unordered_set<int>& blocked_vertices,
  const std::unordered_set<int>& blocked_edges,
  GraphPath& path)
{
  const size_t vertex_capacity = mesh.n_vertices();
  std::vector<double> distance(vertex_capacity, DBL_MAX);
  std::vector<int> previous(vertex_capacity, -1);
  using QueueEntry = std::pair<double, int>;
  std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;
  distance[source.idx()] = 0.0;
  queue.emplace(0.0, source.idx());

  while (!queue.empty())
  {
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
  size_t max_candidates)
{
  std::vector<GraphPath> result;
  GraphPath shortest;
  if (!shortest_path(
    mesh, source, target, blocked_vertices, blocked_edges, shortest))
    return result;
  result.push_back(shortest);

  const size_t edge_count = shortest.vertices.size() > 0 ?
    shortest.vertices.size() - 1 : 0;
  const size_t probes = std::min<size_t>(edge_count, 32);
  for (size_t probe = 0; probe < probes && result.size() < max_candidates; probe++)
  {
    const size_t i = probes == edge_count ? probe :
      (probe * edge_count) / probes;
    const EdgeHandle excluded = edge_between(
      mesh, shortest.vertices[i], shortest.vertices[i + 1]);
    if (!excluded.is_valid())
      continue;
    std::unordered_set<int> trial_blocked_edges = blocked_edges;
    trial_blocked_edges.insert(excluded.idx());
    GraphPath alternative;
    if (!shortest_path(
      mesh, source, target, blocked_vertices,
      trial_blocked_edges, alternative))
      continue;
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

bool construct_rail(
  SMeshT& mesh, const std::vector<VertexHandle>& input_anchors,
  const std::unordered_set<int>& all_anchor_vertices,
  int rail_id)
{
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
    return false;

  std::unordered_set<int> used_vertices;
  std::unordered_set<int> used_edges;
  std::vector<GraphPath> selected(anchors.size());
  size_t search_nodes = 0;
  const size_t max_search_nodes = std::max<size_t>(1024, anchors.size() * 16);

  std::function<bool(size_t)> select_paths = [&](size_t segment)
  {
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
      mesh, source, target, blocked_vertices, used_edges, 8);
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
    }
    return false;
  };

  if (!select_paths(0))
    return false;

  std::vector<VertexHandle> cycle_vertices;
  std::vector<EdgeHandle> cycle_edges;
  if (!validate_cycle(mesh, selected, anchors, cycle_vertices, cycle_edges))
    return false;

  for (VertexHandle vh : cycle_vertices)
    mesh.data(vh).boundary_rail_id = rail_id;
  for (EdgeHandle eh : cycle_edges)
    mesh.data(eh).boundary_rail_id = rail_id;
  return true;
}
}// namespace

bool BoundaryRailBuilder::build()
{
  if (!source || !cage || source->n_vertices() == 0 || cage->n_faces() == 0)
    return false;

  Box source_box;
  for (VertexHandle vh : source->vertices())
    source_box.add(source->point(vh));
  Box cage_box;
  for (VertexHandle vh : cage->vertices())
    cage_box.add(cage->point(vh));
  const double source_scale = std::max((source_box.hi - source_box.lo).length(), 1.0);
  const double cage_scale = std::max((cage_box.hi - cage_box.lo).length(), 1.0);
  const double length_epsilon = source_scale * 1e-12;
  const double intersection_epsilon = cage_scale * 1e-10;
  const double anchor_merge_epsilon = cage_scale * 1e-9;

  std::vector<BoundaryLoop> loops;
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
    if (valid && current == start && loop.halfedges.size() >= 3)
      loops.push_back(std::move(loop));
  }

  if (loops.empty())
  {
    Logger::user_logger->info(
      "boundary rail construction: source mesh has no closed boundary loop.");
    return true;
  }

  RayTriangleBvh ray_tree(*cage);
  std::vector<AnchorRequest> requests;
  std::vector<BoundaryLoop> valid_loops;
  valid_loops.reserve(loops.size());
  for (size_t loop_index = 0; loop_index < loops.size(); loop_index++)
  {
    BoundaryLoop candidate = loops[loop_index];
    std::vector<AnchorRequest> loop_requests;
    loop_requests.reserve(candidate.halfedges.size());
    bool valid = true;
    for (size_t i = 0; i < candidate.halfedges.size(); i++)
    {
      const HalfedgeHandle incoming =
        candidate.halfedges[(i + candidate.halfedges.size() - 1) % candidate.halfedges.size()];
      const HalfedgeHandle outgoing = candidate.halfedges[i];
      Vec3d incoming_outer;
      Vec3d outgoing_outer;
      if (!get_boundary_edge_outer_direction(
          *source, incoming, length_epsilon, incoming_outer) ||
        !get_boundary_edge_outer_direction(
          *source, outgoing, length_epsilon, outgoing_outer))
      {
        valid = false;
        break;
      }

      Vec3d direction = incoming_outer + outgoing_outer;
      // A near-zero average direction is deliberately not given a heuristic
      // fallback in this first implementation.
      if (direction.length() <= length_epsilon)
      {
        valid = false;
        break;
      }
      direction.normalize();
      const VertexHandle source_vertex = source->from_vertex_handle(outgoing);
      const Vec3d origin = source->point(source_vertex);
      RayHit hit;
      if (!ray_tree.first_intersection(
          origin, direction, intersection_epsilon, hit) || !finite_vec(hit.point))
      {
        valid = false;
        break;
      }

      const RayTriangle& triangle = ray_tree.triangle(hit.triangle);
      AnchorRequest request;
      request.loop_index = valid_loops.size();
      request.order_index = i;
      request.point = hit.point;
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
      }
      else
        request.location = AnchorRequest::Location::Face;
      loop_requests.push_back(request);
    }

    if (!valid)
    {
      Logger::user_logger->warn(
        "boundary rail construction skipped source boundary component {} because an anchor ray was undefined or missed the cage.",
        loop_index);
      continue;
    }
    candidate.request_indices.clear();
    for (AnchorRequest& request : loop_requests)
    {
      candidate.request_indices.push_back(requests.size());
      requests.push_back(std::move(request));
    }
    valid_loops.push_back(std::move(candidate));
  }

  if (valid_loops.empty())
    return false;

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
        insertion_ok = false;
        break;
      }
      const VertexHandle inserted = cage->add_vertex(request.point);
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
      if (!barycentric_coordinates(
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
    return false;
  }

  std::unordered_set<int> all_anchor_vertices;
  for (BoundaryLoop& loop : valid_loops)
  {
    for (size_t request_index : loop.request_indices)
    {
      const VertexHandle anchor = requests[request_index].inserted_vertex;
      if (!anchor.is_valid())
        return false;
      loop.anchors.push_back(anchor);
      all_anchor_vertices.insert(anchor.idx());
    }
  }

  size_t built_rails = 0;
  for (BoundaryLoop& loop : valid_loops)
  {
    if (construct_rail(
        *cage, loop.anchors, all_anchor_vertices,
        static_cast<int>(built_rails)))
    {
      built_rails++;
    }
    else
    {
      Logger::user_logger->warn(
        "boundary rail construction could not form a simple closed cage edge loop for source boundary component {}.",
        built_rails);
    }
  }

  size_t rail_vertices = 0;
  size_t rail_edges = 0;
  for (VertexHandle vh : cage->vertices())
    rail_vertices += cage->data(vh).boundary_rail_id != kNoRail;
  for (EdgeHandle eh : cage->edges())
    rail_edges += cage->data(eh).boundary_rail_id != kNoRail;
  double min_quality = DBL_MAX;
  for (FaceHandle fh : cage->faces())
  {
    std::array<Vec3d, 3> points;
    size_t i = 0;
    for (VertexHandle vh : cage->fv_range(fh))
      points[i++] = cage->point(vh);
    const double a = (points[1] - points[0]).length();
    const double b = (points[2] - points[1]).length();
    const double c = (points[0] - points[2]).length();
    const double denominator = a * a + b * b + c * c;
    const double quality = denominator > 0.0 ?
      2.0 * std::sqrt(3.0) *
      (points[1] - points[0]).cross(points[2] - points[0]).length() /
      denominator : 0.0;
    min_quality = std::min(min_quality, quality);
  }
  Logger::user_logger->info(
    "boundary rail construction built {} closed rails with {} anchors, {} rail vertices, and {} rail edges; cage min triangle quality {}.",
    built_rails, requests.size(), rail_vertices, rail_edges,
    min_quality == DBL_MAX ? 0.0 : min_quality);
  return built_rails == valid_loops.size();
}

}// namespace CageInit
}// namespace Cage
