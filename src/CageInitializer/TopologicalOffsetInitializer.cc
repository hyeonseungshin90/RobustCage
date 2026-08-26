#include "TopologicalOffsetInitializer.hh"

#include <algorithm>
#include <array>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "Geometry/Exact/CGALTypes.h"
#include "Utils/logger.hh"

namespace Cage
{
namespace CageInit
{
using namespace Geometry;
using namespace SimpleUtils;

namespace
{
using VertexId = size_t;

struct EdgeKey
{
  VertexId a;
  VertexId b;

  EdgeKey() : a(0), b(0) {}
  EdgeKey(VertexId v0, VertexId v1)
    : a(std::min(v0, v1)), b(std::max(v0, v1))
  {}

  bool operator<(const EdgeKey& rhs) const
  {
    return std::tie(a, b) < std::tie(rhs.a, rhs.b);
  }
};

struct FaceKey
{
  std::array<VertexId, 3> vertices;

  FaceKey() : vertices{ 0, 0, 0 } {}
  FaceKey(VertexId v0, VertexId v1, VertexId v2)
    : vertices{ v0, v1, v2 }
  {
    std::sort(vertices.begin(), vertices.end());
  }

  bool operator<(const FaceKey& rhs) const
  {
    return vertices < rhs.vertices;
  }
};

struct EmbeddedVertex
{
  ExactPoint point;
  VM::VVertexProp prop;
  bool is_offset = false;
};

struct EmbeddedTet
{
  std::array<VertexId, 4> vertices;
  VM::VCellProp prop;
};

bool contains(const EmbeddedTet& tet, VertexId vertex)
{
  return std::find(tet.vertices.begin(), tet.vertices.end(), vertex) !=
    tet.vertices.end();
}

bool contains(const EmbeddedTet& tet, const EdgeKey& edge)
{
  return contains(tet, edge.a) && contains(tet, edge.b);
}

bool contains(const EmbeddedTet& tet, const FaceKey& face)
{
  return contains(tet, face.vertices[0]) &&
    contains(tet, face.vertices[1]) &&
    contains(tet, face.vertices[2]);
}

VertexId oppositeVertex(const EmbeddedTet& tet, const FaceKey& face)
{
  for (VertexId vertex : tet.vertices)
  {
    if (std::find(face.vertices.begin(), face.vertices.end(), vertex) ==
      face.vertices.end())
      return vertex;
  }
  throw std::logic_error("face is not a face of its tetrahedron");
}

std::array<VertexId, 2> oppositeVertices(
  const EmbeddedTet& tet, const EdgeKey& edge)
{
  std::array<VertexId, 2> result;
  size_t count = 0;
  for (VertexId vertex : tet.vertices)
  {
    if (vertex != edge.a && vertex != edge.b)
      result[count++] = vertex;
  }
  if (count != 2)
    throw std::logic_error("edge is not an edge of its tetrahedron");
  return result;
}

std::array<FaceKey, 4> tetFaces(const EmbeddedTet& tet)
{
  const auto& v = tet.vertices;
  return {
    FaceKey(v[1], v[2], v[3]),
    FaceKey(v[0], v[2], v[3]),
    FaceKey(v[0], v[1], v[3]),
    FaceKey(v[0], v[1], v[2])
  };
}

std::array<EdgeKey, 6> tetEdges(const EmbeddedTet& tet)
{
  const auto& v = tet.vertices;
  return {
    EdgeKey(v[0], v[1]), EdgeKey(v[0], v[2]), EdgeKey(v[0], v[3]),
    EdgeKey(v[1], v[2]), EdgeKey(v[1], v[3]), EdgeKey(v[2], v[3])
  };
}

std::array<EdgeKey, 3> faceEdges(const FaceKey& face)
{
  const auto& v = face.vertices;
  return {
    EdgeKey(v[0], v[1]), EdgeKey(v[1], v[2]), EdgeKey(v[0], v[2])
  };
}

ExactPoint averagePoint(
  std::vector<EmbeddedVertex>& vertices,
  const std::vector<VertexId>& ids)
{
  ASSERT(!ids.empty(), "cannot average an empty vertex set");
  Vector_3 sum(ET(0), ET(0), ET(0));
  for (VertexId id : ids)
    sum = sum + vertices[id].point.exactVec3();
  return ExactPoint(sum * (ET(1) / ET(ids.size())));
}

VertexId appendVertex(
  std::vector<EmbeddedVertex>& vertices,
  const ExactPoint& point,
  bool is_offset = false)
{
  EmbeddedVertex vertex;
  vertex.point = point;
  vertex.prop = VM::VVertexProp();
  vertex.is_offset = is_offset;
  vertices.push_back(std::move(vertex));
  return vertices.size() - 1;
}

void splitTetAtCenter(
  const EmbeddedTet& tet,
  VertexId center,
  std::vector<EmbeddedTet>& output)
{
  const auto& v = tet.vertices;
  output.push_back({ { center, v[1], v[2], v[3] }, tet.prop });
  output.push_back({ { v[0], center, v[2], v[3] }, tet.prop });
  output.push_back({ { v[0], v[1], center, v[3] }, tet.prop });
  output.push_back({ { v[0], v[1], v[2], center }, tet.prop });
}

void splitTetsOnFace(
  std::vector<EmbeddedTet>& tets,
  const FaceKey& face,
  VertexId center)
{
  std::vector<EmbeddedTet> result;
  result.reserve(tets.size() + 2);
  for (const EmbeddedTet& tet : tets)
  {
    if (!contains(tet, face))
    {
      result.push_back(tet);
      continue;
    }

    const VertexId opposite = oppositeVertex(tet, face);
    const auto& v = face.vertices;
    result.push_back({ { center, v[0], v[1], opposite }, tet.prop });
    result.push_back({ { center, v[1], v[2], opposite }, tet.prop });
    result.push_back({ { center, v[2], v[0], opposite }, tet.prop });
  }
  tets = std::move(result);
}

void splitTetsOnEdge(
  std::vector<EmbeddedTet>& tets,
  const EdgeKey& edge,
  VertexId split_vertex)
{
  std::vector<EmbeddedTet> result;
  result.reserve(tets.size() + 1);
  for (const EmbeddedTet& tet : tets)
  {
    if (!contains(tet, edge))
    {
      result.push_back(tet);
      continue;
    }

    const auto opposite = oppositeVertices(tet, edge);
    result.push_back(
      { { edge.a, split_vertex, opposite[0], opposite[1] }, tet.prop });
    result.push_back(
      { { split_vertex, edge.b, opposite[0], opposite[1] }, tet.prop });
  }
  tets = std::move(result);
}

void splitTaggedBoundaryTets(
  std::vector<EmbeddedVertex>& vertices,
  std::vector<EmbeddedTet>& tets,
  const std::set<FaceKey>& constraint_faces,
  size_t& split_count)
{
  std::vector<EmbeddedTet> result;
  result.reserve(tets.size());
  for (const EmbeddedTet& tet : tets)
  {
    const auto faces = tetFaces(tet);
    const bool all_faces_constraint = std::all_of(
      faces.begin(), faces.end(),
      [&](const FaceKey& face) { return constraint_faces.count(face) != 0; });
    if (!all_faces_constraint)
    {
      result.push_back(tet);
      continue;
    }

    const VertexId center = appendVertex(
      vertices,
      averagePoint(vertices, {
        tet.vertices[0], tet.vertices[1], tet.vertices[2], tet.vertices[3] }));
    splitTetAtCenter(tet, center, result);
    split_count++;
  }
  tets = std::move(result);
}

void splitTaggedBoundaryFaces(
  std::vector<EmbeddedVertex>& vertices,
  std::vector<EmbeddedTet>& tets,
  const std::set<EdgeKey>& constraint_edges,
  const std::set<FaceKey>& constraint_faces,
  size_t& split_count)
{
  std::set<FaceKey> candidates;
  for (const EmbeddedTet& tet : tets)
  {
    for (const FaceKey& face : tetFaces(tet))
    {
      if (constraint_faces.count(face) != 0)
        continue;
      const auto edges = faceEdges(face);
      if (std::all_of(
        edges.begin(), edges.end(),
        [&](const EdgeKey& edge) { return constraint_edges.count(edge) != 0; }))
        candidates.insert(face);
    }
  }

  std::map<FaceKey, VertexId> centers;
  for (const FaceKey& face : candidates)
  {
    centers.emplace(
      face,
      appendVertex(
        vertices,
        averagePoint(vertices, {
          face.vertices[0], face.vertices[1], face.vertices[2] })));
  }

  std::vector<EmbeddedTet> result;
  for (const EmbeddedTet& tet : tets)
  {
    std::vector<EmbeddedTet> local{ tet };
    std::vector<FaceKey> local_candidates;
    for (const FaceKey& face : tetFaces(tet))
    {
      if (candidates.count(face) != 0)
        local_candidates.push_back(face);
    }
    std::sort(local_candidates.begin(), local_candidates.end());
    for (const FaceKey& face : local_candidates)
      splitTetsOnFace(local, face, centers.at(face));
    result.insert(result.end(), local.begin(), local.end());
  }

  split_count += candidates.size();
  tets = std::move(result);
}

void splitTaggedBoundaryEdges(
  std::vector<EmbeddedVertex>& vertices,
  std::vector<EmbeddedTet>& tets,
  const std::set<EdgeKey>& constraint_edges,
  size_t& split_count)
{
  std::set<EdgeKey> candidates;
  for (const EmbeddedTet& tet : tets)
  {
    for (const EdgeKey& edge : tetEdges(tet))
    {
      if (constraint_edges.count(edge) == 0 &&
        vertices[edge.a].prop.is_constraint &&
        vertices[edge.b].prop.is_constraint)
        candidates.insert(edge);
    }
  }

  std::map<EdgeKey, VertexId> split_vertices;
  for (const EdgeKey& edge : candidates)
  {
    split_vertices.emplace(
      edge,
      appendVertex(
        vertices,
        averagePoint(vertices, { edge.a, edge.b })));
  }

  std::vector<EmbeddedTet> result;
  for (const EmbeddedTet& tet : tets)
  {
    std::vector<EmbeddedTet> local{ tet };
    std::vector<EdgeKey> local_candidates;
    for (const EdgeKey& edge : tetEdges(tet))
    {
      if (candidates.count(edge) != 0)
        local_candidates.push_back(edge);
    }
    std::sort(local_candidates.begin(), local_candidates.end());
    for (const EdgeKey& edge : local_candidates)
      splitTetsOnEdge(local, edge, split_vertices.at(edge));
    result.insert(result.end(), local.begin(), local.end());
  }

  split_count += candidates.size();
  tets = std::move(result);
}

void validateSimplicialEmbedding(
  const std::vector<EmbeddedVertex>& vertices,
  const std::vector<EmbeddedTet>& tets,
  const std::set<EdgeKey>& constraint_edges,
  const std::set<FaceKey>& constraint_faces)
{
  for (size_t tet_index = 0; tet_index < tets.size(); tet_index++)
  {
    const EmbeddedTet& tet = tets[tet_index];
    std::vector<VertexId> input_vertices;
    for (VertexId vertex : tet.vertices)
    {
      if (vertices[vertex].prop.is_constraint)
        input_vertices.push_back(vertex);
    }

    bool valid = input_vertices.size() <= 1;
    if (input_vertices.size() == 2)
      valid = constraint_edges.count(
        EdgeKey(input_vertices[0], input_vertices[1])) != 0;
    else if (input_vertices.size() == 3)
      valid = constraint_faces.count(FaceKey(
        input_vertices[0], input_vertices[1], input_vertices[2])) != 0;
    else if (input_vertices.size() == 4)
      valid = false;

    if (!valid)
    {
      Logger::user_logger->critical(
        "simplicial embedding validation failed at tetrahedron {} with {} input vertices.",
        tet_index, input_vertices.size());
      throw std::logic_error("invalid simplicial embedding");
    }
  }
}

void insertOffset(
  std::vector<EmbeddedVertex>& vertices,
  std::vector<EmbeddedTet>& tets,
  size_t& split_count)
{
  // TopologicalOffsets' current implementation uses scalar values 0 on the
  // input and 1 elsewhere and marches at isovalue 0.25.  Therefore every
  // inserted point lies one quarter of the mixed edge length away from S.
  const ET input_weight = ET(3) / ET(4);
  const ET outside_weight = ET(1) / ET(4);

  std::set<EdgeKey> mixed_edges;
  for (const EmbeddedTet& tet : tets)
  {
    for (const EdgeKey& edge : tetEdges(tet))
    {
      if (vertices[edge.a].prop.is_constraint !=
        vertices[edge.b].prop.is_constraint)
        mixed_edges.insert(edge);
    }
  }

  std::map<EdgeKey, VertexId> offset_vertices;
  for (const EdgeKey& edge : mixed_edges)
  {
    const VertexId input = vertices[edge.a].prop.is_constraint ? edge.a : edge.b;
    const VertexId outside = input == edge.a ? edge.b : edge.a;
    const Vector_3 point =
      vertices[input].point.exactVec3() * input_weight +
      vertices[outside].point.exactVec3() * outside_weight;
    offset_vertices.emplace(
      edge, appendVertex(vertices, ExactPoint(point), true));
  }

  std::vector<EmbeddedTet> result;
  for (const EmbeddedTet& tet : tets)
  {
    std::vector<EmbeddedTet> local{ tet };
    std::vector<EdgeKey> local_edges;
    for (const EdgeKey& edge : tetEdges(tet))
    {
      if (mixed_edges.count(edge) != 0)
        local_edges.push_back(edge);
    }
    std::sort(local_edges.begin(), local_edges.end());
    for (const EdgeKey& edge : local_edges)
      splitTetsOnEdge(local, edge, offset_vertices.at(edge));
    result.insert(result.end(), local.begin(), local.end());
  }

  split_count = mixed_edges.size();
  tets = std::move(result);

  for (EmbeddedTet& tet : tets)
  {
    size_t input_count = 0;
    size_t offset_count = 0;
    size_t other_count = 0;
    for (VertexId vertex : tet.vertices)
    {
      input_count += vertices[vertex].prop.is_constraint;
      offset_count += vertices[vertex].is_offset;
      other_count += !vertices[vertex].prop.is_constraint &&
        !vertices[vertex].is_offset;
    }
    tet.prop.is_adjacent_constraint =
      input_count > 0 && offset_count > 0 && other_count == 0;
  }
}

VM::HalfEdgeHandle orientedHalfEdge(
  VM::VMeshT& mesh,
  VM::EdgeHandle edge,
  VM::VertexHandle from,
  VM::VertexHandle to)
{
  if (mesh.edge(edge).from() == from && mesh.edge(edge).to() == to)
    return halfEdgeHdl(edge, VM::Sequence);
  ASSERT(
    mesh.edge(edge).from() == to && mesh.edge(edge).to() == from,
    "edge endpoints do not match requested halfedge");
  return halfEdgeHdl(edge, VM::Reverse);
}

void rebuildVolumeMesh(
  VM::VMeshT& output,
  const std::vector<EmbeddedVertex>& vertices,
  const std::vector<EmbeddedTet>& tets,
  const std::set<EdgeKey>& constraint_edges,
  const std::set<FaceKey>& constraint_faces)
{
  VM::VMeshT rebuilt;
  rebuilt.vertices.reserve(vertices.size());
  for (const EmbeddedVertex& vertex : vertices)
  {
    const VM::VertexHandle vh = rebuilt.addVertex(vertex.point);
    rebuilt.vertex(vh).prop = vertex.prop;
    rebuilt.vertex(vh).prop.is_on_boundary = false;
  }

  std::map<EdgeKey, VM::EdgeHandle> edge_handles;
  auto ensure_edge = [&](const EdgeKey& key) -> VM::EdgeHandle
  {
    auto found = edge_handles.find(key);
    if (found != edge_handles.end())
      return found->second;
    const VM::EdgeHandle edge = rebuilt.addEdge(
      VM::VertexHandle(key.a), VM::VertexHandle(key.b));
    rebuilt.edge(edge).prop.is_constraint = constraint_edges.count(key) != 0;
    rebuilt.edge(edge).prop.is_on_boundary = false;
    edge_handles.emplace(key, edge);
    return edge;
  };

  std::map<FaceKey, VM::FaceHandle> face_handles;
  auto ensure_face = [&](const FaceKey& key) -> VM::FaceHandle
  {
    auto found = face_handles.find(key);
    if (found != face_handles.end())
      return found->second;

    const VM::VertexHandle v0(key.vertices[0]);
    const VM::VertexHandle v1(key.vertices[1]);
    const VM::VertexHandle v2(key.vertices[2]);
    const VM::EdgeHandle e01 = ensure_edge(EdgeKey(key.vertices[0], key.vertices[1]));
    const VM::EdgeHandle e12 = ensure_edge(EdgeKey(key.vertices[1], key.vertices[2]));
    const VM::EdgeHandle e20 = ensure_edge(EdgeKey(key.vertices[2], key.vertices[0]));
    VM::SubHalfEdges halfedges = {
      orientedHalfEdge(rebuilt, e01, v0, v1),
      orientedHalfEdge(rebuilt, e12, v1, v2),
      orientedHalfEdge(rebuilt, e20, v2, v0)
    };
    const VM::FaceHandle face = rebuilt.addFace(halfedges, false);
    rebuilt.face(face).prop.is_constraint = constraint_faces.count(key) != 0;
    rebuilt.face(face).prop.is_boundary = false;
    face_handles.emplace(key, face);
    return face;
  };

  std::vector<std::array<VM::VertexHandle, 4>> cell_vertices;
  cell_vertices.reserve(tets.size());
  for (size_t tet_index = 0; tet_index < tets.size(); tet_index++)
  {
    const EmbeddedTet& tet = tets[tet_index];
    const Point_3& p0 = vertices[tet.vertices[0]].point.exact();
    const Point_3& p1 = vertices[tet.vertices[1]].point.exact();
    const Point_3& p2 = vertices[tet.vertices[2]].point.exact();
    const Point_3& p3 = vertices[tet.vertices[3]].point.exact();
    if (CGAL::orientation(p0, p1, p2, p3) == CGAL::COPLANAR)
    {
      Logger::user_logger->critical(
        "topological offset construction produced degenerate tetrahedron {}.",
        tet_index);
      throw std::logic_error("degenerate tetrahedron during topological offset insertion");
    }

    const auto faces = tetFaces(tet);
    VM::SubHalfFaces half_faces;
    for (size_t i = 0; i < faces.size(); i++)
    {
      const FaceKey& key = faces[i];
      const VM::FaceHandle face = ensure_face(key);
      const VertexId opposite = oppositeVertex(tet, key);
      const CGAL::Orientation orientation = CGAL::orientation(
        vertices[key.vertices[0]].point.exact(),
        vertices[key.vertices[1]].point.exact(),
        vertices[key.vertices[2]].point.exact(),
        vertices[opposite].point.exact());
      ASSERT(orientation != CGAL::COPLANAR, "degenerate tetrahedron face");
      half_faces[i] = halfFaceHdl(
        face,
        orientation == CGAL::POSITIVE ? VM::Sequence : VM::Reverse);
    }

    const VM::CellHandle cell = rebuilt.addCell(half_faces);
    rebuilt.cell(cell).prop = tet.prop;
    cell_vertices.push_back({
      VM::VertexHandle(tet.vertices[0]), VM::VertexHandle(tet.vertices[1]),
      VM::VertexHandle(tet.vertices[2]), VM::VertexHandle(tet.vertices[3]) });
  }

  for (size_t face_index = 0; face_index < rebuilt.nFaces(); face_index++)
  {
    const VM::FaceHandle face(face_index);
    const size_t adjacent_cells = rebuilt.face(face).nConnCells();
    ASSERT(
      adjacent_cells == 1 || adjacent_cells == 2,
      "rebuilt volume face must have one or two incident tetrahedra");
    if (adjacent_cells != 1)
      continue;

    rebuilt.face(face).prop.is_boundary = true;
    for (VM::EdgeHandle edge : rebuilt.face(face).edges())
    {
      rebuilt.edge(edge).prop.is_on_boundary = true;
      rebuilt.vertex(rebuilt.edge(edge).from()).prop.is_on_boundary = true;
      rebuilt.vertex(rebuilt.edge(edge).to()).prop.is_on_boundary = true;
    }
  }

  rebuilt.updateVertexConnCells(cell_vertices);

  // Every non-input face separating selected and unselected cells must be an
  // offset face.  This is the boundary that retrieveCage() will extract.
  size_t offset_face_count = 0;
  for (size_t face_index = 0; face_index < rebuilt.nFaces(); face_index++)
  {
    const VM::FaceHandle face(face_index);
    if (rebuilt.face(face).prop.is_constraint)
      continue;

    size_t selected_count = 0;
    for (VM::CellHandle cell : rebuilt.face(face).connCells())
      selected_count += rebuilt.cell(cell).prop.is_adjacent_constraint;
    if (selected_count != 1)
      continue;

    for (VM::VertexHandle vertex : rebuilt.findFV(face))
    {
      if (!vertices[vertex.idx()].is_offset)
      {
        Logger::user_logger->critical(
          "offset insertion produced a non-offset vertex on extracted face {}.",
          face_index);
        throw std::logic_error("invalid offset boundary face");
      }
    }
    offset_face_count++;
  }

  if (offset_face_count == 0)
    throw std::logic_error("topological offset insertion produced no boundary faces");

  output = std::move(rebuilt);
  Logger::user_logger->info(
    "topological offset volume: {} vertices, {} edges, {} faces, {} tetrahedra; {} extractable offset faces.",
    output.nVertices(), output.nEdges(), output.nFaces(), output.nCells(),
    offset_face_count);
}

}// namespace

TopologicalOffsetInitializer::TopologicalOffsetInitializer(VM::VMeshT* _mesh)
  : mesh(_mesh)
{}

void TopologicalOffsetInitializer::generate()
{
  ASSERT(mesh != nullptr, "topological offset initializer requires a volume mesh");

  std::vector<EmbeddedVertex> vertices;
  vertices.reserve(mesh->nVertices());
  for (size_t vertex_index = 0; vertex_index < mesh->nVertices(); vertex_index++)
  {
    const VM::VertexHandle vertex(vertex_index);
    EmbeddedVertex embedded_vertex;
    embedded_vertex.point = mesh->exact_point(vertex);
    embedded_vertex.prop = mesh->vertex(vertex).prop;
    vertices.push_back(std::move(embedded_vertex));
  }

  std::set<EdgeKey> constraint_edges;
  for (size_t edge_index = 0; edge_index < mesh->nEdges(); edge_index++)
  {
    const VM::EdgeHandle edge(edge_index);
    if (mesh->edge(edge).prop.is_constraint)
      constraint_edges.emplace(
        mesh->edge(edge).from().idx(), mesh->edge(edge).to().idx());
  }

  std::set<FaceKey> constraint_faces;
  for (size_t face_index = 0; face_index < mesh->nFaces(); face_index++)
  {
    const VM::FaceHandle face(face_index);
    if (!mesh->face(face).prop.is_constraint)
      continue;
    const auto face_vertices = mesh->findFV(face);
    constraint_faces.emplace(
      face_vertices[0].idx(), face_vertices[1].idx(), face_vertices[2].idx());
  }

  std::vector<EmbeddedTet> tets;
  tets.reserve(mesh->nCells());
  for (size_t cell_index = 0; cell_index < mesh->nCells(); cell_index++)
  {
    const VM::CellHandle cell(cell_index);
    if (mesh->deleted(cell))
      continue;
    const auto cell_vertices = mesh->findCV(cell);
    tets.push_back({ {
      size_t(cell_vertices[0].idx()), size_t(cell_vertices[1].idx()),
      size_t(cell_vertices[2].idx()), size_t(cell_vertices[3].idx()) },
      mesh->cell(cell).prop });
  }

  Logger::user_logger->info(
    "begin simplicial embedding: {} vertices, {} tetrahedra, {} constraint edges, {} constraint faces.",
    vertices.size(), tets.size(), constraint_edges.size(), constraint_faces.size());

  size_t split_tets = 0;
  size_t split_faces = 0;
  size_t split_edges = 0;
  splitTaggedBoundaryTets(
    vertices, tets, constraint_faces, split_tets);
  splitTaggedBoundaryFaces(
    vertices, tets, constraint_edges, constraint_faces, split_faces);
  splitTaggedBoundaryEdges(
    vertices, tets, constraint_edges, split_edges);
  validateSimplicialEmbedding(
    vertices, tets, constraint_edges, constraint_faces);

  Logger::user_logger->info(
    "simplicial embedding done: split {} tetrahedra, {} faces, {} edges; {} vertices and {} tetrahedra remain.",
    split_tets, split_faces, split_edges, vertices.size(), tets.size());

  size_t offset_edge_splits = 0;
  insertOffset(vertices, tets, offset_edge_splits);
  Logger::user_logger->info(
    "offset insertion done: split {} input/outside edges; {} vertices and {} tetrahedra remain.",
    offset_edge_splits, vertices.size(), tets.size());

  rebuildVolumeMesh(
    *mesh, vertices, tets, constraint_edges, constraint_faces);
}

}// namespace CageInit
}// namespace Cage
