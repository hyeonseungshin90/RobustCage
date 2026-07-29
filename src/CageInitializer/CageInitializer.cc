#include "CageInitializer.hh"
#include "CageSimplifier/SpaceSearch/DFaceTree.h"
#include <algorithm>
#include <fstream>
#include <iomanip>

namespace Cage
{
namespace CageInit
{

namespace
{
bool hasExactNonAdjacentSelfIntersection(SM::SMeshT* mesh)
{
  Geometry::LightDFaceTree tree(*mesh);

  for (SM::FaceHandle fh : mesh->faces())
  {
    // Adjacent faces meet at a shared edge or vertex by construction. Ignore
    // all faces sharing a vertex with the query face so this query detects
    // only non-adjacent triangle intersections.
    std::vector<int> ignored_faces;
    for (SM::VertexHandle vh : mesh->fv_range(fh))
    {
      for (SM::FaceHandle vfh : mesh->vf_range(vh))
        ignored_faces.push_back(vfh.idx());
    }

    std::sort(ignored_faces.begin(), ignored_faces.end());
    ignored_faces.erase(
      std::unique(ignored_faces.begin(), ignored_faces.end()),
      ignored_faces.end());

    auto fv_it = mesh->fv_begin(fh);
    const SM::VertexHandle v0 = *fv_it++;
    const SM::VertexHandle v1 = *fv_it++;
    const SM::VertexHandle v2 = *fv_it;

    const auto& p0 = mesh->point(v0);
    const auto& p1 = mesh->point(v1);
    const auto& p2 = mesh->point(v2);
    const auto* ep0 = mesh->data(v0).ep.get();
    const auto* ep1 = mesh->data(v1).ep.get();
    const auto* ep2 = mesh->data(v2).ep.get();

    if (tree.do_intersect(
      p0, p1, p2,
      ep0, ep1, ep2,
      std::move(ignored_faces)))
    {
      Logger::user_logger->critical(
        "Phase 1 exact non-adjacent self-intersection detected: query face {}.",
        fh.idx());
      return true;
    }
  }

  Logger::user_logger->info(
    "Phase 1 exact non-adjacent self-intersection check passed.");
  return false;
}

size_t countBoundaryEdges(SM::SMeshT* mesh)
{
  size_t boundary_edges = 0;
  for (SM::EdgeHandle eh : mesh->edges())
  {
    if (mesh->is_boundary(eh))
      boundary_edges++;
  }
  return boundary_edges;
}

void writeVolumeMeshFacesObj(VM::VMeshT* mesh, const std::string& file_name)
{
  std::fstream fout(file_name, std::fstream::out);
  ASSERT(fout.is_open(), "fail to open {}", file_name);

  fout << std::setprecision(17);
  fout << "# volume mesh face dump\n";
  fout << "# vertices " << mesh->nVertices() << "\n";
  fout << "# faces " << mesh->nFaces() << "\n";
  fout << "# cells " << mesh->nCells() << "\n";

  for (size_t vidx = 0; vidx < mesh->nVertices(); vidx++)
  {
    auto& p = mesh->point(VM::VertexHandle(vidx));
    fout << "v " << p.x() << " " << p.y() << " " << p.z() << "\n";
  }

  size_t written_faces = 0;
  for (size_t fidx = 0; fidx < mesh->nFaces(); fidx++)
  {
    VM::FaceHandle fh(fidx);
    if (mesh->face(fh).nConnCells() == 0)
      continue;

    auto fv = mesh->findFV(fh);
    fout << "f "
      << fv[0].idx() + 1 << " "
      << fv[1].idx() + 1 << " "
      << fv[2].idx() + 1 << "\n";
    written_faces++;
  }

  fout.close();
  Logger::user_logger->info("wrote volume mesh face OBJ: {} ({} faces)", file_name, written_faces);
}
}

CageInitializer::CageInitializer()
{
  SMesh = nullptr;
  param = nullptr;
}

CageInitializer::CageInitializer(
  SM::SMeshT* _SMesh, ParamCageInitializer* _param,
  VM::VMeshT* _outVMesh, SM::SMeshT* _outSMesh
)
{
  SMesh = _SMesh;
  param = _param;
  outVMesh = _outVMesh;
  outSMesh = _outSMesh;
}

void CageInitializer::generate()
{
  Logger::user_logger->info("begin generating initial cage.");

  SM::pre_calculate_edge_length(SMesh);
  SM::pre_calculate_face_area(SMesh);
  getBoundingBox();
  const size_t input_boundary_edges = countBoundaryEdges(SMesh);

  // step 1. tetrahedralize space around surface mesh.
  tetrahedralizer = std::make_unique<Tetrahedralizer>(SMesh, &param->paramTetrahedralizer, outVMesh);
  tetrahedralizer->tetrahedralize();
  tetrahedralizer.reset();
  writeVolumeMeshFacesObj(
    outVMesh,
    param->fileOutPath + param->fileName + "_debug_tetrahedralize.obj");

  tetrahedralizationPostProcess();

  if (input_boundary_edges == 0)
  {
    Logger::user_logger->info("separating volume mesh to inside and outside.");
    separateMeshToInOut();
  }
  else
  {
    Logger::user_logger->warn(
      "input mesh is non-watertight ({} boundary edges); "
      "skip inside/outside separation and extract cage from constraint-adjacent tetrahedra.",
      input_boundary_edges);
  }

  // step 2.1. trim tetrahedral mesh, including subdiving and removing tets.
  tetMeshTrimmer = std::make_unique<TetMeshTrimmer>(outVMesh);
  tetMeshTrimmer->trim();
  writeVolumeMeshFacesObj(
    outVMesh,
    param->fileOutPath + param->fileName + "_debug_trim.obj");

  // step 2.2. retrieve cage from tetrahedral mesh.
  retrieveCage(outVMesh, outSMesh);
  if (hasExactNonAdjacentSelfIntersection(outSMesh))
  {
    throw std::logic_error(
      "Phase 1 cage has an exact non-adjacent self-intersection.");
  }

  const std::string retrieve_cage_path =
    param->fileOutPath + param->fileName + "_debug_retrieve_cage.obj";
  OpenMesh::IO::write_mesh(*outSMesh, retrieve_cage_path, OpenMesh::IO::Options::Default, 15);
  Logger::user_logger->info("wrote retrieved cage OBJ: {}", retrieve_cage_path);
  tetMeshTrimmer.reset();
  {Logger::user_logger->info("generating initial cage done!");}
  {Logger::user_logger->info("peak memory used: {} MB", getPeakMegabytesUsed());}
}

///@brief get bounding box for tri_mesh
void CageInitializer::getBoundingBox()
{
  bbox = BoundingBox(SMesh->point(SM::VertexHandle(0)), SMesh->point(SM::VertexHandle(0)));
  for (SM::VertexHandle vh : SMesh->vertices())
    bbox += SMesh->point(vh);
  bbox_diag_length = (bbox.max() - bbox.min()).norm();
}

/// @brief do post-process after tetrahedralization:
/// 1. Detect constraint edges and constraint vertices by constraint faces.
///    (Edges and vertices adjacent to constraint faces are constraits.)
/// 2. Detect tetrahedral mesh's boundary face.
///    (All faces adjacent to exactly one tetrahedron are boundary face.)
/// 3. Detect tetrahedrons that are adjacent to constraint.
///    (Tets adjacent to constraint vertex are marked.)
void CageInitializer::tetrahedralizationPostProcess()
{
  for (auto& face : outVMesh->faces)
  {
    // 1. detect constraint edges and vertices.
    if (face.prop.is_constraint)
    {
      for (auto eh : face.halfedges)
      {
        auto& edge = outVMesh->edge(eh);
        edge.prop.is_constraint = true;
        outVMesh->vertex(edge.from()).prop.is_constraint = true;
        outVMesh->vertex(edge.to()).prop.is_constraint = true;
      }
    }
  }
  // 2. detect tet mesh's boundary
  std::vector<size_t> face_ref_count(outVMesh->nFaces(), 0);
  for (auto& cell : outVMesh->cells)
  {
    for (auto hfh : cell.halffaces)
    {
      face_ref_count[faceHdl(hfh).idx()]++;
    }
  }

  for (size_t fidx = 0;fidx < outVMesh->nFaces();fidx++)
  {
    if (face_ref_count[fidx] == 1)
    {
      auto& face = outVMesh->faces[fidx];
      face.prop.is_boundary = true;
      for (auto eh : face.halfedges)
      {
        auto& edge = outVMesh->edge(eh);
        edge.prop.is_on_boundary = true;
        outVMesh->vertex(edge.from()).prop.is_on_boundary = true;
        outVMesh->vertex(edge.to()).prop.is_on_boundary = true;
      }
    }
  }
  // 3. detect tet adjacent to constraint.
  for (auto& cell : outVMesh->cells)
  {
    for (auto fh : cell.faces())
    {
      if (outVMesh->face(fh).prop.is_constraint)
      {
        cell.prop.is_adjacent_constraint = true;
        goto NEXT_CELL;
      }
      for (auto eh : outVMesh->face(fh).edges())
      {
        if (outVMesh->edge(eh).prop.is_constraint)
        {
          cell.prop.is_adjacent_constraint = true;
          goto NEXT_CELL;
        }
        for (auto vh : outVMesh->edge(eh).vertices)
        {
          if (outVMesh->vertex(vh).prop.is_constraint)
          {
            cell.prop.is_adjacent_constraint = true;
            goto NEXT_CELL;
          }
        }
      }
    }
  NEXT_CELL:;
  }
}

/// @brief separate volume mesh into inside and outside.
void CageInitializer::separateMeshToInOut()
{
  REP_FUNC;
  tetRemover = std::make_unique<TetRemover>(outVMesh);
  std::vector<VM::CellHandle> cellsToDelete;
  cellsToDelete.reserve(outVMesh->nCells() / 2);

  // collect inside tets in outVMesh
  for (int cidx = 0;cidx < outVMesh->nCells();cidx++)
  {
    if (outVMesh->cells[cidx].prop.is_inside)
      cellsToDelete.emplace_back(cidx);
  }
  tetRemover->removeTet(cellsToDelete);
  Logger::user_logger->info("outside cage: {} vertices, {} faces.",
    outVMesh->nVertices(), outVMesh->nFaces());
  tetRemover.reset();
}

/// @brief retrieve boundary of outside/inside mesh.
/// @note won't release memory of volume mesh.
void CageInitializer::retrieveBoundary(VM::VMeshT* vol_mesh, SM::SMeshT* boundary_mesh)
{
  std::vector<SM::VertexHandle> vertex_vol2surf(vol_mesh->nVertices());

  // add vertices on boundary
  for (int vidx = 0;vidx < vol_mesh->nVertices();vidx++)
  {
    if (vol_mesh->vertices[vidx].prop.is_on_boundary && !vol_mesh->vertices[vidx].prop.is_constraint)
    {
      VM::VertexHandle vh(vidx);
      SM::VertexHandle bvh = boundary_mesh->add_vertex(vol_mesh->point(vh));
      if (!vol_mesh->exact_point(vh).isSame())
      {
        boundary_mesh->data(bvh).ep = std::make_unique<ExactPoint>(vol_mesh->exact_point(vh));
      }
      vertex_vol2surf[vidx] = bvh;
    }
  }
  // add faces on boundary
  for (int fidx = 0;fidx < vol_mesh->nFaces();fidx++)
  {
    if (vol_mesh->faces[fidx].prop.is_boundary && !vol_mesh->faces[fidx].prop.is_constraint)
    {
      VM::FaceHandle fh(fidx);
      // find halfface on boundary
      auto& conn_cell = vol_mesh->cell(vol_mesh->face(fh).connCells()[0]);
      ASSERT(vol_mesh->face(fh).nConnCells() == 1, "not boundary face.");
      std::array<VM::VertexHandle, 3> fvh;
      if (std::find(conn_cell.halffaces.begin(), conn_cell.halffaces.end(),
        halfFaceHdl(fh, VM::Sequence)) != conn_cell.halffaces.end())
        fvh = vol_mesh->findHFV(halfFaceHdl(fh, VM::Reverse));
      else
        fvh = vol_mesh->findHFV(halfFaceHdl(fh, VM::Sequence));
      // add face into boundary mesh
      auto surf_fh = boundary_mesh->add_face(
        SM::VertexHandle(vertex_vol2surf[fvh[0].idx()]),
        SM::VertexHandle(vertex_vol2surf[fvh[1].idx()]),
        SM::VertexHandle(vertex_vol2surf[fvh[2].idx()]));
      ASSERT(surf_fh.is_valid(), "adding vol face {} failed", fidx);
    }
  }
  Logger::user_logger->info("retrieve boundary surface: {} vertices and {} faces.",
    boundary_mesh->n_vertices(), boundary_mesh->n_faces());
}

/// @brief retrieve cage from volume mesh. tets non-adjacent to constraint are not removed.
/// @note won't release memory of volume mesh.
void CageInitializer::retrieveCage(VM::VMeshT* vol_mesh, SM::SMeshT* boundary_mesh)
{
  // 1. find faces of cage
  std::vector<VM::HalfFaceHandle> cage_faces;

  for (int fidx = 0;fidx < vol_mesh->nFaces();fidx++)
  {
    if (vol_mesh->faces[fidx].prop.is_constraint)
      continue;

    // count connect cells adjacent to constraint   
    VM::FaceHandle fh(fidx);
    VM::VMeshT::FT& face = vol_mesh->faces[fidx];
    size_t constr_cnt = 0;
    constr_cnt += face.conn_cells[0].isValid() && vol_mesh->cell(face.conn_cells[0]).prop.is_adjacent_constraint;
    constr_cnt += face.conn_cells[1].isValid() && vol_mesh->cell(face.conn_cells[1]).prop.is_adjacent_constraint;

    if (constr_cnt == 1)
    {
      if (face.conn_cells[0].isValid() && vol_mesh->cell(face.conn_cells[0]).prop.is_adjacent_constraint)
        cage_faces.push_back(halfFaceHdl(fh, VM::Reverse));
      else
        cage_faces.push_back(halfFaceHdl(fh, VM::Sequence));
    }
  }
  // 2. insert vertices and faces to surface mesh
  std::vector<SM::VertexHandle> vertex_vol2surf(vol_mesh->nVertices(), SM::SMeshT::InvalidVertexHandle);
  for (VM::HalfFaceHandle hfh : cage_faces)
  {
    std::array<VM::VertexHandle, 3> fvs = vol_mesh->findHFV(hfh);
    for (size_t i = 0;i < 3;i++) if (!vertex_vol2surf[fvs[i].idx()].is_valid())
    {
      vertex_vol2surf[fvs[i].idx()] = boundary_mesh->add_vertex(vol_mesh->point(fvs[i]));
      if (!vol_mesh->exact_point(fvs[i]).isSame())
      {
        boundary_mesh->data(vertex_vol2surf[fvs[i].idx()]).ep =
          std::make_unique<ExactPoint>(vol_mesh->exact_point(fvs[i]));
      }
    }
    auto surf_fh = boundary_mesh->add_face(
      vertex_vol2surf[fvs[0].idx()],
      vertex_vol2surf[fvs[1].idx()],
      vertex_vol2surf[fvs[2].idx()]);
    ASSERT(surf_fh.is_valid(), "adding vol face {} failed", hfh.idx());
  }
  Logger::user_logger->info("retrieve boundary surface: {} vertices and {} faces.",
    boundary_mesh->n_vertices(), boundary_mesh->n_faces());

  // Check possible error from VolumeMesher.
  // Sometimes VolumeMesher will mark some inside tets as outside. Perhaps there are bugs in VolumeMesher's exact arithmetic lib. 
  // You can drop the inside tets to fix this problem, but we report the bug here. Hope it will be fixed.
  /*for (SM::EdgeHandle eh : boundary_mesh->edges())
  {
    if (boundary_mesh->is_boundary(eh))
    {
      Logger::user_logger->critical("find boundary edge in boundary of tetrahedral mesh.");
      throw std::logic_error("find boundary edge in boundary of tetrahedral mesh.");
    }
  }*/
}

}// namespace CageInit
}// namespace Cage 
