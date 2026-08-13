#pragma once

#include "Mesh/SurfaceMesh/SurfaceMeshDefinition.h"

namespace Cage
{
namespace CageInit
{

// Inserts one cage anchor for every vertex of every closed source boundary,
// connects consecutive anchors along cage edges, and labels the resulting
// closed edge loops for Phase 2.
class BoundaryRailBuilder
{
public:
  BoundaryRailBuilder(
    SurfaceMesh::SMeshT* source_mesh,
    SurfaceMesh::SMeshT* cage_mesh)
    : source(source_mesh), cage(cage_mesh)
  {}

  bool build();

private:
  SurfaceMesh::SMeshT* source;
  SurfaceMesh::SMeshT* cage;
};

}// namespace CageInit
}// namespace Cage
