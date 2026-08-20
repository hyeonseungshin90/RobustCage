#pragma once

#include "Mesh/SurfaceMesh/SurfaceMeshDefinition.h"

namespace Cage
{
namespace CageInit
{

// Inserts one cage anchor for every edge of every closed source boundary.
// Each ray starts at the source edge midpoint and follows that edge's outward
// co-normal. Consecutive anchors are then connected along cage edges and the
// resulting closed edge loops are labeled for Phase 2. The source edges keep
// the same labels and co-normals so rail collapses can use their ruled support
// surface instead of shrinking onto cage-edge chords.
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
