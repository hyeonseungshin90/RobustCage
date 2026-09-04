#pragma once

#include "Mesh/SurfaceMesh/SurfaceMeshDefinition.h"

namespace Cage
{
namespace CageInit
{

enum class BoundaryRailAnchorMode
{
  EdgeMidpoint,
  VertexBisector
};

struct BoundaryRailBuildStats
{
  size_t detectedLoopCount = 0;
  size_t rayValidLoopCount = 0;
  size_t builtLoopCount = 0;
  size_t anchorRequestCount = 0;
  size_t railVertexCount = 0;
  size_t railEdgeCount = 0;
  bool successful = false;
};

// Inserts one cage anchor for every element of every closed source boundary.
// The production mode casts from edge midpoints along per-edge outward
// co-normals; the experimental alternative casts from vertices along adjacent
// co-normal bisectors. Consecutive anchors are connected along cage edges and
// the resulting closed loops are labeled for Phase 2. In both modes source
// edges keep their own (unaveraged) co-normals as the rail support metadata.
class BoundaryRailBuilder
{
public:
  BoundaryRailBuilder(
    SurfaceMesh::SMeshT* source_mesh,
    SurfaceMesh::SMeshT* cage_mesh,
    BoundaryRailAnchorMode anchor_mode = BoundaryRailAnchorMode::EdgeMidpoint)
    : source(source_mesh), cage(cage_mesh), mode(anchor_mode)
  {}

  bool build();
  BoundaryRailBuildStats build_with_stats();

private:
  SurfaceMesh::SMeshT* source;
  SurfaceMesh::SMeshT* cage;
  BoundaryRailAnchorMode mode;
};

}// namespace CageInit
}// namespace Cage
