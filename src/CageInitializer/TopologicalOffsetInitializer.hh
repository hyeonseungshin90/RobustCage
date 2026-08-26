#pragma once

#include "Mesh/VolumeMesh/VolumeMeshDefinition.hh"

namespace Cage
{
namespace CageInit
{
namespace VM = VolumeMesh;

// Builds the initial-cage volume used by retrieveCage() with the
// simplicial-embedding and offset-insertion construction from
// Zint et al., "Topological Offsets" (2025).
//
// The implementation is deliberately specialized to RobustCage's embedded
// triangle constraints.  It mirrors TopologicalOffsets' RegularSpace and
// Marching components without introducing a dependency on WMTK.
class TopologicalOffsetInitializer
{
public:
  explicit TopologicalOffsetInitializer(VM::VMeshT* mesh);

  // Replaces mesh with a conforming tetrahedral mesh containing the inserted
  // topological-offset surface.  Cells between the input and offset are tagged
  // is_adjacent_constraint so CageInitializer::retrieveCage() can extract the
  // surface exactly as it does for the legacy path.
  void generate();

private:
  VM::VMeshT* mesh;
};

}// namespace CageInit
}// namespace Cage
