#pragma once
#include "Tetrahedralizer.hh"
#include "TetMeshTrimmer.hh"
#include "TopologicalOffsetInitializer.hh"

namespace Cage
{
namespace CageInit
{
using namespace SimpleUtils;
namespace SM = SurfaceMesh;
namespace VM = VolumeMesh;

// Write a cage as OBJ with round-trip double precision so that it can be read
// back as an initial cage.  OpenMesh's OBJ writer converts points to float,
// which would move the vertices of a cage that a later run resumes from.
// Vertex i of the file is the cage vertex with index i - 1.
bool write_cage_obj(SM::SMeshT& cage, const std::string& path);

class CageInitializer
{
public:
  // input
  SM::SMeshT* SMesh;
  ParamCageInitializer* param;

  // sub components
  std::unique_ptr<Tetrahedralizer> tetrahedralizer;
  std::unique_ptr<TetRemover> tetRemover;
  std::unique_ptr<TetMeshTrimmer> tetMeshTrimmer;
  std::unique_ptr<TopologicalOffsetInitializer> topologicalOffsetInitializer;

  // middle results
  VM::VMeshT* outVMesh;

  // final results
  SM::SMeshT* outSMesh;
public:
  CageInitializer();
  CageInitializer(
    SM::SMeshT* _SMesh, ParamCageInitializer* _param,
    VM::VMeshT* _outVMesh, SM::SMeshT* _outSMesh
  );

  void generate();
  // Skip Phase 1: read the initial cage of an earlier run into outSMesh and
  // apply the Phase 1 exit checks to it.
  void load(const std::string& cage_path);

private:
  BoundingBox bbox;
  double bbox_diag_length;
  void getBoundingBox();

  void tetrahedralizationPostProcess();
  void separateMeshToInOut();
  void retrieveBoundary(VM::VMeshT* vol_mesh, SM::SMeshT* boundary_mesh);
  void retrieveCage(VM::VMeshT* vol_mesh, SM::SMeshT* boundary_mesh);
};
}// namespace CageInit
}// namespace Cage
