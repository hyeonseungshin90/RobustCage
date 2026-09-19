#pragma once
#include "CageInitializer/CageInitializer.hh"
#include "CageInitializer/BoundaryRailBuilder.hh"
#include "CageInitializer/BoundaryRailExporter.hh"
#include "CageSimplifier/CageSimplifier.hh"

namespace Cage
{
using namespace SimpleUtils;
using namespace Geometry;
using namespace CageInit;
using namespace CageSimp;
namespace SM = SurfaceMesh;
namespace VM = VolumeMesh;

class CageGenerator
{
public:
  // input
  std::unique_ptr<SM::SMeshT> originalMesh;
  ParamCageGenerator param;
  // Optional files of an earlier run.  A cage path skips Phase 1; a rail
  // index file (with the cage it refers to) also skips rail construction.
  std::string inputCagePath;
  std::string inputRailPath;

  // sub components
  std::unique_ptr<CageInitializer> cageInitializer;
  std::unique_ptr<CageSimplifier> cageSimplifier;

  // middle results
  std::unique_ptr<VM::VMeshT> VMesh;

  // result
  std::unique_ptr<SM::SMeshT> cage;
public:
  void stageInitialize();
  void stageLoadInitialCage();
  void stageBuildBoundaryRails();
  void stageLoadBoundaryRails();
  void stageExportInitialBoundaryRails();
  void stageSimplify();

  void generate();
};


}// namespace Cage
