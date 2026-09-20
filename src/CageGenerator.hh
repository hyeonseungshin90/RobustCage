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
  // Optional files of an earlier run.  A cage path skips Phase 1 (initial
  // cage); a rail index file (with the cage it refers to) also skips Phase 2
  // (rail construction).
  std::string inputCagePath;
  std::string inputRailPath;

  // sub components
  std::unique_ptr<CageInitializer> cageInitializer;
  std::unique_ptr<CageSimplifier> cageSimplifier;

  // middle results
  std::unique_ptr<VM::VMeshT> VMesh;

  // result
  std::unique_ptr<SM::SMeshT> cage;

  // Computation time of each phase of the last generate(); see PhaseTimer.
  PhaseTimer phase1Timer;
  PhaseTimer phase2Timer;
  PhaseTimer phase3Timer;
public:
  void stageInitialize();
  void stageLoadInitialCage();
  void stageBuildBoundaryRails();
  void stageLoadBoundaryRails();
  void stageExportPhase2();
  void stageSimplify();

  void generate();
  // Logs the phase times of the last generate() and writes them to
  // <model>_timing.json.
  void reportPhaseTimes() const;

  // Output path without extension for one stage of the current cage label,
  // e.g. <out dir>/bunny_phase2_rails, or bunny_phase2_rails_1 for the
  // second nested cage.
  std::string stageOutputPath(const std::string& stage) const;
  // Its file name part, e.g. bunny_phase2_rails.
  std::string stageOutputName(const std::string& stage) const;
};


}// namespace Cage
