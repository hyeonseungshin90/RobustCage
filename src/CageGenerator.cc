#include "CageGenerator.hh"
#include <chrono>
#include "omp.h"

namespace Cage
{
void CageGenerator::stageInitialize()
{
  VMesh = std::make_unique<VM::VMeshT>();
  cage = std::make_unique<SM::SMeshT>();

  cageInitializer = std::make_unique<CageInitializer>(
    originalMesh.get(), &param.paramCageInitializer,
    VMesh.get(), cage.get());

  cageInitializer->generate();
}

void CageGenerator::stageBuildBoundaryRails()
{
  const BoundaryRailAnchorMode anchor_mode =
    param.paramCageSimplifier.boundaryRailAnchorMode == "vertex" ?
      BoundaryRailAnchorMode::VertexBisector :
      BoundaryRailAnchorMode::EdgeMidpoint;
  BoundaryRailBuilder rail_builder(originalMesh.get(), cage.get(), anchor_mode);
  if (!rail_builder.build())
  {
    Logger::user_logger->warn(
      "one or more boundary rails could not be constructed; continuing with the successfully constructed rails.");
  }
}

// Phase 2 rebuilds the rails along with the cage, so the rails as Phase 1 left
// them only exist here.  Write them next to the final rails of the cage so that
// both states of the same rail ids can be compared.
void CageGenerator::stageExportInitialBoundaryRails()
{
  const ParamCageSimplifier& simplifier_param = param.paramCageSimplifier;
  const std::string prefix =
    simplifier_param.fileOutPath + simplifier_param.fileName + "_cage_" +
    std::to_string(simplifier_param.cageLabel) + "_initial_rails";
  const std::string rail_obj_path = prefix + ".obj";
  const std::string rail_txt_path = prefix + ".txt";
  const BoundaryRailExport rail_export = write_boundary_rail_files(
    *cage, rail_obj_path, rail_txt_path,
    "the Phase 1 cage (before Phase 2 simplification)");
  if (rail_export.rail_count > 0)
  {
    Logger::user_logger->info(
      "wrote {} initial boundary rails ({} rail vertices, {} rail edges) to {} and {}.",
      rail_export.rail_count, rail_export.vertex_count, rail_export.edge_count,
      rail_obj_path, rail_txt_path);
  }
  else
  {
    Logger::user_logger->warn(
      "boundary rails were enabled but Phase 1 produced no rail vertex or edge; no initial rail file written.");
  }
}

void CageGenerator::stageSimplify()
{
  cageSimplifier = std::make_unique<CageSimplifier>(
    originalMesh.get(), cage.get(), &param.paramCageSimplifier);

  cageSimplifier->simplify();
}

void CageGenerator::generate()
{
  omp_set_num_threads(12);

  const auto phase1_start = std::chrono::steady_clock::now();
  stageInitialize();
  if (param.paramCageSimplifier.enableBoundaryRails)
  {
    stageBuildBoundaryRails();
    stageExportInitialBoundaryRails();
  }
  const double phase1_seconds = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - phase1_start).count();

  cageInitializer = nullptr;
  VMesh = nullptr;

  const auto phase2_start = std::chrono::steady_clock::now();
  stageSimplify();
  const double phase2_seconds = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - phase2_start).count();

  Logger::user_logger->info(
    "Phase 1 elapsed time: {:.6f} seconds.", phase1_seconds);
  Logger::user_logger->info(
    "Phase 2 elapsed time: {:.6f} seconds.", phase2_seconds);
}

}// namespace Cage
