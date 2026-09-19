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

void CageGenerator::stageLoadInitialCage()
{
  cage = std::make_unique<SM::SMeshT>();

  cageInitializer = std::make_unique<CageInitializer>(
    originalMesh.get(), &param.paramCageInitializer,
    nullptr, cage.get());

  cageInitializer->load(inputCagePath);
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

void CageGenerator::stageLoadBoundaryRails()
{
  Logger::user_logger->info(
    "skip boundary rail construction and load rails: {}", inputRailPath);
  const BoundaryRailExport rail_import =
    read_boundary_rail_files(*cage, *originalMesh, inputRailPath);
  Logger::user_logger->info(
    "loaded {} boundary rails ({} rail vertices, {} rail edges, {} source boundary edges).",
    rail_import.rail_count, rail_import.vertex_count, rail_import.edge_count,
    rail_import.support_edge_count);
  if (rail_import.rail_count == 0)
  {
    Logger::user_logger->warn(
      "boundary rails were enabled but the rail file lists no rail.");
  }
}

// Phase 3 rebuilds the rails along with the cage, so the rails as Phase 2 left
// them only exist here.  Write them next to the final rails of the cage so that
// both states of the same rail ids can be compared.  Anchor insertion and the
// geodesic embedding split the Phase 1 cage, so the rail indices refer to the
// cage written here, not to the retrieved cage; the pair can be passed back
// with --cage/--rails to skip Phases 1 and 2.
void CageGenerator::stageExportInitialBoundaryRails()
{
  const ParamCageSimplifier& simplifier_param = param.paramCageSimplifier;
  const std::string name = simplifier_param.fileName + "_cage_" +
    std::to_string(simplifier_param.cageLabel) + "_initial";
  const std::string prefix = simplifier_param.fileOutPath + name;
  const std::string cage_obj_name = name + ".obj";
  const std::string cage_obj_path = prefix + ".obj";
  if (write_cage_obj(*cage, cage_obj_path))
  {
    Logger::user_logger->info(
      "wrote initial cage with boundary rails inserted ({} vertices, {} faces) to {}.",
      cage->n_vertices(), cage->n_faces(), cage_obj_path);
  }

  const std::string rail_obj_path = prefix + "_rails.obj";
  const std::string rail_txt_path = prefix + "_rails.txt";
  const BoundaryRailExport rail_export = write_boundary_rail_files(
    *cage, rail_obj_path, rail_txt_path, cage_obj_name, originalMesh.get());
  if (rail_export.rail_count > 0)
  {
    Logger::user_logger->info(
      "wrote {} initial boundary rails ({} rail vertices, {} rail edges, {} source boundary edges) to {} and {}.",
      rail_export.rail_count, rail_export.vertex_count, rail_export.edge_count,
      rail_export.support_edge_count, rail_obj_path, rail_txt_path);
  }
  else
  {
    Logger::user_logger->warn(
      "boundary rails were enabled but Phase 2 produced no rail vertex or edge; no initial rail file written.");
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

  const auto seconds_since = [](std::chrono::steady_clock::time_point start)
  {
    return std::chrono::duration<double>(
      std::chrono::steady_clock::now() - start).count();
  };

  // Phase 1: initial cage.
  const auto phase1_start = std::chrono::steady_clock::now();
  if (inputCagePath.empty())
    stageInitialize();
  else
    stageLoadInitialCage();
  const double phase1_seconds = seconds_since(phase1_start);

  // Phase 2: boundary rails.
  const bool rails_enabled = param.paramCageSimplifier.enableBoundaryRails;
  double phase2_seconds = 0.0;
  if (rails_enabled)
  {
    const auto phase2_start = std::chrono::steady_clock::now();
    if (inputRailPath.empty())
      stageBuildBoundaryRails();
    else
      stageLoadBoundaryRails();
    stageExportInitialBoundaryRails();
    phase2_seconds = seconds_since(phase2_start);
  }

  cageInitializer = nullptr;
  VMesh = nullptr;

  // Phase 3: simplification.
  const auto phase3_start = std::chrono::steady_clock::now();
  stageSimplify();
  const double phase3_seconds = seconds_since(phase3_start);

  if (inputCagePath.empty())
  {
    Logger::user_logger->info(
      "Phase 1 elapsed time: {:.6f} seconds.", phase1_seconds);
  }
  else
  {
    Logger::user_logger->info(
      "Phase 1 skipped: initial cage loaded from file in {:.6f} seconds.",
      phase1_seconds);
  }
  if (!rails_enabled)
    Logger::user_logger->info("Phase 2 skipped: boundary rails disabled.");
  else if (inputRailPath.empty())
  {
    Logger::user_logger->info(
      "Phase 2 elapsed time: {:.6f} seconds.", phase2_seconds);
  }
  else
  {
    Logger::user_logger->info(
      "Phase 2 skipped: boundary rails loaded from file in {:.6f} seconds.",
      phase2_seconds);
  }
  Logger::user_logger->info(
    "Phase 3 elapsed time: {:.6f} seconds.", phase3_seconds);
}

}// namespace Cage
