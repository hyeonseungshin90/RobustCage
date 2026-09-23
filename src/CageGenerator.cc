#include "CageGenerator.hh"
#include <fstream>
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

  // This file can be passed back with --cage to skip Phase 1.
  const std::string cage_path = stageOutputPath("phase1_cage") + ".obj";
  if (write_cage_obj(*cage, cage_path))
  {
    Logger::user_logger->info(
      "wrote Phase 1 cage ({} vertices, {} faces) to {}.",
      cage->n_vertices(), cage->n_faces(), cage_path);
  }
}

void CageGenerator::stageLoadInitialCage()
{
  cage = std::make_unique<SM::SMeshT>();

  cageInitializer = std::make_unique<CageInitializer>(
    originalMesh.get(), &param.paramCageInitializer,
    nullptr, cage.get());

  cageInitializer->load(inputCagePath);
}

size_t CageGenerator::stageCleanupCage()
{
  SM::pre_calculate_edge_length(originalMesh.get());
  SM::pre_calculate_face_area(originalMesh.get());
  SM::pre_calculate_edge_length(cage.get());
  SM::pre_calculate_face_area(cage.get());
  CageSimp::init_one_ring_faces(cage.get());
  if (!cage->has_face_normals())
    cage->request_face_normals();
  if (!cage->has_vertex_normals())
    cage->request_vertex_normals();
  cage->update_normals();

  VertexTree vertex_tree(*originalMesh);
  DFaceTree source_tree(*originalMesh);
  FaceGrid source_grid(*originalMesh);
  LightDFaceTree cage_tree(*cage);
  Vec3d lo(DBL_MAX, DBL_MAX, DBL_MAX);
  Vec3d hi(-DBL_MAX, -DBL_MAX, -DBL_MAX);
  for (SM::VertexHandle vh : originalMesh->vertices())
  {
    lo.minimize(originalMesh->point(vh));
    hi.maximize(originalMesh->point(vh));
  }

  CageSimp::DegenerationRemover remover(
    originalMesh.get(), cage.get(), &vertex_tree, &source_tree, &cage_tree,
    &source_grid);
  remover.original_diagonal_length = (hi - lo).norm();
  return remover.perform();
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
// Phase 2 cage written here, not to the Phase 1 cage; the pair can be passed
// back with --cage/--rails to skip Phases 1 and 2.
void CageGenerator::stageExportPhase2()
{
  const std::string cage_obj_name = stageOutputName("phase2_cage") + ".obj";
  const std::string cage_obj_path = stageOutputPath("phase2_cage") + ".obj";
  if (write_cage_obj(*cage, cage_obj_path))
  {
    Logger::user_logger->info(
      "wrote Phase 2 cage with boundary rails inserted ({} vertices, {} faces) to {}.",
      cage->n_vertices(), cage->n_faces(), cage_obj_path);
  }

  const std::string rail_path = stageOutputPath("phase2_rails");
  const std::string rail_obj_path = rail_path + ".obj";
  const std::string rail_txt_path = rail_path + ".txt";
  const BoundaryRailExport rail_export = write_boundary_rail_files(
    *cage, rail_obj_path, rail_txt_path, cage_obj_name, originalMesh.get());
  if (rail_export.rail_count > 0)
  {
    Logger::user_logger->info(
      "wrote {} Phase 2 boundary rails ({} rail vertices, {} rail edges, {} source boundary edges) to {} and {}.",
      rail_export.rail_count, rail_export.vertex_count, rail_export.edge_count,
      rail_export.support_edge_count, rail_obj_path, rail_txt_path);
  }
  else
  {
    Logger::user_logger->warn(
      "boundary rails were enabled but Phase 2 produced no rail vertex or edge; no Phase 2 rail file written.");
  }
}

void CageGenerator::stageSimplify()
{
  cageSimplifier = std::make_unique<CageSimplifier>(
    originalMesh.get(), cage.get(), &param.paramCageSimplifier);

  cageSimplifier->simplify(usesDefaultCleanup());
}

bool CageGenerator::usesDefaultCleanup() const
{
  return param.paramCageInitializer.phase1Mode == "subdivision" &&
    param.paramCageSimplifier.phase3Mode == "fast" &&
    !param.paramCageSimplifier.enableBoundaryRails;
}

void CageGenerator::generate()
{
  omp_set_num_threads(12);
  phase1Timer = PhaseTimer();
  phase2Timer = PhaseTimer();
  phase3Timer = PhaseTimer();

  // Phase 1: initial cage.
  phase1Timer.start();
  if (inputCagePath.empty())
    stageInitialize();
  else
  {
    PhaseTimer::Exclusion exclusion("input");
    stageLoadInitialCage();
  }
  phase1Timer.stop();

  // Keep the original default pipeline's cleanup. Other phase combinations
  // skip it; loaded rails must also retain the cage indices they refer to.
  cleanupBeforePhase2Seconds = 0.0;
  cleanupBeforePhase2Cases = 0;
  if (usesDefaultCleanup() && inputRailPath.empty())
  {
    const auto cleanup_start = PhaseTimer::Clock::now();
    cleanupBeforePhase2Cases = stageCleanupCage();
    cleanupBeforePhase2Seconds = std::chrono::duration<double>(
      PhaseTimer::Clock::now() - cleanup_start).count();
    Logger::user_logger->info(
      "degeneracy cleanup before Phase 2: {} cases in {:.6f} seconds (not counted in the phase time).",
      cleanupBeforePhase2Cases, cleanupBeforePhase2Seconds);
  }
  else if (!usesDefaultCleanup())
  {
    Logger::user_logger->info(
      "degeneracy cleanup before Phase 2 disabled for non-default pipeline.");
  }

  // Phase 2: boundary rails.
  if (param.paramCageSimplifier.enableBoundaryRails)
  {
    phase2Timer.start();
    if (inputRailPath.empty())
      stageBuildBoundaryRails();
    else
    {
      PhaseTimer::Exclusion exclusion("input");
      stageLoadBoundaryRails();
    }
    phase2Timer.stop();
    stageExportPhase2();
  }

  cageInitializer = nullptr;
  VMesh = nullptr;

  // Phase 3: simplification.
  phase3Timer.start();
  stageSimplify();
  phase3Timer.stop();

  reportPhaseTimes();
}

void CageGenerator::reportPhaseTimes() const
{
  const auto excluded_text = [](const PhaseTimer& timer)
  {
    std::string text;
    for (const auto& entry : timer.excludedSeconds())
      text += fmt::format("{} {:.6f} s, ", entry.first, entry.second);
    return text + fmt::format("logging {:.6f} s", timer.loggingSeconds());
  };
  const auto phase_json = [](const PhaseTimer& timer, const char* status)
  {
    boost::json::object excluded;
    for (const auto& entry : timer.excludedSeconds())
      excluded[entry.first] = entry.second;
    excluded["logging"] = timer.loggingSeconds();
    boost::json::object phase;
    phase["status"] = status;
    phase["seconds"] = timer.wasStarted() ? timer.seconds() : 0.0;
    phase["excluded"] = excluded;
    return phase;
  };

  const bool phase1_loaded = !inputCagePath.empty();
  const bool rails_enabled = param.paramCageSimplifier.enableBoundaryRails;
  const bool phase2_loaded = rails_enabled && !inputRailPath.empty();
  const auto input_seconds = [](const PhaseTimer& timer)
  {
    const auto found = timer.excludedSeconds().find("input");
    return found == timer.excludedSeconds().end() ? 0.0 : found->second;
  };

  if (phase1_loaded)
  {
    Logger::user_logger->info(
      "Phase 1 skipped: initial cage loaded from file in {:.6f} seconds.",
      input_seconds(phase1Timer));
  }
  else
  {
    Logger::user_logger->info(
      "Phase 1 elapsed time: {:.6f} seconds; excluded: {}.",
      phase1Timer.seconds(), excluded_text(phase1Timer));
  }
  if (!rails_enabled)
    Logger::user_logger->info("Phase 2 skipped: boundary rails disabled.");
  else if (phase2_loaded)
  {
    Logger::user_logger->info(
      "Phase 2 skipped: boundary rails loaded from file in {:.6f} seconds.",
      input_seconds(phase2Timer));
  }
  else
  {
    Logger::user_logger->info(
      "Phase 2 elapsed time: {:.6f} seconds; excluded: {}.",
      phase2Timer.seconds(), excluded_text(phase2Timer));
  }
  Logger::user_logger->info(
    "Phase 3 elapsed time: {:.6f} seconds; excluded: {}.",
    phase3Timer.seconds(), excluded_text(phase3Timer));

  boost::json::object timing;
  timing["unit"] = "seconds";
  timing["definition"] =
    "computation time: wall time minus log output and the excluded sections "
    "(file input/output and verification the cage does not need)";
  timing["phase1"] = phase_json(phase1Timer, phase1_loaded ? "loaded" : "run");
  timing["phase2"] = phase_json(
    phase2Timer, !rails_enabled ? "disabled" : (phase2_loaded ? "loaded" : "run"));
  timing["phase3"] = phase_json(phase3Timer, "run");
  boost::json::object cleanup;
  boost::json::object before_phase2;
  before_phase2["status"] = !usesDefaultCleanup() ? "disabled" :
    (inputRailPath.empty() ? "run" : "skipped");
  before_phase2["seconds"] = cleanupBeforePhase2Seconds;
  before_phase2["cases"] = cleanupBeforePhase2Cases;
  boost::json::object in_phase3;
  in_phase3["status"] = usesDefaultCleanup() ? "run" : "disabled";
  in_phase3["seconds"] = cageSimplifier ? cageSimplifier->degeneracyCleanupSeconds : 0.0;
  in_phase3["cases"] = cageSimplifier ? cageSimplifier->degeneracyCleanupCases : size_t(0);
  cleanup["before_phase2"] = before_phase2;
  cleanup["phase3"] = in_phase3;
  cleanup["note"] =
    "degeneracy cleanup of the cage; reported on its own and not part of the phase times";
  timing["cleanup"] = cleanup;
  timing["total_seconds"] =
    (phase1_loaded ? 0.0 : phase1Timer.seconds()) +
    (rails_enabled && !phase2_loaded ? phase2Timer.seconds() : 0.0) +
    phase3Timer.seconds();
  const std::string timing_path = stageOutputPath("timing") + ".json";
  std::ofstream timing_file(timing_path.c_str());
  if (timing_file.is_open())
    timing_file << boost::json::serialize(timing) << "\n";
  else
    Logger::user_logger->warn("fail to write phase times: {}", timing_path);
}

std::string CageGenerator::stageOutputName(const std::string& stage) const
{
  // Only nested cages after the first carry their index, so that they do not
  // overwrite each other's files.
  const ParamCageSimplifier& simplifier_param = param.paramCageSimplifier;
  std::string name = simplifier_param.fileName + "_" + stage;
  if (simplifier_param.cageLabel > 0)
    name += "_" + std::to_string(simplifier_param.cageLabel);
  return name;
}

std::string CageGenerator::stageOutputPath(const std::string& stage) const
{
  return param.paramCageSimplifier.fileOutPath + stageOutputName(stage);
}

}// namespace Cage
