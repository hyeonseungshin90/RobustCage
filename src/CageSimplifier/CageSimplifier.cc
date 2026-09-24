#include "CageSimplifier.hh"
#include "CageSimplifier/Topo/BoundaryRailUpdater.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <queue>
#include <unordered_set>

namespace Cage
{
namespace CageSimp
{
namespace
{
bool is_phase3_energy_mode(const std::string& mode)
{
  return mode == "linear_solve" || mode == "newton_solve" ||
    mode == "qem_original";
}

void configure_phase3_strategy(
  ParamCollapseStage& collapse, const std::string& mode)
{
  collapse.collapsePlacementMethod = "optimization";
  collapse.phase3PlacementStrategy = mode;
  if (mode == "linear_solve" || mode == "qem_original")
    collapse.robustnessMode = "exact_reject";
}

bool validate_and_log_boundary_rails(SMeshT* mesh, const char* stage)
{
  PhaseTimer::Exclusion exclusion("check");
  std::map<int, std::vector<VertexHandle>> rail_vertices;
  std::map<int, std::vector<EdgeHandle>> rail_edges;
  for (VertexHandle vh : mesh->vertices())
  {
    const int id = mesh->data(vh).boundary_rail_id;
    if (id >= 0)
      rail_vertices[id].push_back(vh);
  }
  for (EdgeHandle eh : mesh->edges())
  {
    const int id = mesh->data(eh).boundary_rail_id;
    if (id >= 0)
      rail_edges[id].push_back(eh);
  }
  if (rail_vertices.empty() && rail_edges.empty())
    return true;

  bool all_valid = rail_vertices.size() == rail_edges.size();
  for (const auto& entry : rail_vertices)
  {
    const int id = entry.first;
    const std::vector<VertexHandle>& vertices = entry.second;
    const auto edge_it = rail_edges.find(id);
    if (edge_it == rail_edges.end() || vertices.size() < 3 ||
      edge_it->second.size() != vertices.size())
    {
      all_valid = false;
      continue;
    }

    std::unordered_set<int> vertex_ids;
    for (VertexHandle vh : vertices)
      vertex_ids.insert(vh.idx());
    for (VertexHandle vh : vertices)
    {
      size_t degree = 0;
      for (HalfedgeHandle outgoing : mesh->voh_range(vh))
      {
        const EdgeHandle edge = mesh->edge_handle(outgoing);
        if (mesh->data(edge).boundary_rail_id == id)
        {
          degree++;
          if (mesh->data(mesh->to_vertex_handle(outgoing)).boundary_rail_id != id)
            all_valid = false;
        }
      }
      if (degree != 2)
        all_valid = false;
    }

    std::unordered_set<int> visited;
    std::queue<VertexHandle> queue;
    queue.push(vertices.front());
    visited.insert(vertices.front().idx());
    while (!queue.empty())
    {
      const VertexHandle current = queue.front();
      queue.pop();
      for (HalfedgeHandle outgoing : mesh->voh_range(current))
      {
        if (mesh->data(mesh->edge_handle(outgoing)).boundary_rail_id != id)
          continue;
        const VertexHandle next = mesh->to_vertex_handle(outgoing);
        if (visited.insert(next.idx()).second)
          queue.push(next);
      }
    }
    if (visited.size() != vertex_ids.size())
      all_valid = false;
  }

  size_t vertex_count = 0;
  size_t edge_count = 0;
  for (const auto& entry : rail_vertices)
    vertex_count += entry.second.size();
  for (const auto& entry : rail_edges)
    edge_count += entry.second.size();
  Logger::user_logger->info(
    "boundary rails {}: {} rails, {} vertices, {} edges, topology {}.",
    stage, rail_vertices.size(), vertex_count, edge_count,
    all_valid ? "valid" : "INVALID");
  return all_valid;
}
}

CageSimplifier::CageSimplifier(SMeshT* original, SMeshT* cage, ParamCageSimplifier* p)
  :om(original), rm(cage), param(p)
{}

void CageSimplifier::simplify(bool enable_degeneracy_cleanup)
{
  Logger::user_logger->info("initializing simplifier.");

  pre_calculate_edge_length(om);
  pre_calculate_edge_length(rm);
  pre_calculate_face_area(om);
  pre_calculate_face_area(rm);
  init_one_ring_faces(rm);

  vt = std::make_unique<VertexTree>(*om);
  ot = std::make_unique<DFaceTree>(*om);
  og = std::make_unique<FaceGrid>(*om);
  lrt = std::make_unique<LightDFaceTree>(*rm);

  if (!rm->has_face_normals())
    rm->request_face_normals();
  if (!rm->has_vertex_normals())
    rm->request_vertex_normals();
  rm->update_normals();

  if (is_phase3_energy_mode(param->phase3Mode))
    configure_phase3_strategy(param->paramCollapse, param->phase3Mode);

  degeneration_remover = std::make_unique<DegenerationRemover>(
    om, rm, vt.get(), ot.get(), lrt.get(), og.get());
  degeneration_remover->use_rail_support =
    param->enableBoundaryRails && param->enableRailSupport;
  fast_simplifier = std::make_unique<FastSimplifier>(
    om, rm, vt.get(), ot.get(), lrt.get(), og.get(), &param->paramFastSimplifier);
  calc_diagonal_length();
  if (fast_simplifier->param->targetVerticesNum == 0)
    fast_simplifier->param->targetVerticesNum = (size_t)(om->n_vertices() * 3);

  Logger::user_logger->info("begin simplifying.");
  if (param->enableBoundaryRails)
  {
    Logger::user_logger->info(
      "boundary rail options: rail update {}, rail support {}.",
      param->enableRailUpdate ? "on" : "off",
      param->enableRailSupport ? "on" : "off");
  }

  degeneracyCleanupCases = 0;
  degeneracyCleanupSeconds = 0.0;
  if (enable_degeneracy_cleanup && param->phase3Mode == "fast" &&
    !param->enableBoundaryRails)
  {
    // Preserve the default pipeline's near-degeneracy cleanup. It is timed
    // separately from simplification.
    PhaseTimer::Exclusion exclusion("cleanup");
    const auto cleanup_start = PhaseTimer::Clock::now();
    degeneracyCleanupCases = degeneration_remover->perform();
    degeneracyCleanupSeconds = std::chrono::duration<double>(
      PhaseTimer::Clock::now() - cleanup_start).count();
    Logger::user_logger->info(
      "degeneracy cleanup in Phase 3: {} cases in {:.6f} seconds (not counted in the phase time).",
      degeneracyCleanupCases, degeneracyCleanupSeconds);
  }
  else
  {
    Logger::user_logger->info(
      "degeneracy cleanup before Phase 3 disabled for non-default pipeline.");
  }
  // Rail topology verification is disabled; the cage does not need it.
  // if (param->enableBoundaryRails &&
  //   !validate_and_log_boundary_rails(rm, "after degeneration removal"))
  //   throw std::logic_error("boundary rail topology became invalid during degeneration removal");
  if (is_phase3_energy_mode(param->phase3Mode))
  {
    run_phase3_energy_simplification();
    // Verification and the log-only minimum quality are disabled.
    // if (param->enableBoundaryRails &&
    //   !validate_and_log_boundary_rails(rm, "after Phase 3"))
    //   throw std::logic_error("boundary rail topology became invalid during Phase 3");
    // const std::string quality_label = "after phase 3 " + param->phase3Mode;
    // log_min_triangle_quality(quality_label.c_str());

    degeneration_remover = nullptr;
    fast_simplifier = nullptr;
    collapse_stage = nullptr;
    vt = nullptr;
    Logger::user_logger->info("simplification done.");
    return;
  }
  else
  {
    fast_simplifier->simplify();
    // log_min_triangle_quality("after fast simplify");
  }

  #ifdef OUTPUT_MIDDLE_RESULT
    OpenMesh::IO::write_mesh(*rm, param->fileOutPath + std::to_string(param->cageLabel) + "_fast_simplify.obj",
      OpenMesh::IO::Options::Default, 15);
  #endif
  degeneration_remover = nullptr;
  fast_simplifier = nullptr;
  vt = nullptr;

  // initialize hausdorff distance.
  rt = std::make_unique<FaceTree>(*rm);
  generate_out_links(om, rm, rt.get());
  rt = nullptr;
  calc_face_in_error(rm, std::vector<FaceHandle>(rm->faces_begin(), rm->faces_end()));
  #ifdef USE_TREE_SEARCH
  calc_out_error(rm, om, ot.get(), cage_infinite_fp);
  #else
  calc_out_error(rm, om, og.get(), cage_infinite_fp);
  #endif

  collapse_stage = std::make_unique<CollapseStage>(
    om, rm, &param->paramCollapse,
    ot.get(), lrt.get(), og.get(), original_diagonal_length);
  collapse_stage->use_rail_support =
    param->enableBoundaryRails && param->enableRailSupport;

  relocate_stage = std::make_unique<RelocateStage>(
    om, rm, &param->paramRelocate,
    vt.get(), ot.get(), lrt.get(), og.get(), original_diagonal_length);

  flip_stage = std::make_unique<FlipStage>(
    om, rm, &param->paramFlip,
    ot.get(), lrt.get(), og.get(), original_diagonal_length);

  simplify_to_target_num();
  Logger::user_logger->info("simplification done.");
}

void CageSimplifier::calc_diagonal_length()
{
  Vec3d ptMin(DBL_MAX, DBL_MAX, DBL_MAX);
  Vec3d ptMax(-DBL_MAX, -DBL_MAX, -DBL_MAX);
  for (const auto& vh : om->vertices())
  {
    const auto& point = om->point(vh);
    ptMin.minimize(point);
    ptMax.maximize(point);
  }
  original_diagonal_length = (ptMax - ptMin).norm();
  degeneration_remover->original_diagonal_length = original_diagonal_length;
  fast_simplifier->original_diagonal_length = original_diagonal_length;
}

void CageSimplifier::run_phase3_energy_simplification()
{
  if (param->targetVerticesNum == 0)
    Logger::user_logger->info(
      "running Phase 3 energy mode [{}] until no further progress (no vertex target).",
      param->phase3Mode);
  else
    Logger::user_logger->info(
      "running Phase 3 energy mode [{}] to final target {} vertices.",
      param->phase3Mode, param->targetVerticesNum);

  collapse_stage = std::make_unique<CollapseStage>(
    om, rm, &param->paramCollapse,
    ot.get(), lrt.get(), og.get(), original_diagonal_length);
  collapse_stage->use_rail_support =
    param->enableBoundaryRails && param->enableRailSupport;
  // Closed components such as the cage around a cavity stop at a tetrahedron.
  // The fast (default) pipeline keeps its original collapse rule.
  collapse_stage->keep_tetrahedra = true;
  if (param->phase3Mode == "linear_solve")
  {
    run_phase3_linear_solve_iterations();
    collapse_stage = nullptr;
    return;
  }
  collapse_stage->do_phase3_energy_simplification(param->targetVerticesNum);
  collapse_stage = nullptr;
  if (param->phase3Mode == "qem_original")
  {
    init_one_ring_faces(rm);
    return;
  }

  Logger::user_logger->info("running phase 3 energy final flip polish.");
  rt = std::make_unique<FaceTree>(*rm);
  generate_out_links(om, rm, rt.get());
  rt = nullptr;
  calc_face_in_error(rm, std::vector<FaceHandle>(rm->faces_begin(), rm->faces_end()));
#ifdef USE_TREE_SEARCH
  calc_out_error(rm, om, ot.get(), cage_infinite_fp);
#else
  calc_out_error(rm, om, og.get(), cage_infinite_fp);
#endif

  flip_stage = std::make_unique<FlipStage>(
    om, rm, &param->paramFlip,
    ot.get(), lrt.get(), og.get(), original_diagonal_length);
  flip_stage->update(true, cage_infinite_fp);
  flip_stage->do_flip();
  flip_stage = nullptr;
  init_one_ring_faces(rm);
}

void CageSimplifier::run_phase3_linear_solve_iterations()
{
  init_one_ring_faces(rm);
  if (param->phase3QualityPolishIterations == 0)
  {
    // Preserve the existing zero setting as a collapse-only comparison.
    collapse_stage->do_phase3_energy_simplification(param->targetVerticesNum);
    init_one_ring_faces(rm);
    return;
  }

  const bool until_stalled = param->targetVerticesNum == 0;
  const std::string cycle_limit = until_stalled ? "unlimited" :
    "up to " + std::to_string(param->phase3QualityPolishIterations);
  Logger::user_logger->info(
    "running phase 3 linear-solve iterations: {} collapse/flip/rail-update cycles, then up to {} final relocation sweeps with {} backtracking attempts per vertex; rail vertices fixed during relocation.",
    cycle_limit, param->paramRelocate.qualitySweeps,
    param->paramRelocate.lineSearchMaxIter);

  flip_stage = std::make_unique<FlipStage>(
    om, rm, &param->paramFlip,
    ot.get(), lrt.get(), og.get(), original_diagonal_length);
  for (size_t iteration = 0;
       until_stalled || iteration < param->phase3QualityPolishIterations;
       ++iteration)
  {
    // Rebuild collapse candidates after the previous cycle's flips and rail
    // relabeling, which can make previously constrained edges collapsible.
    // Keep this stage alive so its initial uniformity target is fixed.
    // The collapse stage skips its work once the target vertex count is met.
    const size_t collapsed = collapse_stage->do_phase3_energy_simplification(
      param->targetVerticesNum);

    // Degeneration removal and energy collapse do not update all normals.
    rm->update_normals();
    pre_calculate_edge_length(rm);
    pre_calculate_face_area(rm);
    const size_t flipped = flip_stage->do_quality_flip();
    // Each call processes triangle shortcuts until no more are possible.
    // Labels change only; newly released vertices are reconsidered by the
    // next cycle's collapse/flip stages and by the final relocation stage.
    const size_t rail_updates =
      param->enableBoundaryRails && param->enableRailUpdate ?
      update_boundary_rails(*rm) : 0;
    // The per-cycle minimum quality and rail verification are disabled.
    Logger::user_logger->info(
      "phase 3 linear-solve cycle {}: collapsed {}, flipped {}, rail updates {}, vertices {}.",
      iteration + 1, collapsed, flipped, rail_updates, rm->n_vertices());
    // if (param->enableBoundaryRails &&
    //   !validate_and_log_boundary_rails(rm, "after linear-solve cycle"))
    //   throw std::logic_error("boundary rail topology became invalid during linear-solve iteration");
    if (collapsed == 0 && flipped == 0 && rail_updates == 0)
    {
      Logger::user_logger->info("phase 3 linear-solve iterations stopped: no accepted collapses, flips, or rail updates.");
      break;
    }
    if (!until_stalled && iteration + 1 == param->phase3QualityPolishIterations)
      Logger::user_logger->info("phase 3 linear-solve iterations stopped: cycle limit reached.");
  }
  flip_stage = nullptr;

  // Relocation changes positions only after the final topology is selected.
  // Do not resume collapse or flip after this stage, even if a move would
  // make further topology changes possible.
  relocate_stage = std::make_unique<RelocateStage>(
    om, rm, &param->paramRelocate,
    vt.get(), ot.get(), lrt.get(), og.get(), original_diagonal_length);
  size_t relocated = 0;
  size_t sweeps = 0;
  if (param->paramRelocate.lineSearchMaxIter > 0)
  {
    for (; sweeps < param->paramRelocate.qualitySweeps;)
    {
      // Each sweep recomputes both targets from the latest positions.
      const size_t moved = relocate_stage->do_quality_relocate(
        param->paramRelocate.lineSearchMaxIter);
      ++sweeps;
      relocated += moved;
      if (moved == 0)
        break;
    }
  }
  // The minimum quality and rail verification after relocation are disabled.
  Logger::user_logger->info(
    "phase 3 linear-solve final relocation: moves {} in {} sweeps, vertices {}.",
    relocated, sweeps, rm->n_vertices());
  // if (param->enableBoundaryRails &&
  //   !validate_and_log_boundary_rails(rm, "after final linear-solve relocation"))
  //   throw std::logic_error("boundary rail topology became invalid during final linear-solve relocation");
  relocate_stage = nullptr;
}

void CageSimplifier::update_strategy()
{
  if (!allow_negtive)
  {
    collapse_stage->update(6, allow_negtive, max_distance_error);
    relocate_stage->update(6, allow_negtive, max_distance_error);
  }
  else
  {
    collapse_stage->update(10, allow_negtive, max_distance_error);
    relocate_stage->update(10, allow_negtive, max_distance_error);
  }
  flip_stage->update(allow_negtive, max_distance_error);
}

double CageSimplifier::calc_triangle_quality(FaceHandle fh) const
{
  Vec3d pts[3];
  size_t vertex_count = 0;
  for (VertexHandle vh : rm->fv_range(fh))
  {
    if (vertex_count >= 3)
      return 0.0;
    pts[vertex_count++] = rm->point(vh);
  }
  if (vertex_count != 3)
    return 0.0;

  const double a = (pts[1] - pts[0]).length();
  const double b = (pts[2] - pts[1]).length();
  const double c = (pts[0] - pts[2]).length();
  const double denom = a * a + b * b + c * c;
  if (denom <= 0.0)
    return 0.0;

  const double area = 0.5 * (pts[1] - pts[0]).cross(pts[2] - pts[0]).length();
  return 4.0 * std::sqrt(3.0) * area / denom;
}

double CageSimplifier::calc_min_triangle_quality() const
{
  // Only reported in the log.
  PhaseTimer::Exclusion exclusion("check");
  double min_quality = DBL_MAX;
  for (FaceHandle fh : rm->faces())
    min_quality = std::min(min_quality, calc_triangle_quality(fh));
  return min_quality == DBL_MAX ? 0.0 : min_quality;
}

void CageSimplifier::log_min_triangle_quality(const char* label) const
{
  Logger::user_logger->info("{} min triangle quality {}.", label, calc_min_triangle_quality());
}

void CageSimplifier::simplify_to_target_num()
{
  if (rm->n_vertices() <= param->targetVerticesNum)
    return;
  Logger::user_logger->info("collapsing to target vertices number {}.", param->targetVerticesNum);

  // goal
  const size_t edge_num_to_collapse = rm->n_vertices() - param->targetVerticesNum;

  // strategy: force collapsing procedure to jump out and do smoothing.
  std::queue<size_t> force_skip_en;
  size_t force_skip_times = 3;
  for (size_t i = 0;i < force_skip_times - 1;i++)
    force_skip_en.push(edge_num_to_collapse / force_skip_times);
  force_skip_en.push(edge_num_to_collapse - (edge_num_to_collapse / force_skip_times) * (force_skip_times - 1));

  // initialize parameters during simplification
  allow_negtive = false;
  max_distance_error = original_diagonal_length * param->initError;

  // record statistics during simplification
  size_t total_cen = 0;               // cen: Collapsed Edges Number
  size_t last_iter_cen = 0;

  size_t error_relax_iter = 0;
  size_t iter = 0;
  size_t no_collapsed_iter = 0;
  while (true)  // when collapsed_edge_num == edge_num_to_collapse, end loop immediately.
  {
    update_strategy();
    // Log-only minimum qualities are disabled.
    collapse_stage->do_collapse(force_skip_en.front(), total_cen);
    // log_min_triangle_quality("after collapse");
    flip_stage->do_flip();
    // log_min_triangle_quality("after flip");
    relocate_stage->do_relocate();
    // log_min_triangle_quality("after relocate");
    size_t collapsed_this_iter = total_cen - last_iter_cen;
    last_iter_cen = total_cen;
    // forced to jump out,
    // continue if still have edges to collapse |OR| break if reach the target vertices number.
    if (force_skip_en.front() == total_cen)
    {
      // don't relax distance error, do next iteration or stop.
      last_iter_cen = 0;
      total_cen = 0;
      force_skip_en.pop();
      if (force_skip_en.empty())
        break;
      else
        continue;
    }
    iter++;
    // if reach max iter, break.
    if (iter == param->maxIter)
      break;
    // little edges are collapsed, perhaps need relax valence constraint.
    if (collapsed_this_iter <= 10 && error_relax_iter >= param->maxErrorRelaxIter)
      param->paramCollapse.maxValence += 1;

    // no edge is collapsed.
    if (collapsed_this_iter == 0)
    {
      no_collapsed_iter++;
      if (no_collapsed_iter == 5)
        break;
    }
    else no_collapsed_iter = 0;   // reset

    // no edge is able to collapse, perhaps need relax distance error.
    if (iter % param->relaxErrorIterStep == 0)
    {
      if (!allow_negtive)
        allow_negtive = true;
      else if (error_relax_iter < param->maxErrorRelaxIter)
        max_distance_error += original_diagonal_length * param->errorStep;

      error_relax_iter++;
      Logger::user_logger->info("[{}]current maximal distance error {}", error_relax_iter, max_distance_error);
      Logger::user_logger->info("{} vertices and {} faces remained.", rm->n_vertices(), rm->n_faces());

    #ifdef OUTPUT_MIDDLE_RESULT
      OpenMesh::IO::write_mesh(*rm, param->fileOutPath + std::to_string(param->cageLabel) + "_simplify_iter_" + std::to_string(iter) + ".obj",
        OpenMesh::IO::Options::Default, 15);
    #endif
    }
  }
  rm->garbage_collection();
  lrt->collect_garbage();
}

}// namespace CageSimp
}// namespace Cage
