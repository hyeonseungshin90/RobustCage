#pragma once
#include <vector>
#include <string>
#include "boost/json.hpp"

namespace Cage
{
struct ParamLatticePointsGenerator
{
  /// minimal point's distance to surface mesh, nearer points will be ingored.
  double minDistanceFactor;

  /***** used for "Uniform" *****/

  /// Scale mesh's bounding box to a larger one.
  double bboxScale;
  /// Sample n^3 lattice points in bounding box uniformly.
  size_t pointNumAlongAxis;

  /***** used for "Offset" *****/

  /// good approximation: 4, high efficiency: 16.
  double areaThresholdRate;
  /// Multiplier for offset point distance from the surface.
  double offsetLengthScale;
};

struct ParamTetrahedralizer
{
  ParamLatticePointsGenerator paramLatticePointsGenerator;
};

struct ParamCageInitializer
{
  std::string fileOutPath;
  std::string fileName;

  // Phase 1 construction mode:
  // "subdivision" preserves the original two-round 1-to-12 refinement;
  // "topological_offset" uses simplicial embedding and offset insertion.
  std::string phase1Mode = "subdivision";

  ParamTetrahedralizer paramTetrahedralizer;

  boost::json::object serialize()const
  {
    boost::json::object jo;
    jo["phase1Mode"] = phase1Mode;
    return jo;
  }
  void deserialize(const boost::json::object& jo)
  {
    auto phase1_mode_it = jo.find("phase1Mode");
    if (phase1_mode_it != jo.end())
      phase1Mode = std::string(phase1_mode_it->value().as_string().c_str());
  }
};

struct ParamFastSimplifier
{
  std::string fileOutPath;
  std::string fileName;

  // target
  size_t targetVerticesNum;

  // parameters
  size_t collapseIter;
  size_t splitIter;
  size_t equalizeValenceIter;
  size_t smoothIter;
  double enlargeTargetLengthRatio;
};

struct ParamCollapseStage
{
  // constraints
  size_t maxValence;

  // Collapse placement method:
  // "sampling" keeps the original candidate sampling strategy.
  // "optimization" uses local QEM/Newton-style energies to choose the post-collapse vertex position.
  std::string collapsePlacementMethod;
  // Phase 2 placement strategy:
  // "linear_solve": use the QEM-based linear system with quality terms and no Newton.
  // "newton_solve": use Newton placement for every popped edge.
  // "qem_original": pure Garland-Heckbert QEM; non-QEM energies and collision rejection are disabled.
  std::string phase2PlacementStrategy;
  // Newton solver mode: "damped" or "trust_region".
  std::string newtonSolverMode;
  // Robustness mode:
  // "exact_reject": optimize first, then reject invalid final positions.
  // "exact_backtracking": use exact local checks during backtracking.
  std::string robustnessMode;
  // Curvature mode: "none" or "weighted_qem".
  std::string curvatureMode;
  // Uniformity mode: "none", "source", "global".
  std::string uniformityMode;
  // When true, phase2PlacementStrategy "linear_solve" rejects the edge outright (no
  // collapse) if the raw QEM linear-solve point fails the hard validity
  // checks (collision/degenerate/wrinkle), instead of backtracking toward a
  // nearby fallback position. Default false keeps the existing backtracking
  // behavior; this is an opt-in alternative.
  bool phase2LinearSolveCollisionReject;

  size_t newtonMaxIter;
  double newtonGradTol;
  double newtonStepTol;
  double newtonFiniteDiffScale;
  double trustRegionRadiusScale;
  size_t lineSearchMaxIter;
  double qemWeight;
  double triangleQualityWeight;
  double uniformityWeight;

  boost::json::object serialize()const
  {
    boost::json::object jo;
    jo["maxValence"] = maxValence;
    jo["collapsePlacementMethod"] = collapsePlacementMethod;
    jo["phase2PlacementStrategy"] = phase2PlacementStrategy;
    jo["newtonSolverMode"] = newtonSolverMode;
    jo["robustnessMode"] = robustnessMode;
    jo["curvatureMode"] = curvatureMode;
    jo["uniformityMode"] = uniformityMode;
    jo["phase2LinearSolveCollisionReject"] = phase2LinearSolveCollisionReject;
    jo["newtonMaxIter"] = newtonMaxIter;
    jo["newtonGradTol"] = newtonGradTol;
    jo["newtonStepTol"] = newtonStepTol;
    jo["newtonFiniteDiffScale"] = newtonFiniteDiffScale;
    jo["trustRegionRadiusScale"] = trustRegionRadiusScale;
    jo["lineSearchMaxIter"] = lineSearchMaxIter;
    jo["qemWeight"] = qemWeight;
    jo["triangleQualityWeight"] = triangleQualityWeight;
    jo["uniformityWeight"] = uniformityWeight;
    return jo;
  }
  void deserialize(const boost::json::object& jo)
  {
    maxValence = jo.at("maxValence").as_int64();
    auto collapse_placement_method_it = jo.find("collapsePlacementMethod");
    if (collapse_placement_method_it == jo.end())
      collapse_placement_method_it = jo.find("placementMode");
    collapsePlacementMethod =
      collapse_placement_method_it != jo.end() ? std::string(collapse_placement_method_it->value().as_string().c_str()) : "sampling";
    if (collapsePlacementMethod == "newton" || collapsePlacementMethod == "energy")
      collapsePlacementMethod = "optimization";

    auto phase2_placement_strategy_it = jo.find("phase2PlacementStrategy");
    if (phase2_placement_strategy_it == jo.end())
      phase2_placement_strategy_it = jo.find("phase2PlacementMode");
    phase2PlacementStrategy =
      phase2_placement_strategy_it != jo.end() ? std::string(phase2_placement_strategy_it->value().as_string().c_str()) : "linear_solve";

    auto solver_mode_it = jo.find("newtonSolverMode");
    newtonSolverMode = solver_mode_it != jo.end() ? std::string(solver_mode_it->value().as_string().c_str()) : "damped";

    auto robustness_mode_it = jo.find("robustnessMode");
    robustnessMode = robustness_mode_it != jo.end() ? std::string(robustness_mode_it->value().as_string().c_str()) : "exact_backtracking";

    auto curvature_mode_it = jo.find("curvatureMode");
    curvatureMode = curvature_mode_it != jo.end() ? std::string(curvature_mode_it->value().as_string().c_str()) : "none";

    auto uniformity_mode_it = jo.find("uniformityMode");
    uniformityMode = uniformity_mode_it != jo.end() ? std::string(uniformity_mode_it->value().as_string().c_str()) : "none";

    auto phase2_linear_solve_collision_reject_it = jo.find("phase2LinearSolveCollisionReject");
    phase2LinearSolveCollisionReject =
      phase2_linear_solve_collision_reject_it != jo.end() ? phase2_linear_solve_collision_reject_it->value().as_bool() : false;

    auto newton_max_iter_it = jo.find("newtonMaxIter");
    newtonMaxIter = newton_max_iter_it != jo.end() ? newton_max_iter_it->value().as_int64() : 4;

    auto newton_grad_tol_it = jo.find("newtonGradTol");
    newtonGradTol = newton_grad_tol_it != jo.end() ? newton_grad_tol_it->value().as_double() : 1e-8;

    auto newton_step_tol_it = jo.find("newtonStepTol");
    newtonStepTol = newton_step_tol_it != jo.end() ? newton_step_tol_it->value().as_double() : 1e-8;

    auto finite_diff_scale_it = jo.find("newtonFiniteDiffScale");
    newtonFiniteDiffScale = finite_diff_scale_it != jo.end() ? finite_diff_scale_it->value().as_double() : 1e-4;

    auto trust_radius_it = jo.find("trustRegionRadiusScale");
    trustRegionRadiusScale = trust_radius_it != jo.end() ? trust_radius_it->value().as_double() : 0.25;

    auto line_search_max_iter_it = jo.find("lineSearchMaxIter");
    lineSearchMaxIter = line_search_max_iter_it != jo.end() ? line_search_max_iter_it->value().as_int64() : 6;

    auto qem_weight_it = jo.find("qemWeight");
    qemWeight = qem_weight_it != jo.end() ? qem_weight_it->value().as_double() : 1.0;

    auto triangle_quality_weight_it = jo.find("triangleQualityWeight");
    triangleQualityWeight = triangle_quality_weight_it != jo.end() ? triangle_quality_weight_it->value().as_double() : 2.0;

    auto uniformity_weight_it = jo.find("uniformityWeight");
    uniformityWeight = uniformity_weight_it != jo.end() ? uniformity_weight_it->value().as_double() : 1.0;
  }
};

struct ParamRelocateStage
{
  // simplification parameter
  size_t smoothIter;
  // Supported modes: "hausdorff", "triangle_quality_hard".
  std::string priorityMode;

  boost::json::object serialize()const
  {
    boost::json::object jo;
    jo["smoothIter"] = smoothIter;
    jo["priorityMode"] = priorityMode;
    return jo;
  }
  void deserialize(const boost::json::object& jo)
  {
    auto smooth_iter_it = jo.find("smoothIter");
    if (smooth_iter_it != jo.end())
      smoothIter = smooth_iter_it->value().as_int64();

    auto priority_mode_it = jo.find("priorityMode");
    if (priority_mode_it != jo.end())
      priorityMode = std::string(priority_mode_it->value().as_string().c_str());
    else
      priorityMode = "hausdorff";
  }
};

struct ParamFlipStage
{
  // constraints
  size_t maxValence;
  // Supported modes: "valence", "triangle_quality_hard".
  std::string priorityMode;
  // When true, quality-priority flips must also move all four incident
  // vertices toward regular valence 6. Used by the Newton Phase 2 polish.
  bool requireRegularValence;

  boost::json::object serialize()const
  {
    boost::json::object jo;
    jo["maxValence"] = maxValence;
    jo["priorityMode"] = priorityMode;
    jo["requireRegularValence"] = requireRegularValence;
    return jo;
  }
  void deserialize(const boost::json::object& jo)
  {
    maxValence = jo.at("maxValence").as_int64();
    auto priority_mode_it = jo.find("priorityMode");
    if (priority_mode_it != jo.end())
      priorityMode = std::string(priority_mode_it->value().as_string().c_str());
    else
      priorityMode = "valence";

    auto require_regular_valence_it = jo.find("requireRegularValence");
    requireRegularValence =
      require_regular_valence_it != jo.end() ? require_regular_valence_it->value().as_bool() : false;
  }
};

struct ParamCageSimplifier
{
  // target
  size_t targetVerticesNum;
  // Phase 2 simplification mode: "fast" keeps the original FastSimplifier,
  // "linear_solve" uses the QEM-based linear system, "newton_solve" uses
  // Newton placement for every collapse, and "qem_original" uses pure
  // Garland-Heckbert QEM cost/placement only.
  std::string phase2Mode;
  // Build and preserve source-boundary rails on the initial cage.  This is
  // deliberately opt-in so the original/default pipeline is unchanged.
  bool enableBoundaryRails = false;
  // iterations
  size_t maxIter;
  // distance error control
  // error relax iter->Hausdorff Distance: 0->no negtive, 1->initError, 2->initError+errorStep, 3->initError+errorStep*2, ...
  size_t relaxErrorIterStep;
  size_t maxErrorRelaxIter;
  double initError;
  double errorStep;

  // file output
  size_t cageLabel;
  std::string fileOutPath;
  std::string fileName;

  ParamFastSimplifier paramFastSimplifier;
  ParamCollapseStage paramCollapse;
  ParamRelocateStage paramRelocate;
  ParamFlipStage paramFlip;

  boost::json::object serialize()const
  {
    boost::json::object jo;
    jo["maxIter"] = maxIter;
    jo["phase2Mode"] = phase2Mode;
    jo["enableBoundaryRails"] = enableBoundaryRails;
    jo["relaxErrorIterStep"] = relaxErrorIterStep;
    jo["maxErrorRelaxIter"] = maxErrorRelaxIter;
    jo["initError"] = initError;
    jo["errorStep"] = errorStep;
    jo["paramCollapse"] = paramCollapse.serialize();
    jo["paramRelocate"] = paramRelocate.serialize();
    jo["paramFlip"] = paramFlip.serialize();
    return jo;
  }
  void deserialize(const boost::json::object& jo)
  {
    maxIter = jo.at("maxIter").as_int64();
    auto phase2_mode_it = jo.find("phase2Mode");
    phase2Mode = phase2_mode_it != jo.end() ? std::string(phase2_mode_it->value().as_string().c_str()) : "fast";
    auto boundary_rails_it = jo.find("enableBoundaryRails");
    enableBoundaryRails = boundary_rails_it != jo.end() ?
      boundary_rails_it->value().as_bool() : false;
    relaxErrorIterStep = jo.at("relaxErrorIterStep").as_int64();
    maxErrorRelaxIter = jo.at("maxErrorRelaxIter").as_int64();
    initError = jo.at("initError").as_double();
    errorStep = jo.at("errorStep").as_double();
    paramCollapse.deserialize(jo.at("paramCollapse").as_object());
    auto relocate_it = jo.find("paramRelocate");
    if (relocate_it != jo.end())
      paramRelocate.deserialize(relocate_it->value().as_object());
    auto flip_it = jo.find("paramFlip");
    if (flip_it != jo.end())
      paramFlip.deserialize(flip_it->value().as_object());
  }
};

struct ParamCageGenerator
{
  ParamCageInitializer paramCageInitializer;
  ParamCageSimplifier paramCageSimplifier;

  ParamCageGenerator()
  {
    auto& Lpg = paramCageInitializer.paramTetrahedralizer.paramLatticePointsGenerator;
    Lpg.minDistanceFactor = 0.01;
    Lpg.bboxScale = 1.5;
    Lpg.pointNumAlongAxis = 10;
    Lpg.areaThresholdRate = 16;
    Lpg.offsetLengthScale = 1.0;

    auto& simplifier = paramCageSimplifier;
    simplifier.phase2Mode = "fast";
    simplifier.enableBoundaryRails = false;
    simplifier.maxIter = 30;
    simplifier.relaxErrorIterStep = 5;
    simplifier.maxErrorRelaxIter = 4;
    simplifier.initError = 0.005;
    simplifier.errorStep = 0.005;
    simplifier.targetVerticesNum = 0;

    auto& fast = paramCageSimplifier.paramFastSimplifier;
    fast.collapseIter = 3;
    fast.splitIter = 1;
    fast.equalizeValenceIter = 1;
    fast.smoothIter = 3;
    fast.enlargeTargetLengthRatio = 1.1;
    fast.targetVerticesNum = 0;

    auto& collapse = paramCageSimplifier.paramCollapse;
    collapse.maxValence = 8;
    collapse.collapsePlacementMethod = "sampling";
    collapse.phase2PlacementStrategy = "linear_solve";
    collapse.newtonSolverMode = "damped";
    collapse.robustnessMode = "exact_backtracking";
    collapse.curvatureMode = "none";
    collapse.uniformityMode = "none";
    collapse.phase2LinearSolveCollisionReject = false;
    collapse.newtonMaxIter = 4;
    collapse.newtonGradTol = 1e-8;
    collapse.newtonStepTol = 1e-8;
    collapse.newtonFiniteDiffScale = 1e-4;
    collapse.trustRegionRadiusScale = 0.25;
    collapse.lineSearchMaxIter = 6;
    collapse.qemWeight = 1.0;
    collapse.triangleQualityWeight = 2.0;
    collapse.uniformityWeight = 1.0;

    auto& relocate = paramCageSimplifier.paramRelocate;
    relocate.smoothIter = 3;
    relocate.priorityMode = "hausdorff";

    auto& flip = paramCageSimplifier.paramFlip;
    flip.maxValence = 8;
    flip.priorityMode = "valence";
    flip.requireRegularValence = false;
  }

  void setOutputPath(const std::string& outDir, const std::string& outFile)
  {
    paramCageInitializer.fileOutPath = outDir;
    paramCageInitializer.fileName = outFile;
    paramCageSimplifier.fileOutPath = outDir;
    paramCageSimplifier.fileName = outFile;
    paramCageSimplifier.paramFastSimplifier.fileOutPath = outDir;
    paramCageSimplifier.paramFastSimplifier.fileName = outFile;
  }

  void setCageLabel(const size_t label)
  {
    paramCageSimplifier.cageLabel = label;
  }
  void setTargetNumber(const size_t target)
  {
    paramCageSimplifier.targetVerticesNum = target;
  }
  void setFastTargetNumber(const size_t target)
  {
    paramCageSimplifier.paramFastSimplifier.targetVerticesNum = target;
  }

  boost::json::object serialize()const
  {
    boost::json::object jo;
    jo["paramCageInitializer"] = paramCageInitializer.serialize();
    jo["paramCageSimplifier"] = paramCageSimplifier.serialize();
    return jo;
  }
  void deserialize(const boost::json::object& jo)
  {
    auto initializer_it = jo.find("paramCageInitializer");
    if (initializer_it != jo.end())
      paramCageInitializer.deserialize(initializer_it->value().as_object());
    paramCageSimplifier.deserialize(jo.at("paramCageSimplifier").as_object());
  }
};

inline void tag_invoke(boost::json::value_from_tag, boost::json::value& jv, const ParamCageGenerator& param)
{
  auto& jo = jv.emplace_object();
  jo["paramCageGenerator"] = param.serialize();
}

inline ParamCageGenerator tag_invoke(boost::json::value_to_tag<ParamCageGenerator>, boost::json::value& jv)
{
  ParamCageGenerator param;
  param.deserialize(jv.as_object());
  return param;
}
}
