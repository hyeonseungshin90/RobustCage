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

  ParamTetrahedralizer paramTetrahedralizer;
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
  // Supported modes: "hausdorff", "length", "length_quality",
  // "post_edge_length", "post_face_area", "triangle_quality",
  // "post_edge_length_hard", "post_face_area_hard", "triangle_quality_hard".
  std::string priorityMode;
  // Supported length-quality submodes: "weighted", "relative_reject", "absolute_reject", "lexicographic".
  std::string lengthQualitySubMode;
  double lengthQualityWeight;
  double lengthQualityDegradationRatio;
  double lengthQualityMinQuality;

  // Collapse placement method:
  // "sampling" keeps the original candidate sampling strategy.
  // "optimization" uses local QEM/Newton-style energies to choose the post-collapse vertex position.
  std::string collapsePlacementMethod;
  // Phase 2 placement strategy:
  // "adaptive": QEM linear solve first, Newton only for selected hard cases.
  // "linear_only": QEM linear solve only during collapse.
  // "qem": use QEM plus quality/uniformity terms; no Newton.
  // "qem_no_collision": pure Garland-Heckbert QEM only; non-QEM energies and collision rejection are disabled.
  // "final_newton": QEM linear solve during collapse, then one fixed-topology Newton polish.
  // "newton_only": skip QEM linear solve and use Newton placement for every popped edge.
  // "quadratic_surrogate": use the 4x4 quadratic surrogate solve as a separate experimental path.
  std::string phase2PlacementStrategy;
  // Newton solver mode: "damped" or "trust_region".
  std::string newtonSolverMode;
  // Robustness mode:
  // "exact_reject": optimize first, then reject invalid final positions.
  // "exact_backtracking": use exact local checks during backtracking.
  // "ipc_line_search": additionally sample the step path with exact checks.
  std::string robustnessMode;
  // Curvature mode: "none" or "weighted_qem".
  std::string curvatureMode;
  // Uniformity mode: "none", "source", "global".
  std::string uniformityMode;
  // When true, phase2PlacementStrategy "qem" rejects the edge outright (no
  // collapse) if the raw QEM linear-solve point fails the hard validity
  // checks (collision/degenerate/wrinkle), instead of backtracking toward a
  // nearby fallback position. Default false keeps the existing backtracking
  // behavior; this is an opt-in alternative.
  bool phase2QemRejectOnLinearCollision;

  size_t newtonMaxIter;
  double newtonGradTol;
  double newtonStepTol;
  double newtonFiniteDiffScale;
  double trustRegionRadiusScale;
  size_t lineSearchMaxIter;
  size_t lineSearchCcdSamples;
  double phase2NewtonQualityThreshold;
  double phase2NewtonResidualThreshold;
  double phase2NewtonResidualGrowth;
  size_t phase2NewtonFinalRefineCollapses;

  double qemWeight;
  double positionFidelityWeight;
  double triangleQualityWeight;
  double uniformityWeight;

  boost::json::object serialize()const
  {
    boost::json::object jo;
    jo["maxValence"] = maxValence;
    jo["priorityMode"] = priorityMode;
    jo["lengthQualitySubMode"] = lengthQualitySubMode;
    jo["lengthQualityWeight"] = lengthQualityWeight;
    jo["lengthQualityDegradationRatio"] = lengthQualityDegradationRatio;
    jo["lengthQualityMinQuality"] = lengthQualityMinQuality;
    jo["collapsePlacementMethod"] = collapsePlacementMethod;
    jo["phase2PlacementStrategy"] = phase2PlacementStrategy;
    jo["newtonSolverMode"] = newtonSolverMode;
    jo["robustnessMode"] = robustnessMode;
    jo["curvatureMode"] = curvatureMode;
    jo["uniformityMode"] = uniformityMode;
    jo["phase2QemRejectOnLinearCollision"] = phase2QemRejectOnLinearCollision;
    jo["newtonMaxIter"] = newtonMaxIter;
    jo["newtonGradTol"] = newtonGradTol;
    jo["newtonStepTol"] = newtonStepTol;
    jo["newtonFiniteDiffScale"] = newtonFiniteDiffScale;
    jo["trustRegionRadiusScale"] = trustRegionRadiusScale;
    jo["lineSearchMaxIter"] = lineSearchMaxIter;
    jo["lineSearchCcdSamples"] = lineSearchCcdSamples;
    jo["phase2NewtonQualityThreshold"] = phase2NewtonQualityThreshold;
    jo["phase2NewtonResidualThreshold"] = phase2NewtonResidualThreshold;
    jo["phase2NewtonResidualGrowth"] = phase2NewtonResidualGrowth;
    jo["phase2NewtonFinalRefineCollapses"] = phase2NewtonFinalRefineCollapses;
    jo["qemWeight"] = qemWeight;
    jo["positionFidelityWeight"] = positionFidelityWeight;
    jo["triangleQualityWeight"] = triangleQualityWeight;
    jo["uniformityWeight"] = uniformityWeight;
    return jo;
  }
  void deserialize(const boost::json::object& jo)
  {
    maxValence = jo.at("maxValence").as_int64();
    auto priority_mode_it = jo.find("priorityMode");
    if (priority_mode_it != jo.end())
      priorityMode = std::string(priority_mode_it->value().as_string().c_str());
    else
      priorityMode = "hausdorff";

    auto lq_submode_it = jo.find("lengthQualitySubMode");
    if (lq_submode_it != jo.end())
      lengthQualitySubMode = std::string(lq_submode_it->value().as_string().c_str());
    else
      lengthQualitySubMode = "weighted";

    auto lq_weight_it = jo.find("lengthQualityWeight");
    lengthQualityWeight = lq_weight_it != jo.end() ? lq_weight_it->value().as_double() : 5.0;

    auto lq_ratio_it = jo.find("lengthQualityDegradationRatio");
    lengthQualityDegradationRatio = lq_ratio_it != jo.end() ? lq_ratio_it->value().as_double() : 0.5;

    auto lq_min_quality_it = jo.find("lengthQualityMinQuality");
    lengthQualityMinQuality = lq_min_quality_it != jo.end() ? lq_min_quality_it->value().as_double() : 0.05;

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
      phase2_placement_strategy_it != jo.end() ? std::string(phase2_placement_strategy_it->value().as_string().c_str()) : "adaptive";

    auto solver_mode_it = jo.find("newtonSolverMode");
    newtonSolverMode = solver_mode_it != jo.end() ? std::string(solver_mode_it->value().as_string().c_str()) : "damped";

    auto robustness_mode_it = jo.find("robustnessMode");
    robustnessMode = robustness_mode_it != jo.end() ? std::string(robustness_mode_it->value().as_string().c_str()) : "exact_backtracking";

    auto curvature_mode_it = jo.find("curvatureMode");
    curvatureMode = curvature_mode_it != jo.end() ? std::string(curvature_mode_it->value().as_string().c_str()) : "none";

    auto uniformity_mode_it = jo.find("uniformityMode");
    uniformityMode = uniformity_mode_it != jo.end() ? std::string(uniformity_mode_it->value().as_string().c_str()) : "none";

    auto phase2_qem_reject_on_linear_collision_it = jo.find("phase2QemRejectOnLinearCollision");
    phase2QemRejectOnLinearCollision =
      phase2_qem_reject_on_linear_collision_it != jo.end() ? phase2_qem_reject_on_linear_collision_it->value().as_bool() : false;

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

    auto line_search_ccd_samples_it = jo.find("lineSearchCcdSamples");
    lineSearchCcdSamples = line_search_ccd_samples_it != jo.end() ? line_search_ccd_samples_it->value().as_int64() : 4;

    auto phase2_quality_threshold_it = jo.find("phase2NewtonQualityThreshold");
    phase2NewtonQualityThreshold =
      phase2_quality_threshold_it != jo.end() ? phase2_quality_threshold_it->value().as_double() : 0.12;

    auto phase2_residual_threshold_it = jo.find("phase2NewtonResidualThreshold");
    phase2NewtonResidualThreshold =
      phase2_residual_threshold_it != jo.end() ? phase2_residual_threshold_it->value().as_double() : 0.25;

    auto phase2_residual_growth_it = jo.find("phase2NewtonResidualGrowth");
    phase2NewtonResidualGrowth =
      phase2_residual_growth_it != jo.end() ? phase2_residual_growth_it->value().as_double() : 1.5;

    auto phase2_final_refine_it = jo.find("phase2NewtonFinalRefineCollapses");
    phase2NewtonFinalRefineCollapses =
      phase2_final_refine_it != jo.end() ? phase2_final_refine_it->value().as_int64() : 25;

    auto qem_weight_it = jo.find("qemWeight");
    qemWeight = qem_weight_it != jo.end() ? qem_weight_it->value().as_double() : 1.0;

    auto position_fidelity_weight_it = jo.find("positionFidelityWeight");
    positionFidelityWeight =
      position_fidelity_weight_it != jo.end() ? position_fidelity_weight_it->value().as_double() : 1.0;

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
  // "newton" uses adaptive QEM/Newton energy collapses, "linear_only"
  // uses the linear Phase 2 solve without Newton refinement, and "qem"
  // uses QEM plus quality/uniformity terms with hard validity rejection.
  // "qem_no_collision" uses pure Garland-Heckbert QEM cost/placement only.
  std::string phase2Mode;
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
    collapse.priorityMode = "hausdorff";
    collapse.lengthQualitySubMode = "weighted";
    collapse.lengthQualityWeight = 5.0;
    collapse.lengthQualityDegradationRatio = 0.5;
    collapse.lengthQualityMinQuality = 0.1;
    collapse.collapsePlacementMethod = "sampling";
    collapse.phase2PlacementStrategy = "adaptive";
    collapse.newtonSolverMode = "damped";
    collapse.robustnessMode = "exact_backtracking";
    collapse.curvatureMode = "none";
    collapse.uniformityMode = "none";
    collapse.phase2QemRejectOnLinearCollision = false;
    collapse.newtonMaxIter = 4;
    collapse.newtonGradTol = 1e-8;
    collapse.newtonStepTol = 1e-8;
    collapse.newtonFiniteDiffScale = 1e-4;
    collapse.trustRegionRadiusScale = 0.25;
    collapse.lineSearchMaxIter = 6;
    collapse.lineSearchCcdSamples = 4;
    collapse.phase2NewtonQualityThreshold = 0.12;
    collapse.phase2NewtonResidualThreshold = 0.25;
    collapse.phase2NewtonResidualGrowth = 1.5;
    collapse.phase2NewtonFinalRefineCollapses = 25;
    collapse.qemWeight = 1.0;
    collapse.positionFidelityWeight = 1.0;
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
    jo["paramCageSimplifier"] = paramCageSimplifier.serialize();
    return jo;
  }
  void deserialize(const boost::json::object& jo)
  {
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
