#pragma once

#include <cfloat>
#include <queue>

#include "Config.hh"
#include "CageSimplifier/Topo/TopoOperations.h"
#include "CageSimplifier/Topo/EdgeCollapser.h"
#include "CageSimplifier/Topo/VertexRelocater.h"
#include "CageSimplifier/Geom/GeometryCheck.h"

namespace Cage
{
namespace CageSimp
{

class CollapseStage
{
public:
  // input
  ParamCollapseStage* param;
  SMeshT* om;
  SMeshT* rm;
  DFaceTree* ot;
  LightDFaceTree* lrt;
  FaceGrid* og;

  double original_diagonal_length;
  size_t candidate_points_size;
  bool allow_negtive;
  double max_distance_error;
public:
  CollapseStage(
    SMeshT* original, SMeshT* cage, ParamCollapseStage* p,
    DFaceTree* original_tree, LightDFaceTree* remeshing_tree,
    FaceGrid* original_grid,
    double _original_diagonal_length);

  void do_collapse(size_t edge_num_to_collapse, size_t& total_collapsed_edge_num);
  // Replaces legacy Phase 2 by repeatedly collapsing edges until the target
  // vertex count is reached. Each collapse chooses its new vertex position
  // before the topology change is committed.
  void do_phase2_optimization_simplification(size_t target_vertices_num);
  void update(size_t _candidate_points_size, bool _allow_negtive, double _max_distance_error);
private:
  double avg_edge_length;
  double avg_original_edge_length;
  std::vector<size_t> update_states;

  struct CollapseEdgeReward
  {
    EdgeHandle eh;
    size_t state;
    double reward;
    double secondary_reward = 0.0;
    Vec3d new_point;

    CollapseEdgeReward() = default;
    CollapseEdgeReward(EdgeHandle _eh, size_t _state, double _reward, const Vec3d& _new_point, double _secondary_reward = 0.0) :
      eh(_eh), state(_state), reward(_reward), secondary_reward(_secondary_reward), new_point(_new_point)
    {}

    bool operator<(const CollapseEdgeReward& rhs)const
    {
      if (reward != rhs.reward)
        return reward < rhs.reward;
      return secondary_reward < rhs.secondary_reward;
    }
  };
  typedef std::priority_queue<CollapseEdgeReward> CollapseEdgeRewardQueue;
  CollapseEdgeRewardQueue edges_to_collapse;

  struct Phase2EdgeReward
  {
    EdgeHandle eh;
    size_t state;
    double reward;
    Vec3d initial_point;

    Phase2EdgeReward() = default;
    Phase2EdgeReward(EdgeHandle _eh, size_t _state, double _reward, const Vec3d& _initial_point) :
      eh(_eh), state(_state), reward(_reward), initial_point(_initial_point)
    {}

    bool operator<(const Phase2EdgeReward& rhs)const
    {
      return reward < rhs.reward;
    }
  };
  typedef std::priority_queue<Phase2EdgeReward> Phase2EdgeRewardQueue;
  Phase2EdgeRewardQueue phase2_edges_to_collapse;

  enum class EdgeSide
  {
    Inside,
    Outside,
    OnSurface,
    Unknown
  };

  struct EdgeSideStats
  {
    size_t inside = 0;
    size_t outside = 0;
    size_t on_surface = 0;
    size_t unknown = 0;

    size_t total() const { return inside + outside + on_surface + unknown; }
  };

  struct QEMPlane
  {
    Vec3d normal;
    double offset = 0.0;
    double weight = 1.0;
  };

  struct CollapseFanEdge
  {
    Vec3d from;
    Vec3d to;
    Vec3d apex;
  };

  // Local geometric data for evaluating one edge-collapse placement x.
  // This object does not run Newton; it is shared by the QEM proxy solve,
  // fallback candidate checks, and optional Newton refinement.
  struct Phase2PlacementContext
  {
    EdgeHandle edge;
    std::vector<HalfedgeHandle> halfedges;
    std::vector<CollapseFanEdge> fan_edges;
    std::vector<Vec3d> neighbor_points;
    std::vector<QEMPlane> qem_planes;
    std::set<int> ignored_faces;
    Vec3d start_point;
    Vec3d midpoint;
    double local_scale = 1.0;
    double edge_curvature = 0.0;
    double global_target_length = 1.0;
  };

  struct NewtonDerivatives
  {
    double energy = DBL_MAX;
    Vec3d gradient;
    double hessian[3][3] = {};
  };

  struct Phase2PlacementDecision
  {
    Vec3d point;
    double priority_energy = DBL_MAX;
    double min_quality = 1.0;
    double nonlinear_residual = 0.0;
    bool attempted_newton = false;
    bool accepted_newton = false;
    bool newton_failed = false;
    bool used_fallback = false;
  };

  std::vector<Vec3d> generate_candidate_points_for_collapse(EdgeHandle e, EdgeCollapser& edge_collapser);
  bool find_collapse_hausdorff_deviation(
    EdgeHandle eh, double& local_hd_before, double& local_hd_after, Vec3d& new_point, double& priority_score);
  bool refine_collapse_placement_with_newton(
    EdgeHandle eh, EdgeCollapser& edge_collapser, Vec3d& new_point, double& energy,
    const Vec3d* initial_point = nullptr);
  void initialize_collapse_edges_reward();
  void update_after_collapsing(VertexHandle collapsed_center);
  bool try_enqueue_collapse_candidate(
    EdgeHandle eh, size_t state, double local_hd_before, double local_hd_after, const Vec3d& new_point, double priority_score);
  bool is_optimization_collapse_placement_method() const;
  bool is_phase2_adaptive_strategy() const;
  bool is_phase2_linear_only_strategy() const;
  bool is_phase2_final_newton_strategy() const;
  bool is_phase2_newton_only_strategy() const;
  bool is_phase2_quadratic_surrogate_strategy() const;
  bool is_trust_region_solver_mode() const;
  bool is_exact_reject_robustness_mode() const;
  bool is_exact_backtracking_robustness_mode() const;
  bool is_ipc_line_search_robustness_mode() const;
  bool is_length_priority_mode() const;
  bool is_length_quality_priority_mode() const;
  bool is_length_quality_weighted_submode() const;
  bool is_length_quality_relative_reject_submode() const;
  bool is_length_quality_absolute_reject_submode() const;
  bool is_length_quality_lexicographic_submode() const;
  bool is_post_edge_length_priority_mode() const;
  bool is_post_face_area_priority_mode() const;
  bool is_triangle_quality_priority_mode() const;
  bool is_post_metric_priority_mode() const;
  bool is_hard_post_metric_priority_mode() const;
  bool is_length_quality_allowed(double pre_quality, double post_quality) const;
  double calc_length_quality_reward(double normalized_length_score, double quality_delta) const;
  double calc_pre_collapse_metric(EdgeHandle eh) const;
  double calc_post_collapse_metric_score(double pre_metric, SMeshT* local_mesh) const;
  double calc_max_edge_length(SMeshT* mesh) const;
  double calc_max_edge_length(SMeshT* mesh, const std::set<FaceHandle>& faces) const;
  double calc_max_face_area(SMeshT* mesh) const;
  double calc_max_face_area(SMeshT* mesh, const std::set<FaceHandle>& faces) const;
  double calc_triangle_quality(SMeshT* mesh, FaceHandle fh) const;
  double calc_min_triangle_quality(SMeshT* mesh) const;
  double calc_min_triangle_quality(SMeshT* mesh, const std::set<FaceHandle>& faces) const;
  Phase2PlacementContext make_phase2_placement_context(EdgeHandle eh, EdgeCollapser& edge_collapser) const;
  Phase2PlacementContext make_phase2_vertex_relocation_context(VertexHandle vh, VertexRelocater& vertex_relocater) const;
  NewtonDerivatives finite_difference_newton_derivatives(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  NewtonDerivatives approximate_newton_derivatives(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  NewtonDerivatives autodiff_newton_derivatives(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double evaluate_newton_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double evaluate_qem_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double evaluate_original_barrier_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double evaluate_self_barrier_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double evaluate_position_fidelity_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double evaluate_curvature_normal_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double evaluate_triangle_quality_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double evaluate_uniformity_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double evaluate_phase2_proxy_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double evaluate_phase2_refinement_residual(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double calc_phase2_fan_min_quality(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double source_uniformity_target_length(const Phase2PlacementContext& ctx, const Vec3d& p) const;
  double ipc_barrier(double distance, double dhat) const;
  double source_size_at(const Vec3d& p) const;
  bool closest_original_point(const Vec3d& p, Vec3d& closest) const;
  bool collapse_target_valid(EdgeCollapser& edge_collapser, const Vec3d& x) const;
  bool sampled_path_valid(EdgeCollapser& edge_collapser, const Vec3d& from, const Vec3d& to) const;
  bool phase2_placement_satisfies_hard_constraints(
    const Phase2PlacementContext& ctx, EdgeCollapser& edge_collapser, const Vec3d& x) const;
  bool should_refine_phase2_placement_with_newton(
    const Phase2PlacementContext& ctx, const Vec3d& x,
    size_t remaining_vertices, size_t target_vertices_num,
    double min_quality, double nonlinear_residual) const;
  bool choose_phase2_collapse_placement(
    EdgeHandle eh, EdgeCollapser& edge_collapser, const Vec3d& queued_point,
    size_t remaining_vertices, size_t target_vertices_num,
    Phase2PlacementDecision& decision, double& newton_seconds);
  bool solve_phase2_qem_placement(
    const Phase2PlacementContext& ctx, Vec3d& new_point, double& energy) const;
  bool solve_phase2_quadratic_surrogate(
    const Phase2PlacementContext& ctx, Vec3d& new_point, double& energy) const;
  bool compute_phase2_queue_placement_candidate(EdgeHandle eh, Vec3d& new_point, double& energy);
  bool enqueue_phase2_candidate(EdgeHandle eh, size_t state);
  void initialize_phase2_candidates();
  void update_phase2_after_collapsing(VertexHandle collapsed_center);
  bool refine_vertex_relocation_with_newton(
    VertexHandle vh, VertexRelocater& vertex_relocater, Vec3d& new_point, double& energy);
  void do_phase2_final_newton_relocation();
  EdgeSide classify_edge_side(EdgeHandle eh) const;
  EdgeSideStats collect_candidate_edge_side_stats() const;
  void add_edge_side(EdgeSideStats& stats, EdgeSide side) const;
  void log_edge_side_stats(const char* label, const EdgeSideStats& stats) const;


  inline EdgeCollapser new_edge_collapser() { return EdgeCollapser(om, rm, ot, lrt, og); }
  inline VertexRelocater new_vertex_relocater() { return VertexRelocater(om, rm, ot, lrt, nullptr, og); }
};

}// namespace CageSimp
}// namespace Cage
