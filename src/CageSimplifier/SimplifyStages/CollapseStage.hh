#pragma once

#include <cfloat>
#include <queue>
#include "Dense"

#include "Config.hh"
#include "CageSimplifier/Topo/TopoOperations.h"
#include "CageSimplifier/Topo/EdgeCollapser.h"
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
  ~CollapseStage();

  void do_collapse(size_t edge_num_to_collapse, size_t& total_collapsed_edge_num);
  // Replaces legacy Phase 2 by repeatedly collapsing edges until the target
  // vertex count is reached. Each collapse chooses its new vertex position
  // before the topology change is committed.
  void do_phase2_energy_simplification(size_t target_vertices_num);
  void update(size_t _candidate_points_size, bool _allow_negtive, double _max_distance_error);
private:
  double avg_edge_length;
  double avg_source_edge_length;
  bool avg_source_edge_length_initialized;
  double phase2_target_edge_length;
  std::vector<size_t> update_states;

  struct CollapseEdgeReward
  {
    EdgeHandle eh;
    size_t state;
    double reward;
    Vec3d new_point;

    CollapseEdgeReward() = default;
    CollapseEdgeReward(EdgeHandle _eh, size_t _state, double _reward, const Vec3d& _new_point) :
      eh(_eh), state(_state), reward(_reward), new_point(_new_point)
    {}

    bool operator<(const CollapseEdgeReward& rhs)const { return reward < rhs.reward; }
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
  OpenMesh::VPropHandleT<Eigen::Matrix4d> phase2_qem_quadric_prop;
  bool phase2_qem_quadric_prop_added = false;

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
    Eigen::Matrix4d qem_matrix = Eigen::Matrix4d::Zero();
    bool use_qem_matrix = false;
    Vec3d endpoint0;
    Vec3d endpoint1;
    Vec3d tangential_smoothing_point;
    Vec3d midpoint;
    double local_scale = 1.0;
    double edge_curvature = 0.0;
    double avg_cage_edge_length = 1.0;
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

  struct Phase2QueueComponents
  {
    double qem = 0.0;
    double triangle_quality = 0.0;
    double uniformity = 0.0;
  };

  Phase2QueueComponents phase2_queue_component_maxima;

  std::vector<Vec3d> generate_candidate_points_for_collapse(EdgeHandle e, EdgeCollapser& edge_collapser);
  bool find_collapse_hausdorff_deviation(
    EdgeHandle eh, double& local_hd_before, double& local_hd_after, Vec3d& new_point);
  bool refine_collapse_placement_with_newton(
    EdgeHandle eh, EdgeCollapser& edge_collapser, Vec3d& new_point, double& energy,
    const Vec3d* initial_point = nullptr);
  void initialize_collapse_edges_reward();
  void update_after_collapsing(VertexHandle collapsed_center);
  bool is_optimization_collapse_placement_method() const;
  bool is_phase2_qem_based_strategy() const;
  bool is_phase2_linear_solve_strategy() const;
  bool is_phase2_qem_original_strategy() const;
  bool is_phase2_newton_solve_strategy() const;
  bool is_trust_region_solver_mode() const;
  bool is_exact_reject_robustness_mode() const;
  bool is_exact_backtracking_robustness_mode() const;
  Phase2PlacementContext make_phase2_placement_context(EdgeHandle eh, EdgeCollapser& edge_collapser) const;
  NewtonDerivatives finite_difference_newton_derivatives(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  NewtonDerivatives approximate_newton_derivatives(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  NewtonDerivatives autodiff_newton_derivatives(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double evaluate_newton_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double evaluate_qem_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double evaluate_qem_only_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double calc_triangle_quality(SMeshT* mesh, FaceHandle fh) const;
  double evaluate_triangle_quality_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  Eigen::Matrix4d build_phase2_triangle_quality_surrogate_quadric(
    const Phase2PlacementContext& ctx) const;
  double evaluate_phase2_triangle_quality_surrogate(
    const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double evaluate_pre_collapse_triangle_quality_penalty(const Phase2PlacementContext& ctx) const;
  Eigen::Matrix4d build_phase2_uniformity_difference_quadric(
    const Phase2PlacementContext& ctx) const;
  double evaluate_phase2_uniformity_difference(
    const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double evaluate_uniformity_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double evaluate_phase2_proxy_energy(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  Phase2QueueComponents evaluate_phase2_queue_components(
    const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double evaluate_phase2_queue_score(const Phase2QueueComponents& components) const;
  double evaluate_phase2_refinement_residual(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double calc_phase2_fan_min_quality(const Phase2PlacementContext& ctx, const Vec3d& x) const;
  double source_uniformity_target_length(const Phase2PlacementContext& ctx, const Vec3d& sample_point) const;
  double source_edge_length_at(const Vec3d& sample_point) const;
  bool collapse_target_valid(EdgeCollapser& edge_collapser, const Vec3d& x) const;
  bool phase2_placement_satisfies_hard_constraints(
    const Phase2PlacementContext& ctx, EdgeCollapser& edge_collapser, const Vec3d& x) const;
  void set_phase2_edge_collapser_flags(EdgeCollapser& edge_collapser) const;
  bool choose_phase2_collapse_placement(
    EdgeHandle eh, EdgeCollapser& edge_collapser, const Vec3d& queued_point,
    Phase2PlacementDecision& decision, double& newton_seconds);
  bool solve_phase2_quadric_placement(
    const Phase2PlacementContext& ctx, Vec3d& new_point, double& energy) const;
  bool select_phase2_linear_solve_line_search_candidate(
    const Phase2PlacementContext& ctx, EdgeCollapser& edge_collapser,
    const Vec3d& qem_point, Vec3d& selected_point, double& selected_energy) const;
  void initialize_phase2_qem_quadrics();
  void update_phase2_qem_quadric_after_collapse(VertexHandle center_vh, const Eigen::Matrix4d& quadric);
  Eigen::Matrix4d phase2_qem_quadric(VertexHandle vh) const;
  bool select_phase2_feasibility_fallback(
    EdgeHandle eh, const Phase2PlacementContext& ctx, EdgeCollapser& edge_collapser,
    Vec3d& new_point, double& energy, bool qem_only_score = false) const;
  bool compute_phase2_queue_placement_candidate(EdgeHandle eh, Vec3d& new_point, double& energy);
  bool compute_phase2_queue_candidate_data(
    EdgeHandle eh, Vec3d& new_point, Phase2QueueComponents& components);
  bool enqueue_phase2_candidate(EdgeHandle eh, size_t state);
  void refresh_phase2_average_lengths();
  void initialize_phase2_target_edge_length(size_t target_vertices_num);
  void initialize_phase2_candidates();
  void update_phase2_after_collapsing(VertexHandle collapsed_center);
  EdgeSide classify_edge_side(EdgeHandle eh) const;
  EdgeSideStats collect_candidate_edge_side_stats() const;
  void add_edge_side(EdgeSideStats& stats, EdgeSide side) const;
  void log_edge_side_stats(const char* label, const EdgeSideStats& stats) const;


  inline EdgeCollapser new_edge_collapser() { return EdgeCollapser(om, rm, ot, lrt, og); }
};

}// namespace CageSimp
}// namespace Cage
