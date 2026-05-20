#pragma once

#include <queue>

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

  void do_collapse(size_t edge_num_to_collapse, size_t& total_collapsed_edge_num);
  void update(size_t _candidate_points_size, bool _allow_negtive, double _max_distance_error);
private:
  double avg_edge_length;
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

  std::vector<Vec3d> generate_candidate_points_for_collapse(EdgeHandle e, EdgeCollapser& edge_collapser);
  bool find_collapse_hausdorff_deviation(
    EdgeHandle eh, double& local_hd_before, double& local_hd_after, Vec3d& new_point, double& priority_score);
  void initialize_collapse_edges_reward();
  void update_after_collapsing(VertexHandle collapsed_center);
  bool try_enqueue_collapse_candidate(
    EdgeHandle eh, size_t state, double local_hd_before, double local_hd_after, const Vec3d& new_point, double priority_score);
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
  EdgeSide classify_edge_side(EdgeHandle eh) const;
  EdgeSideStats collect_candidate_edge_side_stats() const;
  void add_edge_side(EdgeSideStats& stats, EdgeSide side) const;
  void log_edge_side_stats(const char* label, const EdgeSideStats& stats) const;


  inline EdgeCollapser new_edge_collapser() { return EdgeCollapser(om, rm, ot, lrt, og); }
};

}// namespace CageSimp
}// namespace Cage
