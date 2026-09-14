#pragma once

#include <queue>

#include "Config.hh"
#include "CageSimplifier/Topo/TopoOperations.h"
#include "CageSimplifier/Topo/EdgeFlipper.h"
#include "CageSimplifier/Geom/GeometryCheck.h"

namespace Cage
{
namespace CageSimp
{
class FlipStage
{
public:
  ParamFlipStage* param;
  SMeshT* om;
  SMeshT* rm;
  DFaceTree* ot;
  LightDFaceTree* lrt;
  FaceGrid* og;

  double original_diagonal_length;
  bool allow_negtive;
  double max_distance_error;
public:
  FlipStage(
    SMeshT* original, SMeshT* cage, ParamFlipStage* p,
    DFaceTree* orginal_tree, LightDFaceTree* remeshing_tree,
    FaceGrid* original_grid,
    double _original_diagonal_length);

  void do_flip();
  // Fixed-position quality polish for the linear Phase 2 path. The optional
  // chord-aware mode prioritizes rail-chord removal and forbids creating
  // chords. Returns accepted flips; no Hausdorff sampling links are needed.
  size_t do_quality_flip(bool prioritize_rail_chords = false);
  void update(bool _allow_negtive, double _max_distance_error);
private:
  std::vector<size_t> update_states;

  struct FlipEdgeReward
  {
    EdgeHandle eh;
    size_t state;
    double reward;
    bool removes_rail_chord = false;

    FlipEdgeReward() = default;
    FlipEdgeReward(EdgeHandle _eh, size_t _state, double _reward,
      bool _removes_rail_chord = false) :
      eh(_eh), state(_state), reward(_reward),
      removes_rail_chord(_removes_rail_chord)
    {}

    bool operator<(const FlipEdgeReward& rhs)const
    {
      if (removes_rail_chord != rhs.removes_rail_chord)
        return removes_rail_chord < rhs.removes_rail_chord;
      return reward < rhs.reward;
    }
  };
  typedef std::priority_queue<FlipEdgeReward> FlipEdgeRewardQueue;
  FlipEdgeRewardQueue edges_to_flip;

  void initialize_flip_edges_reward();
  void update_after_flipping(EdgeHandle flipped_edge);
  bool try_enqueue_flip_candidate(EdgeHandle eh, size_t state, double local_hd_before, double local_hd_after);
  bool is_triangle_quality_priority_mode() const;
  bool is_triangle_quality_hard_priority_mode() const;
  bool should_require_regular_valence() const;
  bool is_flip_quality_allowed(EdgeHandle eh) const;
  double calc_flip_quality_delta(EdgeHandle eh) const;
  double calc_pre_flip_quality(EdgeHandle eh) const;
  double calc_post_flip_quality(EdgeHandle eh) const;
  double calc_triangle_quality(FaceHandle fh) const;
  double calc_triangle_quality(const Vec3d& p0, const Vec3d& p1, const Vec3d& p2) const;

  inline EdgeFlipper new_edge_flipper() { return EdgeFlipper(om, rm, ot, lrt, og); }
};
}// namespace CageSimp
}// namespace Cage
