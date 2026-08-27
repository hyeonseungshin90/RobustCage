#include "VertexRounder.hh"

namespace Cage
{
namespace CageInit
{

bool VertexRounder::doRounding(VertexHandle vh)
{
  // constrained point, skip.
  if (VMesh->vertex(vh).prop.is_constraint)
    return false;
  // exact point and float point are same, skip.
  if (VMesh->exact_point(vh).isSame())
    return false;
  // get rounded point
  Point_3 round_ep = VMesh->exact_point(vh).round();
  // for each connected cell, check volume after rounding
  bool all_positive = true;
  for (CellHandle ch : VMesh->getConnCells(vh))
  {
    if (!checkPositiveVolume(ch, vh, round_ep))
    {
      all_positive = false;
      break;
    }
  }
  if (all_positive)
  {
    // do real rounding
    VMesh->exact_point(vh).rounded(VMesh->point(vh));
    return true;
  }
  else return false;
}

size_t VertexRounder::doRounding(const std::vector<VertexHandle>& vrts)
{
  Logger::user_logger->info("rounding vertex.");
  size_t rounded_vertex = 0;
  for (VertexHandle vh : vrts)
  {
    rounded_vertex += doRounding(vh);
  }
  Logger::user_logger->info("rounded {} vertices.", rounded_vertex);
  return rounded_vertex;
}

size_t VertexRounder::doRounding()
{
  Logger::user_logger->info("rounding vertex.");
  size_t rounded_vertex = 0;
  for (size_t vidx = 0;vidx < VMesh->nVertices();vidx++)
  {
    rounded_vertex += doRounding(VertexHandle(vidx));
  }
  Logger::user_logger->info("rounded {} vertices.", rounded_vertex);
  return rounded_vertex;
}

bool VertexRounder::checkPositiveVolume(CellHandle ch, VertexHandle vh, const Point_3& approx_p)
{
  HalfFaceHandle opp_hfh = VMesh->halffaceOppositeVertex(ch, vh);
  std::array<VertexHandle, 3> fvs = VMesh->findHFV(opp_hfh);
  CGAL::Sign ori = CGAL::orientation(
    VMesh->exact_point(fvs[0]).exact(), VMesh->exact_point(fvs[1]).exact(),
    VMesh->exact_point(fvs[2]).exact(), approx_p);
  return ori == CGAL::ON_POSITIVE_SIDE;
}

bool VertexRounder::checkAllTets(const char* stage)
{
  // first we count rounded points.
  size_t rounded_vertices = 0;
  size_t unconstrained_vertices = 0;
  for (size_t vidx = 0;vidx < VMesh->vertices.size();vidx++)
  {
    VertexHandle vh(vidx);
    if (!VMesh->vertex(vh).prop.is_constraint)
    {
      if (VMesh->exact_point(vh).isSame())
        rounded_vertices++;
      unconstrained_vertices++;
    }
  }
  Logger::user_logger->info("[{}] rounded_vertices ratio = {}", stage, (double)rounded_vertices / (double)unconstrained_vertices);
  // second we check tets' volume.  A healthy cell costs a single exact
  // predicate; the other three halffaces are only evaluated to separate a cell
  // that is merely stored with reversed orientation from one that is broken.
  size_t live_cells = 0;
  size_t flipped_cells = 0;      // all four halffaces negative: reversed storage
  size_t degenerate_cells = 0;   // some halfface coplanar: zero volume
  size_t inconsistent_cells = 0; // halffaces disagree: broken orientation
  size_t dumped = 0;
  const size_t max_dump = 10;

  for (size_t cidx = 0;cidx < VMesh->cells.size();cidx++)
  {
    CellHandle ch(cidx);
    if (VMesh->deleted(ch))
      continue;
    live_cells++;

    std::array<CGAL::Sign, 4> ori;
    size_t n_ori = 0;
    for (HalfFaceHandle cf : VMesh->cell(ch).halffaces)
    {
      std::array<VertexHandle, 3> fvs = VMesh->findHFV(cf);
      VertexHandle opp_vh = VMesh->vertexOppositeFace(ch, cf);
      ori[n_ori++] = CGAL::orientation(
        VMesh->exact_point(fvs[0]).exact(), VMesh->exact_point(fvs[1]).exact(),
        VMesh->exact_point(fvs[2]).exact(), VMesh->exact_point(opp_vh).exact());
      if (n_ori == 1 && ori[0] == CGAL::ON_POSITIVE_SIDE)
        break;
    }
    if (n_ori == 1)
      continue;

    size_t n_positive = 0, n_negative = 0, n_zero = 0;
    for (size_t i = 0;i < n_ori;i++)
    {
      if (ori[i] == CGAL::ON_POSITIVE_SIDE) n_positive++;
      else if (ori[i] == CGAL::ON_NEGATIVE_SIDE) n_negative++;
      else n_zero++;
    }

    if (n_zero > 0)
      degenerate_cells++;
    else if (n_negative == n_ori)
      flipped_cells++;
    else
      inconsistent_cells++;

    if (dumped < max_dump)
    {
      dumped++;
      std::array<VertexHandle, 4> cvs = VMesh->findCV(ch);
      Logger::user_logger->warn(
        "[{}] bad tet cell {}: halfface signs +{} -{} 0{}; rounded {}{}{}{}; "
        "v{} ({} {} {}), v{} ({} {} {}), v{} ({} {} {}), v{} ({} {} {}).",
        stage, cidx, n_positive, n_negative, n_zero,
        VMesh->exact_point(cvs[0]).isSame() ? 1 : 0,
        VMesh->exact_point(cvs[1]).isSame() ? 1 : 0,
        VMesh->exact_point(cvs[2]).isSame() ? 1 : 0,
        VMesh->exact_point(cvs[3]).isSame() ? 1 : 0,
        cvs[0].idx(), VMesh->point(cvs[0]).x(), VMesh->point(cvs[0]).y(), VMesh->point(cvs[0]).z(),
        cvs[1].idx(), VMesh->point(cvs[1]).x(), VMesh->point(cvs[1]).y(), VMesh->point(cvs[1]).z(),
        cvs[2].idx(), VMesh->point(cvs[2]).x(), VMesh->point(cvs[2]).y(), VMesh->point(cvs[2]).z(),
        cvs[3].idx(), VMesh->point(cvs[3]).x(), VMesh->point(cvs[3]).y(), VMesh->point(cvs[3]).z());
    }
  }

  Logger::user_logger->info(
    "[{}] tet volume check: {} live cells, {} degenerate, {} inconsistent, {} reversed-orientation.",
    stage, live_cells, degenerate_cells, inconsistent_cells, flipped_cells);

  return degenerate_cells == 0 && inconsistent_cells == 0;
}

}// namespace CageInit
}// namespace Cage