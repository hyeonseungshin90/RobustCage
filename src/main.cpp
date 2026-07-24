#include "CageGenerator.hh"
#include "boost/filesystem.hpp"
#include "boost/algorithm/string.hpp"
#include <cctype>
#include <cfloat>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <sstream>

using namespace Cage;
namespace bf = boost::filesystem;
namespace bj = boost::json;

struct MeshDiagnostics
{
  size_t boundary_edges = 0;
  size_t non_manifold_vertices = 0;
};

MeshDiagnostics analyze_mesh(Cage::SMeshT& mesh)
{
  MeshDiagnostics diagnostics;
  for (auto eh : mesh.edges())
  {
    if (mesh.is_boundary(eh))
      diagnostics.boundary_edges++;
  }
  for (auto vh : mesh.vertices())
  {
    if (!mesh.is_manifold(vh))
      diagnostics.non_manifold_vertices++;
  }
  return diagnostics;
}

bool mesh_valid(Cage::SMeshT& mesh)
{
  return mesh.n_vertices() > 0 && mesh.n_faces() > 0;
}

double calc_triangle_quality(Cage::SM::SMeshT& mesh, Cage::SM::FaceHandle fh)
{
  Cage::SM::Vec3d pts[3];
  size_t vertex_count = 0;
  for (Cage::SM::VertexHandle vh : mesh.fv_range(fh))
  {
    if (vertex_count >= 3)
      return 0.0;
    pts[vertex_count++] = mesh.point(vh);
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

struct TriangleQualityStats
{
  size_t face_count = 0;
  double min_quality = DBL_MAX;
  double max_quality = 0.0;
  double avg_quality = 0.0;
  Cage::SM::FaceHandle min_face = Cage::SM::FaceHandle(-1);
  Cage::SM::FaceHandle max_face = Cage::SM::FaceHandle(-1);
};

TriangleQualityStats calc_triangle_quality_stats(Cage::SM::SMeshT& mesh)
{
  TriangleQualityStats stats;
  for (Cage::SM::FaceHandle fh : mesh.faces())
  {
    const double quality = calc_triangle_quality(mesh, fh);
    if (quality < stats.min_quality)
    {
      stats.min_quality = quality;
      stats.min_face = fh;
    }
    if (quality > stats.max_quality)
    {
      stats.max_quality = quality;
      stats.max_face = fh;
    }
    stats.avg_quality += quality;
    stats.face_count++;
  }

  if (stats.face_count == 0)
  {
    stats.min_quality = 0.0;
    stats.max_quality = 0.0;
  }
  else
    stats.avg_quality /= static_cast<double>(stats.face_count);
  return stats;
}

void log_triangle_quality_stats(Cage::SM::SMeshT& mesh, const std::string& label)
{
  const TriangleQualityStats stats = calc_triangle_quality_stats(mesh);
  Logger::user_logger->info(
    "{} triangle quality: min {} (face {}), max {} (face {}), avg {}, faces {}.",
    label,
    stats.min_quality,
    stats.min_face.idx(),
    stats.max_quality,
    stats.max_face.idx(),
    stats.avg_quality,
    stats.face_count);
}

std::string sanitize_path_component(std::string value)
{
  if (value.empty())
    return "unknown";

  for (char& ch : value)
  {
    const unsigned char uch = static_cast<unsigned char>(ch);
    if (!std::isalnum(uch) && ch != '_' && ch != '-')
      ch = '_';
  }
  return value;
}

std::string format_timestamp(std::time_t time_value)
{
  if (time_value <= 0)
    return "unknown";

  std::tm local_time{};
#ifdef _WIN32
  localtime_s(&local_time, &time_value);
#else
  localtime_r(&time_value, &local_time);
#endif

  std::ostringstream oss;
  oss << std::put_time(&local_time, "%Y%m%d_%H%M%S");
  return oss.str();
}

std::string collapse_mode_label(const Cage::ParamCollapseStage& param)
{
  std::string mode = sanitize_path_component(param.priorityMode);
  if (mode == "length_quality")
    mode += "_" + sanitize_path_component(param.lengthQualitySubMode);
  if (param.collapsePlacementMethod == "optimization" ||
    param.collapsePlacementMethod == "energy" ||
    param.collapsePlacementMethod == "newton")
  {
    mode += "_optimization";
    mode += "_" + sanitize_path_component(param.phase2PlacementStrategy);
    mode += "_" + sanitize_path_component(param.newtonSolverMode);
    mode += "_" + sanitize_path_component(param.robustnessMode);
    if (param.curvatureMode != "none")
      mode += "_curv_" + sanitize_path_component(param.curvatureMode);
    if (param.uniformityMode != "none")
      mode += "_uniform_" + sanitize_path_component(param.uniformityMode);
  }
  return "collapse_" + mode;
}

std::string phase2_mode_label(const Cage::ParamCageSimplifier& param)
{
  return "phase2_" + sanitize_path_component(param.phase2Mode);
}

std::string newton_phase2_detail_label(const Cage::ParamCollapseStage& param)
{
  std::string mode = "newton";
  mode += "_" + sanitize_path_component(param.phase2PlacementStrategy);
  mode += "_" + sanitize_path_component(param.newtonSolverMode);
  mode += "_" + sanitize_path_component(param.robustnessMode);
  if (param.curvatureMode != "none")
    mode += "_curv_" + sanitize_path_component(param.curvatureMode);
  if (param.uniformityMode != "none")
    mode += "_uniform_" + sanitize_path_component(param.uniformityMode);
  return mode;
}

std::string linear_phase2_detail_label(const Cage::ParamCollapseStage& param)
{
  std::string mode = "linear_only";
  if (param.curvatureMode != "none")
    mode += "_curv_" + sanitize_path_component(param.curvatureMode);
  if (param.uniformityMode != "none")
    mode += "_uniform_" + sanitize_path_component(param.uniformityMode);
  return mode;
}

std::string flip_mode_label(const Cage::ParamFlipStage& param)
{
  return "flip_" + sanitize_path_component(param.priorityMode);
}

std::string relocate_mode_label(const Cage::ParamRelocateStage& param)
{
  return "relocate_" + sanitize_path_component(param.priorityMode);
}

std::string build_run_dir_name(const Cage::ParamCageGenerator& param)
{
  const auto& simplifier = param.paramCageSimplifier;
  std::ostringstream oss;
  oss
    << format_timestamp(std::time(nullptr))
    << "__" << phase2_mode_label(simplifier);
  if (simplifier.phase2Mode == "newton")
  {
    oss << "__" << newton_phase2_detail_label(simplifier.paramCollapse);
  }
  else if (simplifier.phase2Mode == "linear_only")
  {
    oss << "__" << linear_phase2_detail_label(simplifier.paramCollapse);
  }
  else if (simplifier.phase2Mode != "qem" &&
    simplifier.phase2Mode != "qem_no_collision")
  {
    oss
      << "__" << collapse_mode_label(simplifier.paramCollapse)
      << "__" << flip_mode_label(simplifier.paramFlip)
      << "__" << relocate_mode_label(simplifier.paramRelocate);
  }
  return oss.str();
}

bf::path create_unique_output_dir(const bf::path& parent_dir, const std::string& dir_name)
{
  for (size_t i = 0;; i++)
  {
    std::ostringstream suffix;
    if (i > 0)
      suffix << "_" << std::setw(3) << std::setfill('0') << i;

    bf::path candidate = parent_dir;
    candidate.append(dir_name + suffix.str());
    if (!bf::exists(candidate))
    {
      bf::create_directory(candidate);
      return candidate;
    }
  }
}

void generate_cages(
  Cage::ParamCageGenerator param,
  bf::path in_model_path,
  bf::path out_data_path,
  std::vector<size_t> target_vn
)
{
  // parse file name
  std::string file_path = in_model_path.string();
  std::vector<std::string> split_file_path;
  boost::split(split_file_path, file_path, boost::is_any_of("/\\"), boost::token_compress_on);
  std::string file_name = split_file_path.back();
  file_name = file_name.substr(0, file_name.find_last_of("."));

  CageGenerator cage_generator;
  cage_generator.param = param;
  cage_generator.originalMesh = std::make_unique<SMeshT>();
  try
  {
    // create output directory
    bf::path input_out_dir = out_data_path;
    input_out_dir.append(file_name);
    if (bf::exists(input_out_dir) && !bf::is_directory(input_out_dir))
      throw logic_error("input output path already exists as a file");
    if (!bf::exists(input_out_dir))
      bf::create_directory(input_out_dir);

    bf::path file_out_dir = create_unique_output_dir(
      input_out_dir,
      build_run_dir_name(param));
    // create log file
    bf::path log_path = file_out_dir;
    log_path.append("log.txt");
    Logger::updateFileLog(true, spdlog::level::trace, log_path.string());
    Logger::user_logger->info("processing {}", file_name);
    Logger::user_logger->info("output directory: {}", file_out_dir.string());
    // read input mesh
    if (!OpenMesh::IO::read_mesh(*cage_generator.originalMesh, file_path))
    {
      Logger::user_logger->warn("fail to read input mesh: {}", file_path);
      throw logic_error("fail to read input mesh");
    }
    // check input
    if (!mesh_valid(*cage_generator.originalMesh))
    {
      Logger::user_logger->warn("invalid mesh: input mesh has no vertices or faces.");
      throw logic_error("invalid mesh");
    }
    const MeshDiagnostics mesh_diagnostics = analyze_mesh(*cage_generator.originalMesh);
    Logger::user_logger->info(
      "input mesh: {} vertices, {} edges, {} faces.",
      cage_generator.originalMesh->n_vertices(),
      cage_generator.originalMesh->n_edges(),
      cage_generator.originalMesh->n_faces());
    if (mesh_diagnostics.non_manifold_vertices > 0)
    {
      Logger::user_logger->warn(
        "input mesh has {} non-manifold vertices; continuing anyway.",
        mesh_diagnostics.non_manifold_vertices);
    }
    if (mesh_diagnostics.boundary_edges > 0)
    {
      Logger::user_logger->warn(
        "input mesh has {} boundary edges (non-watertight); continuing anyway.",
        mesh_diagnostics.boundary_edges);
    }
    // set output dir and filename.
    cage_generator.param.setOutputPath(file_out_dir.string() + "/", file_name);
    // set fast simplification target for cages' construction.
    cage_generator.param.setFastTargetNumber(cage_generator.originalMesh->n_vertices() * 3);

    // do cage or nested cages generation
    for (size_t it = 0;it < target_vn.size();it++)
    {
      cage_generator.param.setCageLabel(it);
      cage_generator.param.setTargetNumber(target_vn[it]);
      cage_generator.generate();
      log_triangle_quality_stats(
        *cage_generator.cage,
        file_name + "_cage_" + std::to_string(it));

      bf::path mesh_out_file = file_out_dir;
      mesh_out_file.append(file_name + "_cage_" + std::to_string(it) + ".obj");
      OpenMesh::IO::write_mesh(*cage_generator.cage, mesh_out_file.string(), OpenMesh::IO::Options::Default, 15);

      *cage_generator.originalMesh = *cage_generator.cage;
    }
  }
  catch (BreakoutExcept be)
  {
    UNUSED(be);
  }
  catch (AssertFailExcept ae)
  {
    UNUSED(ae);
  }
  catch (...)
  {
    Logger::dev_logger->error("unexpected exception, check it!");
  }

  Logger::dev_logger->flush();
  Logger::user_logger->flush();
}

std::string normalize_parameter_token(std::string token)
{
  boost::trim(token);
  boost::replace_all(token, "-", "_");
  return token;
}

void apply_qem_energy_weights(
  Cage::ParamCollapseStage& collapse, bool include_triangle_quality)
{
  collapse.qemWeight = 1.0;
  collapse.selfBarrierWeight = 0.0;
  collapse.originalBarrierWeight = 0.0;
  collapse.positionFidelityWeight = 0.0;
  collapse.curvatureWeight = 0.0;
  collapse.triangleQualityWeight = include_triangle_quality ? 1.0 : 0.0;
  collapse.uniformityWeight = 1.0;
}

void enable_newton_phase2_defaults(Cage::ParamCageGenerator& param)
{
  auto& simplifier = param.paramCageSimplifier;
  auto& collapse = simplifier.paramCollapse;
  simplifier.phase2Mode = "newton";
  collapse.collapsePlacementMethod = "optimization";
  collapse.phase2PlacementStrategy = "adaptive";
  collapse.curvatureMode = "weighted_qem";
  collapse.uniformityMode = "source";
  collapse.selfBarrierWeight = 0.0;
  collapse.originalBarrierWeight = 0.0;
  collapse.positionFidelityWeight = 1.0;
  collapse.triangleQualityWeight = 2.0;
  collapse.uniformityWeight = 1.0;
  collapse.phase2NewtonQualityThreshold = 0.12;
  collapse.phase2NewtonResidualThreshold = 0.25;
  collapse.phase2NewtonResidualGrowth = 1.5;
  collapse.phase2NewtonFinalRefineCollapses = 25;
  param.paramCageSimplifier.paramFlip.priorityMode = "triangle_quality_hard";
  param.paramCageSimplifier.paramFlip.requireRegularValence = false;
}

void enable_linear_only_phase2_defaults(Cage::ParamCageGenerator& param)
{
  enable_newton_phase2_defaults(param);
  auto& simplifier = param.paramCageSimplifier;
  auto& collapse = simplifier.paramCollapse;
  simplifier.phase2Mode = "linear_only";
  collapse.phase2PlacementStrategy = "linear_only";
}

void enable_qem_phase2_defaults(Cage::ParamCageGenerator& param, bool skip_original_collision_check = false)
{
  auto& simplifier = param.paramCageSimplifier;
  auto& collapse = simplifier.paramCollapse;
  simplifier.phase2Mode = skip_original_collision_check ? "qem_no_collision" : "qem";
  collapse.collapsePlacementMethod = "optimization";
  collapse.phase2PlacementStrategy =
    skip_original_collision_check ? "qem_no_collision" : "qem";
  collapse.robustnessMode = "exact_reject";
  if (skip_original_collision_check)
  {
    collapse.curvatureMode = "none";
    collapse.uniformityMode = "none";
    apply_qem_energy_weights(collapse, false);
    collapse.phase2NewtonQualityThreshold = 0.0;
    collapse.phase2NewtonResidualThreshold = 0.0;
    collapse.phase2NewtonFinalRefineCollapses = 0;
  }
  else
  {
    collapse.curvatureMode = "normal_matching";
    collapse.uniformityMode = "source";
    apply_qem_energy_weights(collapse, true);
  }
}

bool apply_parameter_token(const std::string& raw_token, Cage::ParamCageGenerator& param)
{
  const std::string token = normalize_parameter_token(raw_token);
  auto& simplifier = param.paramCageSimplifier;
  auto& collapse = param.paramCageSimplifier.paramCollapse;
  auto& relocate = param.paramCageSimplifier.paramRelocate;
  auto& flip = param.paramCageSimplifier.paramFlip;

  if (token.empty() || token == "default")
    return true;

  if (token == "phase2_newton" || token == "newton_phase2")
  {
    enable_newton_phase2_defaults(param);
    return true;
  }

  if (token == "collapse" || token == "collapse_default" || token == "collapse_hausdorff")
  {
    collapse.priorityMode = "hausdorff";
    return true;
  }
  if (token == "collapse_length" || token == "length")
  {
    collapse.priorityMode = "length";
    return true;
  }
  if (token == "collapse_length_quality" || token == "collapse_length_quality_weighted" ||
    token == "length_quality" || token == "length_quality_weighted")
  {
    collapse.priorityMode = "length_quality";
    collapse.lengthQualitySubMode = "weighted";
    return true;
  }
  if (token == "collapse_length_quality_relative" || token == "length_quality_relative")
  {
    collapse.priorityMode = "length_quality";
    collapse.lengthQualitySubMode = "relative_reject";
    return true;
  }
  if (token == "collapse_length_quality_absolute" || token == "length_quality_absolute")
  {
    collapse.priorityMode = "length_quality";
    collapse.lengthQualitySubMode = "absolute_reject";
    return true;
  }
  if (token == "collapse_length_quality_lexicographic" || token == "length_quality_lexicographic")
  {
    collapse.priorityMode = "length_quality";
    collapse.lengthQualitySubMode = "lexicographic";
    return true;
  }
  if (token == "collapse_post_edge_length" || token == "post_edge_length")
  {
    collapse.priorityMode = "post_edge_length";
    return true;
  }
  if (token == "collapse_post_edge_length_hard" || token == "post_edge_length_hard")
  {
    collapse.priorityMode = "post_edge_length_hard";
    return true;
  }
  if (token == "collapse_post_face_area" || token == "post_face_area")
  {
    collapse.priorityMode = "post_face_area";
    return true;
  }
  if (token == "collapse_post_face_area_hard" || token == "post_face_area_hard")
  {
    collapse.priorityMode = "post_face_area_hard";
    return true;
  }
  if (token == "collapse_triangle_quality" || token == "triangle_quality")
  {
    collapse.priorityMode = "triangle_quality";
    return true;
  }
  if (token == "collapse_triangle_quality_hard" || token == "triangle_quality_hard")
  {
    collapse.priorityMode = "triangle_quality_hard";
    return true;
  }
  if (token == "collapse_sampling" || token == "sampling")
  {
    collapse.collapsePlacementMethod = "sampling";
    return true;
  }
  if (token == "collapse_optimization" || token == "optimization" ||
    token == "collapse_energy" || token == "energy" ||
    token == "collapse_newton" || token == "newton" || token == "placement_newton")
  {
    enable_newton_phase2_defaults(param);
    return true;
  }
  if (token == "phase2_adaptive" || token == "phase2_placement_adaptive" ||
    token == "adaptive")
  {
    enable_newton_phase2_defaults(param);
    collapse.phase2PlacementStrategy = "adaptive";
    return true;
  }
  if (token == "phase2_linear_only" || token == "phase2_placement_linear_only" ||
    token == "linear_only")
  {
    enable_linear_only_phase2_defaults(param);
    return true;
  }
  if (token == "phase2_qem" || token == "qem" ||
    token == "phase2_qem_only" || token == "phase2_placement_qem_only" ||
    token == "qem_only")
  {
    enable_qem_phase2_defaults(param);
    return true;
  }
  if (token == "phase2_qem_no_collision" ||
    token == "qem_no_collision" ||
    token == "phase2_qem_no_intersection" ||
    token == "qem_no_intersection" ||
    token == "phase2_qem_no_intersection_reject" ||
    token == "qem_no_intersection_reject" ||
    token == "phase2_qem_no_collision_reject" ||
    token == "qem_no_collision_reject")
  {
    enable_qem_phase2_defaults(param, true);
    return true;
  }
  if (token == "phase2_final_newton" || token == "phase2_placement_final_newton" ||
    token == "final_newton")
  {
    enable_newton_phase2_defaults(param);
    collapse.phase2PlacementStrategy = "final_newton";
    return true;
  }
  if (token == "phase2_newton_only" || token == "phase2_placement_newton_only" ||
    token == "newton_only")
  {
    enable_newton_phase2_defaults(param);
    collapse.phase2PlacementStrategy = "newton_only";
    apply_qem_energy_weights(collapse, true);
    return true;
  }
  if (token == "phase2_quadratic_surrogate" || token == "phase2_surrogate" ||
    token == "quadratic_surrogate" || token == "surrogate")
  {
    enable_newton_phase2_defaults(param);
    collapse.phase2PlacementStrategy = "quadratic_surrogate";
    return true;
  }
  if (token == "newton_damped" || token == "solver_damped")
  {
    collapse.newtonSolverMode = "damped";
    return true;
  }
  if (token == "newton_trust_region" || token == "trust_region" || token == "solver_trust_region")
  {
    collapse.newtonSolverMode = "trust_region";
    return true;
  }
  if (token == "robust_exact_reject" || token == "exact_reject")
  {
    collapse.robustnessMode = "exact_reject";
    return true;
  }
  if (token == "robust_exact_backtracking" || token == "exact_backtracking")
  {
    collapse.robustnessMode = "exact_backtracking";
    return true;
  }
  if (token == "robust_ipc" || token == "robust_ipc_line_search" ||
    token == "ipc_line_search" || token == "robust_ccd")
  {
    collapse.robustnessMode = "ipc_line_search";
    return true;
  }
  if (token == "curvature_none" || token == "curv_none")
  {
    collapse.curvatureMode = "none";
    return true;
  }
  if (token == "curvature_weighted_qem" || token == "curv_weighted_qem")
  {
    collapse.curvatureMode = "weighted_qem";
    return true;
  }
  if (token == "curvature_normal_matching" || token == "curv_normal_matching")
  {
    collapse.curvatureMode = "normal_matching";
    return true;
  }
  if (token == "uniformity_none" || token == "uniform_none")
  {
    collapse.uniformityMode = "none";
    return true;
  }
  if (token == "uniformity_source" || token == "uniform_source")
  {
    collapse.uniformityMode = "source";
    return true;
  }
  if (token == "uniformity_global" || token == "uniform_global")
  {
    collapse.uniformityMode = "global";
    return true;
  }
  if (token == "flip_valence")
  {
    flip.priorityMode = "valence";
    return true;
  }
  if (token == "flip_triangle_quality_hard")
  {
    flip.priorityMode = "triangle_quality_hard";
    return true;
  }
  if (token == "relocate_hausdorff")
  {
    relocate.priorityMode = "hausdorff";
    return true;
  }
  if (token == "relocate_triangle_quality_hard")
  {
    relocate.priorityMode = "triangle_quality_hard";
    return true;
  }

  Logger::user_logger->error("unknown parameter token: {}", raw_token);
  return false;
}

bool parse_parameter_arg(const std::string& arg_param, Cage::ParamCageGenerator& param)
{
  bf::path json_file_path(arg_param);
  if (bf::is_regular_file(json_file_path))
  {
    fstream json_file;
    json_file.open(json_file_path.string(), fstream::in);
    if (!json_file.is_open())
    {
      Logger::user_logger->error("fail to open json file.");
      return false;
    }

    std::string json_str((std::istreambuf_iterator<char>(json_file)), std::istreambuf_iterator<char>());
    bj::stream_parser sp;
    sp.write(json_str.c_str());
    param.deserialize(sp.release().as_object());
    json_file.close();
    return true;
  }

  std::vector<std::string> parameter_tokens;
  boost::split(parameter_tokens, arg_param, boost::is_any_of("+,"), boost::token_compress_on);
  for (const std::string& token : parameter_tokens)
  {
    if (!apply_parameter_token(token, param))
      return false;
  }
  return true;
}

int main(int argc, char* argv[])
{
  // args tips
  if (argc < 5)
  {
    printf("Need args:\n");
    printf("arg[0]: parameters.\n");
    printf("input \"default\" to set default parameters\n");
    printf("input \"collapse\" to use default collapse priority; combine with flip/relocate presets to change only those stages\n");
    printf("input \"collapse_length\" to set length-based collapse priority\n");
    printf("input \"collapse_length_quality\" to use weighted length + triangle quality priority\n");
    printf("input \"collapse_length_quality_weighted\" to use weighted quality penalty\n");
    printf("input \"collapse_length_quality_relative\" to reject large relative quality drops\n");
    printf("input \"collapse_length_quality_absolute\" to reject triangles below a quality threshold\n");
    printf("input \"collapse_length_quality_lexicographic\" to use length first and quality as tie-break\n");
    printf("input \"collapse_post_edge_length\" to use post-collapse edge length priority\n");
    printf("input \"collapse_post_edge_length_hard\" to require post-collapse max edge length not to increase\n");
    printf("input \"collapse_post_face_area\" to use post-collapse face area priority\n");
    printf("input \"collapse_post_face_area_hard\" to require post-collapse max face area not to increase\n");
    printf("input \"collapse_triangle_quality\" to use triangle quality priority\n");
    printf("input \"collapse_triangle_quality_hard\" to require post-collapse min triangle quality not to decrease\n");
    printf("input \"phase2_newton\" to replace Phase 2 with optimization collapses and default curvature/uniformity energies\n");
    printf("input \"phase2_qem\" to run QEM + quality/uniformity/curvature/dihedral Phase 2 collapses with QEM-first backtracking and hard intersection constraints\n");
    printf("input \"phase2_qem_no_collision\" to run pure Garland-Heckbert QEM without collision rejection\n");
    printf("input \"phase2_quadratic_surrogate\" to use the separate 4x4 quadratic surrogate Phase 2 strategy\n");
    printf("input \"collapse_optimization\" to use the optimization Phase 2 replacement with default curvature/uniformity energies\n");
    printf("input \"newton_damped\" or \"newton_trust_region\" to choose the Newton solver\n");
    printf("input \"robust_exact_reject\", \"robust_exact_backtracking\", or \"robust_ipc\" to choose robustness handling\n");
    printf("input \"curvature_weighted_qem\", \"curvature_normal_matching\", or \"curvature_none\" to choose curvature energy\n");
    printf("input \"uniformity_source\", \"uniformity_global\", or \"uniformity_none\" to choose one uniformity energy\n");
    printf("input \"flip_triangle_quality_hard\" to require post-flip min triangle quality not to decrease\n");
    printf("input \"relocate_triangle_quality_hard\" to require post-relocate min triangle quality not to decrease\n");
    printf("combine presets with '+' or ',', for example \"collapse_length+flip_triangle_quality_hard+relocate_triangle_quality_hard\".\n");
    printf("or a json file to set parameters.\n");
    printf("arg[1]: input model path.\n");
    printf("arg[2]: output dir path.\n");
    printf("(optional)arg[2]: vertices number of nested cage 0.\n");
    printf("(optional)arg[3]: vertices number of nested cage 1.\n");
    printf("(optional)arg[n + 2]: vertices number of nested cage n.\n");
    return 1;
  }

  // logger and random seed.
  Logger::InitLogger(spdlog::level::level_enum::trace);
  srand((unsigned)time(NULL));

  // parse parameters
  std::string arg_param(argv[1]);
  Cage::ParamCageGenerator param;
  if (!parse_parameter_arg(arg_param, param))
    return 1;

  // parse input/output file/directory.
  bf::path in_model_path(argv[2]);
  bf::path out_data_path(argv[3]);

  if (!bf::is_regular_file(in_model_path))
  {
    Logger::user_logger->error("error input file: {}", in_model_path.string());
    return 1;
  }

  if (!bf::exists(out_data_path))
  {
    boost::system::error_code ec;
    if (!bf::create_directories(out_data_path, ec) || ec)
    {
      Logger::user_logger->error("error output directory: {}", out_data_path.string());
      return 1;
    }
  }
  else if (!bf::is_directory(out_data_path))
  {
    Logger::user_logger->error("error output directory: {}", out_data_path.string());
    return 1;
  }

  // parse target vertices numbers.
  int cage_number = argc - 4;
  std::vector<size_t> target_vn;
  try
  {
    for (int i = 0;i < cage_number;i++)
    {
      std::string vn_str(argv[4 + i]);
      size_t vn = std::stol(vn_str);
      target_vn.push_back(vn);
    }
  }
  catch (...)
  {
    Logger::user_logger->error("error in parsing arguments.");
    return 1;
  }

  // generate!
  generate_cages(param, in_model_path, out_data_path, target_vn);
  return 0;
}
