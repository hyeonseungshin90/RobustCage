#include "CageGenerator.hh"
#include "boost/filesystem.hpp"
#include "boost/algorithm/string.hpp"
#include <array>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <cfloat>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <vector>

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

// OpenMesh 8.1 identifies a .stl file as ASCII as soon as its header starts
// with "solid".  Binary STL headers are arbitrary text and are also allowed to
// start with that word, so use the binary STL byte layout before asking the
// reader to parse the stream.
bool has_binary_stl_layout(const bf::path& path)
{
  boost::system::error_code ec;
  const std::uintmax_t file_size = bf::file_size(path, ec);
  if (ec || file_size < 84)
    return false;

  std::ifstream input(path.string().c_str(), std::ios::in | std::ios::binary);
  if (!input.is_open())
    return false;

  input.seekg(80, std::ios::beg);
  std::array<unsigned char, 4> count_bytes{};
  input.read(
    reinterpret_cast<char*>(count_bytes.data()),
    static_cast<std::streamsize>(count_bytes.size()));
  if (!input)
    return false;

  const std::uint32_t triangle_count =
    static_cast<std::uint32_t>(count_bytes[0]) |
    (static_cast<std::uint32_t>(count_bytes[1]) << 8) |
    (static_cast<std::uint32_t>(count_bytes[2]) << 16) |
    (static_cast<std::uint32_t>(count_bytes[3]) << 24);
  const std::uintmax_t expected_size =
    84u + static_cast<std::uintmax_t>(triangle_count) * 50u;
  return expected_size == file_size;
}

bool read_input_mesh(
  Cage::SMeshT& mesh,
  const bf::path& path,
  std::string& format_description)
{
  const std::string extension =
    boost::algorithm::to_lower_copy(path.extension().string());
  const bool is_stl =
    extension == ".stl" || extension == ".stla" || extension == ".stlb";
  if (!is_stl)
  {
    format_description = extension.empty() ? "unknown" : extension.substr(1);
    return OpenMesh::IO::read_mesh(mesh, path.string());
  }

  const bool is_binary = extension == ".stlb" ||
    (extension == ".stl" && has_binary_stl_layout(path));
  format_description = is_binary ? "binary STL" : "ASCII STL";

  // The stream overload obeys Options::Binary directly and therefore avoids
  // OpenMesh's ambiguous filename-based STL detection.
  std::ifstream input(path.string().c_str(), std::ios::in | std::ios::binary);
  if (!input.is_open())
    return false;

  OpenMesh::IO::Options options;
  if (is_binary)
    options += OpenMesh::IO::Options::Binary;
  return OpenMesh::IO::read_mesh(mesh, input, ".stl", options);
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

// Short label of a mode value for the run folder name, e.g.
// topological_offset -> TO.  Values without a fixed label use the initials of
// their words, which also keeps the label a valid path component.
std::string abbreviate_mode(const std::string& value)
{
  static const std::map<std::string, std::string> labels = {
    { "topological_offset", "TO" },
    { "linear_solve", "LS" },
    { "newton_solve", "NS" },
    { "qem_original", "QEM" },
    { "fast", "FS" },
  };
  const auto found = labels.find(value);
  if (found != labels.end())
    return found->second;

  std::string initials;
  bool word_start = true;
  for (char ch : value)
  {
    const unsigned char uch = static_cast<unsigned char>(ch);
    if (!std::isalnum(uch))
    {
      word_start = true;
      continue;
    }
    if (word_start)
      initials += static_cast<char>(std::toupper(uch));
    word_start = false;
  }
  return initials.empty() ? "unknown" : initials;
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

// <timestamp>__<labels joined by '_'>, e.g. 20260919_182543__phase1_TO_phase3_LS.
// COMMAND_ARGUMENTS.md lists the abbreviations.
std::string build_run_dir_name(
  const Cage::ParamCageGenerator& param, bool from_cage, bool from_rails)
{
  const auto& simplifier = param.paramCageSimplifier;
  const auto& collapse = simplifier.paramCollapse;
  std::vector<std::string> labels;
  // The default subdivision Phase 1 is not named.
  if (param.paramCageInitializer.phase1Mode != "subdivision")
    labels.push_back("phase1_" + abbreviate_mode(param.paramCageInitializer.phase1Mode));
  if (from_rails)
    labels.push_back("from_RC");
  else if (from_cage)
    labels.push_back("from_IC");
  if (simplifier.boundaryRailAnchorMode == "compare")
  {
    labels.push_back("BRC");
  }
  else
  {
    labels.push_back("phase3_" + abbreviate_mode(simplifier.phase3Mode));
    if (simplifier.phase3Mode == "newton_solve")
    {
      labels.push_back(abbreviate_mode(collapse.phase3PlacementStrategy));
      labels.push_back(abbreviate_mode(collapse.newtonSolverMode));
      labels.push_back(abbreviate_mode(collapse.robustnessMode));
      if (collapse.curvatureMode != "none")
        labels.push_back("curv" + abbreviate_mode(collapse.curvatureMode));
      if (collapse.uniformityMode != "none")
        labels.push_back("unif" + abbreviate_mode(collapse.uniformityMode));
    }
    else if (simplifier.phase3Mode != "linear_solve" &&
      simplifier.phase3Mode != "qem_original")
    {
      labels.push_back("colH");
      labels.push_back("flip" + abbreviate_mode(simplifier.paramFlip.priorityMode));
      labels.push_back("reloc" + abbreviate_mode(simplifier.paramRelocate.priorityMode));
    }
    // Only the enabled cases are named.
    if (simplifier.enableBoundaryRails && simplifier.enableRailUpdate)
      labels.push_back("RU");
    if (simplifier.enableBoundaryRails && simplifier.enableRailSupport)
      labels.push_back("RS");
  }

  std::ostringstream oss;
  oss << format_timestamp(std::time(nullptr));
  for (size_t i = 0; i < labels.size(); i++)
    oss << (i == 0 ? "__" : "_") << labels[i];
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

std::string csv_cell(std::string value)
{
  size_t position = 0;
  while ((position = value.find('"', position)) != std::string::npos)
  {
    value.insert(position, 1, '"');
    position += 2;
  }
  return '"' + value + '"';
}

std::string boundary_rail_benchmark_winner(
  const Cage::CageInit::BoundaryRailBuildStats& edge,
  const Cage::CageInit::BoundaryRailBuildStats& vertex)
{
  if (edge.builtLoopCount > vertex.builtLoopCount)
    return "edge";
  if (vertex.builtLoopCount > edge.builtLoopCount)
    return "vertex";
  return "tie";
}

bool write_boundary_rail_benchmark_csv(
  const bf::path& path,
  const std::string& model_name,
  const Cage::CageInit::BoundaryRailBuildStats& edge,
  const Cage::CageInit::BoundaryRailBuildStats& vertex,
  double phase1_seconds,
  double edge_seconds,
  double vertex_seconds)
{
  std::ofstream csv(path.string().c_str());
  if (!csv.is_open())
  {
    Logger::user_logger->error(
      "fail to open boundary-rail benchmark CSV: {}", path.string());
    return false;
  }

  csv
    << "model,edge_detected_loops,edge_ray_valid_loops,edge_built_loops,"
    << "edge_anchor_requests,edge_rail_vertices,edge_rail_edges,edge_successful,"
    << "vertex_detected_loops,vertex_ray_valid_loops,vertex_built_loops,"
    << "vertex_anchor_requests,vertex_rail_vertices,vertex_rail_edges,"
    << "vertex_successful,winner,phase1_seconds,edge_seconds,vertex_seconds\n";
  csv << std::setprecision(17)
    << csv_cell(model_name) << ','
    << edge.detectedLoopCount << ','
    << edge.rayValidLoopCount << ','
    << edge.builtLoopCount << ','
    << edge.anchorRequestCount << ','
    << edge.railVertexCount << ','
    << edge.railEdgeCount << ','
    << (edge.successful ? 1 : 0) << ','
    << vertex.detectedLoopCount << ','
    << vertex.rayValidLoopCount << ','
    << vertex.builtLoopCount << ','
    << vertex.anchorRequestCount << ','
    << vertex.railVertexCount << ','
    << vertex.railEdgeCount << ','
    << (vertex.successful ? 1 : 0) << ','
    << boundary_rail_benchmark_winner(edge, vertex) << ','
    << phase1_seconds << ','
    << edge_seconds << ','
    << vertex_seconds << '\n';
  return true;
}

// Phase 2 outputs of both anchor modes, named like the regular Phase 2 cage
// and rails with an _edge/_vertex suffix.  stageInitialize already wrote the
// shared Phase 1 cage.
void write_boundary_rail_benchmark_meshes(
  const CageGenerator& cage_generator,
  Cage::SM::SMeshT& edge_source,
  Cage::SM::SMeshT& edge_cage,
  Cage::SM::SMeshT& vertex_source,
  Cage::SM::SMeshT& vertex_cage)
{
  const auto write_mode = [&](Cage::SM::SMeshT& source, Cage::SM::SMeshT& cage,
    const std::string& mode)
  {
    // Full double precision, so each cage/rail pair can be passed back with
    // --cage/--rails.
    const std::string cage_name =
      cage_generator.stageOutputName("phase2_cage") + "_" + mode + ".obj";
    const std::string cage_path =
      cage_generator.stageOutputPath("phase2_cage") + "_" + mode + ".obj";
    if (!Cage::CageInit::write_cage_obj(cage, cage_path))
    {
      Logger::user_logger->warn(
        "fail to write boundary-rail benchmark cage: {}", cage_path);
    }
    const std::string rail_path =
      cage_generator.stageOutputPath("phase2_rails") + "_" + mode;
    write_boundary_rail_files(
      cage, rail_path + ".obj", rail_path + ".txt", cage_name, &source);
  };

  write_mode(edge_source, edge_cage, "edge");
  write_mode(vertex_source, vertex_cage, "vertex");
}

void run_boundary_rail_anchor_benchmark(
  CageGenerator& cage_generator,
  const bf::path& file_out_dir,
  const std::string& file_name)
{
  // Computation times, as in CageGenerator::generate (see PhaseTimer).
  PhaseTimer phase1_timer;
  phase1_timer.start();
  if (cage_generator.inputCagePath.empty())
    cage_generator.stageInitialize();
  else
  {
    PhaseTimer::Exclusion exclusion("input");
    cage_generator.stageLoadInitialCage();
  }
  phase1_timer.stop();
  const double phase1_seconds = phase1_timer.seconds();

  // Both alternatives start from independent deep copies of exactly the same
  // source and Phase 1 cage.  This isolates the anchor choice from Phase 1 and
  // avoids counting initialization twice.
  Cage::SM::SMeshT edge_source(*cage_generator.originalMesh);
  Cage::SM::SMeshT vertex_source(*cage_generator.originalMesh);
  Cage::SM::SMeshT edge_cage(*cage_generator.cage);
  Cage::SM::SMeshT vertex_cage(*cage_generator.cage);

  Cage::CageInit::BoundaryRailSearchOptions search_options;
  const auto& simplifier_options = cage_generator.param.paramCageSimplifier;
  search_options.candidateLimit = simplifier_options.boundaryRailCandidateLimit;
  search_options.retryCyclicStarts = simplifier_options.boundaryRailRetryCyclicStarts;
  search_options.maxSeconds = simplifier_options.boundaryRailSearchSeconds;

  PhaseTimer edge_timer;
  edge_timer.start();
  Cage::CageInit::BoundaryRailBuilder edge_builder(
    &edge_source, &edge_cage,
    Cage::CageInit::BoundaryRailAnchorMode::EdgeMidpoint, search_options);
  const Cage::CageInit::BoundaryRailBuildStats edge =
    edge_builder.build_with_stats();
  edge_timer.stop();
  const double edge_seconds = edge_timer.seconds();

  PhaseTimer vertex_timer;
  vertex_timer.start();
  Cage::CageInit::BoundaryRailBuilder vertex_builder(
    &vertex_source, &vertex_cage,
    Cage::CageInit::BoundaryRailAnchorMode::VertexBisector, search_options);
  const Cage::CageInit::BoundaryRailBuildStats vertex =
    vertex_builder.build_with_stats();
  vertex_timer.stop();
  const double vertex_seconds = vertex_timer.seconds();

  const std::string winner = boundary_rail_benchmark_winner(edge, vertex);
  Logger::user_logger->info(
    "boundary rail anchor benchmark: edge built {}/{} loops ({} ray-valid); vertex built {}/{} loops ({} ray-valid); winner {}.",
    edge.builtLoopCount, edge.detectedLoopCount, edge.rayValidLoopCount,
    vertex.builtLoopCount, vertex.detectedLoopCount,
    vertex.rayValidLoopCount, winner);
  Logger::user_logger->info(
    "boundary rail anchor benchmark elapsed time: Phase 1 {:.6f}s, Phase 2 edge {:.6f}s, Phase 2 vertex {:.6f}s; Phase 3 skipped.",
    phase1_seconds, edge_seconds, vertex_seconds);

  write_boundary_rail_benchmark_meshes(
    cage_generator, edge_source, edge_cage, vertex_source, vertex_cage);
  bf::path csv_path = file_out_dir;
  csv_path.append("boundary_rail_anchor_benchmark.csv");
  if (write_boundary_rail_benchmark_csv(
      csv_path, file_name, edge, vertex,
      phase1_seconds, edge_seconds, vertex_seconds))
  {
    Logger::user_logger->info(
      "wrote boundary rail anchor benchmark CSV to {}.", csv_path.string());
  }
}

void generate_cages(
  Cage::ParamCageGenerator param,
  bf::path in_model_path,
  bf::path out_data_path,
  std::vector<size_t> target_vn,
  const std::string& input_cage_path,
  const std::string& input_rail_path
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
  cage_generator.inputCagePath = input_cage_path;
  cage_generator.inputRailPath = input_rail_path;
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
      build_run_dir_name(
        param, !input_cage_path.empty(), !input_rail_path.empty()));
    // create log file
    bf::path log_path = file_out_dir;
    log_path.append("log.txt");
    Logger::updateFileLog(true, spdlog::level::trace, log_path.string());
    Logger::user_logger->info("processing {}", file_name);
    Logger::user_logger->info("output directory: {}", file_out_dir.string());
    // read input mesh
    std::string input_format;
    if (!read_input_mesh(*cage_generator.originalMesh, in_model_path, input_format))
    {
      Logger::user_logger->warn("fail to read input mesh: {}", file_path);
      throw logic_error("fail to read input mesh");
    }
    Logger::user_logger->info("input format: {}.", input_format);
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

    if (cage_generator.param.paramCageSimplifier.boundaryRailAnchorMode ==
      "compare")
    {
      if (target_vn.size() > 1)
      {
        Logger::user_logger->warn(
          "phase2_boundary_rail_compare ignores additional nested-cage targets; the benchmark uses one shared Phase 1 cage and skips Phase 3.");
      }
      cage_generator.param.setCageLabel(0);
      cage_generator.param.setTargetNumber(target_vn.front());
      run_boundary_rail_anchor_benchmark(
        cage_generator, file_out_dir, file_name);
      Logger::dev_logger->flush();
      return;
    }

    // do cage or nested cages generation
    for (size_t it = 0;it < target_vn.size();it++)
    {
      cage_generator.param.setCageLabel(it);
      cage_generator.param.setTargetNumber(target_vn[it]);
      cage_generator.generate();
      const std::string cage_name = cage_generator.stageOutputName("phase3_cage");
      log_triangle_quality_stats(*cage_generator.cage, cage_name);

      const std::string cage_path =
        cage_generator.stageOutputPath("phase3_cage") + ".obj";
      // Double precision: OpenMesh's writer would round the cage to float and
      // could move it into the source it was checked against.
      Cage::CageInit::write_cage_obj(*cage_generator.cage, cage_path);

      // export the final boundary rails on their own for separate
      // visualization; the rails Phase 2 handed to Phase 3 were already written
      // as <file>_phase2_rails.obj/.txt.
      const std::string rail_path = cage_generator.stageOutputPath("phase3_rails");
      const Cage::CageInit::BoundaryRailExport rail_export =
        write_boundary_rail_files(
          *cage_generator.cage, rail_path + ".obj", rail_path + ".txt",
          cage_name + ".obj", cage_generator.originalMesh.get());
      if (rail_export.rail_count > 0)
      {
        Logger::user_logger->info(
          "wrote {} Phase 3 boundary rails ({} rail vertices, {} rail edges) to {}.obj and .txt.",
          rail_export.rail_count, rail_export.vertex_count, rail_export.edge_count,
          rail_path);
      }
      else if (cage_generator.param.paramCageSimplifier.enableBoundaryRails)
      {
        Logger::user_logger->warn(
          "boundary rails were enabled but the cage carries no rail vertex or edge; no rail file written.");
      }

      *cage_generator.originalMesh = *cage_generator.cage;
      // The input cage and rails belong to cage 0; each nested cage runs
      // Phases 1 and 2 on the previous cage.
      cage_generator.inputCagePath.clear();
      cage_generator.inputRailPath.clear();
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
  catch (const std::exception& e)
  {
    Logger::user_logger->error("unexpected exception: {}", e.what());
    Logger::dev_logger->error("unexpected exception: {}", e.what());
  }
  catch (...)
  {
    Logger::dev_logger->error("unexpected exception of unknown type, check it!");
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

void apply_plane_energy_weights(
  Cage::ParamCollapseStage& collapse, bool include_triangle_quality)
{
  collapse.planeWeight = 1.0;
  collapse.triangleQualityWeight = include_triangle_quality ? 2.0 : 0.0;
  collapse.uniformityWeight = 4.0;
}

void enable_newton_solve_phase3_defaults(Cage::ParamCageGenerator& param)
{
  auto& simplifier = param.paramCageSimplifier;
  auto& collapse = simplifier.paramCollapse;
  simplifier.phase3Mode = "newton_solve";
  collapse.collapsePlacementMethod = "optimization";
  collapse.phase3PlacementStrategy = "newton_solve";
  collapse.curvatureMode = "weighted_plane";
  collapse.uniformityMode = "source";
  collapse.triangleQualityWeight = 2.0;
  collapse.uniformityWeight = 1.0;
  param.paramCageSimplifier.paramFlip.priorityMode = "triangle_quality_hard";
  param.paramCageSimplifier.paramFlip.requireRegularValence = false;
}

void enable_linear_solve_phase3_defaults(Cage::ParamCageGenerator& param)
{
  auto& simplifier = param.paramCageSimplifier;
  auto& collapse = simplifier.paramCollapse;
  simplifier.phase3Mode = "linear_solve";
  collapse.collapsePlacementMethod = "optimization";
  collapse.phase3PlacementStrategy = "linear_solve";
  collapse.robustnessMode = "exact_reject";
  collapse.curvatureMode = "none";
  collapse.uniformityMode = "global";
  apply_plane_energy_weights(collapse, true);
}

void enable_qem_original_phase3_defaults(Cage::ParamCageGenerator& param)
{
  auto& simplifier = param.paramCageSimplifier;
  auto& collapse = simplifier.paramCollapse;
  simplifier.phase3Mode = "qem_original";
  collapse.collapsePlacementMethod = "optimization";
  collapse.phase3PlacementStrategy = "qem_original";
  collapse.robustnessMode = "exact_reject";
  collapse.curvatureMode = "none";
  collapse.uniformityMode = "none";
  apply_plane_energy_weights(collapse, false);
}

bool apply_parameter_token(const std::string& raw_token, Cage::ParamCageGenerator& param)
{
  const std::string token = normalize_parameter_token(raw_token);
  auto& collapse = param.paramCageSimplifier.paramCollapse;

  if (token.empty() || token == "default")
    return true;

  if (token == "phase1_topological_offset" || token == "topological_offset")
  {
    param.paramCageInitializer.phase1Mode = "topological_offset";
    return true;
  }

  if (token == "phase2_boundary_rail" || token == "boundary_rail")
  {
    param.paramCageSimplifier.enableBoundaryRails = true;
    param.paramCageSimplifier.boundaryRailAnchorMode = "edge";
    return true;
  }
  if (token == "phase2_boundary_rail_vertex" || token == "boundary_rail_vertex")
  {
    param.paramCageSimplifier.enableBoundaryRails = true;
    param.paramCageSimplifier.boundaryRailAnchorMode = "vertex";
    return true;
  }
  if (token == "phase2_boundary_rail_compare" || token == "boundary_rail_compare")
  {
    param.paramCageSimplifier.enableBoundaryRails = true;
    param.paramCageSimplifier.boundaryRailAnchorMode = "compare";
    return true;
  }

  if (token == "phase3_rail_update" || token == "rail_update")
  {
    param.paramCageSimplifier.enableRailUpdate = true;
    return true;
  }
  if (token == "phase3_rail_support" || token == "rail_support")
  {
    param.paramCageSimplifier.enableRailSupport = true;
    return true;
  }

  if (token == "phase3_linear_solve" || token == "linear_solve")
  {
    enable_linear_solve_phase3_defaults(param);
    return true;
  }
  if (token == "phase3_linear_solve_collision_reject" ||
    token == "linear_solve_collision_reject")
  {
    enable_linear_solve_phase3_defaults(param);
    collapse.phase3LinearSolveCollisionReject = true;
    return true;
  }
  if (token == "phase3_qem_original" || token == "qem_original")
  {
    enable_qem_original_phase3_defaults(param);
    return true;
  }
  if (token == "phase3_newton_solve" || token == "newton_solve")
  {
    enable_newton_solve_phase3_defaults(param);
    apply_plane_energy_weights(collapse, true);
    return true;
  }
  Logger::user_logger->error("unknown parameter token: {}", raw_token);
  return false;
}

bool parse_parameter_arg(const std::string& arg_param, Cage::ParamCageGenerator& param)
{
  const auto phase1_mode_is_valid = [&]()
  {
    const std::string& mode = param.paramCageInitializer.phase1Mode;
    if (mode == "subdivision" || mode == "topological_offset")
      return true;
    Logger::user_logger->error(
      "invalid phase1Mode: {} (expected subdivision or topological_offset).",
      mode);
    return false;
  };

  const auto boundary_rail_mode_is_valid = [&]()
  {
    if (!param.paramCageSimplifier.enableBoundaryRails)
      return true;
    const std::string& anchor_mode =
      param.paramCageSimplifier.boundaryRailAnchorMode;
    if (anchor_mode != "edge" && anchor_mode != "vertex" &&
      anchor_mode != "compare")
    {
      Logger::user_logger->error(
        "invalid boundaryRailAnchorMode: {} (expected edge, vertex, or compare).",
        anchor_mode);
      return false;
    }
    if (anchor_mode == "compare")
      return true;
    const std::string& mode = param.paramCageSimplifier.phase3Mode;
    if (mode == "linear_solve" || mode == "newton_solve" ||
      mode == "qem_original")
      return true;
    Logger::user_logger->error(
      "phase2_boundary_rail and phase2_boundary_rail_vertex must be combined with phase3_linear_solve, phase3_linear_solve_collision_reject, phase3_newton_solve, or phase3_qem_original; phase2_boundary_rail_compare skips Phase 3.");
    return false;
  };

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
    return phase1_mode_is_valid() && boundary_rail_mode_is_valid();
  }

  std::vector<std::string> parameter_tokens;
  boost::split(parameter_tokens, arg_param, boost::is_any_of("+,"), boost::token_compress_on);
  for (const std::string& token : parameter_tokens)
  {
    if (!apply_parameter_token(token, param))
      return false;
  }
  return phase1_mode_is_valid() && boundary_rail_mode_is_valid();
}

int main(int argc, char* argv[])
{
  // Options of an earlier run may appear anywhere; the remaining arguments
  // keep their positional meaning.
  std::vector<std::string> args;
  std::string input_cage_path;
  std::string input_rail_path;
  for (int i = 1; i < argc; i++)
  {
    const std::string arg(argv[i]);
    if (arg == "--cage" || arg == "--rails")
    {
      if (i + 1 >= argc)
      {
        printf("missing path after %s\n", arg.c_str());
        return 1;
      }
      (arg == "--cage" ? input_cage_path : input_rail_path) = argv[++i];
    }
    else
      args.push_back(arg);
  }

  // args tips
  if (args.size() < 3)
  {
    printf("Need args:\n");
    printf("arg[0]: parameters.\n");
    printf("input \"default\" to set default parameters\n");
    printf("Phase 1 (initial cage):\n");
    printf("input \"phase1_topological_offset\" to use simplicial embedding and offset insertion\n");
    printf("Phase 2 (boundary rail construction, runs only when requested):\n");
    printf("input \"phase2_boundary_rail\" with a Phase 3 energy preset to build closed boundary rails and preserve them in Phase 3\n");
    printf("input \"phase2_boundary_rail_vertex\" with a Phase 3 energy preset to use vertex-bisector anchor rays\n");
    printf("input \"phase2_boundary_rail_compare\" to compare edge and vertex anchor rays on one shared Phase 1 cage and skip Phase 3\n");
    printf("Phase 3 (simplification):\n");
    printf("input \"phase3_linear_solve\" for repeated linear-solve collapse and quality flip, followed by final Voronoi/surface relocation sweeps\n");
    printf("input \"phase3_linear_solve_collision_reject\" to reject invalid raw linear-solve placements without backtracking\n");
    printf("input \"phase3_newton_solve\" to use Newton placement for every Phase 3 collapse candidate\n");
    printf("input \"phase3_qem_original\" to run original Garland-Heckbert QEM without collision rejection\n");
    printf("input \"phase3_rail_update\" with Phase 2 rails to relabel rails after each linear-solve flip stage (off by default)\n");
    printf("input \"phase3_rail_support\" with Phase 2 rails to project rail-edge collapses onto the source-boundary half-strips (off by default)\n");
    printf("the phaseN_ prefix may be omitted, e.g. \"linear_solve\" or \"boundary_rail\".\n");
    printf("combine presets with '+' or ',' in any order, for example \"phase1_topological_offset+phase2_boundary_rail+phase3_linear_solve\".\n");
    printf("or a json file to set parameters.\n");
    printf("arg[1]: input model path.\n");
    printf("ASCII and binary STL input files are supported.\n");
    printf("arg[2]: output dir path.\n");
    printf("arg[3]: target vertices number of cage 0 (optional for linear_solve).\n");
    printf("(optional)arg[4]: target vertices number of nested cage 1.\n");
    printf("(optional)arg[n + 3]: target vertices number of nested cage n.\n");
    printf("For linear_solve, omit targets or use 0 to repeat collapse/flip (and rail update with phase3_rail_update) until no further progress.\n");
    printf("options:\n");
    printf("--cage <cage.obj>: skip Phase 1 and start from this cage, e.g. <model>_phase1_cage.obj of an earlier run.\n");
    printf("--rails <rails.txt>: with --cage <model>_phase2_cage.obj, also skip Phase 2 and use <model>_phase2_rails.txt.\n");
    return 1;
  }

  // logger and random seed.
  Logger::InitLogger(spdlog::level::level_enum::trace);
  srand((unsigned)time(NULL));

  // parse parameters
  const std::string& arg_param = args[0];
  Cage::ParamCageGenerator param;
  if (!parse_parameter_arg(arg_param, param))
    return 1;

  if (args.size() == 3 && param.paramCageSimplifier.phase3Mode != "linear_solve")
  {
    Logger::user_logger->error(
      "a target vertex count is required; targets may be omitted only for linear_solve.");
    return 1;
  }

  // parse input/output file/directory.
  bf::path in_model_path(args[1]);
  bf::path out_data_path(args[2]);

  if (!bf::is_regular_file(in_model_path))
  {
    Logger::user_logger->error("error input file: {}", in_model_path.string());
    return 1;
  }

  // parse the files of an earlier run.
  if (!input_cage_path.empty() && !bf::is_regular_file(input_cage_path))
  {
    Logger::user_logger->error("error initial cage file: {}", input_cage_path);
    return 1;
  }
  if (!input_rail_path.empty())
  {
    const auto& simplifier = param.paramCageSimplifier;
    if (input_cage_path.empty())
    {
      Logger::user_logger->error(
        "--rails needs --cage with the cage its indices refer to (<model>_phase2_cage.obj).");
      return 1;
    }
    if (!simplifier.enableBoundaryRails ||
      simplifier.boundaryRailAnchorMode == "compare")
    {
      Logger::user_logger->error(
        "--rails needs phase2_boundary_rail or phase2_boundary_rail_vertex; phase2_boundary_rail_compare always builds its rails.");
      return 1;
    }
    // A rail OBJ has no cage indices; its index file is written next to it.
    bf::path rail_path(input_rail_path);
    if (boost::algorithm::iequals(rail_path.extension().string(), ".obj"))
    {
      rail_path.replace_extension(".txt");
      Logger::user_logger->info(
        "--rails {} is a rail OBJ; reading its index file {} instead.",
        input_rail_path, rail_path.string());
      input_rail_path = rail_path.string();
    }
    if (!bf::is_regular_file(rail_path))
    {
      Logger::user_logger->error("error rail file: {}", input_rail_path);
      return 1;
    }
  }
  else if (!input_cage_path.empty() &&
    param.paramCageSimplifier.enableBoundaryRails)
  {
    Logger::user_logger->info(
      "boundary rails will be constructed on the loaded cage; to reuse the rails of a <model>_phase2_cage.obj, pass its rail file with --rails.");
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
  size_t cage_number = args.size() - 3;
  std::vector<size_t> target_vn;
  try
  {
    for (size_t i = 0;i < cage_number;i++)
    {
      std::string vn_str(args[3 + i]);
      size_t vn = std::stol(vn_str);
      target_vn.push_back(vn);
    }
  }
  catch (...)
  {
    Logger::user_logger->error("error in parsing arguments.");
    return 1;
  }

  if (target_vn.empty())
    target_vn.push_back(0);

  // generate!
  generate_cages(
    param, in_model_path, out_data_path, target_vn,
    input_cage_path, input_rail_path);
  return 0;
}
