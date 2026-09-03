#include "CageGenerator.hh"
#include "boost/filesystem.hpp"
#include "boost/algorithm/string.hpp"
#include <cctype>
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

struct BoundaryRailExport
{
  size_t rail_count = 0;
  size_t vertex_count = 0;
  size_t edge_count = 0;
};

// Write the boundary rails of a finished cage into their own files so that the
// rails can be visualized separately from the cage surface.
//   <cage>_rails.obj : rail vertices as OBJ points and rail edges as OBJ line
//                      elements, one object/group per rail id.
//   <cage>_rails.txt : the same rails listed with the 1-based vertex indices of
//                      the cage OBJ, for cross referencing inside the cage.
BoundaryRailExport write_boundary_rail_files(
  Cage::SM::SMeshT& mesh,
  const bf::path& rail_obj_path,
  const bf::path& rail_txt_path,
  const std::string& cage_obj_name)
{
  typedef std::pair<Cage::SM::VertexHandle, Cage::SM::VertexHandle> RailEdge;

  BoundaryRailExport summary;
  std::map<int, std::vector<Cage::SM::VertexHandle>> rail_vertices;
  std::map<int, std::vector<RailEdge>> rail_edges;

  for (Cage::SM::VertexHandle vh : mesh.vertices())
  {
    const int rail_id = mesh.data(vh).boundary_rail_id;
    if (rail_id >= 0)
      rail_vertices[rail_id].push_back(vh);
  }
  for (Cage::SM::EdgeHandle eh : mesh.edges())
  {
    const int rail_id = mesh.data(eh).boundary_rail_id;
    if (rail_id < 0)
      continue;
    const Cage::SM::HalfedgeHandle hh = mesh.halfedge_handle(eh, 0);
    rail_edges[rail_id].push_back(
      RailEdge(mesh.from_vertex_handle(hh), mesh.to_vertex_handle(hh)));
  }

  if (rail_vertices.empty() && rail_edges.empty())
    return summary;

  // Local 1-based numbering of the rail OBJ, keyed by cage vertex index.  Edge
  // endpoints are registered too so that an edge is never dropped even if one
  // of its vertices lost its own rail id.
  std::map<int, size_t> local_index;
  std::vector<Cage::SM::VertexHandle> ordered_vertices;
  const auto register_vertex = [&](Cage::SM::VertexHandle vh)
  {
    if (local_index.emplace(vh.idx(), ordered_vertices.size() + 1).second)
      ordered_vertices.push_back(vh);
  };
  for (const auto& entry : rail_vertices)
  {
    for (Cage::SM::VertexHandle vh : entry.second)
      register_vertex(vh);
  }
  for (const auto& entry : rail_edges)
  {
    for (const RailEdge& edge : entry.second)
    {
      register_vertex(edge.first);
      register_vertex(edge.second);
    }
  }

  std::ofstream rail_obj(rail_obj_path.string().c_str());
  if (!rail_obj.is_open())
  {
    Logger::user_logger->warn(
      "fail to open boundary rail OBJ for writing: {}", rail_obj_path.string());
    return summary;
  }
  rail_obj << std::setprecision(15);
  rail_obj << "# boundary rails extracted from " << cage_obj_name << "\n";
  rail_obj << "# rail edges are OBJ line elements, rail vertices are OBJ points\n";
  rail_obj << "# vertices keep the order of " << rail_txt_path.filename().string() << "\n";
  for (Cage::SM::VertexHandle vh : ordered_vertices)
  {
    const Cage::SM::Vec3d& p = mesh.point(vh);
    rail_obj << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
  }
  for (const auto& entry : rail_edges)
  {
    const int rail_id = entry.first;
    rail_obj << "o rail_" << rail_id << "\n";
    rail_obj << "g rail_" << rail_id << "\n";
    for (const RailEdge& edge : entry.second)
    {
      rail_obj
        << "l " << local_index[edge.first.idx()]
        << " " << local_index[edge.second.idx()] << "\n";
    }
    const auto vertex_it = rail_vertices.find(rail_id);
    if (vertex_it != rail_vertices.end())
    {
      for (Cage::SM::VertexHandle vh : vertex_it->second)
        rail_obj << "p " << local_index[vh.idx()] << "\n";
    }
  }
  rail_obj.close();

  std::ofstream rail_txt(rail_txt_path.string().c_str());
  if (!rail_txt.is_open())
  {
    Logger::user_logger->warn(
      "fail to open boundary rail index file for writing: {}", rail_txt_path.string());
    return summary;
  }
  rail_txt << std::setprecision(15);
  rail_txt << "# boundary rails extracted from " << cage_obj_name << "\n";
  rail_txt << "# cage_vertex is the 1-based OBJ vertex index inside " << cage_obj_name << "\n";
  rail_txt << "# V <rail_id> <cage_vertex> <x> <y> <z>\n";
  rail_txt << "# E <rail_id> <cage_vertex_a> <cage_vertex_b>\n";
  for (const auto& entry : rail_vertices)
  {
    for (Cage::SM::VertexHandle vh : entry.second)
    {
      const Cage::SM::Vec3d& p = mesh.point(vh);
      rail_txt
        << "V " << entry.first << " " << (vh.idx() + 1)
        << " " << p[0] << " " << p[1] << " " << p[2] << "\n";
    }
    summary.vertex_count += entry.second.size();
  }
  for (const auto& entry : rail_edges)
  {
    for (const RailEdge& edge : entry.second)
    {
      rail_txt
        << "E " << entry.first
        << " " << (edge.first.idx() + 1)
        << " " << (edge.second.idx() + 1) << "\n";
    }
    summary.edge_count += entry.second.size();
  }
  rail_txt.close();

  std::set<int> rail_ids;
  for (const auto& entry : rail_vertices)
    rail_ids.insert(entry.first);
  for (const auto& entry : rail_edges)
    rail_ids.insert(entry.first);
  summary.rail_count = rail_ids.size();
  return summary;
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

std::string phase2_mode_label(const Cage::ParamCageSimplifier& param)
{
  return "phase2_" + sanitize_path_component(param.phase2Mode);
}

std::string newton_solve_phase2_detail_label(const Cage::ParamCollapseStage& param)
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
  oss << format_timestamp(std::time(nullptr));
  if (param.paramCageInitializer.phase1Mode != "subdivision")
  {
    oss << "__phase1_" <<
      sanitize_path_component(param.paramCageInitializer.phase1Mode);
  }
  oss << "__" << phase2_mode_label(simplifier);
  if (simplifier.phase2Mode == "newton_solve")
  {
    oss << "__" << newton_solve_phase2_detail_label(simplifier.paramCollapse);
  }
  else if (simplifier.phase2Mode != "linear_solve" &&
    simplifier.phase2Mode != "qem_original")
  {
    oss
      << "__collapse_hausdorff"
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

      // export the boundary rails on their own for separate visualization.
      bf::path rail_obj_file = file_out_dir;
      rail_obj_file.append(file_name + "_cage_" + std::to_string(it) + "_rails.obj");
      bf::path rail_txt_file = file_out_dir;
      rail_txt_file.append(file_name + "_cage_" + std::to_string(it) + "_rails.txt");
      const BoundaryRailExport rail_export = write_boundary_rail_files(
        *cage_generator.cage, rail_obj_file, rail_txt_file,
        mesh_out_file.filename().string());
      if (rail_export.rail_count > 0)
      {
        Logger::user_logger->info(
          "wrote {} boundary rails ({} rail vertices, {} rail edges) to {} and {}.",
          rail_export.rail_count, rail_export.vertex_count, rail_export.edge_count,
          rail_obj_file.string(), rail_txt_file.string());
      }
      else if (cage_generator.param.paramCageSimplifier.enableBoundaryRails)
      {
        Logger::user_logger->warn(
          "boundary rails were enabled but the cage carries no rail vertex or edge; no rail file written.");
      }

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

void apply_qem_energy_weights(
  Cage::ParamCollapseStage& collapse, bool include_triangle_quality)
{
  collapse.qemWeight = 1.0;
  collapse.triangleQualityWeight = include_triangle_quality ? 1.0 : 0.0;
  collapse.uniformityWeight = 10.0;
}

void enable_newton_solve_phase2_defaults(Cage::ParamCageGenerator& param)
{
  auto& simplifier = param.paramCageSimplifier;
  auto& collapse = simplifier.paramCollapse;
  simplifier.phase2Mode = "newton_solve";
  collapse.collapsePlacementMethod = "optimization";
  collapse.phase2PlacementStrategy = "newton_solve";
  collapse.curvatureMode = "weighted_qem";
  collapse.uniformityMode = "source";
  collapse.triangleQualityWeight = 2.0;
  collapse.uniformityWeight = 1.0;
  param.paramCageSimplifier.paramFlip.priorityMode = "triangle_quality_hard";
  param.paramCageSimplifier.paramFlip.requireRegularValence = false;
}

void enable_linear_solve_phase2_defaults(Cage::ParamCageGenerator& param)
{
  auto& simplifier = param.paramCageSimplifier;
  auto& collapse = simplifier.paramCollapse;
  simplifier.phase2Mode = "linear_solve";
  collapse.collapsePlacementMethod = "optimization";
  collapse.phase2PlacementStrategy = "linear_solve";
  collapse.robustnessMode = "exact_reject";
  collapse.curvatureMode = "none";
  collapse.uniformityMode = "global";
  apply_qem_energy_weights(collapse, true);
}

void enable_qem_original_phase2_defaults(Cage::ParamCageGenerator& param)
{
  auto& simplifier = param.paramCageSimplifier;
  auto& collapse = simplifier.paramCollapse;
  simplifier.phase2Mode = "qem_original";
  collapse.collapsePlacementMethod = "optimization";
  collapse.phase2PlacementStrategy = "qem_original";
  collapse.robustnessMode = "exact_reject";
  collapse.curvatureMode = "none";
  collapse.uniformityMode = "none";
  apply_qem_energy_weights(collapse, false);
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

  if (token == "boundary_rail")
  {
    param.paramCageSimplifier.enableBoundaryRails = true;
    return true;
  }

  if (token == "phase2_linear_solve" || token == "linear_solve")
  {
    enable_linear_solve_phase2_defaults(param);
    return true;
  }
  if (token == "phase2_linear_solve_collision_reject" ||
    token == "linear_solve_collision_reject")
  {
    enable_linear_solve_phase2_defaults(param);
    collapse.phase2LinearSolveCollisionReject = true;
    return true;
  }
  if (token == "phase2_qem_original" || token == "qem_original")
  {
    enable_qem_original_phase2_defaults(param);
    return true;
  }
  if (token == "phase2_newton_solve" || token == "newton_solve")
  {
    enable_newton_solve_phase2_defaults(param);
    apply_qem_energy_weights(collapse, true);
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
    const std::string& mode = param.paramCageSimplifier.phase2Mode;
    if (mode == "linear_solve" || mode == "newton_solve" ||
      mode == "qem_original")
      return true;
    Logger::user_logger->error(
      "boundary_rail must be combined with phase2_linear_solve, phase2_linear_solve_collision_reject, phase2_newton_solve, or phase2_qem_original.");
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
  // args tips
  if (argc < 5)
  {
    printf("Need args:\n");
    printf("arg[0]: parameters.\n");
    printf("input \"default\" to set default parameters\n");
    printf("input \"phase1_topological_offset\" to use simplicial embedding and offset insertion in Phase 1\n");
    printf("input \"phase2_linear_solve\" to run linear-solve Phase 2 collapses with Armijo backtracking and hard intersection constraints\n");
    printf("input \"phase2_linear_solve_collision_reject\" to reject invalid raw linear-solve placements without backtracking\n");
    printf("input \"phase2_newton_solve\" to use Newton placement for every Phase 2 collapse candidate\n");
    printf("input \"phase2_qem_original\" to run original Garland-Heckbert QEM without collision rejection\n");
    printf("input \"boundary_rail\" with a Phase 2 energy preset to build and preserve closed boundary rails\n");
    printf("combine presets with '+' or ',', for example \"phase2_linear_solve+boundary_rail\".\n");
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
