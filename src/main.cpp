#include "CageGenerator.hh"
#include "boost/filesystem.hpp"
#include "boost/algorithm/string.hpp"
#include <cfloat>
#include <cmath>

using namespace Cage;
namespace bf = boost::filesystem;
namespace bj = boost::json;

bool mesh_valid(Cage::SMeshT& mesh)
{
  // test close
  /*for (auto eh : mesh.edges())
  {
    if (mesh.is_boundary(eh))
      return false;
  }*/
  // test manifold
  for (auto vh : mesh.vertices())
  {
    if (!mesh.is_manifold(vh))
      return false;
  }
  return true;
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
    bf::path file_out_dir = out_data_path;
    file_out_dir.append(file_name);
    if (!bf::exists(file_out_dir))
      bf::create_directory(file_out_dir);
    // create log file
    bf::path log_path = file_out_dir;
    log_path.append("log.txt");
    Logger::updateFileLog(true, spdlog::level::trace, log_path.string());
    Logger::user_logger->info("processing {}", file_name);
    // read input mesh
    OpenMesh::IO::read_mesh(*cage_generator.originalMesh, file_path);
    // check input
    if (!mesh_valid(*cage_generator.originalMesh))
    {
      Logger::user_logger->warn("invalid mesh.");
      throw logic_error("invalid mesh");
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
  boost::replace_all(token, "_", "-");
  if (token.rfind("default-", 0) == 0)
    token = token.substr(std::string("default-").size());
  return token;
}

bool apply_parameter_token(const std::string& raw_token, Cage::ParamCageGenerator& param)
{
  const std::string token = normalize_parameter_token(raw_token);
  auto& collapse = param.paramCageSimplifier.paramCollapse;
  auto& flip = param.paramCageSimplifier.paramFlip;

  if (token.empty() || token == "default")
    return true;

  if (token == "length")
  {
    collapse.priorityMode = "length";
    return true;
  }
  if (token == "length-quality" || token == "length-quality-weighted")
  {
    collapse.priorityMode = "length_quality";
    collapse.lengthQualitySubMode = "weighted";
    return true;
  }
  if (token == "length-quality-relative")
  {
    collapse.priorityMode = "length_quality";
    collapse.lengthQualitySubMode = "relative_reject";
    return true;
  }
  if (token == "length-quality-absolute")
  {
    collapse.priorityMode = "length_quality";
    collapse.lengthQualitySubMode = "absolute_reject";
    return true;
  }
  if (token == "length-quality-lexicographic")
  {
    collapse.priorityMode = "length_quality";
    collapse.lengthQualitySubMode = "lexicographic";
    return true;
  }
  if (token == "post-edge-length")
  {
    collapse.priorityMode = "post_edge_length";
    return true;
  }
  if (token == "post-edge-length-hard")
  {
    collapse.priorityMode = "post_edge_length_hard";
    return true;
  }
  if (token == "post-face-area")
  {
    collapse.priorityMode = "post_face_area";
    return true;
  }
  if (token == "post-face-area-hard")
  {
    collapse.priorityMode = "post_face_area_hard";
    return true;
  }
  if (token == "triangle-quality")
  {
    collapse.priorityMode = "triangle_quality";
    return true;
  }
  if (token == "triangle-quality-hard")
  {
    collapse.priorityMode = "triangle_quality_hard";
    return true;
  }
  if (token == "flip-valence")
  {
    flip.priorityMode = "valence";
    return true;
  }
  if (token == "flip-triangle-quality-hard")
  {
    flip.priorityMode = "triangle_quality_hard";
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
    printf("input \"default-length\" to set default parameters with length-based collapse priority\n");
    printf("input \"default-length-quality\" to use weighted length + triangle quality priority\n");
    printf("input \"default-length-quality-weighted\" to use weighted quality penalty\n");
    printf("input \"default-length-quality-relative\" to reject large relative quality drops\n");
    printf("input \"default-length-quality-absolute\" to reject triangles below a quality threshold\n");
    printf("input \"default-length-quality-lexicographic\" to use length first and quality as tie-break\n");
    printf("input \"default-post-edge-length\" to use post-collapse edge length priority\n");
    printf("input \"default-post-edge-length-hard\" to require post-collapse max edge length not to increase\n");
    printf("input \"default-post-face-area\" to use post-collapse face area priority\n");
    printf("input \"default-post-face-area-hard\" to require post-collapse max face area not to increase\n");
    printf("input \"default-triangle-quality\" to use triangle quality priority\n");
    printf("input \"default-triangle-quality-hard\" to require post-collapse min triangle quality not to decrease\n");
    printf("input \"flip-triangle-quality-hard\" to require post-flip min triangle quality not to decrease\n");
    printf("combine presets with '+' or ',', for example \"default-length+flip-triangle-quality-hard\".\n");
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
    Logger::user_logger->error("error input file.");
    return 1;
  }

  if (!bf::is_directory(out_data_path))
  {
    Logger::user_logger->error("error output directory.");
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
