#include "CageGenerator.hh"
#include "boost/filesystem.hpp"
#include "boost/algorithm/string.hpp"

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

  if (arg_param == "default-length" || arg_param == "default_length")
  {
    param.paramCageSimplifier.paramCollapse.priorityMode = "length";
  }
  else if (arg_param == "default-length-quality" || arg_param == "default_length_quality")
  {
    param.paramCageSimplifier.paramCollapse.priorityMode = "length_quality";
    param.paramCageSimplifier.paramCollapse.lengthQualitySubMode = "weighted";
  }
  else if (arg_param == "default-length-quality-weighted" || arg_param == "default_length_quality_weighted")
  {
    param.paramCageSimplifier.paramCollapse.priorityMode = "length_quality";
    param.paramCageSimplifier.paramCollapse.lengthQualitySubMode = "weighted";
  }
  else if (arg_param == "default-length-quality-relative" || arg_param == "default_length_quality_relative")
  {
    param.paramCageSimplifier.paramCollapse.priorityMode = "length_quality";
    param.paramCageSimplifier.paramCollapse.lengthQualitySubMode = "relative_reject";
  }
  else if (arg_param == "default-length-quality-absolute" || arg_param == "default_length_quality_absolute")
  {
    param.paramCageSimplifier.paramCollapse.priorityMode = "length_quality";
    param.paramCageSimplifier.paramCollapse.lengthQualitySubMode = "absolute_reject";
  }
  else if (arg_param == "default-length-quality-lexicographic" || arg_param == "default_length_quality_lexicographic")
  {
    param.paramCageSimplifier.paramCollapse.priorityMode = "length_quality";
    param.paramCageSimplifier.paramCollapse.lengthQualitySubMode = "lexicographic";
  }
  else if (arg_param == "default-post-edge-length" || arg_param == "default_post_edge_length")
  {
    param.paramCageSimplifier.paramCollapse.priorityMode = "post_edge_length";
  }
  else if (arg_param == "default-post-edge-length-hard" || arg_param == "default_post_edge_length_hard")
  {
    param.paramCageSimplifier.paramCollapse.priorityMode = "post_edge_length_hard";
  }
  else if (arg_param == "default-post-face-area" || arg_param == "default_post_face_area")
  {
    param.paramCageSimplifier.paramCollapse.priorityMode = "post_face_area";
  }
  else if (arg_param == "default-post-face-area-hard" || arg_param == "default_post_face_area_hard")
  {
    param.paramCageSimplifier.paramCollapse.priorityMode = "post_face_area_hard";
  }
  else if (arg_param == "default-triangle-quality" || arg_param == "default_triangle_quality")
  {
    param.paramCageSimplifier.paramCollapse.priorityMode = "triangle_quality";
  }
  else if (arg_param == "default-triangle-quality-hard" || arg_param == "default_triangle_quality_hard")
  {
    param.paramCageSimplifier.paramCollapse.priorityMode = "triangle_quality_hard";
  }
  else if (arg_param != "default")
  {
    bf::path json_file_path(argv[1]);
    if (bf::is_regular_file(json_file_path))
    {
      fstream json_file;
      json_file.open(json_file_path.string(), fstream::in);
      if (json_file.is_open())
      {
        std::string json_str((std::istreambuf_iterator<char>(json_file)), std::istreambuf_iterator<char>());
        bj::stream_parser sp;
        sp.write(json_str.c_str());
        param.deserialize(sp.release().as_object());
        json_file.close();
      }
      else
      {
        Logger::user_logger->error("fail to open json file.");
        return 1;
      }
    }
    else
    {
      Logger::user_logger->error("wrong json file.");
      return 1;
    }
  }

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
