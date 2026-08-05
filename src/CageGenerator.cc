#include "CageGenerator.hh"
#include <chrono>
#include "omp.h"

namespace Cage
{
void CageGenerator::stageInitialize()
{
  VMesh = std::make_unique<VM::VMeshT>();
  cage = std::make_unique<SM::SMeshT>();

  cageInitializer = std::make_unique<CageInitializer>(
    originalMesh.get(), &param.paramCageInitializer,
    VMesh.get(), cage.get());

  cageInitializer->generate();
}

void CageGenerator::stageSimplify()
{
  cageSimplifier = std::make_unique<CageSimplifier>(
    originalMesh.get(), cage.get(), &param.paramCageSimplifier);

  cageSimplifier->simplify();
}

void CageGenerator::generate()
{
  omp_set_num_threads(12);

  const auto phase1_start = std::chrono::steady_clock::now();
  stageInitialize();
  const double phase1_seconds = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - phase1_start).count();

  cageInitializer = nullptr;
  VMesh = nullptr;

  const auto phase2_start = std::chrono::steady_clock::now();
  stageSimplify();
  const double phase2_seconds = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - phase2_start).count();

  Logger::user_logger->info(
    "Phase 1 elapsed time: {:.6f} seconds.", phase1_seconds);
  Logger::user_logger->info(
    "Phase 2 elapsed time: {:.6f} seconds.", phase2_seconds);
}

}// namespace Cage
