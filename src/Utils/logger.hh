#pragma once
#include <cstdio>
#include <chrono>
#include <iostream>
#include <map>
#include <memory>
#include <string>
// spdlog
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/stopwatch.h"

namespace Cage
{
namespace SimpleUtils
{

/****************************/
/******** Logger ************/
/****************************/

class Logger
{
public:
  static std::shared_ptr<spdlog::logger> dev_logger;
  static std::shared_ptr<spdlog::logger> user_logger;
  static std::shared_ptr<spdlog::stopwatch> sw;

  static void InitLogger(
    spdlog::level::level_enum console_level = spdlog::level::trace,
    bool file_log = false,
    spdlog::level::level_enum file_level = spdlog::level::trace,
    const std::string& file_path = "logs/log.txt"
  );

  static void updateFileLog(bool file_log, spdlog::level::level_enum file_level, const std::string& file_path);

  static void StopProgram();

  // Seconds all log sinks have spent formatting and writing messages so far.
  static double loggingSeconds();
};

/****************************/
/******** Phase time ********/
/****************************/

// Computation time of one pipeline phase: the wall time while it runs, minus
// the time the log sinks spent and minus the sections excluded with
// PhaseTimer::Exclusion (file output and verification that the cage does not
// need).
class PhaseTimer
{
public:
  using Clock = std::chrono::steady_clock;

  // Excludes its lifetime from the running timer, if any, under a label such
  // as "output" or "check".  Nested exclusions count once, under the
  // outermost label.
  class Exclusion
  {
  public:
    explicit Exclusion(const char* label);
    ~Exclusion();
    Exclusion(const Exclusion&) = delete;
    Exclusion& operator=(const Exclusion&) = delete;
  private:
    PhaseTimer* timer;
    const char* label;
    Clock::time_point start;
    double loggingAtStart;
  };

  // A started timer is the one exclusions apply to until it stops.
  void start();
  void stop();
  bool wasStarted() const { return everStarted; }
  // Wall time minus logging and exclusions.
  double seconds() const;
  double loggingSeconds() const { return logging; }
  const std::map<std::string, double>& excludedSeconds() const { return excluded; }

private:
  static PhaseTimer* active;
  static int exclusionDepth;
  bool running = false;
  bool everStarted = false;
  Clock::time_point startTime;
  double loggingAtStart = 0.0;
  double wall = 0.0;
  double logging = 0.0;
  std::map<std::string, double> excluded;
};

/****************************/
/******** Trace *************/
/****************************/

#if defined(TRACE_FUNC)
#define REP_FUNC {Logger::dev_logger->trace(__FUNCTION__);}
#else
#define REP_FUNC
#endif

/*************************/
/******** Time ***********/
/*************************/

#define SW_RESET Logger::sw->reset()
#define SW_COUNT Logger::sw->elapsed().count()

/****************************/
/******** Assert ************/
/****************************/

class AssertFailExcept
{
public:
};

#if defined(CAGE_DEBUG)
#define ASSERT(condition, msg, ...) \
if(!(condition)) {\
  Logger::dev_logger->critical("in file {} line {}: " msg , __FILE__, __LINE__, ##__VA_ARGS__);\
  throw AssertFailExcept();}
#else
#define ASSERT(condition, msg, ...)
#endif

/*****************************/
/******** Breakout ***********/
/*****************************/

class BreakoutExcept
{
public:
};

#define BREAKOUT Logger::dev_logger->debug("break out."); throw BreakoutExcept()
#define UNUSED(var) (void)var


/*****************************/
/********  Memory  ***********/
/*****************************/

/** Copy from VolumeMesher::makePolyhedralMesh.cpp **/
#ifdef _MSC_VER
#include <windows.h>
#include <psapi.h>

// To ensure correct resolution of symbols, add Psapi.lib to TARGETLIBS
// and compile with -DPSAPI_VERSION=1

inline double getPeakMegabytesUsed()
{
  HANDLE hProcess = OpenProcess(PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, FALSE, GetCurrentProcessId());
  if (NULL == hProcess) return 0;

  PROCESS_MEMORY_COUNTERS pmc;
  double mem = 0;
  if (GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc)))
  {
    mem = pmc.PeakWorkingSetSize / 1048576.0;
  }

  CloseHandle(hProcess);
  return mem;
}

inline double getMegabytesUsed()
{
  HANDLE hProcess = OpenProcess(PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, FALSE, GetCurrentProcessId());
  if (NULL == hProcess) return 0;

  PROCESS_MEMORY_COUNTERS pmc;
  double mem = 0;
  if (GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc)))
  {
    mem = pmc.WorkingSetSize / 1048576.0;
  }

  CloseHandle(hProcess);
  return mem;
}
#else
inline double getPeakMegabytesUsed() { return 0.0; }
inline double getMegabytesUsed() { return 0.0; }
#endif
}// namespace SimpleUtils
}// namespace Cage