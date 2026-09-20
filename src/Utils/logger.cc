#include "logger.hh"
#include <atomic>

namespace Cage
{
namespace SimpleUtils
{
std::shared_ptr<spdlog::logger> Logger::dev_logger = nullptr;
std::shared_ptr<spdlog::logger> Logger::user_logger = nullptr;
std::shared_ptr<spdlog::stopwatch> Logger::sw = nullptr;

namespace
{
std::atomic<long long> logging_nanoseconds{ 0 };

// Forwards to a sink and adds the time spent there to logging_nanoseconds,
// so phase timers can leave log output out of the computation time.
class TimedSink final : public spdlog::sinks::sink
{
public:
  explicit TimedSink(std::shared_ptr<spdlog::sinks::sink> sink) :
    inner(std::move(sink))
  {
    set_level(inner->level());
  }

  void log(const spdlog::details::log_msg& msg) override
  {
    const auto start = std::chrono::steady_clock::now();
    inner->log(msg);
    add_since(start);
  }
  void flush() override
  {
    const auto start = std::chrono::steady_clock::now();
    inner->flush();
    add_since(start);
  }
  void set_pattern(const std::string& pattern) override
  {
    inner->set_pattern(pattern);
  }
  void set_formatter(std::unique_ptr<spdlog::formatter> formatter) override
  {
    inner->set_formatter(std::move(formatter));
  }

private:
  static void add_since(std::chrono::steady_clock::time_point start)
  {
    logging_nanoseconds += std::chrono::duration_cast<std::chrono::nanoseconds>(
      std::chrono::steady_clock::now() - start).count();
  }

  std::shared_ptr<spdlog::sinks::sink> inner;
};

std::shared_ptr<spdlog::sinks::sink> timed(std::shared_ptr<spdlog::sinks::sink> sink)
{
  return std::make_shared<TimedSink>(std::move(sink));
}
}

void Logger::InitLogger(
  spdlog::level::level_enum console_level,
  bool file_log,
  spdlog::level::level_enum file_level,
  const std::string& file_path)
{
  auto dev_console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_st>();
  dev_console_sink->set_level(console_level);
  auto user_console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_st>();
  user_console_sink->set_level(console_level);

  if (file_log)
  {
    auto dev_file_sink = std::make_shared<spdlog::sinks::basic_file_sink_st>(file_path, true);
    dev_file_sink->set_level(file_level);
    auto user_file_sink = std::make_shared<spdlog::sinks::basic_file_sink_st>(file_path, true);
    user_file_sink->set_level(file_level);
    dev_logger = std::make_shared<spdlog::logger>(
      spdlog::logger(
        "console and file logger",
        { timed(dev_console_sink), timed(dev_file_sink) }
    ));
    user_logger = std::make_shared<spdlog::logger>(
      spdlog::logger(
        "console and file logger",
        { timed(user_console_sink), timed(user_file_sink) }
    ));
  }
  else
  {
    dev_logger = std::make_shared<spdlog::logger>(
      "console logger", timed(dev_console_sink));
    user_logger = std::make_shared<spdlog::logger>(
      "console logger", timed(user_console_sink));
  }
  dev_logger->set_pattern("%^[%H:%M:%S][dev][%l]%v%$");
  user_logger->set_pattern("%^[%H:%M:%S][user][%l]%v%$");
  dev_logger->set_level(spdlog::level::trace);
  user_logger->set_level(spdlog::level::trace);

  sw = std::make_shared<spdlog::stopwatch>();
  dev_logger->info("create logger.");
  user_logger->info("create logger.");
}

void Logger::updateFileLog(bool file_log, spdlog::level::level_enum file_level, const std::string& file_path)
{
  if (dev_logger->sinks().size() == 2)
    dev_logger->sinks().pop_back();
  if (user_logger->sinks().size() == 2)
    user_logger->sinks().pop_back();

  if (file_log)
  {
    auto dev_file_sink = std::make_shared<spdlog::sinks::basic_file_sink_st>(file_path, true);
    dev_file_sink->set_level(file_level);
    dev_logger->sinks().push_back(timed(dev_file_sink));

    auto user_file_sink = std::make_shared<spdlog::sinks::basic_file_sink_st>(file_path, true);
    user_file_sink->set_level(file_level);
    user_logger->sinks().push_back(timed(user_file_sink));
  }

  dev_logger->set_pattern("%^[%H:%M:%S][dev][%l]%v%$");
  user_logger->set_pattern("%^[%H:%M:%S][user][%l]%v%$");
}

void Logger::StopProgram()
{
  dev_logger->critical("trigger stopping program.");
  exit(1);
}

double Logger::loggingSeconds()
{
  return static_cast<double>(logging_nanoseconds.load()) * 1e-9;
}

PhaseTimer* PhaseTimer::active = nullptr;
int PhaseTimer::exclusionDepth = 0;

void PhaseTimer::start()
{
  if (running)
    return;
  running = true;
  everStarted = true;
  active = this;
  startTime = Clock::now();
  loggingAtStart = Logger::loggingSeconds();
}

void PhaseTimer::stop()
{
  if (!running)
    return;
  wall += std::chrono::duration<double>(Clock::now() - startTime).count();
  logging += Logger::loggingSeconds() - loggingAtStart;
  running = false;
  if (active == this)
    active = nullptr;
}

double PhaseTimer::seconds() const
{
  double total = wall - logging;
  if (running)
  {
    total += std::chrono::duration<double>(Clock::now() - startTime).count() -
      (Logger::loggingSeconds() - loggingAtStart);
  }
  for (const auto& entry : excluded)
    total -= entry.second;
  return total > 0.0 ? total : 0.0;
}

PhaseTimer::Exclusion::Exclusion(const char* label_) :
  timer(exclusionDepth == 0 && active && active->running ? active : nullptr),
  label(label_),
  start(Clock::now()),
  loggingAtStart(Logger::loggingSeconds())
{
  exclusionDepth++;
}

PhaseTimer::Exclusion::~Exclusion()
{
  exclusionDepth--;
  if (!timer)
    return;
  // Logging inside the section is already subtracted as logging.
  const double seconds = std::chrono::duration<double>(Clock::now() - start).count() -
    (Logger::loggingSeconds() - loggingAtStart);
  timer->excluded[label] += seconds > 0.0 ? seconds : 0.0;
}
}// namespace SimpleUtils
}// namespace Cage
