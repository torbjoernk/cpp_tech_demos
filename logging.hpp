#include <iomanip>
#include <string>
using namespace std;

#include <boost/algorithm/string.hpp>

#ifdef WITH_MPI
  #include <mpi.h>
#endif

#define MILLISEC_WIDTH "6"

#ifdef NO_LOGGING
  #define ELPP_DISABLE_LOGGING
#else
  // enable easy logging of STL containers
  #define ELPP_STL_LOGGING
  // disable creation of default log file
  #define ELPP_NO_DEFAULT_LOG_FILE
  // enable passing `--logging-flags` via command line
  #define ELPP_LOGGING_FLAGS_FROM_ARG
  #include "easylogging++.h"
  INITIALIZE_EASYLOGGINGPP
#endif


int get_rank()
{
#ifdef WITH_MPI
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
#else
  return 0;
#endif
}

inline string format_mpi_rank(const char fill = ' ')
{
  ostringstream frmter;
  frmter << std::setw(2) << std::setfill(fill) << get_rank();
  return frmter.str();
}


inline string get_log_file_name()
{
  string log_name = "log";
#ifdef WITH_MPI
  if (log_name.size() > 0) {
    log_name += "_";
  }
  log_name += "mpi-rank-" + format_mpi_rank('0');
#endif
  log_name += ".log";
  return log_name;
}


#ifndef NO_LOGGING
inline void set_global_logging_options(el::Configurations* conf,
                                       const el::Configurations* default_conf = nullptr)
{
  string to_stdout;
  if (default_conf) {
    el::Configurations* default_conf_nc = const_cast<el::Configurations*>(default_conf);
    to_stdout = default_conf_nc->get(el::Level::Info,
                                     el::ConfigurationType::ToStandardOutput)->value();
  } else {
    to_stdout = "true";
  }

  conf->setGlobally(el::ConfigurationType::MillisecondsWidth, MILLISEC_WIDTH);
  conf->setGlobally(el::ConfigurationType::ToStandardOutput, to_stdout);
  conf->setGlobally(el::ConfigurationType::Filename, get_log_file_name());
}
#endif

inline static void add_custom_logger(const string& id)
{
#ifndef NO_LOGGING
  const string TIMESTAMP = "%datetime{%H:%m:%s,%g} ";
  const string LEVEL = "%level";
  const string VLEVEL = "VERB%vlevel";
  const string POSITION = "%fbase:%line";
  const string MESSAGE = "%msg";
#ifdef WITH_MPI
  const string MPI_RANK = ", MPI " + format_mpi_rank();
#else
  const string MPI_RANK = "";
#endif

  const size_t id_length = id.size();
  string id2print = id.substr(0, 6);
  boost::to_upper(id2print);
  if (id_length < 6) {
    id2print.append(6 - id_length, ' ');
  }

  el::Logger* logger = el::Loggers::getLogger(id);
  el::Configurations* conf = logger->configurations();
  const el::Configurations* default_conf = el::Loggers::defaultConfigurations();
  set_global_logging_options(conf, default_conf);

  conf->set(el::Level::Info, el::ConfigurationType::Format,
            TIMESTAMP + "[" + id2print + ", " + LEVEL  + MPI_RANK + "] " + MESSAGE);
  conf->set(el::Level::Debug, el::ConfigurationType::Format,
            TIMESTAMP + "[" + id2print + ", " + LEVEL  + MPI_RANK + "] " + POSITION + " " + MESSAGE);
  conf->set(el::Level::Warning, el::ConfigurationType::Format,
            TIMESTAMP + "[" + id2print + ", " + LEVEL  + MPI_RANK + "] " + MESSAGE);
  conf->set(el::Level::Error, el::ConfigurationType::Format,
            TIMESTAMP + "[" + id2print + ", " + LEVEL  + MPI_RANK + "] " + MESSAGE);
  conf->set(el::Level::Fatal, el::ConfigurationType::Format,
            TIMESTAMP + "[" + id2print + ", " + LEVEL  + MPI_RANK + "] " + POSITION + " " + MESSAGE);
  conf->set(el::Level::Verbose, el::ConfigurationType::Format,
            TIMESTAMP + "[" + id2print + ", " + VLEVEL + MPI_RANK + "] " + MESSAGE);
  el::Loggers::reconfigureLogger(logger, *conf);
#endif
}

inline static void init_log(int argn, char** argv)
{
#ifndef NO_LOGGING
  START_EASYLOGGINGPP(argn, argv);

  el::Loggers::addFlag(el::LoggingFlag::LogDetailedCrashReason);
  el::Loggers::addFlag(el::LoggingFlag::DisableApplicationAbortOnFatalLog);
  el::Loggers::addFlag(el::LoggingFlag::MultiLoggerSupport);
  el::Loggers::addFlag(el::LoggingFlag::CreateLoggerAutomatically);
  el::Loggers::addFlag(el::LoggingFlag::ImmediateFlush);
  el::Loggers::removeFlag(el::LoggingFlag::ColoredTerminalOutput);

  el::Configurations defaultConf;
  defaultConf.setToDefault();

  set_global_logging_options(&defaultConf);

  el::Loggers::setDefaultConfigurations(defaultConf, true);

  add_custom_logger("default");
#endif
}
