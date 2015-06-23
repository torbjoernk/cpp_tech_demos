#include <cassert>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>
using namespace std;

typedef chrono::system_clock Clock;
typedef chrono::nanoseconds ClockResolution;

#include <boost/format.hpp>
boost::format log_fmt;

#include <mpi.h>

#define WITH_MPI
#include "../../logging.hpp"

#ifndef _MPI__SIMPLE_COMM_CONFIG_HPP_
#define _MPI__SIMPLE_COMM_CONFIG_HPP_

#define MAX_ITER                 5
#define BASE_DELAY            1000  // nanoseconds
#define FINE_MULTIPLIER     200000
#define COARSE_MULTIPLIER    10000
#define STATE_MULTIPLIER        10

#define TOTAL_STEPS              4
#define RESIDUAL_TOL             2  // seconds


inline static MPI_Status MPI_Status_factory()
{
  MPI_Status stat;
  stat.MPI_ERROR = MPI_SUCCESS;
  stat.MPI_SOURCE = MPI_ANY_SOURCE;
  stat.MPI_TAG = MPI_ANY_TAG;
  return stat;
}


enum class PState : int {
  // overall state
  CONVERGED        =  0,
  FAILED           =  1,

  // iterating states
  PREDICTING       = 10,
  ITERATING        = 11,

  // coarse level
  PRE_ITER_COARSE  = 20,
  ITER_COARSE      = 21,
  POST_ITER_COARSE = 22,

  // fine level
  PRE_ITER_FINE    = 30,
  ITER_FINE        = 31,
  POST_ITER_FINE   = 32,

  // last
  UNKNOWN          = 99,
};

struct ProcessState
{
  PState state= PState::UNKNOWN;
  int iter = -1;
  double residual = numeric_limits<double>::max();
};

int block_length[3] = {1, 1, 1};
MPI_Aint block_displace[3] = {
  0,
  sizeof(int),
  sizeof(int) + sizeof(int)
};
MPI_Datatype block_types[3] = {
  MPI_INT, MPI_INT, MPI_DOUBLE
};

MPI_Datatype process_state_type;
int process_state_type_size;


#endif // _MPI__SIMPLE_COMM_CONFIG_HPP_
