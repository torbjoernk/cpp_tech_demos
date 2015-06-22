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
#include "../logging.hpp"


inline static MPI_Status MPI_Status_factory()
{
  MPI_Status stat;
  stat.MPI_ERROR = MPI_SUCCESS;
  stat.MPI_SOURCE = MPI_ANY_SOURCE;
  stat.MPI_TAG = MPI_ANY_TAG;
  return stat;
}


#define MAX_ITER                 5
#define BASE_DELAY            1000  // nanoseconds
#define FINE_MULTIPLIER     500000
#define COARSE_MULTIPLIER    10000
#define STATE_MULTIPLIER        10

#define TOTAL_STEPS              4
#define RESIDUAL_TOL             2  // seconds


enum class PState : int {
  // overall state
  CONVERGED        =  0,
  FAILED           =  1,
  // iterating states
  PREDICT          = 10,
  PRE_ITER_COARSE  = 20,
  ITER_COARSE      = 21,
  POST_ITER_COARSE = 22,
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


struct ProcessData
{
  int size = -1;
  int rank = -1;
  bool iam_first;
  bool iam_last;
  int prev = -1;
  int next = -1;

  double coarse_val;
  double fine_val;

  double mpi_start;

  ProcessState state;

  MPI_Status state_stat;
  MPI_Status fine_stat;
  MPI_Status coarse_stat;
  MPI_Request state_req;
  MPI_Request fine_req;
  MPI_Request coarse_req;

  void init(const double start_time) {
    assert(rank >= 0 && size > 0);
    this->mpi_start = start_time;
    this->iam_first = (rank == 0);
    this->iam_last = (rank == size - 1);
    if (!this->iam_first) { this->prev = this->rank - 1; }
    if (!this->iam_last) { this->next = this->rank + 1; }
    this->state_stat = MPI_Status_factory();
    this->fine_stat = MPI_Status_factory();
    this->coarse_stat = MPI_Status_factory();
    this->state_req = MPI_REQUEST_NULL;
    this->fine_req = MPI_REQUEST_NULL;
    this->coarse_req = MPI_REQUEST_NULL;
  }
};

static int fine_tag(const int iter)   { return (iter + 1) * FINE_MULTIPLIER; }
static int coarse_tag(const int iter) { return (iter + 1) * COARSE_MULTIPLIER; }
static int state_tag(const int iter)  { return (iter + 1) * STATE_MULTIPLIER; }


void doing_fine(ProcessData &data, const int iter) {
  int mpi_err = MPI_SUCCESS;

  if (!data.iam_first) {
    mpi_err = MPI_Recv(&(data.fine_val), 1, MPI_DOUBLE, data.prev, fine_tag(iter),
                       MPI_COMM_WORLD, &(data.fine_stat));
    assert(mpi_err == MPI_SUCCESS);
  }

  data.fine_val += (data.state.iter + 1) * FINE_MULTIPLIER + iter * 0.001;

  chrono::time_point<Clock> start, end;
  ClockResolution duration;
  start = Clock::now();
  do {
    end = Clock::now();
    duration = end - start;
  } while(duration.count() < BASE_DELAY * FINE_MULTIPLIER);

  if (!data.iam_last) {
    if (data.fine_req != MPI_REQUEST_NULL) {
      mpi_err = MPI_Wait(&(data.fine_req), &(data.fine_stat));
      assert(mpi_err == MPI_SUCCESS);
    }
    mpi_err = MPI_Isend(&(data.fine_val), 1, MPI_DOUBLE, data.next, fine_tag(iter),
                        MPI_COMM_WORLD, &(data.fine_req));
    assert(mpi_err == MPI_SUCCESS);
  }
}

void doing_coarse(ProcessData &data, const int iter) {
  int mpi_err = MPI_SUCCESS;

  if (!data.iam_first) {
    mpi_err = MPI_Recv(&(data.coarse_val), 1, MPI_DOUBLE, data.prev, coarse_tag(iter),
                       MPI_COMM_WORLD, &(data.coarse_stat));
    assert(mpi_err == MPI_SUCCESS);
  }

  data.coarse_val += (data.state.iter + 1) * COARSE_MULTIPLIER + iter * 0.001;

  chrono::time_point<Clock> start, end;
  ClockResolution duration;
  start = Clock::now();
  do {
    end = Clock::now();
    duration = end - start;
  } while(duration.count() < BASE_DELAY * COARSE_MULTIPLIER);

  if (!data.iam_last) {
    if (data.coarse_req != MPI_REQUEST_NULL) {
      mpi_err = MPI_Wait(&(data.coarse_req), &(data.coarse_stat));
      assert(mpi_err == MPI_SUCCESS);
    }
    mpi_err = MPI_Isend(&(data.coarse_val), 1, MPI_DOUBLE, data.next, coarse_tag(iter),
                        MPI_COMM_WORLD, &(data.coarse_req));
    assert(mpi_err == MPI_SUCCESS);
  }
}


void check_finished(ProcessData &data, const int iter) {
  int mpi_err = MPI_SUCCESS;
  ProcessState other_state;

  if (!data.iam_first) {
    mpi_err = MPI_Recv(&other_state, 1, process_state_type, data.prev, state_tag(iter),
                       MPI_COMM_WORLD, &(data.state_stat));
    assert(mpi_err == MPI_SUCCESS);

    assert(other_state.iter == iter);
  } else {
    other_state.state = PState::CONVERGED;
  }

  double curr_time = MPI_Wtime();
  data.state.residual = curr_time - data.mpi_start;

  if (other_state.state == PState::FAILED) {
      data.state.state = PState::FAILED;
  } else if (other_state.state == PState::CONVERGED && data.state.residual > RESIDUAL_TOL) {
    data.state.state = PState::CONVERGED;
  }

  if (!data.iam_last) {
    if (data.state_req != MPI_REQUEST_NULL) {
      mpi_err = MPI_Cancel(&(data.state_req));
      assert(mpi_err == MPI_SUCCESS);

      mpi_err = MPI_Wait(&(data.state_req), &(data.state_stat));
      assert(mpi_err == MPI_SUCCESS);
    }
    mpi_err = MPI_Isend(&(data.state), 1, process_state_type, data.next, state_tag(iter),
                        MPI_COMM_WORLD, &(data.state_req));
    assert(mpi_err == MPI_SUCCESS);
  }
}


int main(int argn, char** argv) {
//                        iter   residual  coarse     fine
  log_fmt = boost::format("%4.d    %8.4f    %5.1f     %10.1f");

  MPI_Init(&argn, &argv);

  init_log(argn, argv);

  int size = -1;
  int rank = -1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Type_create_struct(3, block_length, block_displace, block_types, &process_state_type);
  MPI_Type_commit(&process_state_type);

  int curr_step_start = 0;

  do {
    double mpi_start = MPI_Wtime();

    int mpi_err = MPI_SUCCESS;

    ProcessData myself;

    myself.rank = rank;
    if (curr_step_start + size - 1 < TOTAL_STEPS) {
      // all ranks fit
      myself.size = size;
    } else {
      if (curr_step_start + rank < TOTAL_STEPS) {
        myself.size = TOTAL_STEPS - curr_step_start;
      } else {
        // this rank hasn't to do anything anymore
        break;
      }
    }
    myself.init(mpi_start);

    for(int iter = 0; myself.state.state > PState::FAILED; ++iter) {
      myself.state.iter = iter;
      doing_coarse(myself, iter);
      doing_fine(myself, iter);
      check_finished(myself, iter);
      log_fmt % iter % myself.state.residual % myself.coarse_val % myself.fine_val;
      LOG(INFO) << log_fmt;
    }

    if (myself.state_req != MPI_REQUEST_NULL) {
      mpi_err = MPI_Wait(&(myself.state_req), &(myself.coarse_stat));
      assert(mpi_err == MPI_SUCCESS);
    }
    if (myself.coarse_req != MPI_REQUEST_NULL) {
      mpi_err = MPI_Wait(&(myself.coarse_req), &(myself.coarse_stat));
      assert(mpi_err == MPI_SUCCESS);
    }
    if (myself.fine_req != MPI_REQUEST_NULL) {
      mpi_err = MPI_Wait(&(myself.fine_req), &(myself.fine_stat));
      assert(mpi_err == MPI_SUCCESS);
    }

    curr_step_start += size;
  } while (curr_step_start < TOTAL_STEPS);

  MPI_Type_free(&process_state_type);
  MPI_Finalize();
}
