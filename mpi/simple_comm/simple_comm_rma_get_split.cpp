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

#include "simple_comm_config.hpp"


enum class PState : int {
  // overall state
  CONVERGED        =  0,
  FAILED           =  1,
  // iterating states
  ITERATING        = 10,
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

  double coarse_val_in, coarse_val, coarse_val_out;
  double fine_val_in, fine_val, fine_val_out;

  double mpi_start;

  ProcessState state;

  MPI_Win coarse_win;
  MPI_Win fine_win;
  MPI_Win state_win;

  void init(const double start_time) {
    assert(rank >= 0 && size > 0);
    this->mpi_start = start_time;
    this->iam_first = (rank == 0);
    this->iam_last = (rank == size - 1);
    if (!this->iam_first) { this->prev = this->rank - 1; }
    if (!this->iam_last) { this->next = this->rank + 1; }
  }

  void create_windows()
  {
    if (this->size > 1) {
      VLOG(6) << "creating windows";
      MPI_Win_create(&fine_val_out, sizeof(double), MPI_DOUBLE, MPI_INFO_NULL, MPI_COMM_WORLD, &fine_win);
      MPI_Win_create(&coarse_val_out, sizeof(double), MPI_DOUBLE, MPI_INFO_NULL, MPI_COMM_WORLD, &coarse_win);
      MPI_Win_create(&state, 1, process_state_type_size, MPI_INFO_NULL, MPI_COMM_WORLD, &state_win);
    } else {
      VLOG(6) << "no need for communication";
    }
  }

  void free_windows()
  {
    if (this->size > 1) {
      VLOG(6) << "freeing windows";
      MPI_Win_free(&fine_win);
      MPI_Win_free(&coarse_win);
      MPI_Win_free(&state_win);
    } else {
      VLOG(6) << "no need for communication";
    }
  }
};

static int fine_tag(const int iter)   { return (iter + 1) * FINE_MULTIPLIER; }
static int coarse_tag(const int iter) { return (iter + 1) * COARSE_MULTIPLIER; }
static int state_tag(const int iter)  { return (iter + 1) * STATE_MULTIPLIER; }


void update_state(ProcessData &data, const PState &new_state, const int iter=-1) {
  int mpi_err = MPI_SUCCESS;

  if (data.size > 1) {
    VLOG(5) << "locking local state window";
    mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, data.rank, 0, data.state_win); assert(mpi_err == MPI_SUCCESS);
  }

  data.state.state = new_state;
  if (iter != -1) {
    data.state.iter = iter;
  }

  if (data.size > 1) {
    mpi_err = MPI_Win_unlock(data.rank, data.state_win); assert(mpi_err == MPI_SUCCESS);
    VLOG(5) << "unlocked local state window";
  }
}


void doing_pre_fine(ProcessData &data) {
  update_state(data, PState::PRE_ITER_FINE);

  int mpi_err = MPI_SUCCESS;

  if (data.size > 1 && !data.iam_first) {
    VLOG(5) << "locking fine window to rank " << data.prev;
    mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, data.prev, 0, data.fine_win); assert(mpi_err == MPI_SUCCESS);

    VLOG(4) << "getting fine value from " << data.prev;
    mpi_err = MPI_Get(&(data.fine_val_in), 1, MPI_DOUBLE, data.prev, 0, 1, MPI_DOUBLE, data.fine_win); assert(mpi_err == MPI_SUCCESS);

    mpi_err = MPI_Win_unlock(data.prev, data.fine_win); assert(mpi_err == MPI_SUCCESS);
    VLOG(5) << "unlocked fine window to rank " << data.prev;
  }
}

void doing_fine(ProcessData &data, const int iter) {
  update_state(data, PState::ITER_FINE);

  VLOG(2) << "start computation";
  data.fine_val = data.fine_val_in + (data.state.iter + 1) * FINE_MULTIPLIER + iter * 0.001;
  VLOG(3) << data.fine_val << " = " << data.fine_val_in << " + " << ((data.rank + 1) * FINE_MULTIPLIER) << " + " << (iter * 0.001);

  chrono::time_point<Clock> start, end;
  ClockResolution duration;
  start = Clock::now();
  do {
    end = Clock::now();
    duration = end - start;
  } while(duration.count() < BASE_DELAY * FINE_MULTIPLIER);
  VLOG(2) << "done computation";
}

void doing_post_fine(ProcessData &data) {
  update_state(data, PState::POST_ITER_FINE);

  int mpi_err = MPI_SUCCESS;

  if (data.size > 1) {
    VLOG(5) << "locking local fine window";
    mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, data.rank, 0, data.fine_win); assert(mpi_err == MPI_SUCCESS);
  }

  data.fine_val_out = data.fine_val;

  if (data.size > 1) {
    mpi_err = MPI_Win_unlock(data.rank, data.fine_win); assert(mpi_err == MPI_SUCCESS);
    VLOG(5) << "unlocked local fine window";
  }
}


void doing_pre_coarse(ProcessData &data) {
  update_state(data, PState::PRE_ITER_COARSE);

  int mpi_err = MPI_SUCCESS;

  if (data.size > 1 && !data.iam_first) {
    VLOG(5) << "locking coarse window to rank " << data.prev;
    mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, data.prev, 0, data.coarse_win); assert(mpi_err == MPI_SUCCESS);

    VLOG(4) << "getting coarse value from " << data.prev;
    mpi_err = MPI_Get(&(data.coarse_val_in), 1, MPI_DOUBLE, data.prev, 0, 1, MPI_DOUBLE, data.coarse_win); assert(mpi_err == MPI_SUCCESS);

    mpi_err = MPI_Win_unlock(data.prev, data.coarse_win); assert(mpi_err == MPI_SUCCESS);
    VLOG(5) << "unlocked coarse window to rank " << data.prev;
  }
}

void doing_coarse(ProcessData &data, const int iter) {
  update_state(data, PState::ITER_COARSE);

  VLOG(2) << "start computation";
  data.coarse_val = data.coarse_val_in + (data.state.iter + 1) * COARSE_MULTIPLIER + iter * 0.001;
  VLOG(3) << data.coarse_val << " = " << data.coarse_val_in << " + " << ((data.rank + 1) * COARSE_MULTIPLIER) << " + " << (iter * 0.001);

  chrono::time_point<Clock> start, end;
  ClockResolution duration;
  start = Clock::now();
  do {
    end = Clock::now();
    duration = end - start;
  } while(duration.count() < BASE_DELAY * COARSE_MULTIPLIER);
  VLOG(2) << "done computation";
}

void doing_post_coarse(ProcessData &data) {
  update_state(data, PState::POST_ITER_COARSE);

  int mpi_err = MPI_SUCCESS;

  if (data.size > 1) {
    VLOG(5) << "locking local coarse window";
    mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, data.rank, 0, data.coarse_win); assert(mpi_err == MPI_SUCCESS);
  }

  data.coarse_val_out = data.coarse_val;

  if (data.size > 1) {
    mpi_err = MPI_Win_unlock(data.rank, data.coarse_win); assert(mpi_err == MPI_SUCCESS);
    VLOG(5) << "unlocked local coarse window";
  }
}


void check_finished(ProcessData &data, const int iter) {
  int mpi_err = MPI_SUCCESS;
  ProcessState other_state;

  if (data.size > 1 && !data.iam_first) {
    VLOG(5) << "locking state window to rank " << data.prev;
    mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, data.prev, 0, data.state_win); assert(mpi_err == MPI_SUCCESS);

    VLOG(4) << "getting state from " << data.prev;
    mpi_err = MPI_Get(&other_state, 1, process_state_type, data.prev, 0, 1, process_state_type, data.state_win); assert(mpi_err == MPI_SUCCESS);

    mpi_err = MPI_Win_unlock(data.prev, data.state_win); assert(mpi_err == MPI_SUCCESS);
    VLOG(5) << "unlocked state window to rank " << data.prev;
  } else {
    other_state.state = PState::CONVERGED;
  }

  if (data.size > 1) {
    VLOG(5) << "locking local state window";
    mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, data.rank, 0, data.state_win); assert(mpi_err == MPI_SUCCESS);
  }

  double curr_time = MPI_Wtime();
  data.state.residual = curr_time - data.mpi_start;

  if (other_state.state == PState::FAILED) {
      data.state.state = PState::FAILED;
  } else if (other_state.state == PState::CONVERGED && data.state.residual > RESIDUAL_TOL) {
    data.state.state = PState::CONVERGED;
  }

  if (data.size > 1) {
    mpi_err = MPI_Win_unlock(data.rank, data.state_win); assert(mpi_err == MPI_SUCCESS);
    VLOG(5) << "unlocked local state window";
  }
}


int main(int argn, char** argv) {
//                         iter      residual   coarse     fine
  log_fmt = boost::format("%4.d      %12.3f      %12.3f      %12.3f");

  MPI_Init(&argn, &argv);

  init_log(argn, argv);

  int size = -1;
  int rank = -1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  assert(size <= TOTAL_STEPS);

  MPI_Type_create_struct(3, block_length, block_displace, block_types, &process_state_type);
  MPI_Type_commit(&process_state_type);
  MPI_Type_size(process_state_type, &process_state_type_size);

  int curr_step_start = 0;
  double initial_value = 1.0;

  do {
    double mpi_start = MPI_Wtime();

    int mpi_err = MPI_SUCCESS;

    ProcessData myself;

    myself.rank = rank;
    int working_size = (curr_step_start + size - 1 < TOTAL_STEPS) ? size : TOTAL_STEPS % size;
    LOG(INFO) << working_size << " processes will work now";

    myself.size = working_size;
    if (rank < working_size) {
      LOG(INFO) << "doing step " << (curr_step_start + rank);
      myself.init(mpi_start);
      myself.create_windows();

      myself.fine_val_in = initial_value;
      myself.coarse_val_in = initial_value / FINE_MULTIPLIER;
      LOG(INFO) << "inital values:\tcoarse=" << std::fixed << std::setprecision(3) << myself.coarse_val_in << "\tfine=" << myself.fine_val_in;

      for(int iter = 0; myself.state.state > PState::FAILED; ++iter) {
        update_state(myself, PState::ITERATING, iter);

        doing_pre_coarse(myself);
        doing_coarse(myself, iter);
        doing_post_coarse(myself);

        doing_pre_fine(myself);
        doing_fine(myself, iter);
        doing_post_fine(myself);

        update_state(myself, PState::ITERATING);

        check_finished(myself, iter);
        log_fmt % iter % myself.state.residual % myself.coarse_val_out % myself.fine_val_out;
        LOG(INFO) << log_fmt;
      }

      myself.free_windows();
    } else {
      // this rank hasn't to do anything anymore
      LOG(WARNING) << "hasn't work anymore";
    }

    VLOG(4) << "broadcasting final value to all";
    mpi_err = MPI_Bcast(&(myself.fine_val_out), 1, MPI_DOUBLE, working_size - 1, MPI_COMM_WORLD); assert(mpi_err == MPI_SUCCESS);
    initial_value = myself.fine_val_out;

    curr_step_start += size;
  } while (curr_step_start < TOTAL_STEPS);

  MPI_Type_free(&process_state_type);
  VLOG(6) << "finalizing";
  MPI_Finalize();
}