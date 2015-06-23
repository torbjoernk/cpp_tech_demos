#include "simple_comm_config.hpp"


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
  MPI_Request fine_req;

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
    this->fine_req = MPI_REQUEST_NULL;
  }
};

static int fine_tag(const int iter)   { return (iter + 1) * FINE_MULTIPLIER; }
static int coarse_tag(const int iter) { return (abs(iter) + 1) * (COARSE_MULTIPLIER / 10); }
static int state_tag(const int iter)  { return (iter + 1) * STATE_MULTIPLIER; }


void doing_fine(ProcessData &data, const int iter) {
  int mpi_err = MPI_SUCCESS;

  if (!data.iam_first && iter > 0) {
    VLOG(5) << "receiving fine data from " << data.prev << " with tag " << fine_tag(iter - 1);
    mpi_err = MPI_Recv(&(data.fine_val), 1, MPI_DOUBLE, data.prev, fine_tag(iter - 1),
                       MPI_COMM_WORLD, &(data.fine_stat));
    assert(mpi_err == MPI_SUCCESS);
  }

  VLOG(2) << "start computation";
  data.fine_val += (data.state.iter + 1) * FINE_MULTIPLIER + iter * 0.001;
  VLOG(3) << data.fine_val << " = " << data.fine_val << " + " << ((data.state.iter + 1) * FINE_MULTIPLIER) << " + " << (iter * 0.001);

  chrono::time_point<Clock> start, end;
  ClockResolution duration;
  start = Clock::now();
  do {
    end = Clock::now();
    duration = end - start;
  } while(duration.count() < BASE_DELAY * FINE_MULTIPLIER);
  VLOG(2) << "done computation";

  if (!data.iam_last) {
    VLOG(5) << "sending fine data to " << data.next << " with tag " << fine_tag(iter);
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

  if (!data.iam_first && iter != - data.rank) {
    VLOG(5) << "receiving coarse data from " << data.prev << " with tag " << coarse_tag(iter);
    mpi_err = MPI_Recv(&(data.coarse_val), 1, MPI_DOUBLE, data.prev, coarse_tag(iter),
                       MPI_COMM_WORLD, &(data.coarse_stat));
    assert(mpi_err == MPI_SUCCESS);
  } else {
    VLOG(5) << "not receiving as I'm first or in my first predict";
  }

  VLOG(2) << "start computation";
  data.coarse_val += (data.state.iter + 1) * COARSE_MULTIPLIER + iter * 0.001;
  VLOG(3) << data.coarse_val << " = " << data.coarse_val << " + " << ((data.state.iter + 1) * COARSE_MULTIPLIER) << " + " << (iter * 0.001);

  chrono::time_point<Clock> start, end;
  ClockResolution duration;
  start = Clock::now();
  do {
    end = Clock::now();
    duration = end - start;
  } while(duration.count() < BASE_DELAY * COARSE_MULTIPLIER);
  VLOG(2) << "done computation";

  if (!data.iam_last) {
    VLOG(5) << "sending coarse data to " << data.next << " with tag " << coarse_tag(iter);
    mpi_err = MPI_Send(&(data.coarse_val), 1, MPI_DOUBLE, data.next, coarse_tag(iter),
                       MPI_COMM_WORLD);
    assert(mpi_err == MPI_SUCCESS);
  }
}


void check_finished(ProcessData &data, const int iter) {
  int mpi_err = MPI_SUCCESS;
  ProcessState other_state;

  if (!data.iam_first) {
    VLOG(5) << "receiving state info from " << data.prev;
    mpi_err = MPI_Recv(&other_state, 1, process_state_type, data.prev, state_tag(iter),
                       MPI_COMM_WORLD, &(data.state_stat));
    assert(mpi_err == MPI_SUCCESS);
  } else {
    other_state.state = PState::CONVERGED;
  }

  double curr_time = MPI_Wtime();
  data.state.residual = curr_time - data.mpi_start;

  if (other_state.state == PState::FAILED) {
      data.state.state = PState::FAILED;
  } else if (other_state.state == PState::CONVERGED) {
    VLOG(2) << "previous converged";
    if (data.state.residual > RESIDUAL_TOL) {
      VLOG(2) << "and I'm done as well";
      data.state.state = PState::CONVERGED;
    } else {
      VLOG(2) << "but I'm not yet done";
    }
  }

  if (!data.iam_last) {
    VLOG(5) << "sending state info to " << data.next;
    mpi_err = MPI_Send(&(data.state), 1, process_state_type, data.next, state_tag(iter),
                       MPI_COMM_WORLD);
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
  double initial_value = 1.0;
  int mpi_err = MPI_SUCCESS;

  do {
    double mpi_start = MPI_Wtime();

    ProcessData myself;

    myself.rank = rank;
    int working_size = (curr_step_start + size - 1 < TOTAL_STEPS) ? size : TOTAL_STEPS % size;
    LOG(INFO) << working_size << " processes will work now";

    myself.size = working_size;
    if (rank < working_size) {
      LOG(INFO) << "doing step " << (curr_step_start + rank);
      myself.init(mpi_start);
      myself.fine_val = initial_value;
      myself.coarse_val = initial_value / FINE_MULTIPLIER;
      LOG(INFO) << "inital values:\tcoarse=" << std::fixed << std::setprecision(3) << myself.coarse_val << "\tfine=" << myself.fine_val;

      for(int iter = 0; myself.state.state > PState::FAILED; ++iter) {
        myself.state.iter = iter;

        if (iter == 0) {
          for (int proc = - myself.rank; proc <= 0; ++proc) {
            VLOG(2) << "predict " << proc;
            doing_coarse(myself, proc);
          }
        } else {
          doing_coarse(myself, iter);
        }

        doing_fine(myself, iter);

        check_finished(myself, iter);
        log_fmt % iter % myself.state.residual % myself.coarse_val % myself.fine_val;
        LOG(INFO) << log_fmt;
      }
    } else {
      // this rank hasn't to do anything anymore
      LOG(WARNING) << "hasn't work anymore";
    }

    if (myself.fine_req != MPI_REQUEST_NULL) {
      mpi_err = MPI_Wait(&(myself.fine_req), &(myself.fine_stat));
      assert(mpi_err == MPI_SUCCESS);
    }

    VLOG(4) << "broadcasting final value to all";
    mpi_err = MPI_Bcast(&(myself.fine_val), 1, MPI_DOUBLE, working_size - 1, MPI_COMM_WORLD);
    assert(mpi_err == MPI_SUCCESS);
    initial_value = myself.fine_val;

    curr_step_start += size;
  } while (curr_step_start < TOTAL_STEPS);

  MPI_Type_free(&process_state_type);
  VLOG(6) << "finalizing";
  MPI_Finalize();
}
