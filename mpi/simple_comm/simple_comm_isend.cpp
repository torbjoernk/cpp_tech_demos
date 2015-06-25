#include "simple_comm_config.hpp"


class NonBlockingSendProcess
  : public Process
{
  public:
    NonBlockingSendProcess(const double start_time)
      : Process(start_time)
    {}
};


class NonBlockingSendCommunicator
  : public Communicator
{
  public:
    MPI_Request fine_req;

    NonBlockingSendCommunicator(const int size)
      : Communicator(size)
    {
      this->fine_req = MPI_REQUEST_NULL;
    }

    void recv_fine(shared_ptr<Process> proc) override
    {
      if (!this->iam_first && proc->state.iter > 0) {
        CVLOG(5, "Communicator") << "receiving fine data from " << this->prev << " with tag "
                                 << fine_tag(proc->state.iter - 1);
        mpi_err = MPI_Recv(&(proc->fine_val), 1, MPI_DOUBLE, this->prev, fine_tag(proc->state.iter - 1),
                           MPI_COMM_WORLD, &(this->fine_stat));
        assert(mpi_err == MPI_SUCCESS);
      }
    }

    void send_fine(shared_ptr<Process> proc) override
    {
      if (!this->iam_last) {
        CVLOG(5, "Communicator") << "sending fine data to " << this->next
                                 << " with tag " << fine_tag(proc->state.iter);
        if (this->fine_req != MPI_REQUEST_NULL) {
          mpi_err = MPI_Wait(&(this->fine_req), &(this->fine_stat));
          assert(mpi_err == MPI_SUCCESS);
        }
        mpi_err = MPI_Isend(&(proc->fine_val), 1, MPI_DOUBLE, this->next, fine_tag(proc->state.iter),
                            MPI_COMM_WORLD, &(this->fine_req));
        assert(mpi_err == MPI_SUCCESS);
      }
    }

    void recv_coarse(shared_ptr<Process> proc) override
    {
      if (!this->iam_first && proc->state.iter != - this->rank) {
        CVLOG(5, "Communicator") << "receiving coarse data from " << this->prev
                                 << " with tag " << coarse_tag(proc->state.iter);
        mpi_err = MPI_Recv(&(proc->coarse_val), 1, MPI_DOUBLE, this->prev, coarse_tag(proc->state.iter),
                           MPI_COMM_WORLD, &(this->coarse_stat));
        assert(mpi_err == MPI_SUCCESS);
      } else {
        CVLOG(5, "Communicator") << "not receiving as I'm first or in my first predict";
      }
    }

    void send_coarse(shared_ptr<Process> proc) override
    {
      if (!this->iam_last) {
        CVLOG(5, "Communicator") << "sending coarse data to " << this->next
                                 << " with tag " << coarse_tag(proc->state.iter);
        mpi_err = MPI_Send(&(proc->coarse_val), 1, MPI_DOUBLE, this->next, coarse_tag(proc->state.iter),
                           MPI_COMM_WORLD);
        assert(mpi_err == MPI_SUCCESS);
      }
    }

    void recv_state(shared_ptr<Process> proc) override
    {
      if (!this->iam_first) {
        CVLOG(5, "Communicator") << "receiving state info from " << this->prev;
        mpi_err = MPI_Recv(&(proc->prev_state), 1, process_state_type, this->prev, state_tag(proc->state.iter),
                           MPI_COMM_WORLD, &(this->state_stat));
        assert(mpi_err == MPI_SUCCESS);
      } else {
        proc->prev_state.state = PState::CONVERGED;
      }
    }

    void send_state(shared_ptr<Process> proc) override
    {
      if (!this->iam_last) {
        CVLOG(5, "Communicator") << "sending state info to " << this->next;
        mpi_err = MPI_Send(&(proc->state), 1, process_state_type, this->next, state_tag(proc->state.iter),
                           MPI_COMM_WORLD);
        assert(mpi_err == MPI_SUCCESS);
      }
    }

    void cleanup() override
    {
      if (this->fine_req != MPI_REQUEST_NULL) {
        mpi_err = MPI_Wait(&(this->fine_req), &(this->fine_stat));
        assert(mpi_err == MPI_SUCCESS);
      }
    }
};


int main(int argn, char** argv) {
  //                       iter    residual    coarse      fine
  log_fmt = boost::format("%4.d    %12.6f      %12.3f      %12.3f");

  MPI_Init(&argn, &argv);

  init_log(argn, argv);
  init_additional_loggers();

  int size = -1;
  int rank = -1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Type_create_struct(3, block_length, block_displace, block_types, &process_state_type);
  MPI_Type_commit(&process_state_type);

  int curr_step_start = 0;
  double initial_value = FINE_MULTIPLIER;

  do {
    double mpi_start = MPI_Wtime();

    int working_size = (curr_step_start + size - 1 < TOTAL_STEPS) ? size : TOTAL_STEPS % size;
    LOG(INFO) << working_size << " processes will work now";

    shared_ptr<NonBlockingSendCommunicator> comm = make_shared<NonBlockingSendCommunicator>(working_size);
    shared_ptr<NonBlockingSendProcess> proc = make_shared<NonBlockingSendProcess>(mpi_start);
    FixedIterationController controll(comm, proc);

//     controll._fine_delay = bind(random_delay, placeholders::_1, BASE_DELAY * FINE_MULTIPLIER, BASE_DELAY * FINE_MULTIPLIER * FINE_DELAY_VARIANCE);
//     controll._coarse_delay = bind(random_delay, placeholders::_1, BASE_DELAY * COARSE_MULTIPLIER, BASE_DELAY * COARSE_MULTIPLIER * COARSE_DELAY_VARIANCE);

//     controll._fine_delay = bind(deminishing_delay, placeholders::_1, BASE_DELAY * FINE_MULTIPLIER, FINE_DEMINISH);
//     controll._coarse_delay = bind(deminishing_delay, placeholders::_1, BASE_DELAY * COARSE_MULTIPLIER, COARSE_DEMINISH);

    controll._fine_delay = bind(specific_delay, placeholders::_1, BASE_DELAY * FINE_MULTIPLIER, comm->rank, true);
    controll._coarse_delay = bind(specific_delay, placeholders::_1, BASE_DELAY * COARSE_MULTIPLIER, comm->rank, false);

    if (rank < working_size) {
      LOG(INFO) << "doing step " << (curr_step_start + rank);
      proc->fine_val = initial_value;
      proc->coarse_val = initial_value / FINE_MULTIPLIER;
      LOG(INFO) << "inital values:\tcoarse=" << std::fixed << std::setprecision(3)
                << proc->coarse_val << "\tfine=" << proc->fine_val;

      controll.run();
    } else {
      // this rank hasn't to do anything anymore
      LOG(WARNING) << "hasn't work anymore";
    }

    initial_value = proc->fine_val;

    curr_step_start += size;
  } while (curr_step_start < TOTAL_STEPS);

  MPI_Type_free(&process_state_type);
  VLOG(6) << "finalizing";
  MPI_Finalize();
}
