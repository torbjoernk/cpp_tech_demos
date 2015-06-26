#include "simple_comm_config.hpp"


class RmaGetProcess
  : public Process
{
  public:
    RmaGetProcess(const double start_time)
      : Process(start_time)
    {}
};


class RmaGetCommunicator
  : public Communicator
{
  public:
    MPI_Win coarse_win;
    MPI_Win fine_win;
    MPI_Win state_win;

    RmaGetCommunicator(const int size)
      : Communicator(size)
    {
      this->coarse_win = MPI_WIN_NULL;
      this->fine_win = MPI_WIN_NULL;
      this->state_win = MPI_WIN_NULL;
    }

    void recv_fine(shared_ptr<Process> proc) override
    {
      if (this->size > 1 && !this->iam_first && proc->state.iter > 0) {
#ifndef NO_LOGGING
        CVLOG(5, "Communicator") << "locking fine window to rank " << this->prev;
#endif
        mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, this->prev, 0, this->fine_win);
        assert(mpi_err == MPI_SUCCESS);

#ifndef NO_LOGGING
        CVLOG(4, "Communicator") << "getting fine value from " << this->prev;
#endif
        mpi_err = MPI_Get(&(proc->fine_val_in), 1, MPI_DOUBLE, this->prev, 0, 1, MPI_DOUBLE, this->fine_win);
        assert(mpi_err == MPI_SUCCESS);

        mpi_err = MPI_Win_unlock(this->prev, this->fine_win);
        assert(mpi_err == MPI_SUCCESS);
#ifndef NO_LOGGING
        CVLOG(5, "Communicator") << "unlocked fine window to rank " << this->prev;
#endif
      }
    }

    void pre_comp_fine(shared_ptr<Process> proc) override
    {
      Communicator::pre_comp_fine(proc);

      if (this->size > 1) {
#ifndef NO_LOGGING
        CVLOG(5, "Communicator") << "locking local fine window";
#endif
        mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, this->rank, 0, this->fine_win);
        assert(mpi_err == MPI_SUCCESS);
      }

      proc->fine_val = proc->fine_val_in;
    }

    void post_comp_fine(shared_ptr<Process> proc) override
    {
      Communicator::post_comp_fine(proc);

      proc->fine_val_out = proc->fine_val;

      if (this->size > 1) {
        mpi_err = MPI_Win_unlock(this->rank, this->fine_win);
        assert(mpi_err == MPI_SUCCESS);
#ifndef NO_LOGGING
        CVLOG(5, "Communicator") << "unlocked local fine window";
#endif
      }
    }

    void recv_coarse(shared_ptr<Process> proc) override
    {
      if (this->size > 1 && !this->iam_first && proc->state.iter != - this->rank) {
#ifndef NO_LOGGING
        CVLOG(5, "Communicator") << "locking coarse window to rank " << this->prev;
#endif
        mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, this->prev, 0, this->coarse_win);
        assert(mpi_err == MPI_SUCCESS);

#ifndef NO_LOGGING
        CVLOG(4, "Communicator") << "getting coarse value from " << this->prev;
#endif
        mpi_err = MPI_Get(&(proc->coarse_val_in), 1, MPI_DOUBLE, this->prev, 0, 1, MPI_DOUBLE, this->coarse_win);
        assert(mpi_err == MPI_SUCCESS);

        mpi_err = MPI_Win_unlock(this->prev, this->coarse_win);
        assert(mpi_err == MPI_SUCCESS);
#ifndef NO_LOGGING
        CVLOG(5, "Communicator") << "unlocked coarse window to rank " << this->prev;
#endif
      }
    }

    void pre_comp_coarse(shared_ptr<Process> proc) override
    {
      Communicator::pre_comp_coarse(proc);

      if (this->size > 1) {
#ifndef NO_LOGGING
        CVLOG(5, "Communicator") << "locking local coarse window";
#endif
        mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, this->rank, 0, this->coarse_win);
        assert(mpi_err == MPI_SUCCESS);
      }

      proc->coarse_val = proc->coarse_val_in;
    }

    void post_comp_coarse(shared_ptr<Process> proc) override
    {
      Communicator::post_comp_coarse(proc);

      proc->coarse_val_out = proc->coarse_val;

      if (this->size > 1) {
        mpi_err = MPI_Win_unlock(this->rank, this->coarse_win);
        assert(mpi_err == MPI_SUCCESS);
#ifndef NO_LOGGING
        CVLOG(5, "Communicator") << "unlocked local coarse window";
#endif
      }
    }

    void recv_state(shared_ptr<Process> proc) override
    {
      if (this->size > 1 && !this->iam_first) {
#ifndef NO_LOGGING
        CVLOG(5, "Communicator") << "locking state window to rank " << this->prev;
#endif
        mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, this->prev, 0, this->state_win);
        assert(mpi_err == MPI_SUCCESS);

#ifndef NO_LOGGING
        CVLOG(4, "Communicator") << "getting state from " << this->prev;
#endif
        mpi_err = MPI_Get(&(proc->prev_state), 1, process_state_type, this->prev, 0, 1,
                          process_state_type, this->state_win);
        assert(mpi_err == MPI_SUCCESS);

        mpi_err = MPI_Win_unlock(this->prev, this->state_win);
        assert(mpi_err == MPI_SUCCESS);
#ifndef NO_LOGGING
        CVLOG(5, "Communicator") << "unlocked state window to rank " << this->prev;
#endif
      } else {
        proc->prev_state.state = PState::CONVERGED;
      }
    }

    void pre_update_state(shared_ptr<Process> proc) override
    {
      if (this->size > 1) {
#ifndef NO_LOGGING
        CVLOG(5, "Communicator") << "locking local state window";
#endif
        mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, this->rank, 0, this->state_win);
        assert(mpi_err == MPI_SUCCESS);
      }
    }

    void post_update_state(shared_ptr<Process> proc) override
    {
      if (this->size > 1) {
        mpi_err = MPI_Win_unlock(this->rank, this->state_win);
        assert(mpi_err == MPI_SUCCESS);
#ifndef NO_LOGGING
        VLOG(5) << "unlocked local state window";
#endif
      }
    }

    void cleanup() override
    {
      if (this->size > 1) {
#ifndef NO_LOGGING
        CVLOG(6, "Communicator") << "freeing windows";
#endif
        MPI_Win_free(&(this->fine_win));
        MPI_Win_free(&(this->coarse_win));
        MPI_Win_free(&(this->state_win));
      }
    }
};


template<typename BaseControll>
class RmaGetController
  : public BaseControll
{
  public:
    RmaGetController(shared_ptr<Communicator> comm, shared_ptr<Process> proc)
      : BaseControll(comm, proc)
    {
      if (this->comm->size > 1) {
#ifndef NO_LOGGING
        CVLOG(6, "Controller") << "creating windows";
#endif
        shared_ptr<RmaGetProcess> _proc = dynamic_pointer_cast<RmaGetProcess>(this->proc);
        shared_ptr<RmaGetCommunicator> _comm = dynamic_pointer_cast<RmaGetCommunicator>(this->comm);
        MPI_Win_create(&(_proc->fine_val_out), sizeof(double), MPI_DOUBLE, MPI_INFO_NULL,
                       MPI_COMM_WORLD, &(_comm->fine_win));
        MPI_Win_create(&(_proc->coarse_val_out), sizeof(double), MPI_DOUBLE, MPI_INFO_NULL,
                       MPI_COMM_WORLD, &(_comm->coarse_win));
        MPI_Win_create(&(_proc->state), 1, process_state_type_size, MPI_INFO_NULL, MPI_COMM_WORLD,
                       &(_comm->state_win));
      }
    }
};


int main(int argn, char** argv) {
#ifndef NO_LOGGING
  //                       iter    residual    coarse      fine
  log_fmt = boost::format("%4.d    %12.6f      %12.3f      %12.3f");
#endif

  MPI_Init(&argn, &argv);

  init_log(argn, argv);
  init_additional_loggers();

  int size = -1;
  int rank = -1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Type_create_struct(3, block_length, block_displace, block_types, &process_state_type);
  MPI_Type_commit(&process_state_type);
  MPI_Type_size(process_state_type, &process_state_type_size);

  int curr_step_start = 0;
  double initial_value = FINE_MULTIPLIER;

  do {
    double mpi_start = MPI_Wtime();

    int working_size = (curr_step_start + size - 1 < TOTAL_STEPS) ? size : TOTAL_STEPS % size;
#ifndef NO_LOGGING
    LOG(INFO) << working_size << " processes will work now";
#endif

    shared_ptr<RmaGetCommunicator> comm = make_shared<RmaGetCommunicator>(working_size);
    shared_ptr<RmaGetProcess> proc = make_shared<RmaGetProcess>(mpi_start);
    RmaGetController<FixedIterationController> controll(comm, proc);

//     controll._fine_delay = bind(random_delay, placeholders::_1, BASE_DELAY * FINE_MULTIPLIER, BASE_DELAY * FINE_MULTIPLIER * FINE_DELAY_VARIANCE);
//     controll._coarse_delay = bind(random_delay, placeholders::_1, BASE_DELAY * COARSE_MULTIPLIER, BASE_DELAY * COARSE_MULTIPLIER * COARSE_DELAY_VARIANCE);

//     controll._fine_delay = bind(deminishing_delay, placeholders::_1, BASE_DELAY * FINE_MULTIPLIER, FINE_DEMINISH);
//     controll._coarse_delay = bind(deminishing_delay, placeholders::_1, BASE_DELAY * COARSE_MULTIPLIER, COARSE_DEMINISH);

    controll._fine_delay = bind(specific_delay, placeholders::_1, BASE_DELAY * FINE_MULTIPLIER, comm->rank, true);
    controll._coarse_delay = bind(specific_delay, placeholders::_1, BASE_DELAY * COARSE_MULTIPLIER, comm->rank, false);

    if (rank < working_size) {
      proc->fine_val = initial_value;
      proc->coarse_val = initial_value / FINE_MULTIPLIER;
#ifndef NO_LOGGING
      LOG(INFO) << "inital values:\tcoarse=" << std::fixed << std::setprecision(3)
                << proc->coarse_val << "\tfine=" << proc->fine_val;
#endif

      controll.run();
    } else {
      // this rank hasn't to do anything anymore
#ifndef NO_LOGGING
      LOG(WARNING) << "hasn't work anymore";
#endif
    }

    initial_value = proc->fine_val;

    curr_step_start += size;
  } while (curr_step_start < TOTAL_STEPS);

  MPI_Type_free(&process_state_type);
#ifndef NO_LOGGING
  VLOG(6) << "finalizing";
#endif
  MPI_Finalize();
}
