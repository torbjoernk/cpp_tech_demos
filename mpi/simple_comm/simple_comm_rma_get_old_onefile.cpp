#include <cassert>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>
#include <functional>
#include <memory>
using namespace std;

typedef chrono::system_clock Clock;
typedef chrono::nanoseconds ClockResolution;

#include <boost/format.hpp>
boost::format log_fmt;

#include <mpi.h>

#define WITH_MPI
#include "../../logging.hpp"

#define MAX_ITER                 5
#define BASE_DELAY            1000  // nanoseconds
#define FINE_MULTIPLIER     200000
#define COARSE_MULTIPLIER    10000
#define STATE_MULTIPLIER        10

#define TOTAL_STEPS              4
#define RESIDUAL_TOL             1  // seconds


inline static void init_additional_loggers()
{
  add_custom_logger("Process");
  add_custom_logger("Communicator");
  add_custom_logger("Controller");
}


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
  double residual = numeric_limits<double>::min();
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


static void default_fine_delay(const int iter=-1)
{
  chrono::time_point<Clock> start, end;
  ClockResolution duration;
  start = Clock::now();
  do {
    end = Clock::now();
    duration = end - start;
  } while(duration.count() < BASE_DELAY * FINE_MULTIPLIER);
}

static void default_coarse_delay(const int iter=-1)
{
  chrono::time_point<Clock> start, end;
  ClockResolution duration;
  start = Clock::now();
  do {
    end = Clock::now();
    duration = end - start;
  } while(duration.count() < BASE_DELAY * COARSE_MULTIPLIER);
}


class Process
{
  public:
    double coarse_val;
    double coarse_val_in, coarse_val_out;
    double fine_val;
    double fine_val_in, fine_val_out;
    int rank;

    double mpi_start;

    ProcessState state;
    ProcessState prev_state;

    virtual void init(const double start_time)
    {
      this->mpi_start = start_time;
    }

    virtual void comp_fine(function<void(int)> delay = bind(default_fine_delay, -1))
    {
      CVLOG(2, "Process") << "start computation";
      this->fine_val += (this->rank + 1) * FINE_MULTIPLIER + this->state.iter * 0.001;
      CVLOG(3, "Process") << this->fine_val << " = " << this->fine_val << " + "
                          << ((this->state.iter + 1) * FINE_MULTIPLIER)
                          << " + " << (this->state.iter * 0.001);

      delay(this->state.iter);
      CVLOG(2, "Process") << "done computation";
    }

    virtual void comp_coarse(function<void(int)> delay = bind(default_coarse_delay, -1))
    {
      CVLOG(2, "Process") << "start computation";
      this->coarse_val += (this->rank + 1) * COARSE_MULTIPLIER + this->state.iter * 0.001;
      CVLOG(3, "Process") << this->coarse_val << " = " << this->coarse_val << " + "
                          << ((this->state.iter + 1) * COARSE_MULTIPLIER)
                          << " + " << (this->state.iter * 0.001);

      delay(this->state.iter);
      CVLOG(2, "Process") << "done computation";
    }

    virtual void check_state()
    {
      double curr_time = MPI_Wtime();
      this->state.residual = curr_time - this->mpi_start;

      if (this->prev_state.state == PState::FAILED) {
          this->state.state = PState::FAILED;
      } else if (this->prev_state.state == PState::CONVERGED) {
        CVLOG(2, "Process") << "previous converged";
        if (this->state.residual > RESIDUAL_TOL) {
          CVLOG(2, "Process") << "and I'm done as well";
          this->state.state = PState::CONVERGED;
        } else {
          CVLOG(2, "Process") << "but I'm not yet done";
        }
      }
    }
};


class Communicator
{
  public:
    int _size;
    int _rank;
    bool iam_first;
    bool iam_last;
    int prev;
    int next;

    MPI_Status state_stat;
    MPI_Status fine_stat;
    MPI_Status coarse_stat;

    int mpi_err = MPI_SUCCESS;

    virtual void init(const int size)
    {
      this->_size = size;
      MPI_Comm_rank(MPI_COMM_WORLD, &(this->_rank));
      assert(_rank >= 0 && size > 0);
      this->iam_first = (this->_rank == 0);
      this->iam_last = (this->_rank == size - 1);
      if (!this->iam_first) { this->prev = this->_rank - 1; }
      if (!this->iam_last) { this->next = this->_rank + 1; }
      this->state_stat = MPI_Status_factory();
      this->fine_stat = MPI_Status_factory();
      this->coarse_stat = MPI_Status_factory();
    }

    virtual void set_pstate(shared_ptr<Process> proc, const PState new_state)
    {
      this->pre_update_state(proc);
      proc->state.state = new_state;
      this->post_update_state(proc);
    }
    virtual void set_iter(shared_ptr<Process> proc, const int iter)
    {
      this->pre_update_state(proc);
      proc->state.iter = iter;
      this->post_update_state(proc);
    }

    virtual void recv_fine(shared_ptr<Process> proc) {}
    virtual void pre_comp_fine(shared_ptr<Process> proc)
    {
      this->set_pstate(proc, PState::PRE_ITER_FINE);
    }
    virtual void post_comp_fine(shared_ptr<Process> proc)
    {
      this->set_pstate(proc, PState::POST_ITER_FINE);
    }
    virtual void send_fine(shared_ptr<Process> proc) {}

    virtual void recv_coarse(shared_ptr<Process> proc) {}
    virtual void pre_comp_coarse(shared_ptr<Process> proc)
    {
      this->set_pstate(proc, PState::PRE_ITER_COARSE);
    }
    virtual void post_comp_coarse(shared_ptr<Process> proc)
    {
      this->set_pstate(proc, PState::POST_ITER_COARSE);
    }
    virtual void send_coarse(shared_ptr<Process> proc) {}

    virtual void recv_state(shared_ptr<Process> proc) {}
    virtual void pre_update_state(shared_ptr<Process> proc) {}
    virtual void post_update_state(shared_ptr<Process> proc) {}
    virtual void send_state(shared_ptr<Process> proc) {}

    virtual void bcast_fine(shared_ptr<Process> proc) {
      CVLOG(4, "Communicator") << "broadcasting final value to all";
      mpi_err = MPI_Bcast(&(proc->fine_val), 1, MPI_DOUBLE, this->_size - 1, MPI_COMM_WORLD);
      assert(mpi_err == MPI_SUCCESS);
    }

    virtual void cleanup() {}
};


template<class CommT, class ProcT>
class Controller
{
  public:
    shared_ptr<CommT> comm;
    shared_ptr<ProcT> proc;

    explicit Controller(shared_ptr<CommT> comm, shared_ptr<ProcT> proc)
      : comm(comm), proc(proc)
    {}

    virtual void init(const int wsize, const double start_time)
    {
      this->comm->init(wsize);
      this->proc->init(start_time);
    }

    virtual void do_fine()
    {
      CVLOG(2, "Controller") << "doing fine";
      this->comm->recv_fine(this->proc);
      this->comm->pre_comp_fine(this->proc);
      this->comm->set_pstate(this->proc, PState::ITER_FINE);
      this->proc->comp_fine();
      this->comm->post_comp_fine(this->proc);
      this->comm->send_fine(this->proc);
    }

    virtual void do_coarse()
    {
      CVLOG(2, "Controller") << "doing coarse";
      this->comm->recv_coarse(this->proc);
      this->comm->pre_comp_coarse(this->proc);
      this->comm->set_pstate(this->proc, PState::ITER_COARSE);
      this->proc->comp_coarse();
      this->comm->post_comp_coarse(this->proc);
      this->comm->send_coarse(this->proc);
    }

    virtual void check_state()
    {
      CVLOG(2, "Controller") << "checking state";
      this->comm->recv_state(this->proc);
      this->comm->pre_update_state(this->proc);
      this->proc->check_state();
      this->comm->post_update_state(this->proc);
      this->comm->send_state(this->proc);
    }

    virtual int run()
    {
      int iter = 0;
      for(; this->proc->state.state > PState::FAILED; ++iter) {
        if (iter == 0) {
          this->comm->set_pstate(this->proc, PState::PREDICTING);
          for (int p = - this->comm->_rank; p <= 0; ++p) {
            CVLOG(2, "Controller") << "predict " << p;
            this->comm->set_iter(this->proc, p);
            this->do_coarse();
            log_fmt % iter % this->proc->state.residual % this->proc->coarse_val % this->proc->fine_val;
            CLOG(INFO, "Controller") << "PREDICT " << log_fmt;
          }
        } else {
          this->comm->set_pstate(this->proc, PState::ITERATING);
          this->comm->set_iter(this->proc, iter);
          this->do_coarse();
        }

        this->do_fine();

        this->check_state();
        log_fmt % iter % this->proc->state.residual % this->proc->coarse_val % this->proc->fine_val;
        CLOG(INFO, "Controller") << "ITERATE " << log_fmt;
      }

      this->comm->cleanup();

      log_fmt % iter % this->proc->state.residual % this->proc->coarse_val % this->proc->fine_val;
      CLOG(INFO, "Controller") << "FINISHED" << log_fmt;

      this->comm->bcast_fine(this->proc);

      return iter;
    }
};


class RmaGetProcess
  : public Process
{};


class RmaGetCommunicator
  : public Communicator
{
  public:
    MPI_Win coarse_win;
    MPI_Win fine_win;
    MPI_Win state_win;

    RmaGetCommunicator()
      : Communicator()
    {
      this->coarse_win = MPI_WIN_NULL;
      this->fine_win = MPI_WIN_NULL;
      this->state_win = MPI_WIN_NULL;
    }

    void recv_fine(shared_ptr<Process> proc) override
    {
      if (this->_size > 1 && !this->iam_first && proc->state.iter > 0) {
        CVLOG(5, "Communicator") << "locking fine window to rank " << this->prev;
        mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, this->prev, 0, this->fine_win);
        assert(mpi_err == MPI_SUCCESS);

        CVLOG(4, "Communicator") << "getting fine value from " << this->prev;
        mpi_err = MPI_Get(&(proc->fine_val_in), 1, MPI_DOUBLE, this->prev, 0, 1, MPI_DOUBLE, this->fine_win);
        assert(mpi_err == MPI_SUCCESS);

        mpi_err = MPI_Win_unlock(this->prev, this->fine_win);
        assert(mpi_err == MPI_SUCCESS);
        CVLOG(5, "Communicator") << "unlocked fine window to rank " << this->prev;
      }
    }

    void pre_comp_fine(shared_ptr<Process> proc) override
    {
      Communicator::pre_comp_fine(proc);

      if (this->_size > 1) {
        CVLOG(5, "Communicator") << "locking local fine window";
        mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, this->_rank, 0, this->fine_win);
        assert(mpi_err == MPI_SUCCESS);
      }

      proc->fine_val = proc->fine_val_in;
    }

    void post_comp_fine(shared_ptr<Process> proc) override
    {
      Communicator::post_comp_fine(proc);

      proc->fine_val_out = proc->fine_val;

      if (this->_size > 1) {
        mpi_err = MPI_Win_unlock(this->_rank, this->fine_win);
        assert(mpi_err == MPI_SUCCESS);
        CVLOG(5, "Communicator") << "unlocked local fine window";
      }
    }

    void recv_coarse(shared_ptr<Process> proc) override
    {
      if (this->_size > 1 && !this->iam_first && proc->state.iter != - this->_rank) {
        CVLOG(5, "Communicator") << "locking coarse window to rank " << this->prev;
        mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, this->prev, 0, this->coarse_win);
        assert(mpi_err == MPI_SUCCESS);

        CVLOG(4, "Communicator") << "getting coarse value from " << this->prev;
        mpi_err = MPI_Get(&(proc->coarse_val_in), 1, MPI_DOUBLE, this->prev, 0, 1, MPI_DOUBLE, this->coarse_win);
        assert(mpi_err == MPI_SUCCESS);

        mpi_err = MPI_Win_unlock(this->prev, this->coarse_win);
        assert(mpi_err == MPI_SUCCESS);
        CVLOG(5, "Communicator") << "unlocked coarse window to rank " << this->prev;
      }
    }

    void pre_comp_coarse(shared_ptr<Process> proc) override
    {
      Communicator::pre_comp_coarse(proc);

      if (this->_size > 1) {
        CVLOG(5, "Communicator") << "locking local coarse window";
        mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, this->_rank, 0, this->coarse_win);
        assert(mpi_err == MPI_SUCCESS);
      }

      proc->coarse_val = proc->coarse_val_in;
    }

    void post_comp_coarse(shared_ptr<Process> proc) override
    {
      Communicator::post_comp_coarse(proc);

      proc->coarse_val_out = proc->coarse_val;

      if (this->_size > 1) {
        mpi_err = MPI_Win_unlock(this->_rank, this->coarse_win);
        assert(mpi_err == MPI_SUCCESS);
        CVLOG(5, "Communicator") << "unlocked local coarse window";
      }
    }

    void recv_state(shared_ptr<Process> proc) override
    {
      if (this->_size > 1 && !this->iam_first) {
        CVLOG(5, "Communicator") << "locking state window to rank " << this->prev;
        mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, this->prev, 0, this->state_win);
        assert(mpi_err == MPI_SUCCESS);

        CVLOG(4, "Communicator") << "getting state from " << this->prev;
        mpi_err = MPI_Get(&(proc->prev_state), 1, process_state_type, this->prev, 0, 1,
                          process_state_type, this->state_win);
        assert(mpi_err == MPI_SUCCESS);

        mpi_err = MPI_Win_unlock(this->prev, this->state_win);
        assert(mpi_err == MPI_SUCCESS);
        CVLOG(5, "Communicator") << "unlocked state window to rank " << this->prev;
      } else {
        proc->prev_state.state = PState::CONVERGED;
      }
    }

    void pre_update_state(shared_ptr<Process> proc) override
    {
      if (this->_size > 1) {
        CVLOG(5, "Communicator") << "locking local state window";
        mpi_err = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, this->_rank, 0, this->state_win);
        assert(mpi_err == MPI_SUCCESS);
      }
    }

    void post_update_state(shared_ptr<Process> proc) override
    {
      if (this->_size > 1) {
        mpi_err = MPI_Win_unlock(this->_rank, this->state_win);
        assert(mpi_err == MPI_SUCCESS);
        VLOG(5) << "unlocked local state window";
      }
    }

    void cleanup() override
    {
      if (this->_size > 1) {
        CVLOG(6, "Communicator") << "freeing windows";
        MPI_Win_free(&(this->fine_win));
        MPI_Win_free(&(this->coarse_win));
        MPI_Win_free(&(this->state_win));
      }
    }
};


class RmaGetController
  : public Controller<RmaGetCommunicator, RmaGetProcess>
{
  public:
    explicit RmaGetController(shared_ptr<RmaGetCommunicator> comm, shared_ptr<RmaGetProcess> proc)
      : Controller<RmaGetCommunicator, RmaGetProcess>(comm, proc)
    {
      if (this->comm->_size > 1) {
        CVLOG(6, "Controller") << "creating windows";
        MPI_Win_create(&(this->proc->fine_val_out), sizeof(double), MPI_DOUBLE, MPI_INFO_NULL,
                       MPI_COMM_WORLD, &(this->comm->fine_win));
        MPI_Win_create(&(this->proc->coarse_val_out), sizeof(double), MPI_DOUBLE, MPI_INFO_NULL,
                       MPI_COMM_WORLD, &(this->comm->coarse_win));
        MPI_Win_create(&(this->proc->state), 1, process_state_type_size, MPI_INFO_NULL, MPI_COMM_WORLD,
                       &(this->comm->state_win));
      }
    }
};


int main(int argn, char** argv) {
  //                       iter    residual    coarse      fine
  log_fmt = boost::format("%4.d    %12.6f      %12.3f      %12.3f");

  MPI_Init(&argn, &argv);

  init_log(argn, argv);

  int size = -1;
  int rank = -1;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Type_create_struct(3, block_length, block_displace, block_types, &process_state_type);
  MPI_Type_commit(&process_state_type);
  MPI_Type_size(process_state_type, &process_state_type_size);

  int curr_step_start = 0;
  double initial_value = FINE_MULTIPLIER;

  shared_ptr<RmaGetCommunicator> comm = make_shared<RmaGetCommunicator>();
  shared_ptr<RmaGetProcess> proc = make_shared<RmaGetProcess>();
  RmaGetController controll(comm, proc);

  do {
    double mpi_start = MPI_Wtime();

    int working_size = (curr_step_start + size - 1 < TOTAL_STEPS) ? size : TOTAL_STEPS % size;
    LOG(INFO) << working_size << " processes will work now";

    controll.init(working_size, mpi_start);

    if (rank < working_size) {
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
