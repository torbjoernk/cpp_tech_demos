#ifndef _MPI__SIMPLE_COMM_CONFIG_HPP_
#define _MPI__SIMPLE_COMM_CONFIG_HPP_

#include <cassert>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>
#include <functional>
#include <memory>
#include <random>
using namespace std;

typedef chrono::system_clock Clock;
typedef chrono::nanoseconds ClockResolution;

#include <boost/format.hpp>
boost::format log_fmt;

#include <mpi.h>

#define WITH_MPI
#include "../../logging.hpp"

#define MAX_ITER                   5
#define BASE_DELAY              1000  // nanoseconds
#define FINE_MULTIPLIER       200000
#define COARSE_MULTIPLIER      10000
#define STATE_MULTIPLIER          10

#define FINE_DELAY_VARIANCE     0.50  // as percent of default delay
#define COARSE_DELAY_VARIANCE   0.25  // as percent of default delay
#define FINE_DEMINISH           0.10  // as percent
#define COARSE_DEMINISH         0.01  // as percent

#define TOTAL_STEPS                4
#define RESIDUAL_TOL               2  // seconds


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


inline static int fine_tag(const int iter)   { return (iter + 1) * FINE_MULTIPLIER; }
inline static int coarse_tag(const int iter) { return (abs(iter) + 1) * (COARSE_MULTIPLIER / 10); }
inline static int state_tag(const int iter)  { return (iter + 1) * STATE_MULTIPLIER; }


inline static void default_delay(const int iter, const long delay)
{
#ifndef NO_LOGGING
  VLOG(2) << "waiting for " << delay << " nanoseconds";
#endif
  chrono::time_point<Clock> start, end;
  ClockResolution duration;
  start = Clock::now();
  do {
    end = Clock::now();
    duration = end - start;
  } while(duration.count() < delay);
}

inline static void random_delay(const int iter, const long base_delay, const long variance)
{
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<long> dis(base_delay - variance, base_delay + variance);
  default_delay(iter, dis(gen));
}

inline static void deminishing_delay(const int iter, const long base_delay, const double deminish)
{
  if (iter >= 0) {
    default_delay(iter, max((long)BASE_DELAY, (long)(base_delay - base_delay * log10(2 * deminish + pow(iter, deminish)))));
  } else {
    default_delay(iter, base_delay);
  }
}

inline static void specific_delay(const int iter, const long base_delay, const int rank, const bool fine)
{
  if (iter == 3 && fine && rank == 2) {
    default_delay(iter, 1.5 * base_delay);
  } else {
    default_delay(iter, base_delay);
  }
}

inline static void default_fine_delay(const int iter=-1)
{
  default_delay(iter, BASE_DELAY * FINE_MULTIPLIER);
}

inline static void default_coarse_delay(const int iter=-1)
{
  default_delay(iter, BASE_DELAY * COARSE_MULTIPLIER);
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

    Process(const double start_time)
      : coarse_val(0.0), fine_val(0.0), rank(-1), mpi_start(start_time),
        coarse_val_in(0.0), coarse_val_out(0.0), fine_val_in(0.0), fine_val_out(0.0)
    {
      MPI_Comm_rank(MPI_COMM_WORLD, &(this->rank));
    }

    virtual void comp_fine(function<void(int)> delay)
    {
#ifndef NO_LOGGING
      CVLOG(2, "Process") << "start computation";
#endif
      this->fine_val += (this->rank + 1) * FINE_MULTIPLIER + this->state.iter * 0.001;
#ifndef NO_LOGGING
      CVLOG(3, "Process") << this->fine_val << " = " << this->fine_val << " + "
                          << ((this->state.iter + 1) * FINE_MULTIPLIER)
                          << " + " << (this->state.iter * 0.001);
#endif

      delay(this->state.iter);
#ifndef NO_LOGGING
      CVLOG(2, "Process") << "done computation";
#endif
    }

    virtual void comp_coarse(function<void(int)> delay)
    {
#ifndef NO_LOGGING
      CVLOG(2, "Process") << "start computation";
#endif
      this->coarse_val += (this->rank + 1) * COARSE_MULTIPLIER + this->state.iter * 0.001;
#ifndef NO_LOGGING
      CVLOG(3, "Process") << this->coarse_val << " = " << this->coarse_val << " + "
                          << ((this->state.iter + 1) * COARSE_MULTIPLIER)
                          << " + " << (this->state.iter * 0.001);
#endif

      delay(this->state.iter);

#ifndef NO_LOGGING
      CVLOG(2, "Process") << "done computation";
#endif
    }

    virtual void check_state()
    {
      double curr_time = MPI_Wtime();
      this->state.residual = curr_time - this->mpi_start;

      if (this->prev_state.state == PState::FAILED) {
          this->state.state = PState::FAILED;
      } else if (this->prev_state.state == PState::CONVERGED) {
#ifndef NO_LOGGING
        CVLOG(2, "Process") << "previous converged";
#endif
        if (this->state.residual > RESIDUAL_TOL) {
#ifndef NO_LOGGING
          CVLOG(2, "Process") << "and I'm done as well";
#endif
          this->state.state = PState::CONVERGED;
        } else {
#ifndef NO_LOGGING
          CVLOG(2, "Process") << "but I'm not yet done";
#endif
        }
      }
    }
};


class Communicator
{
  public:
    int size = -1;
    int rank = -1;
    bool iam_first;
    bool iam_last;
    int prev = -1;
    int next = -1;

    MPI_Status state_stat;
    MPI_Status fine_stat;
    MPI_Status coarse_stat;

    int mpi_err = MPI_SUCCESS;

    Communicator(const int size)
      : size(size)
    {
      MPI_Comm_rank(MPI_COMM_WORLD, &(this->rank));
      assert(rank >= 0 && size > 0);
      this->iam_first = (rank == 0);
      this->iam_last = (rank == size - 1);
      if (!this->iam_first) { this->prev = this->rank - 1; }
      if (!this->iam_last) { this->next = this->rank + 1; }
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
#ifndef NO_LOGGING
      CVLOG(4, "Communicator") << "broadcasting final value from " << (this->size - 1) << " to all";
#endif
      mpi_err = MPI_Bcast(&(proc->fine_val), 1, MPI_DOUBLE, this->size - 1, MPI_COMM_WORLD);
      assert(mpi_err == MPI_SUCCESS);
    }

    virtual void cleanup() {}
};


class Controller
{
  public:
    shared_ptr<Communicator> comm;
    shared_ptr<Process> proc;
    function<void(int)> _fine_delay;
    function<void(int)> _coarse_delay;

    Controller(shared_ptr<Communicator> comm, shared_ptr<Process> proc)
      : comm(comm), proc(proc),
        _fine_delay(bind(default_fine_delay, -1)), _coarse_delay(bind(default_coarse_delay, -1))
    {}

    virtual void do_fine()
    {
#ifndef NO_LOGGING
      CVLOG(2, "Controller") << "doing fine";
#endif
      this->comm->recv_fine(this->proc);
      this->comm->pre_comp_fine(this->proc);
      this->comm->set_pstate(this->proc, PState::ITER_FINE);
      this->proc->comp_fine(this->_fine_delay);
      this->comm->post_comp_fine(this->proc);
      this->comm->send_fine(this->proc);
    }

    virtual void do_coarse()
    {
#ifndef NO_LOGGING
      CVLOG(2, "Controller") << "doing coarse";
#endif
      this->comm->recv_coarse(this->proc);
      this->comm->pre_comp_coarse(this->proc);
      this->comm->set_pstate(this->proc, PState::ITER_COARSE);
      this->proc->comp_coarse(this->_coarse_delay);
      this->comm->post_comp_coarse(this->proc);
      this->comm->send_coarse(this->proc);
    }

    virtual void check_state()
    {
#ifndef NO_LOGGING
      CVLOG(2, "Controller") << "checking state";
#endif
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
          for (int p = - this->comm->rank; p <= 0; ++p) {
#ifndef NO_LOGGING
            CVLOG(2, "Controller") << "predict " << p;
#endif
            this->comm->set_iter(this->proc, p);
            this->do_coarse();
#ifndef NO_LOGGING
            log_fmt % iter % this->proc->state.residual % this->proc->coarse_val % this->proc->fine_val;
            CLOG(INFO, "Controller") << "PREDICT " << log_fmt;
#endif
          }
        } else {
          this->comm->set_pstate(this->proc, PState::ITERATING);
          this->comm->set_iter(this->proc, iter);
          this->do_coarse();
        }

        this->do_fine();

        this->check_state();
#ifndef NO_LOGGING
        log_fmt % iter % this->proc->state.residual % this->proc->coarse_val % this->proc->fine_val;
        CLOG(INFO, "Controller") << "ITERATE " << log_fmt;
#endif
      }

      this->comm->cleanup();

#ifndef NO_LOGGING
      log_fmt % iter % this->proc->state.residual % this->proc->coarse_val % this->proc->fine_val;
      CLOG(INFO, "Controller") << "FINISHED" << log_fmt;
#endif

      this->comm->bcast_fine(this->proc);

      return iter;
    }
};


class FixedIterationController
  : public Controller
{
  public:
    int max_iter = 10;

    FixedIterationController(shared_ptr<Communicator> comm, shared_ptr<Process> proc)
      : Controller(comm, proc)
    {}

    int run() override
    {
      int iter = 0;
      for(; iter < this->max_iter; ++iter) {
        if (iter == 0) {
          this->comm->set_pstate(this->proc, PState::PREDICTING);
          for (int p = - this->comm->rank; p <= 0; ++p) {
#ifndef NO_LOGGING
            CVLOG(2, "Controller") << "predict " << p;
#endif
            this->comm->set_iter(this->proc, p);
            this->do_coarse();
#ifndef NO_LOGGING
            log_fmt % iter % this->proc->state.residual % this->proc->coarse_val % this->proc->fine_val;
            CLOG(INFO, "Controller") << "PREDICT " << log_fmt;
#endif
          }
        } else {
          this->comm->set_pstate(this->proc, PState::ITERATING);
          this->comm->set_iter(this->proc, iter);
          this->do_coarse();
        }

        this->do_fine();

        this->check_state();
#ifndef NO_LOGGING
        log_fmt % iter % this->proc->state.residual % this->proc->coarse_val % this->proc->fine_val;
        CLOG(INFO, "Controller") << "ITERATE " << log_fmt;
#endif
      }

      this->comm->cleanup();

#ifndef NO_LOGGING
      log_fmt % iter % this->proc->state.residual % this->proc->coarse_val % this->proc->fine_val;
      CLOG(INFO, "Controller") << "FINISHED" << log_fmt;
#endif

      this->comm->bcast_fine(this->proc);

      return iter;
    }
};

#endif // _MPI__SIMPLE_COMM_CONFIG_HPP_
