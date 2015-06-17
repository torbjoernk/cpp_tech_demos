#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

#include <mpi.h>


bool check(const int val) { return (val >= 10); }

int compute(const int val, const int rank) { return val + (rank + 1); }


int main(int argn, char** argc)
{
  int size = 0, rank = 0;

  MPI_Init(&argn, &argc);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  bool first = (rank == 0);
  bool last = (rank == size - 1);

  int outgoing_to_rank = rank + 1;
  int incoming_from_rank = rank - 1;

  int outgoing_val = 0;
  int incoming_val = 0;

  size_t iters = 0;

  MPI_Win win;
  int err = MPI_SUCCESS;
  string out_base = to_string(rank) + "/" + to_string(size) + ": ";

  // if not last rank, open window for reading data
  if (!last) {
    err = MPI_Win_create(&outgoing_val, sizeof(int), sizeof(int),
                         MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    assert(err == MPI_SUCCESS);
  } else {
    // no read from last process
    err = MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    assert(err == MPI_SUCCESS);
  }

  for (; !check(incoming_val); ++iters) {
    cout << (out_base + "start iter " + to_string(iters)).c_str() << endl;

    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, incoming_from_rank, 0, win);
    MPI_Get(&incoming_val, 1, MPI_INT,
            incoming_from_rank, 0, 1, MPI_INT,
            win);
    MPI_Win_unlock(incoming_from_rank, win);

    // do the local computation and write into outgoing variable
    MPI_Win_lock(MPI_LOCK_EXCLUSIVE, outgoing_to_rank, 0, win);
    outgoing_val = compute(incoming_val, rank);
    MPI_Win_unlock(outgoing_to_rank, win);
  }

  // finalize RMA windows
  if (!last) {
    err = MPI_Win_unlock(incoming_from_rank, win_read);
    assert(err == MPI_SUCCESS);
  }
  MPI_Win_free(&win_read);
  MPI_Win_free(&win_write);

  cout << (out_base + "done after " + to_string(iters) + " iters: " + to_string(outgoing_val)).c_str() << endl;

  MPI_Finalize();
}
