#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

#include <mpi.h>

int main(int argn, char** argc)
{
  int size = 0, rank = 0;

  MPI_Init(&argn, &argc);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  vector<int> src(size, 0);
  vector<int> dest(size, 0);

  for (size_t i = 0; i < size; ++i) {
    src[i] = i + rank;
  }

  MPI_Win write_win;
  int err = MPI_SUCCESS;
  err = MPI_Win_create(dest.data(), dest.size() * sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &write_win);
  assert(err == MPI_SUCCESS);

  err = MPI_Win_fence(MPI_MODE_NOPRECEDE, write_win);
  assert(err == MPI_SUCCESS);

  err = MPI_Put(src.data(), src.size(), MPI_INT, (rank + 1) % size, 0, src.size(), MPI_INT, write_win);
  assert(err == MPI_SUCCESS);

  err = MPI_Win_fence(0, write_win);
  assert(err == MPI_SUCCESS);

  err = MPI_Win_fence(MPI_MODE_NOSUCCEED, write_win);
  assert(err == MPI_SUCCESS);

  MPI_Win_free(&write_win);

  string out = to_string(rank) + "/" + to_string(size) + ": ";
  for (size_t i = 0; i < size - 1; ++i) {
    out += to_string(dest[i]) + ", ";
  }
  out += to_string(dest[size-1]);
  cout << out.c_str() << endl;

  MPI_Finalize();
}
