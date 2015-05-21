#include <cstdlib>
#include <iostream>
#include <string>
using namespace std;

#include <mpi.h>


inline static void check_mpi_error(const int err_code)
{
  if (err_code != MPI_SUCCESS) {
    char err_str[MPI_MAX_ERROR_STRING];
    int err_len = 0;
    MPI_Error_string(err_code, err_str, &err_len);
    cout << "MPI Error: " + string(err_str, err_len) + " (code=" + to_string(err_code) + ")" << endl;
    MPI_Abort(MPI_COMM_WORLD, err_code);
  }
}

int main(int argn, char** argv)
{
  for(int i = 0; i < argn; ++i) { cout << "i: " << i << " '" << argv[i] << "'" << endl; }
  int num_loops = atoi(argv[1]);

  MPI_Init(&argn, &argv);
  int rank = 0, size = 0;
  int err = MPI_SUCCESS;

  err = MPI_Comm_size(MPI_COMM_WORLD, &size); check_mpi_error(err);
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank); check_mpi_error(err);
  cout << "rank " << rank << " / " << size << ": init done" << endl;
  if (rank == 0) {
    double data = 42.21;
    for (int i = 0; i < num_loops + 1; ++i) {
      cout << "rank " << rank << " / " << size << ": sending " << data << endl;
      err = MPI_Send(&data, 1, MPI_DOUBLE, 1, 0 + i, MPI_COMM_WORLD);
      check_mpi_error(err);
      cout << "rank " << rank << " / " << size << ": done" << endl;
      data *= (i + 1);
    }
  } else {
    double data = 0.0;
    MPI_Status stat;
    for (int i = 0; i < num_loops + 1; ++i) {
      data = 0.0;
      stat.MPI_ERROR = MPI_SUCCESS; stat.MPI_SOURCE = MPI_ANY_SOURCE; stat.MPI_TAG = MPI_ANY_TAG;
      cout << "rank " << rank << " / " << size << ": recieving" << endl;
      err = MPI_Recv(&data, 1, MPI_DOUBLE, 0, 0 + i, MPI_COMM_WORLD, &stat);
      check_mpi_error(err);
      cout << "rank " << rank << " / " << size << ": recieved " << data << endl;
    }
  }
  cout << "rank " << rank << " / " << size << ": finalizing";
  MPI_Finalize();
}
