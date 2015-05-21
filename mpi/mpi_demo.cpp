#include <iostream>
using namespace std;

#include <mpi.h>

#include "datatype.hpp"


int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int size = 0, rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  cout << mpi_demo::mpi_start_log() << "Hello" << endl;

  mpi_demo::MPIDataType* dt = new mpi_demo::MPIDataType(2);
  int dest = (rank < size - 1) ? rank + 1 : 0;
  int src = (rank > 0) ? rank - 1 : size - 1;
  dt->send(dest, 0);
  dt->recv(src, 0);
  delete dt;

  MPI_Finalize();
}
