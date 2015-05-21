#include <iostream>
using namespace std;

#include <mpi.h>

#include "datatype.hpp"


int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  int size_world = -1, rank_world = -1;
  MPI_Comm_size(MPI_COMM_WORLD, &size_world);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_world);

  int num_groups = 2;
  if (argc == 2) {
    num_groups = *argv[1] - '0';
  }
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Comm MPI_COMM_FIRST;
  int color1 = rank_world / (size_world / num_groups);
  int split_err = MPI_Comm_split(MPI_COMM_WORLD, color1, rank_world, &MPI_COMM_FIRST);
  int size_first = -1, rank_first = -1;
  MPI_Comm_size(MPI_COMM_FIRST, &size_first);
  MPI_Comm_rank(MPI_COMM_FIRST, &rank_first);
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Comm MPI_COMM_SECOND;
  int color2 = rank_world % size_first;
  int split_err2 = MPI_Comm_split(MPI_COMM_WORLD, color2, rank_world, &MPI_COMM_SECOND);
  int size_second = -1, rank_second = -1;
  MPI_Comm_size(MPI_COMM_SECOND, &size_second);
  MPI_Comm_rank(MPI_COMM_SECOND, &rank_second);

  int dummy_buff = rank_world;
  if (rank_world == 0) {
    cout << "s(w)\tr(w)\t-\tc(1)\ts(1)\tr(1)\t-\tc(2)\ts(2)\tr(2)" << endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  cout << size_world << "\t" << rank_world << "\t-\t"
       << color1 << "\t" << size_first << "\t" << rank_first << "\t-\t"
       << color2 << "\t" << size_second << "\t" << rank_second << endl;

  MPI_Finalize();
}
