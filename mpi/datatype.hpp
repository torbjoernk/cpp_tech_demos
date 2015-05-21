#include <string>
#include <vector>
using namespace std;

#include <mpi.h>


namespace mpi_demo
{

  static inline string mpi_start_log(int size = -1, int rank = -1)
  {
    if (size == -1) {
      MPI_Comm_size(MPI_COMM_WORLD, &size);
    }
    if (rank == -1) {
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }
    return string("[" + to_string(rank) + "/" + to_string(size) + "]: ");
  }

  class MPIDataType
  {
    protected:
      vector<double> _vec1_src;
      vector<double> _vec2_src;
      vector<double> _vec1_dest;
      vector<double> _vec2_dest;

      int rank = -1;
      int size = -1;

      vector<MPI_Request> send_req;

      static inline string print_vec(vector<double> vec)
      {
        const size_t size = vec.size();
        string out = "[";
        for (size_t i = 0; i < size - 1; ++i) {
          out += to_string(vec[i]) + ", ";
        }
        out += to_string(vec.back()) + "]";
        return out;
      }

    public:
      MPIDataType(const size_t size)
        :   _vec1_src(size), _vec2_src(size)
          , _vec1_dest(size), _vec2_dest(size)
          , send_req(2)
      {
        for (auto req : this->send_req) {
          req = MPI_REQUEST_NULL;
        }
        MPI_Comm_size(MPI_COMM_WORLD, &(this->size));
        MPI_Comm_rank(MPI_COMM_WORLD, &(this->rank));
        for (size_t i = 0; i < size; ++i) {
          this->_vec1_src[i] = this->rank * 100 + i;
          this->_vec2_src[i] = this->rank * 10 + i;
        }
        cout << mpi_start_log() << "datatype with " << size << " elems created" << endl;
      }

      ~MPIDataType()
      {
        cout << mpi_start_log() << "waiting for send to complete" << endl;
        MPI_Status stat_send;
        for (auto req : this->send_req) {
          MPI_Wait(&(req), &stat_send);
        }
        cout << mpi_start_log() << "cleanup. "
                                << "vec1: " << print_vec(_vec1_src) << ", "
                                << "vec2: " << print_vec(_vec2_src) << " --> "
                                << "vec1: " << print_vec(_vec1_dest) << ", "
                                << "vec2: " << print_vec(_vec2_dest) << endl;
      }

      void send(const int dest, const int tag)
      {
        int err = MPI_Isend(this->_vec1_src.data(), sizeof(double) * this->_vec1_src.size(),
                            MPI_CHAR, dest, tag, MPI_COMM_WORLD, &(this->send_req[0]));
        cout << mpi_start_log() << "Isend to " << dest << " with tag " << tag << ": " << err << endl;
        err = MPI_Isend(this->_vec2_src.data(), sizeof(double) * this->_vec2_src.size(),
                        MPI_CHAR, dest, tag + 1, MPI_COMM_WORLD, &(this->send_req[1]));
        cout << mpi_start_log() << "Isend to " << dest << " with tag " << tag + 1 << ": " << err << endl;
      }

      void recv(const int src, const int tag)
      {
        MPI_Status stat;
        int err = MPI_Recv(this->_vec1_dest.data(), sizeof(double) * this->_vec1_dest.size(),
                           MPI_CHAR, src, tag, MPI_COMM_WORLD, &stat);
        cout << mpi_start_log() << "Irecv from " << src << " with tag " << tag << ": " << err << endl;
        err = MPI_Recv(this->_vec2_dest.data(), sizeof(double) * this->_vec2_dest.size(),
                       MPI_CHAR, src, tag + 1, MPI_COMM_WORLD, &stat);
        cout << mpi_start_log() << "Irecv from " << src << " with tag " << tag + 1 << ": " << err << endl;
      }
  };
}  // ::mpi_demo
