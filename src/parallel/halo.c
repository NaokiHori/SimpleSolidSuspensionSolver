#include <mpi.h>
#include "common.h"
#include "parallel.h"


int parallel_communicate_halo_with_ymrank(const parallel_t *parallel, const MPI_Datatype mpi_dtype, const void *sendbuf, void *recvbuf){
  const int sendtag = 0;
  const int recvtag = 0;
  const int ymrank = parallel->ymrank;
  const int yprank = parallel->yprank;
  MPI_Sendrecv(sendbuf, 1, mpi_dtype, yprank, sendtag, recvbuf, 1, mpi_dtype, ymrank, recvtag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  return 0;
}

int parallel_communicate_halo_with_yprank(const parallel_t *parallel, const MPI_Datatype mpi_dtype, const void *sendbuf, void *recvbuf){
  const int sendtag = 0;
  const int recvtag = 0;
  const int ymrank = parallel->ymrank;
  const int yprank = parallel->yprank;
  MPI_Sendrecv(sendbuf, 1, mpi_dtype, ymrank, sendtag, recvbuf, 1, mpi_dtype, yprank, recvtag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  return 0;
}

