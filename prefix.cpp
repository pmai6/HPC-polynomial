/*
 * CX 4220 / CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 1
 *
 *  Prefix Sum and Product function implementations goes here
 *
 */

int prefix_sum(int sum, int rank, int size, MPI_Comm comm) {
  int temp  = 0;
  int key   = 1;
  int i;
  int log_proc = (int)(log(size)*1.44269504089) + 1;
  MPI_Status status;
  for (i = 0; i < log_proc; i++) {
    if (rank <= (size-1-key)) {
      MPI_Send(&sum, 1, MPI_INT, rank+key, 0, comm);
    }
    if(rank >= key) {
      MPI_Recv(&temp, 1, MPI_INT, rank-key, 0, comm, &status);
      sum += temp;
    }
    key = key*2;
    MPI_Barrier(comm);
  }

  return sum;
}
/******************************************************************************
*******************************************************************************
******************************************************************************/
long prefix_product(long prod, int rank, int size, MPI_Comm comm) {
  long temp  = 0;
  int key   = 1;
  int i;
  int log_proc = (int)(log(size)*1.44269504089) + 1;
  MPI_Status status;
  for (i = 0; i < log_proc; i++) {
    if (rank <= (size-1-key)) {
      MPI_Send(&prod, 1, MPI_LONG, rank+key, 0, comm);
    }
    if(rank >= key) {
      MPI_Recv(&temp, 1, MPI_LONG, rank-key, 0, comm, &status);
      prod *= temp;
    }
    key = key*2;
    MPI_Barrier(comm);
  }

  return prod;
}
