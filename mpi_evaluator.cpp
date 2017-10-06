/*
 * CX 4220 / CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 1
 *
 *  MPI polynomial evaluation algorithm function implementations go here
 *
 */

#include "mpi_evaluator.h"
#include "const.h"

// void scatter(const int n, double* scatter_values, int &n_local, double* &local_values, int source_rank, const MPI_Comm comm){
//     //Implementation
//     MPI_Status status;
//     int p; // No. of processes
//     int rank; // rank of current process
//     MPI_Comm_size(comm, &p);
//     MPI_Comm_rank(comm, &rank);
//
//     int i;
//     int dest;
//     int shift;
//     int msg_size;
//
//     if (rank == source_rank) {
//       shift = 0;
//       for ( i = 0; i < p; i++) {
//         if (i == source_rank) {
//           n_local = n/p + n%p;
//           local_values = (double *)(malloc(n_local*sizeof(double)));
//           for (int indx = 0; indx < n_local; indx++ ) {
//             local_values[indx] = scatter_values[shift+indx];
//           }
//           shift += n_local;
//         } else {
//           dest = i;
//           msg_size = n/p;
//           MPI_Send(&scatter_values[shift],msg_size, MPI_DOUBLE,dest,100,comm);
//           shift += msg_size;
//         }
//       }
//     } else {
//       n_local = n/p;
//       local_values = (double *)(malloc(n_local*sizeof(double)));
//       MPI_Recv(&local_values[0],n_local, MPI_DOUBLE,source_rank,100,comm,&status);
//     }
// }

void scatter(const int n, double* scatter_values, int &n_local, double* &local_values, int source_rank, const MPI_Comm comm){
    //Implementation
    int p; // No. of processes
    int rank; // rank of current process
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    int i;
    int dest;
    int shift;
    int msg_size;
    int extra;
    int avg;

    avg = n/p;
    extra = n%p;

    if (rank == source_rank) {
      shift = 0;
      for ( i = 0; i < p; i++ ) {
        dest = i;
        if ( i == source_rank ) {
          msg_size = (dest < extra) ? avg+1 : avg;
          n_local  = msg_size;
          local_values = (double *)(malloc(n_local*sizeof(double)));
          for (int indx = 0; indx < n_local; indx++ ) {
            local_values[indx] = scatter_values[shift+indx];
          }
          shift   += msg_size;
        } else {
          msg_size = (dest < extra) ? avg+1 : avg;
          MPI_Send(&scatter_values[shift],   msg_size, MPI_DOUBLE,  dest,  100, comm);
          shift += msg_size;
        }
      }
    } else {
      n_local = (rank < extra) ? avg+1 : avg;
      local_values = (double *)(malloc(n_local*sizeof(double)));
      MPI_Recv(&local_values[0],  n_local, MPI_DOUBLE,  source_rank, 100, comm, MPI_STATUS_IGNORE);
    }
}

double broadcast(double value, int source_rank, const MPI_Comm comm){
  //Implementation
  int rank;
  MPI_Comm_rank(comm, &rank);
  int size;
  MPI_Comm_size(comm, &size);

  int i;
  int log_proc = parallel_steps(size);
  int flip = 1 << (log_proc-1);
  int mask = flip -1;

  if(rank == source_rank){
    MPI_Send(&value, 1, MPI_DOUBLE,0, 0, comm);
  }
  if(rank == 0){
    MPI_Recv(&value, 1, MPI_DOUBLE, source_rank, 0, comm, MPI_STATUS_IGNORE);
  }
    for (i = 0; i < log_proc; i++) {
      if ((rank & mask) == 0){
        if ((rank & flip) == 0) {
          MPI_Send(&value, 1, MPI_DOUBLE, (rank^flip)%size, 50, comm);
        }else {
          MPI_Recv(&value, 1, MPI_DOUBLE, (rank^flip)%size, 50, comm, MPI_STATUS_IGNORE);
        }
      }
      mask = mask >> 1;
      flip = flip >> 1;
    }

  return value;
}

// implemented in utils library - msg from ashanker9
void parallel_prefix(const int n, const double* values, double* prefix_results, const int OP, const MPI_Comm comm){
    //Implementation

}

double mpi_poly_evaluator(const double x, const int n, const double* constants, const MPI_Comm comm) {

  int p; // No. of processes
  int rank; // rank of current process
  MPI_Comm_size(comm, &p);
  MPI_Comm_rank(comm, &rank);

  double result;
  double v1 = 1.0;
  double factor = 1.0;
  MPI_Status status;
  // power = rank*local_n;

  /*
  * n/(p*p) products calculated at each node
  * product at current node is sent to the right neighbor to achieve following config:
    * at process rank i = 0; factor = 1
    * at other i; factor = x^{(i)*n/p}
  * results in O(n/p) calculations
  */
/******************************************************************************/
// //approach 1
  int powr;
  int avg   = n/p;
  int extra = n%p;
  powr = (rank < extra) ? avg+1 : avg;
  double prod = 1;
  for (int i = 0; i < powr; i++) {
    prod *= x;
  }
  prod = prefix_product(prod, rank, p, comm);
  prod = broadcast(prod, p-1, comm);
  v1 = prod;
  if (rank == 0) {
    MPI_Send(&v1,   1, MPI_DOUBLE,  1,      0, comm);
  } else if (rank == p-1){
    MPI_Recv(&factor,  1, MPI_DOUBLE,  rank-1, 0, comm, &status);
  }else {
    MPI_Recv(&factor,  1, MPI_DOUBLE,  rank-1, 0, comm, &status);
    v1 *= factor;
    MPI_Send(&v1,   1, MPI_DOUBLE,  rank+1, 0, comm);
  }
/******************************************************************************/
// approach 2
/******************************************************************************/
// int power = n;
  // if (rank == 0) {
  //   for ( int i = 0; i < power; i++) {
  //     v1 = v1*x;
  //   }
  //   MPI_Bcast(&v1, 1, MPI_DOUBLE, 0, comm);
  //   MPI_Send(&v1,   1, MPI_DOUBLE,  1,      0, comm);
  // } else if (rank == p-1){
  //   MPI_Bcast(&v1, 1, MPI_DOUBLE, 0, comm);
  //   MPI_Recv(&factor,  1, MPI_DOUBLE,  rank-1, 0, comm, &status);
  // }else {
  //   MPI_Bcast(&v1, 1, MPI_DOUBLE, 0, comm);
  //   MPI_Recv(&factor,  1, MPI_DOUBLE,  rank-1, 0, comm, &status);
  //   v1 *= factor;
  //   MPI_Send(&v1,   1, MPI_DOUBLE,  rank+1, 0, comm);
  // }
/******************************************************************************/

  /*
  * Polynomial is evaluated for x using the sequential algorithm
  * for block of locally avlbl. constants
  */
  result = poly_evaluator(x, n, &constants[0]);
  /*
  * Polynomial is multiplied with prefactor
  */
  result = result*factor;
  /*
  * polynomial chunks evaluated at each node are added using parallel prefix_sum
  * prefix_sum implementation contained in utils.cpp
  */
  result = prefix_sum(result, rank, p, comm);
  result = broadcast(result, p-1, comm);

  return result;
}

int parallel_steps(int size){
  int ans = 0;
  // double s = (double)(size);
  while (size >= 2) {
    size /= 2.;
    ans++;
  }
  ans++;
  return ans;
}
