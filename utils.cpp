/*
 * CX 4220 / CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 1
 *
 *  Utility function implementation
 *
 */

#include <mpi.h>
#include "utils.h"

//------------------- Timer Functions (Do not change) -------------------//
void set_time(double &t, const int rank, MPI_Comm comm){
    if (rank>=0) // Do not call barrier if rank is negative
        MPI_Barrier(comm);
    if (rank <= 0){ //only 1 processor will set the time
        t = MPI_Wtime();
    }
}

double get_duration(double &t_start, double &t_end){
    return (t_end - t_start);
}
//---------------------------------------------------------------------//

/*********************************************************************
 *                 Implement your own functions here                 *
 *********************************************************************/
 double prefix_sum(double sum, int rank, int size, MPI_Comm comm) {
   double temp  = 0;
   int key   = 1;
   int i;
   int log_proc = prefix_steps(size);
   MPI_Status status;
   for (i = 0; i < log_proc; i++) {
     if (rank <= (size-1-key)) {
       MPI_Send(&sum, 1, MPI_DOUBLE, rank+key, 0, comm);
     }
     if(rank >= key) {
       MPI_Recv(&temp, 1, MPI_DOUBLE, rank-key, 0, comm, &status);
       sum += temp;
     }
     key = key*2;
   }
   return sum;
 }
 /******************************************************************************
 *******************************************************************************
 ******************************************************************************/
 double prefix_product(double prod, int rank, int size, MPI_Comm comm) {
   double temp  = 0;
   int key   = 1;
   int i;
   int log_proc = prefix_steps(size);
   MPI_Status status;
   for (i = 0; i < log_proc; i++) {
     if (rank <= (size-1-key)) {
       MPI_Send(&prod, 1, MPI_DOUBLE, rank+key, 0, comm);
     }
     if(rank >= key) {
       MPI_Recv(&temp, 1, MPI_DOUBLE, rank-key, 0, comm, &status);
       prod *= temp;
     }
     key = key*2;
   }
   return prod;
 }
 /******************************************************************************
 *******************************************************************************
 ******************************************************************************/
 int prefix_steps(int size){
   int ans = 0;
   // double s = (double)(size);
   while (size >= 2) {
     size /= 2.;
     ans++;
   }
   ans++;
   return ans;
 }

// ...
