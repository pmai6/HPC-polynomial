/* File: main.c
   Author: Jharrod LaFon
   Date: Spring 2011
   Purpose: Compute the prefix sum of an array
   */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<mpi.h>
#include<math.h>

#define ARRAY_SIZE 30

int prefix_sum(int sum, int rank, int size, MPI_Comm comm);
long prefix_product(long prod, int rank, int size, MPI_Comm comm);
int prefix_steps(int size);

int main( int argc, char *argv[] ) {
  int rank;
  int size;

  if(MPI_Init(&argc, &argv) != MPI_SUCCESS)
  {
      fprintf(stderr, "Unable to initialize MPI!\n");
      return -1;
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if(ARRAY_SIZE % size != 0 && rank == 0)
  {
      fprintf(stderr,"Array size must be multiple of mpi job size.\n");
      return -1;
  }

  MPI_Status status;

  int block_size = ARRAY_SIZE/size;

  int * array = (int *) malloc(sizeof(int) * ARRAY_SIZE);
  int * chunk = (int *) malloc(sizeof(int) * block_size);

// initialize a dummy array filled with random values
  int i = 0;
  int total_sum = 0;
  long total_prod = 1;
  for(i = 0; i < ARRAY_SIZE; i++)
  {
      // array[i] = rand() % 1024;
      array[i]    =   i+1;
      total_sum   +=  array[i];
      total_prod  *=  (long)(array[i]);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Scatter(array, block_size, MPI_INT, chunk, block_size, MPI_INT, 0, MPI_COMM_WORLD);

  long prod = 1;
  for( i = 0; i < block_size; i++ )
    prod *= (long)(chunk[i]);

  int sum = 0;
  for( i = 0; i < block_size; i++ )
    sum += chunk[i];


  sum = prefix_sum(sum, rank, size, MPI_COMM_WORLD);
  prod = prefix_product(prod, rank, size, MPI_COMM_WORLD);
  if(rank == size-1)
  {
    fprintf(stderr,"Total: %d\n", sum);
    fprintf(stderr,"Correct Sum: %d\n", total_sum);
    fprintf(stderr,"Total: %ld\n", total_prod);
    fprintf(stderr,"Correct Product: %ld\n", prod);
  }

  free(array);
  free(chunk);
  MPI_Finalize();
  return 0;
}

int prefix_sum(int sum, int rank, int size, MPI_Comm comm) {
  int temp  = 0;
  int key   = 1;
  int i;
  // int log_proc = (int)(log(size)*1.44269504089) + 1;
  int log_proc = prefix_steps(size);
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

long prefix_product(long prod, int rank, int size, MPI_Comm comm) {
  long temp  = 0;
  int key   = 1;
  int i;
  // int log_proc = (int)(log(size)*1.44269504089) + 1;
  int log_proc = prefix_steps(size);
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
