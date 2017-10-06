/*
 * CX 4220 / CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 1
 *
 *  Utility function definitions
 *
 */

/*
 * File:   utils.h
 * Author: samindaw
 * Modified by: ashanker9
 * Created on February 17, 2017, 11:05 PM
 * Modified on Feb 21, 2017, 3.23 PM
 * Prefix Sum and Product functions added
 */

#ifndef UTILS_H
#define UTILS_H

#include <mpi.h>
#include <time.h>

//------------------- Timer Functions (Do not change) -------------------//
/**
 * Update to current time
 * @param t         Time saved in this data structure
 * @param rank      Rank of the processor requesting time
 * @param comm      MPI communication object
 */
void set_time(double &t, const int rank, MPI_Comm comm);

/**
 * Find the time durations and return the the result in seconds
 * @param t_start   Duration start time
 * @param t_end     Duration end time
 * @return          Duration in seconds
 */
double get_duration(double &t_start, double &t_end);
//---------------------------------------------------------------------//

/*********************************************************************
 *                  Prefix Sum and Product fucntions                  *
 *********************************************************************/

/**
* Calculates prefix sum and returns appropriate values at each node
* @param sum : the locally calculated sum at current node
* @param rank : rank of current node
* @param size : total no. of worker processes
* @return : sum upto the current node
*/
double prefix_sum(double sum, int rank, int size, MPI_Comm comm);

/**
* Calculates prefix product and returns appropriate values at each node
* @param sum : the locally calculated sum at current node
* @param rank : rank of current node
* @param size : total no. of worker processes
* @return : product upto the current node
*/

double prefix_product(double prod, int rank, int size, MPI_Comm comm);
/**
* Calculates required number of steps for prefix sum and product
* @param size : total no. of worker processes
* @return : no. of steps required
*/
int prefix_steps(int size);

// ...

#endif /* UTILS_H */
