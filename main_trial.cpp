/*
 * CX 4220 / CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 1
 *
 *  Main executor for poly-eval
 *
 */

/*
 * File:   main_serial.cpp
 * Author: Apaar Shanker
 *
 * Created on February 17, 2017, 2:29 PM
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
// #include "mpi_evaluator.h"
#include "evaluator.h"
#include "const.h"
#include "io.h"
#include "utils.h"

using namespace std;

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    int p, rank;
    const MPI_Comm comm = MPI_COMM_WORLD;

    //get the rank and number of processors
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    // Get the input data
    int n, m;
    vector<double> global_constants;
    vector<double> x;
    if (rank == 0) {
        setup(argc, argv, n, global_constants, m, x);
    }

    double poly_evaltime = 0.0;
    double t_start, t_end;
    
    if (rank == 0) {
            //------------------- run the serial version ---------------------//
        for(int i=0; i<m; i++){
            set_time(t_start, -1, comm);
            double result = poly_evaluator(x[i], n, &global_constants[0]); //eval
            set_time(t_end, -1, comm);
            poly_evaltime = get_duration(t_start, t_end);

            print_output(x[i], result, -1, poly_evaltime, 0, 0);            //print result
        }
            //----------------------------------------------------------------//
    } else {
        printf("Hello, World!\n");
    }
    MPI_Finalize();
    return (EXIT_SUCCESS);
}
