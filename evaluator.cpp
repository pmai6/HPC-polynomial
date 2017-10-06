/*
 * CX 4220 / CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 1
 *
 *  Serial polynomial evaluation algorithm function implementations goes here
 *
 */

double poly_evaluator(const double x, const int n, const double* constants){
    //Implementation
    int i;
    double ans = 0.0;
    ans = constants[n-1];
    for (i = n-1; i > 0; i--) {
      ans = x*ans + constants[i-1];
    }
    return ans;
}
