#include <mpi.h>
#ifndef PIVOT_H_

void pivotMatrix(double *A, int *pivot, int b, int m);
void pivotFullRow(double *A, double *indices, int n, int m, int b, int locBlockNum, MPI_Comm comm);

#endif


