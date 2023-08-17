#include <mpi.h>
#ifndef TSLU_H_
int tournPivot(double *locA, int n, int b, int locBlockNum, double **indices, MPI_Comm comm);
#endif
