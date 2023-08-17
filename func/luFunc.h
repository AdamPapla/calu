#include <mpi.h>
#ifndef MATRIX_H_

void unpivotedLUA(double *locA, int blockNum, int locBlockNum, int block, int m, int b, int rank, int nprocs, MPI_Comm comm);
void finalLUA(double *A, int blockNum, int locBlockNum, int block, int m, int b, int rem);

#endif
