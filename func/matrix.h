#include <mpi.h>
#ifndef MATRIX_H_

void fillPanel(double *locA, double *locPanel, int panelNum, int locBlockNum, int block, int m, int b, int rank, int rem);
void extractL(double *A, int blockNum, int locBlockNum, int m, int b, int rem, double *L);
void gatherMatrix(double *locA, int n, int m, int b, int block, int rank, int nprocs, MPI_Comm comm, double *A);
void gatherPivot(int *ipvt, int n, int b, int block, int rank, int nprocs, MPI_Comm comm, int *ipvtFull);
void triangSub(double *A, int n, int m, double *b, double *x);
#endif
