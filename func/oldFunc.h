#include <mpi.h>
#ifndef OLDFUNC_H_
void distributeMatrix(double *A, double *locA, int n, int m, int b, int myid, int nprocs, int s, int e);
void extractU12(double *locA, double *blockU, int m, int b, int blockNum, int locBlockNum);
void insertU12(double *locA, double *blockU, int m, int b, int blockNum, int locBlockNum);
void pivotTransform(int **pivot, int n, int b);
void readArray(const char *filename, double *a, int rows, int columns);
#endif
