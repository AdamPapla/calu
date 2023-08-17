#include <mpi.h>
#ifndef MATRIXIO_H_

void read_chunk_from_file(const char* filename, double* A, int m, int chunkStart, int chunkEnd, int nprocs, int rank, MPI_Comm comm);
void readArray(const char *filename, double *a, int rows, int columns);
void print_mat(double *A, int n, int m);
void writeToFile(const char* filename, double *A, int n, int m);
void printFull(double *A, int n, int m, int nprocs, int rank, MPI_Comm comm);
void intPrintMat(int *A, int m, int n);
void intPrintFull(int *A, int n, int m, int nprocs, int rank, MPI_Comm comm);
#endif
