#include <mpi.h>
#ifndef	RANDOM_GEN_H_

void createRandomMatrix(double *A, int n, int m);
void writeToFileMPI(const char* filename, double* A, int rows, int cols, MPI_Comm comm);
void generateRMFile(char *fileA, char *fileb, int n, int m, int txt);

#endif
