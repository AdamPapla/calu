#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mkl_types.h>
#include <mkl_cblas.h>
#include <mkl.h>
#include <mpi.h>

#include "matrixIO.h"

void createRandomMatrix(double *A, int n, int m){
    // Seed the random number generator
    srand(1);

    // Generate and write the matrix elements to the file
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            // Generate a random double between 0 and 1
            A[i*m+j] = (double)rand() / RAND_MAX;
		}
	}
}
void writeToFileMPI(const char* filename, double* A, int rows, int cols, MPI_Comm comm) {
    MPI_Offset filesize = rows * cols * sizeof(double);

    // Open the file collectively
    MPI_File file;
    MPI_File_open(comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);

    // Set the file view for each process to write its chunk of data
    MPI_Offset displacement = 0; // Offset from the beginning of the file
    MPI_Datatype filetype;
    MPI_Type_create_subarray(2, (int[]){rows, cols}, (int[]){rows, cols}, (int[]){0, 0}, MPI_ORDER_C, MPI_DOUBLE, &filetype);
    MPI_Type_commit(&filetype);
    MPI_File_set_view(file, displacement, MPI_DOUBLE, filetype, "native", MPI_INFO_NULL);

    // Write the data to the file
    MPI_File_write_all(file, A, rows * cols, MPI_DOUBLE, MPI_STATUS_IGNORE);

    // Close the file
    MPI_File_close(&file);
}

void generateRMFile(char *fileA, char *fileb, int n, int m, int txt){
	int rank; int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		
	double *A = malloc(n*m*sizeof(double));
	double *x = malloc(m*sizeof(double));
	double *b = malloc(n*sizeof(double));

	createRandomMatrix(A, n, m);
	for (int i=0; i<m; i++){
		x[i] = 2.0;
	}
	cblas_dgemv(CblasRowMajor, CblasNoTrans, n, m, 1.0, A, m, x, 1, 1.0, b, 1);
	writeToFileMPI(fileA, A, n, m, MPI_COMM_WORLD);
	writeToFileMPI(fileb, b, n, 1, MPI_COMM_WORLD);
	if (txt == 1){
		writeToFile("A.txt", A, n, m);
		writeToFile("b.txt", b, n, 1);
	}
	free(A);
	free(b);
	free(x);
}
