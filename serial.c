#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <mkl_types.h>
#include <mkl_cblas.h>
#include <mkl.h>
#include <string.h>
#include <time.h>
#include "func/matrixIO.h"
int main(int argc, char * argv[]){
	MPI_Init(&argc, &argv);
	double elapsedTime;

	int n = 10000;		// Number of rows in matrix A
	int m = 10000;				// Number of columns in matrix A
		// Read in matrix on root process and distribute	


	double * A;
	double * bx;
	A = malloc(n*m*sizeof(double));
	readArray("A.txt", A, n, m);
	bx = malloc(n*sizeof(double));
	readArray("b.txt", bx, n, 1);
	int info;
	double start = MPI_Wtime();
	int *ipvt = malloc(n*sizeof(int));
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, A, n, ipvt);
	double end = MPI_Wtime();
	elapsedTime = end - start;
	printf("\nElapsed time = %lf\n", elapsedTime);
	/*LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', n, 1, A, n, ipvt, bx, 1);
	double sum = 0.0;
	for (int i=0; i<m; i++){
		sum += bx[i];
	}
	sum /= m;
	printf("\n\naverage value in solution x = %lf\n", sum);*/
	free(A);
	free(bx);
	MPI_Finalize();
}
	

	
