/**
* @file matrix.c
* @author Adam Graham
* @brief Create functions to allocate memory for matrices and print them
* @date 2023-3-5
*/


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "func/calcBlock.h"
#include "func/matrixIO.h"
#include "func/matrix.h"
#include "func/randomGen.h"
#include "func/calu.h"

int main(int argc, char * argv[]){
	MPI_Init(&argc, &argv);
	double startTime,endTime;	// Used for timings
	double elapsedTime;
	int n = 10000;		// Number of rows in matrix A
	int m = 10000;				// Number of columns in matrix A
	int b = 200;				// Panel width


	int nprocs;
	int rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int s = 0;
	int e = 0;
	calcBlock(n, nprocs, b, rank, &s, &e);		// Start and end of this processors chunk of matrix
	int block = e-s;
	double *locA = malloc(block*m*sizeof(double));		
	double *locb = malloc(block*sizeof(double));	
	// Read in matrix on root process and distribute	
	read_chunk_from_file("A.bin", locA, m, s, e, nprocs, rank, MPI_COMM_WORLD); 
	read_chunk_from_file("b.bin", locb, 1, s, e, nprocs, rank, MPI_COMM_WORLD); 


	MPI_Barrier(MPI_COMM_WORLD);
	startTime = MPI_Wtime();

	calu(locA, locb, n, m, b);

	MPI_Barrier(MPI_COMM_WORLD);
	endTime = MPI_Wtime();

	elapsedTime = endTime - startTime;
	if (rank == 0)
		printf("\nTime for CALU: %lf seconds\n", elapsedTime);
	double *A;
	double *bx;
	if (rank == 0){
		A = malloc(m*n*sizeof(double));
		bx = malloc(n*sizeof(double));
	}
	MPI_Barrier(MPI_COMM_WORLD);
	gatherMatrix(locA, n, m, b, block, rank, nprocs, MPI_COMM_WORLD, A);
	gatherMatrix(locb, n, 1, b, block, rank, nprocs, MPI_COMM_WORLD, bx);
	MPI_Barrier(MPI_COMM_WORLD);	

	if (rank == 0){
			
		double *x = malloc(n*sizeof(double));
		triangSub(A, n, m, bx, x);
		double sum = 0.0;
		for (int i=0; i<m; i++){
			sum += x[i];
		}
		sum /= m;
		printf("\n\naverage value in solution x = %lf\n", sum);
		writeToFile("x.txt", x, m, 1);
		free(x);
		free(A);
		free(bx);
	}

	// Freeing allocated memory
	free(locA);
	free(locb);
	MPI_Finalize();

}	

	
