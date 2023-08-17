/**
* @file oldFunc.c
* @author Adam Graham
* @brief Functions in previous implementations which were made redundant by
* 		 optimisations referenced in the report
* @date 2023-8-15
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "matrix.h"
#include "decomp1d.h"
#include "matrixIO.h"

void distributeMatrix(double *A, double *locA, int n, int m, int b, int myid, int nprocs, int s, int e){
	if (myid == 0){		// Loop through processes on root, calculate chunk and send
		for (int i=1; i<nprocs; i++){
			int sDest, eDest;
			calcBlock(n, nprocs, b, i, &sDest, &eDest);
			int blockDest = eDest - sDest;
			MPI_Send(A+sDest*m, blockDest*m, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
		memcpy(locA, A, (e-s)*m*sizeof(double));
	}
	else{
		int block = e-s;
		MPI_Recv(locA, block*m, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}

void pivotTransform(int **pivot, int n, int b){
	printf("\nEntering pivot transform\n");
	int *newPivot = malloc(n*sizeof(int));
	
	for (int i=0; i<n; i++){
		newPivot[i] = i;
	}
	for (int i=b-1; i>=0; i--){
		int temp = newPivot[(*pivot)[i]];
		newPivot[(*pivot)[i]] = newPivot[i];
		newPivot[i] = temp;
	}
	for (int i=0; i<b; i++)
		(*pivot)[i] = newPivot[i];	
	free(newPivot);
	
}


void extractU12(double *locA, double *blockU, int m, int b, int blockNum, int locBlockNum){
	int rowLength = m-(blockNum + 1)*b;
	int startingRowIndex = locBlockNum*b * m;
	int startingColIndex = (blockNum + 1)*b;	// +1 to account for the LU block not included at the start
	for (int i=0; i<b; i++){
		for (int j=0; j<rowLength; j++){
			blockU[i*rowLength + j] = locA[startingRowIndex + i*m + startingColIndex + j];
		}
	}
}

void insertU12(double *locA, double *blockU, int m, int b, int blockNum, int locBlockNum){
	int rowLength = m-(blockNum + 1)*b;
	int startingRowIndex = locBlockNum*b * m;
	int startingColIndex = (blockNum + 1)*b;	// +1 to account for the LU block not included at the start
	for (int i=0; i<b; i++){
		for (int j=0; j<rowLength; j++){
			locA[startingRowIndex + i*m + startingColIndex + j] = blockU[i*rowLength + j];
		}
	}
}

void extractL21(double *locA, double *L21, int b, int m, int block, int blockNum, int locBlockNum, int rank){
	int startIndex = (rank==0) ? (locBlockNum+1) * b * m + blockNum*b : blockNum * b;
	int colLength = (rank == 0) ? block - b : block;
	for (int i=0; i<colLength; i++){
		for (int j=0; j<b; j++){
			L21[i*b + j] = locA[startIndex + i*m + j];
		}
	}
}

void extractA22(double *locA, double *A22, int b, int m, int block, int blockNum, int locBlockNum, int rank){
	int startIndex = (rank == 0) ? (locBlockNum+1) * b * m + (blockNum+1) * b : (blockNum+1) * b;
	int rowLength = m - (blockNum+1)*b;
	int colLength = (rank == 0) ? block - b : block;
	for (int i=0; i < colLength; i++){
		for (int j=0; j<rowLength; j++){
			A22[i*rowLength + j] = locA[startIndex + i*m + j];
		}
	}
}

void inputA22(double *locA, double *A22, int b, int m, int block, int blockNum, int locBlockNum, int rank){
	int startIndex = (rank == 0) ? (locBlockNum+1) * b * m + (blockNum+1) * b : (blockNum+1) * b;
	int rowLength = m - (blockNum+1)*b;
	int colLength = (rank == 0) ? block - b : block;
	for (int i=0; i < colLength; i++){
		for (int j=0; j<rowLength; j++){
			locA[startIndex + i*m + j] = A22[i*rowLength + j]; 
		}
	}
}



