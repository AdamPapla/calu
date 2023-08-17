/**
* @file pivot.c
* @author Adam Graham
* @brief Create functions to pivot matrix
* @date 2023-15-7
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "matrix.h"
#include "calcBlock.h"
#include "pivot.h"

/**
 * @brief	Pivot a matrix locally, with pivot array being applied iteratively from 
 * 			bottom to top
 * @param[in]	A		Matrix to be pivoted
 * @param[in] 	pivot	Array holding the pivot indices, which should be applied
 * 						iteratively in reverse to A. For detail, see report.
 * @param[in]	b		Length of pivot array.
 * @param[in]	m		Number of columns in m
 * @param[out]	A		Pivoted matrix A
 */
void pivotMatrix(double *A, int *pivot, int b, int m){
	double *temp = malloc(m*sizeof(double));
	for (int i=b-1; i>=0; i--){
		for (int j=0; j<m; j++){
			temp[j] = A[i*m + j]; 	// Store ith row in temp
			A[i*m + j] = A[pivot[i]*m + j]; // Store A[pivot[i],j] in A[i,j]
			A[pivot[i]*m + j] = temp[j];	// Complete swap
		}
	}
	free(temp);
}
/**
 * @brief	Pivot a matrix that is spread out over processes in a communicator.
 * 			Handles the cases relevant in recursion where matrix to be pivoted
 * 			essentially starts at a specified number of rows into A. 
 * 			By design, CALU never has to pivot from one arbitrary process to another.
 * 			Pivoting always occurs between root and the other processes, as the 
 * 			top b rows being pivoted always belong to the root process.
 * 			Locates on which process the row to be pivoted to root resides and 
 * 			handles all MPI communication necessary to perform the pivot.
 * 			Uses memcpy for local pivoting.
 * @param[in]	A			Local chunk of matrix to be pivoted
 * @param[in] 	indices		Array returned from ::tslu function specifying which indices
 * 							must be pivoted to the top. Again, these should be performed
 * 							iteratively in reverse order. 
 * @param[in]	n			Number of rows remaining to be processed in the global matrix.
 * @param[in]	m			Number of columns in global matrix.
 * @param[in]	b			Block size for CALU.
 * @param[in]	locBlockNum	Number of blocks already processed on the root processor.
 * @param[in]	comm		MPI communicator with currently active processes.
 * @param[out]	A			Local chunk of the pivoted matrix 
 */
void pivotFullRow(double *A, double *indices, int n, int m, int b, int locBlockNum, MPI_Comm comm){
	int nprocs;
	int rank;
	MPI_Comm_size(comm, &nprocs);
	MPI_Comm_rank(comm, &rank);

	// For rank 0
	
	if (rank == 0){

		// Calculate the starting index of all processors

		int *starts = malloc((nprocs+1)*sizeof(int));	// +1 so we can add the end index
		for (int i=0; i<nprocs; i++){
			int end; 
			calcBlock(n+locBlockNum*b, nprocs, b, i, &starts[i], &end);
		}
		// Account for locBlockNum
		for (int i=0; i<nprocs; i++){
			if (i!= 0)
				starts[i] -= locBlockNum*b;
		}
		starts[nprocs] = n;
		double *temp = malloc(m * sizeof(double));

		// Decision tree to decide whether to memcopy this row or recieve and send to other process

		// Calculate process which owns row to be pivoted
		for (int i=b-1; i>= 0; i--){	// Loop through indices array
			int destRank = -1;
			for (int j=0; j<nprocs; j++){		// loop through array of starts
				if (indices[i] < starts[j+1]){
					destRank = j;
					break;
				}			
		 	}

			// Send/recv or memcpy
			if (destRank != 0){
				MPI_Send(A + locBlockNum*b*m + i*m, m, MPI_DOUBLE, destRank, 0, comm);
				MPI_Recv(A + locBlockNum*b*m + i*m, m, MPI_DOUBLE, destRank, 0, comm, MPI_STATUS_IGNORE);
			}
			else {
				int startIndex = locBlockNum*m*b + indices[i]*m;
				// Use memcpy to swap rows locally
				memcpy(temp, A + locBlockNum*b*m + i*m, m*sizeof(double));	// Copy ith row of A to temp
				memcpy(A + locBlockNum*b*m + i*m, A + startIndex, m*sizeof(double));	// Copy indices[i]th row of A to ith row of A
				memcpy(A + startIndex, temp, m*sizeof(double)); 	// Copy temp to indices[i]th row of A
			}
		}
		free(temp);
	}

	// For all other ranks 
	
	else {
		int start = 0;
		int end = 0;
		// Calculate the start and end row index of this processors chunk
		calcBlock(n + locBlockNum*b, nprocs, b, rank, &start, &end);
		// Necessary offset to account for rows already processed- details in report
		start -= locBlockNum*b;
		end -= locBlockNum*b;
		double *temp = malloc(m*sizeof(double));
		for (int i=b-1; i>=0; i--){
			if (indices[i] >= start && indices[i] < end){
				int index = indices[i] - start;
				// Copy row to be sent to temp then recieve + send temp
				memcpy(temp, A + index*m, m*sizeof(double));
				MPI_Recv(A + index*m, m, MPI_DOUBLE, 0, 0, comm, MPI_STATUS_IGNORE);
				MPI_Send(temp, m, MPI_DOUBLE, 0, 0, comm);
			}
			
		}
		free(temp);
	}
}
