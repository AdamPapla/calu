/**
* @file	luFunc.c 
* @author Adam Graham
* @brief Functions to perform in place LU factorisation both in parallel and serial
* @date 2023-8-2
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include "matrix.h"
#include "calcBlock.h"
#include "matrixIO.h"
#include "luFunc.h"
/**
 * @brief	Performs an unpivoted parallel LU factorisation on a panel of the matrix to complete TSLU.
 * @param[in]	A			Array containing the processes local chunk of the overall matrix to be factorised.
 * @param[in]	blockNum	The number of blocks of b columns and rows already processed in the CALU recursion.
 * 							Necessary for controlling column offsets of local chunk.
 * @param[in] 	locBlockNum	Number of blocks of b columns and rows already processed on the currently active processors
 * 							local chunk of A. Necessary for controlling row offset in local array. 
 * @param[in]	block 		The number of rows remaining unprocessed in the current local chunk of the matrix. Dictates the
 * 							number of rows locally in the panel to be factorised via TSLU.
 * @param[in]	m			The number of columns in A.
 * @param[in] 	b			The block size of the CALU recursion.
 * @param[in] 	rank		Rank of processor in communicator.
 * @param[in] 	nprocs		Number of processes in communicator.
 * @param[in]	comm		Communicator of all active processes.
 * @param[out] 	A   		Updated local chunk of the matrix, with the panel at position specified by 
 * 							blockNum and locBlockNum factorised in place into L and U
 */
void unpivotedLUA(double *A, int blockNum, int locBlockNum, int block, int m, int b, int rank, int nprocs, MPI_Comm comm){
	// On root, account for row offset from rows already processed locally. On all processes account for global column offset
	int startIndex = (rank==0) ? locBlockNum * b * m + blockNum*b : blockNum * b;
	int start; int end;
	// Must broadcast U11 from root to all processes for LU factorisation. See report
	double *U11 = malloc(b*b* sizeof(double));

	if (rank == 0){
		for (int k=0; k<b; k++){
			for (int i= k+1; i<block; i++){
				// Divide row below diagonal by diagonal element
				A[startIndex + i*m + k] /= A[startIndex + k*m + k];
				for (int j=k+1; j<b; j++){
					// Update trailing matrix
					A[startIndex + i*m + j] -= A[startIndex + i*m + k] * A[startIndex + k*m + j];
				}
			}
		}
		// Construct U from A
		for (int i=0; i<b; i++)
			for (int j=0; j<b; j++)
				U11[i*b + j] = (j>= i) ? A[startIndex + i*m + j] : 0;
	}
	// Broadcast U
	MPI_Bcast(U11, b*b, MPI_DOUBLE, 0, comm);
	// Use U to update L on other processes
	if (rank != 0){
		for (int k=0; k<b; k++){
			for (int i=0; i<block; i++){
				A[startIndex + i*m + k] /= U11[k*b + k];
				for (int j=k+1; j<b; j++){
					A[startIndex + i*m + j] -= A[startIndex + i*m + k] * U11[k*b+j];
				}
			}
		}
	}	
	free(U11);
}

/**
 * @brief	Performs local unpivoted LU on the last part of the matrix A. If m divides b, this is a block of 
 * 			b columns. Else it is a block of m%b columns.
 * @param[in]	A			Array containing the processes local chunk of the overall matrix to be factorised.
 * @param[in]	blockNum	The number of blocks of b columns and rows already processed in the CALU recursion.
 * 							Necessary for controlling column offsets of local chunk.
 * @param[in] 	locBlockNum	Number of blocks of b columns and rows already processed on the currently active processors
 * 							local chunk of A. Necessary for controlling row offset in local array. 
 * @param[in]	block 		The number of rows remaining unprocessed in the current local chunk of the matrix. Dictates the
 * 							number of rows locally in the panel to be factorised via TSLU.
 * @param[in]	m			The number of columns in A.
 * @param[in] 	b			The block size of the CALU recursion.
 * @param[out] 	A   		Updated local chunk of the matrix, with the final panel of the matrix factorised in place
 * 			 				into L and U, completing CALU.
 */

void finalLUA(double *A, int blockNum, int locBlockNum, int block, int m, int b, int rem){
	int startIndex = locBlockNum*b * m + blockNum*b;
	for (int k=0; k<rem; k++){
		for (int i=k+1; i<block; i++){
			A[startIndex + i*m + k] /= A[startIndex + k*m + k];	// Divide by diagonal elements
			for (int j=k+1; j<rem; j++){
				A[startIndex + i*m + j] -= A[startIndex + i*m + k] * A[startIndex + k*m + j];
			}
		}
	}
}


