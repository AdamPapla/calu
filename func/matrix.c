/**
* @file matrix.c
* @author Adam Graham
* @brief Create functions to allocate memory for matrices and print them
* @date 2023-3-5
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "matrix.h"
#include "calcBlock.h"
#include "matrixIO.h"

/**
 * @brief 	Copies the next panel to be factorised with TSLU into a separate
 * 			array for performing tournament pivoting. 
 * @param[in]	A		Array containing the processes local chunk of the overall matrix to be factorised.
 * @param[in]	blockNum	The number of blocks of b columns and rows already processed in the CALU recursion.
 * 							Necessary for controlling column offsets of local chunk.
 * @param[in] 	locBlockNum	Number of blocks of b columns and rows already processed on the currently active processors
 * 							local chunk A. Necessary for controlling row offset in local array. 
 * @param[in]	block 		The number of rows remaining unprocessed in the current local chunk A. Dictates the
 * 							number of rows locally in the panel to be factorised via TSLU.
 * @param[in]	m			The number of columns in A.
 * @param[in] 	b			The block size of the CALU recursion.
 * @param[in] 	rank		Rank of processor in communicator
 * @param[in]	rem			Number of columns in panel. Always b except in the last step of the recursion, where it 
 * 							will have value m%b to account for the last columns which are not enough for a full b chunk.
 * @param[out]	locPanel	Array containing this processors chunk of the panel to be factorised using TSLU.
 */
void fillPanel(double *locA, double *locPanel, int blockNum, int locBlockNum, int block, int m, int b, int rank, int rem){
	// If root, must account for rows already processed. Account for columns already processed on all processes.
	int startIndex = (rank == 0) ? locBlockNum*b * m + blockNum*b : blockNum*b;

	for (int i=0; i<block; i++){
		for (int j=0; j<rem; j++){
			locPanel[i*rem + j] = locA[startIndex + i*m + j];
		}
	}
}

/**
 * @brief 	Extracts a copy of L_11 from the in place position in the A matrix into an array containing it 
 * 			in its full lower triangular form. This is used for the A12 update.
 * @param[in]	A			Array containing the processes local chunk of the overall matrix to be factorised.
 * @param[in]	blockNum	The number of blocks of b columns and rows already processed in the CALU recursion.
 * 							Necessary for controlling column offsets of local chunk.
 * @param[in] 	locBlockNum	Number of blocks of b columns and rows already processed on the currently active processors
 * 							local chunk A. Necessary for controlling row offset in local array. 
 * @param[in]	m			The number of columns in A.
 * @param[in] 	b			The block size of the CALU recursion.
 * 							will have value m%b to account for the last columns which are not enough for a full b chunk.
 * @param[in]	rem			The size of the L array. Will be b in all but the final step where m>n
 * @param[out]	L			Array containing the L11 matrix in full lower triangular form.
 */

void extractL(double *A, int blockNum, int locBlockNum, int m, int b, int rem, double *L){
	// Start index of L11
	int startIndex = locBlockNum*b * m + blockNum*b;
	for (int i=0; i<rem; i++){
		for (int j=0; j<i; j++){
			L[i*rem + j] = A[startIndex + i*m + j];
		}
		// Set upper triangular entries to 0
		for (int j=i+1; j<rem; j++){
			L[i*rem+j] = 0;
		}
		// Set diagonal entries to 1
		L[i*rem + i] = 1;
	}
}

/**
 * @brief	Gathers the full matrix onto the root process.
 * @param[in]	locA		Local chunk of the matrix A
 * @param[in]	n			The number of rows in A.
 * @param[in]	m			The number of columns in A.
 * @param[in] 	b			The block size of the CALU recursion.
 * @param[in]	block 		Number of rows in local chunk
 * @param[in] 	rank		Rank of processor in communicator.
 * @param[in] 	nprocs		Number of processes in communicator.
 * @param[in]	comm		Communicator of all active processes.
 * @param[out] 	A   		Full matrix A, gathered from all processes in communicator.
 */
void gatherMatrix(double *locA, int n, int m, int b, int block, int rank, int nprocs, MPI_Comm comm, double *A){
	int s, e;
	if (rank == 0){
		int currentIndex = 0;		// Keeps track of index that next process will begin copyign to
		memcpy(A, locA, block*m*sizeof(double));	// Copy local array to the start of global array
		currentIndex += block*m;	// Update current index to receive from first process
		for (int i=1; i<nprocs; i++){
			calcBlock(n, nprocs, b, i, &s, &e);		// Calculate number of rows to receive
			int sourceBlock = e-s;
			MPI_Recv(A+currentIndex, sourceBlock*m, MPI_DOUBLE, i, 0, comm, MPI_STATUS_IGNORE);
			currentIndex += sourceBlock*m;		// Update current index to receive from next process
		}
	}
	// Send from all nonroot processes
	else{
		MPI_Send(locA, block*m, MPI_DOUBLE, 0, 0, comm);
	}
}
	
/**
 * @brief	Gathers the full pivot array onto the root process.
 * @param[in]	ipvt		Local chunk of the array of pivots
 * @param[in]	n			The number of rows in pivot array.
 * @param[in] 	b			The block size of the CALU recursion.
 * @param[in]	block 		The number of rows in the local
 * @param[in] 	rank		Rank of processor in communicator.
 * @param[in] 	nprocs		Number of processes in communicator.
 * @param[in]	comm		Communicator of all active processes.
 * @param[out] 	ipvtFull	Full pivot array, gathered from all processes
 */
void gatherPivot(int *ipvt, int n, int b, int block, int rank, int nprocs, MPI_Comm comm, int *ipvtFull){
	int s, e;
	if (rank == 0){
		int currentIndex = 0;
		memcpy(ipvtFull, ipvt, block*sizeof(int));
		currentIndex += block;
		for (int i=1; i<nprocs; i++){
			calcBlock(n, nprocs, b, i, &s, &e);
			int sourceBlock = e-s;
			MPI_Recv(ipvtFull+currentIndex, sourceBlock, MPI_INT, i, 0, comm, MPI_STATUS_IGNORE);
			currentIndex += sourceBlock;
		}
	}
	else{
		MPI_Send(ipvt, block, MPI_INT, 0, 0, comm);
	}
}
		
/**
 * @brief	Back substitution to solve a system of linear equations given the answer vector
 * 			and the LU factorisation of A
 * @param[in]	A		Array containing the full LU factorisation of A (stored in place)
 * @param[in]	n		Number of rows in A
 * @param[in]	m		Number of columns in A
 * @param[in]	b		Array containing answer vector (pivoted along with LU)
 * @param[out]	x		Solution vector 
 */
void triangSub(double *A, int n, int m, double *b, double *x){
	// Forward substitiution: solve Lc = b
	double *c = malloc(n*sizeof(double));
	for (int i=0; i<n; i++){
		c[i] = b[i];
		for (int j=0; j<i; j++){
			c[i] -= A[i*m+j]*c[j];
		}
	}
	// Now solve Ux=c
	for (int i=n-1; i>=0; i--){
		x[i] = c[i];
		for (int j=i+1; j<m; j++){
			x[i] -= A[i*m + j] * x[j];
		}
		x[i] /= A[i*m + i];
	}
	free(c);
}


