#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <mkl_types.h>
#include <mkl_cblas.h>
#include <mkl.h>
#include <string.h>
#include "calcBlock.h"
#include "matrix.h"
#include "matrixIO.h"
#include "pivot.h"
#include "luFunc.h"
/**
 * @brief 	Function to perform parallel TSLU. Allocates everything but locA internally. 
 * @param[in] 	locA 		Pointer to local chunk of matrix
 * @param[in] 	n 			Number of rows in overall matrix that remains to be processed
 * @param[in] 	b 			Number of columns of panel
 * @param[in]	locBlockNum	Number of blocks of b rows already processed on root
 * @param[in]	comm		Communicator with all active processes
 * @param[out] 	indices 	Pointer to array holding the b global indices that should be 
 * 						   	pivoted to the top. Note these indices start at the 
 * 						   	first row of the matrix that remains to be processed.
 * 						   	Fully processed rows are excluded from this indexing
 */
int tournPivot(double *locA, int n, int b, int locBlockNum, double **indices, MPI_Comm comm){
	int nprocs;
	int rank;
	int myInfo=0;
	MPI_Comm_size(comm, &nprocs);
	MPI_Comm_rank(comm, &rank);
	int activeProcs = nprocs; 
	

	// Calculate start and end of this processors chunk of matrix
	

	int s = 0;
	int e = 0;
	// Must calculate block as if it's the first panel of the process then shift backward
	calcBlock(n+locBlockNum*b, nprocs, b, rank, &s, &e);	
	if (rank != 0){
		s -= locBlockNum * b;
	}
	e -= locBlockNum * b;
	int block = e-s;
	
	
	// Allocate to store pivots and initialise indices array
	

	*indices = malloc(block*sizeof(double));
	for (int i=0; i<block; i++)
		(*indices)[i] = i+s; 	// Indices of elements in this procs locA matrix
	int *pivot = malloc(b*sizeof(int));
	double *PLU = malloc(block*b*sizeof(double));
	

	// Variables to control binomial tree
	

	int oddProcess = 0;	// Will be set to true if process is odd
	int isActive = 0;
	int wasActive = 0;
	int iter = 1;


	// Enter binomial tree
	

	while (activeProcs > 1){
		
		isActive = pow(2,iter);	// Active if receiving
		wasActive = pow(2,iter-1);	// Was active if sending then terminates


		// Perform local LU to find pivots and apply them (unless process is odd- then this has already been done)


		if (oddProcess == 0){
			for (int i=0; i<block*b; i++)
				PLU[i] = locA[i];	// Pass in copy of A to avoid overwriting A with L and U
			int ret = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, block, b, PLU, b, pivot); 	// Perform local LU
			for (int i=0; i<b; i++)
				pivot[i] -= 1;	// Account for pivot indexing starting at 1 in return
			
			pivotMatrix(locA, pivot, b, b);		// Pivot A and indices arrays
			pivotMatrix((*indices), pivot, b, 1);
		}


		// If process is receiving allocate and store source processes best b rows under local best b rows


		if (rank%isActive == 0){
			int sourceRank = rank+pow(2,iter-1);
		
			// Reallocate locA and indices to hold 2b*b and 2b, ready to recieve
			block = 2*b;
			if (iter == 1){
				locA = realloc(locA, (block*b)*sizeof(double));
				PLU = realloc(PLU, block*b*sizeof(double));
				*indices = realloc(*indices, block*sizeof(double)); 
			}
			
			// If process is odd, perform necessary updates. No recieves
			if (sourceRank >= nprocs){
				oddProcess = 1;
				activeProcs = activeProcs/2 + 1;	// We know activeProcs is odd so can do this
			}
			
			// Else if process is even, recieve from source process and update
			else if (sourceRank < nprocs) {
				MPI_Recv(locA + b*b, b*b, MPI_DOUBLE, sourceRank, 0, comm, MPI_STATUS_IGNORE);
				MPI_Recv((*indices)+b, b, MPI_DOUBLE, sourceRank, 0, comm, MPI_STATUS_IGNORE);
				activeProcs = (activeProcs%2 == 0) ? activeProcs/2 : activeProcs/2 + 1;
			}

		}


		// If process is even, send then terminate
		

		else if (rank%wasActive == 0){
			int targetRank = rank-pow(2,iter-1); // Rank this proc will send to
			MPI_Send(locA, b*b, MPI_DOUBLE, targetRank, 0, comm);
			MPI_Send((*indices), b, MPI_DOUBLE, targetRank, 0, comm);
			activeProcs = 1;
		}
		iter++;
		
	}


	// Final step of TSLU - pivot result of binomial tree
	

	if (rank==0){
		for (int i=0; i<block*b; i++)
			PLU[i] = locA[i];
		LAPACKE_dgetrf(LAPACK_ROW_MAJOR, block, b, PLU, b, pivot); 	// Perform local LU
		for (int i=0; i<b; i++)
			pivot[i] -= 1;	// Account for pivot indexing starting at 1 in return
		pivotMatrix(locA, pivot, b, b);
		pivotMatrix((*indices), pivot, b, 1);
	}

	free(pivot);
	free(PLU);
	return 0;
}


	
