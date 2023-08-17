/**
* @file calu.c
* @author Adam Graham
* @brief 	A communication avoiding LU factorisation algorithm which performs an in place 
* 			factorisation in parallel through MPI. CBLAS routines are used for local matrix-matrix
* 			and matrix-vector operations. Works for any arbitrary m x n matrix on any number of processors
* @date 2023-8-15
*/

// Libraries used throughout CALU
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

// User defined function files
#include "calcBlock.h"
#include "matrixIO.h"
#include "matrix.h"
#include "pivot.h"
#include "luFunc.h"
#include "tslu.h"
#include "calu.h"

/**
 * @brief	Function for performing parallel Communication Avoiding LU Factorisation. Performs 
 * 			factorisation in place, with low communication and low memory footprint. Works
 * 			on any arbitrary m x n matrix and mantains stability in all cases tested.
 * 			CALU is a block factorisation algorithm, proceeding recursively on the matrix in blocks
 * 			of b rows and columns. 
 * 			Note there may be some special case matrices which have stability concerns, as is the 
 * 			case with gaussian elimination.
 * 			Note also that MPI should be initialised before using this function through calling 
 * 			MPI_Init in the main file.
 * @param[in]	locA		Local chunk of the matrix A to be factorised. Overall matrix A should be 
 * 							read in using parallel MPI IO/distributed from root before passing to this 
 * 							function
 * @param[in]	locb		Effectively keeps track of pivot information. Can be passed in in multiple 
 * 							forms.
 * 							For solving a system of linear equations, pass in the rhs vector b (with
 * 							Ax=b being the system).
 * 							For more general purposes, pass in an array containing the row indices of A
 * 							i.e. [0:n-1] where n is the number of rows in A. This will hold information
 * 							on the pivots applied to A throughout the algorithm
 * @param[in] 	n			Number of rows in matrix A (and vector b)
 * @param[in]	m			Number of columns in matrix A (and rows in solution vector x for system of
 * 							linear equations
 * @param[in]	b			Block size of the algorithm
 * @param[out]	locA		Local chunk of global matrix A holding U in the upper triangular (including diagonal) 
 * 							section and L in the lower triangular section (when constructing L, diagonal should 
 * 							be set =1) such that LU = PA i.e. L*U is equal to A, where the rows of A have been 
 * 							permuted. Note locA holds only some number of rows of the overall matrix A.
 * @param[out]	locb		Local chunk of vector b, which contains the pivot information. 
 * 							If passed in as an rhs vector b of a system of linear equations, it will be 
 * 							pivoted in the same fashion as A. In this form, one can directly use A and b 
 * 							to solve for x by back substitution with no additional pivoting.
 * 							If passed in as an array containg [0:n-1], it will output the new order of the row
 * 							indices of A. For example, if A=[[1,3,5], [4,6,8], [2,1,3]] on input, and 
 * 							the output L and U are s.t L*U=[[4,6,8], [2,1,3], [1,3,5]], we will have
 * 							b=[1,2,0], the new order of the rows.  
 */


void calu(double *locA, double *locb, int n, int m, int b){
	int nprocs;		// Number of processes in MPI_COMM_WORLD
	int rank;		// Rank of current process
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
		printf("\n\nComputing CALU of %d x %d matrix with block size %d on %d processors\n\n", n, m, b, nprocs);


	// Calculating local block size and initialising variables to control recursion


	int s = 0;
	int e = 0;
	calcBlock(n, nprocs, b, rank, &s, &e);		// Start and end of this processors chunk of matrix
	int block = e-s;
	int blockNumber = 0;
	int locBlockNum = 0;



	// Setting up communicator and variables to determine when it must be redefined
	

	MPI_Group oldProcs, newProcs;
	MPI_Comm comm, newComm;
	int activeProcs = nprocs;
	int rootEnd, junk;
	comm = MPI_COMM_WORLD;
	calcBlock(n, nprocs, b, 0, &junk, &rootEnd);		// Calculating end of initial root


	// Recursion loop
	

	double *indices;	// Pivot information
	double *locPanel;	// Holds panel for TSLU
	MPI_Datatype sendU12;

	while ((blockNumber+1)*b < n && (blockNumber+1)*b < m){



		// TSLU of current panel in recursion
		


		// Extract panel
	 	locPanel = malloc(block*b*sizeof(double));	// Allocate memory for panel
		fillPanel(locA, locPanel, blockNumber, locBlockNum, block, m, b, rank, b);	// Fill the nth panel of the processor

		// Perform TSLU on panel 
		tournPivot(locPanel, n-blockNumber*b, b, locBlockNum, &indices, comm);
		// Broadcast indices to all processes 
		if (rank != 0)
			indices = malloc(b*sizeof(double));
		MPI_Bcast(indices, b, MPI_DOUBLE, 0, comm);

		// Perform pivots on full matrix
		pivotFullRow(locA, indices, n-blockNumber*b, m, b, locBlockNum, comm);
		pivotFullRow(locb, indices, n-blockNumber*b, 1, b, locBlockNum, comm);

		// Perform unpivoted parallel LU factorisation on now pivoted panel (performed in place)
		unpivotedLUA(locA, blockNumber, locBlockNum, block, m, b, rank, nprocs, comm);



		// Update of U12 and trailing matrix A22
	


		double *L;
		double *U12;
		double *L21;
		double *A22;
		int rowLength = m - (blockNumber+1) * b;	// Number of cols in U12 and A22

		// Define datatype for broadcasting noncontiguous U12 matrix 
		MPI_Type_vector(b, rowLength, m, MPI_DOUBLE, &sendU12);
		MPI_Type_commit(&sendU12);
		
		// On root, update U, broadcast it and update its chunk of the trailing matrix
		if (rank == 0){
			// Indices of the start of relevant blocks
			int startIndexA12 = locBlockNum*b * m + (blockNumber+1)*b;
			int startIndexA22 = (locBlockNum+1)*b * m + (blockNumber+1)*b;
			int startIndexL21 = (locBlockNum+1)*b * m + blockNumber*b;
			int colLength = block - b;
			
			// Extract L11
			L = malloc(b*b*sizeof(double));
			extractL(locA, blockNumber, locBlockNum, m, b, b, L);

			// Update U using dtrsm according to U12 = inverse(L11) * A12
			cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, b, rowLength, 1.0, L, b, locA+startIndexA12, m);

			// Broadcast U12 to all processes using defined MPI vector type
			MPI_Bcast(locA + startIndexA12, 1, sendU12, 0, comm);

			// Perform update on roots chunk of trailing matrix according to A22 = A22 - L21 * U12
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, colLength, rowLength, b, -1, locA + startIndexL21, m, locA + startIndexA12, m, 1, locA+startIndexA22, m);

			free(L);
		}

		// On nonroot active processes, receive U12 from root and update trailing matrix
		else {
			// Local start index of A22 and L21. Everything is acted on in place in row major order, so these blocks of A are naturally discontiguous
			int startIndexA22 = (blockNumber+1) * b;
			int startIndexL21 = blockNumber * b;
			
			// Recieve U from top process and perform updates
			U12 = malloc(b * rowLength * sizeof(double)); 
			MPI_Bcast(U12, rowLength*b, MPI_DOUBLE, 0, comm);

			// Perform update on trailing matrix according to A22 = A22 - L21 * U12
			cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, block, rowLength, b, -1, locA + startIndexL21, m, U12, rowLength, 1, locA+startIndexA22, m);

			free(U12);
		}
		

		
		// Recursion variable and communicator updates




		// Recursion variable updates
		if (rank == 0){
			block -= b;		// On root, the remaining rows to be processed decreases by b
		}
		blockNumber++;		// Another block has been processed, both locally and globally. 
		locBlockNum++;

		// Redefining communicator if at end of roots block of rows
		MPI_Barrier(comm);
		if (rootEnd == locBlockNum*b){
			MPI_Comm_group(comm, &oldProcs);		// Create group from current communicator
			int exclRank[1] = {0};
			MPI_Group_excl(oldProcs, 1, exclRank, &newProcs);		// Exclude root process from communicator
			MPI_Comm_create(comm, newProcs, &newComm);		// Create new comm from this group
			comm = newComm;		// Update comm and oldProcs
			oldProcs = newProcs;
			if (rank == 0){		// Terminate on now excluded rank
				activeProcs = 1;
				rank = -1;
				break;
			}
			MPI_Comm_rank(comm, &rank);	// Update ranks to those on new communicator
			activeProcs--;
			calcBlock(n-blockNumber*b, activeProcs, b, 0, &junk, &rootEnd);		// Calculating end of new root's block of rows
			locBlockNum = 0;	// Reset local number of blocks processed
			block = e-s;
		}
		MPI_Type_free(&sendU12);
		free(indices);
		free(locPanel);
	}


	// Performing final step of recursion. Final block will have <=b rows or columns
	

	if (rank == 0){
		// Calculating how many rows/columns are remaining (depending on which is terminating factor)
		int rem = 0;
		if (n >= m)
			rem = (m%b == 0) ? b : m%b;		// Number of columns of final block 
		else
		 	rem = (n%b == 0) ? b : n%b;		// Number of rows of final block

		// Find pivots locally for this block of rows using dgetrf (local LU factorisation)

		locPanel = malloc(block * rem * sizeof(double));
		fillPanel(locA, locPanel, blockNumber, locBlockNum, block, m, b, rank, rem);	// Fill the nth panel of the processor
		int * ipvt = malloc(rem * sizeof(int));
		LAPACKE_dgetrf(LAPACK_ROW_MAJOR, block, rem, locPanel, rem, ipvt); 	// Perform local LU
		for (int i=0; i<rem; i++)
			ipvt[i] -= 1;

		// Pivot matrix with these pivots

		pivotMatrix(locA + locBlockNum*b*m, ipvt, rem, m);
		pivotMatrix(locb + locBlockNum*b, ipvt, rem, 1);

		// Perform local unpivoted LU

		finalLUA(locA, blockNumber, locBlockNum, block, m, b, rem);

		
		// If matrix has more columns than rows, must perform a final U update
		if (n < m){
			int startIndexA12 = locBlockNum*b * m + (blockNumber*b + rem);
			int rowLength = m - (blockNumber*b + rem);
			// Allocate for and extract L
			double *L = malloc(rem*rem*sizeof(double));
			extractL(locA, blockNumber, locBlockNum, m, b, rem, L);

			// Update U using dtrsm according to U12 = inverse(L11) * A12
			cblas_dtrsm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, rem, rowLength, 1.0, L, rem, locA+startIndexA12, m);
			free(L);
		}
		
		free(locPanel);
		free(ipvt);
	}
	
}
