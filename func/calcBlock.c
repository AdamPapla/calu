/**
* @file calcBlock.c
* @author Adam Graham
* @brief File containing function to load balance a matrix in chunks of b rows
* 		 for the CALU algorithm
* @date 2023-10-6
*/


#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<math.h>
/**
* @brief Calculates start and end indices of the local chunk of a matrix
* 		 ensuring load balancing with minimal unit being b rows.
* @param[in] n length of list
* @param[in] p number of processors
* @param[in] b number of 
* @param[in] rank rank of processor
* @param[out] *s pointer to integer that will hold starting index
* @param[out] *e pointer to integer that will hold ending index
*/

void calcBlock(int n, int p, int b, int rank, int *s, int *e)
{
	int nOriginal = n;
	if (n%b != 0)	// If n does not divide b, must add so that it does for sake of load balancing
		n += b - n%b;
	int bChunks = n/b;	// Number of chunks of b in matrix (n should always divide b)
	int base_chunk = bChunks/p;	// The minimum number of b chunks per processor
	int chunk = base_chunk;
	int rem = bChunks % p;	 // How many procs will need an additional b chunk
	int start; int end;
	if (rank < rem){
		chunk += 1;
		start = rank + base_chunk*rank; 	// If not all additional elements have been added, rank many of them will have been
	}
	else{
		start = rem + base_chunk * rank;	// If all have been added, can simply add remainder to basechunk*rank
	}
	end = start + chunk;
	start *= b;
	end *= b;
	if (rank == p-1 && nOriginal%b != 0)	// Subtract the rows initially added (as they dont exist)
		end -= b-nOriginal%b;
	*s = start;
	*e = end;
}
