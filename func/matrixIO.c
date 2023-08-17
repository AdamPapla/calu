/**
* @file matrixIO.c
* @author Adam Graham
* @brief Functions for input/output of matrices to files/stdout
* @date 10-7-2023
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include "matrix.h"
#include "calcBlock.h"
#include "matrixIO.h"

/**
 * @brief 	Use MPI IO to read the correct rows of the overall matrix from a binary file
 * @param[in]	filename	Name of binary file containing matrix
 * @param[out]	A			Matrix to store local chunk of matrix read from file
 * @param[in]	m			Number of columns in matrix
 * @param[in]	chunkStart	Start row index of process' chunk of the matrix
 * @param[in]	chunkEnd	End row index of process' chunk of the matrix
 * @param[in]	nprocs		Number of processes in communicator
 * @param[in]	comm		Communicator
 */
void read_chunk_from_file(const char* filename, double* A, int m, int chunkStart, int chunkEnd, int nprocs, int rank, MPI_Comm comm) {
    int num_rows = chunkEnd - chunkStart;
    int num_elements = num_rows * m;
    // Open the file collectively
    MPI_File file;
    MPI_File_open(comm, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

    // Set the file view for each process to read its chunk of data
    MPI_Offset displacement = chunkStart * m * sizeof(double);
    MPI_File_set_view(file, displacement, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
	MPI_Status status;
    MPI_File_read_all(file, A, num_elements, MPI_DOUBLE, &status);
    // Close the file
    MPI_File_close(&file);

}

/**
 * @brief	Reads an array from file
 * @param[in]	filename		Name of file containing array of doubles
 * @param[in]	rows			Number of rows in matrix being read in
 * @param[in]	columns			Number of columns in matrix being read in
 * @param[out]	a				Array containing matrix
 */
void readArray(const char *filename, double *a, int rows, int columns) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error opening file: %s\n", filename);
        return;
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            if (fscanf(file, "%lf", &a[i * columns + j]) != 1) {
                printf("Error reading data from file.\n");
                fclose(file);
                return;
            }
        }
    }

    fclose(file);
	return;
}

/**
 * @brief	Write an array to a file
 * @param[in]	filename	Filename to store the array in
 * @param[in]	A			Array to write to file
 * @param[in]	n			Number of rows in A
 * @param[in]	m			Number of columns in A
 */
void writeToFile(const char* filename, double *A, int n, int m) {
    // Open the file for writing
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening the file.\n");
        return;
    }


    // Generate and write the matrix elements to the file
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            // Write the matrix entry to the file
            fprintf(file, "%lf ", A[i*m+j]);
        }
        fprintf(file, "\n");
    }

    // Close the file
    fclose(file);
}

/**
 * @brief	Print matrix of doubles to stdout
 * @param[in]	A		Array holding matrix
 * @param[in]	n		Number of rows in matrix
 * @param[in]	m		Number of columns in matrix
 */
void print_mat(double *A, int n, int m){
	for (int i=0; i<n; i++){
		for (int j=0; j<m; j++){
			printf("%4.1lf ", A[i*m+j]);
		}
		printf("\n");
	}
}
/**
 * @brief	Print matrix of doubles that has been distributed among processes
 * 			to stdout
 * @param[in]	A		Local chunk of matrix to be printed
 * @param[in]	n		Number of rows in local chunk
 * @param[in]	m		Number of columns in local chunk
 * @param[in]	nprocs	Number of processes in communicator
 * @param[in]	rank	Rank or process in communicator
 */
void printFull(double *A, int n, int m, int nprocs, int rank, MPI_Comm comm){
	for (int i=0; i<nprocs; i++){
		if (rank == i)
			print_mat(A, n, m);
		MPI_Barrier(comm);
	}
}
/**
 * @brief	Print matrix of integers to stdout
 * @param[in]	A		Array holding matrix
 * @param[in]	n		Number of rows in matrix
 * @param[in]	m		Number of columns in matrix
 */

void intPrintMat(int *A, int m, int n){
	for (int i=0; i<m; i++){
		for (int j=0; j<n; j++){
			printf("%d ", A[i*n+j]);
		}
		printf("\n");
	}
}
/**
 * @brief	Print matrix of integers that has been distributed among processes
 * 			to stdout
 * @param[in]	A		Local chunk of matrix to be printed
 * @param[in]	n		Number of rows in local chunk
 * @param[in]	m		Number of columns in local chunk
 * @param[in]	nprocs	Number of processes in communicator
 * @param[in]	rank	Rank or process in communicator
 */

void intPrintFull(int *A, int n, int m, int nprocs, int rank, MPI_Comm comm){
	for (int i=0; i<nprocs; i++){
		if (rank == i)
			intPrintMat(A, n, m);
		MPI_Barrier(comm);
	}
}

