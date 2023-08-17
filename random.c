#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mkl_types.h>
#include <mkl_cblas.h>
#include <mkl.h>
#include <mpi.h>

#include "func/randomGen.h"
#include "func/matrixIO.h"

int main(int argc, char * argv[]){
	MPI_Init(&argc, &argv);

	int n = 10000;
	int m = 10000;

	generateRMFile("A.bin", "b.bin", n, m, 1);

	MPI_Finalize();
	return 0;

}
