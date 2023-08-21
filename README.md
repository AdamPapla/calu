# Communication Avoiding LU Factorisation (CALU) on a 1d process grid  
Adam Graham  
17-08-23  
Based on the paper "CALU: A Communication Optimal LU Factorization Algorithm" by L.Grigori et. al.
****

## Overview  
This algorithm provides an implementation of communication avoiding LU factorisation, capable of 
performing a stable factorisation of any arbitrary n x m matrix distributed among p processes.  
MPI is used for inter process communication. The algorithm implemented is a block LU factorisation
algorithm which; 
- factorises a panel of b columns using TSLU with tournament pivoting
- updates the corresponding b rows
- updates the Schurs complement (trailing matrix)
- Proceeds recursively on this trailing matrix  
The tournament pivoting scheme massively reduces communication from that in parallel Gaussian
Elimination with Partial Pivoting (GEPP), while mantaining the high level of stability GEPP offers.  
CBLAS functions are used for local operations where possible, owing to their high level of optimisation.  
For a detailed report on the implementation and results of this algorithm, see the report.  
For a description of all functions used in the algorithm, see the Doxygen documentation.  
## Usage
The __main.c__ file provided shows an example of how the CALU algorithm can be used. Note that for CALU, it is 
assumed that MPI has been initialised with `MPI_Init`, and that the matrix has already been appropriately 
distributed among processes (or read in in parallel) according to the load balancing given by the `calcBlock` function.  
For a list of arguments and returns, see the Doxygen documentation. A brief outline is provided below  

    void calu(double *locA, double *locb, int n, int m, int b)

- `locA`: An array of doubles containing the process' chunk of the A matrix to be factorised
- `locb`: An array of doubles containing either the rhs vector of a system of linear equations
or an indexed array for [0:n-1] that will track row pivoting
- `n`: Number of rows in matrix A
- `m`: Number of columns in matrix A
- `b`: Block size for algorithms recursion  
Returns an in place LU factorisation, where each process' local chunk now holds a section of the LU matrix. U is stored in
diagonal and upper triangular entries of the global matrix. L is stored in lower triangular entries of global matrix.  
The __random.c__ file can be readily used to generate binary files of random matrices A and b such that Ax=b. x is a constant vector of
user defined value. Allows for LU factorisation to be quickly verified.
The __main.c__ file can be readily used to test CALU on these generated matrices. However, the user must ensure that the 
_n_ and _m_ values of the generated matrices specified in __random.c__ match those specified in **main.c**.
## Function file Descriptions
Function files are found in the __func__ directory. Doxygen documentation is provided, as well as in file commenting
- __calcBlock.c__ Calculates load balancing of rows of matrix in blocks of b between processes. Accounts
for cases where n does not divide b.
- __matrixIO.c__ Contains a variety of functions involving IO, including MPI IO for reading/writing to files
in parallel
- __matrix.c__ Contains miscellaneous functions concerning matrices, eg. extracting blocks of the matrix.
- __pivot.c__ Contains functions concerning the pivoting of matrix rows, both in serial and in parallel.
- __luFunc.c__ Contains functions concerning the serial and parallel unpivoted LU factorisation of matrices
- __tslu.c__ Defines the tournament pivoting function, which determines the best b pivot rows of a panel
by performing a binomial tree reduction, comparing two process' best b rows at each step.
- __calu.c__ Defines the full CALU algorithm. Performs TSLU on panels of the matrix, updates relevant blocks and 
handles all elements of recursion, both intra process recursion and inter process recursion. Performs the 
factorisation in place. For further details, see the report.
- __randomGen.c__ Contains functions for generating binary files containing matrices A and b such that Ax=b for some
constant vector x. This allows the LU factorisation to be quickly validated by using it to solve the system and 
checking that x has the desired constant entries. 
- __oldFunc.c__ Functions used at one point in the implementation, but optimisations rendered them unnecessary.
Included as reference to certain sections of the report.
 


