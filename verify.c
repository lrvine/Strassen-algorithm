#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "strassen.h"

// Use to verify the correctness and efficiency comparison
int main ( int argc, char ** argv )
{
	// decide if we want to verify answer
	int verification = 1;
	// input square matrix's size N
	int N = 512;

	// assign input square matrix's size N
	if( argc >= 2 )
		N = atoi(argv[1]);

	// decide whether we want to run sequential matrix multiplication for verification.
	if( argc >= 3 )
		verification = atoi(argv[2]);

	// allocate two dimentional array
	double **A = (double **)malloc(sizeof(double*) * N);// input matrix A
	double **B = (double **)malloc(sizeof(double*) * N);// input matrix B
	double **C = (double **)malloc(sizeof(double*) * N);// answer matrix C
	double **T = (double **)malloc(sizeof(double*) * N);// Store the correct answer for verifying the correctness of C

	if( (A==NULL) || (B==NULL) || (C==NULL) || (T==NULL) ){ 
		printf("malloc failed ! Out of memory !\n");
		exit(EXIT_FAILURE);
	}

	// allocate in this fasion to avoid memory fragmentation
	double *trueA = (double*)malloc(sizeof(double) * N * N);
	double *trueB = (double*)malloc(sizeof(double) * N * N);
	double *trueC = (double*)malloc(sizeof(double) * N * N);
	double *trueT = (double*)malloc(sizeof(double) * N * N);

	if( (trueA==NULL) || (trueB==NULL) || (trueC==NULL) || (trueT==NULL) ){ 
		printf("malloc failed ! Out of memory !\n");
		exit(EXIT_FAILURE);
	}

	for(int i=0; i<N; ++i) {
    		A[i] = trueA;
    		B[i] = trueB;
    		C[i] = trueC;
    		T[i] = trueT;
    		trueA += N;
    		trueB += N;
    		trueC += N;
    		trueT += N;
	}

	// assign value into matrices A and B for testing
	for(int i=0; i<N; ++i) {
		for(int j=0; j<N; ++j) {
			A[i][j]=i;
			B[i][j]=j;
		}
	}	

	clock_t begin;
        clock_t end;	
        double traditional_time_spent;// store calculation time for traditional sequential method	
        double fastMatrixMultiplication_time_spent;// store calculation time for this implementation	

	if( verification ){
	// Calculate correct answer and check time spent 
		begin = clock();
		sequentialMatrixMult(N,A,B,T); //C=S(A,B)
		end = clock();
		traditional_time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	}

	// Calculate answer by new implementation and check time spent
	begin = clock();
	Strassen(N,A,B,C);  //C=S(A,B)
	end = clock();
	fastMatrixMultiplication_time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

	//If we want to verify the answer, compare each element in the matrices.
	if( verification ){
		for(int i=0; i<N; ++i) {
			for(int j=0; j<N; ++j) {
				if( T[i][j] != C[i][j] ){ 
					printf("ERROR !! Calculation result DID NOT match! Correct is %F, and we got %F\n",C[i][j],T[i][j]);
					exit(EXIT_FAILURE);
				}
			}
		}	
		printf( "\nTraditional method would spend %F seconds for %d X %d matrices.\n",traditional_time_spent,N,N);
		printf( "\nThis implementation by Strassen Algorithm spent %F seconds for %d X %d matrices.\n\n",fastMatrixMultiplication_time_spent,N,N);
	}else{
		printf( "\nThis implementation by Strassen Algorithm spent %F seconds for %d X %d matrices.\n\n",fastMatrixMultiplication_time_spent,N,N);
	}
	//free memory
	free(A[0]);
	free(B[0]);
	free(C[0]);
	free(T[0]);
	free(A);
	free(B);
	free(C);
	free(T);
	return 0;
}

// Traditional sequential matrix mutiplication
// Used for correctness verification and efficiency comparison.
// Need a copy here to avoid compiler warning for inline function
void inline sequentialMatrixMult(int n, double ** A, double ** B, double ** C) {
	for (int i = 0 ; i < n ; i++ ) {
		for (int j = 0 ; j < n ; j++)  {
			C[i][j] = 0.0;
			for (int k = 0 ; k < n ; k++ )
				C[i][j] += A[i][k] * B[k][j];
		}
	}		
}
