#include "strassen.h"
#include <stdio.h>
#include <stdlib.h>

//#define CHECK_DATA_RANGE

#if __STDC_VERSION__ >= 199901L

#ifdef CHECK_DATA_RANGE
#include <fenv.h>
#include <math.h>
#endif

#endif

// Traditional sequential matrix mutiplication
// Used when matrix size is smaller than 16 to increase efficiency for Strassen
// Algorithm.
void inline sequentialMatrixMult(int n, double **A, double **B, double **C) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      C[i][j] = 0.0;
      for (int k = 0; k < n; k++) {
        C[i][j] += A[i][k] * B[k][j];
#ifdef CHECK_DATA_RANGE
        if (isinf(C[i][j]))
          printf("%s WARNING : DATA IS INFINITY ! \n", __func__);
        if (isnan(C[i][j]))
          printf("%s WARNING : DATA IS Not-A-Number(NaN)! \n", __func__);
#endif
      }
    }
  }
}

// Add two square matrices
void inline matrixAdd(int n, double **A, double **B, double **C) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      C[i][j] = A[i][j] + B[i][j];
#ifdef CHECK_DATA_RANGE
      if (isinf(C[i][j]))
        printf("%s WARNING : DATA IS INFINITY ! \n", __func__);
      if (isnan(C[i][j]))
        printf("%s WARNING : DATA IS Not-A-Number(NaN)! \n", __func__);
#endif
    }
  }
}

// Substract two square matrices
void inline matrixSub(int n, double **A, double **B, double **C) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      C[i][j] = A[i][j] - B[i][j];
#ifdef CHECK_DATA_RANGE
      if (isinf(C[i][j]))
        printf("%s WARNING : DATA IS INFINITY ! \n", __func__);
      if (isnan(C[i][j]))
        printf("%s WARNING : DATA IS Not-A-Number(NaN)! \n", __func__);
#endif
    }
  }
}

// Get sub matrix from a matrix
void inline get_sub_matrix(double **matrix, double **subMatrix, int row_size,
                           int column_size, int row_offset, int column_offset) {
  for (int i = 0; i < row_size; i++) {
    for (int j = 0; j < column_size; j++) {
      subMatrix[i][j] = matrix[i + row_offset][j + column_offset];
    }
  }
}

// Form a matrix from four sub-matrices (C11 C12 C21 C22), and avoid copying
// last row and column if we padded zero for odd N.
void four_blocks_to_matrix(int N, int padding, double **matrix, double **C11,
                           double **C12, double **C21, double **C22) {
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) matrix[i][j] = C11[i][j];

  for (int i = 0; i < N; i++)
    for (int j = 0; j < (N - padding); j++) matrix[i][N + j] = C12[i][j];

  for (int i = 0; i < (N - padding); i++)
    for (int j = 0; j < N; j++) matrix[N + i][j] = C21[i][j];

  for (int i = 0; i < (N - padding); i++)
    for (int j = 0; j < (N - padding); j++) matrix[N + i][N + j] = C22[i][j];
}

// Strassen Algorithm
void Strassen(int n, double **A, double **B, double **C) {
  // Use sequential mutiplication when matrix size is smaller than 16 to
  // increase efficiency. The best threshold might vary between different
  // platforms.
  if (n <= 16) {
    sequentialMatrixMult(n, A, B, C);
    return;
  }

  // Check if n is odd. If yes, we need padding and set padding as 1.
  int padding = n % 2;
  // Make sure N is even.
  int N = n / 2 + padding;

  // Allocate memory for splitting blocks in Strassen Algorithm
  // Cannot move allocation memory to another function for clean code because it
  // will encounter Illegal instruction error
  double **A11 = (double **)malloc(sizeof(double *) * N);
  double **A12 = (double **)malloc(sizeof(double *) * N);
  double **A21 = (double **)malloc(sizeof(double *) * N);
  double **A22 = (double **)malloc(sizeof(double *) * N);
  double **B11 = (double **)malloc(sizeof(double *) * N);
  double **B12 = (double **)malloc(sizeof(double *) * N);
  double **B21 = (double **)malloc(sizeof(double *) * N);
  double **B22 = (double **)malloc(sizeof(double *) * N);
  double **C11 = (double **)malloc(sizeof(double *) * N);
  double **C12 = (double **)malloc(sizeof(double *) * N);
  double **C21 = (double **)malloc(sizeof(double *) * N);
  double **C22 = (double **)malloc(sizeof(double *) * N);

  if ((A11 == NULL) || (A12 == NULL) || (A21 == NULL) || (A22 == NULL) ||
      (B11 == NULL) || (B12 == NULL) || (B21 == NULL) || (B22 == NULL) ||
      (C11 == NULL) || (C12 == NULL) || (C21 == NULL) || (C22 == NULL)) {
    printf("malloc failed ! Out of memory !\n");
    exit(EXIT_FAILURE);
  }

  // Allocate two dimentional array in this fashion to avoid memory
  // fragmentation Use calloc to initialize as 0 for padding
  double *trueA11 = (double *)calloc(N * N, sizeof(double));
  double *trueA12 = (double *)calloc(N * N, sizeof(double));
  double *trueA21 = (double *)calloc(N * N, sizeof(double));
  double *trueA22 = (double *)calloc(N * N, sizeof(double));
  double *trueB11 = (double *)calloc(N * N, sizeof(double));
  double *trueB12 = (double *)calloc(N * N, sizeof(double));
  double *trueB21 = (double *)calloc(N * N, sizeof(double));
  double *trueB22 = (double *)calloc(N * N, sizeof(double));
  double *trueC11 = (double *)calloc(N * N, sizeof(double));
  double *trueC12 = (double *)calloc(N * N, sizeof(double));
  double *trueC21 = (double *)calloc(N * N, sizeof(double));
  double *trueC22 = (double *)calloc(N * N, sizeof(double));

  if ((trueA11 == NULL) || (trueA12 == NULL) || (trueA21 == NULL) ||
      (trueA22 == NULL) || (trueB11 == NULL) || (trueB12 == NULL) ||
      (trueB21 == NULL) || (trueB22 == NULL) || (trueC11 == NULL) ||
      (trueC12 == NULL) || (trueC21 == NULL) || (trueC22 == NULL)) {
    printf("calloc failed ! Out of memory !\n");
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < N; ++i) {
    A11[i] = trueA11;
    A12[i] = trueA12;
    A21[i] = trueA21;
    A22[i] = trueA22;
    B11[i] = trueB11;
    B12[i] = trueB12;
    B21[i] = trueB21;
    B22[i] = trueB22;
    C11[i] = trueC11;
    C12[i] = trueC12;
    C21[i] = trueC21;
    C22[i] = trueC22;
    trueA11 += N;
    trueA12 += N;
    trueA21 += N;
    trueA22 += N;
    trueB11 += N;
    trueB12 += N;
    trueB21 += N;
    trueB22 += N;
    trueC11 += N;
    trueC12 += N;
    trueC21 += N;
    trueC22 += N;
  }

  // Copy original matrix A to local four sub-matrix, and avoid copying the last
  // row and column if we use padding
  get_sub_matrix(A, A11, N, N, 0, 0);
  get_sub_matrix(A, A12, N, (N - padding), 0, N);
  get_sub_matrix(A, A21, (N - padding), N, N, 0);
  get_sub_matrix(A, A22, (N - padding), (N - padding), N, N);
  // Copy original matrix B to local four sub-matrix, and avoid copying the last
  // row and column if we use padding
  get_sub_matrix(B, B11, N, N, 0, 0);
  get_sub_matrix(B, B12, N, (N - padding), 0, N);
  get_sub_matrix(B, B21, (N - padding), N, N, 0);
  get_sub_matrix(B, B22, (N - padding), (N - padding), N, N);

#ifdef FE_ALL_EXCEPT
  // clears all of the exception flags
  feclearexcept(FE_ALL_EXCEPT);
#endif

  // Mutiply matrix A with B and save the result at C by Strassen Algorithm.
  // Refer to "Memory efficient scheduling of Strassen-Winograd’s matrix
  // multiplication algorithm" Table 3: IP schedule for operation C <- A×B in
  // place http://lig-membres.imag.fr/pernet/Publications/fp05-dumas.pdf
  matrixSub(N, A11, A21, C11);
  matrixAdd(N, A21, A22, A21);
  matrixSub(N, B12, B11, C22);
  matrixSub(N, B22, B12, B12);
  Strassen(N, C11, B12, C21);
  matrixSub(N, A21, A11, C12);
  Strassen(N, A11, B11, C11);
  matrixSub(N, B22, C22, B11);
  Strassen(N, A21, C22, A11);
  matrixSub(N, B11, B21, C22);
  Strassen(N, A22, C22, A21);

  matrixSub(N, A12, C12, A22);
  Strassen(N, C12, B11, C22);
  matrixAdd(N, C11, C22, C22);
  Strassen(N, A12, B21, C12);
  matrixAdd(N, C11, C12, C11);
  matrixAdd(N, C22, A11, C12);
  matrixAdd(N, C22, C21, C22);
  matrixSub(N, C22, A21, C21);
  matrixAdd(N, C22, A11, C22);
  Strassen(N, A22, B22, A12);
  matrixAdd(N, C12, A12, C12);

#if defined(FE_OVERFLOW) && defined(FE_UNDERFLOW)
  // Test whether the exception flags overflow and underflow are currently set.
  int raised = fetestexcept(FE_OVERFLOW | FE_UNDERFLOW);
  // Print warning message when overflow or underflow occurred.
  if (raised & FE_OVERFLOW) {
    printf(" Overflow detected !\n");
  }
  if (raised & FE_UNDERFLOW) {
    printf(" Underflow detected !\n");
  }
  // NOTE: This might not work due to compiler and optimization. BUG :
  // https://llvm.org/bugs/show_bug.cgi?id=6050
#endif

  // Put the calculation result of four sub-matrices (C11 C12 C21 C22) into one
  // matrix C Avoid copying the last row and column if we use padding
  four_blocks_to_matrix(N, padding, C, C11, C12, C21, C22);

  // free memory
  // Cannot move free memory to another function for clean code because it will
  // encounter Illegal instruction error
  free(A11[0]);
  free(A12[0]);
  free(A21[0]);
  free(A22[0]);
  free(B11[0]);
  free(B12[0]);
  free(B21[0]);
  free(B22[0]);
  free(C11[0]);
  free(C12[0]);
  free(C21[0]);
  free(C22[0]);
  free(A11);
  free(A12);
  free(A21);
  free(A22);
  free(B11);
  free(B12);
  free(B21);
  free(B22);
  free(C11);
  free(C12);
  free(C21);
  free(C22);
}
