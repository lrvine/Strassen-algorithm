#ifndef STRASSEN_H
#define STRASSEN_H

void inline sequentialMatrixMult(int n, double ** A, double ** B, double ** C);
void inline matrixAdd(int n, double ** A, double ** B, double ** C);
void inline matrixSub(int n, double ** A, double ** B, double ** C);
void inline get_sub_matrix( double ** matrix, double ** subMatrix, int row_size,int column_size, int row_offset, int column_offset );
void four_blocks_to_matrix( int N, int padding, double ** matrix, double ** C11, double ** C12, double ** C21, double ** C22 );
void Strassen(int n, double** A, double** B, double** C);
void fastMatrixMultiplication(int n, double** A, double** B, double** C);

#endif
