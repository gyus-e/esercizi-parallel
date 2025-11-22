#include "matmatblock.h"
#include "matmat.h"

/*
A: N1 x N2
B: N2 x N3
C: N1 x N3

N1: rows in A and C
dbA: rows in A and C per block.

N2: columns in A and rows in B (the shared dimension)
dbB: columns in A and rows in B per block (the number of partial products).

N3: columns in B and C
dbC: columns in B and C per block.
*/

void matmatblock_ijk(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                     int N1, int N2, int N3, int dbA, int dbB, int dbC) {
  const unsigned int num_blocks_rows = N1 / dbA;
  const unsigned int num_blocks_products = N2 / dbB;
  const unsigned int num_blocks_columns = N3 / dbC;

  unsigned int row_A, col_B, curr;
  unsigned int idxA, idxB, idxC;
  unsigned int ii, jj, kk;
  for (ii = 0; ii < num_blocks_rows; ii++) {
    row_A = ii * dbA;
    for (jj = 0; jj < num_blocks_columns; jj++) {
      col_B = jj * dbC;
      idxC = row_A * ldC + col_B;
      for (kk = 0; kk < num_blocks_products; kk++) {
        curr = kk * dbB;
        idxA = row_A * ldA + curr;
        idxB = curr * ldB + col_B;
        matmatijk(ldA, ldB, ldC, &A[idxA], &B[idxB], &C[idxC], dbA, dbB, dbC);
      }
    }
  }
}

void matmatblock_ikj(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                     int N1, int N2, int N3, int dbA, int dbB, int dbC) {
  const unsigned int num_blocks_rows = N1 / dbA;
  const unsigned int num_blocks_products = N2 / dbB;
  const unsigned int num_blocks_columns = N3 / dbC;

  unsigned int row_A, col_B, curr;
  unsigned int idxA, idxB, idxC;
  unsigned int ii, jj, kk;
  for (ii = 0; ii < num_blocks_rows; ii++) {
    row_A = ii * dbA;
    for (kk = 0; kk < num_blocks_products; kk++) {
      curr = kk * dbB;
      idxA = row_A * ldA + curr;
      for (jj = 0; jj < num_blocks_columns; jj++) {
        col_B = jj * dbC;
        idxB = curr * ldB + col_B;
        idxC = row_A * ldC + col_B;
        matmatikj(ldA, ldB, ldC, &A[idxA], &B[idxB], &C[idxC], dbA, dbB, dbC);
      }
    }
  }
}