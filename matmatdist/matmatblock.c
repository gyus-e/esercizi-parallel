#include "matmatblock.h"
#include "matmat.h"

void matmatblock_ikj(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                     int N1, int N2, int N3, int dbA, int dbB, int dbC) {
  const unsigned int num_blocks_A = N1 / dbA;
  const unsigned int num_blocks_B = N2 / dbB;
  const unsigned int num_blocks_C = N3 / dbC;

  unsigned int row, col, tmp;
  unsigned int idxA, idxB, idxC;
  unsigned int i, j, k;
  for (i = 0; i < num_blocks_A; i++) {
    row = i * dbA;
    for (k = 0; k < num_blocks_B; k++) {
      tmp = k * dbB;
      idxA = row * ldA + tmp;
      for (j = 0; j < num_blocks_C; j++) {
        col = j * dbC;
        idxC = row * ldC + col;
        idxB = tmp * ldB + col;
        matmatikj(ldA, ldB, ldC, &A[idxA], &B[idxB], &C[idxC], dbA, dbB, dbC);
      }
    }
  }
}

void matmatblock_ijk(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                     int N1, int N2, int N3, int dbA, int dbB, int dbC) {
  const unsigned int num_blocks_A = N1 / dbA;
  const unsigned int num_blocks_B = N2 / dbB;
  const unsigned int num_blocks_C = N3 / dbC;

  unsigned int row, col, tmp;
  unsigned int idxA, idxB, idxC;
  unsigned int i, j, k;
  for (i = 0; i < num_blocks_A; i++) {
    row = i * dbA;
    for (j = 0; j < num_blocks_C; j++) {
      col = j * dbC;
      idxC = row * ldC + col;
      for (k = 0; k < num_blocks_B; k++) {
        tmp = k * dbB;
        idxA = row * ldA + tmp;
        idxB = tmp * ldB + col;
        matmatijk(ldA, ldB, ldC, &A[idxA], &B[idxB], &C[idxC], dbA, dbB, dbC);
      }
    }
  }
}
