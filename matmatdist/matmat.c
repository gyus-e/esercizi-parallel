#include "matmat.h"

// A: N1 x N2
// B: N2 x N3
// C: N1 x N3

// C = C + A*B
// for i=1..N1, j=1..N3
// C[i][j] = C[i][j] + sum k=1..N2 (A[i][k] * B[k][j])

void matmatblockikj(int ldA, int ldB, int ldC, double *A, double *B, double *C,
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

void matmatblock(int ldA, int ldB, int ldC, double *A, double *B, double *C,
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

void matmatijk(int ldA, int ldB, int ldC, double *A, double *B, double *C,
               int N1, int N2, int N3) {
  unsigned int i, j, k;
  for (i = 0; i < N1; i++) {
    for (j = 0; j < N3; j++) {
      for (k = 0; k < N2; k++) {
        C[i * ldC + j] += A[i * ldA + k] * B[k * ldB + j];
      }
    }
  }
}

void matmatikj(int ldA, int ldB, int ldC, double *A, double *B, double *C,
               int N1, int N2, int N3) {
  unsigned int i, j, k;
  for (i = 0; i < N1; i++) {
    for (k = 0; k < N2; k++) {
      for (j = 0; j < N3; j++) {
        C[i * ldC + j] += A[i * ldA + k] * B[k * ldB + j];
      }
    }
  }
}

void matmatjik(int ldA, int ldB, int ldC, double *A, double *B, double *C,
               int N1, int N2, int N3) {
  unsigned int i, j, k;
  for (j = 0; j < N3; j++) {
    for (i = 0; i < N1; i++) {
      for (k = 0; k < N2; k++) {
        C[i * ldC + j] += A[i * ldA + k] * B[k * ldB + j];
      }
    }
  }
}

void matmatjki(int ldA, int ldB, int ldC, double *A, double *B, double *C,
               int N1, int N2, int N3) {
  unsigned int i, j, k;
  for (j = 0; j < N3; j++) {
    for (k = 0; k < N2; k++) {
      for (i = 0; i < N1; i++) {
        C[i * ldC + j] += A[i * ldA + k] * B[k * ldB + j];
      }
    }
  }
}

void matmatkij(int ldA, int ldB, int ldC, double *A, double *B, double *C,
               int N1, int N2, int N3) {
  unsigned int i, j, k;
  for (k = 0; k < N2; k++) {
    for (i = 0; i < N1; i++) {
      for (j = 0; j < N3; j++) {
        C[i * ldC + j] += A[i * ldA + k] * B[k * ldB + j];
      }
    }
  }
}

void matmatkji(int ldA, int ldB, int ldC, double *A, double *B, double *C,
               int N1, int N2, int N3) {
  unsigned int i, j, k;
  for (k = 0; k < N2; k++) {
    for (j = 0; j < N3; j++) {
      for (i = 0; i < N1; i++) {
        C[i * ldC + j] += A[i * ldA + k] * B[k * ldB + j];
      }
    }
  }
}