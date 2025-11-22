#include "matmat.h"

// A: N1 x N2
// B: N2 x N3
// C: N1 x N3

// C = C + A*B
// for i=1..N1, j=1..N3
// C[i][j] = C[i][j] + sum k=1..N2 (A[i][k] * B[k][j])

void matmatblock(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                 int N1, int N2, int N3, int dbA, int dbB, int dbC) {
  const unsigned int blocksA = N1 / dbA;
  const unsigned int blocksB = N2 / dbB;
  const unsigned int blocksC = N3 / dbC;
  unsigned int idxA, idxB, idxC;
  unsigned int i, j, k;
  unsigned int ii, jj, kk;
  for (i = 0; i < blocksA; i++) {
    for (j = 0; j < blocksC; j++) {
      ii = i * dbA;
      jj = j * dbC;
      idxC = ii * ldC + jj;
      for (k = 0; k < blocksB; k++) {
        kk = k * dbB;
        idxA = ii * ldA + kk;
        idxB = kk * ldB + jj;
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