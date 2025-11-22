#include "matmat.h"

// A: N1 x N2
// B: N2 x N3
// C: N1 x N3

// C = C + A*B
// for i=1..N1, j=1..N3
// C[i][j] = C[i][j] + sum k=1..N2 (A[i][k] * B[k][j])

void matmatblock(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                 int N1, int N2, int N3, int dbA, int dbB, int dbC) {
  double *startA, *startB, *startC;
  int ii, jj, kk;
  for (ii = 0; ii < N1; ii += dbA) {
    for (jj = 0; jj < N3; jj += dbB) {
      for (kk = 0; kk < N2; kk += dbC) {
        startA = A + ii * ldA + kk;
        startB = B + kk * ldB + jj;
        startC = C + ii * ldC + jj;
        matmatijk(ldA, ldB, ldC, startA, startB, startC, N1 / dbA, N2 / dbB,
                  N3 / dbC);
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