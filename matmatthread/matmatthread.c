#include <omp.h>

#define min(a, b) (((a) < (b)) ? (a) : (b))

void matmatijk(int ldA, int ldB, int ldC, double *A, double *B, double *C,
               int N1, int N2, int N3);
void matmatikj(int ldA, int ldB, int ldC, double *A, double *B, double *C,
               int N1, int N2, int N3);
void matmatjik(int ldA, int ldB, int ldC, double *A, double *B, double *C,
               int N1, int N2, int N3);
void matmatjki(int ldA, int ldB, int ldC, double *A, double *B, double *C,
               int N1, int N2, int N3);
void matmatkij(int ldA, int ldB, int ldC, double *A, double *B, double *C,
               int N1, int N2, int N3);
void matmatkji(int ldA, int ldB, int ldC, double *A, double *B, double *C,
               int N1, int N2, int N3);
void matmatblock(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                 int N1, int N2, int N3, int dbA, int dbB, int dbC);
void matmatthread(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                  int N1, int N2, int N3, int dbA, int dbB, int dbC, int NTROW,
                  int NTCOL);

/*
  In una griglia di thread NTROW x NTCOL,
  Calcolare dimensione delle sottomatrici per ogni thread.
  Il thread con identificativo (IDi, Idj)
  calcola sottomatrice C(IDi, IDj) di dimensione myN1 x myN3
  prendendo sottomatrice A(IDi, 0) di dimensione myN1 x N2
  e sottomatrice di B(0, IDj) di dimensione N2 x myN3.

  ldA, ldB, ldC - Leading dimension delle tre matrici
  A, B, C - Puntatori al primo elemento delle tre matrici
  N1 - Numero di righe di A e C
  N2 - Numero di colonne di A e righe di B
  N3 - Numero di colonne di B e C
  dbA - Righe di A e C contenute in un blocco di cache
  dbB - Colonne di A e righe di B contenute in un blocco di cache
  dbC - Colonne di B e C contenute in un blocco di cache
  NTROW - Numero di righe della griglia di thread
  NTCOL - Numero di colonne della griglia di thread
*/
void matmatthread(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                  int N1, int N2, int N3, int dbA, int dbB, int dbC, int NTROW,
                  int NTCOL) {

  const unsigned int NT = NTROW * NTCOL;
  const unsigned int myN1 = N1 / NTROW;
  const unsigned int myN3 = N3 / NTCOL;

  unsigned int myID, IDi, IDj;
  unsigned int start_row, start_col;

  omp_set_num_threads(NT);

  #pragma omp parallel private(myID, IDi, IDj, start_row, start_col)
  {
    myID = omp_get_thread_num();
    IDi = myID / NTCOL;
    IDj = myID % NTCOL;

    start_row = myN1 * IDi;
    start_col = myN3 * IDj;

    matmatblock(ldA, ldB, ldC, &A[start_row * ldA], &B[start_col],
                &C[start_row * ldC + start_col], myN1, N2, myN3, dbA, dbB, dbC);
  }
}

/*
  Calcolare numero di sottomatrici della stessa dimensione dei blocchi di cache.
  Per ogni sottomatrice C(ii, jj) di dimensione dbA x dbC,
  effetturare prodotto righe per colonne di
  sottomatrice A(ii, 0) di dimensione dbA x N2
  sottomatrice B(0, jj) di dimensione N2 x dbC
  suddividendole a loro volta in
  sottomatrice A(ii, kk) di dimensione dbA x dbB
  sottomatrice B(kk, jj) di dimensione dbB x dbC
  effettuando per ognuna di esse il prodotto righe per colonne
  e sommando il risultato alla sottomatrice C(ii, jj).

  ldA, ldB, ldC - Leading dimension delle tre matrici
  A, B, C - Puntatori al primo elemento delle tre matrici
  N1 - Numero di righe di A e C
  N2 - Numero di colonne di A e righe di B
  N3 - Numero di colonne di B e C
  dbA - Righe di A e C contenute in un blocco di cache
  dbB - Colonne di A e righe di B contenute in un blocco di cache
  dbC - Colonne di B e C contenute in un blocco di cache
*/
void matmatblock(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                 int N1, int N2, int N3, int dbA, int dbB, int dbC) {
  dbA = min(N1, dbA);
  dbB = min(N2, dbB);
  dbC = min(N3, dbC);

  const unsigned int num_submatrixes_A = N1 / dbA;
  const unsigned int num_submatrixes_B = N3 / dbC;
  const unsigned int num_subsubmatrixes = N2 / dbB;

  unsigned int row_A, col_B, curr_subsubmatrix;
  unsigned int idxA, idxB, idxC;
  unsigned int ii, jj, kk;
  for (ii = 0; ii < num_submatrixes_A; ii++) {
    row_A = ii * dbA;
    for (jj = 0; jj < num_submatrixes_B; jj++) {
      col_B = jj * dbC;
      idxC = row_A * ldC + col_B;
      for (kk = 0; kk < num_subsubmatrixes; kk++) {
        curr_subsubmatrix = kk * dbB;
        idxA = row_A * ldA + curr_subsubmatrix;
        idxB = curr_subsubmatrix * ldB + col_B;
        matmatikj(ldA, ldB, ldC, &A[idxA], &B[idxB], &C[idxC], dbA, dbB, dbC);
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
