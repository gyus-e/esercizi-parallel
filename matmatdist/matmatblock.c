#include "matmatblock.h"
#include "matmat.h"

/**
  Calcolare numero di sottomatrici della stessa dimensione dei blocchi di cache.
  Per ogni sottomatrice C(ii, jj) di dimensione dbA x dbC,
  effetturare prodotto righe per colonne di
  sottomatrice A(ii, 0) di dimensione dbA x N2
  sottomatrice B(0, jj) di dimensione N2 x dbC
  suddividendole a loro volta in
  sottomatrice A(ii, kk) di dimensione dbA x dbB
  sottomatrice B(ii, kk) di dimensione dbB x dbC
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
