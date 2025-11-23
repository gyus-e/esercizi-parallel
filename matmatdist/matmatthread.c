#include "matmatthread.h"
#include "matmatblock.h"
#include <omp.h>

/**
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

  const unsigned int myN1 = N1 / NTROW;
  const unsigned int myN3 = N3 / NTCOL;

  const unsigned int NT = NTROW * NTCOL;
  omp_set_num_threads(NT);

#pragma omp parallel firstprivate(myN1, myN3)
  {
    const unsigned int myID = omp_get_thread_num();
    const unsigned int IDi = myID / NTCOL;
    const unsigned int IDj = myID % NTCOL;

    const unsigned int start_row = myN1 * IDi;
    const unsigned int start_col = myN3 * IDj;

    matmatblock(ldA, ldB, ldC, &A[start_row * ldA], &B[start_col],
                &C[start_row * ldC + start_col], myN1, N2, myN3, dbA, dbB, dbC);
  }
}
