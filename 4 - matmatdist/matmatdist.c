#include <mpi.h>
#include <omp.h>
#include <stdlib.h>

static void copy_submatrix(double *src, double *dst, int rows, int cols,
                           int ldsrc, int lddst);
inline static int mcm(const unsigned int A, const unsigned int B);
inline static int MCD(const unsigned int A, const unsigned int B);
static int euclid(unsigned int A, unsigned int B);

void matmatikj(int ldA, int ldB, int ldC, double *A, double *B, double *C,
               int N1, int N2, int N3);
void matmatblock(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                 int N1, int N2, int N3, int dbA, int dbB, int dbC);
void matmatthread(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                  int N1, int N2, int N3, int dbA, int dbB, int dbC, int NTROW,
                  int NTCOL);
void matmatdist(MPI_Comm Gridcom, int LDA, int LDB, int LDC, double *A,
                double *B, double *C, int N1, int N2, int N3, int DB1, int DB2,
                int DB3, int NTrow, int NTcol);

/*
  Moltiplicazione distribuita di matrici in una griglia di processi K1 x K3.
  K2 è il minimo comune multiplo di K1 e K3.

  La matrice A viene suddivisa in K1 x K2 blocchi di dimensione (N1/K1) x (N2/K2), 
  la matrice B viene suddivisa in K2 x K3 blocchi di dimensione (N2/K2) x (N3/K3), 
  la matrice C viene suddivisa in K1 x K3 blocchi di dimensione (N1/K1) x (N3/K3).

  Le sottomatrici vengono distribuite ai processi della griglia in modo ciclico.
  Al processo con coordinate (i, j) vengono assegnati i seguenti blocchi:
  - tutti i blocchi A(i, k) dove k % K3 == j;
  - tutti i blocchi B(k, j) dove k % K1 == i;
  - il blocco C(i, j).

  Per calcolare il blocco C(i, j), il processo (i, j) ha bisogno di tutti
  i blocchi A(i, k) e B(k, j), per k = 0, ..., K2-1. Per ogni k,
  - il processo che possiede il blocco A(i, k) lo invia a tutti i processi
    della riga x della griglia,
  - il processo che possiede il blocco B(k, j) lo invia a tutti i processi
    della colonna y della griglia. 
  - Ogni processo calcola il prodotto dei blocchi ricevuti 
    e somma il risultato al proprio blocco C(i, j).

  Gridcom - Comunicatore della griglia di processi
  LDA, LDB, LDC - Leading dimension delle tre matrici
  A, B, C - Puntatori al primo elemento delle tre matrici
  N1 - Numero di righe di A e C
  N2 - Numero di colonne di A e righe di B
  N3 - Numero di colonne di B e C
  DB1 - Righe di A e C contenute in un blocco di cache
  DB2 - Colonne di A e righe di B contenute in un blocco di cache
  DB3 - Colonne di B e C contenute in un blocco di cache
  NTrow - Numero di righe della griglia di thread
  NTcol - Numero di colonne della griglia di thread
*/
void matmatdist(MPI_Comm Gridcom, int LDA, int LDB, int LDC, double *A,
                double *B, double *C, int N1, int N2, int N3, int DB1, int DB2,
                int DB3, int NTrow, int NTcol) {
  int myid, NP, dims[2], periods[2], mycoords[2];
  MPI_Comm rowcom, colcom;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &NP);
  MPI_Cart_get(Gridcom, 2, dims, periods, mycoords);
  MPI_Cart_sub(Gridcom, (int[]){0, 1}, &rowcom);
  MPI_Cart_sub(Gridcom, (int[]){1, 0}, &colcom);

  // dimensioni della griglia di processi
  const unsigned int K1 = dims[0];
  const unsigned int K2 = mcm(dims[0], dims[1]);
  const unsigned int K3 = dims[1];

  // dimensioni dei blocchi di matrici
  const unsigned int localN1 = N1 / K1;
  const unsigned int localN2 = N2 / K2;
  const unsigned int localN3 = N3 / K3;

  const unsigned int A_block_dim = localN1 * localN2;
  const unsigned int B_block_dim = localN2 * localN3;

  double *A_block = (double *)malloc(sizeof(double) * A_block_dim);
  double *B_block = (double *)malloc(sizeof(double) * B_block_dim);

  unsigned int row_sender_id, col_sender_id;
  unsigned int i, k, j;
  unsigned int itr;
  for (itr = 0; itr < K2; itr++) {
    row_sender_id = itr % dims[0];
    col_sender_id = itr % dims[1];

    i = row_sender_id * localN1;
    k = itr * localN2;
    j = col_sender_id * localN3;

    if (mycoords[0] == row_sender_id) {
      copy_submatrix(&A[i * LDA + k], A_block, localN1, localN2,
                     LDA, localN2);
    }
    MPI_Bcast(A_block, A_block_dim, MPI_DOUBLE, col_sender_id, rowcom);

    if (mycoords[1] == col_sender_id) {
      copy_submatrix(&B[k * LDB + j], B_block, localN2, localN3,
                     LDB, localN3);
    }
    MPI_Bcast(B_block, B_block_dim, MPI_DOUBLE, row_sender_id, colcom);

    matmatthread(localN2, localN3, LDC, A_block, B_block, C, localN1, localN2, localN3, DB1, DB2,
                 DB3, NTrow, NTcol);
  }

  free(A_block);
  free(B_block);
}

/*
  Moltiplicazione parallela di matrici in una griglia di thread NTROW x NTCOL,
  Calcola le dimensioni delle sottomatrici da assegnare a ogni thread.
  Il thread con identificativo (IDi, IDj)
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
  Calcola il numero di sottomatrici della stessa dimensione dei blocchi di
  cache. Per ogni sottomatrice C(ii, jj) di dimensione dbA x dbC, effettua il
  prodotto righe per colonne di:
  - sottomatrice A(ii, 0) di dimensione dbA x N2,
  - sottomatrice B(0, jj) di dimensione N2 x dbC,
  suddividendole a loro volta in:
  - sottomatrice A(ii, kk) di dimensione dbA x dbB,
  - sottomatrice B(kk, jj) di dimensione dbB x dbC,
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
  dbA = MCD(N1, dbA);
  dbB = MCD(N2, dbB);
  dbC = MCD(N3, dbC);

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

static void copy_submatrix(double *src, double *dst, int Nrows, int Ncols,
                           int ldsrc, int lddst) {
  unsigned int i, j;
  for (i = 0; i < Nrows; i++) {
    for (j = 0; j < Ncols; j++) {
      dst[i * lddst + j] = src[i * ldsrc + j];
    }
  }
}

inline static int mcm(const unsigned int A, const unsigned int B) {
  return (A * B) / euclid(A, B);
}

inline static int MCD(const unsigned int A, const unsigned int B) {
  return euclid(A, B);
}

static int euclid(unsigned int A, unsigned int B) {
  int Q, R;
  while (B != 0) {
    Q = A / B;
    R = A % B;
    A = B;
    B = R;
  }
  return A;
}