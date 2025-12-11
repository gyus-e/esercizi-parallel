#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

static void copy_matrix_row(double *src, double *dst, int rows, int cols, int ldsrc);
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

void matmatdist(MPI_Comm Gridcom, // comunicatore della griglia di processi
                int LDA, int LDB, int LDC,       // leading dimension
                double *A, double *B, double *C, // matrici
                int N1, int N2, int N3,          // dimensioni matrici
                int DB1, int DB2, int DB3,       // dimensioni blocchi
                int NTrow, int NTcol) {          // config. thread
  int myid, NP, dims[2], periods[2], coords[2];
  MPI_Comm rowcom, colcom;

  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &NP);
  MPI_Cart_get(Gridcom, 2, dims, periods, coords);
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

  // dimensioni dei buffer per le comunicazioni
  const unsigned int A_col_dim = localN1 * localN2;
  const unsigned int B_row_dim = localN2 * localN3;

  // buffer per le comunicazioni
  double *Acol = (double *)malloc(sizeof(double) * A_col_dim);
  double *Brow = (double *)malloc(sizeof(double) * B_row_dim);

  unsigned int k, row_sender_id, col_sender_id;
  for (k = 0; k < K2; k++) {
    row_sender_id = k % dims[0];
    col_sender_id = k % dims[1];

    if (coords[0] == row_sender_id) {
      copy_matrix_row(A, Acol, localN1, localN2, LDA);
    }
    MPI_Bcast(Acol, A_col_dim, MPI_DOUBLE, col_sender_id, rowcom);
    
    if (coords[1] == col_sender_id) {
      copy_matrix_row(B, Brow, localN2, localN3, LDB);
    }
    MPI_Bcast(Brow, B_row_dim, MPI_DOUBLE, row_sender_id, colcom);

    matmatthread(LDA, LDB, LDC, A, B, C, localN1, localN2, localN3, DB1, DB2,
                 DB3, NTrow, NTcol);
  }

  free(Acol);
  free(Brow);
}

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

void matmatblock(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                 int N1, int N2, int N3, int dbA, int dbB, int dbC) {
  dbA = MCD(N1, dbA);
  dbB = MCD(N2, dbB);
  dbC = MCD(N3, dbC);

  const unsigned int num_submatrixes_A = dbA > 0 ? N1 / dbA : 1;
  const unsigned int num_submatrixes_B = dbC > 0 ? N3 / dbC : 1;
  const unsigned int num_subsubmatrixes = dbB > 0 ? N2 / dbB : 1;

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

static void copy_matrix_row(double *src, double *dst, int rows, int cols, int ldsrc) {
  unsigned int i, j, dst_idx = 0;
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      dst[dst_idx++] = src[i * ldsrc + j];
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