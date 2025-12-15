#include <mpi.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>

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

void matmatdist(MPI_Comm Gridcom, int LDA, int LDB, int LDC, double *A,
                double *B, double *C, int N1, int N2, int N3, int DB1, int DB2,
                int DB3, int NTrow, int NTcol) {
  int my_id, grid_dims[2], grid_periods[2], my_coords[2];
  const int row_dir[2] = {0, 1};
  const int col_dir[2] = {1, 0};
  MPI_Comm rowcom, colcom;

  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Cart_get(Gridcom, 2, grid_dims, grid_periods, my_coords);
  MPI_Cart_sub(Gridcom, row_dir, &rowcom);
  MPI_Cart_sub(Gridcom, col_dir, &colcom);

  const unsigned int grid_dims_mcm = mcm(grid_dims[0], grid_dims[1]);

  const unsigned int blockN1 = N1 / grid_dims[0];
  const unsigned int blockN2 = N2 / grid_dims_mcm;
  const unsigned int blockN3 = N3 / grid_dims[1];

  const unsigned int A_block_dim = blockN1 * blockN2;
  const unsigned int B_block_dim = blockN2 * blockN3;

  double *A_block = (double *)malloc(sizeof(double) * A_block_dim);
  double *B_block = (double *)malloc(sizeof(double) * B_block_dim);

  int k, c, r, i;

  for (k = 0; k < grid_dims_mcm; k++) {
    r = k % grid_dims[0];
    c = k % grid_dims[1];

    if (my_coords[0] == r) {
      for (int i = 0; i < blockN2; i++) {
        memcpy(&B_block[i * blockN3], &B[i * LDB], sizeof(double) * blockN3);
      }
    }
    
    if (my_coords[1] == c) {
      for (i = 0; i < blockN1; i++) {
        memcpy(&A_block[i * blockN2], &A[i * LDA], sizeof(double) * blockN2);
      }
    }

    MPI_Bcast(B_block, B_block_dim, MPI_DOUBLE, r, colcom);
    MPI_Bcast(A_block, A_block_dim, MPI_DOUBLE, c, rowcom);

    matmatthread(blockN2, blockN3, LDC, A_block, B_block, C,
                 blockN1, blockN2, blockN3, DB1, DB2, DB3, NTrow, NTcol);
  }

  free(A_block);
  free(B_block);
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