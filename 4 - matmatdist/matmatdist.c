#include <mpi.h>
#include <omp.h>

#define min(a, b) (((a) < (b)) ? (a) : (b))

void matmatikj(int ldA, int ldB, int ldC, double *A, double *B, double *C,
               int N1, int N2, int N3);
void matmatblock(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                 int N1, int N2, int N3, int dbA, int dbB, int dbC);
void matmatthread(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                  int N1, int N2, int N3, int dbA, int dbB, int dbC, int NTROW,
                  int NTCOL);
void matmatdist(MPI_Comm Gridcom, 
                int LDA, int LDB, int LDC, 
                double *A, double *B, double *C, 
                int N1, int N2, int N3, 
                int DB1, int DB2, int DB3, 
                int NTrow, int Ntcol);

void matmatdist(MPI_Comm Gridcom,                // comunicatore della griglia di processi
                int LDA, int LDB, int LDC,       // leading dimension
                double *A, double *B, double *C, // matrici
                int N1, int N2, int N3,          // dimensioni matrici
                int DB1, int DB2, int DB3,       // dimensioni blocchi
                int NTrow, int Ntcol) {          // config. thread
  return;
  
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
