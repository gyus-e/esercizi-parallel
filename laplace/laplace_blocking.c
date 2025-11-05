#include <mpi.h>

float get_B_ij(float up, float down, float left, float right) {
  return (up + down + left + right) * 0.25;
}

void init_first_row(float *A, float *B, float *daprev, int N, int LD) {
  int j;
  float up, down, left, right;

  for (j = 1; j < N - 1; j++) {
    up = daprev[j];
    down = A[1 * LD + j];
    left = A[0 * LD + (j - 1)];
    right = A[0 * LD + (j + 1)];

    B[0 * LD + j] = get_B_ij(up, down, left, right);
  }
}

void init_row_i(int i, float *A, float *B, int N, int LD) {
  int j;
  float up, down, left, right;

  for (j = 1; j < N - 1; j++) {
    up = A[(i - 1) * LD + j];
    down = A[(i + 1) * LD + j];
    left = A[i * LD + (j - 1)];
    right = A[i * LD + (j + 1)];

    B[i * LD + j] = get_B_ij(up, down, left, right);
  }
}

void init_last_row(float *A, float *B, float *danext, int rows_per_proc, int N,
                   int LD) {
  int j;
  float up, down, left, right;

  for (j = 1; j < N - 1; j++) {
    up = A[(rows_per_proc - 2) * LD + j];
    down = danext[j];
    left = A[(rows_per_proc - 1) * LD + (j - 1)];
    right = A[(rows_per_proc - 1) * LD + (j + 1)];

    B[(rows_per_proc - 1) * LD + j] = get_B_ij(up, down, left, right);
  }
}

void copy_rows(int start_row, int end_row, float *A, float *B, int N, int LD) {
  int i, j;
  for (i = start_row; i <= end_row; i++) {
    for (j = 1; j < N - 1; j++) {
      A[i * LD + j] = B[i * LD + j];
    }
  }
}

void laplace(float *A, float *B, float *daprev, float *danext, int N, int LD,
             int Niter) {
  int i, j, k;
  int nproc, myid;
  int rows_per_proc, start_row, end_row;
  MPI_Status status;

  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  rows_per_proc = N / nproc;
  start_row = 1;
  end_row = rows_per_proc - 2;

  for (k = 0; k < Niter; k++) {
    if (myid > 0) {
      MPI_Send(&A[0 * LD], N, MPI_FLOAT, myid - 1, 0, MPI_COMM_WORLD);
      MPI_Recv(daprev, N, MPI_FLOAT, myid - 1, 0, MPI_COMM_WORLD, &status);
      init_first_row(A, B, daprev, N, LD);
      start_row = 0;
    }

    for (i = 1; i < rows_per_proc - 1; i++) {
      init_row_i(i, A, B, N, LD);
    }

    if (myid < nproc - 1) {
      MPI_Send(&A[(rows_per_proc - 1) * LD], N, MPI_FLOAT, myid + 1, 0,
               MPI_COMM_WORLD);
      MPI_Recv(danext, N, MPI_FLOAT, myid + 1, 0, MPI_COMM_WORLD, &status);
      init_last_row(A, B, danext, rows_per_proc, N, LD);
      end_row = rows_per_proc - 1;
    }

    copy_rows(start_row, end_row, A, B, N, LD);
  }
}