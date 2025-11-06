#include <mpi.h>

void laplace(float *A, float *B, float *daprev, float *danext, int N, int LD,
             int Niter);

void laplace_nb(float *A, float *B, float *daprev, float *danext, int N, int LD,
                int Niter);

static inline void init_first_row(float *A, float *B, float *daprev, int N,
                                  int LD);

static inline void init_last_row(float *A, float *B, float *danext,
                                 int rows_per_proc, int N, int LD);

static inline void init_row_i(int i, float *A, float *B, int N, int LD);

static inline void copy_rows(int start_row, int end_row, float *A, float *B,
                             int N, int LD);

static inline float get_B_ij(float up, float down, float left, float right);

void laplace(float *A, float *B, float *daprev, float *danext, int N, int LD,
             int Niter) {
  int i, k;
  int nproc, myid;
  int rows_per_proc, start_row, end_row;
  int idx;
  MPI_Status status;

  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  rows_per_proc = N / nproc;
  start_row = 1;
  end_row = rows_per_proc - 2;

  for (k = 0; k < Niter; k++) {
    if (myid > 0) {
      if (myid % 2 == 1) {
        MPI_Send(&A[0], N, MPI_FLOAT, myid - 1, 0, MPI_COMM_WORLD);
        MPI_Recv(daprev, N, MPI_FLOAT, myid - 1, 0, MPI_COMM_WORLD, &status);
      } else {
        MPI_Recv(daprev, N, MPI_FLOAT, myid - 1, 0, MPI_COMM_WORLD, &status);
        MPI_Send(&A[0], N, MPI_FLOAT, myid - 1, 0, MPI_COMM_WORLD);
      }
    }

    if (myid < nproc - 1) {
      idx = (rows_per_proc - 1) * LD;
      if (myid % 2 == 1) {
        MPI_Send(&A[idx], N, MPI_FLOAT, myid + 1, 0, MPI_COMM_WORLD);
        MPI_Recv(danext, N, MPI_FLOAT, myid + 1, 0, MPI_COMM_WORLD, &status);
      } else {
        MPI_Recv(danext, N, MPI_FLOAT, myid + 1, 0, MPI_COMM_WORLD, &status);
        MPI_Send(&A[idx], N, MPI_FLOAT, myid + 1, 0, MPI_COMM_WORLD);
      }
    }

    for (i = 1; i < rows_per_proc - 1; i++) {
      init_row_i(i, A, B, N, LD);
    }

    if (myid > 0) {
      start_row = 0;
      init_first_row(A, B, daprev, N, LD);
    }
    if (myid < nproc - 1) {
      end_row = rows_per_proc - 1;
      init_last_row(A, B, danext, rows_per_proc, N, LD);
    }

    copy_rows(start_row, end_row, A, B, N, LD);
  }
}

void laplace_nb(float *A, float *B, float *daprev, float *danext, int N, int LD,
                int Niter) {
  int i, k;
  int nproc, myid;
  int rows_per_proc, start_row, end_row;
  MPI_Status status_send_first, status_send_last;
  MPI_Request send_first, send_last, recv_first, recv_last;

  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  rows_per_proc = N / nproc;
  start_row = 1;
  end_row = rows_per_proc - 2;

  for (k = 0; k < Niter; k++) {
    if (myid > 0) {
      MPI_Isend(&A[0], N, MPI_FLOAT, myid - 1, 0, MPI_COMM_WORLD, &send_first);
      MPI_Irecv(daprev, N, MPI_FLOAT, myid - 1, 0, MPI_COMM_WORLD, &recv_first);
    }

    if (myid < nproc - 1) {
      MPI_Isend(&A[(rows_per_proc - 1) * LD], N, MPI_FLOAT, myid + 1, 0,
                MPI_COMM_WORLD, &send_last);
      MPI_Irecv(danext, N, MPI_FLOAT, myid + 1, 0, MPI_COMM_WORLD, &recv_last);
    }

    for (i = 1; i < rows_per_proc - 1; i++) {
      init_row_i(i, A, B, N, LD);
    }

    if (myid > 0) {
      start_row = 0;
      MPI_Wait(&recv_first, &status_send_first);
      init_first_row(A, B, daprev, N, LD);
    }
    if (myid < nproc - 1) {
      end_row = rows_per_proc - 1;
      MPI_Wait(&recv_last, &status_send_last);
      init_last_row(A, B, danext, rows_per_proc, N, LD);
    }

    if (myid > 0) {
      MPI_Wait(&send_first, &status_send_first);
    }
    if (myid < nproc - 1) {
      MPI_Wait(&send_last, &status_send_last);
    }
    copy_rows(start_row, end_row, A, B, N, LD);
  }
}

static inline void init_first_row(float *A, float *B, float *daprev, int N,
                                  int LD) {
  int j;
  for (j = 1; j < N - 1; j++) {
    B[j] = get_B_ij(daprev[j], A[LD + j], A[j - 1], A[j + 1]);
  }
}

static inline void init_last_row(float *A, float *B, float *danext,
                                 int rows_per_proc, int N, int LD) {
  int j;
  int idx = (rows_per_proc - 1) * LD;
  int prev = idx - LD;
  for (j = 1; j < N - 1; j++) {
    B[(rows_per_proc - 1) * LD + j] =
        get_B_ij(A[prev + j], danext[j], A[idx + (j - 1)], A[idx + (j + 1)]);
  }
}

static inline void init_row_i(int i, float *A, float *B, int N, int LD) {
  int j;
  for (j = 1; j < N - 1; j++) {
    B[i * LD + j] = get_B_ij(A[(i - 1) * LD + j], A[(i + 1) * LD + j],
                             A[i * LD + (j - 1)], A[i * LD + (j + 1)]);
  }
}

static inline void copy_rows(int start_row, int end_row, float *A, float *B,
                             int N, int LD) {
  int i, j;
  for (i = start_row; i <= end_row; i++) {
    for (j = 1; j < N - 1; j++) {
      A[i * LD + j] = B[i * LD + j];
    }
  }
}

static inline float get_B_ij(float up, float down, float left, float right) {
  return (up + down + left + right) * 0.25;
}