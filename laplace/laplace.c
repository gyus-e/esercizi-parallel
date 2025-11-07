#include <mpi.h>

void laplace(float *A, float *B, float *daprev, float *danext, int N, int LD,
             int Niter);

void laplace_nb(float *A, float *B, float *daprev, float *danext, int N, int LD,
                int Niter);

void laplace(float *A, float *B, float *daprev, float *danext, int N, int LD,
             int Niter) {
  int i, j, k;
  int nproc, myid;
  int curr_row, prev, next;
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
      curr_row = i * LD;
      prev = curr_row - LD;
      next = curr_row + LD;
      for (j = 1; j < N - 1; j++) {
        B[curr_row + j] =
            0.25 * (A[prev + j] + A[next + j] + A[curr_row + (j - 1)] +
                    A[curr_row + (j + 1)]);
      }
    }

    if (myid > 0) {
      start_row = 0;
      for (j = 1; j < N - 1; j++) {
        B[j] = 0.25 * (daprev[j] + A[LD + j] + A[j - 1] + A[j + 1]);
      }
    }
    if (myid < nproc - 1) {
      end_row = rows_per_proc - 1;
      curr_row = (rows_per_proc - 1) * LD;
      prev = curr_row - LD;
      for (j = 1; j < N - 1; j++) {
        B[curr_row + j] =
            0.25 * (A[prev + j] + danext[j] + A[curr_row + (j - 1)] +
                    A[curr_row + (j + 1)]);
      }
    }

    for (i = start_row; i <= end_row; i++) {
      for (j = 1; j < N - 1; j++) {
        A[i * LD + j] = B[i * LD + j];
      }
    }
  }
}

void laplace_nb(float *A, float *B, float *daprev, float *danext, int N, int LD,
                int Niter) {
  int i, j, k;
  int nproc, myid;
  int curr_row, prev, next;
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
      MPI_Irecv(daprev, N, MPI_FLOAT, myid - 1, 1, MPI_COMM_WORLD, &recv_first);
    }

    if (myid < nproc - 1) {
      MPI_Isend(&A[(rows_per_proc - 1) * LD], N, MPI_FLOAT, myid + 1, 1,
                MPI_COMM_WORLD, &send_last);
      MPI_Irecv(danext, N, MPI_FLOAT, myid + 1, 0, MPI_COMM_WORLD, &recv_last);
    }

    for (i = 1; i < rows_per_proc - 1; i++) {
      curr_row = i * LD;
      prev = curr_row - LD;
      next = curr_row + LD;
      for (j = 1; j < N - 1; j++) {
        B[curr_row + j] =
            0.25 * (A[prev + j] + A[next + j] + A[curr_row + (j - 1)] +
                    A[curr_row + (j + 1)]);
      }
    }

    if (myid > 0) {
      start_row = 0;
      MPI_Wait(&recv_first, &status_send_first);
      for (j = 1; j < N - 1; j++) {
        B[j] = 0.25 * (daprev[j] + A[LD + j] + A[j - 1] + A[j + 1]);
      }
    }
    if (myid < nproc - 1) {
      end_row = rows_per_proc - 1;
      curr_row = (rows_per_proc - 1) * LD;
      prev = curr_row - LD;
      MPI_Wait(&recv_last, &status_send_last);
      for (j = 1; j < N - 1; j++) {
        B[curr_row + j] =
            0.25 * (A[prev + j] + danext[j] + A[curr_row + (j - 1)] +
                    A[curr_row + (j + 1)]);
      }
    }

    if (myid > 0) {
      MPI_Wait(&send_first, &status_send_first);
    }
    if (myid < nproc - 1) {
      MPI_Wait(&send_last, &status_send_last);
    }

    for (i = start_row; i <= end_row; i++) {
      for (j = 1; j < N - 1; j++) {
        A[i * LD + j] = B[i * LD + j];
      }
    }
  }
}
