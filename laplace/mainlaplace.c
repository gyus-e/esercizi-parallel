#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int main(int argc, char *argv[]) {
  int nproc, myid, prev, next;
  int N, i, j, ifirst, iter, Niter, LD, local_row;
  float *A, *Anew, *daprev, *danext;
  MPI_Status status;
  double get_cur_time(), t1, t2;
  void laplace(float *, float *, float *, float *, int, int, int);

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  printf("hello from %d di %d processi \n", myid, nproc);
  sleep(1);

  N = 400;
  Niter = 8000;
  LD = 500;
  A = (float *)malloc(500 * 500 * sizeof(float));
  Anew = (float *)malloc(500 * 500 * sizeof(float));
  daprev = (float *)malloc(500 * sizeof(float));
  danext = (float *)malloc(500 * sizeof(float));

  // inizializzazione matrice

  for (i = 0; i < N / nproc; i++) { // tutta la matrice locale = 0
    for (j = 0; j < N; j++) {
      A[i * LD + j] = 0.;
    }
  }
  if (myid == 0)
    for (j = 0; j < N; j++)
      A[0 * LD + j] = j; // prima riga matrice del proc id=0  da 0 a 390

  if (myid == nproc - 1)
    for (j = 0; j < N; j++)
      A[(N / nproc - 1) * LD + j] =
          N - 1 - j; // ultima riga matrice del proc id=nproc-1 da 390 a 0

  ifirst = myid * N / nproc;
  for (i = 0; i < N / nproc; i++) {
    A[i * LD + 0] =
        ifirst + i; // bordo sinistro da ifirst a ilast-1 in ogni proc
    A[i * LD + N - 1] = N - 1 - A[i * LD + 0]; // A[i][0] + A[i][N-1] = 0 sempre
  }

  if (myid == 0)
    printf("\n esecuzione con N = %d  e %d iterazioni\n\n", N, Niter);

  sleep(1);

  printf("processo %d ha righe da %d a %d \n", myid, ifirst,
         ifirst + N / nproc - 1);

  t1 = get_cur_time();

  laplace(A, Anew, daprev, danext, N, LD, Niter);

  t2 = get_cur_time();

  if (myid == 0)
    printf("con %d processi, il tempo e' %f\n", nproc, t2 - t1);

  sleep(1);

  if (myid == 0) {
    local_row = 1;
    printf("prima  %d -->   A[1][1]=%f  A[1][398]=%f  \n", myid, A[1 * LD + 1],
           A[1 * LD + 398]);
  }
  if (myid == nproc / 2 - 1) {
    local_row = 199 - ifirst;
    printf("centro %d -->   A[199][199]=%f  A[199][200]=%f  \n", myid,
           A[local_row * LD + 199], A[local_row * LD + 200]);
  }
  if (myid == nproc / 2) {
    local_row = 200 - ifirst;
    printf("centro %d -->   A[200][199]=%f  A[200][200]=%f  \n", myid,
           A[local_row * LD + 199], A[local_row * LD + 200]);
  }
  if (myid == nproc - 1) {
    local_row = 398 - ifirst;
    printf("ultima %d -->   A[398][1]=%f  A[398][398]=%f  \n", myid,
           A[local_row * LD + 1], A[local_row * LD + 398]);
  }

  MPI_Finalize();
}
