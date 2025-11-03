#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 2

int main(int argc, char **argv) {
  int rank = 0;
  int size = 0;
  MPI_Status status;
  float message[N];
  int i = 0;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0) {

    srand(time(NULL));

    printf("Process %d will send values", rank);
    for (i = 0; i < N; i++) {
      message[i] = (float)(rand() % 100);
      printf(" %.2f", message[i]);
    }
    printf("\n");
    MPI_Send(message, N, MPI_FLOAT, 1, 10, MPI_COMM_WORLD);

  } else if (rank == 1) {

    MPI_Recv(message, N, MPI_FLOAT, 0, 10, MPI_COMM_WORLD, &status);
    printf("Process %d received values", rank);
    for (i = 0; i < N; i++) {
      printf(" %.2f", message[i]);
    }
    printf("\n");
  }

  MPI_Finalize();
  return 0;
}
