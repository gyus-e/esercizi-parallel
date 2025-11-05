#include <mpi.h>
#include <stdio.h>

int somma(int *A, int N);

int main(int argc, char **argv) {
  int A[500];
  MPI_Init(&argc, &argv);

  int sum = somma(A, 10);
  printf("Sum = %d\n", sum);

  MPI_Finalize();
  return 0;
}

int somma(int *A, int N) {
  MPI_Status status;
  int rank, size;
  int i, start, end;
  int next, prev;
  int token;
  int sum = 0;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  start = rank * N;
  end = (rank + 1) * N;
  for (i = start; i < end; i++) {
    A[i] = i;
    sum += A[i];
  }

  next = (rank + size + 1) % size;
  prev = (rank + size - 1) % size;
  token = sum;
  for (i = 0; i < size - 1; i++) {
    if (i % 2 == 0) {
      MPI_Send(&token, 1, MPI_INT, next, 10, MPI_COMM_WORLD);
      MPI_Recv(&token, 1, MPI_INT, prev, 10, MPI_COMM_WORLD, &status);
    } else {
      MPI_Recv(&token, 1, MPI_INT, prev, 10, MPI_COMM_WORLD, &status);
      MPI_Send(&token, 1, MPI_INT, next, 10, MPI_COMM_WORLD);
    }
    sum += token;
  }
}

return sum;
}
