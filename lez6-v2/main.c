#include "maxsum.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define N_arr_size 7
#define NT_arr_size 4

int main(int argc, char *argv[]) {
  int N_arr[] = {800, 1600, 2400, 3200, 4000, 4800, 5400};
  int NT_arr[] = {1, 2, 4, 8};
  int n_idx = 0, nt_idx = 0;

  int N = 0;
  int NT = 0;

  double *A = NULL;
  int i = 0, j = 0;

  double max_sum = 0.0;

  double start_time = 0.0;
  double end_time = 0.0;
  double time = 0.0;

  printf("\n");

  for (n_idx = 0; n_idx < N_arr_size; n_idx++) {

    printf(
        "\n************************************************************\n\n");

    N = N_arr[n_idx];
    printf("Size of matrix: %d x %d\n\n", N, N);

    A = (double *)malloc(N * N * sizeof(double));
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        A[i * N + j] = (double)(rand() % 100);
      }
    }

    for (nt_idx = 0; nt_idx < NT_arr_size; nt_idx++) {
      printf("------------------------------\n");

      NT = NT_arr[nt_idx];
      printf("Number of threads: %d\n", NT);

      start_time = omp_get_wtime();

      max_sum = maxsum(N, N, A, NT);

      end_time = omp_get_wtime();

      time = end_time - start_time;
      printf("Time for maxsum: %f seconds\n", time);
      printf("Max sum: %f\n", max_sum);
    }

    free(A);
  }

  printf("\n");
  return 0;
}
