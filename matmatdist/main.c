#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "matmat.h"

double get_cur_time();

void init_matrix(double *M, int LD, int rows, int cols) {
  srand(time(NULL));
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      M[i * LD + j] = (double)(rand() % 1000) * 0.1;
    }
  }
}

void benchmark_matmat_func(matmat_func func, const char *func_name, int ldA,
                           int ldB, int ldC, double *A, double *B, double *C,
                           int N1, int N2, int N3) {
  long long num_ops = 2LL * N1 * N2 * N3;
  double time, gflops;
  double start, end;

  start = get_cur_time();
  func(ldA, ldB, ldC, A, B, C, N1, N2, N3);
  end = get_cur_time();

  time = end - start;
  gflops = (double)num_ops / (time * 1e9);

  printf("%s: %f seconds, %f GFLOPS\n", func_name, time, gflops);
}

void benchmark_matmat_block_func(matmat_block_func func, const char *func_name,
                                 int ldA, int ldB, int ldC, double *A,
                                 double *B, double *C, int N1, int N2, int N3,
                                 int dbA, int dbB, int dbC) {
  long long num_ops = 2LL * N1 * N2 * N3;
  double time, gflops;
  double start, end;

  start = get_cur_time();
  func(ldA, ldB, ldC, A, B, C, N1, N2, N3, dbA, dbB, dbC);
  end = get_cur_time();

  time = end - start;
  gflops = (double)num_ops / (time * 1e9);

  printf("%s: %f seconds, %f GFLOPS\n", func_name, time, gflops);
}

int main() {
  matmat_func functions[] = {matmatijk, matmatikj, matmatjik,
                             matmatjki, matmatkij, matmatkji};
  char *function_names[] = {"matmatijk", "matmatikj", "matmatjik",
                            "matmatjki", "matmatkij", "matmatkji"};

  double *A, *B, *C;
  int N, LD;
  int L = 64; // block size, must divide N
  for (N = 256; N <= 1536; N += 256) {
    printf("Matrix size N=%d\n", N);

    LD = N + 16;
    A = (double *)malloc(LD * N * sizeof(double));
    B = (double *)malloc(LD * N * sizeof(double));
    C = (double *)malloc(LD * N * sizeof(double));

    init_matrix(A, LD, N, N);
    init_matrix(B, LD, N, N);

    init_matrix(C, LD, N, N);
    benchmark_matmat_block_func(matmatblock, "matmatblock", LD, LD, LD, A, B, C,
                                N, N, N, L, L, L);

    for (int i = 0; i < 6; i++) {
      if (i == 3 || i == 5)
        continue; // Skip matmatjki and matmatkji for brevity
      init_matrix(C, LD, N, N);
      benchmark_matmat_func(functions[i], function_names[i], LD, LD, LD, A, B,
                            C, N, N, N);
    }
    free(A);
    free(B);
    free(C);
  }

  return 0;
}