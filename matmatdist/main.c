#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "matmat.h"
#include "matmatblock.h"
#include "matmatthread.h"

#define DEFAULT_TOLERANCE 1e-3
#define MAX_THREADS 6
#define MAX_THREAD_ROWS 2
#define MAX_THREAD_COLS 4
#define BLOCK_SIZE 256
#define LEADING_DIM 2100

double get_cur_time();

void init_matrix_rand(double *M, const unsigned int LD, const unsigned int rows,
                      const unsigned int cols) {
  unsigned int i, j;
  srand(time(NULL));
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      M[i * LD + j] = (double)(rand() % 1000) * 0.1;
    }
  }
}

void init_matrix_zero(double *M, const unsigned int LD, const unsigned int rows,
                      const unsigned int cols) {
  unsigned int i, j;
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      M[i * LD + j] = 0.0;
    }
  }
}

void copy_matrix(double *dest, const double *src, const unsigned int LD,
                 const unsigned int rows, const unsigned int cols) {
  unsigned int i, j;
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      dest[i * LD + j] = src[i * LD + j];
    }
  }
}

int compare_matrices(const double *A, const double *B, const unsigned int LD,
                     const unsigned int rows, const unsigned int cols,
                     const double tol) {
  unsigned int i, j;
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      double diff = A[i * LD + j] - B[i * LD + j];
      if (diff < -tol || diff > tol) {
        fprintf(stderr, "Mismatch at (%u, %u): %f != %f\n", i, j, A[i * LD + j],
                B[i * LD + j]);
        return 0;
      }
    }
  }
  return 1;
}

void benchmark_matmat_func(matmat_func func, const char *func_name, int ldA,
                           int ldB, int ldC, double *A, double *B, double *C,
                           int N1, int N2, int N3) {
  const unsigned long long num_ops = 2LL * N1 * N2 * N3;
  double time, gflops;
  double start, end;

  start = get_cur_time();
  func(ldA, ldB, ldC, A, B, C, N1, N2, N3);
  end = get_cur_time();

  time = end - start;
  gflops = (double)num_ops / (time * 1e9);

  printf("%s: %f seconds, %f GFLOPS\n", func_name, time, gflops);
}

void benchmark_matmatblock_func(matmatblock_func func, const char *func_name,
                                int ldA, int ldB, int ldC, double *A, double *B,
                                double *C, int N1, int N2, int N3, int dbA,
                                int dbB, int dbC) {
  const unsigned long long num_ops = 2LL * N1 * N2 * N3;
  double time, gflops;
  double start, end;

  start = get_cur_time();
  func(ldA, ldB, ldC, A, B, C, N1, N2, N3, dbA, dbB, dbC);
  end = get_cur_time();

  time = end - start;
  gflops = (double)num_ops / (time * 1e9);

  printf("%s: %f seconds, %f GFLOPS\n", func_name, time, gflops);
}

void benchmark_matmatthread_func(matmatthread_func func, const char *func_name,
                                 int ldA, int ldB, int ldC, double *A,
                                 double *B, double *C, int N1, int N2, int N3,
                                 int dbA, int dbB, int dbC, int NTROW,
                                 int NTCOL) {
  const unsigned long long num_ops = 2LL * N1 * N2 * N3;
  double time, gflops;
  double start, end;

  start = get_cur_time();
  func(ldA, ldB, ldC, A, B, C, N1, N2, N3, dbA, dbB, dbC, NTROW, NTCOL);
  end = get_cur_time();

  time = end - start;
  gflops = (double)num_ops / (time * 1e9);

  printf("%s: %f seconds, %f GFLOPS\n", func_name, time, gflops);
}

void verify_correctness(const double *A, const double *B, const unsigned int LD,
                        const unsigned int N, const double tol,
                        const char *function_a, const char *function_b) {
  if (!compare_matrices(A, B, LD, N, N, tol)) {
    fprintf(stderr, "Error: Result of %s does not match %s!\n", function_a,
            function_b);
    exit(EXIT_FAILURE);
  }
}

int main() {
  matmat_func functions[] = {matmatijk, matmatikj, matmatjik,
                             matmatjki, matmatkij, matmatkji};
  const char *function_names[] = {"matmatijk", "matmatikj", "matmatjik",
                                  "matmatjki", "matmatkij", "matmatkji"};
  const int num_functions = 6;

  matmatblock_func block_functions[] = {matmatblock_ijk, matmatblock_ikj};
  const char *block_function_names[] = {"matmatblock_ijk", "matmatblock_ikj"};
  const int num_block_functions = 2;

  const char *thread_function_name = "matmatthread";

  const unsigned int L = BLOCK_SIZE;
  const unsigned int LD = LEADING_DIM;
  unsigned int N;
  double *A, *B, *C, *C_ref;

  unsigned int NTROW, NTCOL;

  printf("Using block size L=%u and leading dimension LD=%u\n", L, LD);

  for (N = 256; N <= 2048; N += 256) {
    printf("Matrix size N=%u\n", N);

    A = (double *)malloc(LD * N * sizeof(double));
    B = (double *)malloc(LD * N * sizeof(double));
    C = (double *)malloc(LD * N * sizeof(double));
    C_ref = (double *)malloc(LD * N * sizeof(double));

    init_matrix_rand(A, LD, N, N);
    init_matrix_rand(B, LD, N, N);
    init_matrix_zero(C, LD, N, N);

    benchmark_matmat_func(functions[0], function_names[0], LD, LD, LD, A, B, C,
                          N, N, N);
    copy_matrix(C_ref, C, LD, N, N);

    for (int i = 1; i < num_functions; i++) {
      // Skip all but matmatikj (fastest) for brevity
      if (i != 1)
        break;

      // Skip matmatjki and matmatkji (slowest) for brevity
      // if (i == 3 || i == 5)
      //   continue;

      init_matrix_zero(C, LD, N, N);
      benchmark_matmat_func(functions[i], function_names[i], LD, LD, LD, A, B,
                            C, N, N, N);
      verify_correctness(C_ref, C, LD, N, DEFAULT_TOLERANCE, function_names[0],
                         function_names[i]);
    }

    for (int i = 0; i < num_block_functions; i++) {
      init_matrix_zero(C, LD, N, N);
      benchmark_matmatblock_func(block_functions[i], block_function_names[i],
                                 LD, LD, LD, A, B, C, N, N, N, L, L, L);
      verify_correctness(C_ref, C, LD, N, DEFAULT_TOLERANCE, function_names[0],
                         block_function_names[i]);
    }

    for (NTROW = 1; NTROW <= MAX_THREAD_ROWS; NTROW *= 2) {
      for (NTCOL = NTROW; NTCOL <= NTROW * 2; NTCOL *= 2) {
        if (NTROW * NTCOL >= MAX_THREADS)
          break;

        printf("Num threads: %d\n", NTROW * NTCOL);

        init_matrix_zero(C, LD, N, N);
        benchmark_matmatthread_func(matmatthread, thread_function_name, LD, LD,
                                    LD, A, B, C, N, N, N, L, L, L, NTROW,
                                    NTCOL);
        verify_correctness(C_ref, C, LD, N, DEFAULT_TOLERANCE,
                           function_names[0], thread_function_name);
      }
    }

    free(A);
    free(B);
    free(C);
    free(C_ref);

    printf("\n");
  }
  return 0;
}
