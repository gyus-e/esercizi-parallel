#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "matmatthread.h"

#define DEFAULT_TOLERANCE 1e-3
#define MAX_THREADS 6
#define MAX_THREAD_ROWS 2
#define MAX_THREAD_COLS 4
#define BLOCK_SIZE 256
#define LEADING_DIM 2100

typedef void (*matmat_func)(int, int, int, double *, double *, double *, int,
                            int, int);

typedef void (*matmatblock_func)(int, int, int, double *, double *, double *,
                                 int, int, int, int, int, int);

typedef void (*matmatthread_func)(int, int, int, double *, double *, double *,
                                  int, int, int, int, int, int, int, int);

double get_cur_time();
void init_matrix_rand(double *M, const unsigned int LD, const unsigned int rows,
                      const unsigned int cols);
void init_matrix_zero(double *M, const unsigned int LD, const unsigned int rows,
                      const unsigned int cols);
void copy_matrix(double *dest, const double *src, const unsigned int LD,
                 const unsigned int rows, const unsigned int cols);
int compare_matrices(const double *C_ref, const double *C,
                     const unsigned int LD, const unsigned int rows,
                     const unsigned int cols, const double tol);
void benchmark_matmat_func(matmat_func func, const char *func_name, int ldA,
                           int ldB, int ldC, double *A, double *B, double *C,
                           int N1, int N2, int N3);
void benchmark_matmatblock_func(matmatblock_func func, const char *func_name,
                                int ldA, int ldB, int ldC, double *A, double *B,
                                double *C, int N1, int N2, int N3, int dbA,
                                int dbB, int dbC);
void benchmark_matmatthread_func(matmatthread_func func, const char *func_name,
                                 int ldA, int ldB, int ldC, double *A,
                                 double *B, double *C, int N1, int N2, int N3,
                                 int dbA, int dbB, int dbC, int NTROW,
                                 int NTCOL);
void verify_correctness(const double *C_ref, const double *C,
                        const unsigned int LD, const unsigned int N,
                        const double tol, const char *ref_function_name,
                        const char *function_name);

int main() {
  matmat_func functions[] = {matmatijk, matmatikj, matmatjik,
                             matmatjki, matmatkij, matmatkji};
  const char *function_names[] = {"matmatijk", "matmatikj", "matmatjik",
                                  "matmatjki", "matmatkij", "matmatkji"};
  const int num_functions = 6;

  const unsigned int L = BLOCK_SIZE;
  const unsigned int LD = LEADING_DIM;
  unsigned int N;
  double *A, *B, *C, *C_ref;

  unsigned int NTROW, NTCOL;

  printf("* Block size L=%u\n* Leading dimension LD=%u\n", L, LD);

  for (N = 256; N <= 2048; N *= 2) {

    printf("*********************************************************\n");
    printf("* Matrix size N=%u\n", N);
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

    init_matrix_zero(C, LD, N, N);
    benchmark_matmatblock_func(matmatblock, "matmatblock", LD, LD, LD, A, B, C,
                               N, N, N, L, L, L);
    verify_correctness(C_ref, C, LD, N, DEFAULT_TOLERANCE, function_names[0],
                       "matmatblock");

    for (NTROW = 1; NTROW <= MAX_THREAD_ROWS; NTROW *= 2) {
      for (NTCOL = NTROW; NTCOL <= NTROW * 2; NTCOL *= 2) {
        if (NTROW * NTCOL > MAX_THREADS)
          break;

        printf("\n* Threads grid: %dx%d\n", NTROW, NTCOL);

        init_matrix_zero(C, LD, N, N);
        benchmark_matmatthread_func(matmatthread, "matmatthread", LD, LD, LD, A,
                                    B, C, N, N, N, L, L, L, NTROW, NTCOL);
        verify_correctness(C_ref, C, LD, N, DEFAULT_TOLERANCE,
                           function_names[0], "matmatthread");
      }
    }

    free(A);
    free(B);
    free(C);
    free(C_ref);
  }
  return 0;
}

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

int compare_matrices(const double *C_ref, const double *C,
                     const unsigned int LD, const unsigned int rows,
                     const unsigned int cols, const double tol) {
  unsigned int i, j;
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      double diff = C_ref[i * LD + j] - C[i * LD + j];
      if (diff < -tol || diff > tol) {
        fprintf(stderr, "Mismatch:\n C_ref(%u, %u) = %f\n C(%u, %u) = %f\n", i,
                j, C_ref[i * LD + j], i, j, C[i * LD + j]);
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

  printf("%s: \t %f seconds, \t %f GFLOPS\n", func_name, time, gflops);
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

  printf("%s: \t %f seconds, \t %f GFLOPS\n", func_name, time, gflops);
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

  printf("%s: \t %f seconds, \t %f GFLOPS\n", func_name, time, gflops);
}

void verify_correctness(const double *C_ref, const double *C,
                        const unsigned int LD, const unsigned int N,
                        const double tol, const char *ref_function_name,
                        const char *function_name) {
  if (!compare_matrices(C_ref, C, LD, N, N, tol)) {
    fprintf(stderr, "Error:\n Result of %s does not match %s.\n",
            ref_function_name, function_name);
    // exit(EXIT_FAILURE);
  }
}
