#include <math.h>
#include <omp.h>

double maxsum(int N, int LD, double *A, int NT) {
  double max_sum_of_roots = 0;

  double thread_max_sum_of_roots = 0;
  double sum_of_roots_row_i = 0;
  int thread_id = 0;
  int thread_N = N / NT;
  int start_row = 0;
  int end_row = 0;
  int i = 0, j = 0;

  omp_set_num_threads(NT);
  #pragma omp parallel firstprivate(thread_max_sum_of_roots, sum_of_roots_row_i, \
                                      thread_id, thread_N, start_row, end_row, i, j)
  {
    thread_id = omp_get_thread_num();
    start_row = thread_N * thread_id;
    end_row = (thread_id == NT - 1 ? N : start_row + thread_N);

    for (i = start_row; i < end_row; i++) {
      sum_of_roots_row_i = 0;
      for (j = 0; j < N; j++) {
        sum_of_roots_row_i += sqrt(A[i * LD + j]);
      }

      if (thread_max_sum_of_roots < sum_of_roots_row_i) {
        thread_max_sum_of_roots = sum_of_roots_row_i;
      }
    }
    #pragma omp critical
    {
      if (max_sum_of_roots < thread_max_sum_of_roots) {
        max_sum_of_roots = thread_max_sum_of_roots;
      }
    }
  }

  return max_sum_of_roots;
}
