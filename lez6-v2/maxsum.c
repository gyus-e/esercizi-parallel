#include "maxsum.h"

#include <math.h>
#include <omp.h>

double maxsum(int N, int LD, double *A, int NT) {
  double max_sum_of_roots = 0;

  double thread_max_sum_of_roots = 0;
  double sum_of_roots_row_i = 0;
  int i = 0, j = 0;

  omp_set_num_threads(NT);
  #pragma omp parallel firstprivate(thread_max_sum_of_roots, sum_of_roots_row_i, \
                                      i, j)
  {
    #pragma omp for
    for (i = 0; i < N; i++) {
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
