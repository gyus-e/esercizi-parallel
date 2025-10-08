#include "maxsum.h"

#include <stdio.h>
#include <math.h>
#include <omp.h>

double maxsum(int N, int LD, double *A, int NT)
{
	double max_sum_of_roots = 0;
	int row_with_max_sum = 0;

	omp_set_num_threads(NT);

	#pragma omp parallel
	{
		double thread_max_sum_of_roots = 0;
		int thread_row_with_max_sum = 0;

		#pragma omp for
		for (int i = 0; i < N; i++)
		{
			double sum_of_roots_row_i = 0;

			for (int j = 0; j < N; j++)
			{
				sum_of_roots_row_i += sqrt(A[i * LD + j]);
			}

			if (thread_max_sum_of_roots < sum_of_roots_row_i)
			{
				thread_max_sum_of_roots = sum_of_roots_row_i;
				thread_row_with_max_sum = i;
			}
		}
		#pragma omp critical
		{
			if (max_sum_of_roots < thread_max_sum_of_roots)
			{
				max_sum_of_roots = thread_max_sum_of_roots;
				row_with_max_sum = thread_row_with_max_sum;
			}
		}
	}
	printf("\n");
	printf("max sum of roots = %f\n", max_sum_of_roots);
	printf("in row: %d\n", row_with_max_sum);
	return max_sum_of_roots;
}
