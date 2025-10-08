#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "maxsum.h"

#define NT_arr_size 4
#define N_arr_size 7

int main(int argc, char *argv[])
{
	int NT_arr[] = {1, 2, 4, 8};
	int N_arr[] = {800, 1600, 2400, 3200, 4000, 4800, 5400};
	double **A;
	// srand(time(NULL));

	printf("\n");

	for (int n_idx = 0; n_idx < N_arr_size; n_idx++)
	{

		printf("\n************************************************************\n\n");

		int N = N_arr[n_idx];
		printf("Size of matrix: %d x %d\n\n", N, N);

		A = (double **)malloc(N * sizeof(double *));
		for (int i = 0; i < N; i++)
		{
			A[i] = (double *)malloc(N * sizeof(double));
			for (int j = 0; j < N; j++)
			{
				A[i][j] = (double)(rand() % 100);
			}
		}

		for (int nt_idx = 0; nt_idx < NT_arr_size; nt_idx++)
		{
			printf("------------------------------\n");

			int NT = NT_arr[nt_idx];
			printf("Number of threads: %d\n", NT);

			double start_time = omp_get_wtime();

			maxsum(N, N, A, NT);

			double end_time = omp_get_wtime();

			double time = end_time - start_time;
			printf("Time for maxsum: %.6f seconds\n", time);
		}

		for (int i = 0; i < N; i++)
		{
			free(A[i]);
		}
		free(A);
	}

	printf("\n");
	return 0;
}
