#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
int main(){
   double get_cur_time();
   int NT, N, i, j, LD;
   double MAX, *A ;
   double t1, t2, save;
   double maxsum(int , int , double* , int);
  	
   for (N = 800; N <= 3200; N*=2) {
      printf("===============\n");
      printf(" DIMENSIONE DELLA MATRICE = %d \n", N);
      
      LD = N;
      A = (double*)malloc(sizeof(double)*LD*LD);

      for (i = 0; i < N; i++) {
         for (j = 0; j < N; j++) {
            A[i*LD+j] = (rand()%100);
         }
      }

      for (NT = 1; NT <=8; NT = NT*2){

         printf("===============\n");
         printf(" NUMERO THREAD = %d \n", NT);

         t1 = get_cur_time();
         MAX = maxsum(N, LD, A, NT);
         t2 = get_cur_time();

         if (NT == 1) save = t2-t1;
         printf("il massimo della somma dei moduli con N = %d e' %f \n", N,  MAX);
         printf("il tempo totale e' %e , lo speedup = %f , l'efficienza = %f \n", t2-t1, save/(t2-t1), save/(t2-t1)/NT);
      
         }

      free(A);
   }
}




