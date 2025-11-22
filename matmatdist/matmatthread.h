#ifndef MATMATTHREAD_H
#define MATMATTHREAD_H

typedef void (*matmatthread_func)(int, int, int, double *, double *, double *,
                                  int, int, int, int, int, int, int, int);

void matmatthread(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                  int N1, int N2, int N3, int dbA, int dbB, int dbC, int NTROW,
                  int NTCOL);

#endif // MATMATTHREAD_H