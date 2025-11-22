#ifndef MATMATBLOCK_H
#define MATMATBLOCK_H

typedef void (*matmatblock_func)(int, int, int, double *, double *, double *,
                                 int, int, int, int, int, int);

void matmatblock_ikj(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                     int N1, int N2, int N3, int dbA, int dbB, int dbC);

void matmatblock_ijk(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                     int N1, int N2, int N3, int dbA, int dbB, int dbC);

#endif /* MATMATBLOCK_H */
