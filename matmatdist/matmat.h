#ifndef MATMAT_H
#define MATMAT_H

typedef void (*matmat_func)(int, int, int, double *, double *, double *, int,
                            int, int);

typedef void (*matmat_block_func)(int, int, int, double *, double *, double *,
                                  int, int, int, int, int, int);

void matmatblockikj(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                    int N1, int N2, int N3, int dbA, int dbB, int dbC);

void matmatblock(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                 int N1, int N2, int N3, int dbA, int dbB, int dbC);

void matmatijk(int ldA, int ldB, int ldC, double *A, double *B, double *C,
               int N1, int N2, int N3);
void matmatikj(int ldA, int ldB, int ldC, double *A, double *B, double *C,
               int N1, int N2, int N3);
void matmatjik(int ldA, int ldB, int ldC, double *A, double *B, double *C,
               int N1, int N2, int N3);
void matmatjki(int ldA, int ldB, int ldC, double *A, double *B, double *C,
               int N1, int N2, int N3);
void matmatkij(int ldA, int ldB, int ldC, double *A, double *B, double *C,
               int N1, int N2, int N3);
void matmatkji(int ldA, int ldB, int ldC, double *A, double *B, double *C,
               int N1, int N2, int N3);

#endif /* MATMAT_H */