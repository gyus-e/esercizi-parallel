#ifndef MATMATBLOCK_H
#define MATMATBLOCK_H

typedef void (*matmat_block_func)(int, int, int, double *, double *, double *,
                                  int, int, int, int, int, int);

void matmatblockikj(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                    int N1, int N2, int N3, int dbA, int dbB, int dbC);

void matmatblock(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                 int N1, int N2, int N3, int dbA, int dbB, int dbC);

#endif /* MATMATBLOCK_H */
