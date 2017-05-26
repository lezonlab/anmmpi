#ifndef BLZSOLVE_H
#define BLZSOLVE_H
int blzSolveSparse(double *A, int *IRN, int *JCN, int ne, double *VAL, double **VEC, int nn, int lo_eig, int hi_eig, int nvbset);
#endif
