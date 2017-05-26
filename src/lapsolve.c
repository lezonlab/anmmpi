/*****************************************************************************/
/*                                                                           */
/*        Finds eigenvalues and eigenvectors of a matrix using LAPACK        */
/*                                                                           */
/*                       Created 03/03/2017 by Tim Lezon                     */
/*                                                                           */
/*****************************************************************************/
/*
 * MM[1..lda][1..lda] is the matrix to be decomposed.
 * VAL[1..num_eigs] contains eigenvalues on output
 * VEC[1..nn][1..num_eigs] contains eigenvectors on output
 * lda is the size of the array allocated for MM. This is not necessarily the
 *     size of the matrix that is being decomposed: lda >= nn
 * nn is the size of the matrix to be decomposed. The program expects that 
 *     MM[1..nn][1..nn] is a real, symmetric matrix.
 * lo_eig is the index of the lowest requested eigenvalue (lo_eig>0)
 * hi_eig is the index of the greatest requested eigenvalue (hi_eig<=nn)
 * num_eigs = (hi_eig - lo_eig + 1)
 */
#include "lapsolve.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int lapsolve(double **MM, double *VAL, double **VEC, long lda, long nn, long lo_eig, long hi_eig)
{
  long num_eigs;
  int i,j,k;

  /* LAPACK variables */
  char *jobz="Vectors",*uplo="Upper";
  int len_jobz=1,len_uplo=1;
  double *WORK;
  double wkopt=0.0;
  long lwork=0, info=0;

  /* dsyevx-specific */
  char *range="Inbetween";
  int len_range=1;
  double vl=0.0,vu=0.0;
  int *IWORK,*IFAIL;
  double *Z;
  double abstol = 0.0;

  if(lo_eig<1){
    fprintf(stderr,"\n*ERROR: In 'lapsolve', invalid lo_eig value (%ld)\n",lo_eig);
    return -1;}
  if(hi_eig>nn){
    fprintf(stderr,"\n*ERROR: In 'lapsolve', invalid hi_eig value (%ld)\n",hi_eig);
    return -1;}
  if(hi_eig < lo_eig){
    fprintf(stderr,"\n*ERROR: In 'lapsolve', invalid range [%ld, %ld]\n",lo_eig,hi_eig);
    return -1;}

  /* NOTE: Eigenvalues will be computed most accurately when ABSTOL is
          set to twice the underflow threshold 2*DLAMCH('S'), not zero. */
  abstol = 2.0*dlamch_("S");
  num_eigs = hi_eig - lo_eig + 1;

  Z = calloc(nn*num_eigs,sizeof(double));
  if(!Z){
    fprintf(stderr,"\nfailure allocating Z in lapsolve\n\n");
    return -1;}
  IWORK = malloc( 5*nn*sizeof(int) );
  if(!IWORK){
    fprintf(stderr,"\nfailure allocating IWORK in lapsolve\n\n");
    return -1;}
  IFAIL = malloc( nn*sizeof(int) );
  if(!IFAIL){
    fprintf(stderr,"\nfailure allocating IFAIL in lapsolve\n\n");
    return -1;}

  /* Call 'dsyevx' once with lwork=-1 to find the optimal block size. 
     Then allocate memory for the WORK array. */
  lwork = -1;
  dsyevx_(jobz, range, uplo, &nn, &MM[1][1], &lda, &vl, &vu, &lo_eig, &hi_eig, &abstol, &num_eigs, &VAL[1], Z, &nn, &wkopt, &lwork, IWORK, IFAIL, &info );

  lwork = (long)wkopt;
  WORK = malloc( lwork*sizeof(double) );
  if(!WORK){
    fprintf(stderr,"\nfailure allocating WORK in lapsolve\n\n");
    return -1;}

  /* Solve eigenproblem */
  dsyevx_(jobz, range, uplo, &nn, &MM[1][1], &lda, &vl, &vu, &lo_eig, &hi_eig, &abstol, &num_eigs, &VAL[1], Z, &nn, WORK, &lwork, IWORK, IFAIL, &info );	

  /* Put the solution into the output matrix */
  for(j=1; j<=num_eigs; j++)
    for(i=1; i<=nn; i++)
      VEC[i][j] = Z[nn*(j-1)+i-1];

  free(WORK);
  free(IWORK);
  free(IFAIL);
  free(Z);

  return 0;
}
