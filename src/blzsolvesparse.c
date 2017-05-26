/*****************************************************************************/
/*                                                                           */
/*  C version of the FORTRAN BLZPACK driver drvgp1.c.  This code can be      */
/*  directly incorporated into a C program, enabling the efficient           */
/*  decomposition of large sparse matrices.  It requires the BLZPACK         */
/*  libraries, which--at the time that this file's creation--can be          */
/*  downloaded from http://crd.lbl.gov/~osni.                                */
/*                                                                           */
/*  The code diagonalizes a sparse symmetric matrix and returns the          */
/*  eigenvalues and eigenvectors in double-precision arrays.  All memory     */
/*  is dynamically allocated at run time, so the program will only use       */
/*  as much memory as it needs for the task at hand.                         */
/*                                                                           */
/*  Because both BLZPACK and the HSL libraries are FORTRAN, some measures    */
/*  had to be taken to get them to run in C.  To call a FORTRAN routine      */
/*  from C, all 2D arrays should be stored as zero-offset 1D arrays.  The    */
/*  routine name is given in lowercase and appended with an underscore,      */
/*  and everything but arrays are passed to functions as pointers.  Also,    */
/*  note that FORTRAN uses column-major matrix indexing and unit-offset.     */
/*                                                                           */
/*  Example:  If A is a matrix and B is a scalar,                            */
/*                                                                           */
/*  FORTRAN                                 C                                */
/*                                                                           */
/*  DOUBLE PRECISION A(N,M)                 double A[n*m];                   */
/*  INTEGER B                               int b;                           */
/*  CALL FOO(A,B)                           foo_(A,&b);                      */
/*  A(I,J)                                  A[n*(J-1)+I-1]                   */
/*                                                                           */
/*                                                                           */
/* ------------------------------------------------------------------------- */
/*                                                                           */
/*                       CREATED 02/16/11 BY TIM LEZON                       */
/*                                                                           */
/*                     Overhauled 01/20/2016 by Tim Lezon                    */
/*                                                                           */
/*****************************************************************************/
/* 
   TO COMPILE: 
   $ gfortran -c MA47.f
   $ gcc -c blzsolvesparse.c -O3
   $ gcc blzdecomp_free.o blzsolvesparse.o MA47.o -o blzdecomp -lblzpack -llapack -lblas 
*/
#include<stdio.h>  
#include<math.h>  
#include<stdlib.h>

/* Calculates the bottom 'nreig' eigenvalues and eigenvectors of the symmetric matrix 'A[0..(ne-1)]'. 
   The non-zero element 'A[i]' has indices given by 'IRN[i]' and 'JCN[i]'. */
int blzSolveSparse(double *A, int *IRN, int *JCN, int ne, double *VAL, double **VEC, int nn, int lo_eig, int hi_eig, int nvbset)
{
  int ok=0;
  int nreig = hi_eig;
  int num_eigs = hi_eig - lo_eig + 1;

  /* Parameters for BLZDRD */
  int leig,listor,lni,lrstor,ncuv;//nvbset;
  int maxiwi;
  int i,j,k;
  int nneig=0,nvopu=0,lflag=0;
  double sigma=0.0;
  int *ISTOR;
  double *RSTOR,*U,*V,*X,*EIG;

  /* Parameters used in MA47 */
  int *INFO;
  int *ICNTL,*IW1,*IW2,*KEEP;
  long maxiw1,maxa,maxne;
  double *RINFO,*CNTL,*W,*C;
  int ione=1;

  /*****************************************************************************/
  /*                                                                           */
  /*    Part 1: Allocate memory and initialize numerical parameters.           */
  /*                                                                           */
  /*****************************************************************************/
  

  /* --------------------- MA47 initialization --------------------- 
     Initiate MA47 parameter values.  'maxa' is the size of array 'A'; 
     'maxne' is the maximum number of unique entries in a symmetric 
     nn-by-nn matrix; 'maxiw1' is the size of IW1, which is used as 
     workspace in MA47AD and needs to be 1.2*(2*ne +5*nn +4) or larger. */
  lni = nn;                                    /* First dimension of arrays U, V and X. lni >= nn */
  maxa = (long)2*nn*(nn+1);
  maxne = (long)maxa/4;
  maxiw1 = (long)2*(2*ne + 5*nn + 4);
  CNTL = calloc(2,sizeof(double));             /* Used in MA47ID, MA47BD */
  ICNTL = calloc(7,sizeof(int));               /* Used in MA47ID, MA47AD, MA47BD, MA47CD */
  IW1 = calloc(maxiw1,sizeof(int));            /* Used in MA47AD, MA47BD */
  KEEP = calloc(maxne+lni*5+2,sizeof(int));    /* Used in MA47AD, MA47BD */
  RINFO = calloc(4,sizeof(double));            /* Used in MA47AD, MA47BD */
  INFO = calloc(24,sizeof(int));               /* Used in MA47AD, MA47BD */
  IW2 = calloc(2*lni+2,sizeof(int));           /* Used in MA47BD, MA47CD */
  W = calloc(nn,sizeof(double));               /* Used in MA47CD */

  /* Set default values for MA47 routines */
  while(ok==0){
    ok = 1;
    ma47id_(CNTL,ICNTL);
    ma47ad_(&nn, &ne, IRN, JCN, IW1, &maxiw1, KEEP, ICNTL, RINFO, INFO);
    if(INFO[0]!=0){
      fprintf(stderr,"\n*ERROR: MA47AD, INFO(1) = %d\n\n",INFO[0]);
      return -1;}

    /* Increase the sizes of some arrays empirically to avoid errors in ma47bd:
       maxa >= INFO(6) after ma47ad (parameter LA in MA47BD)
       maxiw1 >= INFO(7) after ma47ad (parameter LIW in MA47BD) */
    if(INFO[5] > maxa){
      ok = 0;
      maxa = INFO[5];
      maxne = (long)maxa/4;
      free(KEEP);
      KEEP = calloc(maxne+lni*5+2,sizeof(int));
    }
    if(INFO[6] > maxiw1){
      ok = 0;
      maxiw1 = INFO[6];
      free(IW1);
      IW1 = calloc(maxiw1,sizeof(int));
    }
  }
  C = calloc(maxa,sizeof(double));
  /* ------------------- End of MA47 initialization ------------------ */



  /* ------------------- Calculate BLZPACK parameter values ------------------ */
  leig = 2*nreig+10;                /* 1st dimension of EIG & 2nd dimension of X. leig>=min(nreig+10,2*nreig) */
  //nvbset = 6;                       /* Number of vectors in a block. 6=blocks, 0=auto */
  ncuv = nvbset==0 ? 3 : nvbset;    /* 2nd dimension of arrays U and V. */
  listor = 2400;
  lrstor = 40*nn;
  U = calloc(lni*ncuv,sizeof(double));
  V = calloc(lni*ncuv,sizeof(double));        /* Used in MA47CD */
  X = calloc(lni*leig,sizeof(double));
  EIG = calloc(2*leig,sizeof(double));
  RSTOR = calloc(lrstor+5,sizeof(double));
  ISTOR = calloc(listor+17,sizeof(int));
  if(!RSTOR || !ISTOR || !U || !V || !X || !EIG || !INFO || !IW1 || 
     !IW2 || !KEEP || !RINFO || !W || !CNTL || !ICNTL || !C){
    fprintf(stderr,"\n * ERROR * Allocation failure in blzSolveSparse\n\n");
    fprintf(stderr,"nn = %d\nlni = %d\nleig = %d\nlistor = %d\nlrstor = %d\nmaxa = %ld\nmaxne = %ld\nmaxiw1 = %ld\n",nn,lni,leig,listor,lrstor,maxa,maxne,maxiw1);
    return -1;}
  

  /* Call blzdrd once to check the size of ISTOR and RSTOR */
  ISTOR[0] = nn;
  ISTOR[1] = lni;
  ISTOR[2] = nreig;
  ISTOR[3] = leig;
  ISTOR[4] = nvbset;
  ISTOR[8] = ISTOR[9] = ISTOR[10] = 1;
  ISTOR[12] = 6;

  blzdrd_(ISTOR, RSTOR, &sigma, &nneig, U, V, &lflag, &nvopu, EIG, X);

  if(ISTOR[14] > listor){
    listor = ISTOR[14];
    free(ISTOR);
    ISTOR = calloc(listor+17,sizeof(int)); 
    if(!ISTOR){
      fprintf(stderr,"\n * ERROR * ISTOR Allocation failure in blzSolveSparse\n\n");
      return -1;}
  }
  if((int)RSTOR[3] > lrstor){
    lrstor = (int)RSTOR[3];
    free(RSTOR);
    RSTOR = calloc(lrstor+5,sizeof(double));
    if(!RSTOR){
      fprintf(stderr,"\n * ERROR * RSTOR Allocation failure in blzSolveSparse\n\n");
      return -1;}
  }

  /* Assign values to BLZPACK inputs 'ISTOR' and 'RSTOR' */
  ISTOR[0] = nn;       /* Number of active rows of U, V and X. Equal to dimension of A in sequential mode. */
  ISTOR[1] = lni;
  ISTOR[2] = nreig;    /* Number of required eigenvalues and eigenvectors. */
  ISTOR[3] = leig;
  ISTOR[4] = nvbset;
  ISTOR[5] = 0;        /* Maximum number of steps per run (0=auto) */
  ISTOR[6] = 0;        /* Number of input starting vectors in V */
  ISTOR[7] = 0;        /* Number of input eigenpairs in EIG and X */
  ISTOR[8] = 1;        /* Problem type (0=standard,1=generalized) */
  ISTOR[9] = 1;        /* Spectrum slicing */
  ISTOR[10] = 1;       /* Purify (Correct for singularities in B) */
  ISTOR[11] = 0;       /* Print level: 0=Nothing; 1=errors & warnings; 2=1 + eigenvalues; 
			  3=2 + statistics; 4=3 + run-by-run; 5=4 + eigenvectors; 6=5+step-by-step */
  ISTOR[12] = 6;       /* Output file (0 || 6 = stdout; fort.k for all other k:0<k<100) */
  ISTOR[13] = 0;       /* LCOMM: Required for parallel mode. */
  ISTOR[14] = listor;  /* Size of ISTOR */
  RSTOR[0] = 0.0;      /* Lower limit of interval (Ignored for standard problem) */
  RSTOR[1] = 0.0;      /* Upper limit of interval (Ignored for standard problem) */
  RSTOR[2] = 0.0;      /* Convergence threshold */
  RSTOR[3] = (double)lrstor;   /* Size of RSTOR */



  /*****************************************************************************/
  /*                                                                           */
  /*    Part 2: Decompose                                                      */
  /*                                                                           */
  /*****************************************************************************/

  /* ------ *** This is where the work gets done *** ----- */
  lflag=0;
  do{
    blzdrd_(ISTOR, RSTOR, &sigma, &nneig, U, V, &lflag, &nvopu, EIG, X);

    if(lflag<0){
      fprintf(stderr,"\n *** ERROR: Abnormal exit of 'blzdrd_' with error code %d\n",ISTOR[15]);
      return -1;}

    /* Given U, solve C*V=U */
    else if(lflag==1){
      for(i=0; i<nvopu; i++){
	dcopy_(&nn, &U[lni*i], &ione, &V[lni*i], &ione);
	ma47cd_(&nn, C, &maxa, IW1, &maxiw1, W, &V[lni*i], IW2, ICNTL);
      }
    }


    /* Given U, find V=B*U */
    else if(lflag==2){
      for(i=0;i<nvopu;i++)
	for(j=0;j<nn;j++)
	  V[lni*i+j] = U[lni*i+j];
    }

    /* Given sigma, find C = A - sigma*B */
    else if(lflag==3){
      for(i=0; i<ne; i++){ 
	C[i] = A[i];
	if(IRN[i]==JCN[i])
	  C[i] -= sigma;
      }

      /* Factor C = LDL' */
      ma47bd_(&nn, &ne, JCN, C, &maxa, IW1, &maxiw1, KEEP, CNTL, ICNTL, IW2, RINFO, INFO);

      /* These are optional error messages for the MA47BD routine */
      if(INFO[0]!=0){
	fprintf(stderr,"\n* Error: ma47bd, INFO[0]!=0\n\n");
	return -1;}
      if(INFO[23]>0){
	fprintf(stderr,"\n* Error: ma47bd, INFO[23]>0\n\n");
	return -1;}
      nneig=INFO[22];
    }

  }while(lflag!=0);
  /* ----------- ***** End of Diagonalization Block ***** ---------- */

  /* Print BLZPACK errors */
  if(ISTOR[15]!=0)
    fprintf(stderr,"ERRORS: %d\n",ISTOR[15]);
  if(ISTOR[16]>=8192)
    fprintf(stderr,"** INCOMPLETE BASIS!! ***\n\nWARNINGS: %d\n\n",ISTOR[16]);


  /* Store eigenvalues in 'VAL' and eigenvectors in 'VEC' */
  for(j=0; j<num_eigs; j++){
    VAL[j+1] = EIG[ lo_eig+j-1 ];
    for(i=0; i<nn; i++)
      VEC[i+1][j+1] = X[ lni*(lo_eig+j-1) + i ];
  }
  /*
  for(j=0; j<nreig; j++){
    VAL[j+1] = EIG[j];
    for(i=0; i<nn; i++)
      VEC[i+1][j+1] = X[lni*j+i];
  }
  */

  /* Free memory */
  free(RSTOR);
  free(ISTOR);
  free(U);
  free(V);
  free(X);
  free(EIG);
  free(INFO);
  free(IW2);
  free(KEEP);
  free(RINFO);
  free(W);
  free(CNTL);
  free(ICNTL);
  free(C);
  free(IW1);

  return 0;
}
