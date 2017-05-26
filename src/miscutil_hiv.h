#ifndef MISCUTIL_HIV_H
#define MISCUTIL_HIV_H
#include<math.h>

/* ------------------ Numerical Recipes Macros ------------------ */

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

typedef struct {int **IDX; double *X; } dSparse_Matrix;

void SortIntMatrixByCol(int **RA,int rlo,int rhi,int clo,int chi,int idx);
int SortSparseMatrixByIndex(dSparse_Matrix *MM,int num_elements,int sort_index);
void CopySparseMatrix(dSparse_Matrix *source_matrix,dSparse_Matrix *target_matrix,int lo,int hi);
void AssignIndexStretches(int *stretches,dSparse_Matrix *PP,int elm,int n,int idx);

char **cmatrix(long nrl,long nrh,long ncl,long nch);
void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch);
int *ivector(long nl, long nh);
void free_ivector(int *v, long nl, long nh);
long *lvector(long nl, long nh);
void free_lvector(long *v, long nl, long nh);
double *dvector(long nl, long nh);
void free_dvector(double *v, long nl, long nh);
int **imatrix(long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,long ndl, long ndh);
int nintelemstr(char *text);
int hintelemstr(char *text,int *H,int nn);
double **submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch, long newrl, long newcl);
void free_submatrix(double **b, long nrl, long nrh, long ncl, long nch);

/* ------------------ Numerical Recipes Routines ------------------ */
double dpythag(double a, double b);
void dsvdcmp(double **a, int m, int n, double w[], double **v);
void deigsrt(double d[], double **v, int n);
void nrerror(char error_text[]);

#endif
