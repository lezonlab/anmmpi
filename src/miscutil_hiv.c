#include "miscutil_hiv.h"
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>


/* "SortSparseMatrixByIndex" sorts a dSparse_Matrix containing 'num_elements' elements by the index 'sort_index' */
/* Replaces 'dsort_PP2' */
int SortSparseMatrixByIndex(dSparse_Matrix *MM,int num_elements,int sort_index)
{
  double x;
  int i,ir,j,l,hi,i1,i2,ndx;
  long rra,*ra;

  if(num_elements < 2) return;
  if(sort_index<1 || sort_index>2){
    fprintf(stderr,"SortSparseMatrixByIndex: Bad index value %d\n\n",sort_index);
    return -1;}
  ndx = sort_index==1 ? 2 : 1;

  /* Create a vector to index the elements of MM */
  hi = 0;
  for(i=1; i<=num_elements; i++)
    if(MM->IDX[i][ndx]>hi) hi = MM->IDX[i][ndx];
  ra = lvector(1,num_elements);
  for(i=1; i<=num_elements; i++) ra[i] = (long)hi*(MM->IDX[i][sort_index]-1) + MM->IDX[i][ndx];

  /* Sort */
  l = (num_elements >> 1)+1;
  ir = num_elements;
  for(;;){
    if(l > 1){
      rra = ra[--l];
      i1 = MM->IDX[l][sort_index];
      i2 = MM->IDX[l][ndx];
      x = MM->X[l];
    }
    else {
      rra = ra[ir];
      i1 = MM->IDX[ir][sort_index];
      i2 = MM->IDX[ir][ndx];
      x = MM->X[ir];
      ra[ir] = ra[1];
      MM->IDX[ir][sort_index] = MM->IDX[1][sort_index];
      MM->IDX[ir][ndx] = MM->IDX[1][ndx];
      MM->X[ir] = MM->X[1];
      if (--ir == 1) {
	ra[1] = rra;
	MM->IDX[1][sort_index] = i1;
	MM->IDX[1][ndx] = i2;
	MM->X[1] = x;
	break;
      }
    }
    i = l;
    j = l + l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
	ra[i] = ra[j];
	MM->IDX[i][sort_index] = MM->IDX[j][sort_index];
	MM->IDX[i][ndx] = MM->IDX[j][ndx];
	MM->X[i] = MM->X[j];
	i = j;
	j <<= 1;
      } else j = ir + 1;
    }
    ra[i] = rra;
    MM->IDX[i][sort_index] = i1;
    MM->IDX[i][ndx] = i2;
    MM->X[i] = x;
  }
  free_lvector(ra,1,num_elements);

  /* Check for proper sorting */
  for(i=1; i<num_elements; i++){
    if(MM->IDX[i+1][sort_index] < MM->IDX[i][sort_index] || (MM->IDX[i+1][sort_index]==MM->IDX[i][sort_index] && MM->IDX[i+1][ndx] < MM->IDX[i][ndx])){
      fprintf(stderr,"\n\nImproper sort in SortSparseMatrixByIndex:\n%d:\t%d\t%d\n%d:\t%d\t%d\n\n",i,MM->IDX[i][sort_index],MM->IDX[i][ndx],i+1,MM->IDX[i+1][sort_index],MM->IDX[i+1][ndx]);
      return -1;
    }
  }
  return 0;
}


/* 'CopySparseMatrix" copies elements 'lo' through 'hi' of dSparse_Matrix 'source_matrix' to dSparse_Matrix 'target_matrix' */
/* Replaces "copy_dsparse" */
void CopySparseMatrix(dSparse_Matrix *source_matrix,dSparse_Matrix *target_matrix,int lo,int hi)
{
  int i;

  for(i=lo; i<=hi; i++){
    target_matrix->IDX[i][1] = source_matrix->IDX[i][1];
    target_matrix->IDX[i][2] = source_matrix->IDX[i][2];
    target_matrix->X[i] = source_matrix->X[i];
  }
}


/* "SortIntMatrixByCol" sorts the rows of the integer matrix RA[rlo..rhi][clo..chi] according to column 'idx'. */
/* WARNING: This isn't tested in its current format! */
void SortIntMatrixByCol(int **RA,int rlo,int rhi,int clo,int chi,int idx)
{
  int *RRA;
  int i,ir,j,l,k;
  int row=rhi-rlo+1;
  int col=chi-clo+1;

  if(row<2) return;
  RRA = ivector(1,col);
  l = (row >> 1) + 1;
  ir = row;
  for(;;){
    if(l > 1){
      l--;
      for(k=1; k<=col; k++) RRA[k] = RA[l-1+rlo][k-1+clo];
    } 
    else {
      for(k=1; k<=col; k++){
	RRA[k] = RA[ir-1+rlo][k-1+clo];
	RA[ir-1+rlo][k-1+clo] = RA[rlo][k-1+clo];
      }
      if (--ir == 1) {
	for(k=1; k<=col; k++) RA[rlo][k-1+clo] = RRA[k];
	break;
      }
    }
    i = l;
    j = l+l;
    while (j <= ir) {
      if (j < ir && RA[j-1+rlo][idx] < RA[j+rlo][idx]) j++;
      if (RRA[idx] < RA[j-1+rlo][idx]) {
	for(k=1; k<=col; k++) RA[i-1+rlo][k-1+clo] = RA[j-1+rlo][k-1+clo];
	i = j;
	j <<= 1;
      } else j = ir + 1;
    }
    for(k=1; k<=col; k++) RA[i-1+rlo][k-1+clo] = RRA[k];
  }
  free_ivector(RRA,1,col);
}


/* "AssignIndexStretches" populates the n-component vector 'stretches': Given the 'elm'-element 
   sparse matrix 'PP', sorted by index 'idx', the function populates 'stretches' such that 
   for all j: stretches[i] <= j < stretches[i+1], PP->IDX[j][idx]=i */
/* Replaces 'init_bst' */
void AssignIndexStretches(int *stretches,dSparse_Matrix *PP,int elm,int n,int idx)
{
  int i;

  for(i=1; i<n; i++) stretches[i] = 0;
  for(i=elm; i>0; i--) stretches[PP->IDX[i][idx]] = i;
  stretches[n] = elm + 1;
  for(i=n-1; i>0; i--)
    if(stretches[i]==0)
      stretches[i] = stretches[i+1];
}


/* allocate a char matrix with subscript range m[nrl..nrh][ncl..nch] */
char **cmatrix(long nrl,long nrh,long ncl,long nch)
{
  long nrow=nrh-nrl+1;
  long ncol=nch-ncl+2; /* +2 because the last element is '\0' */
  long i;
  static char **m;

  /* allocate pointers to rows */
  m = malloc(nrow*sizeof(char*));
  if(!m){
    fprintf(stderr,"\nallocation failure 1 in cmatrix\n\n");
    return NULL;}
  m-=nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = malloc(nrow*ncol*sizeof(char));
  if(!m[nrl]){
    fprintf(stderr,"\nallocation failure 2 in cmatrix\n\n");
    return NULL;}
  m[nrl]-=ncl;

  for(i=nrl+1; i<=nrh; i++) m[i] = m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}


/* free a char matrix allocated by cmatrix() */
void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch)
{
  free(m[nrl]+ncl);
  free(m+nrl);
}

/* allocate an int vector with subscript range v[nl..nh] */
int *ivector(long nl, long nh)
{
  int *v;

  v = malloc( (nh-nl+1)*sizeof(int) );
  if (!v) {
    fprintf(stderr,"\nAllocation failure in ivector\n\n");
    return NULL;}
  return v-nl;
}

/* free an int vector allocated with ivector() */
void free_ivector(int *v, long nl, long nh) {free(v+nl);}

/* allocate a long vector with subscript range v[nl..nh] */
long *lvector(long nl, long nh)
{
  long *v;

  v = malloc( (nh-nl+1)*sizeof(long) );
  if (!v) {
    fprintf(stderr,"\nAllocation failure in lvector\n\n");
    return NULL;}
  return v-nl;
}

/* free an int vector allocated with lvector() */
void free_lvector(long *v, long nl, long nh) {free(v+nl);}

/* allocate a double vector with subscript range v[nl..nh] */
double *dvector(long nl, long nh)
{
  double *v;

  v = malloc( (nh-nl+1)*sizeof(double) );
  if (!v) {
    fprintf(stderr,"\nAllocation failure in dvector\n\n");
    return NULL;}
  return v-nl;
}

/* free an int vector allocated with dvector() */
void free_dvector(double *v, long nl, long nh) {free(v+nl);}

/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
int **imatrix(long nrl, long nrh, long ncl, long nch)
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  int **m;

  /* allocate pointers to rows */
  m = malloc(nrow*sizeof(int*));
  if (!m){
    fprintf(stderr,"\nAllocation failure 1 in imatrix()\n\n");
    return NULL;}
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = malloc( nrow*ncol*sizeof(int) );
  if (!m[nrl]){
    fprintf(stderr,"\nAllocation failure 2 in imatrix()\n\n");
    return NULL;}
  m[nrl] -= ncl;

  for(i=nrl+1; i<=nrh; i++) m[i] = m[i-1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

/* free an int matrix allocated by imatrix() */
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
{
  free(m[nrl]+ncl);
  free(m+nrl);
}

/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /* allocate pointers to rows */
  m = malloc(nrow*sizeof(double*));
  if (!m){
    fprintf(stderr,"\nAllocation failure 1 in dmatrix()\n\n");
    return NULL;}
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = malloc( nrow*ncol*sizeof(double) );
  if (!m[nrl]){
    fprintf(stderr,"\nAllocation failure 2 in dmatrix()\n\n");
    return NULL;}
  m[nrl] -= ncl;

  for(i=nrl+1; i<=nrh; i++) m[i] = m[i-1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

/* free an int matrix allocated by imatrix() */
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
  free(m[nrl]+ncl);
  free(m+nrl);
}




/* ------------------ Numerical Recipes Routines ------------------ */
double dpythag(double a, double b)
{
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+DSQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+DSQR(absa/absb)));
}

void dsvdcmp(double **a, int m, int n, double w[], double **v)
{
  double dpythag(double a, double b);
  int flag,i,its,j,jj,k,l,nm;
  double anorm,c,f,g,h,s,scale,x,y,z,*rv1;
  static int maxits=100;

  rv1=dvector(1,n);
  g=scale=anorm=0.0;
  for (i=1;i<=n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m) {
      for (k=i;k<=m;k++) scale += fabs(a[k][i]);
      if (scale) {
	for (k=i;k<=m;k++) {
	  a[k][i] /= scale;
	  s += a[k][i]*a[k][i];
	}
	f=a[i][i];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][i]=f-g;
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
	  f=s/h;
	  for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
	}
	for (k=i;k<=m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i <= m && i != n) {
      for (k=l;k<=n;k++) scale += fabs(a[i][k]);
      if (scale) {
	for (k=l;k<=n;k++) {
	  a[i][k] /= scale;
	  s += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][l]=f-g;
	for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
	for (j=l;j<=m;j++) {
	  for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
	  for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
	}
	for (k=l;k<=n;k++) a[i][k] *= scale;
      }
    }
    anorm=DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n;i>=1;i--) {
    if (i < n) {
      if (g) {
	for (j=l;j<=n;j++) v[j][i]=(a[i][j]/a[i][l])/g;
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
	  for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
	}
      }
      for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=IMIN(m,n);i>=1;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j<=n;j++) a[i][j]=0.0;
    if (g) {
      g=1.0/g;
      for (j=l;j<=n;j++) {
	for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
	f=(s/a[i][i])*g;
	for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
      }
      for (j=i;j<=m;j++) a[j][i] *= g;
    } else for (j=i;j<=m;j++) a[j][i]=0.0;
    ++a[i][i];
  }
  for (k=n;k>=1;k--) {
    for (its=1;its<=maxits;its++) {
      flag=1;
      for (l=k;l>=1;l--) {
	nm=l-1;
	if ((double)(fabs(rv1[l])+anorm) == anorm) {
	  flag=0;
	  break;
	}
	if ((double)(fabs(w[nm])+anorm) == anorm) break;
      }
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<=k;i++) {
	  f=s*rv1[i];
	  rv1[i]=c*rv1[i];
	  if ((double)(fabs(f)+anorm) == anorm) break;
	  g=w[i];
	  h=dpythag(f,g);
	  w[i]=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=1;j<=m;j++) {
	    y=a[j][nm];
	    z=a[j][i];
	    a[j][nm]=y*c+z*s;
	    a[j][i]=z*c-y*s;
	  }
	}
      }
      z=w[k];
      if (l == k) {
	if (z < 0.0) {
	  w[k] = -z;
	  for (j=1;j<=n;j++) v[j][k] = -v[j][k];
	}
	break;
      }
      if (its == maxits) nrerror("no convergence in many dsvdcmp iterations");
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=dpythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=dpythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g = g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=1;jj<=n;jj++) {
	  x=v[jj][j];
	  z=v[jj][i];
	  v[jj][j]=x*c+z*s;
	  v[jj][i]=z*c-x*s;
	}
	z=dpythag(f,h);
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	for (jj=1;jj<=m;jj++) {
	  y=a[jj][j];
	  z=a[jj][i];
	  a[jj][j]=y*c+z*s;
	  a[jj][i]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  free_dvector(rv1,1,n);
}

void deigsrt(double d[], double **v, int n)
{
  int k,j,i;
  double p;
  
  for (i=1;i<n;i++) {
    p=d[k=i];
    for (j=i+1;j<=n;j++)
      if (d[j] >= p) p=d[k=j];
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
      for (j=1;j<=n;j++) {
	p=v[j][i];
	v[j][i]=v[j][k];
	v[j][k]=p;
      }
    }
  }
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}



/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
  long i,j;
  long nrow=nrh-nrl+1;
  long ncol=nch-ncl+1;
  long ndep=ndh-ndl+1;
  double ***t;

  /* allocate pointers to pointers to rows */
  t = malloc( nrow*sizeof(double**) );
  if (!t) nrerror("allocation failure 1 in d3tensor()");
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl] = malloc( nrow*ncol*sizeof(double*) );
  if (!t[nrl]) nrerror("allocation failure 2 in d3tensor()");
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl] = malloc( nrow*ncol*ndep*sizeof(double) );
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in d3tensor()");
  t[nrl][ncl] -= ndl;

  for(j=ncl+1; j<=nch; j++) t[nrl][j] = t[nrl][j-1] + ndep;
  for(i=nrl+1; i<=nrh; i++){
    t[i] = t[i-1] + ncol;
    t[i][ncl] = t[i-1][ncl] + ncol*ndep;
    for(j=ncl+1; j<=nch; j++) t[i][j] = t[i][j-1] + ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}

/* free a double d3tensor allocated by d3tensor() */
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,long ndl, long ndh)
/* free a double d3tensor allocated by d3tensor() */
{
  free(t[nrl][ncl]+ndl);
  free(t[nrl]+ncl);
  free(t+nrl);
}


/* "nintelemstr" returns the number of integer elements 
   that are specified by the input string */ 
int nintelemstr(char *text)
{
  char *str1,*str2,*stp1,*stp2,*boo,*dup;
  int tlo,thi,n=0;

  /* First copy the string */
  dup=(char *)calloc(strlen(text)+1,sizeof(char));
  strcpy(dup,text);

  /* Token-ize string using commas */
  str1=strtok_r(dup,",",&stp1);
  while(str1!=NULL){

    /* Token-ize the comma-tokenized strings using dashes */
    str2=strtok_r(str1,"-",&stp2);
    while(str2!=NULL){
      tlo=strtol(str2,&boo,10);
      str2=strtok_r(NULL,"-",&stp2);

      /* If there was a dash, add the whole range */
      if(str2!=NULL){
	thi=strtol(str2,&boo,10);
	n+=thi-tlo;	
      }

      /* Otherwise add one for the single element */
      else n++;
    }
    str1=strtok_r(NULL,",",&stp1);
  }
  free(dup);
  return n;
}

/* "hintelemstr" places the integer elements specified 
   by the input string into an integer array */ 
int hintelemstr(char *text,int *H,int nn)
{
  char *str1,*str2,*stp1,*stp2,*boo,*dup;
  int tlo,thi,i=0,j;

  /* First copy the string */
  dup=(char *)calloc(strlen(text)+1,sizeof(char));
  strcpy(dup,text);

  /* Token-ize string using commas */
  str1=strtok_r(dup,",",&stp1);
  while(str1!=NULL){

    /* Token-ize the comma-tokenized strings using dashes */
    str2=strtok_r(str1,"-",&stp2);
    while(str2!=NULL){
      tlo=strtol(str2,&boo,10);
      str2=strtok_r(NULL,"-",&stp2);

      /* If there was a dash, add the whole range */
      if(str2!=NULL){
	thi=strtol(str2,&boo,10);
	for(j=tlo;j<thi;j++)
	  H[++i]=j;
      }

      /* Otherwise add one for the single element */
      else 
	H[++i]=tlo;
    }
    str1=strtok_r(NULL,",",&stp1);
  }
  free(dup);
}

/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
double **submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch, long newrl, long newcl)
{
  long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
  double **m;

  /* allocate array of pointers to rows */
  m = malloc( nrow*sizeof(double*) );
  if(!m){ 
    fprintf(stderr,"\nAllocation failure in submatrix()\n\n");
    return NULL;}
  m -= newrl;

  /* set pointers to rows */
  for(i=oldrl,j=newrl; i<=oldrh; i++,j++) m[j] = a[i] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_submatrix(double **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
  free( b+nrl );
}
