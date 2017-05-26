#include "rtbutil_hiv.h"
#include "miscutil_hiv.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* "ReadBlockfile" reads information from a block file. */
int ReadBlockfile(const char *file,char ***pdbfile,int *num_pdbs,Rigid_Block **blocks,long *num_blocks)
{
  FILE *data;
  Rigid_Block *temp_blocks;
  char buffer[RTB_BUFF_LENGTH];
  char **pdbid;
  char **ptmp,header[128],line[128];
  int nwd,i,j,k;

  /* 
     -- Skip over leading whitespace
     -- Handle everything after '#' as a comment
     -- Look for keywords
     -- Read data
  */

  /* 
     TODO: Allow wildcards so blocks don't have to be specified for every PDB
  */

  if( (data=fopen(file,"r") )==NULL ){
    printf("\nReadBlockfile: Unable to open %s\n\n",file);
    return -1;}

  /* -------- First pass counts number of PDB files and blocks -------- */
  *num_pdbs = *num_blocks = 0;
  while(fgets(buffer,RTB_BUFF_LENGTH,data)!=NULL){
    /* Skip leading whitespace and read the heading */
    i = j = 0;
    while(isspace(buffer[i])) i++;
    while(!isspace(buffer[i]) && buffer[i]!='\0') header[j++] = buffer[i++];
    header[j] = '\0';
    if(!strncmp(header,"PDB",strlen(header)) && strlen(header)==3) (*num_pdbs)++;
    else if(!strncmp(header,"BLOCK",strlen(header)) && strlen(header)==5) (*num_blocks)++;
  }

  /* Allocate memory */
  ptmp = cmatrix(1,(*num_pdbs),0,RTB_WORD_LENGTH);
  pdbid = cmatrix(1,*num_pdbs,0,RTB_WORD_LENGTH);
  temp_blocks = malloc( ((*num_blocks)+1)*sizeof(Rigid_Block) );
  rewind(data);

  /* -------- Second pass reads the stuff -------- */
  *num_pdbs = *num_blocks = 0;
  while(fgets(buffer,RTB_BUFF_LENGTH,data)!=NULL){
    i = j = 0;
    while(isspace(buffer[i])) i++;
    while(!isspace(buffer[i]) && buffer[i]!='\0')  header[j++] = buffer[i++];
    header[j]='\0';
    if(!strncmp(header,"PDB",strlen(header)) && strlen(header)==3){
      (*num_pdbs)++;
      while(isspace(buffer[i])) i++;
      j=0;
      while(!isspace(buffer[i]) && buffer[i]!='\0') ptmp[(*num_pdbs)][j++] = buffer[i++];
      ptmp[(*num_pdbs)][j]='\0';
      while(isspace(buffer[i])) i++;
      j=0;
      while(!isspace(buffer[i]) && buffer[i]!='\0') pdbid[*num_pdbs][j++] = buffer[i++];
      pdbid[*num_pdbs][j]='\0';
    }


    /* Blocks */
    /* BLOCK blknum [pdbid] lores lochain lonum hires hichain hinum */
    if(!strncmp(header,"BLOCK",strlen(header)) && strlen(header)==5){
      (*num_blocks)++;

      /* Count the number of words in the block */
      /* TODO: Put this into its own function! */
      while(isspace(buffer[i])) i++;
      j=i;
      nwd=0;
      while(buffer[j]!='\0'){
	if(!isspace(buffer[j]) && isspace(buffer[j-1])) nwd++;
	j++;
      }

      /* Block number */
      j=0;
      while(!isspace(buffer[i]) && buffer[i]!='\0') line[j++]=buffer[i++];
      line[j]='\0';
      sscanf(line,"%d",&temp_blocks[*num_blocks].blknum);
      
      /* PDB id */
      if(nwd!=8 && nwd!=7){
	printf("\nReadBlockfile: Bad line:\n\n%s\n",buffer);
	return -1;}
      else if(nwd==8){
	j=0;
	while(isspace(buffer[i])) i++;
	while(!isspace(buffer[i]) && buffer[i]!='\0') line[j++] = buffer[i++];
	line[j]='\0';
	for(k=1; k<=(*num_pdbs); k++)
	  if(!strncmp(line,pdbid[k],strlen(line))) temp_blocks[*num_blocks].pdbid = k;
      }
      else temp_blocks[*num_blocks].pdbid = -1;

      /* lores, lochain, lonum */
      j=0;
      while(isspace(buffer[i])) i++;
      while(!isspace(buffer[i]) && buffer[i]!='\0' && j<4) line[j++] = buffer[i++];
      line[j]='\0';
      //strncpy(temp_blocks[*num_blocks].LORES,line,3);
      strcpy(temp_blocks[*num_blocks].LORES,line);
      j=0;
      while(isspace(buffer[i])) i++;
      while(!isspace(buffer[i]) && buffer[i]!='\0') line[j++] = buffer[i++];
      line[j]='\0';
      temp_blocks[*num_blocks].lochain = line[0];
      j=0;
      while(isspace(buffer[i])) i++;
      while(!isspace(buffer[i]) && buffer[i]!='\0') line[j++] = buffer[i++];
      line[j]='\0';
      sscanf(line, "%d", &temp_blocks[*num_blocks].lonum);

      /* hires, hichain, hinum */
      j=0;
      while(isspace(buffer[i])) i++;
      while(!isspace(buffer[i]) && buffer[i]!='\0' && j<4) line[j++] = buffer[i++];
      line[j]='\0';
      //strncpy(temp_blocks[*num_blocks].HIRES,line,3);
      strcpy(temp_blocks[*num_blocks].HIRES,line);
      j=0;
      while(isspace(buffer[i])) i++;
      while(!isspace(buffer[i]) && buffer[i]!='\0') line[j++] = buffer[i++];
      line[j]='\0';
      temp_blocks[*num_blocks].hichain = line[0];
      j=0;
      while(isspace(buffer[i])) i++;
      while(!isspace(buffer[i]) && buffer[i]!='\0') line[j++] = buffer[i++];
      line[j]='\0';
      sscanf(line,"%d",&temp_blocks[*num_blocks].hinum);
    }
  }

  *pdbfile = ptmp;
  *blocks = temp_blocks;

  free_cmatrix(pdbid, 1, (*num_pdbs), 0, RTB_WORD_LENGTH);
  fclose(data);

  return 0;
}


/* "AssignRigidBlocks" assigns each residue to a rigid block and returns the total number of rigid blocks. 
   'pdb' points to an array of 'num_pdbs' PDB_File structures. The array 'blocks' contains 'num_input_blocks' 
   block definitions to be assigned to the residues in the PDB file. Some of these may be overlapping. 
   On return, 'max_block_size' contains the number of residues in the largest block found. Block IDs are 
   assigned beginning with 1 and extending as high as necessary. Unspecified residues are assigned 
   to block 0. */
int AssignRigidBlocks(PDB_File *pdb, Rigid_Block *blocks, int num_pdbs, long num_input_blocks, long *max_block_size)
{
  char ch;
  int *unique_blocks;
  int **blockmap;
  int num_unique_blocks = 0;
  int high_block_index = 0;
  int kwn = 0;
  int ok,rr,i,j,k,p;

  /* Make a list of unique block IDs */
  unique_blocks = ivector(1, num_input_blocks);
  for(i=1; i<=num_input_blocks; i++) unique_blocks[i] = 0;
  for(i=1; i<=num_input_blocks; i++){
    ok = 1;
    for(j=1; j<=num_unique_blocks; j++)
      if(blocks[i].blknum==unique_blocks[j]){
	ok=0;
	break;
      }
    if(ok==1) unique_blocks[++num_unique_blocks] = blocks[i].blknum;
  }

  /* Assigns each residue to a block:  blockmap[i][1] contains the integer ID of 
     block i as provided in the .blk file.  blockmap[i][2] contains the sequential 
     ID of block i that will be used internally. */
  blockmap = imatrix(1, num_unique_blocks, 1, 2);
  for(i=1; i<=num_unique_blocks; i++) blockmap[i][1] = blockmap[i][2] = 0;
  for(p=1; p<=num_pdbs; p++){
    for(i=1; i<=pdb[p].num_residues; i++){
      rr = pdb[p].atom[i].resnum;
      ch = pdb[p].atom[i].chain;
      ok = 0;
      for(j=1; j<=num_input_blocks; j++){
	if( blocks[j].pdbid==p && 
	    (
	     (blocks[j].lochain<ch && ch<blocks[j].hichain) ||
	     (blocks[j].lochain==ch && ch<blocks[j].hichain && rr>=blocks[j].lonum) ||
	     (blocks[j].lochain<ch && ch==blocks[j].hichain && rr<=blocks[j].hinum) ||
	     (blocks[j].lochain==ch && blocks[j].hichain==ch && 
	      rr>=blocks[j].lonum && rr<=blocks[j].hinum)
	     )
	    ){
	  for(k=1; k<=kwn; k++)
	    if(blocks[j].blknum==blockmap[k][1]){
	      ok=1;
	      pdb[p].atom[i].model = blockmap[k][2];
	      break;
	    }
	  if(ok==0){
	    blockmap[++kwn][1] = blocks[j].blknum;
	    blockmap[kwn][2] = pdb[p].atom[i].model = ++high_block_index;
	    ok = 1;
	  }
	  break;
	}
      }
      if(ok==0) pdb[p].atom[i].model = 0;
    }
  }
  free_ivector(unique_blocks, 1, num_input_blocks);
  free_imatrix(blockmap, 1, num_unique_blocks, 1, num_unique_blocks);


  /* Find the size of the largest block */
  unique_blocks = ivector(1, high_block_index);
  for(i=1; i<=high_block_index; i++) unique_blocks[i] = 0;
  for(p=1; p<=num_pdbs; p++)
    for(i=1; i<=pdb[p].num_residues; i++) 
      if(pdb[p].atom[i].model!=0)
	unique_blocks[pdb[p].atom[i].model]++;
  (*max_block_size) = 0;
  for(i=1; i<=high_block_index; i++)
    if(unique_blocks[i]>(*max_block_size)) *max_block_size = unique_blocks[i];
  free_ivector(unique_blocks, 1, high_block_index);

  return high_block_index;
}



/* "CalcProjectionMatrix" calculates the projection matrix from the full residue space to rigid block space. 
   It stores the projection matrix in sparse format in 'projection_matrix' and stores the centers of mass of 
   the blocks in the matrix 'block_center_of_mass'. Returns the number of elements in 'projection_matrix'. */
/* NOTE: Replaces 'dblock_projections3' */
int CalcProjectionMatrix(dSparse_Matrix *projection_matrix,PDB_File *pdb,double **block_center_of_mass,int num_pdbs,long num_blocks,long max_block_size)
{
  double **X,**I,**IC,*W,**A,**ISQT;
  double x,tr,dd,df;
  int runtot=0;
  int *IDX,nbp,b,i,j,k,q,ii,jj,aa,bb,num_projection_elements;

  /* Initialize local arrays */
  num_projection_elements = 0;
  X = dmatrix(1,max_block_size,1,3);
  IDX = ivector(1,max_block_size);
  I = dmatrix(1,3,1,3);
  IC = dmatrix(1,3,1,3);
  W = dvector(1,3);
  A = dmatrix(1,3,1,3);
  ISQT = dmatrix(1,3,1,3);


  /* NOTE:  This is an inefficient (num_blocks*num_residues) loop! */ 
  /* Cycle through blocks */
  for(b=1; b<=num_blocks; b++){

    /* Clear matrices */
    for(j=1; j<=3; j++){
      block_center_of_mass[b][j] = 0.0;
      for(i=1; i<=3; i++) I[i][j] = 0.0;
      for(i=1; i<=max_block_size; i++) X[i][j] = 0.0;
    }

    /* Store values for current block.  'runtot' keeps a running total of 
       the number of residues (1...tres) */
    nbp = runtot = 0;
    for(q=1; q<=num_pdbs; q++){ 
      if(q>1) runtot += pdb[q-1].num_residues;
      for(i=1; i<=pdb[q].num_residues; i++){
	if(pdb[q].atom[i].model==b){
	  IDX[++nbp] = runtot+i;
	  for(j=1; j<=3; j++){
	    x = (double)pdb[q].atom[i].X[j-1];
	    X[nbp][j] = x;
	    block_center_of_mass[b][j] += x;
	  }   
	}
      }
    }

    /* Translate block CM to origin */
    for(j=1; j<=3; j++) block_center_of_mass[b][j]/=(double)nbp;
    for(i=1; i<=nbp; i++)
      for(j=1; j<=3; j++)
	X[i][j] -= block_center_of_mass[b][j];

    /* Calculate inertia tensor */
    for(k=1; k<=nbp; k++){
      dd=0.0;
      for(j=1; j<=3; j++){
	df = X[k][j];
	dd += df*df;
      }
      for(i=1; i<=3; i++){
	I[i][i] += (dd-X[k][i]*X[k][i]);
	for(j=i+1; j<=3; j++){
	  I[i][j] -= X[k][i]*X[k][j];
	  I[j][i] = I[i][j];
	}
      }
    }

    /* Diagonalize inertia tensor */
    for(i=1; i<=3; i++)
      for(j=1; j<=3; j++)
	IC[i][j] = I[i][j];
    dsvdcmp(IC,3,3,W,A);
    /* NOTE: Must this be a right-handed coordinate system? */

    /* Find its square root */
    for(i=1; i<=3; i++)
      for(j=1; j<=3; j++){
	dd=0.0;
	for(k=1; k<=3; k++)
	  dd += A[i][k]*A[j][k]/sqrt(W[k]);
	ISQT[i][j] = dd;
      }

    /* Update projection_matrix with the rigid motions of the block */
    tr = 1.0/sqrt((double)nbp);
    for(i=1; i<=nbp; i++){

      /* Translations: 3*(IDX[i]-1)+1 = x-coordinate of residue IDX[i];
	 6*(b-1)+1 = x-coordinate of block b */
      for(j=1; j<=3; j++){
	num_projection_elements++;
	projection_matrix->IDX[num_projection_elements][1] = 3*(IDX[i]-1)+j;
	projection_matrix->IDX[num_projection_elements][2] = 6*(b-1)+j;
	projection_matrix->X[num_projection_elements] = tr;
      }

      /* Rotations */
      if(nbp>1){
	for(ii=1; ii<=3; ii++){
	  for(jj=1; jj<=3; jj++){
	    if(jj==1) {aa=2; bb=3;}
	    else if(jj==2) {aa=3; bb=1;}
	    else {aa=1; bb=2;}
	    dd = ISQT[ii][aa]*X[i][bb] - ISQT[ii][bb]*X[i][aa];
	    num_projection_elements++;
	    projection_matrix->IDX[num_projection_elements][1] = 3*(IDX[i]-1)+jj;
	    projection_matrix->IDX[num_projection_elements][2] = 6*(b-1)+3+ii;
	    projection_matrix->X[num_projection_elements] = dd;
	  }
	}
      }
    }
  }
  free_dmatrix(X,1,max_block_size,1,3);
  free_ivector(IDX,1,max_block_size);
  free_dmatrix(I,1,3,1,3);
  free_dmatrix(IC,1,3,1,3);
  free_dvector(W,1,3);
  free_dmatrix(A,1,3,1,3);
  free_dmatrix(ISQT,1,3,1,3);

  return num_projection_elements;
}

/* "CondensePrjMtx" condenses a projection matrix by offsetting the column 
   elements if there are blocks with fewer than six degrees of freedom. */
void CondensePrjMtx(dSparse_Matrix *projection_matrix, long num_projection_elements)
{
  int *I1,*I2,max=0,i,j=0;

  for(i=1; i<=num_projection_elements; i++)
    if(projection_matrix->IDX[i][2] > max)
      max = projection_matrix->IDX[i][2];
  I1 = ivector(1,max);
  I2 = ivector(1,max);
  for(i=1; i<=max; i++) I1[i] = 0;
  for(i=1; i<=num_projection_elements; i++) I1[projection_matrix->IDX[i][2]] = projection_matrix->IDX[i][2];
  for(i=1; i<=max; i++){
    if(I1[i]!=0) j++;
    I2[i]=j;
  }
  for(i=1; i<=num_projection_elements; i++)
    if(projection_matrix->X[i] != 0.0)
      projection_matrix->IDX[i][2] = I2[projection_matrix->IDX[i][2]];
  free_ivector(I1,1,max);
  free_ivector(I2,1,max);
}


/* "PrintProjectionMatrix" creates a copy of the input projection matrix, condenses it to eliminate empty
   columns, prints it and destroys it. The input matrix should remain unchanged. */
int PrintProjectionMatrix(char *outfile, dSparse_Matrix *projection_matrix, long num_projection_elements, long total_residues)
{
  FILE *data;
  dSparse_Matrix *proj_copy;
  int i;

  proj_copy = malloc(sizeof(dSparse_Matrix));
  proj_copy->IDX = imatrix(1, num_projection_elements, 1, 2);
  proj_copy->X = dvector(1, num_projection_elements);
  CopySparseMatrix(projection_matrix, proj_copy, 1, num_projection_elements);

  if(num_projection_elements < 12*total_residues){
    //printf("Condensing projection matrix to eliminate %ld columns\n", 
    //	      (12*total_residues-num_projection_elements)/3);
    CondensePrjMtx(proj_copy, num_projection_elements); }
  if((data=fopen(outfile,"w"))==NULL){
    fprintf(stderr,"\nUnable to open %s\n\n",outfile);
    return -1;}
  printf("Writing projection matrix to %s\n",outfile);
  for(i=1; i<=num_projection_elements; i++) 
    fprintf(data,"%8d%8d% 25.15e\n", proj_copy->IDX[i][1], proj_copy->IDX[i][2], proj_copy->X[i]);
  fclose(data);

  free_imatrix(proj_copy->IDX, 1, num_projection_elements, 1, 2);
  free_dvector(proj_copy->X, 1, num_projection_elements);
  free(proj_copy);
  return 0;
}


/* "FindBlockContactsCypa" generates a contact map for rigid blocks, taking CypA binding into account. Block 
   definitions are taken from 'model' values for atoms. On exit, each nonzero element 'block_contacts[i][j]' 
   contains a unique index for the contact between blocks i and j. If blocks are not in contact, 
   block_contacts[i][j]=0. The mass centers of the blocks are stored in 'block_center_of_mass'. In addition 
   to finding contacts between blocks with residues within 'g_anm_cutoff' distance from each other, this function 
   assigns contacts to blocks joined by CypA binding. The 'num_blocks'-dimensional square matrix 'cypa_contacts' 
   contains these contacts across blocks. The function returns the total number of block-block contacts, 
   which is also the largest element of 'block_contacts'. */  
int FindBlockContactsCypa(int **block_contacts, PDB_File *PDB, int **blockmap, int *blockmap_index, 
			  double **block_center_of_mass, int num_pdbs, int num_blocks, 
			  int **chain_starts, int **cypa_contacts, int num_chains)
{
  double *DMAX;
  double cutsq=g_anm_cutoff*g_anm_cutoff;
  int nc,i,j,k,p;
  int blk;
  int a,b,p1,p2,r1,r2,ok;
  double df,dd;

  /* Get the maximum radius of each block */
  DMAX = dvector(1,num_blocks);
  for(i=1; i<=num_blocks; i++) DMAX[i] = 0.0;
  j = 0;
  for(p=1; p<=num_pdbs; p++){
    for(i=1; i<=PDB[p].num_residues; i++){
      dd = 0.0;
      blk = PDB[p].atom[i].model;
      for(k=1; k<=3; k++){
	df = PDB[p].atom[i].X[k-1] - block_center_of_mass[blk][k];
	dd += df*df;
      }
      if(dd>DMAX[blk]) DMAX[blk] = dd;
    }
  }
  for(i=1; i<=num_blocks; i++) DMAX[i] = sqrt(DMAX[i]);

  /* Tentatively join blocks that could possibly touch: 
     Distance between mass centers is less than sum of radii + cutoff */
  for(i=1; i<=num_blocks; i++){
    block_contacts[i][i] = 1;
    for(j=i+1; j<=num_blocks; j++){
      dd = 0.0;
      for(k=1; k<=3; k++){
	df = block_center_of_mass[i][k] - block_center_of_mass[j][k];
	dd += df*df;
      }
      if( sqrt(dd) <= DMAX[i] + DMAX[j] + g_anm_cutoff ) block_contacts[i][j] = block_contacts[j][i] = 1;
      else block_contacts[i][j] = block_contacts[j][i] = 0;
    }
  }
  free_dvector(DMAX,1,num_blocks);

  /* Refine matrix by explicitly checking for contacts */
  for(a=1; a<=num_blocks; a++){
    for(b=a+1; b<=num_blocks; b++){
      if(cypa_contacts[a][b]!=0) block_contacts[a][b] = block_contacts[b][a] = 1;
      else if(block_contacts[a][b]!=0){
	ok=0;
	for(i=blockmap_index[a]; i<blockmap_index[a+1]; i++){
	  p1 = blockmap[i][2];
	  r1 = blockmap[i][3];
	  for(j=blockmap_index[b]; j<blockmap_index[b+1]; j++){
	    p2 = blockmap[j][2];
	    r2 = blockmap[j][3];
	    dd = 0.0;
	    for(k=0; k<3; k++){
	      df = PDB[p1].atom[r1].X[k] - PDB[p2].atom[r2].X[k];
	      dd+=df*df;
	    }
	    if(dd<=cutsq){
	      ok=1;
	      break;
	    }
	  }
	  if(ok==1) break;
	}
	if(ok==0) block_contacts[a][b] = block_contacts[b][a] = 0;
      }
    }
  }

  /* Give a unique index to each contacting pair */
  nc=0;
  for(i=1; i<=num_blocks; i++)
    for(j=i; j<=num_blocks; j++)
      if(block_contacts[i][j]!=0) block_contacts[i][j] = block_contacts[j][i] = ++nc;
  return nc;
}


/* "AddFullHessRowsToBlockHessCypa" updates the array HT with products of the projection matrix and a subset 
   of rows of the full Hessian matrix. 'HT[1..num_block_contacts][1..6][1..6]' is an array of 6-by-6 arrays, 
   each corresponding to interations between two blocks, 'hessian_superrow[1..3*num_residues][1..3]' is the 
   transpose of a Hessian superrow, etc. The superrows corresponding to residues in the range [lores,hires] 
   is multiplied by the projection matrix, and the 'nko' residues with resids in 'KO' are simultaneously 
   perturbed. Finally, the chains connected by CypA are joined: 'cypa_contacts[1..num_blocks][1..num_blocks]' 
   contains the contact map between blocks connected with CypA. */
void AddFullHessRowsToBlockHessCypa(double ***HT, int num_block_contacts, double **hessian_superrow, 
				    int **block_contacts, dSparse_Matrix *projections_sorted_by_1,
				    dSparse_Matrix *projections_sorted_by_2, int num_projection_elements,
				    int **blockmap, int *projection_index_1, int *projection_index_2,
				    int *blockmap_index, int num_blocks, PDB_File *PDB, int num_pdbs, int lores, 
				    int hires, int total_residues, int *KO, int nko, int **cypa_contacts)
{
  int pdb,res;
  int bi,bj,ti,tj,sb;
  int q,q1,q2;
  int i,j,k,p,ii;

  /* Zero HT */
  for(i=1; i<=num_block_contacts; i++)
    for(j=1; j<=6; j++)
      for(k=j; k<=6; k++)
	HT[i][j][k] = HT[i][k][j] = 0;

  /* Find the pdb file and residue number of lores */
  pdb = 1;
  j = PDB[pdb].num_residues;
  while(j < lores && pdb<=num_pdbs) j += PDB[++pdb].num_residues;
  res = lores - j + PDB[pdb].num_residues;

  /* Include residues from lores to hires */
  for(ii=lores; ii<=hires; ii++){
    if(res>PDB[pdb].num_residues) { res=1; pdb++; }
    if(PDB[pdb].atom[res].model!=0){

      CalcHessianSuperrowCypa(hessian_superrow, block_contacts, blockmap, blockmap_index, 
			      num_blocks, PDB, num_pdbs, ii, pdb, res, KO, nko, cypa_contacts);

      q1 = projection_index_1[3*(ii-1)+2];
      q2 = projection_index_1[3*(ii-1)+3];
      for(k=projection_index_1[3*(ii-1)+1]; k<projection_index_1[3*ii+1]; k++){
	if(k<q1) q=1;
	else if(k<q2) q=2;
	else q=3;
	i = projections_sorted_by_1->IDX[k][2];
	bi = (i-1)/6+1;
	ti = i-6*(bi-1);
	for(p=projection_index_2[i]; p<=num_projection_elements; p++){
	  j = projections_sorted_by_2->IDX[p][2];
	  bj = (j-1)/6+1;
	  sb = block_contacts[bi][bj];
	  if(sb!=0){
	    tj=j-6*(bj-1);
	    HT[sb][ti][tj] += projections_sorted_by_1->X[k] * projections_sorted_by_2->X[p] * hessian_superrow[projections_sorted_by_2->IDX[p][1]][q];
	  }
	  else p = projection_index_2[++j] - 1;
	}
      }
    }
    res++;
  }/* <------ End of 'for(ii=1...*/
  return;
}

/* "CalcHessianSuperrowCypa" calculates the 'row'-th super-row of the Hessian, perturbing those involving 
   the 'nko' residues with resids listed in 'KO', and accounting for inter-block CypA binding. The 
   super-row corresponds to atom 'r1' of the 'p1'-th PDB object. Blocks connected by CypA are joined:
   'cypa_contacts[1..num_blocks][1..num_blocks]' contains the contact map between blocks connected with CypA. */
void CalcHessianSuperrowCypa(double **hessian_superrow, int **block_contacts, int **blockmap, 
			     int *blockmap_index, int num_blocks, PDB_File *PDB, int num_pdbs, 
			     int row, int p1, int r1, int *KO, int nko, int **cypa_contacts)
{
  double DX[3],csq=g_anm_cutoff*g_anm_cutoff,df,dd;
  int a,b,p2,r2;
  int i,j,k,q,jj,ii;
 

  /* Clear the diagonal super-element */
  for(i=1; i<=3; i++)
    for(j=1; j<=3; j++)
      hessian_superrow[3*(row-1)+i][j] = 0.0;

  /* Calculate submatrices using contacting blocks only */
  a = PDB[p1].atom[r1].model;
  for(b=1; b<=num_blocks; b++){
    if(block_contacts[a][b]!=0){
      for(ii=blockmap_index[b]; ii<blockmap_index[b+1]; ii++){
	p2 = blockmap[ii][2];
	r2 = blockmap[ii][3];
	jj = blockmap[ii][1];
	if(jj!=row){
	  dd=0.0;
	  for(k=0; k<3; k++){
	    DX[k] = (double)PDB[p1].atom[r1].X[k] - PDB[p2].atom[r2].X[k];
	    dd+=(DX[k]*DX[k]);
	  }
	  /* Either the distance is within cutoff, or these are CypA binding sites */
	  if(dd<csq || PDB[p1].atom[r1].resnum>=88 && PDB[p1].atom[r1].resnum<=90 && 
	     PDB[p2].atom[r2].resnum>=88 && PDB[p2].atom[r2].resnum<=90 && cypa_contacts[a][b]!=0){
	    for(i=1; i<=3; i++){
	      for(j=i; j<=3; j++){
		df = g_anm_gamma*DX[i-1]*DX[j-1]/dd;

		/* The perturbation */
		for(q=1; q<=nko; q++)
		  if(PDB[p1].atom[r1].resnum==KO[q] || PDB[p2].atom[r2].resnum==KO[q]){ 
		    df*=g_anm_eta;
		    break;
		  }

		/* Off-diagonal super-elements */
		hessian_superrow[3*(jj-1)+i][j] = hessian_superrow[3*(jj-1)+j][i] = -df;
	      
		/* Diagonal super-elements */
		hessian_superrow[3*(row-1)+i][j] += df;
		if(i != j) hessian_superrow[3*(row-1)+j][i] += df;
	      }
	    }
	  }
	  /* Zero the elements that are lacking contacts */
	  else
	    for(i=1; i<=3; i++)
	      for(j=1; j<=3; j++)
		hessian_superrow[3*(jj-1)+i][j] = hessian_superrow[3*(jj-1)+j][i] = 0.0;
	}
      }
    }
  }
}


/* "ConvertSuperblocksToMatrix" transfers the block Hessian from 
   the tensor 'superblocks' into the array 'block_hessian' */
int ConvertSuperblocksToMatrix(double **block_hessian, double ***superblocks, int **block_contacts, int num_blocks)
{
  int *I1,*I2,i,j,p,sb,ii,jj,max,a,b,imx;
 

  max = 6*num_blocks;
  I1 = ivector(1, max);
  I2 = ivector(1, max);

  /* I1[i]==i iff there is a non-zero element in column i (removes zeroes that are caused by single-node blocks) */
  for(i=1; i<=max; i++){ 
    I1[i] = I2[i] = 0;
    for(j=i; j<=max; j++) block_hessian[i][j] = block_hessian[j][i] = 0.0;
  }
  for(ii=1; ii<=num_blocks; ii++){
    for(jj=ii; jj<=num_blocks; jj++){
      sb = block_contacts[ii][jj];
      if(sb!=0){
	for(i=1; i<=6; i++){
	  p = jj==ii ? i : 1;
	  for(j=p; j<=6; j++)
	    if(fabs(superblocks[sb][i][j]) > DBL_EPSILON){
	      I1[6*(jj-1)+j] = 6*(jj-1) + j;
	    }
	}
      }
    }
  }


  /* If I1[i]!=0, then I2[i] is a sequential index */
  imx = 0;
  for(i=1; i<=max; i++){
    if(I1[i]!=0) imx++;
    I2[i] = imx;
  }
  for(ii=1; ii<=num_blocks; ii++){
    for(jj=ii; jj<=num_blocks; jj++){
      sb = block_contacts[ii][jj];
      if(sb!=0){
	for(i=1; i<=6; i++){
	  p = jj==ii ? i : 1;
	  for(j=p; j<=6; j++)
	    if(fabs(superblocks[sb][i][j]) > DBL_EPSILON){
	      a = I2[6*(ii-1)+i];
	      b = I2[6*(jj-1)+j];
	      block_hessian[a][b] = block_hessian[b][a] = superblocks[sb][i][j];
	    }
	}
      }
    }
  }
  
  free_ivector(I1,1,max);
  free_ivector(I2,1,max);
  return imx;
}


int PrintHessianMatrix(char *outfile, double **hessian, long dim)
{
  FILE *data;
  int i,j;

  if((data=fopen(outfile,"w"))==NULL){
    fprintf(stderr,"\nERROR: Unable to open %s\n",outfile);
    return -1;}
  for(i=1; i<=dim; i++)
    for(j=i; j<=dim; j++)
      if(fabs(hessian[i][j]) > 1.0e-10)
	fprintf(data,"%8d%8d% 20.10e\n", i, j, hessian[i][j]);
  fclose(data);
  return 0;
}


int PrintSparseHessianMatrix(char *outfile, double *A, int *RIDX, int *CIDX, long num_hess_elements)
{
  FILE *data;
  int i;

  if((data=fopen(outfile,"w"))==NULL){
    printf("\nUnable to open %s\n\n",outfile);
    return -1;}
  for(i=0; i<num_hess_elements; i++) fprintf(data,"%8d%8d% 20.10e\n", RIDX[i], CIDX[i], A[i]);
  fclose(data);
  return 0;
}


/* "BlockVecs2Full" projects eigenvectors from the block space back into the all-residue space. 
   'projection_matrix' is a sparse 'fdim' by 'bdim' matrix with 'num_projection_elements' elements.
   Other variables are: block_vectors[1..bdim][1..num_modes], full_vectors[1..fdim][1..num_modes]. */
void BlockVecs2Full(dSparse_Matrix *projection_matrix, double **block_vectors, double **full_vectors,
		    long num_projection_elements, long bdim, long fdim, int num_modes)
{
  double *MAG;
  double dd;
  int kold;
  int i,j,k,p;

  if(num_projection_elements < 4*fdim){
    printf("Condensing projection matrix to eliminate %ld columns\n",
  	    (4*fdim - num_projection_elements)/3);
    CondensePrjMtx(projection_matrix, num_projection_elements); }

  /* Calculate the magnitude of each vector.  This may be superfluous... */
  MAG = dvector(1, num_modes);
  for(i=1; i<=num_modes; i++) MAG[i] = 0.0;
  k = kold = 1;
  for(i=1; i<=fdim; i++){
    for(j=1; j<=num_modes; j++){
      dd = 0.0;
      k = kold;
      while(k<=num_projection_elements && projection_matrix->IDX[k][1]==i){
	p = projection_matrix->IDX[k][2];
	dd += projection_matrix->X[k] * block_vectors[p][j];
	k++;
      }
      MAG[j] += (dd*dd);
    }
    kold = k;
  }
  for(i=1; i<=num_modes; i++) MAG[i] = sqrt(MAG[i]);


  /* Store the normalized vectors in full_vectors */
  k = kold = 1;
  for(i=1; i<=fdim; i++){
    for(j=1; j<=num_modes; j++){
      dd = 0.0;
      k = kold;
      while(k<=num_projection_elements && projection_matrix->IDX[k][1]==i){
	p = projection_matrix->IDX[k][2];
	dd += projection_matrix->X[k] * block_vectors[p][j];
	k++;
      }
      full_vectors[i][j] = dd/MAG[j];
    }
    kold = k;
  }
  free_dvector(MAG, 1, num_modes);
}


/* "ConvertSuperblocksToSparseMatrix" transfers the block Hessian from 
   the tensor 'superblocks' into the dSparse_Matrix 'block_hessian' */
int ConvertSuperblocksToSparseMatrix(double **A, int **RIDX, int **CIDX, double ***superblocks, int **block_contacts, int num_blocks, long *num_elements)
{
  int *I1,*I2,i,j,p,sb,ii,jj,max,a,b,imx;
  long index = 0;
  
  *num_elements = 0;
  max = 6*num_blocks;
  I1 = ivector(1, max);
  I2 = ivector(1, max);

  /* I1[i]==i iff there is a non-zero element in column i (removes zeroes that are caused by single-node blocks) */
  for(i=1; i<=max; i++) I1[i] = I2[i] = 0;
  for(ii=1; ii<=num_blocks; ii++){
    for(jj=ii; jj<=num_blocks; jj++){
      sb = block_contacts[ii][jj];
      if(sb!=0){
	for(i=1; i<=6; i++){
	  p = jj==ii ? i : 1;
	  for(j=p; j<=6; j++)
	    if(fabs(superblocks[sb][i][j]) > DBL_EPSILON){
	      I1[6*(jj-1)+j] = 6*(jj-1) + j;
	      (*num_elements)++;
	    }
	}
      }
    }
  }

  /* If I1[i]!=0, then I2[i] is a sequential index */
  imx = 0;
  for(i=1; i<=max; i++){
    if(I1[i]!=0) imx++;
    I2[i] = imx;
  }

  *A = dvector(0, (*num_elements)-1);
  *RIDX = ivector(0, (*num_elements)-1);
  *CIDX = ivector(0, (*num_elements)-1);

  for(ii=1; ii<=num_blocks; ii++){
    for(jj=ii; jj<=num_blocks; jj++){
      sb = block_contacts[ii][jj];
      if(sb!=0){
	for(i=1; i<=6; i++){
	  p = jj==ii ? i : 1;
	  for(j=p; j<=6; j++)
	    if(fabs(superblocks[sb][i][j]) > DBL_EPSILON){
	      a = I2[6*(ii-1)+i];
	      b = I2[6*(jj-1)+j];
	      (*A)[index] = superblocks[sb][i][j];
	      (*RIDX)[index] = a;
	      (*CIDX)[index] = b;
	      index++;
	    }
	}
      }
    }
  }
  
  free_ivector(I1,1,max);
  free_ivector(I2,1,max);
  return imx;
}

/*
int PrintEigenSystemMulti(char *outfile, dSparse_Matrix *projection_matrix, double **block_vectors, 
			   double *eigenvalues, int num_projection_elements, int bdim, int total_residues, 
			   int num_modes)
{
  FILE *data;
  double **full_vectors;
  int fdim = 3*total_residues;
  int loeig = 1;
  int i;
  
  full_vectors = dmatrix(1, fdim, 1, num_modes+6);

  blockVecs2Full(projection_matrix, block_vectors, full_vectors, num_projection_elements, bdim, fdim, num_modes+6);

  for(i=1; i<=num_modes; i++)
    if(eigenvalues[i]<1.0e-6) loeig++;
  if(loeig!=7) printf("\n\t*** WARNING: %d ZERO EIGENVALUES DETECTED!!! ***\n\n",loeig-1);
  outfile[(int)(strlen(outfile)-strlen(strrchr(outfile,'.')))+1]='\0';
  strcat(outfile,"val");
  printf("Writing eigenvalues to %s\n",outfile);
  if((data=fopen(outfile,"w"))==NULL){
    printf("\nUnable to open %s\n\n",outfile);
    return 1;}
  for(i=loeig; i<=num_modes+loeig-1; i++) fprintf(data, "%e\n", eigenvalues[i]);
  fclose(data);
  //fprintPerturbVecsMulti(PDB,pdbfile,VEC,npdb,fdim,loeig,CL.neig+loeig-1,knockout);


  free_dmatrix(full_vectors, 1, fdim, 1, num_modes);
}
*/
