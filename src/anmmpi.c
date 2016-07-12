/*
  Licensed under the MIT License (MIT)

  Copyright (c) 2008-2015 Timothy Lezon

  This file is part of the ANMMPI software package

  Permission is hereby granted, free of charge, to any person obtaining a 
  copy of this software and associated documentation files (the "Software"), 
  to deal in the Software without restriction, including without limitation 
  the rights to use, copy, modify, merge, publish, distribute, sublicense, 
  and/or sell copies of the Software, and to permit persons to whom the 
  Software is furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
  DEALINGS IN THE SOFTWARE.


  ------------------------------------------------------------------------
  Author: Tim Lezon

  This is a partially parallel version of the Anisotropic Network Model (ANM)
  suitable for RTB calculations. This program generates a block Hessian matrix 
  that can be decomposed elsewhere. It further enables residue perturbations,
  as applied to the HIV-1 capsid and described by Bergman & Lezon [REF]. 

  Compile with mpicc.

  Run:
  $ mpirun -np 4 anmmpi file.pdb [OPTIONS]
  
  Allowed flags:
     -c      Specify cutoff distance
     -s      Specify perturbation scaling factor
     -n      Specify number of eigenvalue/eigenvector pairs to calculate
     -rtb    Assume input is RTB definitions, instead of PDB file
     -ko     Specify knockouts
     -b      Specify CypA binding site (HIV-1 specific)
     -p      Suppress normal output and print projection matrix
     -h      Suppress normal output and print Hessian matrix
*/

#include "nrstuff.h"
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<stddef.h>
#include<mpi.h>
#define PDB_MAX_LINE 90
#define RTB_WORD_LENGTH 256      /* Size of a reasonable file name */
#define RTB_BUFF_LENGTH 1024     /* Buffer size */


/* Global variables: ANM parameters */
double g_anm_cutoff = 15.0;
double g_anm_gamma = 1.0;
double g_anm_eta = 10.0;
int g_bound_resid = 89;
int g_num_eigs = 20;
const double g_anm_cutoff_default = 15.0;
const double g_anm_gamma_default = 1.0;
const double g_anm_eta_default = 10.0;


/* Structures */
typedef struct {int blknum; char LORES[4]; char lochain; int lonum; char HIRES[4]; char hichain; int hinum; int pdbid; } Rigid_Block;
typedef struct {char HEAD[7]; int atmnum; char ATOM[5]; char RES[4]; char chain; int resnum; float X[3]; float beta; char ELEMENT[3]; int model;} Atom_Line;
typedef struct {char **HEADER; Atom_Line *atom; int num_residues; int num_headers;} PDB_File;
typedef struct {int **IDX; double *X; } dSparse_Matrix;

/* Structure to hold command-line arguments */
typedef struct {char *infile; int is_blocks; int *KO; int nko; int print_prj_mtx; int add_cypa; int print_hessian_mtx;} Local_CL_Opt;

void ReadCommandLine(int argc, char *argv[], Local_CL_Opt *CL);
void CondensePrjMtx(dSparse_Matrix *projection_matrix, int num_projection_elements);
int AssignChainStarts(PDB_File *pdb, int num_pdbs, int **chain_starts, double **chain_directions, 
		      double **mer_centers, int num_chains);
int AssignCypaContacts(PDB_File *pdb, int num_pdbs, int **chain_starts, double **chain_directions, 
		       int **cypa_contacts, int num_chains, int num_blocks, double **mer_centers, int num_mers);
void PrintCypaLines(char *filename, PDB_File *pdb, int **cypa_contacts, 
		    int **blockmap, int *blockmap_index, int num_blocks);
int hintelemstr(char *text,int *H,int nn);
int nintelemstr(char *text);
void SortIntMatrixByCol(int **RA,int rlo,int rhi,int clo,int chi,int idx);
void SortSparseMatrixByIndex(dSparse_Matrix *MM,int num_elements,int sort_index);
void CopySparseMatrix(dSparse_Matrix *source_matrix,dSparse_Matrix *target_matrix,int lo,int hi);
void AssignIndexStretches(int *stretches,dSparse_Matrix *PP,int elm,int n,int idx);
void PdbFileSize(char *file,int *num_headers,int *num_models,int *num_chains,int *num_ca,int *num_atoms);
void AllocatePdbFile(PDB_File *PDB,int num_headers,int num_residues);
int ReadPdbFile(char *file,PDB_File *PDB,int num_headers,int num_residues);
void free_pdb(PDB_File *PDB);
int AssignRigidBlocks(PDB_File *pdb,Rigid_Block *blocks,int num_pdbs,int num_input_blocks,int *max_block_size);
int CalcProjectionMatrix(dSparse_Matrix *projection_matrix,PDB_File *pdb,double **block_center_of_mass,int num_pdbs,int num_blocks,int max_block_size);
void ReadBlockfile(char *file,char ***pdbfile,int *num_pdbs,Rigid_Block **blocks,int *num_blocks);
int ConvertSuperblocksToMatrix(double **block_hessian,double ***superblocks,int **block_contacts,int num_blocks);
int FindBlockContactsCypa(int **block_contacts,PDB_File *PDB,int **blockmap,int *blockmap_index,double **block_center_of_mass,int num_pdbs,int num_blocks,double cut,int **chain_starts,int **cypa_contacts,int num_chains);
void AddFullHessRowsToBlockHessCypa(double ***HT,int num_block_contacts,double **hessian_superrow,int **block_contacts,dSparse_Matrix *projections_sorted_by_1,dSparse_Matrix *projections_sorted_by_2,int num_projection_elements,int **blockmap,int *projection_index_1,int *projection_index_2,int *blockmap_index,int num_blocks,PDB_File *PDB,int num_pdbs,int lores,int hires,int total_residues,int *KO,int nko,int **cypa_contacts);
void CalcHessianSuperrowCypa(double **hessian_superrow,int **block_contacts,int **blockmap,int *blockmap_index,int num_blocks,PDB_File *PDB,int num_pdbs,int row,int p1,int r1,int *KO,int nko,int **cypa_contacts);

int main(int argc,char *argv[])
{
  /* Command-line variables */
  Local_CL_Opt CL;
  double t0,tdiff,tlast;
  FILE *data;
  PDB_File *PDB;
  Rigid_Block *blocks;
  dSparse_Matrix projection_matrix;
  dSparse_Matrix *projections_sorted_by_1,*projections_sorted_by_2;
  double **hessian_superrow,***superblocks,***superblocks_partial;
  int **block_contacts,*projection_index_1,*projection_index_2;
  int **blockmap,*blockmap_index;
  double **HH;
  double **block_center_of_mass;
  double **mer_centers;
  double *A,*VAL,**VEC;
  int *RIDX,*CIDX;
  char **pdbfilenames;
  char *outfile;
  char *bunk;
  double dd;
  int num_residues,total_residues,num_input_blocks;
  int max_block_size,num_blocks,bdim;
  int num_pdbs,num_headers,num_models;
  int num_projection_elements;
  int f,i,j,k,p;
  int kold;
  int num_block_contacts;
  int num_chains=0,num_mers=0;
  int **chain_starts;
  double **chain_directions;
  int **cypa_contacts;
  int num_bound_cypas=0;
  int ne;

  /* MPI stuff */
  MPI_Status status;
  int nprocs,rank;
  int lores,hires;
  int avgload;
  int send_tag=0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  /* Formalities */
  ReadCommandLine(argc,argv,&CL);
  if(rank==0){
    t0 = MPI_Wtime();
    if(CL.is_blocks==1) fprintf(stderr,"Blockfile: %s\n",CL.infile);
    else fprintf(stderr,"PDB file: %s\n",CL.infile);
    fprintf(stderr,"Cutoff: %.3f\n",g_anm_cutoff);
    if(CL.nko!=0){ 
      fprintf(stderr,"knockouts:");
      for(i=1; i<=CL.nko; i++) fprintf(stderr,"\t%d",CL.KO[i]);
      fprintf(stderr,"\nPerturbation scaling is %e\n",g_anm_eta);
    }
    else fprintf(stderr,"No knockouts\n");
    if(CL.add_cypa==1) fprintf(stderr,"Binding chains at residue %d\n",g_bound_resid);
  }
  if(CL.is_blocks==1) ReadBlockfile(argv[1],&pdbfilenames,&num_pdbs,&blocks,&num_input_blocks);
  else{
    pdbfilenames = cmatrix(1,1,0,RTB_WORD_LENGTH);
    strncpy(pdbfilenames[1],argv[1],strlen(argv[1]));
    num_pdbs = 1;
  }
  if(rank==0){
    outfile = (char *)calloc((size_t) (strlen(argv[1])+99),sizeof(char));
    bunk = (char *)calloc(99,sizeof(char));
    strncpy(outfile,argv[1],strlen(argv[1])-strlen(strrchr(argv[1],'.')));
    strcat(outfile,".something");
  }

  /* Read the PDB files */
  PDB = (PDB_File *)malloc( (size_t) (1+num_pdbs)*sizeof(PDB_File));
  if(!PDB){
    fprintf(stderr,"\nPID% 3d: Allocation failure in PDB\n\n",rank);
    exit(1);}
  for(f=1; f<=num_pdbs; f++){
    PdbFileSize(pdbfilenames[f],&num_headers,&num_models,&j,&num_residues,&i);
    num_chains += j;
    if(num_models!=1){
      fprintf(stderr,"\n%s: %d models (needs to be one)!\n\n",pdbfilenames[f],num_models);
      exit(1);}
    AllocatePdbFile(&(PDB[f]),0,num_residues);
    ReadPdbFile(pdbfilenames[f],&PDB[f],0,num_residues);
  }

  /* Assign each residue to a block.  Blocks are given sequential indices starting at 1.  */
  if(CL.is_blocks==1){ 
    num_blocks = AssignRigidBlocks(PDB,blocks,num_pdbs,num_input_blocks,&max_block_size);
    free(blocks);}
  else{
    num_blocks = num_residues;
    max_block_size = 1;
    for(i=1; i<=num_residues; i++) PDB[1].atom[i].model = i;
  }
  if(rank==0){
    fprintf(stderr,"%d total rigid blocks in structure\n",num_blocks);
    fprintf(stderr,"Largest block contains %d residues\n",max_block_size); }

  /* Find the projection matrix.  The 3 coordinates of each residue project onto the 6 coordinates of 
     exactly one block.  The translational component of the projection will always be a diagonal 
     matrix, so there are 12 required components of the projection matrix for each residue.  The 
     coordinates of the block mass centers are stored in block_center_of_mass for later use. */
  total_residues = 0;
  for(i=1; i<=num_pdbs; i++) total_residues += PDB[i].num_residues;
  if(rank==0) fprintf(stderr,"Total number of residues: %d\n",total_residues);
  projection_matrix.IDX = imatrix(1,12*total_residues,1,2);
  projection_matrix.X = dvector(1,12*total_residues);
  block_center_of_mass = dmatrix(1,num_blocks,1,3);
  num_projection_elements = CalcProjectionMatrix(&projection_matrix, PDB, block_center_of_mass, num_pdbs, 
						 num_blocks, max_block_size);
  if( num_projection_elements > 12*total_residues ){
    fprintf(stderr,"\nUnexpected size for projection matrix.\n");
    fprintf(stderr,"This may result from an unassigned residue.\n\n");
    exit(1);}
  if(rank==0) fprintf(stderr,"%d non-zero elements in the projection matrix\n",num_projection_elements);
  SortSparseMatrixByIndex(&projection_matrix,num_projection_elements,1);
  projections_sorted_by_1 = &projection_matrix; 

  /* Print the projection matrix to file if requested */
  if(CL.print_prj_mtx==1){
    if(rank==0){
      if(num_projection_elements < 12*total_residues){
	fprintf(stderr,"Condensing projection matrix to eliminate %d columns\n",
		(12*total_residues-num_projection_elements)/3);
	CondensePrjMtx(&projection_matrix, num_projection_elements); }
      outfile[(int)(strlen(outfile)-strlen(strrchr(outfile,'.')))+1]='\0';
      strcat(outfile,"prj");
      if((data=fopen(outfile,"w"))==NULL){
	fprintf(stderr,"\n%s: Unable to open %s\n\n",argv[0],outfile);
	exit(1);}
      fprintf(stderr,"Writing projection matrix to %s\n",outfile);
      for(i=1; i<=num_projection_elements; i++) 
	fprintf(data,"%8d%8d% 25.15e\n", projection_matrix.IDX[i][1], 
		projection_matrix.IDX[i][2], projection_matrix.X[i]);
      tdiff = MPI_Wtime() - t0;
      fprintf(stderr,"Finished in %e seconds\n",tdiff);
      fclose(data);
      free(outfile);
      free(bunk);
    }
    free_dmatrix(block_center_of_mass, 1, num_blocks, 1, 3);
    free_dvector(projection_matrix.X,1,12*total_residues);
    free_imatrix(projection_matrix.IDX,1,12*total_residues,1,2);
    for(f=1; f<=num_pdbs; f++) free_pdb(&(PDB[f]));
    free(PDB);
    free(CL.infile);
    MPI_Finalize();
    return 0;
  }

  /* ------------- Find out which CypA binding sites should be connected ---------- */
  /* 
     A. Assign a unique chain ID to each chain, and a unique mer ID to each hexa/pentamer
     B. Determine which chains are in which mers by distance between intra-hexamer residues.
     -- Now each block belongs to a single chain, and each chain belongs to a single mer.
     C. Calculate a local orientation vector for each chain. This is just the vector pointing
     -- from the mer CM to the chain CM. 
     D. Use local orientation vectors and positions of CypA loops to determine which chains 
     -- are connected by CypA binding. This information needs to be passed to the function
     -- that calculates the Hessian matrix.
  */
  /* Assign a unique chain ID to each chain. These will be sequential: 1...num_chains */
  /* Chain 'i' starts at PDB[chain_starts[i][1]].atom[chain_starts[i][2]] */
  /* Chain 'i' is in mer chain_starts[i][3] */
  chain_starts = imatrix(1,num_chains,1,3);
  chain_directions = dmatrix(1,num_chains,0,2);
  mer_centers = dmatrix(1,num_chains,0,2);
  cypa_contacts = imatrix(1,num_blocks,1,num_blocks);
  if(CL.add_cypa==1){
    num_mers = AssignChainStarts(PDB, num_pdbs, chain_starts, chain_directions, mer_centers, num_chains);
    if(rank==0) fprintf(stderr,"Found %d hexa/pentamers\n",num_mers);
    num_bound_cypas = AssignCypaContacts(PDB, num_pdbs, chain_starts, chain_directions, cypa_contacts, 
					 num_chains, num_blocks, mer_centers, num_mers);
  }
  else
    for(i=1; i<=num_blocks; i++)
      for(j=1; j<=num_blocks; j++)
	cypa_contacts[i][j] = 0;
  if(rank==0 && num_bound_cypas!=0) fprintf(stderr,"%d CypA molecules bound\n",num_bound_cypas);
  

  /* -------------------- Calculate the block Hessian ------------------------ */
  /*  blockmap is a mapping between blocks and residues:
      blockmap[i][1] is global residue index; blockmap[i][2] is PDB index; 
      blockmap[i][3] is resnum; blockmap[i][4] is block 
      blockmap gets sorted by the block index, then blockmap_index stores the beginning
      positions of each of the blocks in blockmap: 
      Block 'i' is rows blockmap_index[i] to blockmap_index[i+1]-1 of blockmap.
      projection_index_1 and projection_index_2 are defined such that 
      for all j: projection_index_1[i] <= j < projection_index_1[i+1], projections_sorted_by_1->IDX[j][1] = i 
      for all j: projection_index_2[i] <= j < projection_index_2[i+1], projections_sorted_by_2->IDX[j][2] = i
  */
  block_contacts = imatrix(0, num_blocks, 0, num_blocks);
  blockmap = imatrix(1, total_residues, 1, 4);
  blockmap_index = ivector(1, num_blocks+1);
  projections_sorted_by_2 = (dSparse_Matrix *)malloc((size_t)sizeof(dSparse_Matrix));
  projections_sorted_by_2->IDX = imatrix(1, num_projection_elements, 1, 2);
  projections_sorted_by_2->X = dvector(1, num_projection_elements);
  CopySparseMatrix(projections_sorted_by_1, projections_sorted_by_2, 1, num_projection_elements);
  SortSparseMatrixByIndex(projections_sorted_by_2, num_projection_elements, 2);
  projection_index_1 = ivector(1, 3*total_residues+1);
  projection_index_2 = ivector(1, 6*num_blocks+1);

  k=1;
  for(p=1; p<=num_pdbs; p++)
    for(i=1; i<=PDB[p].num_residues; i++){
      blockmap[k][1] = k;
      blockmap[k][2] = p;
      blockmap[k][3] = i;
      blockmap[k][4] = PDB[p].atom[i].model;
      k++;
    }
  SortIntMatrixByCol(blockmap, 1, total_residues, 1, 4, 4);
  i = k = 1;
  while(blockmap[k][4] < 1) k++;
  blockmap_index[i] = k;
  while(i<=num_blocks){
    while( k<=total_residues && blockmap[k][4]==blockmap[blockmap_index[i]][4] ) k++;
    blockmap_index[++i] = k;
  }
  AssignIndexStretches(projection_index_1, projections_sorted_by_1, 
		       num_projection_elements, 3*total_residues+1, 1);
  AssignIndexStretches(projection_index_2, projections_sorted_by_2, 
		       num_projection_elements, 6*num_blocks+1, 2);
  if(rank==0){ 
    fprintf(stderr,"Finding block contacts...");
    tlast = MPI_Wtime();
  }

  num_block_contacts = FindBlockContactsCypa(block_contacts, PDB, blockmap, blockmap_index, 
					     block_center_of_mass, num_pdbs, num_blocks, g_anm_cutoff, 
					     chain_starts, cypa_contacts, num_chains);

  /* Create a TCL file showing CypA as cylinders */
  if(rank==0){ 
    tdiff = MPI_Wtime() - tlast;
    fprintf(stderr,"Finished in %e seconds\n",tdiff);
    fprintf(stderr,"%d block contacts found\n",num_block_contacts);
    if(CL.add_cypa==1){
      outfile[(int)(strlen(outfile)-strlen(strrchr(outfile,'.')))+1]='\0';
      strcat(outfile,"tcl");
      PrintCypaLines(outfile, PDB, cypa_contacts, blockmap, blockmap_index, num_blocks);
    }
  }

  /* ---------------------- Parallel --------------------- */
  /* Everything up to this point has been done on all processors simultaneously. It's a waste of compute 
     power, but it's an easy way to bypass the need to broadcast a bunch of easily calculated values. 
     Here we will actually divide things up by process. */
  if(rank==0){ 
    fprintf(stderr,"Calculating block Hessian...");
    tlast = MPI_Wtime();
  }
  superblocks_partial = d3tensor(1,num_block_contacts,1,6,1,6);
  hessian_superrow = dmatrix(1,3*total_residues,1,3);
  avgload = total_residues/nprocs;
  lores = rank*avgload + 1;
  hires = (rank + 1)*avgload;
  if (rank==nprocs - 1) hires = total_residues;
  AddFullHessRowsToBlockHessCypa(superblocks_partial, num_block_contacts, hessian_superrow, block_contacts,
				 projections_sorted_by_1, projections_sorted_by_2, num_projection_elements,
				 blockmap, projection_index_1, projection_index_2, blockmap_index, num_blocks,
				 PDB, num_pdbs, lores, hires, total_residues, CL.KO, CL.nko, cypa_contacts);

  /* Divide the matrix calculation among processes. */
  if(rank==0){
    superblocks = d3tensor(1,num_block_contacts,1,6,1,6);
    for(i=1; i<=num_block_contacts; i++)
      for(j=1; j<=6; j++)
	for(k=1; k<=6; k++)
	  superblocks[i][j][k] = superblocks_partial[i][j][k];
    for(p=1; p<nprocs; p++){
      MPI_Recv(&superblocks_partial[1][1][1], 6*6*num_block_contacts, MPI_DOUBLE,
	       MPI_ANY_SOURCE, send_tag, MPI_COMM_WORLD, &status);
      for(i=1; i<=num_block_contacts; i++)
	for(j=1; j<=6; j++)
	  for(k=1; k<=6; k++)
	    superblocks[i][j][k] += superblocks_partial[i][j][k];
    }
  }
  else MPI_Send(&superblocks_partial[1][1][1], 6*6*num_block_contacts, MPI_DOUBLE, 0, send_tag, MPI_COMM_WORLD);

  if(rank==0){
    HH = dmatrix(1,6*num_blocks,1,6*num_blocks);
    bdim = ConvertSuperblocksToMatrix(HH,superblocks,block_contacts,num_blocks);
    free_d3tensor(superblocks,1,num_block_contacts,1,6,1,6);
    tdiff = MPI_Wtime() - tlast;
    fprintf(stderr,"Finished in %e seconds\n",tdiff);
  }

  free_imatrix(block_contacts, 0, num_blocks, 0, num_blocks);
  free_imatrix(blockmap, 1, total_residues, 1, 4);
  free_ivector(blockmap_index, 1, num_blocks+1);
  free_imatrix(projections_sorted_by_2->IDX, 1, num_projection_elements, 1, 2);
  free_dvector(projections_sorted_by_2->X, 1, num_projection_elements);
  free(projections_sorted_by_2);
  free_ivector(projection_index_1, 1, 3*total_residues+1);
  free_ivector(projection_index_2, 1, 6*num_blocks+1);
  free_d3tensor(superblocks_partial, 1, num_block_contacts, 1, 6, 1, 6);
  free_dmatrix(hessian_superrow, 1, 3*total_residues, 1, 3);
  /* ------------------------------- End of parallel part --------------------------------- */


  /* Print block Hessian to file, if requested */ 
  if(CL.print_hessian_mtx==1){
    if(rank==0){
      outfile[(int)(strlen(outfile)-strlen(strrchr(outfile,'.')))+1]='\0';
      strcat(outfile,"blockhessian");
      fprintf(stderr,"Writing block Hessian to %s\n",outfile);
      if((data=fopen(outfile,"w"))==NULL){
	fprintf(stderr,"\n%s: Unable to open %s\n\n",argv[0],outfile);
	exit(1);}
      for(i=1; i<=bdim; i++)
	for(j=i; j<=bdim; j++)
	  if(fabs(HH[i][j]) > 1.0e-10)
	    fprintf(data,"%8d%8d% 20.10e\n",i,j,HH[i][j]);
      fclose(data);
      tdiff = MPI_Wtime() - t0;
      fprintf(stderr,"Elapsed time: %e seconds\n",tdiff);
      free(outfile);
      free(bunk);
    }
    free_cmatrix(pdbfilenames,1,num_pdbs,0,RTB_WORD_LENGTH);
    for(f=1;f<=num_pdbs;f++) free_pdb(&(PDB[f]));
    free(PDB);
    free(CL.infile);
    free_imatrix(projection_matrix.IDX,1,12*total_residues,1,2);
    free_dvector(projection_matrix.X,1,12*total_residues);
    free_dmatrix(block_center_of_mass,1,num_blocks,1,3);
    free_imatrix(chain_starts,1,num_chains,1,3);
    free_dmatrix(chain_directions,1,num_chains,0,2);
    free_dmatrix(mer_centers,1,num_chains,0,2);
    free_imatrix(cypa_contacts,1,num_blocks,1,num_blocks);
    MPI_Finalize();
    return 0;
  }

  /* Decompose the Hessian */
  else if(rank==0){
    ne = 0;
    for(i=1; i<=bdim; i++)
      for(j=i; j<=bdim; j++)
	if(fabs(HH[i][j]) > 1.0e-10) ne++;
    A = dvector(0,ne-1);
    RIDX = ivector(0,ne-1);
    CIDX = ivector(0,ne-1);
    k = 0;
    for(i=1; i<=bdim; i++)
      for(j=i; j<=bdim; j++)
	if(fabs(HH[i][j]) > 1.0e-10){
	  RIDX[k] = i;
	  CIDX[k] = j;
	  A[k] = HH[i][j];
	  k++;
	}
    free_dmatrix(HH, 1, 6*num_blocks, 1, 6*num_blocks);
    VAL = dvector(1,g_num_eigs);
    VEC = dmatrix(1,bdim,1,g_num_eigs);
    blzSolveSparse(A, RIDX, CIDX, ne, VAL, VEC, bdim, g_num_eigs);
    outfile[(int)(strlen(outfile)-strlen(strrchr(outfile,'.')))+1]='\0';
    strcat(outfile,"val");
    data = fopen(outfile,"w");
    for(i=1; i<=g_num_eigs; i++) fprintf(data,"% 16.7e\n",VAL[i]);
    fclose(data);
    outfile[(int)(strlen(outfile)-strlen(strrchr(outfile,'.')))+1]='\0';
    strcat(outfile,"vec");
    data = fopen(outfile,"w");

    /* Print the block vectors */
    /*
    for(i=1; i<=bdim; i++){
      for(j=1; j<=g_num_eigs; j++) fprintf(data,"% 16.7e",VEC[i][j]);
      fprintf(data,"\n");
    }
    */

    /* Print the all-residue vectors */
    k = kold = 1;
    for(i=1; i<=total_residues; i++){
      for(j=1; j<=g_num_eigs; j++){
	dd = 0.0;
	k = kold;
	while(k<num_projection_elements && projection_matrix.IDX[k][1]==i){
	  p = projection_matrix.IDX[k][2];
	  dd += projection_matrix.X[k]*VEC[p][j];
	  k++;
	}
	fprintf(data,"% 16.7e",dd);
      }
      kold = k;
      fprintf(data,"\n");
    }
    fclose(data);

    /* Check time */
    tdiff = MPI_Wtime() - t0;
    fprintf(stderr,"Elapsed time: %e seconds\n",tdiff);
    free(outfile);
    free(bunk);
    free_dvector(A,0,ne-1);
    free_ivector(RIDX,0,ne-1);
    free_ivector(CIDX,0,ne-1);
    free_dvector(VAL,1,g_num_eigs);
    free_dmatrix(VEC,1,bdim,1,g_num_eigs);
  }

  /* ------------------------------ Free! ------------------------------ */
  free_cmatrix(pdbfilenames,1,num_pdbs,0,RTB_WORD_LENGTH);
  for(f=1;f<=num_pdbs;f++) free_pdb(&(PDB[f]));
  free(PDB);
  free(CL.infile);
  free_imatrix(projection_matrix.IDX,1,12*total_residues,1,2);
  free_dvector(projection_matrix.X,1,12*total_residues);
  free_dmatrix(block_center_of_mass,1,num_blocks,1,3);
  free_imatrix(chain_starts,1,num_chains,1,3);
  free_dmatrix(chain_directions,1,num_chains,0,2);
  free_dmatrix(mer_centers,1,num_chains,0,2);
  free_imatrix(cypa_contacts,1,num_blocks,1,num_blocks);

  /* ------------------------------  END  ------------------------------ */ 
  MPI_Finalize();  
  return 0;
}


/* "ReadCommandLine" reads the command line */
void ReadCommandLine(int argc,char *argv[],Local_CL_Opt *CL)
{
  int ok=0;
  int i;

  /* Allowed flags:
     -c      Specify cutoff distance
     -s      Specify perturbation scaling factor
     -n      Specify number of eigenvalue/eigenvector pairs to calculate
     -rtb    Assume input is RTB definitions, instead of PDB file
     -ko     Specify knockouts
     -l      Specify CypA binding site (HIV-1 specific)
     -p      Suppress normal output and print projection matrix
     -h      Suppress normal output and print Hessian matrix
  */

  if(argc > 1){
    ok = 1;

    /* Assume first argument is filename */
    CL->infile = (char *)malloc((size_t) (((int)strlen(argv[1])+1)*sizeof(char)));
    strcpy(CL->infile,argv[1]);
    CL->nko = 0;
    CL->print_prj_mtx = 0;
    CL->add_cypa = 0;
    CL->print_hessian_mtx = 0;
    CL->is_blocks = 0;

    /* Go through arguments */
    i=2;
    while(i<argc){
      if(strncmp(argv[i],"-c",2)==0) sscanf(argv[++i],"%lf",&g_anm_cutoff);      /* Cutoff */
      else if(strncmp(argv[i],"-ko",3)==0){
	CL->nko = nintelemstr(argv[++i]);
	CL->KO = ivector(1,CL->nko);
	hintelemstr(argv[i],CL->KO,CL->nko);
      }
      else if(strncmp(argv[i],"-p",2)==0) CL->print_prj_mtx = 1;                 /* Print projection matrix */
      else if(strncmp(argv[i],"-h",2)==0) CL->print_hessian_mtx = 1;             /* Print Hessian matrix */
      else if(strncmp(argv[i],"-rtb",4)==0) CL->is_blocks = 1;                /* Block File */
      else if(strncmp(argv[i],"-s",2)==0) sscanf(argv[++i],"%lf",&g_anm_eta);    /* Scale factor */
      else if(strncmp(argv[i],"-n",2)==0) sscanf(argv[++i],"%d",&g_num_eigs);    /* Number of eigenvalues */
      else if(strncmp(argv[i],"-l",2)==0){                                       /* CypA binding  */
	sscanf(argv[++i],"%d",&g_bound_resid);
	CL->add_cypa = 1;
      }
      else{
	fprintf(stderr,"\n%s: Unknown argument: %s\n\n",argv[0],argv[i]);
	ok=0;
	break;
      }
      i++;	
    }  /* <-- End of 'while(i<argc)...' */
  }

  /* Print error message if things are not all right */
  if(ok!=1){
    fprintf(stderr,"\nUsage:\nmpirun -np N %s file.dat [OPTIONS]\n\n",argv[0]);
    fprintf(stderr,"\tfile.dat -- PDB file\n\n");
    fprintf(stderr,"OPTIONS:\n");
    fprintf(stderr,"\t-c CUT\tSet ANM cutoff distance to CUT (def %.3f)\n",g_anm_cutoff);
    fprintf(stderr,"\t-s SCL\tSet perturbation scaling factor to SCL (def %.3f)\n",g_anm_eta);
    fprintf(stderr,"\t-n NUM\tCalculate NUM eigenvalue/eigenvector pairs (def %d)\n",g_num_eigs);
    fprintf(stderr,"\t-rtb\tAssume input is RTB definitions, instead of PDB file\n");
    fprintf(stderr,"\t-ko\tSpecify indices of residues to perturb\n");
    fprintf(stderr,"\t\tAccepts integers separated by commas (no white space!),\n");
    fprintf(stderr,"\t\tas well as ranges indicated with dashes\n");
    fprintf(stderr,"\t-l RES\tBind neighboring chains at RES (HIV-specific)\n");
    fprintf(stderr,"\t-p\tSuppress normal output and print projection matrix\n");
    fprintf(stderr,"\t-h\tSuppress normal output and print Hessian matrix\n");
    fprintf(stderr,"\nOutput:\nCalculates eigenvalues and eigenvectors\n\n");
    exit(1);
  }
}


/* "init_bmap" initializes the mapping from residues to blocks */
void init_bmap(int **blockmap,PDB_File *PDB,int num_pdbs)
{
  int i,p,n=0;

  for(p=1; p<=num_pdbs; p++){
    for(i=1; i<=PDB[p].num_residues; i++){
      n++;
      blockmap[n][1] = n;
      blockmap[n][2] = PDB[p].atom[i].model;
    }
  }
}

/* "CondensePrjMtx" condenses a projection matrix by offsetting the column 
   elements if there are blocks with fewer than six degrees of freedom. */
void CondensePrjMtx(dSparse_Matrix *projection_matrix,int num_projection_elements)
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

/* "AssignChainStarts" puts the PDB number and atom number of the beginning of each chain into the array 
   'chain_starts': Chain 'i' starts at PDB[chain_starts[i][1]].atom[chain_starts[i][2]] and is a member 
   of mer chain_starts[i][3]. Row 'i' of 'chain_directions' gives the vector pointing from the center of 
   mer chain_starts[i][3] to the center of chain i. */
int AssignChainStarts(PDB_File *pdb, int num_pdbs, int **chain_starts, 
		      double **chain_directions, double **mer_centers, int num_chains)
{
  char chain_old='\0';
  double dd,df;
  int num_mers;
  int num_parts_in_whole=0;
  int i,j=1,k,p;

  /* Zero chain directions */
  for(i=1; i<=num_chains; i++)
    for(k=0; k<3; k++)
      chain_directions[i][k] = 0.0;

  /* Find out where each chain begins, and find its center of mass. */
  chain_starts[j][1] = 1;
  chain_starts[j][2] = 1;
  chain_starts[j][3] = 0;
  chain_old = pdb[1].atom[1].chain;
  for(p=1; p<=num_pdbs; p++){
    for(i=1; i<=pdb[p].num_residues; i++){
      if(pdb[p].atom[i].chain != chain_old){

	/* Residue 89 */
	for(k=0; k<3; k++) 
	  chain_directions[j][k] = pdb[ chain_starts[j][1]].atom[chain_starts[j][2] + g_bound_resid - 1 ].X[k];
	j++;
	chain_starts[j][1] = p;
	chain_starts[j][2] = i;
	chain_starts[j][3] = 0;
	chain_old = pdb[p].atom[i].chain;
      }
    }
  }
  /* Account for the last chain */
  for(k=0; k<3; k++) 
    chain_directions[j][k] = pdb[ chain_starts[j][1]].atom[chain_starts[j][2] + g_bound_resid - 1 ].X[k];
  if(j != num_chains){
    fprintf(stderr,"\n * ERROR in 'AssignChainStarts'\n\n");
    exit(1);}

  /* We will say that if ARG18 residues on two chains are within 
     20 Angstroms of each other, the chains are in the same mer. */
  num_mers = 0;
  for(i=1; i<=num_chains; i++){
    if(chain_starts[i][3]==0) chain_starts[i][3] = ++num_mers;
    for(j=i+1; j<=num_chains; j++){
      if(chain_starts[j][3]==0){
	dd = 0.0;
	for(k=0; k<3; k++){
	  df = pdb[chain_starts[j][1]].atom[chain_starts[j][2]+17].X[k] 
	    - pdb[chain_starts[i][1]].atom[chain_starts[i][2]+17].X[k];
	  dd += df*df;
	}
	if(dd <= 400.0) chain_starts[j][3] = chain_starts[i][3];
      }
    }
  }


  /* Find the center of each mer, and the vectors pointing from 
     the mer center to the center of each chain it contains. */
  for(i=1; i<=num_mers; i++){

    /* Get the center of mass of each mer */
    num_parts_in_whole = 0;
    for(k=0; k<3; k++) mer_centers[i][k] = 0.0;
    for(j=1; j<=num_chains; j++)
      if( chain_starts[j][3] == i ){
	for(k=0; k<3; k++) mer_centers[i][k] += chain_directions[j][k];
	num_parts_in_whole++;
      }
    for(k=0; k<3; k++) mer_centers[i][k] /= (double)num_parts_in_whole;

    /* Get the relative direction of each chain as a unit vector */
    for(j=1; j<=num_chains; j++)
      if( chain_starts[j][3] == i ){
	dd = 0.0;
	for(k=0; k<3; k++){
	  chain_directions[j][k] -= mer_centers[i][k];
	  dd += chain_directions[j][k]*chain_directions[j][k];
	}
	dd = sqrt(dd);
	for(k=0; k<3; k++) chain_directions[j][k] /= dd;
      }
  }
  return num_mers;
}


/* "AssignCypaContacts" populates the matrix cypa_contacts[1..num_blocks][1..num_blocks] with integers 
   indicating which blocks are joined by CypA molecules. It returns the number of inter-CA CypA contacts 
   formed. The CypA binding loop is residues 88-90. Here it will be represented by one residue only. */
int AssignCypaContacts(PDB_File *pdb, int num_pdbs, int **chain_starts, double **chain_directions,
		       int **cypa_contacts, int num_chains, int num_blocks, double **mer_centers, int num_mers)
{
  double mer_vector,cypa_vector;
  double loop_loop_dist,mer_mer_dist;
  double chain_chain_angle,mer_cypa_angle,chain_cypa_angle;
  int *num_bound;
  int *block_for_chain;
  int num_contacts = 0;
  int i,j,k,p;

  /* num_bound[i] is # of CypA bound to chain i;
     block_for_chain[i] is index of block corresponding to residue 'g_bound_resid' of chain i. */
  num_bound = ivector(1,num_chains);
  block_for_chain = ivector(1,num_chains);
  for(i=1; i<=num_chains; i++) num_bound[i] = block_for_chain[i] = 0;
  for(i=1; i<=num_blocks; i++)
    for(j=i; j<=num_blocks; j++) cypa_contacts[i][j] = cypa_contacts[j][i] = 0;
  j = 0;
  for(p=1; p<=num_pdbs; p++)
    for(i=1; i<=pdb[p].num_residues; i++)
      if(pdb[p].atom[i].resnum == g_bound_resid) block_for_chain[++j] = pdb[p].atom[i].model;

  /* Each chain can only have one contact. */
  for(i=1; i<=num_chains; i++){
    if(num_bound[i]==0){
      for(j=i+1; j<=num_chains; j++){
	if(num_bound[j]==0){
	  loop_loop_dist = chain_chain_angle = mer_mer_dist = mer_cypa_angle = chain_cypa_angle = 0.0;
	  for(k=0; k<3; k++){

	    /* Angle between chains */
	    chain_chain_angle += chain_directions[i][k]*chain_directions[j][k];

	    /* Vector in direction of CypA */
	    cypa_vector = pdb[ chain_starts[j][1]].atom[chain_starts[j][2] + g_bound_resid - 1 ].X[k] 
	      - pdb[ chain_starts[i][1]].atom[chain_starts[i][2] + g_bound_resid - 1 ].X[k];
	    loop_loop_dist += cypa_vector*cypa_vector;

	    /* Vector between mer centers */
	    mer_vector = mer_centers[chain_starts[j][3]][k] - mer_centers[chain_starts[i][3]][k];
	    mer_mer_dist += mer_vector*mer_vector;
	    mer_cypa_angle += mer_vector*cypa_vector;
	    chain_cypa_angle += chain_directions[i][k]*cypa_vector;
	  }
	  loop_loop_dist = sqrt(loop_loop_dist);
	  mer_mer_dist = sqrt(mer_mer_dist);
	  chain_chain_angle = 180.0*acos(chain_chain_angle)/M_PI;
	  mer_cypa_angle /= (mer_mer_dist*loop_loop_dist);
	  mer_cypa_angle = 180.0*acos(mer_cypa_angle)/M_PI;
	  chain_cypa_angle /= loop_loop_dist;
	  chain_cypa_angle = 180.0*acos(chain_cypa_angle)/M_PI;

	  if(loop_loop_dist < 47.5 && chain_chain_angle>140.0){ /* Mock binding */
	    cypa_contacts[block_for_chain[i]][block_for_chain[j]] = 1;
	    cypa_contacts[block_for_chain[j]][block_for_chain[i]] = 1;
	    num_bound[i]++;
	    num_bound[j]++;
	    num_contacts++;
	    break;
	  }
	} /* <-- if(num_bound[j]==0)... */
      }   /* <-- for(j=i+1... */
    }     /* <--if(num_bound[i]==0)... */
  }
  free_ivector(num_bound,1,num_chains);
  free_ivector(block_for_chain,1,num_chains);
  return num_contacts;
}

/* "PrintCypaLines" generates a TCL file that shows CypA molecules as cylinders. */
void PrintCypaLines(char *filename, PDB_File *pdb, int **cypa_contacts, 
		    int **blockmap, int *blockmap_index, int num_blocks)
{
  FILE *data;
  int color = 0;
  int i,j,k;

  if((data=fopen(filename,"w"))==NULL){
    fprintf(stderr,"\nPrintCypaLines: Unable to open '%s'\n",filename);
    exit(1);}

  /* This prints cylinders */
  fprintf(data,"draw color %d\n",color);
  for(i=1; i<=num_blocks; i++)
    for(j=i+1; j<=num_blocks; j++)
      if(cypa_contacts[i][j] != 0){
	for(k=blockmap_index[i]; k<blockmap_index[i+1]; k++)
	  if(pdb[blockmap[k][2]].atom[blockmap[k][3]].resnum == g_bound_resid){
	    fprintf(data,"draw cylinder {%.3f %.3f %.3f} ", pdb[blockmap[k][2]].atom[blockmap[k][3]].X[0], 
		    pdb[blockmap[k][2]].atom[blockmap[k][3]].X[1], pdb[blockmap[k][2]].atom[blockmap[k][3]].X[2]);
	    break;
	  }
	for(k=blockmap_index[j]; k<blockmap_index[j+1]; k++)
	  if(pdb[blockmap[k][2]].atom[blockmap[k][3]].resnum == g_bound_resid){
	    fprintf(data,"{%.3f %.3f %.3f} radius 5.1 resolution 10\n", 
		    pdb[blockmap[k][2]].atom[blockmap[k][3]].X[0], pdb[blockmap[k][2]].atom[blockmap[k][3]].X[1], 
		    pdb[blockmap[k][2]].atom[blockmap[k][3]].X[2]);
	    break;
	  }
      }
  fclose(data);
  return;
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

/* "SortSparseMatrixByIndex" sorts a dSparse_Matrix containing 'num_elements' elements by the index 'sort_index' */
/* Replaces 'dsort_PP2' */
void SortSparseMatrixByIndex(dSparse_Matrix *MM,int num_elements,int sort_index)
{
  double x;
  int i,ir,j,l,hi,i1,i2,ndx;
  long rra,*ra;

  if(num_elements < 2) return;
  if(sort_index<1 || sort_index>2){
    fprintf(stderr,"SortSparseMatriByIndex: Bad index value %d\n\n",sort_index);
    exit(1);}
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
      exit(1);
    }
  }
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


/* "PdbFileSize" estimates the size of a PDB file, so memory can be allocated. On return, 
   'num_headers' contains the number HEADERs, 'num_models' contains the number of unique MODELs, 
   'num_residues' contains the number of CA atoms, 'num_atoms' contains the total number of atoms, */
void PdbFileSize(char *file,int *num_headers,int *num_models,int *num_chains,int *num_residues,int *num_atoms)
{
  FILE *data;
  char HED[7],ATM[5],calt,chn='\0',c;
  int i;

  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nPdbFileSize: unable to open %s\n\n",file);
    exit(1);}
  (*num_headers) = (*num_models) = (*num_chains) = (*num_residues) = (*num_atoms) = 0;
  for(;;){

    /* Cols 1-6 are field name */
    for(i=0; i<6; i++) HED[i] = getc(data);

    /* Check for ATOM lines */
    if(!strncmp(HED,"ATOM  ",6)){
      for(i=7; i<=12; i++) c = getc(data);            /* Skip atom numbers */
      for(i=13; i<=16; i++) ATM[i-13] = getc(data);   /* Find atom name */
      calt = getc(data);                              /* Find alternative location */
      if(calt=='A' || calt==' '){
	(*num_atoms)++;
	if(strstr(ATM,"CA")!=NULL) (*num_residues)++;
      }
      for(i=18; i<=22; i++) c = getc(data);           /* Skip ahead to chain id */
      if(c!=chn){
	(*num_chains)++;
	chn = c;
      }
    }

    /* Check for lines to be included in the header */
    else if(!strncmp(HED,"HEADER",6) || !strncmp(HED,"TITLE ",6) || 
	    !strncmp(HED,"COMPND",6) || !strncmp(HED,"SEQRES",6) || 
	    !strncmp(HED,"HELIX ",6) || !strncmp(HED,"SHEET ",6) || 
	    !strncmp(HED,"TURN  ",6))
      (*num_headers)++;

    /* Check for models */
    else if(!strncmp(HED,"MODEL ",6)) (*num_models)++;

    else if(!strncmp(HED,"END   ",6) || feof(data))  /* Identify END of file */
      break;
    do{
      c = getc(data);
    }while(c!='\n');
  }
  fclose(data);
  if(*num_models==0) *num_models = 1;
}


/* "AllocatePdbFile" allocates memory for a PDB_File structure: 'HEADER' is a character 
   array with rows (1,num_headers) and columns (0,PDB_MAX_LINE-1); Calpha is an array of 
   Atom_Line structures with elements (1,num_residues). */
void AllocatePdbFile(PDB_File *PDB,int num_headers,int num_residues)
{

  /* Allocate memory for PDB.HEADER */
  /* NOTE (TRL 07/10/2015):  For some reason, cmatrix causes a memory leak if 
     the first row index is 1.  Here zero is used to prevent this. */
  PDB->HEADER = cmatrix(0,num_headers,0,PDB_MAX_LINE-1);

  /* Allocate memory for PDB.atom */
  PDB->atom=malloc((size_t)((num_residues+2)*sizeof(Atom_Line)));
  if(!PDB->atom){
    fprintf(stderr,"\nAllocatePdbFile: fail to allocate atom\n\n");
    exit(1);}
  PDB->num_residues = num_residues;
}


/* "ReadPdbFile" reads HEADER and CA information from a PDB file into a PDB_File object */
int ReadPdbFile(char *file,PDB_File *PDB,int num_headers,int num_residues)
{
  FILE *data;
  char LINE[PDB_MAX_LINE],HED[7],ATM[5],HOLD[9],c,calt;
  int ca,hd,mdl,i;

  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nReadPdbFile: Unable to open %s\n\n",file);
    exit(1);}
  PDB->num_residues = num_residues;
  ca = hd = mdl = 1;
  for(;;){
    i = 0;
    do{                                  /* Skip remaining columns */
      c = getc(data);
      LINE[i++] = c;
    }while(c!='\n' && !feof(data));
    LINE[i] = '\0';
    for(i=0;i<=5;i++) HED[i] = LINE[i];  /* Cols 1-6 are field name */
    HED[6] = '\0';

    /* ---------- Check for ATOM lines ---------- */
    if(!strncmp(HED,"ATOM  ",6)){
      for(i=0; i<=5; i++) PDB->atom[ca].HEAD[i] = LINE[i];
      PDB->atom[ca].HEAD[6] = '\0';

      for(i=12;i<=15;i++) ATM[i-12] = LINE[i];        /* Cols 13-16 are atom name*/
      ATM[4] = '\0';
      calt = LINE[16];                               /* Col 17 is alt. location */

      /* Keep atom if it is alpha carbon */
      if(strstr(ATM,"CA")!=NULL && (calt=='A' || calt==' ')){
	for(i=6; i<=10; i++) HOLD[i-6] = LINE[i];      /* Cols 7-11 are atom # */
	HOLD[5] = '\0';
	sscanf(HOLD,"%d",&PDB->atom[ca].atmnum);
	for(i=12; i<=15; i++) PDB->atom[ca].ATOM[i-12] = LINE[i];    /* Cols 13-16 are atom name*/
	PDB->atom[ca].ATOM[4] = '\0';
	for(i=17; i<=19; i++) PDB->atom[ca].RES[i-17] = LINE[i];   /* Cols 18-20 are res. name */
	PDB->atom[ca].RES[3] = '\0';
	PDB->atom[ca].chain = LINE[21];       /* Col 22 is chain name */

	/* Assign a model number */
	PDB->atom[ca].model = mdl;
	for(i=22; i<=25; i++) HOLD[i-22] = LINE[i];        /* Cols 23-26 are residue # */
	HOLD[4] = '\0';
	sscanf(HOLD,"%d",&PDB->atom[ca].resnum);
	for(i=30; i<=37; i++) HOLD[i-30] = LINE[i];        /* Cols 31-38 are x-coordinate */
	HOLD[8] = '\0';
	sscanf(HOLD,"%f",&PDB->atom[ca].X[0]);
	for(i=38; i<=45; i++) HOLD[i-38] = LINE[i];        /* Cols 39-46 are y-coordinate */
	HOLD[8] = '\0';
	sscanf(HOLD,"%f",&PDB->atom[ca].X[1]);
	for(i=46; i<=53; i++) HOLD[i-46] = LINE[i];        /* Cols 47-54 are z-coordinate */
	HOLD[8] = '\0';
	sscanf(HOLD,"%f",&PDB->atom[ca].X[2]);
	for(i=60; i<=65; i++) HOLD[i-60] = LINE[i];        /* Columns 61-66 are beta */
	HOLD[6] = '\0';
	sscanf(HOLD,"%f",&PDB->atom[ca].beta);             /* Columns 77-78 are element */
	for(i=76; i<=77; i++) PDB->atom[ca].ELEMENT[i-76] = LINE[i]=='\n' ? '\0' : LINE[i];
	PDB->atom[ca].ELEMENT[2] = '\0';
	if(++ca > num_residues){
	  fclose(data);
	  return 0;
	}
      }
    }/* <---- End of 'if(!strncmp(HED,"ATOM  ",6)){... */


    /* ---------- Check for lines to be included in the header ---------- */
    else if((!strncmp(HED,"HEADER",6) || !strncmp(HED,"TITLE ",6) || 
	    !strncmp(HED,"COMPND",6) || !strncmp(HED,"SEQRES",6) ||  
	    !strncmp(HED,"HELIX ",6) || !strncmp(HED,"SHEET ",6) || 
	    !strncmp(HED,"TURN  ",6)) && hd<=num_headers){
      sprintf(PDB->HEADER[hd],"%s",LINE);
      hd++;
    }

    /* ---------- Check for MODEL lines ----------- */
    else if(!strncmp(HED,"ENDMDL",6)) mdl++;

    /* ---------- Check for end of file ---------- */
    else if(!strncmp(HED,"END",3) || feof(data)){
      fclose(data);
      return 0;
    }
  }
}


/* "free_pdb" frees memory from a PDB_File structure. */
void free_pdb(PDB_File *PDB)
{
  free_cmatrix(PDB->HEADER,0,PDB->num_headers,0,PDB_MAX_LINE-1);
  free(PDB->atom);
}


/* "AssignRigidBlocks" assigns each residue to a rigid block and returns the total number of rigid blocks. 
   'pdb' points to an array of 'num_pdbs' PDB_File structures. The array 'blocks' contains 'num_input_blocks' 
   block definitions to be assigned to the residues in the PDB file. Some of these may be overlapping. 
   On return, 'max_block_size' contains the number of residues in the largest block found. Block IDs are 
   assigned beginning with 1 and extending as high as necessary. Unspecified residues are assigned 
   to block 0. */
/* This function replaces 'assign_rigid_blocks2' */
int AssignRigidBlocks(PDB_File *pdb,Rigid_Block *blocks,int num_pdbs,int num_input_blocks,int *max_block_size)
{
  char ch;
  int *unique_blocks,**blockmap,num_unique_blocks,high_block_index,ok,kwn,rr,i,j,k,p;

  /* Make a list of unique block IDs */
  unique_blocks = ivector(1,num_input_blocks);
  for(i=1; i<=num_input_blocks; i++) unique_blocks[i] = 0;
  num_unique_blocks = 0;
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
  blockmap = imatrix(1,num_unique_blocks,1,2);
  high_block_index = kwn = 0;
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
  free_ivector(unique_blocks,1,num_input_blocks);
  free_imatrix(blockmap,1,num_unique_blocks,1,num_unique_blocks);


  /* Find the size of the largest block */
  unique_blocks=ivector(1,high_block_index);
  for(i=1; i<=high_block_index; i++) unique_blocks[i]=0;
  for(p=1; p<=num_pdbs; p++)
    for(i=1; i<=pdb[p].num_residues; i++) 
      if(pdb[p].atom[i].model!=0)
	unique_blocks[pdb[p].atom[i].model]++;
  (*max_block_size) = 0;
  for(i=1; i<=high_block_index; i++)
    if(unique_blocks[i]>(*max_block_size)) *max_block_size = unique_blocks[i];
  free_ivector(unique_blocks,1,high_block_index);

  return high_block_index;
}


/* "CalcProjectionMatrix" calculates the projection matrix from the full residue space to rigid block space. 
   It stores the projection matrix in sparse format in 'projection_matrix' and stores the centers of mass of 
   the blocks in the matrix 'block_center_of_mass'. Returns the number of elements in 'projection_matrix'. */
/* NOTE: Replaces 'dblock_projections3' */
int CalcProjectionMatrix(dSparse_Matrix *projection_matrix,PDB_File *pdb,double **block_center_of_mass,int num_pdbs,int num_blocks,int max_block_size)
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
      for(i=1; i<=pdb[q].num_residues;i++){
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


/* "ReadBlockfile" reads information from a block file. */
void ReadBlockfile(char *file,char ***pdbfile,int *num_pdbs,Rigid_Block **blocks,int *num_blocks)
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

  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nReadBlockfile: Unable to open %s\n\n",file);
    exit(1);}

  /* -------- First pass counts number of PDB files and blocks -------- */
  *num_pdbs = *num_blocks=0;
  while(fgets(buffer,RTB_BUFF_LENGTH,data)!=NULL){
    /* Skip leading whitespace and read the heading */
    i = j = 0;
    while(isspace(buffer[i])) i++;
    while(!isspace(buffer[i]) && buffer[i]!='\0') header[j++] = buffer[i++];
    header[j]='\0';
    if(!strncmp(header,"PDB",strlen(header)) && strlen(header)==3) (*num_pdbs)++;
    else if(!strncmp(header,"BLOCK",strlen(header)) && strlen(header)==5)  (*num_blocks)++;
  }

  /* Allocate memory */
  ptmp = cmatrix(1,(*num_pdbs),0,RTB_WORD_LENGTH);
  pdbid = cmatrix(1,*num_pdbs,0,RTB_WORD_LENGTH);
  temp_blocks = (Rigid_Block *)malloc((size_t)((*num_blocks)+1)*sizeof(Rigid_Block));
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
	fprintf(stderr,"\nread_blockfile: bad line:\n\n%s\n",buffer);
	exit(1);}
      else if(nwd==8){
	j=0;
	while(isspace(buffer[i])) i++;
	while(!isspace(buffer[i]) && buffer[i]!='\0') line[j++] = buffer[i++];
	line[j]='\0';
	for(k=1; k<=(*num_pdbs); k++)
	  if(!strncmp(line,pdbid[k],strlen(line))) temp_blocks[*num_blocks].pdbid = k;
      }
      else temp_blocks[*num_blocks].pdbid=-1;

      /* lores, lochain, lonum */
      j=0;
      while(isspace(buffer[i])) i++;
      while(!isspace(buffer[i]) && buffer[i]!='\0') line[j++] = buffer[i++];
      line[j]='\0';
      strncpy(temp_blocks[*num_blocks].LORES,line,3);
      j=0;
      while(isspace(buffer[i])) i++;
      while(!isspace(buffer[i]) && buffer[i]!='\0') line[j++] = buffer[i++];
      line[j]='\0';
      temp_blocks[*num_blocks].lochain=line[0];
      j=0;
      while(isspace(buffer[i])) i++;
      while(!isspace(buffer[i]) && buffer[i]!='\0') line[j++] = buffer[i++];
      line[j]='\0';
      sscanf(line,"%d",&temp_blocks[*num_blocks].lonum);

      /* hires, hichain, hinum */
      j=0;
      while(isspace(buffer[i])) i++;
      while(!isspace(buffer[i]) && buffer[i]!='\0') line[j++] = buffer[i++];
      line[j]='\0';
      strncpy(temp_blocks[*num_blocks].HIRES,line,3);
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

  free_cmatrix(pdbid,1,(*num_pdbs),0,RTB_WORD_LENGTH);
  fclose(data);
}


/* "ConvertSuperblocksToMatrix" transfers the block Hessian from 
   the tensor 'superblocks' into the array 'block_hessian' */
int ConvertSuperblocksToMatrix(double **block_hessian, double ***superblocks, int **block_contacts, int num_blocks)
{
  int *I1,*I2,i,j,p,sb,ii,jj,max,a,b,imx;
  
  max = 6*num_blocks;
  I1 = ivector(1,max);
  I2 = ivector(1,max);

  /* I1[i]==i iff there is a non-zero element in column i (removes zeroes that are caused by single-node blocks) */
  for(i=1; i<=max; i++){ 
    I1[i] = I2[i] = 0;
    for(j=i; j<=max; j++) block_hessian[i][j] = block_hessian[j][i] = 0.0;
  }
  for(ii=1; ii<=num_blocks; ii++){
    for(i=1; i<=6; i++){
      for(jj=ii; jj<=num_blocks; jj++){
	sb = block_contacts[ii][jj];
	if(sb!=0){
	  p = jj==ii ? i : 1;
	  for(j=p; j<=6; j++)
	    if(superblocks[sb][i][j]!=0)
	      I1[6*(jj-1)+j] = 6*(jj-1) + j;
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
    for(i=1; i<=6; i++){
      for(jj=ii; jj<=num_blocks; jj++){
	sb = block_contacts[ii][jj];
	if(sb!=0){
	  p = jj==ii ? i : 1;
	  for(j=p; j<=6; j++)
	    if(superblocks[sb][i][j]!=0){
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


/* "FindBlockContactsCypa" generates a contact map for rigid blocks, taking CypA binding into account. Block 
   definitions are taken from 'model' values for atoms. On exit, each nonzero element 'block_contacts[i][j]' 
   contains a unique index for the contact between blocks i and j. If blocks are not in contact, 
   block_contacts[i][j]=0. The mass centers of the blocks are stored in 'block_center_of_mass'. In addition 
   to finding contacts between blocks with residues within 'cut' distance from each other, this function 
   assigns contacts to blocks joined by CypA binding. The 'num_blocks'-dimensional square matrix 'cypa_contacts' 
   contains these contacts across blocks. The function returns the total number of block-block contacts, 
   which is also the largest element of 'block_contacts'. */  
int FindBlockContactsCypa(int **block_contacts, PDB_File *PDB, int **blockmap, int *blockmap_index, 
			  double **block_center_of_mass, int num_pdbs, int num_blocks, double cut, 
			  int **chain_starts, int **cypa_contacts, int num_chains)
{
  double *DMAX;
  double cutsq=cut*cut;
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
      if( sqrt(dd) <= DMAX[i] + DMAX[j] + cut ) block_contacts[i][j] = block_contacts[j][i] = 1;
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

