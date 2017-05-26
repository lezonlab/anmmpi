#include "run_anmmpi.h"
#include "pdbutil_hiv.h"
#include "miscutil_hiv.h"
#include "rtbutil_hiv.h"
#include "blzsolvesparse.h"
#include "lapsolve.h"
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<mpi.h>

const int g_bound_resid_default = 89;
const double g_anm_cutoff_default = 15.0;
const double g_anm_gamma_default = 1.0;
const int g_num_modes_default = 20;
const double g_anm_eta_default = 10.0;

double g_anm_cutoff = 15.0;
double g_anm_gamma = 1.0;
double g_anm_eta = 10.0;
int g_bound_resid = 89;
int nprocs = 1;
int rank = 0;

char *NameOutputFile(const Local_CL_Opt *CL)
{
  char *outfile;
  char *bunk;

  outfile = calloc( (strlen(CL->infile)+99),sizeof(char) );
  bunk = calloc(99, sizeof(char));
  strncpy( outfile, CL->infile, strlen(CL->infile)-strlen(strrchr(CL->infile,'.')) );
  if (CL->nko==1){
    sprintf(bunk,"-%03d",CL->KO[1]);
    strcat(outfile,bunk);
  }
  strcat(outfile,".something");
  free(bunk);
  return outfile;
}


int run_anmmpi(Local_CL_Opt *CL)
{
  Rigid_Block *blocks;
  PDB_File *PDB;
  dSparse_Matrix projection_matrix;
  dSparse_Matrix *projections_sorted_by_1,*projections_sorted_by_2;
  char **pdbfilenames;
  double ***superblocks_partial;
  double **hessian_superrow;
  double **block_center_of_mass;
  int **block_contacts;
  int **blockmap;
  int *blockmap_index;
  int *projection_index_1;
  int *projection_index_2;
  int **chain_starts;
  int **cypa_contacts;
  int num_mers;
  int num_bound_cypas = 0;
  int num_block_contacts = 0;
  long total_residues = 0;
  long total_chains = 0;
  long num_blocks = 0;
  long num_input_blocks = 0;
  long max_block_size = 0;
  long num_projection_elements = 0;
  long lores = 1;
  long hires = 0;
  long avgload = 0;
  int num_pdbs = 1;
  int i,j,k,p,f;

  /* These are only required for master process */
  FILE *data;
  char *outfile;
  double ***superblocks;
  double **HH;
  double **VEC;
  double **WVEC;
  double *WVAL;
  double *A;
  int *RIDX, *CIDX;
  long fdim;
  long bdim = 0;
  long lo_mode=7, hi_mode, num_modes;
  long num_hess_elements = 0;
  int knockout;
  int nvbset = 0; /* Tells BLZPACK to automatically group vectors */

  /* MPI stuff */
  MPI_Status status;
  double t0,tdiff,tlast;
  int send_tag=0;


  if(rank==0){ 
    t0 = MPI_Wtime();
    printf("\n ***** run_anm *****\n");
  }
  if(CL->is_blocks==1){
    ReadBlockfile(CL->infile, &pdbfilenames, &num_pdbs, &blocks, &num_input_blocks);
    if(rank==0){
      nvbset = 6;
      printf("Blockfile: %s\n",CL->infile);
      printf("%ld blocks\n",num_input_blocks);
      printf("%d PDB files:\n",num_pdbs);
      for(f=1; f<=num_pdbs; f++) printf("\t'%s'\n", pdbfilenames[f]);
    }
  }
  else{
    if(rank==0) printf("PDB file: '%s'\n",CL->infile);
    pdbfilenames = cmatrix(1, 1, 0, RTB_WORD_LENGTH);
    strncpy(pdbfilenames[1], CL->infile, RTB_WORD_LENGTH);
    num_pdbs = 1;
  }
  if(rank==0) printf("Cutoff: %.3f\n",g_anm_cutoff);


  /* Name the output file */
  if(rank==0) outfile = NameOutputFile(CL);
  
  /* Read the PDB files: PDB[1, num_pdbs] */
  PDB = ReadMultiplePdbFiles(pdbfilenames, num_pdbs, &total_chains, &total_residues);
  if(PDB==NULL){
    if(rank==0) fprintf(stderr,"\nERROR: Failed to read PDB files\n");
    return -1;}

  /* Assign each residue to a block.  Blocks are given sequential indices starting at 1.  */
  if(CL->is_blocks==1){
    num_blocks = AssignRigidBlocks(PDB, blocks, num_pdbs, num_input_blocks, &max_block_size);
    free(blocks);}
  else{
    num_blocks = PDB[1].num_residues;
    max_block_size = 1;
    for(i=1; i<=PDB[1].num_residues; i++) PDB[1].atom[i].model = i;
  }
  if(rank==0){
    printf("Total number of residues: %ld\n",total_residues);
    printf("%ld total rigid blocks in structure\n",num_blocks);
    printf("Largest block contains %ld residue%c\n",max_block_size, max_block_size==1 ? ' ' : 's');
  }
  /* TODO: If max_block_size==1, simply run ANM without RTB */


  /* Find the projection matrix.  The 3 coordinates of each residue project onto the 6 coordinates of 
     exactly one block.  The translational component of the projection will always be a diagonal 
     matrix, so there are 12 required components of the projection matrix for each residue.  The 
     coordinates of the block mass centers are stored in block_center_of_mass for later use. */
  /* TODO: Stuff this into a function */
  projection_matrix.IDX = imatrix(1, 12*total_residues, 1, 2);
  projection_matrix.X = dvector(1, 12*total_residues);
  block_center_of_mass = dmatrix(1, num_blocks, 1, 3);
  num_projection_elements = CalcProjectionMatrix(&projection_matrix, PDB, block_center_of_mass, num_pdbs, 
						 num_blocks, max_block_size);
  if( num_projection_elements > 12*total_residues ){
    if(rank==0){
      fprintf(stderr,"\nERROR: Unexpected size for projection matrix.\n");
      fprintf(stderr,"This may result from an unassigned residue.\n");}
    return -1;}
  if(rank==0) printf("%ld non-zero elements in the projection matrix\n",num_projection_elements);
  if(SortSparseMatrixByIndex(&projection_matrix, num_projection_elements, 1) != 0){
    if(rank==0) fprintf(stderr,"\nERROR: Sorting error\n");
    return -1;}


  /* Print the projection matrix to file if requested. If we're only printing 
     the projection matrix, we can stop here; otherwise, carry on.  */
  if(CL->print_prj_mtx==1){
    if(rank==0){
      outfile[(int)(strlen(outfile)-strlen(strrchr(outfile,'.')))+1]='\0';
      strcat(outfile,"prj");
      if(PrintProjectionMatrix(outfile, &projection_matrix, num_projection_elements, total_residues) != 0){
	fprintf(stderr,"\nERROR: Failed to print projection matrix\n");
	return -1;}
    }
    if(CL->print_block_hessian==0 && CL->print_vectors==0){
      if(rank==0) free(outfile);
      free_cmatrix(pdbfilenames, 1, 1, 0, RTB_WORD_LENGTH);
      for(f=1; f<=num_pdbs; f++) free_pdb(&(PDB[f]), 1, 1);
      free(PDB+1);
      free_imatrix(projection_matrix.IDX, 1, 12*total_residues, 1, 2);
      free_dvector(projection_matrix.X, 1, 12*total_residues);
      free_dmatrix(block_center_of_mass, 1, num_blocks, 1, 3);
      return 0;
    }
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
  /* Assign a unique chain ID to each chain. These will be sequential: 1...total_chains */
  /* Chain 'i' starts at PDB[chain_starts[i][1]].atom[chain_starts[i][2]] */
  /* Chain 'i' is in mer chain_starts[i][3] */
  chain_starts = imatrix(1, total_chains, 1, 3);
  cypa_contacts = imatrix(1, num_blocks, 1, num_blocks);
  if(CL->add_cypa==1)
    num_bound_cypas = AssignCypaContacts(PDB, num_pdbs, chain_starts, cypa_contacts, total_chains, num_blocks, &num_mers);
  else
    for(i=1; i<=num_blocks; i++)
      for(j=1; j<=num_blocks; j++)
	cypa_contacts[i][j] = 0;
  if(num_bound_cypas!=0 && rank==0) printf("%d CypA molecules bound\n", num_bound_cypas);

  /* -------------------- Calculate the block Hessian ------------------------ */
  /*  blockmap is a mapping between blocks and residues:
      blockmap[i][1] is global residue index; blockmap[i][2] is PDB index; 
      blockmap[i][3] is resnum; blockmap[i][4] is block 
      blockmap gets sorted by the block index, then blockmap_index stores the beginning
      positions of each of the blocks in blockmap: 
      Block 'i' is rows blockmap_index[i] to blockmap_index[i+1]-1 of blockmap.
  */
  block_contacts = imatrix(0, num_blocks, 0, num_blocks);
  blockmap = imatrix(1, total_residues, 1, 4);
  blockmap_index = ivector(1, num_blocks+1);
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
  num_block_contacts = FindBlockContactsCypa(block_contacts, PDB, blockmap, blockmap_index, 
					     block_center_of_mass, num_pdbs, num_blocks, 
					     chain_starts, cypa_contacts, total_chains);


  /*
    Make 2 copies of the projection matrix: One sorted by each index. Mark where indices change.
    projection_index_1 and projection_index_2 are defined such that 
    for all j: projection_index_1[i] <= j < projection_index_1[i+1], projections_sorted_by_1->IDX[j][1] = i 
    for all j: projection_index_2[i] <= j < projection_index_2[i+1], projections_sorted_by_2->IDX[j][2] = i
    TODO: Eliminate the extra allocation and matrix copy. Do it with pointers, instead.
  */
  projections_sorted_by_1 = &projection_matrix;
  projections_sorted_by_2 = malloc(sizeof(dSparse_Matrix));
  projections_sorted_by_2->IDX = imatrix(1, num_projection_elements, 1, 2);
  projections_sorted_by_2->X = dvector(1, num_projection_elements);
  projection_index_1 = ivector(1, 3*total_residues+1);
  projection_index_2 = ivector(1, 6*num_blocks+1);
  CopySparseMatrix(projections_sorted_by_1, projections_sorted_by_2, 1, num_projection_elements);
  if(SortSparseMatrixByIndex(projections_sorted_by_2, num_projection_elements, 2) != 0){
    if(rank==0) fprintf(stderr,"\nERROR: Sorting error\n");
    return -1;}
  AssignIndexStretches(projection_index_1, projections_sorted_by_1, 
		       num_projection_elements, 3*total_residues+1, 1);
  AssignIndexStretches(projection_index_2, projections_sorted_by_2, 
		       num_projection_elements, 6*num_blocks+1, 2);


  /* ---------------------- Parallel --------------------- */
  /* Hessian calculation */
  if(rank==0){
    printf("Calculating block Hessian...");
    tlast = MPI_Wtime();}
  superblocks_partial = d3tensor(1, num_block_contacts, 1, 6, 1, 6);
  hessian_superrow = dmatrix(1, 3*total_residues, 1, 3);
  avgload = total_residues/nprocs;
  lores = rank*avgload + 1;
  hires = (rank + 1)*avgload;
  if (rank==nprocs - 1) hires = total_residues;
  AddFullHessRowsToBlockHessCypa(superblocks_partial, num_block_contacts, hessian_superrow, block_contacts,
				 projections_sorted_by_1, projections_sorted_by_2, num_projection_elements,
				 blockmap, projection_index_1, projection_index_2, blockmap_index, num_blocks,
				 PDB, num_pdbs, lores, hires, total_residues, CL->KO, CL->nko, cypa_contacts);

  /* Divide the matrix calculation among processes. */
  if(rank==0){
    superblocks = d3tensor(1, num_block_contacts, 1, 6, 1, 6);
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
    if(CL->is_blzpack==0){
      HH = dmatrix(1, 6*num_blocks, 1, 6*num_blocks);
      bdim = ConvertSuperblocksToMatrix(HH, superblocks, block_contacts, num_blocks);
    }
    else{
      bdim = ConvertSuperblocksToSparseMatrix(&A, &RIDX, &CIDX, superblocks, block_contacts, num_blocks, &num_hess_elements);
      printf("%ld non-zero elements in the block Hessian\n", num_hess_elements);
    }
    free_d3tensor(superblocks, 1, num_block_contacts, 1, 6, 1, 6);
    tdiff = MPI_Wtime() - tlast;
    printf("Finished in %e seconds\n",tdiff);
  }
  free_imatrix(block_contacts, 0, num_blocks, 0, num_blocks);
  free_imatrix(blockmap, 1, total_residues, 1, 4);
  free_ivector(blockmap_index, 1, num_blocks+1);
  free_imatrix(projections_sorted_by_2->IDX, 1, num_projection_elements, 1, 2);
  free_dvector(projections_sorted_by_2->X, 1, num_projection_elements);
  free(projections_sorted_by_2);
  free_ivector(projection_index_1, 1, 3*total_residues+1);
  free_ivector(projection_index_2, 1, 6*num_blocks+1);
  free_dmatrix(hessian_superrow, 1, 3*total_residues, 1, 3);
  free_d3tensor(superblocks_partial, 1, num_block_contacts, 1, 6, 1, 6);
  /* ------------------------------- End of parallel part --------------------------------- */

  /* NOTE: If we're only printing the Hessian (not calculating modes), we can stop here; otherwise,
     carry on. */
  /* Print block Hessian to file */
  if(rank==0){
    if(CL->print_block_hessian==1){
      outfile[(int)(strlen(outfile)-strlen(strrchr(outfile,'.')))+1]='\0';
      if(CL->is_blocks==1){
	strcat(outfile,"blockhessian");
	printf("Writing block Hessian to %s\n",outfile);
      }
      else{
	strcat(outfile,"sparsehessian");
	printf("Writing sparse Hessian to %s\n",outfile);
      }

      if(CL->is_blzpack==0){
	if(PrintHessianMatrix(outfile, HH, bdim) !=0 ){
	  fprintf(stderr,"\nERROR: Unable to print Hessian matrix\n");
	  return -1;}
      }
      else{
	if(PrintSparseHessianMatrix(outfile, A, RIDX, CIDX, num_hess_elements) != 0){
	  fprintf(stderr,"\nERROR: Unable to print Hessian matrix\n");
	  return -1;}
      }

    }
    if(CL->print_vectors==1){
      num_modes = CL->num_modes;
      if(lo_mode + num_modes - 1 > bdim) num_modes = bdim - lo_mode + 1;
      hi_mode = lo_mode + num_modes - 1;
      WVEC = dmatrix(1, bdim, 1, num_modes);
      WVAL = dvector(1, num_modes);

      /* Decompose */
      /* LAPACK */
      if(CL->is_blzpack==0){
	printf("Calculating %ld modes with LAPACK\n",num_modes);
	if(lapsolve(HH, WVAL, WVEC, 6*num_blocks, bdim, lo_mode, hi_mode) != 0){
	  fprintf(stderr,"\nERROR: Failed to decompose Hessian\n");
	  return -1;}
      }
      else{
	/* BLZPACK */
	printf("Calculating %ld modes with BLZPACK\n",num_modes);
	printf("lo_mode = %ld\nhi_mode = %ld\nbdim = %ld\n",lo_mode,hi_mode,bdim);
	if(blzSolveSparse(A, RIDX, CIDX, num_hess_elements, WVAL, WVEC, bdim, lo_mode, hi_mode, nvbset) != 0){
	  fprintf(stderr,"\nERROR: Failed to decompose Hessian\n");
	  return -1;}
      }
      tdiff = MPI_Wtime() - tlast;
      printf("Finished in %e seconds\n",tdiff);

      /* Back-project modes into full space */
      fdim = 3*total_residues;
      VEC = dmatrix(1, fdim, 1, num_modes);
      BlockVecs2Full(&projection_matrix, WVEC, VEC, num_projection_elements, bdim, fdim, num_modes);
      outfile[(int)(strlen(outfile)-strlen(strrchr(outfile,'.')))+1]='\0';
      strcat(outfile,"val");
      printf("Writing eigenvalues to %s\n",outfile);

      /* This is not implemented yet... */
      //PrintEigenSystemMulti(outfile, &projection_matrix, WVEC, WVAL, num_projection_elements, bdim, fdim, num_modes);
    
      if((data=fopen(outfile,"w"))==NULL){
	printf("\nUnable to open %s\n\n",outfile);
	return -1;}
      for(i=1; i<=num_modes; i++) fprintf(data,"% 16.7e\n",WVAL[i]);
      fclose(data);
    
      /* Print the eigenvectors to separate files if there are multiple PDBs specified in the 
	 block file, or to a single file in the current directory otherwise. */
      knockout = CL->nko==0 ? 0 : CL->KO[1];
      if(fprintPerturbVecsMulti(PDB, pdbfilenames, VEC, num_pdbs, fdim, num_modes, knockout) != 0){
	fprintf(stderr,"\nERROR: Failed to print modes\n\n");
	return -1;}

      free_dmatrix(VEC, 1, fdim, 1, num_modes);
      free_dmatrix(WVEC, 1, bdim, 1, num_modes);
      free_dvector(WVAL, 1, num_modes);
    }

    free(outfile);
    if(CL->is_blzpack==0) free_dmatrix(HH, 1, 6*num_blocks, 1, 6*num_blocks);
    else{
      free_dvector(A, 0, num_hess_elements-1);
      free_ivector(RIDX, 0, num_hess_elements-1);
      free_ivector(CIDX, 0, num_hess_elements-1);}
  }


  /* Clean up */
  free_cmatrix(pdbfilenames, 1, 1, 0, RTB_WORD_LENGTH);
  for(f=1; f<=num_pdbs; f++) free_pdb(&(PDB[f]), 1, 1);
  free(PDB+1);
  free_imatrix(projection_matrix.IDX, 1, 12*total_residues, 1, 2);
  free_dvector(projection_matrix.X, 1, 12*total_residues);
  free_dmatrix(block_center_of_mass, 1, num_blocks, 1, 3);
  free_imatrix(chain_starts, 1, total_chains, 1, 3);
  free_imatrix(cypa_contacts, 1, num_blocks, 1, num_blocks);
  if(rank==0){
    tdiff = MPI_Wtime() - t0;
    printf("Elapsed time: %e seconds\n",tdiff);
  }  
  return 0;
}
