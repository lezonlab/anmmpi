#ifndef RTBUTIL_HIV_H
#define RTBUTIL_HIV_H
#include "pdbutil_hiv.h"
#include "miscutil_hiv.h"

#define RTB_WORD_LENGTH 256      /* Size of a reasonable file name */
#define RTB_BUFF_LENGTH 1024     /* Buffer size */

extern double g_anm_cutoff;
extern double g_anm_gamma;
extern double g_anm_eta;

typedef struct {
  int blknum; 
  char LORES[4]; 
  char lochain; 
  int lonum; 
  char HIRES[4]; 
  char hichain; 
  int hinum; 
  int pdbid; 
} Rigid_Block;

int ReadBlockfile(const char *file,char ***pdbfile,int *num_pdbs,Rigid_Block **blocks,long *num_blocks);
int AssignRigidBlocks(PDB_File *pdb,Rigid_Block *blocks,int num_pdbs,long num_input_blocks,long *max_block_size);
int CalcProjectionMatrix(dSparse_Matrix *projection_matrix,PDB_File *pdb,double **block_center_of_mass,int num_pdbs,long num_blocks,long max_block_size);
void CondensePrjMtx(dSparse_Matrix *projection_matrix, long num_projection_elements);
int PrintProjectionMatrix(char *outfile, dSparse_Matrix *projection_matrix, long num_projection_elements, long total_residues);
int FindBlockContactsCypa(int **block_contacts, PDB_File *PDB, int **blockmap, int *blockmap_index, 
			  double **block_center_of_mass, int num_pdbs, int num_blocks, 
			  int **chain_starts, int **cypa_contacts, int num_chains);
void AddFullHessRowsToBlockHessCypa(double ***HT, int num_block_contacts, double **hessian_superrow, 
				    int **block_contacts, dSparse_Matrix *projections_sorted_by_1,
				    dSparse_Matrix *projections_sorted_by_2, int num_projection_elements,
				    int **blockmap, int *projection_index_1, int *projection_index_2,
				    int *blockmap_index, int num_blocks, PDB_File *PDB, int num_pdbs, int lores, 
				    int hires, int total_residues, int *KO, int nko, int **cypa_contacts);
void CalcHessianSuperrowCypa(double **hessian_superrow, int **block_contacts, int **blockmap, 
			     int *blockmap_index, int num_blocks, PDB_File *PDB, int num_pdbs, 
			     int row, int p1, int r1, int *KO, int nko, int **cypa_contacts);
int ConvertSuperblocksToMatrix(double **block_hessian, double ***superblocks, int **block_contacts, int num_blocks);
int PrintHessianMatrix(char *outfile, double **hessian, long dim);
int PrintSparseHessianMatrix(char *outfile, double *A, int *RIDX, int *CIDX, long num_hess_elements);
void BlockVecs2Full(dSparse_Matrix *projection_matrix, double **block_vectors, double **full_vectors,
		    long num_projection_elements, long bdim, long fdim, int num_modes);
int ConvertSuperblocksToSparseMatrix(double **A, int **RIDX, int **CIDX, double ***superblocks, int **block_contacts, int num_blocks, long *num_elements);

#endif
