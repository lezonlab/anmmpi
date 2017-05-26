#ifndef CAPSIDUTIL_HIV_H
#define CAPSIDUTIL_HIV_H
#include "pdbutil_hiv.h"

extern int g_bound_resid;

int AssignChainStarts(PDB_File *pdb, int num_pdbs, int **chain_starts, double **chain_directions, 
		      double **mer_centers, int num_chains);
int AssignCypaContacts(PDB_File *pdb, int num_pdbs, int **chain_starts, int **cypa_contacts, 
		       int num_chains, int num_blocks, int *num_mers);

int AssignCypaContacts_old(PDB_File *pdb, int num_pdbs, int **chain_starts, double **chain_directions, 
		       int **cypa_contacts, int num_chains, int num_blocks, double **mer_centers, int num_mers);
#endif
