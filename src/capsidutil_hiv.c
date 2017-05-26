#include "capsidutil_hiv.h"
#include "miscutil_hiv.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>


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
    return -1;}

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
int AssignCypaContacts(PDB_File *pdb, int num_pdbs, int **chain_starts, int **cypa_contacts, 
		       int num_chains, int num_blocks, int *num_mers)
{
  double **chain_directions;
  double **mer_centers;
  double mer_vector,cypa_vector;
  double loop_loop_dist,mer_mer_dist;
  double chain_chain_angle,mer_cypa_angle,chain_cypa_angle;
  int *num_bound;
  int *block_for_chain;
  int num_contacts = 0;
  int i,j,k,p;

  chain_directions = dmatrix(1, num_chains, 0, 2);
  mer_centers = dmatrix(1, num_chains, 0, 2);

  *num_mers = AssignChainStarts(pdb, num_pdbs, chain_starts, chain_directions, mer_centers, num_chains);

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
  free_dmatrix(chain_directions, 1, num_chains, 0, 2);
  free_dmatrix(mer_centers, 1, num_chains, 0, 2);
  free_ivector(num_bound, 1, num_chains);
  free_ivector(block_for_chain, 1, num_chains);
  return num_contacts;
}


/* "AssignCypaContacts" populates the matrix cypa_contacts[1..num_blocks][1..num_blocks] with integers 
   indicating which blocks are joined by CypA molecules. It returns the number of inter-CA CypA contacts 
   formed. The CypA binding loop is residues 88-90. Here it will be represented by one residue only. */
int AssignCypaContacts_old(PDB_File *pdb, int num_pdbs, int **chain_starts, double **chain_directions,
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

