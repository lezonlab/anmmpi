#include "pdbutil_hiv.h"
#include "miscutil_hiv.h"
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>


/* "AllocatePdbFile" allocates memory for a PDB_File structure: 'HEADER' is a character 
   array with rows (lo_header, hi_header) and columns (0,PDB_MAX_LINE); Calpha is an array of 
   Atom_Line structures with elements (lo_residue, high_residue). */
int AllocatePdbFile(PDB_File *PDB, long lo_header, long hi_header, long lo_residue, long hi_residue)
{
  PDB->num_residues = hi_residue - lo_residue + 1;
  PDB->num_headers = hi_header - lo_header + 1;

  /* Allocate memory for PDB.HEADER */
  PDB->HEADER = cmatrix(lo_header, lo_header+PDB->num_headers-1, 0, PDB_MAX_LINE);

  /* Allocate memory for PDB.atom */
  PDB->atom = malloc( (PDB->num_residues)*sizeof(Atom_Line) );
  if(!PDB->atom){
    fprintf(stderr,"\nERROR: Failure to allocate PDB file\n\n");
    return -1;}
  PDB->atom-=lo_residue;
  return 0;
}


/* "free_pdb" frees memory from a PDB_File structure. */
void free_pdb(PDB_File *PDB, long lo_header, long lo_residue)
{
  free_cmatrix(PDB->HEADER, lo_header, lo_header+PDB->num_headers-1, 0, PDB_MAX_LINE);
  free(PDB->atom+lo_residue);
}


/* "PdbFileSize" estimates the size of a PDB file, so memory can be allocated. On return, 
   'num_headers' contains the number HEADERs, 'num_models' contains the number of unique MODELs, 
   'num_residues' contains the number of CA atoms, 'num_atoms' contains the total number of atoms, */
int PdbFileSize(const char *file,int *num_headers,int *num_models,int *num_chains,long *num_residues,long *num_atoms)
{
  FILE *data;
  char HED[7],ATM[5],calt,chn='\0',c;
  int i;

  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nERROR: Unable to open '%s'\n\n",file);
    return -1;}
  (*num_headers) = (*num_models) = (*num_chains) = (*num_residues) = (*num_atoms) = 0;
  for(;;){

    /* Cols 1-6 are field name */
    for(i=0; i<6; i++) HED[i] = getc(data);

    /* Check for ATOM lines */
    if(!strncmp(HED,"ATOM  ",6)){
      for(i=7; i<=12; i++) c = getc(data);            /* Skip atom numbers */
      for(i=13; i<=16; i++) ATM[i-13] = getc(data);   /* Find atom name */
      ATM[4] = '\0';
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
  return 0;
}


/* "ReadPdbFile" reads HEADER and CA information from a PDB file into a PDB_File object */
int ReadPdbFile(const char *file, PDB_File *PDB, int num_headers, long num_residues)
{
  FILE *data;
  char LINE[PDB_MAX_LINE],HED[7],ATM[5],HOLD[9],c,calt;
  int ca,hd,mdl,i;

  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nERROR: Unable to open %s\n\n",file);
    return 1;}
  PDB->num_residues = num_residues;
  ca = hd = mdl = 1;

  for(;;){
    i = 0;
    do{                                  /* Skip remaining columns */
      c = getc(data);
      LINE[i++] = c;
    }while(c!='\n' && !feof(data));
    LINE[i] = '\0';
    for(i=0; i<=5; i++) HED[i] = LINE[i];  /* Cols 1-6 are field name */
    HED[6] = '\0';

    /* ---------- Check for ATOM lines ---------- */
    if(!strncmp(HED,"ATOM  ",6)){
      for(i=0; i<=5; i++) PDB->atom[ca].HEAD[i] = LINE[i];
      PDB->atom[ca].HEAD[6] = '\0';

      for(i=12; i<=15; i++) ATM[i-12] = LINE[i];        /* Cols 13-16 are atom name */
      ATM[4] = '\0';
      calt = LINE[16];                               /* Col 17 is alt. location */

      /* Keep atom if it is alpha carbon */
      if(strstr(ATM,"CA")!=NULL && (calt=='A' || calt==' ')){
	for(i=6; i<=10; i++) HOLD[i-6] = LINE[i];      /* Cols 7-11 are atom # */
	HOLD[5] = '\0';
	sscanf(HOLD,"%d",&PDB->atom[ca].atmnum);
	for(i=12; i<=15; i++) PDB->atom[ca].ATOM[i-12] = LINE[i];    /* Cols 13-16 are atom name */
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


PDB_File *ReadMultiplePdbFiles(char **filenames, int num_pdbs, long *total_chains, long *total_residues)
{
  PDB_File *PDB;
  long num_residues = 0;
  long num_atoms = 0;
  int num_headers = 0;
  int num_models = 1;
  int num_chains = 0;
  int f;

  *total_chains = 0;
  PDB = malloc( (num_pdbs)*sizeof(PDB_File));
  if(!PDB){
    fprintf(stderr,"\nERROR: Allocation error in ReadMultiplePdbFiles\n\n");
    return;}
  PDB--;
  for(f=1; f<=num_pdbs; f++){
    if(PdbFileSize(filenames[f], &num_headers, &num_models, &num_chains, &num_residues, &num_atoms)!=0)
      return NULL;
    (*total_chains) += num_chains;
    (*total_residues) += num_residues;
    if(num_models!=1){
      fprintf(stderr,"\n%s: %d models (needs to be one)!\n\n", filenames[f], num_models);
      return NULL;}
    if(AllocatePdbFile(&(PDB[f]), 1, num_headers, 1, num_residues) != 0) return NULL;
    if(ReadPdbFile(filenames[f], &PDB[f], 0, num_residues) != 0) return NULL;
  }
  return PDB;
}


/* "fprintPerturbVecsMulti" breaks up the input eigenvectors 
   by PDB and spreads them among multiple files.  It handles 
   single-residue perturbations. */
int fprintPerturbVecsMulti(PDB_File *PDB, char **pdbfilenames, double **VEC,
			    int num_pdbs, int fdim, int num_modes, int knockout)
{
  FILE *data;
  char *vecfile;
  char *bunk;
  int row=0;
  int f,i,j,q;

  vecfile = (char *)calloc(256,sizeof(char));
  bunk = (char *)calloc(99,sizeof(char));
  for(f=1; f<=num_pdbs; f++){
    strncpy(vecfile, pdbfilenames[f], strlen(pdbfilenames[f])-strlen(strrchr(pdbfilenames[f],'.')));
    vecfile[(int)(strlen(pdbfilenames[f])-strlen(strrchr(pdbfilenames[f],'.')))]='\0';
    sprintf(bunk,"-%03d",knockout);
    strcat(vecfile,bunk);
    strcat(vecfile,".vec");
    printf("%s\n",vecfile);
    if((data=fopen(vecfile,"w"))==NULL){
      fprintf(stderr,"\nERROR: Unable to open %s\n\n",vecfile);
      return -1;}
    for(i=1; i<=PDB[f].num_residues; i++){
      for(q=1; q<=3; q++){
	row++;
	for(j=1; j<num_modes; j++)
	  fprintf(data,"%1.7e\t",PDB[f].atom[i].model==0 ? 0.0 : VEC[row][j]);
	fprintf(data,"%1.7e\n",PDB[f].atom[i].model==0 ? 0.0 : VEC[row][num_modes]);
      }
    }
    fclose(data);
  }
  free(vecfile);
  free(bunk);
  return 0;
}
