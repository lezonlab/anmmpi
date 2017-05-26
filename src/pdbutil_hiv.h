#ifndef PDBUTIL_HIV_H
#define PDBUTIL_HIV_H

#define PDB_MAX_LINE 90

typedef struct {
  char HEAD[7]; 
  int atmnum; 
  char ATOM[5]; 
  char RES[4]; 
  char chain; 
  int resnum; 
  float X[3]; 
  float beta; 
  char ELEMENT[3]; 
  int model;
} Atom_Line;

typedef struct {
  char **HEADER; 
  Atom_Line *atom; 
  int num_residues; 
  int num_headers;
} PDB_File;

int AllocatePdbFile(PDB_File *PDB, long lo_header, long hi_header, long lo_residue, long high_residue);
void free_pdb(PDB_File *PDB, long lo_header, long lo_residue);
int PdbFileSize(const char *file,int *num_headers,int *num_models,int *num_chains,long *num_ca,long *num_atoms);
int ReadPdbFile(const char *file, PDB_File *PDB, int num_headers, long num_residues);
PDB_File *ReadMultiplePdbFiles(char **filenames, int num_pdbs, long *total_chains, long *total_residues);
int fprintPerturbVecsMulti(PDB_File *PDB, char **pdbfilenames, double **VEC,
			    int num_pdbs, int fdim, int num_modes, int knockout);

#endif
