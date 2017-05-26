#ifndef READCOMMANDLINE_H
#define READCOMMANDLINE_H

extern const int g_bound_resid_default;
extern const double g_anm_cutoff_default;
extern const double g_anm_gamma_default;
extern const int g_num_modes_default;
extern const double g_anm_eta_default;
extern double g_anm_cutoff;
extern double g_anm_eta;
extern int g_bound_resid;
extern double g_anm_gamma;

/* Structure to hold command-line arguments */
typedef struct {
  const char *infile; 
  int is_blocks;
  int num_modes; 
  int *KO; 
  int nko; 
  int print_prj_mtx; 
  int print_block_hessian;
  int print_vectors;
  int add_cypa;
  int is_blzpack;
} Local_CL_Opt;

int ReadCommandLine(int argc, char *argv[], Local_CL_Opt *CL);

#endif
