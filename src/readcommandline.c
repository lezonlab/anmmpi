#include "readcommandline.h"
#include "miscutil_hiv.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>


const char *get_filename_ext(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}


/* "ReadCommandLine" reads the command line */
int ReadCommandLine(int argc, char *argv[], Local_CL_Opt *CL)
{
  int ok=0;
  int i;

  /* Allowed flags:
     -c      Specify cutoff distance
     -n      Specify number of modes
     -s      Specify perturbation scaling factor
     -rtb    Assume input is RTB definitions, instead of PDB file
     -ko     Specify knockouts
     -l      Specify CypA binding site (HIV-1 specific)
     -p      Print projection matrix
     -h      Print block Hessian
     --lapack Use LAPACK instead of BLZPACK for decomposition
  */

  if(argc > 1){
    ok = 1;

    CL->infile = argv[1];
    CL->is_blocks = 0;
    CL->num_modes = g_num_modes_default;
    CL->nko = 0;
    CL->print_prj_mtx = 0;
    CL->print_block_hessian = 0;
    CL->print_vectors = 1;
    CL->add_cypa = 0;
    CL->is_blzpack = 1;

    /* Assume first argument is filename */
    if(strcmp(get_filename_ext(CL->infile),"pdb")==0)
      CL->is_blocks = 0;
    else if(strcmp(get_filename_ext(CL->infile),"blk")==0)
      CL->is_blocks = 1;
    else{ 
      printf("Unrecognized file extension: %s\n", get_filename_ext(CL->infile));
      return -1;}

    /* Go through arguments */
    i=2;
    while(i<argc){
      if(strncmp(argv[i],"-c",2)==0){ 
	sscanf(argv[++i],"%lf",&g_anm_cutoff);
	if(g_anm_cutoff < 0.0){
	  printf("\nInvalid cutoff value: %lf\n", g_anm_cutoff);
	  ok=0;}
      }
      else if(strncmp(argv[i],"-n",2)==0){
	sscanf(argv[++i],"%d",&CL->num_modes);
	if(CL->num_modes < 0){
	  printf("\nInvalid number of modes: %d\n", CL->num_modes);
	  ok=0;}
      }
      else if(strncmp(argv[i],"-k",2)==0) sscanf(argv[++i],"%lf",&g_anm_gamma);
      else if(strncmp(argv[i],"-s",2)==0) sscanf(argv[++i], "%lf", &g_anm_eta);
      else if(strncmp(argv[i],"-ko",3)==0){
	CL->nko = nintelemstr(argv[++i]);
	CL->KO = ivector(1, CL->nko);
	hintelemstr(argv[i], CL->KO, CL->nko);
      }
      else if(strncmp(argv[i],"-p",2)==0) CL->print_prj_mtx = 1;
      else if(strncmp(argv[i],"-h",2)==0) CL->print_block_hessian = 1;
      else if(strncmp(argv[i],"-q",2)==0) CL->print_vectors = 0;
      //else if(strncmp(argv[i],"-rtb",4)==0) CL->is_blocks = 1;                   /* Force blocks */
      //else if(strncmp(argv[i],"-l",2)==0){                                       /* CypA binding  */
      //sscanf(argv[++i],"%d",&g_bound_resid);
      //CL->add_cypa = 1;
      //}
      else if(strncmp(argv[i],"--lapack",8)==0) CL->is_blzpack = 0;
      else{
	printf("\n%s: Unknown argument: %s\n\n",argv[0],argv[i]);
	ok=0;
	break;
      }
      i++;	
    }  /* <-- End of 'while(i<argc)...' */
  }

  /* Print error message if things are not all right */
  if(ok!=1){
    printf("\nUsage:\n%s file.dat [OPTIONS]\n\n",argv[0]);
    printf("\tfile.dat -- PDB file\n\n");
    printf("OPTIONS:\n");
    printf("\t-c CUT\tSet ANM cutoff distance to CUT (def %.3f)\n",g_anm_cutoff_default);
    printf("\t-n\tSpecify ANM force constant (default: %.3lf)\n",g_anm_gamma_default);
    printf("\t-n\tSpecify number of modes (default: %d)\n",g_num_modes_default);
    printf("\t-s SCL\tSet perturbation scaling factor to SCL (def %.3f)\n",g_anm_eta_default);
    //printf("\t-rtb\tAssume input is RTB definitions, instead of PDB file\n");
    printf("\t-ko\tSpecify indices of residues to perturb\n");
    printf("\t\tAccepts integers separated by commas (no white space!),\n");
    printf("\t\tas well as ranges indicated with dashes\n");
    //printf("\t-l RES\tBind neighboring chains at RES (HIV-specific)\n");
    printf("\t-p\tPrint projection matrix to file\n");
    printf("\t-h\tPrint block Hessian matrix to file\n");
    printf("\t-q\tSuppress calculation of eigenmodes\n");
    printf("\t--lapack\tDecompose Hessian using LAPACK instead of BLZPACK\n");
    printf("\nOutput:\nCalculates ANM modes.\n\n");
    return -1;
  }
  return 0;
}
