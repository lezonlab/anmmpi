/* This is the pure CLI version */
#include "run_anmmpi.h"
#include <mpi.h>
#include<stdio.h>

int nprocs, rank;


int main(int argc, char *argv[])
{
  Local_CL_Opt params;
  int flag = 1;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if(ReadCommandLine(argc, argv, &params) != 0){
    MPI_Finalize();
    return 1;
  }
  flag = run_anmmpi(&params);
  if(params.nko > 0) free_ivector(params.KO, 1, params.nko);
  MPI_Finalize();
  return 0;
}
