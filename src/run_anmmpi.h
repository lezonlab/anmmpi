#ifndef RUN_ANMMPI_H
#define RUN_ANMMPI_H
#include "readcommandline.h"
extern int rank;
extern int nprocs;
int run_anmmpi(Local_CL_Opt *CL);
#endif
