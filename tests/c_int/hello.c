/* Basic example of libnegf call from C */

#include "libnegf.h"
#include <stdio.h>
#include <string.h>

int main()
{

  int handler[NEGF_HSIZE];
  char realmat[7] = "HR.dat";
  char imagmat[7] = "HI.dat";
  struct lnparams params;
  int surfend[2] = {60,80};
  int contend[2] = {80,100};
  int plend[1] = {60};
  int cblk[2] = {1,1};
  
  printf("Initializing libNEGF \n");
  negf_init_session(handler);
  negf_init(handler);
  negf_read_hs(handler, &realmat[0], &imagmat[0], 0);
  negf_set_s_id(handler, 100);
  negf_init_structure(handler, 2, &contend[0], &surfend[0], 1, &plend[0], &cblk[0]);

  //Set parameters  
  negf_get_params(handler, &params);
  params.emin = -3.0;
  params.emax = 3.0;
  params.estep = 0.01;
  params.wght = 3.0;
  negf_set_params(handler, &params);
  
  //Run calculation and write result to file
  negf_solve_landauer(handler);
  negf_write_tunneling_and_dos(handler);

  //Release library
  negf_destruct_libnegf(handler);
  negf_destruct_session(handler);
  printf("Done \n");

}
