/* Basic example of libnegf call from C */

#include "libnegf.h"
#include <stdio.h>
#include <string.h>

int main()
{

  int handler[NEGF_HSIZE];
  int *hand = &handler[0];
  char realmat[7] = "HR.dat";
  char imagmat[7] = "HI.dat";
  struct lnparams params;
  int surfend[2] = {60,80};
  int contend[2] = {80,100};
  int plend[1] = {60};
  int cblk[2] = {1,1};
  
  printf("Initializing libNEGF \n");
  negf_init_session(hand);
  negf_init(hand);
  negf_read_hs(hand, &realmat[0], &imagmat[0], 0);
  negf_set_s_id(hand, 100);
  negf_init_structure(hand, 2, &contend[0], &surfend[0], 1, &plend[0], &cblk[0]);

  //Set parameters  
  negf_get_params(hand, &params);
  params.emin = -3.0;
  params.emax = 3.0;
  params.estep = 0.01;
  params.wght = 3.0;
  negf_set_params(hand, &params);
  
  //Run calculation and write result to file
  negf_solve_landauer(hand);
  negf_write_tunneling_and_dos(hand);

  //Release library
  negf_destruct_libnegf(handler);
  negf_destruct_session(handler);
  printf("Done \n");

}
