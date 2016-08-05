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
  negf_read_input(handler);
  negf_solve_landauer(handler);
  negf_write_tunneling_and_dos(handler);
  negf_destruct_libnegf(handler);
  negf_destruct_session(handler);
  printf("Done \n");

}
