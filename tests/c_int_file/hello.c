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
  negf_read_input(hand);
  negf_solve_landauer(hand);
  negf_write_tunneling_and_dos(hand);
  negf_destruct_libnegf(hand);
  negf_destruct_session(hand);
  printf("Done \n");

}
