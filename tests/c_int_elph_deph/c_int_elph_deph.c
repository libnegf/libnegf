/* Basic example of calculation of current with dephasing with C interface */

#include "libnegf.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

int main()
{
  struct lnparams params;
  int handler[NEGF_HSIZE];
  int *hand = &handler[0];
  char realmat[7] = "HR.dat";
  char imagmat[7] = "HI.dat";
  int surfstart[2] = {61,81};
  int surfend[2] = {60,80};
  int contend[2] = {80,100};
  int plend[6] = {10, 20, 30, 40, 50, 60};
  int cblk[2] = {6, 1};
  double currents[1];
  double coupling[60];
  int leadpairs;
  int i;
  int foo[1] = {0};

  printf("Initializing libNEGF \n");
  negf_init_session(hand);
  negf_init(hand);
  negf_read_hs(hand, &realmat[0], &imagmat[0], 0);
  negf_set_s_id(hand, 100, 1);
  negf_init_contacts(hand, 2);
  negf_init_structure(hand, 2, &surfstart[0], &surfend[0], &contend[0], 6, &plend[0], &cblk[0]);

  //Set parameters
  negf_get_params(hand, &params);
  params.emin = -2.0;
  params.emax = 2.0;
  params.estep = 0.1;
  params.kbt_t[0] = 0.001;
  params.kbt_t[1] = 0.001;
  params.mu[0] = 0.5;
  params.mu[1] = -0.5;
  params.verbose = 100;
  negf_set_params(hand, &params);

  // Set the dephasing model as simple diagonal dephasing.
  for (i = 0; i < 60; ++i)
  {
    coupling[i] = 0.05;
  }
  negf_set_elph_dephasing(hand, &coupling[0], 60, 10);

  // Calculate the current.
  negf_solve_landauer(hand);
  negf_get_currents(hand, &leadpairs, &currents[0], 1);

  printf("Current: %f\n",currents[0]);
  //Release library
  negf_destruct_libnegf(handler);
  negf_destruct_session(handler);
  negf_mem_stats(handler);
  printf("Done \n");

  if (currents[0] > 1.95 || currents[0] < 1.90) {
    printf("Error in current value not between 1.90 and 1.95 \n");
    return 1;
  }
  else
  {
    return 0;
  }  
}
