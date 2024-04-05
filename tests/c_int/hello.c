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
  int surfstart[2] = {61,81};
  int surfend[2] = {60,80};
  int contend[2] = {80,100};
  int ldos_start[1] = {1};
  int ldos_end[1] = {1};
  int plend[1] = {60};
  int cblk[2] = {1,1};
  double *tr;
  int tr_shape[2];

  printf("Initializing libNEGF \n");
  negf_init_session(hand);
  negf_init(hand);
  negf_read_hs(hand, &realmat[0], &imagmat[0], 0);
  negf_set_s_id(hand, 100, 1);
  negf_init_contacts(hand, 2);
  negf_init_structure(hand, 2, &surfstart[0], &surfend[0], &contend[0], 1, &plend[0], &cblk[0]);

  //Set parameters
  negf_get_params(hand, &params);
  params.emin = -3.0;
  params.emax = 3.0;
  params.estep = 0.01;
  params.kbt_t[0] = 0.001;
  params.kbt_t[1] = 0.001;
  //use default for Np_n = (20,20); n_kt=10; n_poles=3; 
  // delta = 1e-4; dos_delta = 1e-4; g_spin = 2;
  negf_set_params(hand, &params);

  //Set ldos parameters
  negf_init_ldos(hand, 1);
  negf_set_ldos_intervals(hand, 1, &ldos_start[0], &ldos_end[0]);

  //Run calculation and write result to file
  //Also get a reference to transmission, just to try
  negf_solve_landauer(hand);
  negf_associate_transmission(hand, tr_shape, &tr);
  printf("shape %d %d \n",tr_shape[0], tr_shape[1]);
  printf("tr %e %e \n",tr[0], tr[1]);
  negf_write_tunneling_and_dos(hand);

  //Release library
  negf_destruct_libnegf(handler);
  negf_destruct_session(handler);
  printf("Done \n");

}
