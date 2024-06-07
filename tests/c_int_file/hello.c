/* Basic example of libnegf call from C */

#include "libnegf.h"
#include <stdio.h>
#include <string.h>

int main()
{
  int handler[NEGF_HSIZE];
  
  printf("Initializing libNEGF \n");
  negf_init_session(handler);
  negf_init(handler);
  negf_init_contacts(handler, 2);
  negf_read_input(handler);
  negf_solve_landauer(handler);
  negf_write_tunneling_and_dos(handler);
  negf_destruct_libnegf(handler);
  negf_destruct_session(handler);
  printf("Done \n");

}
