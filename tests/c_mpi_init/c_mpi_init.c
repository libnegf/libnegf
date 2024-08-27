
#include "libnegf.h"
#include <stdio.h>
#include <string.h>
#include <mpi.h>

int main()
{
  int handler[NEGF_HSIZE];
  int *hand = &handler[0];
  int ierr;

  ierr = MPI_Init(NULL, NULL);
  if (ierr != 0)
  {
     printf("Error in mpi_init \n");
     return 1;  
  }

  printf("Initializing libNEGF \n");
  negf_init_session(hand);
  negf_init(hand);

  MPI_Fint global_comm_f = MPI_Comm_c2f(MPI_COMM_WORLD);
  negf_set_mpi_fcomm(hand, global_comm_f);
  MPI_Fint cart_comm, k_comm, en_comm;
  negf_cartesian_init(hand, global_comm_f, 1, &cart_comm, &k_comm, &en_comm);

  //Release library
  negf_destruct_libnegf(hand);
  negf_destruct_session(hand);
  printf("Done \n");

  MPI_Finalize();

  return 0;
}
