module mpi_globals

#ifdef MPI
  
  include 'mpif.h'                                   !#MPI#
  
  ! MPI COMMON VARIABLES
  INTEGER, SAVE ::  mpi_comm, id, numprocs, ierr
  LOGICAL, SAVE ::  id0
  
#else
  
  logical, parameter ::  id0 = .true.
  integer, parameter :: id = 0, numprocs = 1
         
#endif


end module mpi_globals
