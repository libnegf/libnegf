!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Green's functions.        !
!! Copyright (C) 2012                                                       !
!!                                                                          ! 
!! This file is part of libNEGF: a library for                              !
!! Non Equilibrium Green's Function calculation                             ! 
!!                                                                          !
!! Developers: Alessandro Pecchia, Gabriele Penazzi                         !
!! Former Conctributors: Luca Latessa, Aldo Di Carlo                        !
!!                                                                          !
!! libNEGF is free software: you can redistribute it and/or modify          !
!! it under the terms of the GNU Lesse General Public License as published  !
!! by the Free Software Foundation, either version 3 of the License, or     !
!! (at your option) any later version.                                      !
!!                                                                          !
!!  You should have received a copy of the GNU Lesser General Public        !
!!  License along with libNEGF.  If not, see                                !
!!  <http://www.gnu.org/licenses/>.                                         !  
!!--------------------------------------------------------------------------!


module mpi_globals

#:if defined("MPI")
  
  use libmpifx_module
  
  INTEGER, SAVE ::  mpi_comm
  INTEGER, SAVE ::  numprocs
  INTEGER, SAVE ::  id
  LOGICAL, SAVE ::  id0

  contains
    
    subroutine negf_mpi_init(energyComm, ioMaster)
      type(mpifx_comm) :: energyComm
      logical, optional :: ioMaster
      
      id =energyComm%rank
      numprocs = energyComm%size
      
      if (present(ioMaster)) then
        id0 = ioMaster 
      else
        id0 = (id == 0)
      end if    
      !print*, 'INIT MPI-NEGF ON',numprocs,'NODES' 
      !print*, 'CPU',id,'READY'
      !print*, 'PRINTING CPU:',id0
      
    end subroutine negf_mpi_init

#:else
  
  logical, parameter ::  id0 = .true.
  integer, parameter :: id = 0, numprocs = 1
         
#:endif


end module mpi_globals
