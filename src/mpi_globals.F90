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

#:endif

  INTEGER, SAVE ::  mpi_comm
  INTEGER, SAVE ::  numprocs = 1
  INTEGER, SAVE ::  id = 0
  LOGICAL, SAVE ::  id0 = .true.

#:if defined("MPI")

  contains

    subroutine negf_mpi_init(energyComm, ioProc)
      type(mpifx_comm) :: energyComm
      logical, optional :: ioProc

      id =energyComm%rank
      numprocs = energyComm%size

      if (present(ioProc)) then
        id0 = ioProc
      else
        id0 = (id == 0)
      end if

    end subroutine negf_mpi_init

    ! Initialize a 2D cartesian grid
    subroutine negf_cart_init(inComm, nk, cartComm, energyComm, kComm)
      type(mpifx_comm), intent(in) :: inComm
      integer, intent(in) :: nk
      type(mpifx_comm), intent(out) :: cartComm
      type(mpifx_comm), intent(out) :: energyComm
      type(mpifx_comm), intent(out) :: kComm

      integer :: outComm
      integer :: ndims = 2
      integer :: dims(2)
      integer :: period(2) = 0
      integer :: remain_dims(2)
      integer :: nE
      integer :: reorder = 1
      integer :: mpierr

      if (mod(inComm%size,nk) /=0 ) then
        stop "Error in cart_init: cannot build a 2D cartesian grid with incompatible sizes"
      end if
      nE = inComm%size/nk
      dims(1)=nk; dims(2)=nE

      call MPI_CART_CREATE(inComm%id, ndims, dims, periods, reorder, outComm, mpierr)
      call cartComm%init(outComm, mpierr)

      remain_dims = (/0, 1/)
      call MPI_CART_SUB(cartComm%id, remain_dims, outComm, mpierr)
      call energyComm%init(outComm, mpierr)

      remain_dims = (/1, 0/)
      call MPI_CART_SUB(cartComm%id, remain_dims, outComm, mpierr)
      call kComm%init(outComm, mpierr)

    end subroutine negf_cart_init


#:endif

end module mpi_globals
