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

  use mpi, only : mpi_cart_create, mpi_cart_sub
  use libmpifx_module, only : mpifx_comm

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
      type(mpifx_comm), intent(out) :: energyComm
      type(mpifx_comm), intent(out) :: kComm

      type(mpifx_comm) :: cartComm
      integer :: outComm
      integer :: ndims
      integer :: dims(2)
      logical :: periods(2)
      logical :: remain_dims(2)
      integer :: nE
      logical :: reorder
      integer :: mpierr

      ndims = 2
      periods(:) = .false.
      reorder = .true.

      if (mod(inComm%size,nk) /=0 ) then
        stop "Error in cart_init: cannot build a 2D cartesian grid with incompatible sizes"
      end if
      nE = inComm%size/nk
      dims(1)=nk; dims(2)=nE

      call MPI_CART_CREATE(inComm%id, ndims, dims, periods, reorder, outComm, mpierr)
      call cartComm%init(outComm, mpierr)

      remain_dims(:) = [.false., .true.]
      call MPI_CART_SUB(cartComm%id, remain_dims, outComm, mpierr)
      call energyComm%init(outComm, mpierr)

      remain_dims(:) = [.true., .false.]
      call MPI_CART_SUB(cartComm%id, remain_dims, outComm, mpierr)
      call kComm%init(outComm, mpierr)

    end subroutine negf_cart_init


#:endif

end module mpi_globals
