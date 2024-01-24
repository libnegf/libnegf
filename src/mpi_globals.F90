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
  use mpi
  use libmpifx_module, only : mpifx_comm
  private
#:endif

  integer, public ::  numprocs = 1
  integer, public ::  id = 0
  logical, public ::  id0 = .true.

#:if defined("MPI")
  public :: negf_mpi_init
  public :: negf_cart_init
  public :: check_cart_comm

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
    !
    ! Order: dim 1: k; dim 2: E
    ! It must be periodic in k for all-to-all communications
    ! CAVEAT:
    ! All processes MUST have the same number of points in K and E
    ! For E it is used to compute where E +/- wq are located
    ! For K it is used to compute where another q is placed
    !
    subroutine negf_cart_init(inComm, nk, cartComm, energyComm, kComm, bareCartComm, barekComm)
      !> Input communicator
      type(mpifx_comm), intent(in) :: inComm
      !> Number of processors for k
      integer, intent(in) :: nk
      !> Output 2D cartesian communicator
      type(mpifx_comm) :: cartComm
      !> Output communicator for the energy sub-grid
      type(mpifx_comm), intent(out) :: energyComm
      !> Output communicator for the k sub-grid
      type(mpifx_comm), intent(out) :: kComm

      !> Output communicators of type int for TiberCAD
      integer, intent(out), optional :: bareCartComm, barekComm

      integer :: outComm
      integer :: ndims = 2
      integer :: dims(2)
      logical :: periods(2) = .false.
      logical :: remain_dims(2)
      integer :: nE
      logical :: reorder = .true.
      integer :: mpierr

      if (mod(inComm%size,nk) /=0 ) then
        stop "Error in cart_init: cannot build a 2D cartesian grid with incompatible sizes"
      end if

      nE = inComm%size/nk
      dims(1)=nk; dims(2)=nE
      periods(1) = .true.

      call MPI_CART_CREATE(inComm%id, ndims, dims, periods, reorder, outComm, mpierr)
      call cartComm%init(outComm, mpierr)
      ! Global master id=0 node as writing node
      id0 = (cartComm%rank == 0)
      if (present(bareCartComm)) bareCartComm = outComm

      ! Extract sub-communicators
      remain_dims(:) = [.false., .true.]
      call MPI_CART_SUB(cartComm%id, remain_dims, outComm, mpierr)
      call energyComm%init(outComm, mpierr)
      id = energyComm%rank
      numprocs = energyComm%size

      remain_dims(:) = [.true., .false.]
      call MPI_CART_SUB(cartComm%id, remain_dims, outComm, mpierr)
      call kComm%init(outComm, mpierr)
      if (present(barekComm)) barekComm = outComm

    end subroutine negf_cart_init

    subroutine check_cart_comm(cartComm, mpierror)
      !> Input 2d cartesian communicator
      type(mpifx_comm), intent(in) :: cartComm
      !> output error
      integer, intent(out) :: mpierror

      integer :: coords(2)

      mpierror = 0

      call MPI_Cart_coords(cartComm%id, 0, 2, coords, mpierror)

    end subroutine check_cart_comm

#:endif

end module mpi_globals
