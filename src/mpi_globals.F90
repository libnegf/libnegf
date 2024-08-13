!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Greens functions.         !
!! Copyright (C) 2012 - 2026                                                !
!!                                                                          !
!! This file is part of libNEGF: a library for                              !
!! Non Equilibrium Green's Functions calculations                           !
!!                                                                          !
!! Developers: Alessandro Pecchia, Daniele Soccodato                        !
!! Former Contributors: Gabriele Penazzi, Luca Latessa, Aldo Di Carlo       !
!!                                                                          !
!! libNEGF is free software: you can redistribute and/or modify it          !
!! under the terms of the GNU Lesser General Public License as published    !
!! by the Free Software Foundation, either version 3 of the License, or     !
!! (at your option) any later version.                                      !
!!                                                                          !
!!  You should have received a copy of the GNU Lesser General Public        !
!!  License along with libNEGF.  If not, see                                !
!!  <http://www.gnu.org/licenses/>.                                         !
!!--------------------------------------------------------------------------!


module mpi_globals

#:if defined("MPI")
  use mpi_f08
  use omp_lib, only : omp_get_max_threads, omp_get_num_threads
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
  public :: MPI_Comm

  contains

    ! Purpose: 
    ! Set the global variables 'numprocs' and proc 'id' 
    ! within the energy communicator    
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

    ! Initialize a 2D cartesian grid.
    ! INPUT: global communicator and size of k communicator
    ! OUPUT: the cartesian communicator, the energy- and k- communicators
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
      type(MPI_Comm), intent(out), optional :: bareCartComm, barekComm

      type(MPI_Comm) :: outComm
      integer :: ndims = 2
      integer :: dims(2)
      logical :: periods(2) = .false.
      logical :: remain_dims(2)
      integer :: nE
      logical :: reorder = .true.
      integer :: mpierr

      if (mod(inComm%size,nk) /=0 ) then
        error stop "Error in cart_init: cannot build a 2D cartesian grid with incompatible sizes"
      end if

      !call check_omp_mpi(inComm, mpierr)

      nE = inComm%size/nk
      dims(1)=nk; dims(2)=nE
      periods(1) = .true.

      call MPI_CART_CREATE(inComm%comm, ndims, dims, periods, reorder, outComm, mpierr)
      call cartComm%init(outComm, mpierr)
      ! Global master id=0 node as writing node
      id0 = (cartComm%rank == 0)
      if (present(bareCartComm)) bareCartComm = outComm

      ! Extract sub-communicators
      remain_dims(:) = [.false., .true.]
      call MPI_CART_SUB(cartComm%comm, remain_dims, outComm, mpierr)
      call energyComm%init(outComm, mpierr)
      id = energyComm%rank
      numprocs = energyComm%size

      remain_dims(:) = [.true., .false.]
      call MPI_CART_SUB(cartComm%comm, remain_dims, outComm, mpierr)
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

      call MPI_Cart_coords(cartComm%comm, 0, 2, coords, mpierror)

    end subroutine check_cart_comm

   ! subroutine check_omp_mpi(comm, mpierror)
   !   type(mpifx_comm), intent(in) :: comm
   !   integer, intent(out) :: mpierror

   !   integer :: newcomm, info, nprocs_shared, phys_cores, num_threads

   !   call MPI_COMM_SPLIT_TYPE(comm%id, MPI_COMM_TYPE_SHARED, 1, info, newcomm, mpierror)
   !   if (mpierror /= 0) then
   !      stop "ERROR in MPI_COMM_SPLIT_TYPE"
   !   end if
   !   call MPI_COMM_SIZE(newcomm, nprocs_shared, mpierror)

   !   phys_cores = get_num_cores()
   !   num_threads = omp_get_max_threads()

   !   if (comm%rank==0) then
   !     print*, "Number of physical cores:", phys_cores
   !     print*, "Number of processors on same shared memory:",nprocs_shared
   !     print*, "Maximum number of OMP threads:",num_threads
   !   end if
   !   if (num_threads*nprocs_shared > phys_cores) then
   !     call omp_set_max_threads(phys_cores/nprocs_shared)
   !   end if

   ! end subroutine check_omp_mpi

#:endif

end module mpi_globals
