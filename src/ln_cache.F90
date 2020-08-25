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

module ln_cache

  use lib_param
  use mat_def, only: z_DNS

  implicit none

  public :: add_surface_green_to_cache, retrieve_surface_green_from_cache

  private

  !> Linked list to store surface green in memory
  type :: TCachedSurfaceGreen
    type(TCachedSurfaceGreen), pointer :: next => null()
    type(z_DNS) :: surface_green
    integer :: kpoint = -1
    integer :: energy_point = -1
    integer :: spin = -1
    integer :: contact = -1
  end type

  type(TCachedSurfaceGreen), target, allocatable, save :: cache_start

contains

  subroutine add_surface_green_to_cache(surface_green, contact, nkp, pnt, nsp)
    type(z_DNS) :: surface_green
    integer :: nkp
    integer :: pnt
    integer :: nsp
    integer :: contact

    type(TCachedSurfaceGreen), pointer :: p

    write(*,*) "Add ", nkp, pnt, nsp
    if (.not.allocated(cache_start)) then
        write(*,*) "Here"
      allocate(cache_start)
      p => cache_start
    else
      p => cache_start
      do
        !while (associated(p%next))
        ! If the point is already present, substitute and exit
        if (p%kpoint .eq. nkp .and. &
        & p%energy_point .eq. pnt .and. &
        & p%spin .eq. nsp .and. p%contact .eq. contact) then

          p%surface_green = surface_green
          write (*, *) 'Substitute ', nkp, pnt, nsp
          return

        else if (.not.associated(p%next)) then
            write (*, *) 'From scratch ', nkp, pnt, nsp
            allocate (p%next)
            p => p%next
            exit
        end if
        p => p%next
      end do

    end if
    p%surface_green = surface_green
    p%kpoint = nkp
    p%energy_point = pnt
    p%spin = nsp
    p%contact = contact

  end subroutine

  subroutine retrieve_surface_green_from_cache(surface_green, contact, nkp, pnt, nsp)
    type(z_DNS) :: surface_green
    integer :: nkp
    integer :: pnt
    integer :: nsp
    integer :: contact

    type(TCachedSurfaceGreen), pointer :: p

    write(*,*) "Retrieve ", nkp, pnt, nsp
    if (.not. allocated(cache_start)) then
      error stop "No entry in surface green cache"
    else
      p => cache_start
    end if

    do while (associated(p))
        write(*,*) "Points B ", p%kpoint, p%energy_point, p%spin
        if (p%kpoint .eq. nkp .and. &
        & p%energy_point .eq. pnt .and. &
        & p%spin .eq. nsp .and. p%contact .eq. contact) then
            write(*,*) "found during loop B"
            surface_green = p%surface_green
            return
        end if
        p => p%next
    end do
    write(*,*) "Points C ", p%kpoint, p%energy_point, p%spin

    error stop "Entry not found in surface green cache"

      end subroutine

      end module
