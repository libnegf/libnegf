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

  use mat_def, only: z_DNS

  implicit none

  public :: TSurfaceGreenCache

  private

  !> Linked list to store surface green in memory
  type :: TSurfaceGreenCacheEntry
    type(TSurfaceGreenCacheEntry), pointer :: next => null()
    type(z_DNS), allocatable :: surface_green
    integer :: kpoint = -1
    integer :: energy_point = -1
    integer :: spin = -1
    integer :: contact = -1
  end type

  type :: TSurfaceGreenCache
    type(TSurfaceGreenCacheEntry), pointer :: first => null()
  contains
    procedure :: add => add_surface_green_to_cache
    procedure :: retrieve => retrieve_surface_green_from_cache
    procedure :: destroy => destroy_surface_green_cache
  end type

contains

  subroutine add_surface_green_to_cache(this, surface_green, contact, nkp, pnt, nsp)
    class(TSurfaceGreenCache) :: this
    type(z_DNS):: surface_green
    integer :: nkp
    integer :: pnt
    integer :: nsp
    integer :: contact

    type(TSurfaceGreenCacheEntry), pointer :: p

    if (.not. associated(this%first)) then
      allocate (this%first)
      p => this%first
    else
      p => this%first
      do
        !while (associated(p%next))
        ! If the point is already present, substitute and exit
        if (p%kpoint .eq. nkp .and. &
        & p%energy_point .eq. pnt .and. &
        & p%spin .eq. nsp .and. p%contact .eq. contact) then

          p%surface_green = surface_green
          return

        else if (.not. associated(p%next)) then
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

  subroutine retrieve_surface_green_from_cache(this, surface_green, contact, nkp, pnt, nsp)
    class(TSurfaceGreenCache) :: this
    type(z_DNS) :: surface_green
    integer :: nkp
    integer :: pnt
    integer :: nsp
    integer :: contact

    type(TSurfaceGreenCacheEntry), pointer :: p

    if (.not. associated(this%first)) then
      error stop "No entry in surface green cache"
    else
      p => this%first
    end if

    do while (associated(p))
      if (p%kpoint .eq. nkp .and. &
      & p%energy_point .eq. pnt .and. &
      & p%spin .eq. nsp .and. p%contact .eq. contact) then
        surface_green = p%surface_green
        return
      end if
      p => p%next
    end do

    error stop "Entry not found in surface green cache"

  end subroutine

  subroutine destroy_surface_green_cache(this)
    class(TSurfaceGreenCache) :: this

    type(TSurfaceGreenCacheEntry), pointer :: p, previous

    if (.not.associated(this%first)) then
      return
    end if

    p => this%first
    do while (associated(p))
      previous => p
      p => p%next
      deallocate (previous)
    end do

  end subroutine

end module
