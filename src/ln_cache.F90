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
  use ln_precision
  use globals
  use outmatrix, only: outmat_c, inmat_c

  implicit none

  public :: TSurfaceGreenCache, TSurfaceGreenCacheMem, TSurfaceGreenCacheDisk
  public :: TSurfaceGreenCacheDummy

  private

  type, abstract :: TSurfaceGreenCache
  contains
    procedure(abst_add_surface_green_to_cache), deferred :: add
    procedure(abst_retrieve_surface_green_from_cache), deferred :: retrieve
    procedure(abst_destroy_surface_green_cache), deferred :: destroy
    procedure(abst_is_cached), deferred :: is_cached
  end type TSurfaceGreenCache

  abstract interface
    subroutine abst_add_surface_green_to_cache(this, surface_green, contact, nkp, pnt, nsp)
      import :: TSurfaceGreenCache
      import :: z_DNS
      class(TSurfaceGreenCache) :: this
      type(z_DNS):: surface_green
      integer :: nkp
      integer :: pnt
      integer :: nsp
      integer :: contact
    end subroutine

    subroutine abst_retrieve_surface_green_from_cache(this, surface_green, contact, nkp, pnt, nsp)
      import :: TSurfaceGreenCache
      import :: z_DNS
      class(TSurfaceGreenCache) :: this
      type(z_DNS) :: surface_green
      integer :: nkp
      integer :: pnt
      integer :: nsp
      integer :: contact

    end subroutine

    function abst_is_cached(this, contact, nkp, pnt, nsp) result(val)
      import :: TSurfaceGreenCache
      !> The bound class
      class(TSurfaceGreenCache) :: this
      !> The k point index
      integer :: nkp
      !> The energy point index
      integer :: pnt
      !> The spin index
      integer :: nsp
      !> The contact index
      integer :: contact
      !> Whether the corresponding surface green function is cached
      logical :: val
    end function

    subroutine abst_destroy_surface_green_cache(this)
      import :: TSurfaceGreenCache
      class(TSurfaceGreenCache) :: this
    end subroutine
  end interface

  !> Linked list to store surface green in memory
  type :: TSurfaceGreenCacheEntry
    type(TSurfaceGreenCacheEntry), pointer :: next => null()
    type(z_DNS), allocatable :: surface_green
    integer :: kpoint = -1
    integer :: energy_point = -1
    integer :: spin = -1
    integer :: contact = -1
  end type

  type, extends(TSurfaceGreenCache) :: TSurfaceGreenCacheMem
    type(TSurfaceGreenCacheEntry), pointer :: first => null()
  contains
    procedure :: add => mem_add
    procedure :: retrieve => mem_retrieve
    procedure :: destroy => mem_destroy
    procedure :: is_cached => mem_is_cached
  end type

  type, extends(TSurfaceGreenCache) :: TSurfaceGreenCacheDisk
    character(len=LST) :: scratch_path
  contains
    procedure :: add => disk_add
    procedure :: retrieve => disk_retrieve
    procedure :: destroy => disk_destroy
    procedure :: is_cached => disk_is_cached
  end type

  type, extends(TSurfaceGreenCache) :: TSurfaceGreenCacheDummy
  contains
    procedure :: add => dummy_add
    procedure :: retrieve => dummy_retrieve
    procedure :: destroy => dummy_destroy
    procedure :: is_cached => dummy_is_cached
  end type

contains

  !! Definitions for TSurfaceGreenCacheMem

  subroutine mem_add(this, surface_green, contact, nkp, pnt, nsp)
    class(TSurfaceGreenCacheMem) :: this
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

  subroutine mem_retrieve(this, surface_green, contact, nkp, pnt, nsp)
    class(TSurfaceGreenCacheMem) :: this
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

    error stop "Cannot retrieve surface green function"

  end subroutine

  function mem_is_cached(this, contact, nkp, pnt, nsp) result(val)
    class(TSurfaceGreenCacheMem) :: this
    integer :: nkp
    integer :: pnt
    integer :: nsp
    integer :: contact
    logical :: val

    type(TSurfaceGreenCacheEntry), pointer :: p

    if (.not. associated(this%first)) then
      val = .false.
      return
    else
      p => this%first
    end if

    do while (associated(p))
      if (p%kpoint .eq. nkp .and. &
      & p%energy_point .eq. pnt .and. &
      & p%spin .eq. nsp .and. p%contact .eq. contact) then
        val = .true.
        return
      end if
      p => p%next
    end do

    val = .false.

  end function

  subroutine mem_destroy(this)
    class(TSurfaceGreenCacheMem) :: this

    type(TSurfaceGreenCacheEntry), pointer :: p, previous

    if (.not. associated(this%first)) then
      return
    end if

    p => this%first
    do while (associated(p))
      previous => p
      p => p%next
      deallocate (previous)
    end do

  end subroutine

  !! Definitions for TSurfaceGreenCacheDisk

  subroutine disk_add(this, surface_green, contact, nkp, pnt, nsp)
    class(TSurfaceGreenCacheDisk) :: this
    type(z_DNS):: surface_green
    integer :: nkp
    integer :: pnt
    integer :: nsp
    integer :: contact

    real(kind=dp) :: dens
    character(2) :: ofcont
    character(1) :: ofspin
    character(10) :: ofkpnt
    character(10) :: ofpnt
    character(LST) :: filename
    integer :: file_unit

    call disk_indices_to_filename(filename, contact, nkp, pnt, nsp)
    open (newunit=file_unit, file=trim(this%scratch_path)//filename, form='UNFORMATTED')
    call outmat_c(file_unit, .false., surface_green%val, surface_green%nrow, surface_green%ncol)
    close (file_unit)

  end subroutine

  subroutine disk_retrieve(this, surface_green, contact, nkp, pnt, nsp)
    class(TSurfaceGreenCacheDisk) :: this
    !> Input surface green to be retrieved: It needs to be already allocated with
    !> the correct number of rows and columns.
    type(z_DNS) :: surface_green
    integer :: nkp
    integer :: pnt
    integer :: nsp
    integer :: contact

    real(kind=dp) :: dens
    character(2) :: ofcont
    character(1) :: ofspin
    character(10) :: ofkpnt
    character(10) :: ofpnt
    character(LST) :: filename
    logical :: file_exists
    integer :: file_unit

    call disk_indices_to_filename(filename, contact, nkp, pnt, nsp)
    inquire (file=trim(this%scratch_path)//filename, EXIST=file_exists)

    if (.not. file_exists) then
      error stop "Cannot retrieve surface green function from disk: file not found"
    end if

    open (newunit=file_unit, file=trim(this%scratch_path)//filename, form='UNFORMATTED')
    call inmat_c(file_unit, .false., surface_green%val, surface_green%nrow, surface_green%ncol)
    close (file_unit)

  end subroutine

  function disk_is_cached(this, contact, nkp, pnt, nsp) result(val)
    class(TSurfaceGreenCacheDisk) :: this
    integer :: nkp
    integer :: pnt
    integer :: nsp
    integer :: contact
    logical :: val

    character(LST) :: filename
    logical :: file_exists
    integer :: file_unit

    call disk_indices_to_filename(filename, contact, nkp, pnt, nsp)
    inquire (file=trim(this%scratch_path)//filename, EXIST=file_exists)

    if (file_exists) then
      val = .true.
    else
      val = .false.
    end if

  end function

  subroutine disk_destroy(this)
    class(TSurfaceGreenCacheDisk) :: this
    ! Empty. Defined only for interface.
  end subroutine

  subroutine disk_indices_to_filename(filename, contact, nkp, pnt, nsp)
    character(LST), intent(out) :: filename
    integer :: nkp
    integer :: pnt
    integer :: nsp
    integer :: contact

    character(2) :: ofcont
    character(1) :: ofspin
    character(10) :: ofkpnt
    character(10) :: ofpnt

    write (ofcont, '(i2.2)') contact
    if (nkp .le. 99) write (ofkpnt, '(i2.2)') nkp
    if (nkp .gt. 99) write (ofkpnt, '(i3.3)') nkp
    if (nkp .gt. 999) write (ofkpnt, '(i4.4)') nkp
    if (nkp .gt. 9999) stop 'ERROR: too many k-points (> 9999)'
    if (pnt .le. 999) write (ofpnt, '(i3.3)') pnt
    if (pnt .gt. 999) write (ofpnt, '(i4.4)') pnt
    if (pnt .gt. 9999) write (ofpnt, '(i5.5)') pnt
    if (pnt .gt. 99999) stop 'ERROR: too many contour points (> 99999)'
    if (nsp .eq. 1) ofspin = 'u'
    if (nsp .eq. 2) ofspin = 'd'
    filename = 'GS'//ofspin//ofcont//'_'//trim(ofkpnt)//'_'//trim(ofpnt)//'.dat'

  end subroutine

  !! Definitions for TSurfaceGreenCacheDummy

  subroutine dummy_add(this, surface_green, contact, nkp, pnt, nsp)
    class(TSurfaceGreenCacheDummy) :: this
    type(z_DNS):: surface_green
    integer :: nkp
    integer :: pnt
    integer :: nsp
    integer :: contact

    ! Dummy operation

  end subroutine

  subroutine dummy_retrieve(this, surface_green, contact, nkp, pnt, nsp)
    class(TSurfaceGreenCacheDummy) :: this
    !> Input surface green to be retrieved: It needs to be already allocated with
    !> the correct number of rows and columns.
    type(z_DNS) :: surface_green
    integer :: nkp
    integer :: pnt
    integer :: nsp
    integer :: contact

    ! Dummy operation

  end subroutine

  function dummy_is_cached(this, contact, nkp, pnt, nsp) result(val)
    class(TSurfaceGreenCacheDummy) :: this
    integer :: nkp
    integer :: pnt
    integer :: nsp
    integer :: contact
    logical :: val

    val = .false.

  end function

  subroutine dummy_destroy(this)
    class(TSurfaceGreenCacheDummy) :: this
  end subroutine
end module
