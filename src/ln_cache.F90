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

module ln_cache
  use iso_c_binding
  use ln_precision, only: dp
  use mat_def, only: z_DNS, c_DNS, x_DNS => z_DNS, destroy, assignment(=)
  use globals
  use outmatrix, only: outmat_c, inmat_c
  use ln_messages, only: error_msg

  implicit none
  private

  public :: TMatrixCache
  public :: TMatrixCacheMem, TMatrixCacheDisk
  public :: TMatrixCacheDummy
  public :: TMatLabel
  public :: print_label
  public :: get_string_label

  !> Default point of labels and cache entries: callers that don't set it never match a set one
  real(dp), parameter :: UNSET_POINT = huge(1.0_dp)

  !> Tolerance to decide whether a cached matrix belongs to the requested point
  real(dp), parameter :: EPS_POINT = 1.0e-10_dp

  type TMatLabel
    integer :: kpoint = -1
    integer :: energy_point = -1
    integer :: spin = -1
    integer :: row_block = -1
    integer :: col_block = -1
    !> Energy and k-vector of the matrix, checked on lookup.
    !> Callers that don't set them keep the defaults, which match each other.
    !> Only the real part of energy is UNSET_POINT: with both parts set, abs() of the
    !> difference to a set energy would overflow, while (huge - huge) is exactly zero.
    complex(dp) :: energy = (UNSET_POINT, 0.0_dp)
    real(dp) :: kvec(3) = UNSET_POINT
  end type TMatLabel

  type, abstract :: TMatrixCache
  contains
    procedure(abst_add_matrix_dp), deferred :: add_dp
    procedure(abst_add_matrix_sp), deferred :: add_sp
    generic :: add => add_sp, add_dp
    procedure(abst_retrieve_matrix_dp), deferred :: retrieve_dp
    procedure(abst_retrieve_matrix_sp), deferred :: retrieve_sp
    generic :: retrieve => retrieve_sp, retrieve_dp
    procedure(abst_retrieve_matrix_pointer), deferred :: retrieve_pointer
    procedure(abst_retrieve_matrix_loc), deferred :: retrieve_loc
    procedure(abst_destroy_matrix), deferred :: destroy
    procedure(abst_clear_matrix), deferred :: clear
    procedure(abst_is_cached), deferred :: is_cached
    procedure(abst_list_cache), deferred :: list_cache
  end type TMatrixCache

  abstract interface
    subroutine abst_add_matrix_dp(this, matrix, label)
      import :: TMatrixCache
      import :: z_DNS
      import :: TMatLabel
      class(TMatrixCache) :: this
      type(z_DNS):: matrix
      type(TMatLabel) :: label
    end subroutine

    subroutine abst_add_matrix_sp(this, matrix, label)
      import :: TMatrixCache
      import :: c_DNS
      import :: TMatLabel
      class(TMatrixCache) :: this
      type(c_DNS):: matrix
      type(TMatLabel) :: label
    end subroutine

    subroutine abst_retrieve_matrix_dp(this, matrix, label)
      import :: TMatrixCache
      import :: z_DNS
      import :: TMatLabel
      class(TMatrixCache) :: this
      type(z_DNS) :: matrix
      type(TMatLabel) :: label
    end subroutine

    subroutine abst_retrieve_matrix_sp(this, matrix, label)
      import :: TMatrixCache
      import :: c_DNS
      import :: TMatLabel
      class(TMatrixCache) :: this
      type(c_DNS) :: matrix
      type(TMatLabel) :: label
    end subroutine

    subroutine abst_retrieve_matrix_pointer(this, pmatrix, label)
      import :: x_DNS
      import :: TMatrixCache
      import :: TMatLabel
      class(TMatrixCache) :: this
      type(x_DNS), pointer :: pmatrix
      type(TMatLabel) :: label
    end subroutine

    function abst_retrieve_matrix_loc(this, label) result(pmatrix)
      use iso_c_binding
      import :: TMatrixCache
      import :: TMatLabel
      class(TMatrixCache) :: this
      type(TMatLabel) :: label
      type(C_PTR) :: pmatrix
    end function

    function abst_is_cached(this, label) result(val)
      import :: TMatrixCache
      import :: TMatLabel
      !> The bound class
      class(TMatrixCache) :: this
      !> matrix label identifier
      type(TMatLabel) :: label
      logical :: val
    end function

    subroutine abst_list_cache(this)
      import :: TMatrixCache
      import :: TMatLabel
      !> The bound class
      class(TMatrixCache) :: this
    end subroutine

    subroutine abst_destroy_matrix(this)
      import :: TMatrixCache
      class(TMatrixCache) :: this
    end subroutine

    subroutine abst_clear_matrix(this)
      import :: TMatrixCache
      class(TMatrixCache) :: this
    end subroutine
  end interface

  !> Linked list to store matrix in memory
  !type :: TMatrixCacheEntry
  !  type(TMatrixCacheEntry), pointer :: next => null()
  !  type(z_DNS),  pointer :: matrix => null()
  !  type(TMatLabel) :: mat_label
  !end type

  !> element to store matrix in memory - single or double
  type :: TMatrixCacheEntry
    type(x_DNS),  pointer :: matrix => null()
    !> Point the matrix was computed at (see TMatLabel)
    complex(dp) :: energy = (UNSET_POINT, 0.0_dp)
    real(dp) :: kvec(3) = UNSET_POINT
  end type

  !> Mem derived class
  !> Internal storage within an array. Blocks stored per diagonals.
  !> Diagonal, Superdiagonal, Subdiagonal.
  type, extends(TMatrixCache) :: TMatrixCacheMem
    type(TMatrixCacheEntry), allocatable :: MatArray(:,:,:,:)
    integer :: Nblocks = 0
    logical :: isInitialized = .false.
    character(len=LST) :: tagname
  contains
    procedure :: init
    procedure :: checks
    procedure :: in_range
    procedure :: add_sp => mem_add_sp
    procedure :: add_dp => mem_add_dp
    procedure :: retrieve_sp => mem_retrieve_sp
    procedure :: retrieve_dp => mem_retrieve_dp
    procedure :: retrieve_pointer => mem_retrieve_pointer
    procedure :: retrieve_loc => mem_retrieve_loc
    procedure :: destroy => mem_destroy
    procedure :: clear => mem_clear
    procedure :: is_cached => mem_is_cached
    procedure :: list_cache => mem_list_cache
  end type

  !> Disk derived class
  type, extends(TMatrixCache) :: TMatrixCacheDisk
    character(len=LST) :: scratch_path
  contains
    procedure :: add_sp => disk_add_sp
    procedure :: add_dp => disk_add_dp
    procedure :: retrieve_dp => disk_retrieve_dp
    procedure :: retrieve_sp => disk_retrieve_sp
    procedure :: retrieve_pointer => disk_retrieve_pointer
    procedure :: retrieve_loc => disk_retrieve_loc
    procedure :: destroy => disk_destroy
    procedure :: clear => disk_clear
    procedure :: is_cached => disk_is_cached
    procedure :: list_cache => disk_list_cache
  end type

  !> Dummy cache when no caching is performed
  type, extends(TMatrixCache) :: TMatrixCacheDummy
  contains
    procedure :: add_sp => dummy_add_sp
    procedure :: add_dp => dummy_add_dp
    procedure :: retrieve_dp => dummy_retrieve_dp
    procedure :: retrieve_sp => dummy_retrieve_sp
    procedure :: retrieve_pointer => dummy_retrieve_pointer
    procedure :: retrieve_loc => dummy_retrieve_loc
    procedure :: destroy => dummy_destroy
    procedure :: clear => dummy_clear
    procedure :: is_cached => dummy_is_cached
    procedure :: list_cache => dummy_list_cache
  end type

contains

  subroutine print_label(label)
    type(TMatLabel) :: label

    print*,'k=',label%kpoint,'E=',label%energy_point,'s=',label%spin, &
          & 'i=',label%row_block,'j=',label%col_block

  end subroutine print_label

  !! Definitions for TMatrixCacheMem
  subroutine init(this, NE, Nk, NPL, Ndiag, Nspin)
    class(TMatrixCacheMem) :: this
    integer, intent(in) :: Nk, NE, NPL, Ndiag, Nspin

    integer :: ierr
    character(len=256) :: errmsg

    if (allocated(this%MatArray)) then
      if (size(this%MatArray) /= NE*Nk*NPL*Ndiag*Nspin) then
        call this%destroy()
      end if
    end if
    if (.not.allocated(this%MatArray)) then
      allocate(this%MatArray(NE,Nk,NPL*Ndiag,Nspin), stat=ierr, errmsg=errmsg)
      if (ierr /= 0) then
        call error_msg( "Allocation error of MatArray in init "//trim(this%tagname)//": "&
            & //trim(errmsg))
      end if
    end if
    this%Nblocks = NPL
    this%isInitialized = .true.
  end subroutine init

  subroutine checks(this, label, bl)
    class(TMatrixCacheMem) :: this
    type(TMatLabel), intent(in) :: label
    integer, intent(out) :: bl

    if (.not.allocated(this%MatArray)) then
      call error_msg( "MatArray not allocated for "//trim(this%tagname))
    end if
    if (label%energy_point < 1 .or. label%energy_point > size(this%MatArray,1)) then
       print*,trim(this%tagname)
       call print_label(label)
       call error_msg( "energy_point out of range")
    end if
    if (label%kpoint < 1 .or. label%kpoint > size(this%MatArray,2)) then
       print*,trim(this%tagname)
       call print_label(label)
       call error_msg( "k_point out of range")
    end if
    if (label%spin < 1 .or. label%spin > size(this%MatArray,4)) then
       print*,trim(this%tagname)
       call print_label(label)
       call error_msg( "spin out of range")
    end if

    bl = block_index(this, label)
    if (bl < 1 .or. bl > size(this%MatArray,3)) then
       print*,trim(this%tagname)
       call print_label(label)
       print*, "bl=",bl
       call error_msg( "block out of range")
    end if
  end subroutine checks

  !> Like checks, but returns .false. instead of stopping when the slot is out of range
  function in_range(this, label, bl) result(val)
    class(TMatrixCacheMem) :: this
    type(TMatLabel), intent(in) :: label
    integer, intent(out) :: bl
    logical :: val

    bl = block_index(this, label)
    val = allocated(this%MatArray)
    if (.not.val) return
    val = label%energy_point >= 1 .and. label%energy_point <= size(this%MatArray,1) &
        & .and. label%kpoint >= 1 .and. label%kpoint <= size(this%MatArray,2) &
        & .and. bl >= 1 .and. bl <= size(this%MatArray,3) &
        & .and. label%spin >= 1 .and. label%spin <= size(this%MatArray,4)
  end function in_range

  !> Index of block (row_block, col_block) in the third dimension of MatArray
  pure function block_index(this, label) result(bl)
    class(TMatrixCacheMem), intent(in) :: this
    type(TMatLabel), intent(in) :: label
    integer :: bl

    integer :: r, c

    r = label%row_block
    c = label%col_block
    if (c >= r) then
       bl = (c-r)*this%Nblocks + c
    else
       bl = 2*(r-c)*this%Nblocks + r
    end if
  end function block_index

  !> Whether a matrix computed at (energy, kvec) belongs to the point of label
  pure function same_point(energy, kvec, label) result(val)
    complex(dp), intent(in) :: energy
    real(dp), intent(in) :: kvec(3)
    type(TMatLabel), intent(in) :: label
    logical :: val

    val = abs(energy - label%energy) < EPS_POINT .and. all(abs(kvec - label%kvec) < EPS_POINT)
  end function same_point

  subroutine mem_add_dp(this, matrix, label)
    class(TMatrixCacheMem) :: this
    type(z_DNS):: matrix
    type(TMatLabel) :: label

    integer :: b

    ! blocks r,r first. r,r+1 after. r+1,r
    ! 11, 22, 33, 44, 12, 23, 34, 21, 32, 43
    !  1,  2,  3,  4,  6,  7,  8, 10, 11, 12
    call checks(this, label, b)
    associate(p => this%MatArray(label%energy_point, label%kpoint, b, label%spin))
      if (associated(p%matrix)) then
         call destroy(p%matrix)
      else
         allocate(p%matrix)
      end if
      p%matrix = matrix
      p%energy = label%energy
      p%kvec = label%kvec
    end associate

  end subroutine mem_add_dp

  subroutine mem_add_sp(this, matrix, label)
    class(TMatrixCacheMem) :: this
    type(c_DNS):: matrix
    type(TMatLabel) :: label

    integer :: b

    ! blocks r,r first. r,r+1 after. r+1,r
    ! 11, 22, 33, 44, 12, 23, 34, 21, 32, 43
    !  1,  2,  3,  4,  6,  7,  8, 10, 11, 12
    call checks(this, label, b)
    associate(p => this%MatArray(label%energy_point, label%kpoint, b, label%spin))
      if (associated(p%matrix)) then
         call destroy(p%matrix)
      else
         allocate(p%matrix)
      end if
      p%matrix = matrix
      p%energy = label%energy
      p%kvec = label%kvec
    end associate

  end subroutine mem_add_sp


  !> Retrieve with copy - double precision output
  subroutine mem_retrieve_dp(this, matrix, label)
    class(TMatrixCacheMem) :: this
    type(z_DNS) :: matrix
    type(TMatLabel) :: label

    integer :: b

    call checks(this, label, b)
    associate(p => this%MatArray(label%energy_point, label%kpoint, b, label%spin))
      if (associated(p%matrix)) then
        call destroy(matrix)
        matrix = p%matrix
      else
        print*, 'Retrieve matrix for '//trim(this%tagname)
        call print_label(label)
        call error_msg( "Cannot retrieve matrix")
      end if
    end associate

  end subroutine mem_retrieve_dp

  !> Retrieve with copy - single precision output
  subroutine mem_retrieve_sp(this, matrix, label)
    class(TMatrixCacheMem) :: this
    type(c_DNS) :: matrix
    type(TMatLabel) :: label

    integer :: b

    call checks(this, label, b)
    associate(p => this%MatArray(label%energy_point, label%kpoint, b, label%spin))
      if (associated(p%matrix)) then
        call destroy(matrix)
        matrix = p%matrix
      else
        print*, 'Retrieve matrix for '//trim(this%tagname)
        call print_label(label)
        call error_msg( "Cannot retrieve matrix")
      end if
    end associate

  end subroutine mem_retrieve_sp


  !> Retrieve memory pointer - single precision output
  subroutine mem_retrieve_pointer(this, pmatrix, label)
    class(TMatrixCacheMem) :: this
    type(TMatLabel) :: label
    type(x_DNS), pointer :: pmatrix

    integer :: b

    call checks(this, label, b)
    associate(p => this%MatArray(label%energy_point, label%kpoint, b, label%spin))
      if (associated(p%matrix)) then
        pmatrix => p%matrix
      else
        print*, 'Retrieve matrix for '//trim(this%tagname)
        call print_label(label)
        call error_msg( "Cannot retrieve matrix")
      end if
    end associate

  end subroutine mem_retrieve_pointer


  function mem_retrieve_loc(this, label) result(pmatrix)
    class(TMatrixCacheMem) :: this
    type(TMatLabel) :: label
    type(C_PTR) :: pmatrix

    integer :: b

    call checks(this, label, b)
    associate(p => this%MatArray(label%energy_point, label%kpoint, b, label%spin))
      if (associated(p%matrix)) then
        pmatrix = c_loc(p%matrix%val)
      else
        print*, 'Retrieve matrix for '//trim(this%tagname)
        call print_label(label)
        call error_msg( "Cannot retrieve matrix")
      end if
    end associate

  end function mem_retrieve_loc

  function mem_is_cached(this, label) result(val)
    class(TMatrixCacheMem) :: this
    type(TMatLabel) :: label
    logical :: val

    integer :: b

    val = this%in_range(label, b)
    if (.not.val) return
    associate(p => this%MatArray(label%energy_point, label%kpoint, b, label%spin))
      val = associated(p%matrix)
      if (val) val = same_point(p%energy, p%kvec, label)
    end associate

  end function mem_is_cached

  subroutine mem_list_cache(this)
    class(TMatrixCacheMem) :: this

  end subroutine mem_list_cache

  ! this%first => p%matrix     p%matrix      p%matrix
  !               p%label      p%label       p%label
  !               p%next   =>  p%next    =>  p%next
  subroutine mem_destroy(this)
    class(TMatrixCacheMem) :: this
    integer :: ii, jj, kk, ll, total_size

    if (allocated(this%MatArray)) then
       total_size = size(this%MatArray)

       do ll = 1, size(this%MatArray,4)
         do kk = 1, size(this%MatArray,3)
           do jj = 1, size(this%MatArray,2)
             do ii = 1, size(this%MatArray,1)
                if (associated(this%MatArray(ii,jj,kk,ll)%matrix)) then
                  call destroy(this%MatArray(ii,jj,kk,ll)%matrix)
                  deallocate(this%MatArray(ii,jj,kk,ll)%matrix)
                end if
             end do
           end do
         end do
       end do
       deallocate(this%MatArray)
       this%isInitialized = .false.
    end if
  end subroutine mem_destroy


  subroutine mem_clear(this)
    class(TMatrixCacheMem) :: this
    integer :: ii, jj, kk, ll, total_size

    if (allocated(this%MatArray)) then
       total_size = size(this%MatArray)

       !$omp parallel do private(ii,jj,kk,ll) collapse(4) if(total_size > 1000)
       do ll = 1, size(this%MatArray,4)
         do kk = 1, size(this%MatArray,3)
           do jj = 1, size(this%MatArray,2)
             do ii = 1, size(this%MatArray,1)

                if (associated(this%MatArray(ii,jj,kk,ll)%matrix)) then
                   ! Azzeramento generico senza dipendenze da _dp o _lp
                   this%MatArray(ii,jj,kk,ll)%matrix%val = &
                      cmplx(0.0, 0.0, kind=kind(this%MatArray(ii,jj,kk,ll)%matrix%val))
                end if

             end do
           end do
         end do
       end do
       !$omp end parallel do
    end if
  end subroutine mem_clear

  !! Definitions for TMatrixCacheDisk

  subroutine disk_add_dp(this, matrix, label)
    class(TMatrixCacheDisk) :: this
    type(z_DNS):: matrix
    type(TMatLabel) :: label

    character(LST) :: filename
    integer :: file_unit

    call disk_indices_to_filename(filename, label)
    open (newunit=file_unit, file=trim(this%scratch_path)//filename, form='UNFORMATTED')
    ! Header record with the point, checked by disk_is_cached
    write (file_unit) label%energy, label%kvec
    call outmat_c(file_unit, .false., matrix%val, matrix%nrow, matrix%ncol)
    close (file_unit)

  end subroutine

  subroutine disk_add_sp(this, matrix, label)
    class(TMatrixCacheDisk) :: this
    type(c_DNS):: matrix
    type(TMatLabel) :: label

    character(LST) :: filename
    integer :: file_unit

    call disk_indices_to_filename(filename, label)
    open (newunit=file_unit, file=trim(this%scratch_path)//filename, form='UNFORMATTED')
    ! call outmat_c(file_unit, .false., matrix%val, matrix%nrow, matrix%ncol)
    close (file_unit)

  end subroutine


  subroutine disk_retrieve_dp(this, matrix, label)
    class(TMatrixCacheDisk) :: this
    !> Input matrix to be retrieved: It needs to be already allocated with
    !> the correct number of rows and columns.
    type(z_DNS) :: matrix
    type(TMatLabel) :: label

    character(LST) :: filename
    logical :: file_exists
    integer :: file_unit

    call disk_indices_to_filename(filename, label)
    inquire (file=trim(this%scratch_path)//filename, EXIST=file_exists)

    if (.not. file_exists) then
      call error_msg( "Cannot retrieve matrix from disk: file not found")
    end if

    open (newunit=file_unit, file=trim(this%scratch_path)//filename, form='UNFORMATTED')
    ! Skip the header record
    read (file_unit)
    call inmat_c(file_unit, .false., matrix%val, matrix%nrow, matrix%ncol)
    close (file_unit)

  end subroutine disk_retrieve_dp

  subroutine disk_retrieve_sp(this, matrix, label)
    class(TMatrixCacheDisk) :: this
    !> Input matrix to be retrieved: It needs to be already allocated with
    !> the correct number of rows and columns.
    type(c_DNS) :: matrix
    type(TMatLabel) :: label

  end subroutine disk_retrieve_sp

  ! This routine gives error because a pointer cannot be associated
  subroutine disk_retrieve_pointer(this, pmatrix, label)
    class(TMatrixCacheDisk) :: this
    type(TMatLabel) :: label
    type(x_DNS), pointer :: pmatrix
    pmatrix => null()
    call error_msg( "cannot retrieve pointer from disk cache")
  end subroutine disk_retrieve_pointer

  ! This routine gives error because a pointer cannot be associated
  function disk_retrieve_loc(this, label) result(pmatrix)
    class(TMatrixCacheDisk) :: this
    type(TMatLabel) :: label
    type(C_PTR) :: pmatrix

    pmatrix = C_NULL_PTR
    call error_msg( "cannot retrieve pointer from disk cache")

  end function disk_retrieve_loc


  function disk_is_cached(this, label) result(val)
    class(TMatrixCacheDisk) :: this
    type(TMatLabel) :: label
    logical :: val

    character(LST) :: filename
    logical :: file_exists
    integer :: file_unit, ierr
    complex(dp) :: energy
    real(dp) :: kvec(3)

    call disk_indices_to_filename(filename, label)
    inquire (file=trim(this%scratch_path)//filename, EXIST=file_exists)

    val = .false.
    if (.not.file_exists) return

    ! Files without a valid header (e.g. from older versions) count as not cached
    open (newunit=file_unit, file=trim(this%scratch_path)//filename, form='UNFORMATTED', &
        & status='OLD', action='READ', iostat=ierr)
    if (ierr /= 0) return
    read (file_unit, iostat=ierr) energy, kvec
    close (file_unit)
    if (ierr /= 0) return
    val = same_point(energy, kvec, label)

  end function

  subroutine disk_list_cache(this)
    class(TMatrixCacheDisk) :: this

    ! no good. Probably the object should keep track of
    ! the cached matrices with a mapping label-> filename
  end subroutine


  subroutine disk_destroy(this)
    class(TMatrixCacheDisk) :: this
    ! Empty. Defined only for interface.
  end subroutine

  subroutine disk_clear(this)
    class(TMatrixCacheDisk) :: this
    ! Empty. Defined only for interface.
  end subroutine

  subroutine disk_indices_to_filename(filename, label)
    character(LST), intent(out) :: filename
    type(TMatLabel) :: label

    character(2) :: ofrow
    character(2) :: ofcol
    character(1) :: ofspin
    character(10) :: ofkpnt
    character(10) :: ofEpnt

    write (ofrow, '(i2.2)') label%row_block
    write (ofcol, '(i2.2)') label%col_block

    call get_string_label(label%kpoint, ofkpnt)
    call get_string_label(label%energy_point, ofEpnt)
    if (label%spin .eq. 1) ofspin = 'u'
    if (label%spin .eq. 2) ofspin = 'd'
    filename = 'Mat'//ofspin//ofrow//'_'//ofcol//'_'//trim(ofkpnt)//'_'//trim(ofEpnt)//'.dat'

  end subroutine

  !! Definitions for TMatrixCacheDummy

  subroutine dummy_add_dp(this, matrix, label)
    class(TMatrixCacheDummy) :: this
    type(z_DNS):: matrix
    type(TMatLabel) :: label
    ! Dummy operation
  end subroutine

  subroutine dummy_add_sp(this, matrix, label)
    class(TMatrixCacheDummy) :: this
    type(c_DNS):: matrix
    type(TMatLabel) :: label
    ! Dummy operation
  end subroutine

  subroutine dummy_retrieve_dp(this, matrix, label)
    class(TMatrixCacheDummy) :: this
    !> Input matrix to be retrieved: It needs to be already allocated with
    !> the correct number of rows and columns.
    type(z_DNS) :: matrix
    type(TMatLabel) :: label

  end subroutine

  subroutine dummy_retrieve_sp(this, matrix, label)
    class(TMatrixCacheDummy) :: this
    !> Input matrix to be retrieved: It needs to be already allocated with
    !> the correct number of rows and columns.
    type(c_DNS) :: matrix
    type(TMatLabel) :: label

  end subroutine

  ! This routine gives error because a pointer cannot be associated
  subroutine dummy_retrieve_pointer(this, pmatrix, label)
    class(TMatrixCacheDummy) :: this
    type(TMatLabel) :: label
    type(x_DNS), pointer :: pmatrix
    ! Dummy operation
    pmatrix => null()
  end subroutine dummy_retrieve_pointer


  ! This routine gives error because a pointer cannot be associated
  function dummy_retrieve_loc(this, label) result(pmatrix)
    class(TMatrixCacheDummy) :: this
    type(TMatLabel) :: label
    type(C_PTR) :: pmatrix
    ! Dummy operation
    pmatrix = C_NULL_PTR

  end function dummy_retrieve_loc


  function dummy_is_cached(this, label) result(val)
    class(TMatrixCacheDummy) :: this
    type(TMatLabel) :: label
    logical :: val

    val = .false.

  end function

  subroutine dummy_list_cache(this)
    class(TMatrixCacheDummy) :: this
  end subroutine

  subroutine dummy_destroy(this)
    class(TMatrixCacheDummy) :: this
  end subroutine

  subroutine dummy_clear(this)
    class(TMatrixCacheDummy) :: this
    ! Empty. Defined only for interface.
  end subroutine


  ! utility to obtain a 0-padded integer label
  subroutine get_string_label(int_label, str_label)
    integer, intent(in) :: int_label
    character(*), intent(inout) :: str_label
    if (int_label .le. 99) write (str_label, '(i2.2)') int_label
    if (int_label .gt. 99) write (str_label, '(i3.3)') int_label
    if (int_label .gt. 999) write (str_label, '(i4.4)') int_label
    if (int_label .gt. 9999) write (str_label, '(i5.5)') int_label
    if (int_label .gt. 99999) write (str_label, '(i6.6)') int_label
    if (int_label .gt. 999999) write (str_label, '(i7.7)') int_label
    if (int_label .gt. 9999999) call error_msg( 'ERROR: label too large (> 9999999)')
  end subroutine get_string_label

end module
