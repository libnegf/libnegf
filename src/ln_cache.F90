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
  use iso_c_binding
  use mat_def, only: z_DNS, destroy, assignment(=)
  use ln_precision
  use globals
  use outmatrix, only: outmat_c, inmat_c

  implicit none
  private

  public :: TMatrixCache
  public :: TMatrixCacheMem, TMatrixCacheDisk
  public :: TMatrixCacheDummy
  public :: TMatLabel
  public :: print_label
  public :: get_string_label

  type TMatLabel
    integer :: kpoint = -1
    integer :: energy_point = -1
    integer :: spin = -1
    integer :: row_block = -1
    integer :: col_block = -1
  end type TMatLabel

  interface operator (==)
    module procedure labeleq
  end interface

  type, abstract :: TMatrixCache
  contains
    procedure(abst_add_matrix), deferred :: add
    procedure(abst_retrieve_matrix), deferred :: retrieve
    procedure(abst_retrieve_matrix_pointer), deferred :: retrieve_pointer
    procedure(abst_retrieve_matrix_loc), deferred :: retrieve_loc
    procedure(abst_destroy_matrix), deferred :: destroy
    procedure(abst_is_cached), deferred :: is_cached
    procedure(abst_list_cache), deferred :: list_cache
  end type TMatrixCache

  abstract interface
    subroutine abst_add_matrix(this, matrix, label)
      import :: TMatrixCache
      import :: z_DNS
      import :: TMatLabel
      class(TMatrixCache) :: this
      type(z_DNS):: matrix
      type(TMatLabel) :: label
    end subroutine

    subroutine abst_retrieve_matrix(this, matrix, label)
      import :: TMatrixCache
      import :: z_DNS
      import :: TMatLabel
      class(TMatrixCache) :: this
      type(z_DNS) :: matrix
      type(TMatLabel) :: label
    end subroutine

    function abst_retrieve_matrix_pointer(this, label) result(pmatrix)
      import :: z_DNS
      import :: TMatrixCache
      import :: TMatLabel
      class(TMatrixCache) :: this
      type(TMatLabel) :: label
      type(z_DNS), pointer :: pmatrix
    end function 

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
  end interface

  !> Linked list to store matrix in memory
  type :: TMatrixCacheEntry
    type(TMatrixCacheEntry), pointer :: next => null()
    type(z_DNS), allocatable :: matrix
    type(TMatLabel) :: mat_label
  end type

  !> Mem derived class
  type, extends(TMatrixCache) :: TMatrixCacheMem
    type(TMatrixCacheEntry), pointer :: first => null()
    character(len=LST) :: tagname 
  contains
    procedure :: add => mem_add
    procedure :: retrieve => mem_retrieve
    procedure :: retrieve_pointer => mem_retrieve_pointer
    procedure :: retrieve_loc => mem_retrieve_loc
    procedure :: destroy => mem_destroy
    procedure :: is_cached => mem_is_cached
    procedure :: list_cache => mem_list_cache
  end type

  !> Disk derived class
  type, extends(TMatrixCache) :: TMatrixCacheDisk
    character(len=LST) :: scratch_path
  contains
    procedure :: add => disk_add
    procedure :: retrieve => disk_retrieve
    procedure :: retrieve_pointer => disk_retrieve_pointer
    procedure :: retrieve_loc => disk_retrieve_loc
    procedure :: destroy => disk_destroy
    procedure :: is_cached => disk_is_cached
    procedure :: list_cache => disk_list_cache
  end type

  !> Dummy cache when no caching is performed
  type, extends(TMatrixCache) :: TMatrixCacheDummy
  contains
    procedure :: add => dummy_add
    procedure :: retrieve => dummy_retrieve
    procedure :: retrieve_pointer => dummy_retrieve_pointer
    procedure :: retrieve_loc => dummy_retrieve_loc
    procedure :: destroy => dummy_destroy
    procedure :: is_cached => dummy_is_cached
    procedure :: list_cache => dummy_list_cache
  end type

contains

  function labeleq(label1, label2) result(var)
     type(TMatLabel), intent(in) :: label1, label2
     logical :: var
     var = .false.
     if (label1%kpoint        == label2%kpoint        .and. &
       & label1%spin          == label2%spin          .and. &
       & label1%row_block     == label2%row_block     .and. &
       & label1%col_block     == label2%col_block     .and. &
       & label1%energy_point  == label2%energy_point ) then
          var = .true.
     end if
  end function labeleq

  subroutine print_label(label)
    type(TMatLabel) :: label

    print*,'k=',label%kpoint,'E=',label%energy_point,'s=',label%spin, &
          & 'i=',label%row_block,'j=',label%col_block
      
  end subroutine print_label

  !! Definitions for TMatrixCacheMem

  subroutine mem_add(this, matrix, label)
    class(TMatrixCacheMem) :: this
    type(z_DNS):: matrix
    type(TMatLabel) :: label

    type(TMatrixCacheEntry), pointer :: p

    if (.not. associated(this%first)) then
      allocate (this%first)
      p => this%first
    else
      ! check if a label is already present.
      do while (associated(p))
        if (p%mat_label .eq. label) then
          call destroy(p%matrix)
          p%matrix = matrix
          return
        end if
        p => p%next
      end do
      ! If the matrix is not found, add the new matrix
      ! Always add at the beginning, because it is faster.
      allocate(p)
      p%next => this%first
      this%first => p
    end if
    allocate(p%matrix)
    p%matrix = matrix
    p%mat_label = label
  end subroutine

  subroutine mem_retrieve(this, matrix, label)
    class(TMatrixCacheMem) :: this
    type(z_DNS) :: matrix
    type(TMatLabel) :: label

    type(TMatrixCacheEntry), pointer :: p

    if (.not. associated(this%first)) then
      print*, 'Retrieve matrix for '//trim(this%tagname)    
      call print_label(label)    
      error stop "Internal error: no entry in matrix cache"
    else
      p => this%first
    end if
    do while (associated(p))
      if (p%mat_label .eq. label) then
        call destroy(matrix)
        matrix = p%matrix
        return
      end if
      p => p%next
    end do

    print*, 'Retrieve matrix for '//trim(this%tagname)    
    call print_label(label)    
    error stop "Cannot retrieve matrix"

  end subroutine
  
  function mem_retrieve_pointer(this, label) result(pmatrix)
    class(TMatrixCacheMem) :: this
    type(TMatLabel) :: label
    type(z_DNS), pointer :: pmatrix  
  
    type(TMatrixCacheEntry), pointer :: p

    if (.not. associated(this%first)) then
      print*, trim(this%tagname)    
      error stop "No entry in matrix cache"
    else
      p => this%first
    end if
    do while (associated(p))
      if (p%mat_label .eq. label) then
        pmatrix => p%matrix
        return 
      end if
      p => p%next
    end do

    print*, 'Retrieve pointer for '//trim(this%tagname)    
    call print_label(label)    
    error stop "Cannot retrieve matrix in cache"

  end function 
  
  function mem_retrieve_loc(this, label) result(pmatrix)
    class(TMatrixCacheMem) :: this
    type(TMatLabel) :: label
    type(C_PTR) :: pmatrix
    
    type(TMatrixCacheEntry), pointer :: p

    if (.not. associated(this%first)) then
      error stop "No entry in matrix cache"
    else
      p => this%first
    end if

    do while (associated(p))
      if (p%mat_label .eq. label) then
        pmatrix = c_loc(p%matrix%val)
        return 
      end if
      p => p%next
    end do

    print*, 'Retrieve loc for '//trim(this%tagname)    
    call print_label(label)
    error stop "Cannot retrieve matrix in cache"

  end function 
  
  function mem_is_cached(this, label) result(val)
    class(TMatrixCacheMem) :: this
    type(TMatLabel) :: label
    logical :: val

    type(TMatrixCacheEntry), pointer :: p

    if (.not. associated(this%first)) then
      val = .false.
      return
    else
      p => this%first
    end if

    do while (associated(p))
      if (p%mat_label .eq. label) then
        val = .true.
        return
      end if
      p => p%next
    end do

    val = .false.

  end function
    
  subroutine mem_list_cache(this)
    class(TMatrixCacheMem) :: this
    
    type(TMatLabel) :: label
    type(TMatrixCacheEntry), pointer :: p
    type(z_DNS), pointer :: pmatrix  
    
    if (.not. associated(this%first)) then
      return
    end if
    print*,'-------------------------------------'
    print*,'List cache:'
    p => this%first
    do while (associated(p))
      call print_label(p%mat_label) 
      pmatrix => p%matrix
      print*,'matrix size:',size(pmatrix%val)
      p => p%next
    end do
    print*,'-------------------------------------'

  end subroutine

  ! this%first => p%matrix     p%matrix      p%matrix
  !               p%label      p%label       p%label 
  !               p%next   =>  p%next    =>  p%next
  subroutine mem_destroy(this)
    class(TMatrixCacheMem) :: this

    type(TMatrixCacheEntry), pointer :: p, previous
    if (.not. associated(this%first)) then
      return
    end if

    p => this%first
    do while (associated(p))
      call destroy(p%matrix)
      deallocate(p%matrix)
      previous => p
      p => p%next
      deallocate(previous)
    end do
    this%first => null()

  end subroutine


  !! Definitions for TMatrixCacheDisk

  subroutine disk_add(this, matrix, label)
    class(TMatrixCacheDisk) :: this
    type(z_DNS):: matrix
    type(TMatLabel) :: label

    real(kind=dp) :: dens
    character(2) :: ofcont
    character(1) :: ofspin
    character(10) :: ofkpnt
    character(10) :: ofpnt
    character(LST) :: filename
    integer :: file_unit

    call disk_indices_to_filename(filename, label)
    open (newunit=file_unit, file=trim(this%scratch_path)//filename, form='UNFORMATTED')
    call outmat_c(file_unit, .false., matrix%val, matrix%nrow, matrix%ncol)
    close (file_unit)

  end subroutine

  subroutine disk_retrieve(this, matrix, label)
    class(TMatrixCacheDisk) :: this
    !> Input matrix to be retrieved: It needs to be already allocated with
    !> the correct number of rows and columns.
    type(z_DNS) :: matrix
    type(TMatLabel) :: label

    real(kind=dp) :: dens
    character(LST) :: filename
    logical :: file_exists
    integer :: file_unit

    call disk_indices_to_filename(filename, label)
    inquire (file=trim(this%scratch_path)//filename, EXIST=file_exists)

    if (.not. file_exists) then
      error stop "Cannot retrieve matrix from disk: file not found"
    end if

    open (newunit=file_unit, file=trim(this%scratch_path)//filename, form='UNFORMATTED')
    call inmat_c(file_unit, .false., matrix%val, matrix%nrow, matrix%ncol)
    close (file_unit)

  end subroutine
  
  ! This routine gives error because a pointer cannot be associated
  function disk_retrieve_pointer(this, label) result(pmatrix)
    class(TMatrixCacheDisk) :: this
    type(TMatLabel) :: label
    type(z_DNS), pointer :: pmatrix

    pmatrix => null() 
    error stop "cannot retrieve pointer from disk cache"

  end function disk_retrieve_pointer
  
  ! This routine gives error because a pointer cannot be associated
  function disk_retrieve_loc(this, label) result(pmatrix)
    class(TMatrixCacheDisk) :: this
    type(TMatLabel) :: label
    type(C_PTR) :: pmatrix

    pmatrix = C_NULL_PTR
    error stop "cannot retrieve pointer from disk cache"

  end function disk_retrieve_loc


  function disk_is_cached(this, label) result(val)
    class(TMatrixCacheDisk) :: this
    type(TMatLabel) :: label
    logical :: val

    character(LST) :: filename
    logical :: file_exists
    integer :: file_unit

    call disk_indices_to_filename(filename, label)
    inquire (file=trim(this%scratch_path)//filename, EXIST=file_exists)

    if (file_exists) then
      val = .true.
    else
      val = .false.
    end if

  end function
  
  subroutine disk_list_cache(this)
    class(TMatrixCacheDisk) :: this
    
    type(TMatLabel) :: label
    
    ! no good. Probably the object should keep track of
    ! the cached matrices with a mapping label-> filename
  end subroutine


  subroutine disk_destroy(this)
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

  subroutine dummy_add(this, matrix, label)
    class(TMatrixCacheDummy) :: this
    type(z_DNS):: matrix
    type(TMatLabel) :: label

    ! Dummy operation

  end subroutine

  subroutine dummy_retrieve(this, matrix, label)
    class(TMatrixCacheDummy) :: this
    !> Input matrix to be retrieved: It needs to be already allocated with
    !> the correct number of rows and columns.
    type(z_DNS) :: matrix
    type(TMatLabel) :: label

    ! Dummy operation

  end subroutine
  
  ! This routine gives error because a pointer cannot be associated
  function dummy_retrieve_pointer(this, label) result(pmatrix)
    class(TMatrixCacheDummy) :: this
    type(TMatLabel) :: label
    type(z_DNS), pointer :: pmatrix
    ! Dummy operation
    pmatrix => null() 

  end function dummy_retrieve_pointer
  
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
    if (int_label .gt. 9999999) stop 'ERROR: label too large (> 9999999)'
  end subroutine get_string_label

end module
