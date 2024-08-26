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
  use mat_def, only: z_DNS, c_DNS, x_DNS => z_DNS, destroy, assignment(=)
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
    procedure :: add_sp => mem_add_sp
    procedure :: add_dp => mem_add_dp
    procedure :: retrieve_sp => mem_retrieve_sp 
    procedure :: retrieve_dp => mem_retrieve_dp
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
    procedure :: add_sp => disk_add_sp
    procedure :: add_dp => disk_add_dp
    procedure :: retrieve_dp => disk_retrieve_dp
    procedure :: retrieve_sp => disk_retrieve_sp
    procedure :: retrieve_pointer => disk_retrieve_pointer
    procedure :: retrieve_loc => disk_retrieve_loc
    procedure :: destroy => disk_destroy
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

    if (allocated(this%MatArray)) then
      if (size(this%MatArray) /= NE*Nk*NPL*Ndiag*Nspin) then   
        call this%destroy()
      end if  
    end if  
    if (.not.allocated(this%MatArray)) then
      allocate(this%MatArray(NE,Nk,NPL*Ndiag,Nspin), stat=ierr)
      if (ierr /= 0) then    
        error stop "Allocation error of MatArray in init "//trim(this%tagname)
      end if   
    end if
    this%Nblocks = NPL    
    this%isInitialized = .true.
  end subroutine init

  subroutine checks(this, label, bl)
    class(TMatrixCacheMem) :: this
    type(TMatLabel), intent(in) :: label
    integer, intent(out) :: bl

    integer :: r, c
    
    if (.not.allocated(this%MatArray)) then
      error stop "MatArray not allocated for "//trim(this%tagname)
    end if
    if (label%energy_point < 1 .or. label%energy_point > size(this%MatArray,1)) then
       print*,trim(this%tagname)   
       call print_label(label)
       error stop "energy_point out of range"
    end if
    if (label%kpoint < 1 .or. label%kpoint > size(this%MatArray,2)) then
       print*,trim(this%tagname)   
       call print_label(label)
       error stop "k_point out of range"
    end if
    if (label%spin < 1 .or. label%spin > size(this%MatArray,4)) then
       print*,trim(this%tagname)   
       call print_label(label)
       error stop "spin out of range"
    end if
        
    r = label%row_block
    c = label%col_block
    if (c >= r) then
       bl = (c-r)*this%Nblocks + c
    else
       bl = 2*(r-c)*this%Nblocks + r   
    end if
    if (bl < 1 .or. bl > size(this%MatArray,3)) then
       print*,trim(this%tagname)   
       call print_label(label)
       print*, "bl=",bl
       error stop "block out of range"
    end if
  end subroutine checks

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
        error stop "Cannot retrieve matrix"
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
        error stop "Cannot retrieve matrix"
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
        error stop "Cannot retrieve matrix"
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
        error stop "Cannot retrieve matrix"
      end if
    end associate

  end function mem_retrieve_loc

  function mem_is_cached(this, label) result(val)
    class(TMatrixCacheMem) :: this
    type(TMatLabel) :: label
    logical :: val

    integer :: b
    
    call checks(this, label, b) 
    associate(p => this%MatArray(label%energy_point, label%kpoint, b, label%spin))
      if (associated(p%matrix)) then
         val = .true. 
      else      
         val = .false. 
      end if
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

    integer :: ii, jj, kk, ll
    
    if (allocated(this%MatArray)) then
      do ii = 1, size(this%MatArray,1)
        do jj = 1, size(this%MatArray,2)
          do kk = 1, size(this%MatArray,3)
            do ll = 1, size(this%MatArray,4)
               if (associated(this%MatArray(ii,jj,kk,ll)%matrix)) then
                  call destroy(this%MatArray(ii,jj,kk,ll)%matrix)
               end if
            end do
          end do 
        end do
      end do
      deallocate(this%MatArray)
      this%isInitialized = .false.
    end if
  end subroutine mem_destroy
   

  !! Definitions for TMatrixCacheDisk

  subroutine disk_add_dp(this, matrix, label)
    class(TMatrixCacheDisk) :: this
    type(z_DNS):: matrix
    type(TMatLabel) :: label

    character(LST) :: filename
    integer :: file_unit

    call disk_indices_to_filename(filename, label)
    open (newunit=file_unit, file=trim(this%scratch_path)//filename, form='UNFORMATTED')
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
      error stop "Cannot retrieve matrix from disk: file not found"
    end if

    open (newunit=file_unit, file=trim(this%scratch_path)//filename, form='UNFORMATTED')
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
    error stop "cannot retrieve pointer from disk cache"
  end subroutine disk_retrieve_pointer

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
    if (int_label .gt. 9999999) error stop 'ERROR: label too large (> 9999999)'
  end subroutine get_string_label

end module
