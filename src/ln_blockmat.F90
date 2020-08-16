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

module ln_blockmat
  use ln_precision
  use ln_allocation
  use mat_def, only: z_CSR, z_DNS, create, destroy
  use sparsekit_drv, only: extract
  implicit none
  private

  public :: TSquareBlockZDns
  public :: create_blockmat, destroy_blockmat, create_tridiagonal_blockmat, &
            subtract_from_block

!> A square block matrix with blocks of type z_DNS.
  type TSquareBlockZDns
    type(z_DNS), dimension(:, :), allocatable :: blocks
    integer :: nblocks = 0
    integer, dimension(:), allocatable :: blockind
  end type TSquareBlockZDns

  !> Empty block matrix constructor
  interface create_blockmat
    module procedure create_blockmat_zdns
  end interface

  !> Destructor
  interface destroy_blockmat
    module procedure destroy_blockmat_zdns
  end interface

    !> Initialized tridiagonal block matrix constructor
  interface create_tridiagonal_blockmat
    module procedure create_zeroed_tridiagonal_zdns
    module procedure create_tridiagonal_from_zcsr_zdns
  end interface

  !> Subtract a matrix from a specific block.
  interface subtract_from_block
    module procedure subtract_from_block_zdns
  end interface

contains

  !> Create a block matrix. The individual blocks are not allocated here.
  subroutine create_blockmat_zdns(this, blockind)
    !> The matrix to be created.
    type(TSquareBlockZDns), intent(out) :: this
    !> The array with the starting index of each diagonal block. The last element
    !> must contain the end of the last block + 1.
    integer, dimension(:), intent(in) :: blockind

    if (allocated(this%blocks)) then
      call destroy_blockmat(this)
    end if

    this%nblocks = size(blockind) - 1
    allocate (this%blocks(this%nblocks, this%nblocks))
    this%blockind = blockind

  end subroutine create_blockmat_zdns

  !> Destroy a block matrix.
  subroutine destroy_blockmat_zdns(this)
    !> The matrix to be created.
    type(TSquareBlockZDns), intent(inout) :: this

    if (allocated(this%blocks)) then
      deallocate (this%blocks)
    end if
    if (allocated(this%blockind)) then
      deallocate (this%blockind)
    end if
    this%nblocks = 0

  end subroutine destroy_blockmat_zdns

  !> Initialize a block to a zero matrix.
  subroutine init_block_to_zero(this, i, j)
    !> The block matrix.
    type(TSquareBlockZDns), intent(inout) :: this
    !> The block row index.
    integer, intent(in) :: i
    !> The coumn block index.
    integer, intent(in) :: j

    integer nrow, ncol

    nrow = this%blockind(i + 1) - this%blockind(i)
    ncol = this%blockind(i + 1) - this%blockind(i)

    call create(this%blocks(i, j), nrow, ncol)
    this%blocks(i, j)%val = (0.d0, 0.d0)

  end subroutine init_block_to_zero

  !> Initialize a tridiagonal block matrix with the same blocks
  !> as a sample one, but filled with zeros.
  subroutine create_zeroed_tridiagonal_zdns(this, blockind)
    !> The matrix to initialize
    type(TSquareBlockZDns), intent(out) :: this
    !> The array with the starting index of each diagonal block.
    integer, dimension(:), intent(in) :: blockind

    integer :: nbl, j

    call create_blockmat(this, blockind)

    nbl = this%nblocks

    call init_block_to_zero(this, 1, 1)
    do j = 2, nbl - 1
      call init_block_to_zero(this, j - 1, j)
      call init_block_to_zero(this, j, j)
      call init_block_to_zero(this, j, j - 1)
    end do
    if (nbl .gt. 1) then
      call init_block_to_zero(this, nbl, nbl)
    end if

  end subroutine create_zeroed_tridiagonal_zdns

  !> Convert a CSR matrix to a tridiagonal block square. Entries out of the
  !> tridiagonal blocks are ignored. The block matrix is constructed internally.
  subroutine create_tridiagonal_from_zcsr_zdns(this, csr_matrix, indices)
    !> The output block matrix.
    type(TSquareBlockZDns), intent(inout) :: this
    !> The CSR matrix to copy values from.
    type(z_CSR), intent(in) :: csr_matrix
    !> The starting index of each block. The last index must be the final row + 1.
    integer, dimension(:), allocatable, intent(in) :: indices

    integer :: i, nbl

    call create_blockmat(this, indices)

    do i = 1, this%nblocks
      call extract(csr_matrix, indices(i), indices(i + 1) - 1, indices(i), &
          & indices(i + 1) - 1, this%blocks(i, i))
    end do

    do i = 2, this%nblocks
      call extract(csr_matrix, indices(i - 1), indices(i) - 1, indices(i), &
          & indices(i + 1) - 1, this%blocks(i - 1, i))
      call extract(csr_matrix, indices(i), indices(i + 1) - 1, indices(i - 1), &
          & indices(i) - 1, this%blocks(i, i - 1))
    end do

  end subroutine create_tridiagonal_from_zcsr_zdns

  !> Subtract a dense matrix from a block, e.g. block_matrix(row, column) -= dns_matrix.
  subroutine subtract_from_block_zdns(this, dns_matrix, row, column)
    !> The block matrix to subtract from.
    type(TSquareBlockZDns), intent(inout) :: this
    !> The dense matrix to subtract.
    type(z_DNS), intent(in) :: dns_matrix
    !> The row index of the block.
    integer, intent(in) :: row
    !> The column index of the block.
    integer, intent(in) :: column

    this%blocks(row, column)%val = this%blocks(row, column)%val - dns_matrix%val

  end subroutine subtract_from_block_zdns

  !TODO: FINISH ME
  !> Convert a block matrix to a CSR folowing a given non-zero values pattern
  !> Note: this has been probably been used only on tri-diagonal matrices.
  subroutine blockmat_to_masked_csr(G, indices, P, Gcsr)

    type(z_DNS), dimension(:,:) :: G
    integer, dimension(:), intent(in) :: indices
    type(z_CSR) :: Gcsr
    type(z_CSR) :: P, G_sp

    integer :: nbl, oldx, row, col, iy, ix, x, y, ii, jj, nrows

    nbl = size(indices)
    nrows = indices(nbl) - 1

    !create Gcsr with same pattern of P
    call create(Gcsr,P%nrow,P%ncol,P%nnz)
    Gcsr%rowpnt = P%rowpnt
    Gcsr%colind = P%colind
    Gcsr%nzval = (0.0_dp, 0.0_dp)

    associate(indblk=>indices)
    !Cycle upon all rows
    x = 1
    do ii = 1, nrows
      !Search block x containing row ii
      oldx = x
      if (oldx.EQ.nbl) THEN
        x = oldx
      ELSE
        do ix = oldx, oldx+1
          if ( (ii.GE.indblk(ix)).AND.(ii.LT.indblk(ix+1)) ) x = ix
        end do
      endif

      !Offset: row is the index for separate blocks
      row = ii - indblk(x) + 1

      !Cycle upon columns of Gcsr (which has been ALREADY MASKED by S)
      do jj = Gcsr%rowpnt(ii), Gcsr%rowpnt(ii+1) -1
        if (Gcsr%colind(jj).gt.nrows) CYCLE
        !Choose which block column we're dealing with
        y = 0
        if (x.eq.1) then
          if ( (Gcsr%colind(jj).GE.indblk(x)).AND.(Gcsr%colind(jj).LT.indblk(x + 1)) ) then
            y = 1
          ELSEif ( (Gcsr%colind(jj).GE.indblk(x + 1)).AND.(Gcsr%colind(jj).LT.indblk(x + 2)) ) then
            y = 2
          endif
        elseif (x.eq.nbl) then
          if ( (Gcsr%colind(jj).GE.indblk(x)).AND.(Gcsr%colind(jj).LT.indblk(x + 1)) ) then
            y = nbl
          ELSEif ( (Gcsr%colind(jj).GE.indblk(x - 1)).AND.(Gcsr%colind(jj).LT.indblk(x)) ) then
            y = nbl - 1
          endif
        else
          do iy = x-1, x+1
            if ( (Gcsr%colind(jj).GE.indblk(iy)).AND.(Gcsr%colind(jj).LT.indblk(iy + 1)) ) y = iy
          end do
        endif

        if (y.EQ.0) THEN
          write(*,*)
          write(*,*) 'ERROR in blk2csr: probably wrong PL size', x
          write(*,*) 'row',ii,nrows,Gcsr%colind(jj)
          write(*,*) 'block indeces:',indblk(1:nbl)
          STOP
        endif

        col = Gcsr%colind(jj) - indblk(y) + 1

        if (allocated(G(x,y)%val)) THEN
          Gcsr%nzval(jj) = Gcsr%nzval(jj) + G(x,y)%val(row,col)
        endif

      end do

    end do
    end associate

  end subroutine blockmat_to_masked_csr

end module ln_blockmat
