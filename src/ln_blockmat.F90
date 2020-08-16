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
            csr_to_tridiagonal_blockmat, subtract_from_block

!> A square block matrix with blocks of type z_DNS.
  type TSquareBlockZDns
    type(z_DNS), dimension(:, :), allocatable :: blocks
    integer :: nrow = 0
  end type TSquareBlockZDns

  interface create_blockmat
    module procedure create_blockmat_zdns
  end interface

  interface destroy_blockmat
    module procedure destroy_blockmat_zdns
  end interface

  interface create_tridiagonal_blockmat
    module procedure create_tridiagonal_blockmat_zdns
  end interface

  interface csr_to_tridiagonal_blockmat
    module procedure zcsr_to_tridiagonal_blockmat_zdns
  end interface

  interface subtract_from_block
    module procedure subtract_from_block_zdns
  end interface

contains

  !> Create a block matrix. The individual blocks are not allocated here.
  subroutine create_blockmat_zdns(this, nrow)
    !> The matrix to be created.
    type(TSquareBlockZDns), intent(out) :: this
    !> The number of block rows.
    integer, intent(in) :: nrow

    if (allocated(this%blocks)) then
      call destroy_blockmat(this)
    end if

    allocate (this%blocks(nrow, nrow))
    this%nrow = nrow

  end subroutine create_blockmat_zdns

  !> Destroy a block matrix.
  subroutine destroy_blockmat_zdns(this)
    !> The matrix to be created.
    type(TSquareBlockZDns), intent(inout) :: this

    if (allocated(this%blocks)) then
      deallocate (this%blocks)
    end if
    this%nrow = 0

  end subroutine destroy_blockmat_zdns

  !> Initialize a tridiagonal block matrix with the same blocks
  !> as a sample one, but filled with zeros.
  subroutine create_tridiagonal_blockmat_zdns(this, sample_matrix)
    !> The matrix which provides the tridiagonal blocks size
    type(TSquareBlockZDns), intent(in) :: sample_matrix
    !> The matrix to initialize
    type(TSquareBlockZDns), intent(out) :: this

    integer :: nbl, j

    call create_blockmat(this, sample_matrix%nrow)

    nbl = sample_matrix%nrow

    associate (bl => sample_matrix%blocks)

      call create(this%blocks(1, 1), bl(1, 1)%nrow, bl(1, 1)%ncol)
      this%blocks(1, 1)%val = (0.0_dp, 0.0_dp)
      do j = 2, nbl - 1
        call create(this%blocks(j - 1, j), bl(j - 1, j)%nrow, bl(j - 1, j)%ncol)
        this%blocks(j - 1, j)%val = (0.0_dp, 0.0_dp)
        call create(this%blocks(j, j), bl(j, j)%nrow, bl(j, j)%ncol)
        this%blocks(j, j)%val = (0.0_dp, 0.0_dp)
        call create(this%blocks(j, j - 1), bl(j, j - 1)%nrow, bl(j, j - 1)%ncol)
        this%blocks(j, j - 1)%val = (0.0_dp, 0.0_dp)
      end do
      if (nbl .gt. 1) then
        call create(this%blocks(nbl, nbl), bl(nbl, nbl)%nrow, bl(nbl, nbl)%ncol)
        this%blocks(nbl, nbl)%val = (0.0_dp, 0.0_dp)
        call create(this%blocks(nbl - 1, nbl), bl(nbl - 1, nbl)%nrow, bl(nbl - 1, nbl)%ncol)
        this%blocks(nbl - 1, nbl)%val = (0.0_dp, 0.0_dp)
        call create(this%blocks(nbl, nbl - 1), bl(nbl, nbl - 1)%nrow, bl(nbl, nbl - 1)%ncol)
        this%blocks(nbl, nbl - 1)%val = (0.0_dp, 0.0_dp)
      end if

    end associate

  end subroutine create_tridiagonal_blockmat_zdns

  !> Convert a CSR matrix to a tridiagonal block square. Entries out of the
  !> tridiagonal blocks are ignored. The block matrix is constructed internally.
  subroutine zcsr_to_tridiagonal_blockmat_zdns(csr_matrix, indices, block_matrix)
    !> The CSR matrix to copy values from.
    type(z_CSR), intent(in) :: csr_matrix
    !> The starting index of each box. The last index must be the final row + 1.
    integer, dimension(:), allocatable, intent(in) :: indices
    !> The output block matrix.
    type(TSquareBlockZDns), intent(out) :: block_matrix

    integer :: i, nbl

    nbl = size(indices) - 1

    call create_blockmat(block_matrix, nbl)

    do i = 1, nbl
      call extract(csr_matrix, indices(i), indices(i + 1) - 1, indices(i), &
          & indices(i + 1) - 1, block_matrix%blocks(i, i))
    end do

    do i = 2, nbl
      call extract(csr_matrix, indices(i - 1), indices(i) - 1, indices(i), &
          & indices(i + 1) - 1, block_matrix%blocks(i - 1, i))
      call extract(csr_matrix, indices(i), indices(i + 1) - 1, indices(i - 1), &
          & indices(i) - 1, block_matrix%blocks(i, i - 1))
    end do

  end subroutine zcsr_to_tridiagonal_blockmat_zdns

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
