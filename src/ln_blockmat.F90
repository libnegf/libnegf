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
  public :: create_blockmat, destroy_blockmat, init_tridiagonal_blockmat, &
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

  interface init_tridiagonal_blockmat
    module procedure init_tridiagonal_blockmat_zdns
  end interface

  interface csr_to_tridiagonal_blockmat
    module procedure zcsr_to_tridiagonal_blockmat_zdns
  end interface

  interface subtract_from_block
    module procedure subtract_from_block_zdns
  end interface

contains

!> Create a block matrix. The individual blocks are not allocated here.
  subroutine create_blockmat_zdns(matrix, nrow)
    !> The matrix to be created.
    type(TSquareBlockZDns), intent(out) :: matrix
    !> The number of block rows.
    integer, intent(in) :: nrow

    if (allocated(matrix%blocks)) then
      call destroy_blockmat(matrix)
    end if

    allocate (matrix%blocks(nrow, nrow))
    matrix%nrow = nrow

  end subroutine create_blockmat_zdns

!> Destroy a block matrix.
  subroutine destroy_blockmat_zdns(matrix)
    !> The matrix to be created.
    type(TSquareBlockZDns), intent(out) :: matrix

    if (allocated(matrix%blocks)) then
      deallocate (matrix%blocks)
    end if
    matrix%nrow = 0

  end subroutine destroy_blockmat_zdns

!> Initialize a tridiagonal block matrix with the same blocks
!> as a sample one, but filled with zeros.
  subroutine init_tridiagonal_blockmat_zdns(sample_matrix, matrix)
    !> The matrix which provides the tridiagonal blocks size
    type(TSquareBlockZDns), intent(in) :: sample_matrix
    !> The matrix to initialize
    type(TSquareBlockZDns), intent(out) :: matrix

    integer :: nbl, j

    call create_blockmat(matrix, sample_matrix%nrow)

    nbl = sample_matrix%nrow

    associate (bl => sample_matrix%blocks)

      call create(matrix%blocks(1, 1), bl(1, 1)%nrow, bl(1, 1)%ncol)
      matrix%blocks(1, 1)%val = (0.0_dp, 0.0_dp)
      do j = 2, nbl - 1
        call create(matrix%blocks(j - 1, j), bl(j - 1, j)%nrow, bl(j - 1, j)%ncol)
        matrix%blocks(j - 1, j)%val = (0.0_dp, 0.0_dp)
        call create(matrix%blocks(j, j), bl(j, j)%nrow, bl(j, j)%ncol)
        matrix%blocks(j, j)%val = (0.0_dp, 0.0_dp)
        call create(matrix%blocks(j, j - 1), bl(j, j - 1)%nrow, bl(j, j - 1)%ncol)
        matrix%blocks(j, j - 1)%val = (0.0_dp, 0.0_dp)
      end do
      if (nbl .gt. 1) then
        call create(matrix%blocks(nbl, nbl), bl(nbl, nbl)%nrow, bl(nbl, nbl)%ncol)
        matrix%blocks(nbl, nbl)%val = (0.0_dp, 0.0_dp)
        call create(matrix%blocks(nbl - 1, nbl), bl(nbl - 1, nbl)%nrow, bl(nbl - 1, nbl)%ncol)
        matrix%blocks(nbl - 1, nbl)%val = (0.0_dp, 0.0_dp)
        call create(matrix%blocks(nbl, nbl - 1), bl(nbl, nbl - 1)%nrow, bl(nbl, nbl - 1)%ncol)
        matrix%blocks(nbl, nbl - 1)%val = (0.0_dp, 0.0_dp)
      end if

    end associate

  end subroutine init_tridiagonal_blockmat_zdns

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

    call create_blockmat(block_matrix, csr_matrix%nrow)

    nbl = size(indices) - 1

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
  subroutine subtract_from_block_zdns(block_matrix, dns_matrix, row, column)
    !> The block matrix to subtract from.
    type(TSquareBlockZDns), intent(inout) :: block_matrix
    !> The dense matrix to subtract.
    type(z_DNS), intent(in) :: dns_matrix
    !> The row index of the block.
    integer, intent(in) :: row
    !> The column index of the block.
    integer, intent(in) :: column

    block_matrix%blocks(row, column)%val = block_matrix%blocks(row, column)%val - dns_matrix%val

  end subroutine subtract_from_block_zdns

end module ln_blockmat
