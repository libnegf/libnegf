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
  use mat_def, only: z_DNS, create, destroy
  implicit none
  private

  public :: TSquareBlockZDns
  public :: create_blockmat, destroy_blockmat, init_tridiagonal_blockmat

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

contains

!> Create a block matrix. The individual blocks are not allocated here.
  subroutine create_blockmat_zdns(matrix, nrow)
    !> The matrix to be created.
    type(TSquareBlockZDns), intent(out) :: matrix
    !> The number of block rows.
    integer, intent(in) :: nrow

    allocate (matrix%blocks(nrow, nrow))
    matrix%nrow = nrow

  end subroutine create_blockmat_zdns

!> Destroy a block matrix.
  subroutine destroy_blockmat_zdns(matrix)
    !> The matrix to be created.
    type(TSquareBlockZDns), intent(out) :: matrix

    if (allocated(matrix%blocks)) then
      deallocate (matrix%blocks)
    endif
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

    call destroy_blockmat(matrix)
    call create_blockmat(matrix, sample_matrix%nrow)

    nbl = sample_matrix%nrow

    call create(matrix%blocks(1, 1), sample_matrix%blocks(1, 1)%nrow, sample_matrix%blocks(1, 1)%ncol)
    matrix%blocks(1, 1)%val = (0.0_dp, 0.0_dp)
    do j = 2, nbl - 1
      call create(matrix%blocks(j - 1, j), sample_matrix%blocks(j - 1, j)%nrow, sample_matrix%blocks(j - 1, j)%ncol)
      matrix%blocks(j - 1, j)%val = (0.0_dp, 0.0_dp)
      call create(matrix%blocks(j, j), sample_matrix%blocks(j, j)%nrow, sample_matrix%blocks(j, j)%ncol)
      matrix%blocks(j, j)%val = (0.0_dp, 0.0_dp)
      call create(matrix%blocks(j, j - 1), sample_matrix%blocks(j, j - 1)%nrow, sample_matrix%blocks(j, j - 1)%ncol)
      matrix%blocks(j, j - 1)%val = (0.0_dp, 0.0_dp)
    end do
    if (nbl .gt. 1) then
      call create(matrix%blocks(nbl, nbl), sample_matrix%blocks(nbl, nbl)%nrow, sample_matrix%blocks(nbl, nbl)%ncol)
      matrix%blocks(nbl, nbl)%val = (0.0_dp, 0.0_dp)
      call create(matrix%blocks(nbl - 1, nbl), sample_matrix%blocks(nbl - 1, nbl)%nrow, sample_matrix%blocks(nbl - 1, nbl)%ncol)
      matrix%blocks(nbl - 1, nbl)%val = (0.0_dp, 0.0_dp)
      call create(matrix%blocks(nbl, nbl - 1), sample_matrix%blocks(nbl, nbl - 1)%nrow, sample_matrix%blocks(nbl, nbl - 1)%ncol)
      matrix%blocks(nbl, nbl - 1)%val = (0.0_dp, 0.0_dp)
    endif

  end subroutine init_tridiagonal_blockmat_zdns

end module ln_blockmat
