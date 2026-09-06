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


Module block_3_diag
  use ln_precision
  use ln_allocation
  use ln_structure, only: Tstruct_Info
  use mat_def
  use iso_c_binding
  use mpi_globals
  use sparsekit_drv, only: extract
  implicit none
  private

public :: init_tridiag_blk
public :: destroy_tridiag_blk
public :: blk2dns
public :: dns2blk
public :: destroy_full_blk
public :: csr2blk3diag

interface csr2blk3diag
   module procedure zcsr2blk3diag
end interface

interface destroy_tridiag_blk
   module procedure destroy_tridiag_blk_2d
   module procedure destroy_tridiag_blk_3d
end interface

contains

subroutine destroy_tridiag_blk_2d(M, str)
   type(z_DNS), dimension(:,:), allocatable :: M
   character(len=*), intent(in), optional :: str

   integer :: i, nbl

   if (.not.allocated(M)) return

   nbl=size(M,1)

   do i=1,nbl
     if (allocated(M(i,i)%val)) then
       call destroy(M(i,i))
     end if
   end do
   do i=2,nbl
     if (allocated(M(i-1,i)%val)) then
       call destroy(M(i-1,i))
     end if
     if (allocated(M(i,i-1)%val)) then
       call destroy(M(i,i-1))
     end if
   end do

 end subroutine destroy_tridiag_blk_2d

 subroutine destroy_tridiag_blk_3d(M, str)
   type(z_DNS), dimension(:,:,:), allocatable :: M
   character(len=*), intent(in), optional :: str

   integer :: i, k, nbl, nK

   if (.not.allocated(M)) return

   nbl=size(M,2)
   nK=size(M,1)

   do k = 1,nK
      do i=1,nbl
      if (allocated(M(k,i,i)%val)) then
         call destroy(M(k,i,i))
      end if
      end do
      do i=2,nbl
      if (allocated(M(k,i-1,i)%val)) then
         call destroy(M(k,i-1,i))
      end if
      if (allocated(M(k,i,i-1)%val)) then
         call destroy(M(k,i,i-1))
      end if
      end do
   end do

 end subroutine destroy_tridiag_blk_3d

 subroutine destroy_full_blk(M, str)
   type(z_DNS), dimension(:,:), allocatable :: M
   character(len=*), intent(in), optional :: str

   integer :: i, j, nbl

   if (.not.allocated(M)) return

   nbl=size(M,1)

   do i=1,nbl
     do j=1,nbl
       if (allocated(M(i,j)%val)) then
         call destroy(M(i,j))
       end if
     end do
   end do

 end subroutine destroy_full_blk

 subroutine init_tridiag_blk(Matrix,S)
   type(z_DNS), dimension(:,:) :: Matrix,S

   integer :: nbl, j

   nbl = SIZE(Matrix,1)

#:if defined("GPU")
   call createAll(Matrix(1,1),S(1,1)%nrow,S(1,1)%ncol)
   Matrix(1,1)%val=(0.0_dp,0.0_dp)
   if (nbl.gt.1) then
      do j=2,nbl
        call createAll(Matrix(j-1,j),S(j-1,j)%nrow,S(j-1,j)%ncol)
        Matrix(j-1,j)%val=(0.0_dp,0.0_dp)
        call createAll(Matrix(j,j),S(j,j)%nrow,S(j,j)%ncol)
        Matrix(j,j)%val=(0.0_dp,0.0_dp)
        call createAll(Matrix(j,j-1),S(j,j-1)%nrow,S(j,j-1)%ncol)
        Matrix(j,j-1)%val=(0.0_dp,0.0_dp)
      end do
   endif
#:else
   ! if (allocated(Matrix(1,1)%val)) call destroy(Matrix(1,1))
   call create(Matrix(1,1),S(1,1)%nrow,S(1,1)%ncol)
   Matrix(1,1)%val=(0.0_dp,0.0_dp)
   if (nbl.gt.1) then
      do j=2,nbl
      !   if (allocated(Matrix(j-1,j)%val)) call destroy(Matrix(j-1,j))
        call create(Matrix(j-1,j),S(j-1,j)%nrow,S(j-1,j)%ncol)
        Matrix(j-1,j)%val=(0.0_dp,0.0_dp)
      !   if (allocated(Matrix(j,j)%val)) call destroy(Matrix(j,j))
        call create(Matrix(j,j),S(j,j)%nrow,S(j,j)%ncol)
        Matrix(j,j)%val=(0.0_dp,0.0_dp)
      !   if (allocated(Matrix(j,j-1)%val)) call destroy(Matrix(j,j-1))
        call create(Matrix(j,j-1),S(j,j-1)%nrow,S(j,j-1)%ncol)
        Matrix(j,j-1)%val=(0.0_dp,0.0_dp)
      end do
   endif
#:endif
 end subroutine init_tridiag_blk

 ! TRIDIAG blk-sparse to dense conversion
 ! CAREFUL! Does not work for generic blk matrices,
 ! but only if the dimensions of Gdns are the same as str%central_dim = nrow(H) = ncol(H)
 subroutine blk2dns(G,str,Gdns)
   type(z_DNS), dimension(:,:), intent(in) :: G
   type(Tstruct_info), intent(in) :: str
   type(z_DNS) :: Gdns

   integer :: n, m, ii, jj, pl_start1, pl_end1, pl_start2, pl_end2, nbl, nrows

   nbl = str%num_PLs
   nrows = str%central_dim
   call create(Gdns,nrows,nrows)

   do n = 1, nbl
     do m= 1, nbl
        pl_start1 = str%mat_PL_start(n)
        pl_end1 = str%mat_PL_end(n)
        pl_start2 = str%mat_PL_start(m)
        pl_end2 = str%mat_PL_end(m)

        if (allocated(G(n,m)%val)) then
          do ii = 1, pl_end1 - pl_start1 + 1
            do jj = 1, pl_end2 - pl_start2 + 1
                  ! if (id0) print*, "DEBUG: check2.5"
               Gdns%val(pl_start1 + ii - 1,pl_start2 + jj - 1) = G(n,m)%val(ii,jj)
            end do
          end do
        else
         do ii = 1, pl_end1 - pl_start1 + 1
            do jj = 1, pl_end2 - pl_start2 + 1
               Gdns%val(pl_start1 + ii - 1,pl_start2 + jj - 1) = 0.0_dp
            end do
          end do
        end if

     end do
   end do

 end subroutine blk2dns

 ! dense to blk-sparse conversion
 ! CAREFUL! Does not work for generic tri-diagonal blk matrices,
 ! but only if the dimensions of G are the same as str%central_dim = nrow(H) = ncol(H)
 ! Also, it does not work with tridiagonal blocks, but only with fully allocated blk matrices
 subroutine dns2blk(Gdns,str,Gblk)
   type(z_DNS), intent(in) :: Gdns
   type(Tstruct_info), intent(in) :: str
   type(z_DNS), allocatable, dimension(:,:) :: Gblk

   integer :: n, m, ii, jj, pl_start1, pl_end1, pl_start2, pl_end2, nbl, nrows

   nbl = str%num_PLs
   nrows = str%central_dim
   allocate(Gblk(nbl,nbl))
   ! call create(Gdns,nrows,nrows)

   do n = 1, nbl
     do m= 1, nbl
        pl_start1 = str%mat_PL_start(n)
        pl_end1 = str%mat_PL_end(n)
        pl_start2 = str%mat_PL_start(m)
        pl_end2 = str%mat_PL_end(m)

        if (allocated(Gblk(n,m)%val)) call destroy(Gblk(n,m))
        call create(Gblk(n,m), pl_end1-pl_start1+1, pl_end2-pl_start2+1)

        do ii = 1, pl_end1 - pl_start1 + 1
          do jj = 1, pl_end2 - pl_start2 + 1
            Gblk(n,m)%val(ii,jj) = Gdns%val(pl_start1 + ii - 1,pl_start2 + jj - 1)
          end do
        end do
     end do
   end do
 end subroutine dns2blk

 !------------------------------------------------------------
 !> Convert a CSR in a block dense matrix, only taking the
 !  diagonal and sub/over diag
 subroutine zcsr2blk3diag(Acsr,Ablk,indblk,str)
   integer :: i
   type(z_CSR) :: Acsr
   integer :: nbl
   type(z_DNS), DIMENSION(:,:) :: Ablk
   integer, DIMENSION(:) :: indblk
   character(len=*),intent(in),optional :: str

   character(1) :: si1, si2
   nbl = size(Ablk,1)

   do i=1,nbl
     write(si1,'(I1)') i
     CALL extract(Acsr,indblk(i),indblk(i+1)-1,indblk(i),indblk(i+1)-1,Ablk(i,i),str//"("//si1//","//si1//")")

   end do

   do i=2,nbl

     write(si1,'(I1)') i
     write(si2,'(I1)') i-1
     CALL extract(Acsr,indblk(i-1),indblk(i)-1,indblk(i),indblk(i+1)-1,Ablk(i-1,i),str//"("//si2//","//si1//")")
     CALL extract(Acsr,indblk(i),indblk(i+1)-1,indblk(i-1),indblk(i)-1,Ablk(i,i-1),str//"("//si1//","//si2//")")

   end do

 end subroutine zcsr2blk3diag

end module block_3_diag
