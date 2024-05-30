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

module iterative_gpu
  use ln_precision
  use ln_constants, only : pi
  use ln_allocation
  use mat_def
  use sparsekit_drv, only : zspectral
  use ln_structure, only : TStruct_Info
  use lib_param, only : MAXNCONT, Tnegf, intarray, cusolverDnHandle, cublasHandle
  use interactions, only : TInteraction, TInteractionList, TInteractionNode
  use mpi_globals, only : id, numprocs, id0
  use clock
  use cudautils
  use, intrinsic :: ieee_arithmetic

  implicit none
  private

  public :: calculate_gsmr_blocks
  public :: calculate_Gr_tridiag_blocks
  public :: calculate_Gn_tridiag_blocks
  public :: calculate_single_transmission_2_contacts
  public :: calculate_single_transmission_N_contacts
  public :: check_convergence_trid
  public :: check_convergence_vec

  interface calculate_gsmr_blocks
     module procedure  calculate_gsmr_blocks_sp
     module procedure  calculate_gsmr_blocks_dp
  end interface calculate_gsmr_blocks

  interface calculate_Gr_tridiag_blocks
     module procedure calculate_Gr_tridiag_blocks_sp
     module procedure calculate_Gr_tridiag_blocks_dp
  end interface calculate_Gr_tridiag_blocks

  interface calculate_single_transmission_2_contacts
     module procedure calculate_single_transmission_2_contacts_sp
     module procedure calculate_single_transmission_2_contacts_dp
  end interface calculate_single_transmission_2_contacts

  interface calculate_single_transmission_N_contacts
     module procedure calculate_single_transmission_N_contacts_sp
     module procedure calculate_single_transmission_N_contacts_dp
  end interface calculate_single_transmission_N_contacts

  interface get_tun_mask
     module procedure get_tun_mask_sp
     module procedure get_tun_mask_dp
  end interface get_tun_mask

contains

  subroutine calculate_gsmr_blocks_sp(negf,ESH,sbl,ebl,gsmr,keepall)

    !In/Out
    type(c_DNS), dimension(:), intent(inout) :: gsmr
    type(Tnegf), intent(in) :: negf
    type(c_DNS), dimension(:,:), intent(inout) :: ESH
    integer, intent(in) :: sbl,ebl                 ! start block, end block
    logical, intent(in), optional :: keepall

    !Work
    type(CublasHandle) :: hh
    type(CusolverDnHandle) :: hhsol
    complex(sp), parameter :: one = (1.0_sp, 0.0_sp)
    complex(sp), parameter :: mone = (-1.0_sp, 0.0_sp)
    complex(sp), parameter :: zero = (0.0_sp, 0.0_sp)
    type(c_DNS) :: work1, work2
    integer :: nrow, ncol
    integer :: i, nbl, istat
    real(sp) :: sum1
    logical :: keep

    if (sbl.lt.ebl) return

    nbl = size(ESH,1)
    if (nbl.eq.1) return

    keep = .true.
    if (present(keepall)) then
       keep = keepall
    end if

    nrow=ESH(sbl,sbl)%nrow

    hh = negf%hcublas
    hhsol = negf%hcusolver

    call createAll(gsmr(sbl),nrow,nrow)
    call copyToGPU(ESH(sbl,sbl))
    call inverse_gpu(hh, hhsol, ESH(sbl,sbl), gsmr(sbl), istat)
    call deleteGPU(ESH(sbl,sbl))

    do i=sbl-1,ebl,-1

       call createAll(work1, ESH(i,i+1)%nrow, gsmr(i+1)%ncol)
       call copyToGPU(ESH(i,i+1))
       call matmul_gpu(hh, one, ESH(i,i+1), gsmr(i+1), zero, work1)
       call deleteGPU(ESH(i,i+1))

       if (.not.keep) then
          call destroyAll(gsmr(i+1))
       end if

       call createAll(work2, work1%nrow, ESH(i+1,i)%ncol)
       call copyToGPU(ESH(i+1,i))
       call matmul_gpu(hh, one, work1, ESH(i+1,i), zero, work2)
       call deleteGPU(ESH(i+1,i))

       call destroyAll(work1)

       call createAll(work1,  ESH(i,i)%nrow, ESH(i,i)%ncol)
       call copyToGPU(ESH(i,i))
       call matsum_gpu(hh, one, ESH(i,i), mone, work2, work1)
       call deleteGPU(ESH(i,i))

       call destroyAll(work2)

       call createAll(gsmr(i), work1%nrow, work1%ncol)
       call inverse_gpu(hh, hhsol, work1, gsmr(i), istat)
       call destroyAll(work1)

    end do

  end subroutine calculate_gsmr_blocks_sp

  subroutine calculate_gsmr_blocks_dp(negf,ESH,sbl,ebl,gsmr,keepall)

    !In/Out
    type(z_DNS), dimension(:), intent(inout) :: gsmr
    type(Tnegf), intent(in) :: negf
    type(z_DNS), dimension(:,:), intent(inout) :: ESH
    integer, intent(in) :: sbl,ebl                 ! start block, end block
    logical, intent(in), optional :: keepall

    !Work
    type(CublasHandle) :: hh
    type(CusolverDnHandle) :: hhsol
    complex(dp), parameter :: one = (1.0_dp, 0.0_dp)
    complex(dp), parameter :: mone = (-1.0_dp, 0.0_dp)
    complex(dp), parameter :: zero = (0.0_dp, 0.0_dp)
    type(z_DNS) :: work1, work2
    integer :: nrow, ncol
    integer :: i, nbl, istat
    real(dp) :: sum1
    logical :: keep

    if (sbl.lt.ebl) return

    nbl = size(ESH,1)
    if (nbl.eq.1) return

    keep = .true.
    if (present(keepall)) then
       keep = keepall
    end if


    nrow=ESH(sbl,sbl)%nrow

    hh = negf%hcublas
    hhsol = negf%hcusolver

    call createAll(gsmr(sbl),nrow,nrow)
    call copyToGPU(ESH(sbl,sbl))
    call inverse_gpu(hh, hhsol, ESH(sbl,sbl), gsmr(sbl), istat)
    call deleteGPU(ESH(sbl,sbl))

    do i=sbl-1,ebl,-1

       call createAll(work1, ESH(i,i+1)%nrow, gsmr(i+1)%ncol)
       call copyToGPU(ESH(i,i+1))
       call matmul_gpu(hh, one, ESH(i,i+1), gsmr(i+1), zero, work1)
       call deleteGPU(ESH(i,i+1))

       if (.not.keep) then
          call destroyAll(gsmr(i+1))
       end if

       call createAll(work2, work1%nrow, ESH(i+1,i)%ncol)
       call copyToGPU(ESH(i+1,i))
       call matmul_gpu(hh, one, work1, ESH(i+1,i), zero, work2)
       call deleteGPU(ESH(i+1,i))

       call destroyAll(work1)

       call createAll(work1,  ESH(i,i)%nrow, ESH(i,i)%ncol)
       call copyToGPU(ESH(i,i))
       call matsum_gpu(hh, one, ESH(i,i), mone, work2, work1)
       call deleteGPU(ESH(i,i))

       call destroyAll(work2)

       call createAll(gsmr(i), work1%nrow, work1%ncol)
       call inverse_gpu(hh, hhsol, work1, gsmr(i), istat)
       call destroyAll(work1)

    end do

  end subroutine calculate_gsmr_blocks_dp

  subroutine calculate_Gr_tridiag_blocks_sp(negf,ESH,gsmr,Gr,sbl,ebl)
    !In/Out
    type(c_DNS), dimension(:,:), intent(inout) :: Gr
    type(c_DNS), dimension(:), intent(in) :: gsmr
    type(Tnegf), intent(in) :: negf
    type(c_DNS), dimension(:,:), intent(inout) :: ESH
    integer, intent(in) :: sbl
    integer, intent(in), optional :: ebl

    !Work
    type(CublasHandle) :: hh
    type(CusolverDnHandle) ::hhsol
    integer :: i,nrow,nbl
    type(c_DNS), target :: work1, work2, work3
    complex(sp), parameter :: one = (1.0_sp, 0.0_sp)
    complex(sp), parameter :: mone = (-1.0_sp, 0.0_sp)
    complex(sp), parameter :: zero = (0.0_sp, 0.0_sp)
    integer :: istat

    nbl = size(ESH,1)
    hh = negf%hcublas
    hhsol = negf%hcusolver

    if (sbl.gt.nbl) return
    if (sbl.lt.1) return

    if (.not.present(ebl)) then
       if (nbl.eq.1) then
          call createAll(Gr(sbl,sbl), ESH(sbl,sbl)%nrow, ESH(sbl,sbl)%ncol)
          call copyToGPU(ESH(sbl,sbl))
          call inverse_gpu(hh, hhsol, ESH(sbl,sbl), Gr(sbl,sbl), istat)
          call deleteGPU(ESH(sbl,sbl))
       else
          call createAll(work1, ESH(sbl,sbl)%nrow, ESH(sbl,sbl)%ncol)
          call copyToGPU(ESH(sbl,sbl))
          call copy_mat_gpu(hh, ESH(sbl,sbl), work1)
          call deleteGPU(ESH(sbl,sbl))

          if (sbl+1.le.nbl) then
             call createAll(work2, ESH(sbl,sbl+1)%nrow, gsmr(sbl+1)%ncol)
             call copyToGPU(ESH(sbl,sbl+1))
             call matmul_gpu(hh, one, ESH(sbl,sbl+1), gsmr(sbl+1), zero, work2)
             call deleteGPU(ESH(sbl,sbl+1))

             call createAll(work3, work2%nrow, ESH(sbl+1,sbl)%ncol)
             call copyToGPU(ESH(sbl+1,sbl))
             call matmul_gpu(hh, one, work2, ESH(sbl+1,sbl), zero, work3)
             call deleteGPU(ESH(sbl+1,sbl))

             call copyToGPU(ESH(sbl,sbl))
             call matsum_gpu(hh, one, ESH(sbl,sbl), mone, work3, work1)
             call deleteGPU(ESH(sbl,sbl))

             call destroyAll(work2)
             call destroyAll(work3)
          end if
          if (sbl-1.ge.1) then
             stop "Error: Gr_tridiag requires gsml"
          end if

          call createAll(Gr(sbl,sbl), work1%nrow, work1%ncol)
          call inverse_gpu(hh, hhsol, work1, Gr(sbl,sbl), istat)
          call destroyAll(work1)
       endif
       return
    endif
    !***
    !Diagonal, Subdiagonal and Superdiagonal blocks
    !***
    if ((ebl.ge.sbl).and.(ebl.gt.1).and.(sbl.gt.1)) THEN
       do i=sbl,ebl,1
          call createAll(work1, gsmr(i)%nrow, ESH(i,i-1)%ncol)
          call createAll(Gr(i,i-1), work1%nrow, Gr(i-1,i-1)%ncol)
          call copyToGPU(ESH(i,i-1))
          call matmul_gpu(hh, one, gsmr(i), ESH(i,i-1), zero, work1)
          call deleteGPU(ESH(i,i-1))
          call matmul_gpu(hh, mone, work1, Gr(i-1,i-1), zero, Gr(i,i-1))

          call destroyAll(work1)

          call createAll(work2, ESH(i-1,i)%nrow, gsmr(i)%ncol)
          call createAll(Gr(i-1,i), Gr(i-1,i-1)%nrow, work2%ncol)
          call copyToGPU(ESH(i-1,i))
          call matmul_gpu(hh, one, ESH(i-1,i), gsmr(i), zero, work2)
          call deleteGPU(ESH(i-1,i))
          call matmul_gpu(hh, mone, Gr(i-1,i-1), work2, zero, Gr(i-1,i))

          call createAll(work1, Gr(i,i-1)%nrow, work2%ncol)
          call matmul_gpu(hh, mone, Gr(i,i-1), work2, zero, work1)

          call destroyAll(work2)
          call createAll(Gr(i,i), gsmr(i)%nrow, gsmr(i)%ncol)
          call matsum_gpu(hh, one, gsmr(i), one, work1, Gr(i,i))

          call destroyAll(work1)
       end do
    else
       stop "Error: Gr_tridiag requires gsml"
    endif

  end subroutine calculate_Gr_tridiag_blocks_sp

  subroutine calculate_Gr_tridiag_blocks_dp(negf,ESH,gsmr,Gr,sbl,ebl)
    !In/Out
    type(z_DNS), dimension(:,:), intent(inout) :: Gr
    type(z_DNS), dimension(:), intent(in) :: gsmr
    type(Tnegf), intent(in) :: negf
    type(z_DNS), dimension(:,:), intent(inout) :: ESH
    integer, intent(in) :: sbl
    integer, intent(in), optional :: ebl

    !Work
    type(CublasHandle) :: hh
    type(CusolverDnHandle) :: hhsol
    integer :: i, nrow, nbl, istat
    type(z_DNS), target :: work1, work2, work3
    complex(dp), parameter :: one = (1.0_dp, 0.0_dp)
    complex(dp), parameter :: mone = (-1.0_dp, 0.0_dp)
    complex(dp), parameter :: zero = (0.0_dp, 0.0_dp)
    real(dp) :: summ

    nbl = size(ESH,1)
    hh = negf%hcublas
    hhsol = negf%hcusolver

    if (sbl.gt.nbl) return
    if (sbl.lt.1) return

    if (.not.present(ebl)) then
       if (nbl.eq.1) then
          call createAll(Gr(sbl,sbl), ESH(sbl,sbl)%nrow, ESH(sbl,sbl)%ncol)
          call copyToGPU(ESH(sbl,sbl))
          call inverse_gpu(hh, hhsol, ESH(sbl,sbl), Gr(sbl,sbl), istat)
          call deleteGPU(ESH(sbl,sbl))
       else
          call createAll(work1, ESH(sbl,sbl)%nrow, ESH(sbl,sbl)%ncol)
          call copyToGPU(ESH(sbl,sbl))
          call copy_mat_gpu(hh, ESH(sbl,sbl), work1)
          call deleteGPU(ESH(sbl,sbl))

          if (sbl+1.le.nbl) then
             call createAll(work2, ESH(sbl,sbl+1)%nrow, gsmr(sbl+1)%ncol)
             call copyToGPU(ESH(sbl,sbl+1))
             call matmul_gpu(hh, one, ESH(sbl,sbl+1), gsmr(sbl+1), zero, work2)
             call deleteGPU(ESH(sbl,sbl+1))

             call createAll(work3, work2%nrow, ESH(sbl+1,sbl)%ncol)
             call copyToGPU(ESH(sbl+1,sbl))
             call matmul_gpu(hh, one, work2, ESH(sbl+1,sbl), zero, work3)
             call deleteGPU(ESH(sbl+1,sbl))

             call copyToGPU(ESH(sbl,sbl))
             call matsum_gpu(hh, one, ESH(sbl,sbl), mone, work3, work1)
             call deleteGPU(ESH(sbl,sbl))

             call destroyAll(work2)
             call destroyAll(work3)
          end if
          if (sbl-1.ge.1) then
             stop "Error: Gr_tridiag requires gsml"
          end if

          call createAll(Gr(sbl,sbl), work1%nrow, work1%ncol)
          call inverse_gpu(hh, hhsol, work1, Gr(sbl,sbl), istat)
          call destroyAll(work1)
       endif
       return
    endif
    !***
    !Diagonal, Subdiagonal and Superdiagonal blocks
    !***
    if ((ebl.ge.sbl).and.(ebl.gt.1).and.(sbl.gt.1)) THEN
       do i=sbl,ebl,1
          call createAll(work1, gsmr(i)%nrow, ESH(i,i-1)%ncol)
          call createAll(Gr(i,i-1), work1%nrow, Gr(i-1,i-1)%ncol)
          call copyToGPU(ESH(i,i-1))
          call matmul_gpu(hh, one, gsmr(i), ESH(i,i-1), zero, work1)
          call deleteGPU(ESH(i,i-1))
          call matmul_gpu(hh, mone, work1, Gr(i-1,i-1), zero, Gr(i,i-1))

          call destroyAll(work1)

          call createAll(work2, ESH(i-1,i)%nrow, gsmr(i)%ncol)
          call createAll(Gr(i-1,i), Gr(i-1,i-1)%nrow, work2%ncol)
          call copyToGPU(ESH(i-1,i))
          call matmul_gpu(hh, one, ESH(i-1,i), gsmr(i), zero, work2)
          call deleteGPU(ESH(i-1,i))
          call matmul_gpu(hh, mone, Gr(i-1,i-1), work2, zero, Gr(i-1,i))

          call createAll(work1, Gr(i,i-1)%nrow, work2%ncol)
          call matmul_gpu(hh, mone, Gr(i,i-1), work2, zero, work1)

          call destroyAll(work2)
          call createAll(Gr(i,i), gsmr(i)%nrow, gsmr(i)%ncol)
          call matsum_gpu(hh, one, gsmr(i), one, work1, Gr(i,i))
          call destroyAll(work1)
       end do
    else
       stop "Error: Gr_tridiag requires gsml"
    endif

  end subroutine calculate_Gr_tridiag_blocks_dp

  subroutine check_convergence_trid(negf,T,nome,gpu)
    type(Tnegf), intent(in) :: negf
    type(z_DNS), dimension(:,:), intent(in) :: T
    character(3), intent(in) :: nome
    logical, intent(in) :: gpu

    integer :: nbl, i
    type(CublasHandle) :: hh
    real(dp) :: summ
    character(1) :: ci, cim1

    nbl = size(T,1)
    hh = negf%hcublas
    if (gpu) then
       !write(*,*) '~-~-~-~-',nome,' check convergence GPU: ~-~-~-~-'
       call checksum(hh, T(1,1), nome//'(1,1)')

       do i= 2,nbl-1
          write(ci,'(i1)') i
          write(cim1,'(i1)') i-1
          call checksum(hh, T(i,i), nome//'('//ci//','//ci//')')
          call checksum(hh, T(i-1,i), nome//'('//cim1//','//ci//')')
          call checksum(hh, T(i,i-1), nome//'('//ci//','//cim1//')')
       end do
       !write(*,*) '~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-'
    else
       !write(*,*) '~-~-~-~-',nome,' check convergence CPU: ~-~-~-~-'
       summ = sum(ABS(T(1,1)%val))
       if (ieee_is_nan(summ)) write(*,*) 'CPU:   ',nome,'(',1,1,')=', summ

       do i= 2,nbl-1
          summ = sum(ABS(T(i,i)%val))
          if (ieee_is_nan(summ)) write(*,*) 'CPU:   ',nome,'(',i,i,')=', summ
          summ = sum(ABS(T(i-1,i)%val))
          if (ieee_is_nan(summ)) write(*,*) 'CPU:   ',nome,'(',i-1,i,')=', summ
          summ = sum(ABS(T(i,i-1)%val))
          if (ieee_is_nan(summ)) write(*,*) 'CPU:   ',nome,'(',i,i-1,')=', summ
       end do
       !write(*,*) '~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-'
    endif
  end subroutine check_convergence_trid

  subroutine check_convergence_vec(negf,T,nome,gpu)
    type(Tnegf), intent(in) :: negf
    type(z_DNS), dimension(:), intent(in) :: T
    character(4), intent(in) :: nome
    logical, intent(in) :: gpu

    integer :: nbl, i
    type(CublasHandle) :: hh
    real(dp) :: summ
    character(1) :: ci

    nbl = size(T)
    hh = negf%hcublas
    if (gpu) then
       !write(*,*) '~-~-~-~-',nome,' check convergence GPU: ~-~-~-~-'

       do i= 2,nbl
          write(ci,'(i1)') i
          call checksum(hh, T(i), nome//'('//ci//')')
       end do
       !write(*,*) '~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-'
    else
       !write(*,*) '~-~-~-~-',nome,' check convergence CPU: ~-~-~-~-'

       do i= 2,nbl
          summ = sum(ABS(T(i)%val))
          if (ieee_is_nan(summ)) write(*,*) 'CPU:   ',nome,'(',i,')=', summ
       end do
       !write(*,*) '~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-'
    endif

  end subroutine check_convergence_vec

  subroutine calculate_Gn_tridiag_blocks(negf,ESH,SelfEneR,frm,ref,struct,gsmr,Gr,Gn)
    type(TNegf), intent(in) :: negf
    type(z_DNS), dimension(:,:), intent(inout) :: ESH, Gr
    type(z_DNS), dimension(:), intent(in) :: SelfEneR, gsmr
    real(dp), dimension(:), intent(in) :: frm
    integer, intent(in) :: ref
    type(Tstruct_info), intent(in) :: struct
    type(z_DNS), dimension(:,:), intent(inout) :: Gn

    !Work
    type(CublasHandle) :: hh
    complex(dp) :: one = (1.0_dp, 0.0_dp)
    complex(dp) :: zero = (0.0_dp, 0.0_dp)
    complex(dp) :: minusone = (-1.0_dp, 0.0_dp)
    type(z_DNS), dimension(:,:), allocatable :: Sigma_n
    type(z_DNS) :: work1, Gam
    complex(dp) :: frmdiff
    integer :: i, j
    integer :: nbl, ncont, cb

    ncont = struct%num_conts
    nbl = struct%num_PLs
    hh = negf%hcublas

    !build Sigma_n from SelfEneR
    call allocate_blk_dns(Sigma_n, nbl)
    call init_tridiag_blk(Sigma_n, ESH)

    ! Add interaction self-energies
    call add_sigma_n(negf, Sigma_n)

    call copy_trid_toGPU(Sigma_n)
    ! Add contact self-energies
    do j=1,ncont
      frmdiff = cmplx(frm(j) - frm(ref),0.0_dp,dp)
      if (j.NE.ref .AND. ABS(frmdiff).GT.EPS) THEN
        cb=struct%cblk(j) ! block corresponding to contact j
        call zspectral(SelfEneR(j),SelfEneR(j),0,Gam)
        call copyToGPU(Gam)
        call matsum_gpu(hh, one, Sigma_n(cb,cb), frmdiff, Gam, Sigma_n(cb,cb))
        call destroyAll(Gam)
      endif
    end do

    call calculate_sigma_n()

    call createAll(work1,Sigma_n(1,1)%nrow, Gr(1,1)%nrow)
    call matmul_gpu(hh, one, Sigma_n(1,1), Gr(1,1), zero, work1, 'dag_2nd')
    call matmul_gpu(hh, one, Gr(1,1), work1, zero, Gn(1,1))
    call destroyAll(work1)

    if (nbl .eq. 1) then
      call destroy_tridiag_blk(Sigma_n)
      call deallocate_blk_dns(Sigma_n)
      return
    endif

    !Explicit formulae:
    !Gn(i+1,i) = gsmr(i+1)*[Sigma(i+1,i)Ga(i,i) + Sigma(i+1,i+1)Ga(i+1,i) - Tr(i+1,i)Gn(i,i)]
    !Gn(i,i+1) = [Gr(i,i)Sigma(i,i+1) + Gr(i,i+1)Sigma(i+1,i+1) - Gn(i,i)Ta(i,i+1)] * gsma(i+1)
    !Use Hermitian property of Gn:
    !Gn(i,i+1) = Gn(i+1,i)^dag
    !Gn(i+1,i+1) = gsmr(i+1) * [Sigma(i+1,i)Ga(i,i+1) + Sigma(i+1,i+1)Ga(i+1,i+1) - Tr(i+1,i)Gn(i,i+1)]
    !Implementation exploits cumulative sum of prealloc_mult, C = C + A*B

    do i = 1, nbl-1

        call createAll(work1, Sigma_n(i+1,i)%nrow, Gr(i,i)%nrow)
        call matmul_gpu(hh, one, Sigma_n(i+1,i), Gr(i,i), zero, work1, 'dag_2nd')
        call matmul_gpu(hh, one, Sigma_n(i+1,i+1), Gr(i,i+1), one, work1, 'dag_2nd')

        call copyToGPU(ESH(i+1,i))
        call matmul_gpu(hh, minusOne, ESH(i+1,i), Gn(i,i), one, work1)
        call deleteGPU(ESH(i+1,i))

        call matmul_gpu(hh, one, gsmr(i+1), work1, zero, Gn(i+1,i))
        call destroyAll(work1)

        call dagger_gpu(hh, Gn(i+1,i),Gn(i,i+1))

        call createAll(work1, Sigma_n(i+1,i)%nrow, Gr(i+1,i)%nrow)
        call matmul_gpu(hh, one,Sigma_n(i+1,i), Gr(i+1,i), zero, work1, 'dag_2nd')

        call matmul_gpu(hh, one,Sigma_n(i+1,i+1), Gr(i+1,i+1), one, work1, 'dag_2nd')

        call copyToGPU(ESH(i+1,i))
        call matmul_gpu(hh, minusOne, ESH(i+1,i), Gn(i,i+1), one, work1)
        call deleteGPU(ESH(i+1,i))


        call matmul_gpu(hh, one,gsmr(i+1), work1, zero, Gn(i+1,i+1))
        call destroyAll(work1)

    end do

    call destroy_tridiag_blk(Sigma_n)
    call deallocate_blk_dns(Sigma_n)

    contains
    ! Recursive calculation of Sigma_n:
    ! gns(i+1) = gsmr(i+1) Sigma(i+1,i+1) gsmr(i+1)^dag
    ! Sigma(i,i) = Sigma(i,i) + Tr(i,i+1) gns(i+1) Ta(i+1,i)
    !                         - Tr(i,i+1) gsmr(i+1) Sigma(i+1,i)
    !                         - Sigma(i,i+1) gsmr^dag(i+1) Ta(i+1,i)]
    !
    subroutine calculate_sigma_n()
      !Work
      type(z_DNS) :: work, gns

      ! if nbl = 1 => Sigma_n(1,1) is ready
      if (nbl.eq.1) return
      !g^n(nbl) = gsmr(nbl) Sigma(nbl,nbl) gsma(nbl)
      call createAll(work1, gsmr(nbl)%nrow, Sigma_n(nbl,nbl)%ncol)
      call createAll(gns, work1%nrow, gsmr(nbl)%nrow)
      call matmul_gpu(hh, one, gsmr(nbl), Sigma_n(nbl,nbl), zero, work1)
      call matmul_gpu(hh, one, work1, gsmr(nbl), zero, gns, 'dag_2nd')
      call destroyAll(work1)

      do i = nbl-1, 1, -1
        !work1 = Tr(i,i+1) gns(i+1) Ta(i+1,i)
        ! Tr(i,i+1) = ESH(i,i+1);  Ta(i+1,i) = ESH(i,i+1)^dag
        call createAll(work, ESH(i,i+1)%nrow, gns%ncol)


        call copyToGPU(ESH(i,i+1))
        call matmul_gpu(hh, one, ESH(i,i+1), gns, zero, work)
        call matmul_gpu(hh, one, work, ESH(i,i+1), one, Sigma_n(i,i), 'dag_2nd')
        call deleteGPU(ESH(i,i+1))

        call destroyAll(work)
        call destroyAll(gns)

        !work2 = Sigma(i,i+1) gsmr^dag(i+1) Ta(i+1,i)
        call createAll(work, Sigma_n(i,i+1)%nrow, gsmr(i+1)%ncol)

        call matmul_gpu(hh, one, Sigma_n(i,i+1), gsmr(i+1), zero, work, 'dag_2nd')
        call copyToGPU(ESH(i,i+1))
        call matmul_gpu(hh, minusOne, work, ESH(i,i+1), one, Sigma_n(i,i), 'dag_2nd')
        call deleteGPU(ESH(i,i+1))

        call destroyAll(work)

        !work3 = ESH(i,i+1) gsmr(i+1) Sigma(i+1,i)
        call createAll(work, ESH(i,i+1)%nrow, gsmr(i+1)%ncol)

        call copyToGPU(ESH(i,i+1))
        call matmul_gpu(hh, one, ESH(i,i+1), gsmr(i+1), zero, work)
        call deleteGPU(ESH(i,i+1))
        call matmul_gpu(hh, minusOne, work, Sigma_n(i+1,i), one, Sigma_n(i,i))

        call destroyAll(work)

        if (i > 1) then
          !gns(i) = gsmr(i) * Sigma_n(i,i) * gsmr^dag(i)
          call createAll(work, gsmr(i)%nrow, Sigma_n(i,i)%ncol)
          call createAll(gns, work%nrow, gsmr(i)%nrow)

          call matmul_gpu(hh, one, gsmr(i), Sigma_n(i,i), zero, work)
          call matmul_gpu(hh, one, work, gsmr(i), zero, gns, 'dag_2nd')

          call destroyAll(work)
        end if

      end do

    end subroutine calculate_sigma_n

  end subroutine calculate_Gn_tridiag_blocks

  !--------------------------------------------------------------------------
  ! Add all Self-enrgies to Sigma_n
  subroutine add_sigma_n(negf, sigma_n)
    class(TNegf) :: negf
    type(z_DNS) :: sigma_n(:,:)

    type(TInteractionNode), pointer :: it
    it => negf%interactList%first
    
    do while (associated(it))
      call it%inter%add_sigma_n(sigma_n, negf%iEloc, negf%iKloc, negf%spin)
      it => it%next
    end do
  end subroutine add_sigma_n



  subroutine calculate_single_transmission_2_contacts_sp(negf,ni,nf,ESH,SelfEneR,cblk,tun_proj,Gr,TUN)
    type(TNegf), intent(in) :: negf
    integer, intent(in) :: ni,nf
    type(c_DNS), intent(in) :: SelfEneR(MAXNCONT)
    type(c_DNS), intent(in) :: ESH(:,:)
    integer, intent(in) :: cblk(:)
    type(intarray), intent(in) :: tun_proj
    type(c_DNS), dimension(:,:), intent(in) :: Gr
    real(sp), intent(out) :: TUN

    !Work variables
    type(CublasHandle) :: hh
    Integer :: ct1, bl1
    logical, dimension(:), allocatable :: tun_mask
    Type(c_DNS) :: work1, work2, GAM1_dns, TRS, AA
    complex(sp), parameter :: j = (0.0_sp,1.0_sp)  ! CMPX unity
    complex(sp), parameter :: mj = (0.0_sp,-1.0_sp)
    complex(sp), parameter :: one = (1.0_sp, 0.0_sp)
    complex(sp), parameter :: mone = (-1.0_sp, 0.0_sp)
    complex(sp), parameter :: zero = (0.0_sp, 0.0_sp)

    if (size(cblk).gt.2) then
       write(*,*) "ERROR: calculate_single_transmission_2_contacts is valid only for 2 contacts"
       TUN = 0.0_sp
       return
    endif

    hh = negf%hcublas

    !Arrange contacts in way that order between first and second is always the
    !same (always ct1 < ct2)

    if (cblk(ni).lt.cblk(nf)) then
       ct1=ni;
    else
       ct1=nf;
    endif

    bl1=cblk(ct1);

    ! Computes the Gamma matrices
    call createAll(GAM1_dns, SelfEneR(ct1)%nrow, SelfEneR(ct1)%ncol)
    call spectral_gpu(hh, SelfEneR(ct1), GAM1_dns)

    ! Work to compute transmission matrix (Gamma G Gamma G)
    call createAll(work1, GAM1_dns%nrow, Gr(bl1,bl1)%ncol)
    call createAll(work2, work1%nrow, GAM1_dns%ncol)
    call matmul_gpu(hh, one, GAM1_dns, Gr(bl1,bl1), zero, work1)
    call matmul_gpu(hh, one, work1, GAM1_dns, zero, work2)

    call destroyAll(work1)

    call createAll(work1, work2%nrow, Gr(bl1,bl1)%nrow)
    call matmul_gpu(hh, one, work2, Gr(bl1,bl1), zero, work1,'dag_2nd')
    call destroyAll(work2)

    call createAll(AA, Gr(bl1,bl1)%ncol, Gr(bl1,bl1)%nrow)
    call matsum_gpu(hh, j, Gr(bl1,bl1), mj, Gr(bl1,bl1), AA, 'dag_2nd')

    call createAll(work2, GAM1_dns%nrow, AA%ncol)
    call matmul_gpu(hh, one, GAM1_dns, AA, zero, work2)
    call destroyAll(GAM1_dns)
    call destroyAll(AA)

    call createAll(TRS, work1%nrow, work1%ncol)
    call matsum_gpu(hh, one, work2, mone, work1, TRS)

    !call get_tun_mask(ESH, bl1, tun_proj, tun_mask)
    !call trace_gpu(hh, TRS, TUN, tun_mask)
    !call log_deallocate(tun_mask)
    call trace_gpu(hh, TRS, TUN)

    call destroyAll(TRS)
    call destroyAll(work1)
    call destroyAll(work2)

  end subroutine calculate_single_transmission_2_contacts_sp

  subroutine calculate_single_transmission_2_contacts_dp(negf,ni,nf,ESH,SelfEneR,cblk,tun_proj,Gr,TUN)
    type(TNegf), intent(in) :: negf
    integer, intent(in) :: ni,nf
    type(z_DNS), intent(in) :: SelfEneR(MAXNCONT)
    type(z_DNS), intent(in) :: ESH(:,:)
    integer, intent(in) :: cblk(:)
    type(intarray), intent(in) :: tun_proj
    type(z_DNS), dimension(:,:), intent(in) :: Gr
    real(dp), intent(out) :: TUN

    !Work variables
    type(CublasHandle) :: hh
    Integer :: ct1, bl1
    logical, dimension(:), allocatable :: tun_mask
    Type(z_DNS) :: work1, work2, GAM1_dns, TRS, AA
    complex(dp), parameter ::    j = (0.0_dp,1.0_dp)  ! CMPX unity
    complex(dp), parameter ::    mj = (0.0_dp,-1.0_dp)
    complex(dp), parameter :: one = (1.0_dp, 0.0_dp)
    complex(dp), parameter :: mone = (-1.0_dp, 0.0_dp)
    complex(dp), parameter :: zero = (0.0_dp, 0.0_dp)
    real(dp) :: summ

    if (size(cblk).gt.2) then
       write(*,*) "ERROR: calculate_single_transmission_2_contacts is valid only for 2 contacts"
       TUN = 0.0_dp
       return
    endif

    hh = negf%hcublas

    !Arrange contacts in way that order between first and second is always the
    !same (always ct1 < ct2)

    if (cblk(ni).lt.cblk(nf)) then
       ct1=ni;
    else
       ct1=nf;
    endif

    bl1=cblk(ct1);

    ! Computes the Gamma matrices
    call createAll(GAM1_dns, SelfEneR(ct1)%nrow, SelfEneR(ct1)%ncol)
    call spectral_gpu(hh, SelfEneR(ct1), GAM1_dns)

    ! Work to compute transmission matrix (Gamma G Gamma G)
    call createAll(work1, GAM1_dns%nrow, Gr(bl1,bl1)%ncol)
    call createAll(work2, work1%nrow, GAM1_dns%ncol)
    call matmul_gpu(hh, one, GAM1_dns, Gr(bl1,bl1), zero, work1)
    call matmul_gpu(hh, one, work1, GAM1_dns, zero, work2)

    call destroyAll(work1)

    call createAll(work1, work2%nrow, Gr(bl1,bl1)%nrow)
    call matmul_gpu(hh, one, work2, Gr(bl1,bl1), zero, work1,'dag_2nd')
    call destroyAll(work2)

    call createAll(AA, Gr(bl1,bl1)%ncol, Gr(bl1,bl1)%nrow)
    call matsum_gpu(hh, j, Gr(bl1,bl1), mj, Gr(bl1,bl1), AA, 'dag_2nd')

    call createAll(work2, GAM1_dns%nrow, AA%ncol)
    call matmul_gpu(hh, one, GAM1_dns, AA, zero, work2)
    call destroyAll(GAM1_dns)
    call destroyAll(AA)

    call createAll(TRS, work1%nrow, work1%ncol)
    call matsum_gpu(hh, one, work2, mone, work1, TRS)

    !call get_tun_mask(ESH, bl1, tun_proj, tun_mask)
    !call trace_gpu(hh, TRS, TUN, tun_mask)
    call trace_gpu(hh, TRS, TUN)
    !call log_deallocate(tun_mask)

    call destroyAll(TRS)
    call destroyAll(work1)
    call destroyAll(work2)

  end subroutine calculate_single_transmission_2_contacts_dp

  subroutine calculate_single_transmission_N_contacts_sp(negf,ni,nf,ESH,SelfEneR,cblk,tun_proj,gsmr,Gr,TUN)
    type(TNegf), intent(in) :: negf
    integer, intent(in) :: ni,nf
    type(c_DNS), intent(in) :: SelfEneR(MAXNCONT)
    type(c_DNS), intent(inout) :: ESH(:,:)
    type(c_DNS), dimension(:),intent(in) :: gsmr
    type(c_DNS), dimension(:,:),intent(inout) :: Gr
    integer, intent(in) :: cblk(:)
    type(intarray), intent(in) :: tun_proj
    real(sp), intent(out) :: TUN

    !Work variables
    type(CublasHandle) :: hh
    Integer :: ct1, ct2, bl1, bl2, i, nbl
    logical, dimension(:), allocatable :: tun_mask
    Type(c_DNS) :: work1, work2, GAM1_dns, GAM2_dns, TRS
    Real(sp) :: max
    complex(sp), parameter :: one = (1.0_sp, 0.0_sp)
    complex(sp), parameter :: mone = (-1.0_sp, 0.0_sp)
    complex(sp), parameter :: zero = (0.0_sp, 0.0_sp)

    !Arrange contacts in way that order between first and second is always the
    !same (always ct1 < ct2)

    if (cblk(ni).lt.cblk(nf)) then
       ct1=ni;ct2=nf;
    else
       ct1=nf;ct2=ni;
    endif

    bl1=cblk(ct1); bl2=cblk(ct2);
    nbl = size(cblk)
    hh = negf%hcublas
    ! in this way nt1 < nt2 by construction
    if ( nbl.gt.1 .and. (bl2-bl1).gt.1) then

       ! Compute column-blocks of Gr(i,bl1) up to i=bl2
       ! Gr(i,bl1) = -gr(i) T(i,i-1) Gr(i-1,bl1)
       do i = bl1+1, bl2
          !Checks whether previous block is non null.
          !If so next block is also null => TUN = 0
          call copyFromGPU(Gr(i-1,bl1))
          max=maxval(abs(Gr(i-1,bl1)%val))

          if (max.lt.EPS) then
             TUN = EPS*EPS !for log plots
             !Destroy also the block adjecent to diagonal since
             !this is not deallocated anymore in calling subroutine
             if (i.gt.(bl1+1)) call destroyAll(Gr(i-1,bl1))
             return
          endif

          !Checks whether block has been created, if not do it
          if (.not.allocated(Gr(i,bl1)%val)) then

             call createAll(work1, gsmr(i)%nrow, ESH(i,i-1)%ncol)
             call createAll(Gr(i,bl1), work1%nrow, Gr(i-1,bl1)%ncol)
             call copyToGPU(ESH(i,i-1))
             call matmul_gpu(hh, mone, gsmr(i), ESH(i,i-1), zero, work1)
             call deleteGPU(ESH(i,i-1))

             call matmul_gpu(hh, one, work1, Gr(i-1,bl1) ,zero, Gr(i,bl1))
             call destroyAll(work1)

          endif

          ! avoid destroying blocks closer to diagonal
          if (i.gt.(bl1+2)) call destroyAll(Gr(i-1,bl1))
       end do

    endif
    ! Computes the Gamma matrices
    call createAll(GAM1_dns, SelfEneR(ct1)%nrow, SelfEneR(ct1)%ncol)
    call createAll(GAM2_dns, SelfEneR(ct2)%nrow, SelfEneR(ct2)%ncol)
    call spectral_gpu(hh, SelfEneR(ct1), GAM1_dns)
    call spectral_gpu(hh, SelfEneR(ct2), GAM2_dns)

    ! Work to compute transmission matrix (Gamma2 Gr Gamma1 Ga)
    call createAll(work1, Gr(bl2,bl1)%nrow, GAM1_dns%ncol)
    call matmul_gpu(hh, one, Gr(bl2,bl1), GAM1_dns, zero, work1)

    call createAll(work2, GAM2_dns%nrow, work1%ncol)
    call matmul_gpu(hh, one, GAM2_dns, work1, zero, work2)

    call destroyAll(work1)
    call destroyAll(GAM2_dns)
    call destroyAll(GAM1_dns)

    call createAll(TRS, work2%nrow, Gr(bl2,bl1)%nrow)

    call matmul_gpu(hh, one, work2, Gr(bl2,bl1), zero, TRS, 'dag_2nd')
    call destroyAll(work2)
    if (bl2.gt.bl1+1) call destroyAll(Gr(bl2,bl1))

    !call get_tun_mask(ESH, bl2, tun_proj, tun_mask)
    !call trace_gpu(hh, TRS, TUN, tun_mask)
    call trace_gpu(hh, TRS, TUN)
    !call log_deallocate(tun_mask)

    call destroyAll(TRS)

  end subroutine calculate_single_transmission_N_contacts_sp

  subroutine calculate_single_transmission_N_contacts_dp(negf,ni,nf,ESH,SelfEneR,cblk,tun_proj,gsmr,Gr,TUN)
    type(TNegf), intent(in) :: negf
    integer, intent(in) :: ni,nf
    type(z_DNS), intent(in) :: SelfEneR(MAXNCONT)
    type(z_DNS), intent(inout) :: ESH(:,:)
    type(z_DNS), dimension(:),intent(in) :: gsmr
    type(z_DNS), dimension(:,:),intent(inout) :: Gr
    integer, intent(in) :: cblk(:)
    type(intarray), intent(in) :: tun_proj
    real(dp), intent(out) :: TUN

    !Work variables
    type(CublasHandle) :: hh
    Integer :: ct1, ct2, bl1, bl2, i, nbl
    logical, dimension(:), allocatable :: tun_mask
    Type(z_DNS) :: work1, work2, GAM1_dns, GAM2_dns, TRS
    real(dp) :: max
    complex(dp), parameter :: one = (1.0_dp, 0.0_dp)
    complex(dp), parameter :: mone = (-1.0_dp, 0.0_dp)
    complex(dp), parameter :: zero = (0.0_dp, 0.0_dp)

    !Arrange contacts in way that order between first and second is always the
    !same (always ct1 < ct2)

    if (cblk(ni).lt.cblk(nf)) then
       ct1=ni;ct2=nf;
    else
       ct1=nf;ct2=ni;
    endif

    bl1=cblk(ct1); bl2=cblk(ct2);
    nbl = size(cblk)
    hh = negf%hcublas
    ! in this way nt1 < nt2 by construction
    if ( nbl.gt.1 .and. (bl2-bl1).gt.1) then

       ! Compute column-blocks of Gr(i,bl1) up to i=bl2
       ! Gr(i,bl1) = -gr(i) T(i,i-1) Gr(i-1,bl1)
       do i = bl1+1, bl2
          !Checks whether previous block is non null.
          !If so next block is also null => TUN = 0
          call copyFromGPU(Gr(i-1,bl1))
          max=maxval(abs(Gr(i-1,bl1)%val))

          if (max.lt.EPS) then
             TUN = EPS*EPS !for log plots
             !Destroy also the block adjecent to diagonal since
             !this is not deallocated anymore in calling subroutine
             if (i.gt.(bl1+1)) call destroyAll(Gr(i-1,bl1))
             return
          endif

          !Checks whether block has been created, if not do it
          if (.not.allocated(Gr(i,bl1)%val)) then

             call createAll(work1, gsmr(i)%nrow, ESH(i,i-1)%ncol)
             call createAll(Gr(i,bl1), work1%nrow, Gr(i-1,bl1)%ncol)
             call copyToGPU(ESH(i,i-1))
             call matmul_gpu(hh, mone, gsmr(i), ESH(i,i-1), zero, work1)
             call deleteGPU(ESH(i,i-1))

             call matmul_gpu(hh, one, work1, Gr(i-1,bl1) ,zero, Gr(i,bl1))
             call destroyAll(work1)

          endif

          ! avoid destroying blocks closer to diagonal
          if (i.gt.(bl1+2)) call destroyAll(Gr(i-1,bl1))
       end do

    endif
    ! Computes the Gamma matrices
    call createAll(GAM1_dns, SelfEneR(ct1)%nrow, SelfEneR(ct1)%ncol)
    call createAll(GAM2_dns, SelfEneR(ct2)%nrow, SelfEneR(ct2)%ncol)
    call spectral_gpu(hh, SelfEneR(ct1), GAM1_dns)
    call spectral_gpu(hh, SelfEneR(ct2), GAM2_dns)

    ! Work to compute transmission matrix (Gamma2 Gr Gamma1 Ga)
    call createAll(work1, Gr(bl2,bl1)%nrow, GAM1_dns%ncol)
    call matmul_gpu(hh, one, Gr(bl2,bl1), GAM1_dns, zero, work1)

    call createAll(work2, GAM2_dns%nrow, work1%ncol)
    call matmul_gpu(hh, one, GAM2_dns, work1, zero, work2)

    call destroyAll(work1)
    call destroyAll(GAM2_dns)
    call destroyAll(GAM1_dns)

    call createAll(TRS, work2%nrow, Gr(bl2,bl1)%nrow)

    call matmul_gpu(hh, one, work2, Gr(bl2,bl1), zero, TRS, 'dag_2nd')
    call destroyAll(work2)
    if (bl2.gt.bl1+1) call destroyAll(Gr(bl2,bl1))

    !call get_tun_mask(ESH, bl2, tun_proj, tun_mask)
    !call trace_gpu(hh, TRS, TUN, tun_mask)
    call trace_gpu(hh, TRS, TUN)
    !call log_deallocate(tun_mask)

    call destroyAll(TRS)

  end subroutine calculate_single_transmission_N_contacts_dp

! Based on projection indices build a logical mask just on contact block
  subroutine get_tun_mask_sp(ESH,nbl,tun_proj,tun_mask)
    Type(c_DNS), intent(in) :: ESH(:,:)
    integer, intent(in) :: nbl
    type(intarray), intent(in) :: tun_proj
    logical, intent(out), allocatable :: tun_mask(:)

    integer :: ii, istart, iend, ind

    call log_allocate(tun_mask, ESH(nbl,nbl)%nrow)

    if (allocated(tun_proj%indexes)) then
       tun_mask = .false.

       ! set the start/end indices of nbl
       ! NB: istart has offset -1 to avoid +/-1 operations
       istart = 0
       do ii = 1, nbl-1
          istart = istart + ESH(ii,ii)%nrow
       end do
       iend = istart + ESH(nbl,nbl)%nrow + 1

       ! select the indices in tun_proj

       do ii = 1, size(tun_proj%indexes)
          ind = tun_proj%indexes(ii)
          if (ind > istart .and. ind < iend) then
             tun_mask(ind - istart) = .true.
          end if
       end do
    else
       tun_mask = .true.
    end if

  end subroutine get_tun_mask_sp

  subroutine get_tun_mask_dp(ESH,nbl,tun_proj,tun_mask)
    Type(z_DNS), intent(in) :: ESH(:,:)
    integer, intent(in) :: nbl
    type(intarray), intent(in) :: tun_proj
    logical, intent(out), allocatable :: tun_mask(:)

    integer :: ii, istart, iend, ind

    call log_allocate(tun_mask, ESH(nbl,nbl)%nrow)

    if (allocated(tun_proj%indexes)) then
       tun_mask = .false.

       ! set the start/end indices of nbl
       ! NB: istart has offset -1 to avoid +/-1 operations
       istart = 0
       do ii = 1, nbl-1
          istart = istart + ESH(ii,ii)%nrow
       end do
       iend = istart + ESH(nbl,nbl)%nrow + 1

       ! select the indices in tun_proj

       do ii = 1, size(tun_proj%indexes)
          ind = tun_proj%indexes(ii)
          if (ind > istart .and. ind < iend) then
             tun_mask(ind - istart) = .true.
          end if
       end do
    else
       tun_mask = .true.
    end if

  end subroutine get_tun_mask_dp

!---------------------------------------------------

  subroutine init_tridiag_blk(Matrix,S)
    type(z_DNS), dimension(:,:) :: Matrix,S

    integer :: nbl, j

    nbl = SIZE(Matrix,1)

    call create(Matrix(1,1),S(1,1)%nrow,S(1,1)%ncol)
    Matrix(1,1)%val=(0.0_dp,0.0_dp)
    do j=2,nbl-1
       call create(Matrix(j-1,j),S(j-1,j)%nrow,S(j-1,j)%ncol)
       Matrix(j-1,j)%val=(0.0_dp,0.0_dp)

       call create(Matrix(j,j),S(j,j)%nrow,S(j,j)%ncol)
       Matrix(j,j)%val=(0.0_dp,0.0_dp)

       call create(Matrix(j,j-1),S(j,j-1)%nrow,S(j,j-1)%ncol)
       Matrix(j,j-1)%val=(0.0_dp,0.0_dp)
    end do
    if (nbl.gt.1) then
       call create(Matrix(nbl,nbl),S(nbl,nbl)%nrow,S(nbl,nbl)%ncol)
       Matrix(nbl,nbl)%val=(0.0_dp,0.0_dp)

       call create(Matrix(nbl-1,nbl),S(nbl-1,nbl)%nrow,S(nbl-1,nbl)%ncol)
       Matrix(nbl-1,nbl)%val=(0.0_dp,0.0_dp)

       call create(Matrix(nbl,nbl-1),S(nbl,nbl-1)%nrow,S(nbl,nbl-1)%ncol)
       Matrix(nbl,nbl-1)%val=(0.0_dp,0.0_dp)
    endif

  end subroutine init_tridiag_blk

!---------------------------------------------------

  subroutine destroy_tridiag_blk(M)
    type(z_DNS), dimension(:,:), allocatable :: M

    integer :: i, nbl

    if (.not.allocated(M)) return

    nbl=size(M,1)

    do i=1,nbl
      if (allocated(M(i,i)%val)) then
        call destroyAll(M(i,i))
      end if
    end do
    do i=2,nbl
      if (allocated(M(i-1,i)%val)) then
        call destroyAll(M(i-1,i))
      end if
      if (allocated(M(i,i-1)%val)) then
        call destroyAll(M(i,i-1))
      end if
    end do

  end subroutine destroy_tridiag_blk

!---------------------------------------------------

  subroutine allocate_blk_dns(blkM,nbl)
    type(z_DNS), dimension(:,:), allocatable :: blkM
    integer :: nbl, ierr

    allocate(blkM(nbl,nbl),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate block-Matrix'

  end subroutine allocate_blk_dns

  !---------------------------------------------------

  subroutine deallocate_blk_dns(blkM)
    type(z_DNS), dimension(:,:), allocatable :: blkM
    integer :: ierr

    deallocate(blkM,stat=ierr)
    if (ierr.ne.0) stop 'DEALLOCATION ERROR: could not deallocate block-Matrix'

  end subroutine deallocate_blk_dns

end module iterative_gpu
