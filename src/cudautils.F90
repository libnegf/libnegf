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

module cudautils
   use ln_precision
   use mat_def
   use lib_param
   use iso_c_binding
   use, intrinsic :: ieee_arithmetic
   implicit none
   private

   public :: createGPU
   public :: copyToGPU
   public :: copyFromGPU
   public :: deleteGPU
   public :: createAll
   public :: destroyAll

   public :: copy_trid_toGPU
   public :: copy_trid_toHOST
   public :: copy_vdns_toGPU
   public :: delete_vdns_fromGPU
   public :: delete_trid_fromGPU

   public :: matmul_gpu
   public :: inverse_gpu
   public :: matsum_gpu
   public :: kernelsum_gpu
   public :: spectral_gpu
   public :: init_gpu
   public :: trace_gpu
   public :: copy_mat_gpu
   public :: asum_gpu
   public :: dagger_gpu

   public :: checksum

   interface createGPU
      module procedure createGPU_sp
      module procedure createGPU_dp
   end interface createGPU

   interface deleteGPU
      module procedure deleteGPU_sp
      module procedure deleteGPU_dp
   end interface deleteGPU

   interface createAll
      module procedure createAll_sp
      module procedure createAll_dp
   end interface createAll

   interface destroyAll
      module procedure destroyAll_sp
      module procedure destroyAll_dp
   end interface destroyAll

   interface copyToGPU
      module procedure copyToGPU_sp
      module procedure copyToGPU_dp
   end interface copyToGPU

   interface copyFromGPU
      module procedure copyFromGPU_sp
      module procedure copyFromGPU_dp
   end interface copyFromGPU



   interface copy_trid_toGPU
      module procedure copy_trid_toGPU_sp
      module procedure copy_trid_toGPU_dp
   end interface copy_trid_toGPU

   interface copy_trid_toHOST
      module procedure copy_trid_toHOST_sp
      module procedure copy_trid_toHOST_dp
   end interface copy_trid_toHOST

   interface delete_vdns_fromGPU
      module procedure delete_vdns_fromGPU_sp
      module procedure delete_vdns_fromGPU_dp
   end interface delete_vdns_fromGPU

   interface copy_vdns_toGPU
      module procedure copy_vdns_toGPU_sp
      module procedure copy_vdns_toGPU_dp
   end interface copy_vdns_toGPU

   interface delete_trid_fromGPU
      module procedure delete_trid_fromGPU_sp
      module procedure delete_trid_fromGPU_dp
   end interface delete_trid_fromGPU



   interface matmul_gpu
      module procedure matmul_gpu_sp
      module procedure matmul_gpu_dp
   end interface matmul_gpu

   interface inverse_gpu
      module procedure inverse_gpu_sp
      module procedure inverse_gpu_dp
   end interface inverse_gpu

   interface matsum_gpu
      module procedure matsum_gpu_sp
      module procedure matsum_gpu_dp
   end interface matsum_gpu

   interface kernelsum_gpu
      module procedure kernelsum_gpu_sp
      module procedure kernelsum_gpu_dp
   end interface kernelsum_gpu

   interface spectral_gpu
      module procedure spectral_gpu_sp
      module procedure spectral_gpu_dp
   end interface spectral_gpu

   interface init_gpu
      module procedure init_gpu_sp
      module procedure init_gpu_dp
   end interface init_gpu

   interface trace_gpu
      module procedure trace_gpu_sp
      module procedure trace_gpu_dp
   end interface trace_gpu

   interface copy_mat_gpu
      module procedure copy_mat_gpu_sp
      module procedure copy_mat_gpu_dp
   end interface copy_mat_gpu

   interface asum_gpu
      module procedure asum_gpu_sp
      module procedure asum_gpu_dp
   end interface asum_gpu

   interface dagger_gpu
      module procedure dagger_gpu_sp
      module procedure dagger_gpu_dp
   end interface dagger_gpu

   integer, parameter :: REAL_SIZE=4
   integer, parameter :: COMPLEX_SIZE=2*REAL_SIZE
   integer, parameter :: DOUBLE_SIZE=8
   integer, parameter :: DOUBLE_COMPLEX_SIZE=2*DOUBLE_SIZE
   ! Notes: type(c_ptr) -> void**
   !        type(c_ptr), value -> void*

   interface
     integer(c_int) function cu_createMat(d_A, siz) bind(C, name='cu_createMat')
       use iso_c_binding
       ! not by value since pointer has to be initialized, hence pass its reference
       type(c_ptr) :: d_A
       integer(c_int), value :: siz
     end function cu_createMat

     integer(c_int) function cu_copyMatD2H(h_A, d_A, siz) bind(C, name='cu_copyMatD2H')
       use iso_c_binding
       type(c_ptr), value :: h_A
       type(c_ptr), value :: d_A
       integer(c_int), value :: siz
     end function cu_copyMatD2H

     integer(c_int) function cu_copyMatH2D(h_A, d_A, siz) bind(C, name='cu_copyMatH2D')
       use iso_c_binding
       type(c_ptr), value :: h_A
       type(c_ptr), value :: d_A
       integer(c_int), value :: siz
     end function cu_copyMatH2D

     integer(c_int) function cu_deleteMat(d_A) bind(C, name='cu_deleteMat')
       use iso_c_binding
       type(c_ptr), value :: d_A
     end function cu_deleteMat

     ! C = alpha*A*B + beta*C
     integer(c_int) function cu_CmultMat(hcublas, m, n, k, alpha, d_A, d_B, beta, d_C, dagger) &
                  &   bind(C, name='cu_CmultMat')
       use iso_c_binding
       import cublasHandle
       type(cublasHandle), value :: hcublas
       integer(c_int), value :: m
       integer(c_int), value :: n
       integer(c_int), value :: k
       complex(c_float_complex) :: alpha
       type(c_ptr), value :: d_A
       type(c_ptr), value :: d_B
       complex(c_float_complex) :: beta
       type(c_ptr), value :: d_C
       integer(c_int), value :: dagger
     end function cu_CmultMat

     integer(c_int) function cu_ZmultMat(hcublas, m, n, k, alpha, d_A, d_B, beta, d_C, dagger) &
                  &   bind(C, name='cu_ZmultMat')
       use iso_c_binding
       import cublasHandle
       type(cublasHandle), value :: hcublas
       integer(c_int), value :: m
       integer(c_int), value :: n
       integer(c_int), value :: k
       complex(c_double_complex) :: alpha
       type(c_ptr), value :: d_A
       type(c_ptr), value :: d_B
       complex(c_double_complex) :: beta
       type(c_ptr), value :: d_C
       integer(c_int), value :: dagger
     end function cu_ZmultMat

     integer(c_int) function cu_Cinverse(hcublas, hcusolver, d_A, d_Ainv, N) &
                  &   bind(C, name='cu_Cinverse')
       use iso_c_binding
       import cublasHandle
       import cusolverDnHandle
       type(cusolverDnHandle), value :: hcusolver
       type(cublasHandle), value :: hcublas
       integer(c_int), value :: N
       type(c_ptr), value :: d_A
       type(c_ptr), value :: d_Ainv
     end function cu_Cinverse

     integer(c_int) function cu_Zinverse(hcublas, hcusolver, d_A, d_Ainv, N) &
                  &   bind(C, name='cu_Zinverse')
       use iso_c_binding
       import cublasHandle
       import cusolverDnHandle
       type(cusolverDnHandle), value :: hcusolver
       type(cublasHandle), value :: hcublas
       integer(c_int), value :: N
       type(c_ptr), value :: d_A
       type(c_ptr), value :: d_Ainv
     end function cu_Zinverse

     integer(c_int) function cu_Ckernelsum(d_C,alpha,d_A,beta,d_B,msize) &
                  &   bind(C, name='cu_Ckernelsum')
       use iso_c_binding
       type(c_ptr), value :: d_C
       complex(c_float_complex) :: alpha
       complex(c_float_complex) :: beta
       type(c_ptr), value :: d_A
       type(c_ptr), value :: d_B
       integer(c_int), value :: msize
     end function cu_Ckernelsum

     integer(c_int) function cu_Zkernelsum(d_C,alpha,d_A,beta,d_B,msize) &
                  &   bind(C, name='cu_Zkernelsum')
       use iso_c_binding
       type(c_ptr), value :: d_C
       complex(c_double_complex) :: alpha
       complex(c_double_complex) :: beta
       type(c_ptr), value :: d_A
       type(c_ptr), value :: d_B
       integer(c_int), value :: msize
     end function cu_Zkernelsum

     integer(c_int) function cu_Cmatsum(hcublas, m, n, alpha, d_A, beta, d_B, d_C, dagger) &
                  &   bind(C, name='cu_Cmatsum')
       use iso_c_binding
       import cublasHandle
       type(cublasHandle), value :: hcublas
       integer(c_int), value :: m
       integer(c_int), value :: n
       integer(c_int), value :: dagger
       complex(c_float_complex) :: alpha
       type(c_ptr), value :: d_A
       type(c_ptr), value :: d_B
       complex(c_float_complex) :: beta
       type(c_ptr), value :: d_C
     end function cu_Cmatsum

     integer(c_int) function cu_Zmatsum(hcublas, m, n, alpha, d_A, beta, d_B, d_C, dagger) &
                  &   bind(C, name='cu_Zmatsum')
       use iso_c_binding
       import cublasHandle
       type(cublasHandle), value :: hcublas
       integer(c_int), value :: m
       integer(c_int), value :: n
       integer(c_int), value :: dagger
       complex(c_double_complex) :: alpha
       type(c_ptr), value :: d_A
       type(c_ptr), value :: d_B
       complex(c_double_complex) :: beta
       type(c_ptr), value :: d_C
     end function cu_Zmatsum

     integer(c_int) function cu_Cinitmat(d_A, nrow) &
                  &   bind(C, name='cu_Cinitmat')
       use iso_c_binding
       type(c_ptr), value :: d_A
       integer(c_int), value :: nrow
     end function cu_Cinitmat

     integer(c_int) function cu_Zinitmat(d_A, nrow) &
                  &   bind(C, name='cu_Zinitmat')
       use iso_c_binding
       type(c_ptr), value :: d_A
       integer(c_int), value :: nrow
     end function cu_Zinitmat

     real(c_float) function cu_Ctrace(hcublas, d_A, nrow, h_tun, mask_present) &
                  &   bind(C, name='cu_Ctrace')
       use iso_c_binding
       import cublasHandle
       type(cublasHandle), value :: hcublas
       type(c_ptr), value :: d_A
       type(c_ptr), value :: h_tun
       integer(c_int), value :: nrow
       integer(c_int), value :: mask_present
     end function cu_Ctrace

     real(c_double) function cu_Ztrace(hcublas, d_A, nrow, h_tun, mask_present) &
                  &   bind(C, name='cu_Ztrace')
       use iso_c_binding
       import cublasHandle
       type(cublasHandle), value :: hcublas
       type(c_ptr), value :: d_A
       type(c_ptr), value :: h_tun
       integer(c_int), value :: nrow
       integer(c_int), value :: mask_present
     end function cu_Ztrace

     integer(c_int) function cu_Cmatcopy(hcublas, d_A, d_B, msize) &
                  &   bind(C, name='cu_Cmatcopy')
       use iso_c_binding
       import cublasHandle
       type(cublasHandle), value :: hcublas
       type(c_ptr), value :: d_A
       type(c_ptr), value :: d_B
       integer(c_int), value :: msize
     end function cu_Cmatcopy

     integer(c_int) function cu_Zmatcopy(hcublas, d_A, d_B, msize) &
                  &   bind(C, name='cu_Zmatcopy')
       use iso_c_binding
       import cublasHandle
       type(cublasHandle), value :: hcublas
       type(c_ptr), value :: d_A
       type(c_ptr), value :: d_B
       integer(c_int), value :: msize
     end function cu_Zmatcopy

     integer(c_int) function cu_Casum(hcublas, d_A, summ, N) &
                  &   bind(C, name='cu_Casum')
       use iso_c_binding
       import cublasHandle
       type(cublasHandle), value :: hcublas
       type(c_ptr), value :: d_A
       real(c_float) :: summ
       integer(c_int), value :: N
     end function cu_Casum

     integer(c_int) function cu_Zasum(hcublas, d_A, summ, N) &
                  &   bind(C, name='cu_Zasum')
       use iso_c_binding
       import cublasHandle
       type(cublasHandle), value :: hcublas
       type(c_ptr), value :: d_A
       real(c_double) :: summ
       integer(c_int), value :: N
     end function cu_Zasum

end interface


 contains

!~-~-~-~-~-~-~-~-~-~-~-~ DATA MOVEMENT ROUTINES  ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
   !-~-~-~-~ Single precision  ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
   subroutine createGPU_sp(A)
     type(c_DNS), intent(in) :: A

     integer :: err
     err = cu_createMat(A%d_addr, size(A%val)*COMPLEX_SIZE)
   end subroutine createGPU_sp

   subroutine copyToGPU_sp(A)
     type(c_DNS), intent(in), target :: A
     integer :: err
     call createGPU(A)
     err = cu_copyMatH2D(c_loc(A%val), A%d_addr, size(A%val)*COMPLEX_SIZE)
   end subroutine copyToGPU_sp

   subroutine copyFromGPU_sp(A)
     type(c_DNS), intent(in), target :: A
     integer :: err
     err = cu_copyMatD2H(c_loc(A%val), A%d_addr, size(A%val)*COMPLEX_SIZE)
   end subroutine copyFromGPU_sp

   subroutine deleteGPU_sp(A)
     type(c_DNS), intent(in) :: A
     integer :: err
     err = cu_deleteMat(A%d_addr)
   end subroutine deleteGPU_sp

   subroutine createAll_sp(A, nrow, ncol)
     type(c_DNS) :: A
     integer, intent(in) :: nrow, ncol

     call create(A, nrow, ncol)
     call createGPU(A)
   end subroutine createAll_sp

   subroutine destroyAll_sp(A)
     type(c_DNS) :: A
     call deleteGPU(A)
     call destroy(A)
   end subroutine destroyAll_sp

   !-~-~-~-~ Double precision  ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
   subroutine createGPU_dp(A)
     type(z_DNS), intent(in) :: A

     integer :: err
     err = cu_createMat(A%d_addr, size(A%val)*DOUBLE_COMPLEX_SIZE)
   end subroutine createGPU_dp

   subroutine copyToGPU_dp(A)
     type(z_DNS), intent(in), target :: A
     integer :: err
     call createGPU(A)
     err = cu_copyMatH2D(c_loc(A%val), A%d_addr, size(A%val)*DOUBLE_COMPLEX_SIZE)
   end subroutine copyToGPU_dp

   subroutine copyFromGPU_dp(A)
     type(z_DNS), intent(in), target :: A
     integer :: err
     err = cu_copyMatD2H(c_loc(A%val), A%d_addr, size(A%val)*DOUBLE_COMPLEX_SIZE)
   end subroutine copyFromGPU_dp

   subroutine deleteGPU_dp(A)
     type(z_DNS), intent(in) :: A
     integer :: err
     err = cu_deleteMat(A%d_addr)
   end subroutine deleteGPU_dp

   subroutine createAll_dp(A, nrow, ncol)
     type(z_DNS) :: A
     integer, intent(in) :: nrow, ncol

     call create(A, nrow, ncol)
     call createGPU(A)
   end subroutine createAll_dp

   subroutine destroyAll_dp(A)
     type(z_DNS) :: A
     call deleteGPU(A)
     call destroy(A)
   end subroutine destroyAll_dp

!~-~-~-~-~-~-~-~-~-~-~-~ MATRIX COMPUTATION ROUTINES  ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

   !-~-~-~-~ Single precision  ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
   ! C = alpha*A*B + beta*C
   subroutine matmul_gpu_sp(hcublas, alpha, A, B, beta, C, dagger)
     type(cublasHandle), intent(in) :: hcublas
     complex(sp), intent(in) :: alpha
     type(c_DNS), intent(in) :: A
     type(c_DNS), intent(in) :: B
     complex(sp), intent(in) :: beta
     type(c_DNS), intent(inout) :: C
     character(*), intent(in), optional :: dagger

     integer :: istat

     if (.not.present(dagger)) then
     istat = cu_CmultMat(hcublas, C%nrow, C%ncol, A%ncol, alpha, A%d_addr, &
             & B%d_addr, beta, C%d_addr, 0)
     else
       select case(dagger)    
       case('dag_1st')
         istat = cu_CmultMat(hcublas, C%nrow, C%ncol, B%nrow, alpha, A%d_addr, &
               & B%d_addr, beta, C%d_addr, 1)
       case('dag_2nd')
         istat = cu_CmultMat(hcublas, C%nrow, C%ncol, A%ncol, alpha, A%d_addr, &
               & B%d_addr, beta, C%d_addr, 2)
       case default  
         error stop 'Error in matmul_gpu'    
       end select
     endif

   end subroutine matmul_gpu_sp

   subroutine inverse_gpu_sp(hcublas, hcusolver, A, Ainv, err)
     type(cublasHandle) :: hcublas
     type(cusolverDnHandle) :: hcusolver
     type(c_DNS), intent(in) :: A
     type(c_DNS), intent(inout) :: Ainv
     integer, intent(out) :: err

     integer :: istat

     call init_gpu_sp(Ainv)

     istat = cu_Cinverse(hcublas,hcusolver, A%d_addr, Ainv%d_addr, A%nrow)
     err = istat

   end subroutine inverse_gpu_sp

   subroutine kernelsum_gpu_sp(C,alpha,A,beta,B)
     type(c_DNS), intent(inout) :: C
     type(c_DNS), intent(in) :: A
     type(c_DNS), intent(in) :: B
     complex(sp), intent(in) :: alpha
     complex(sp), intent(in) :: beta

     integer :: istat

     istat = cu_Ckernelsum(C%d_addr, alpha, A%d_addr, beta, B%d_addr, size(A%val))

   end subroutine kernelsum_gpu_sp

   subroutine matsum_gpu_sp(hcublas, alpha, A, beta, B, C, dagger)
     type(cublasHandle), intent(in) :: hcublas
     complex(sp), intent(in) :: alpha
     type(c_DNS), intent(in) :: A
     type(c_DNS), intent(in) :: B
     complex(sp), intent(in) :: beta
     type(c_DNS), intent(inout) :: C
     character(*), intent(in), optional :: dagger

     integer :: istat

     if (.not.present(dagger)) then
       istat = cu_Cmatsum(hcublas, C%nrow, C%ncol, alpha, A%d_addr, beta, B%d_addr, C%d_addr, 0)
     else
       if (dagger == 'dag_1st') then
         istat = cu_Cmatsum(hcublas, C%nrow, C%ncol, alpha, A%d_addr, beta, B%d_addr, C%d_addr, 1)
       endif
       if (dagger == 'dag_2nd') then
         istat = cu_Cmatsum(hcublas, C%nrow, C%ncol, alpha, A%d_addr, beta, B%d_addr, C%d_addr, 2)
       endif
     endif

   end subroutine matsum_gpu_sp

   subroutine spectral_gpu_sp(hcublas, G_in, G_out)
      type(cublasHandle), intent(in) :: hcublas
      type(c_DNS), intent(in) :: G_in
      type(c_DNS), intent(inout) :: G_out

      complex(sp) :: alpha, beta

      alpha = cmplx(0.0, 1.0, sp)
      beta = cmplx(0.0, -1.0, sp)

      call matsum_gpu_sp(hcublas, alpha, G_in, beta, G_in, G_out, 'dag_2nd')

   end subroutine spectral_gpu_sp

   subroutine init_gpu_sp(A)
      type(c_DNS), intent(inout) :: A

      integer :: istat

      istat = cu_Cinitmat(A%d_addr, A%nrow)

   end subroutine init_gpu_sp

   subroutine trace_gpu_sp(hcublas, A, trace, tun_mask)
      type(cublasHandle), intent(in) :: hcublas
      type(c_DNS), intent(inout) :: A
      real(sp), intent(out) :: trace
      logical, intent(in), optional, target :: tun_mask(:)

      type(c_ptr) :: dummy

      if (.not.present(tun_mask)) then
         trace = cu_Ctrace(hcublas, A%d_addr, A%nrow, dummy, 0)
      else
         trace = cu_Ctrace(hcublas, A%d_addr, A%nrow, c_loc(tun_mask), 1)
      endif

   end subroutine trace_gpu_sp

   subroutine copy_mat_gpu_sp(hcublas, A, Acopy)
      type(cublasHandle), intent(in) :: hcublas
      type(c_DNS), intent(inout) :: A
      type(c_DNS), intent(inout) :: Acopy

      integer :: istat

      istat = cu_Cmatcopy(hcublas, A%d_addr, Acopy%d_addr, A%nrow*A%ncol)

   end subroutine copy_mat_gpu_sp

   subroutine asum_gpu_sp(hcublas, A, summ)
      type(cublasHandle), intent(in) :: hcublas
      type(c_DNS), intent(in) :: A
      real(sp), intent(out) :: summ

      integer :: istat

      istat = cu_Casum(hcublas, A%d_addr, summ, A%nrow*A%ncol)

   end subroutine asum_gpu_sp

   subroutine dagger_gpu_sp(hcublas, G, G_dag)
      type(cublasHandle), intent(in) :: hcublas
      type(c_DNS), intent(in) :: G
      type(c_DNS), intent(inout) :: G_dag

      complex(sp) :: alpha, beta

      alpha = cmplx(0.0, 0.0, sp)
      beta = cmplx(1.0, 0.0, sp)

      call matsum_gpu_sp(hcublas, alpha, G, beta, G, G_dag, 'dag_2nd')

  end subroutine dagger_gpu_sp

   !-~-~-~-~ Double precision  ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
   ! C = alpha AB + beta C
   subroutine matmul_gpu_dp(hcublas, alpha, A, B, beta, C, dagger)
     type(cublasHandle), intent(in) :: hcublas
     complex(dp), intent(in) :: alpha
     type(z_DNS), intent(in) :: A
     type(z_DNS), intent(in) :: B
     complex(dp), intent(in) :: beta
     type(z_DNS), intent(inout) :: C
     character(*), intent(in), optional :: dagger

     integer :: istat

     if (.not.present(dagger)) then
       istat = cu_ZmultMat(hcublas, C%nrow, C%ncol, A%ncol, alpha, A%d_addr, &
             & B%d_addr, beta, C%d_addr, 0)
     else
       select case(dagger)    
       case('dag_1st')
         istat = cu_ZmultMat(hcublas, C%nrow, C%ncol, B%nrow, alpha, A%d_addr, &
               & B%d_addr, beta, C%d_addr, 1)
       case('dag_2nd')
         istat = cu_ZmultMat(hcublas, C%nrow, C%ncol, A%ncol, alpha, A%d_addr, &
               & B%d_addr, beta, C%d_addr, 2)
       case default  
         error stop 'Error in matmul_gpu'    
       end select
     endif

   end subroutine matmul_gpu_dp

   subroutine inverse_gpu_dp(hcublas, hcusolver, A, Ainv, err)
     type(cublasHandle) :: hcublas
     type(cusolverDnHandle) :: hcusolver
     type(z_DNS), intent(in) :: A
     type(z_DNS), intent(inout) :: Ainv
     integer, intent(out) :: err

     integer :: istat

     call init_gpu_dp(Ainv)

     istat = cu_Zinverse(hcublas,hcusolver, A%d_addr, Ainv%d_addr, A%nrow)
     err = istat

   end subroutine inverse_gpu_dp

   subroutine kernelsum_gpu_dp(C,alpha,A,beta,B)
     type(z_DNS), intent(inout) :: C
     type(z_DNS), intent(in) :: A
     type(z_DNS), intent(in) :: B
     complex(dp), intent(in) :: alpha
     complex(dp), intent(in) :: beta

     integer :: istat

     istat = cu_Zkernelsum(C%d_addr, alpha, A%d_addr, beta, B%d_addr, size(A%val))

   end subroutine kernelsum_gpu_dp

   subroutine matsum_gpu_dp(hcublas, alpha, A, beta, B, C, dagger)
     type(cublasHandle), intent(in) :: hcublas
     complex(dp), intent(in) :: alpha
     type(z_DNS), intent(in) :: A
     type(z_DNS), intent(in) :: B
     complex(dp), intent(in) :: beta
     type(z_DNS), intent(inout) :: C
     character(*), intent(in), optional :: dagger

     integer :: istat

     if (.not.present(dagger)) then
       istat = cu_Zmatsum(hcublas, C%nrow, C%ncol, alpha, A%d_addr, beta, B%d_addr, C%d_addr, 0)
     else
       if (dagger == 'dag_1st') then
         istat = cu_Zmatsum(hcublas, C%nrow, C%ncol, alpha, A%d_addr, beta, B%d_addr, C%d_addr, 1)
       endif
       if (dagger == 'dag_2nd') then
         istat = cu_Zmatsum(hcublas, C%nrow, C%ncol, alpha, A%d_addr, beta, B%d_addr, C%d_addr, 2)
       endif
     endif

   end subroutine matsum_gpu_dp

   subroutine spectral_gpu_dp(hcublas, G_in, G_out)
      type(cublasHandle), intent(in) :: hcublas
      type(z_DNS), intent(in) :: G_in
      type(z_DNS), intent(inout) :: G_out

      complex(dp) :: alpha, beta

      alpha = cmplx(0.0, 1.0, dp)
      beta = cmplx(0.0, -1.0, dp)

      call matsum_gpu_dp(hcublas, alpha, G_in, beta, G_in, G_out, 'dag_2nd')

   end subroutine spectral_gpu_dp

   subroutine init_gpu_dp(A)
      type(z_DNS), intent(inout) :: A

      integer :: istat

      istat = cu_Zinitmat(A%d_addr, A%nrow)

   end subroutine init_gpu_dp

   subroutine trace_gpu_dp(hcublas, A, trace, tun_mask)
      type(cublasHandle), intent(in) :: hcublas
      type(z_DNS), intent(inout) :: A
      real(dp), intent(out) :: trace
      logical, intent(in), optional, target :: tun_mask(:)

      type(c_ptr) :: dummy

      if (.not.present(tun_mask)) then
         trace = cu_Ztrace(hcublas, A%d_addr, A%nrow, dummy, 0)
      else
         trace = cu_Ztrace(hcublas, A%d_addr, A%nrow, c_loc(tun_mask), 1)
      endif

   end subroutine trace_gpu_dp

   subroutine copy_mat_gpu_dp(hcublas, A, Acopy)
      type(cublasHandle), intent(in) :: hcublas
      type(z_DNS), intent(inout) :: A
      type(z_DNS), intent(inout) :: Acopy

      integer :: istat

      istat = cu_Zmatcopy(hcublas, A%d_addr, Acopy%d_addr, A%nrow*A%ncol)

   end subroutine copy_mat_gpu_dp

   subroutine asum_gpu_dp(hcublas, A, summ)
      type(cublasHandle), intent(in) :: hcublas
      type(z_DNS), intent(in) :: A
      real(dp), intent(out) :: summ

      integer :: istat

      istat = cu_Zasum(hcublas, A%d_addr, summ, A%nrow*A%ncol)

   end subroutine asum_gpu_dp

   subroutine dagger_gpu_dp(hcublas, G, G_dag)
      type(cublasHandle), intent(in) :: hcublas
      type(z_DNS), intent(in) :: G
      type(z_DNS), intent(inout) :: G_dag

      complex(dp) :: alpha, beta

      alpha = cmplx(0.0, 0.0, dp)
      beta = cmplx(1.0, 0.0, dp)

      call matsum_gpu_dp(hcublas, alpha, G, beta, G, G_dag, 'dag_2nd')

  end subroutine dagger_gpu_dp

!~-~-~-~-~-~-~-~-~-~-~-~ MATRIX MOVEMENTS ROUTINES  ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
   !-~-~-~-~ Single precision  ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
  subroutine copy_trid_toGPU_sp(M)
    type(c_DNS), dimension(:,:), intent(in) :: M
    integer :: ii, nbl

    nbl = size(M,1)
    call copyToGPU(M(1,1))
    do ii=2,nbl

       call copyToGPU(M(ii,ii))
       call copyToGPU(M(ii-1,ii))
       call copyToGPU(M(ii,ii-1))

    end do
  end subroutine copy_trid_toGPU_sp

  subroutine delete_trid_fromGPU_sp(M)
    type(c_DNS), dimension(:,:), intent(inout) :: M
    integer :: ii, nbl

    nbl = size(M,1)
    call deleteGPU(M(1,1))
    do ii=2,nbl
       call deleteGPU(M(ii,ii))
       call deleteGPU(M(ii-1,ii))
       call deleteGPU(M(ii,ii-1))
    end do
  end subroutine delete_trid_fromGPU_sp

  subroutine copy_trid_toHOST_sp(M)
    type(c_DNS), dimension(:,:), intent(inout) :: M
    integer :: ii, nbl

    nbl = size(M,1)
    call copyFromGPU(M(1,1))
    do ii=2,nbl
       call copyFromGPU(M(ii,ii))
       call copyFromGPU(M(ii-1,ii))
       call copyFromGPU(M(ii,ii-1))
    end do
  end subroutine copy_trid_toHOST_sp

  subroutine copy_vdns_toGPU_sp(V)
    type(c_DNS), dimension(:), intent(in) :: V
    integer :: ii, nbl

    nbl = size(V)
    do ii=1,nbl
       if(allocated(V(ii)%val)) then
          call copyToGPU(V(ii))
       endif
    end do
  end subroutine copy_vdns_toGPU_sp

  subroutine delete_vdns_fromGPU_sp(V)
    type(c_DNS), dimension(:) :: V
    integer :: ii, nbl

    nbl = size(V)
    do ii=1,nbl
       if(allocated(V(ii)%val)) then
          call deleteGPU(V(ii))
       endif
    end do
  end subroutine delete_vdns_fromGPU_sp

  !-~-~-~-~ Double precision  ~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
  subroutine copy_trid_toGPU_dp(M)
    type(z_DNS), dimension(:,:), intent(in) :: M
    integer :: ii, nbl

    nbl = size(M,1)
    do ii=2,nbl

       call copyToGPU(M(ii,ii-1))
       call copyToGPU(M(ii-1,ii))
       call copyToGPU(M(ii,ii))

    end do
    call copyToGPU(M(1,1))
  end subroutine copy_trid_toGPU_dp

  subroutine delete_trid_fromGPU_dp(M)
    type(z_DNS), dimension(:,:), intent(inout) :: M
    integer :: ii, nbl

    nbl = size(M,1)
    call deleteGPU(M(1,1))
    do ii=2,nbl
       call deleteGPU(M(ii,ii))
       call deleteGPU(M(ii-1,ii))
       call deleteGPU(M(ii,ii-1))
    end do
  end subroutine delete_trid_fromGPU_dp

  subroutine copy_trid_toHOST_dp(M)
    type(z_DNS), dimension(:,:), intent(inout) :: M
    integer :: ii, nbl

    nbl = size(M,1)
    call copyFromGPU(M(1,1))
    do ii=2,nbl
       call copyFromGPU(M(ii,ii))
       call copyFromGPU(M(ii-1,ii))
       call copyFromGPU(M(ii,ii-1))
    end do
  end subroutine copy_trid_toHOST_dp

  subroutine copy_vdns_toGPU_dp(V)
    type(z_DNS), dimension(:), intent(in) :: V
    integer :: ii, nbl

    nbl = size(V)
    do ii=1,nbl
       if(allocated(V(ii)%val)) then
          call copyToGPU(V(ii))
       endif
    end do
  end subroutine copy_vdns_toGPU_dp

  subroutine delete_vdns_fromGPU_dp(V)
    type(z_DNS), dimension(:) :: V
    integer :: ii, nbl

    nbl = size(V)
    do ii=1,nbl
       if(allocated(V(ii)%val)) then
          call deleteGPU(V(ii))
       endif
    end do
  end subroutine delete_vdns_fromGPU_dp

  subroutine checksum(hcublas, A, nome)
    type(cublasHandle), intent(in) :: hcublas
    type(z_DNS), intent(in) :: A
    character(*), intent(in) :: nome

    real(dp) :: summ

    call asum_gpu(hcublas, A, summ)

    if (ieee_is_nan(summ)) then
       write(*,*) 'GPU:   ',trim(nome),summ
    end if

  end subroutine checksum
end module cudautils
