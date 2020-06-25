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


!****************************************************
!                                                   |
!Collection of subroutines for matrix inversion     |
!                                                   |
!****************************************************

Module inversions

use ln_precision
use ln_allocation 
use mat_def
use sparsekit_drv
private

INTEGER :: t1_i,t2_I,cr_i,cm_i,t1_ii,t2_ii,cr_ii,cm_ii
LOGICAL, PARAMETER :: timing=.FALSE.

! SPARSKIT iterative solvers removed (Alex)
!public :: ZMSKINVP_LA,  ZMSKINVP_MA, zINVP_MA
public :: zINV_LAPACK
#:if defined("__SUPERLU")
public :: zINV_LU
#:endif
#:if defined("__PARDISO")
public :: zINV_PARDISO
#:endif

public :: compGreen ! wrapper to different type of computations
public :: inverse, block2Green, block3Green
public :: zlin

interface compGreen
  module procedure compGreen_dns
  module procedure compGreen_arr
end interface

interface inverse
   module procedure cinv, zinv, rinv
end interface


contains

  !--------------------------------------------------------------------------
  subroutine compGreen_dns(G,A,n)
    Type(z_DNS) :: A, G
    Integer :: n

    Integer :: sel, iter

    sel = 1 
    !if(A%nrow.gt.100) sel = 2 
    !if(A%nrow.gt.1200) sel = 3     

    select case(sel)
    case(1)
       call inverse(G%val,A%val,n)
    case(2)
       iter = 1     
       call block2Green(G%val,A%val,n,iter)
    case(3)
       call block3Green(G%val,A%val,n)
    end select

  end subroutine compGreen_dns
  !--------------------------------------------------------------------------
  subroutine compGreen_arr(G,A,n)
    complex(dp), dimension(:,:) :: A, G
    Integer :: n

    Integer :: sel, iter

    sel = 1 

    select case(sel)
    case(1)
       call inverse(G,A,n)
    case(2)
       iter = 1     
       call block2Green(G,A,n,iter)
    case(3)
       call block3Green(G,A,n)
    end select

  end subroutine compGreen_arr


!******************************************************
!                                                     | 
!  PGMRES based inversion subroutine with masking     |
!                                                     |
!******************************************************

!subroutine ZMSKINVP_LA(A_csr, M_csc, INV_csc)

!*************************************************************************
!Calculate the inverted of A_csr matrix in CSR format masked             |
!with M_csc matrix in CSC format. Result is INV_csc in CSC format.       |
!                                                                        |
!Input:                                                                  |
!A_csr: Matrix to be inverted in Real CSR format                         |
!M_csc: Mask matrix in csc format                                        |
!Output:                                                                 |
!INV_csc: Inverted masked matrix                                         | 
!                                                                        | 
!Note for allocation: INV_csc has the same dimension of M_csc            |
!*************************************************************************

!integer :: i, j, index_start, index_end, row_index, nzval_address  
!integer :: ierr, LU_iwk, iout
!type(z_CSR) :: A_csr
!complex(kind=dp), DIMENSION(:), ALLOCATABLE :: y 
!complex(kind=dp), DIMENSION(:), ALLOCATABLE :: x 
!type(z_MSR) :: LU_msr
!integer, DIMENSION(:), ALLOCATABLE :: LU_levs 
!integer, DIMENSION(:), ALLOCATABLE :: LU_ju 
!complex(kind=dp), DIMENSION(:), ALLOCATABLE :: vv 
!type(z_CSC) :: M_csc, INV_csc

!Define parameters for PGMRES and FILTER calls
!im = krylov subspace dimension for PGMRES solver
!maxits = maximum iterations number allowed in PGMRES solver
!eps = maximum allowed error in PGMRES solver
!
!integer :: im=15, maxits=35
!real(kind=dp), PARAMETER :: eps=1e-16
!
!
!write(*,*) "Allocation of LU_msr, LU_ju e LU_levs for preconditioning"
!  
!  LU_iwk=A_csr%nnz+1
!  call create(LU_msr,A_csr%nrow,A_csr%ncol,LU_iwk)
!  call log_allocate(LU_ju,A_csr%nrow)
!  call log_allocate(LU_levs,LU_iwk) 
!  !End Allocation
!    
!write(*,*)  "Call LU0 preconditioning for pgmres solver"
!  call  ziluk_st(A_csr, 0, LU_msr, LU_ju, LU_levs, LU_iwk) 
!  
!   call log_deallocate(LU_levs)
!  
!  nzval_address=1
!  INV_csc%rowind=M_csc%rowind
!  INV_csc%colpnt=M_csc%colpnt
!
!write(*,*) "Allocations and initializations for PGMRES solver"
!write(*,*) im*A_csr%nrow+1
!  iout=0
!
!  call log_allocate(vv,im*A_csr%nrow+1) 
! call log_allocate(x,A_csr%nrow)
! call log_allocate(y,A_csr%nrow)  
! 
! x(:)=(0.d0,0.d0) !Solution initial guess
!
! do i=1,A_csr%nrow
!
!    y(:)=(0.d0,0.d0)
!    y(i)=(1.d0,0.d0)
!
!    !Call iterative approssimative solver PGMRES
!    call zpgmres(A_csr%nrow, im, y, x, vv, eps, maxits, iout, A_csr%nzval, A_csr%colind, &
!          A_csr%rowpnt, LU_msr%nzval, LU_msr%index, LU_ju, ierr)
!
!    if (ierr.ne.0) then
!       write(*,*)'Error number ',ierr,' in pgmres linear solver'
!    endif
!
!    index_start=M_csc%colpnt(i)
!    index_end=M_csc%colpnt(i+1)
!
!    do j=index_start,index_end-1
!
!       row_index=M_csc%rowind(j)
!       INV_csc%nzval(nzval_address)=x(row_index)
!       nzval_address=nzval_address+1
!    enddo
!
! enddo
! 
!write(*,*) "Work arrays and structures deallocation"
! call destroy(LU_msr)
! call log_deallocate(vv)
! call log_deallocate(LU_ju)
! call log_deallocate(x)
!  call log_deallocate(y)

 
!END subroutine ZMSKINVP_LA
 
 
!*************************************************************************
!                                                                        | 
!  PGMRES based inversion subroutine with masking and more arguments     |
!                                                                        |
!*************************************************************************
!
!subroutine ZMSKINVP_MA(A_csr, M_csc, INV_csc, im, maxits, eps)
!
!
!!*************************************************************************
!!Calculate the inverted of A_csr matrix in CSR format masked             |
!!with M_csc matrix in CSC format. Result is INV_csc in CSC format.       |
!!                                                                        |
!!Input:                                                                  |
!!A_csr: Matrix to be inverted in Real CSR format                         |
 !M_csc: Mask matrix in csc format                                        |
!!im: krylov subspace dimension for PGMRES solver                         |
!!maxits: maximum iterations number for PGMRES solver                     |
!!eps: error allowed in PGMRES solver                                     |
!!Output:                                                                 |
!!INV_csc: Inverted masked matrix                                         |
!!*************************************************************************
!
!integer :: i, j, index_start, index_end, row_index, nzval_address  
!integer :: ierr, LU_iwk, iout
!type(z_CSR) :: A_csr
!complex(kind=dp), DIMENSION(:), ALLOCATABLE :: y 
!complex(kind=dp), DIMENSION(:), ALLOCATABLE :: x 
!type(z_MSR) :: LU_msr
!integer, DIMENSION(:), ALLOCATABLE :: LU_levs 
!integer, DIMENSION(:), ALLOCATABLE :: LU_ju 
!complex(kind=dp), DIMENSION(:), ALLOCATABLE :: vv 
!type(z_CSC) :: M_csc, INV_csc
!
!integer :: im, maxits
!real(kind=dp) :: eps
!
!
!write(*,*) "Allocation of LU_msr, LU_ju e LU_levs for preconditioning"
!  
!  LU_iwk=A_csr%nnz+1
!  call create(LU_msr,A_csr%nrow,A_csr%ncol,LU_iwk)
!  call log_allocate(LU_ju,A_csr%nrow)
!  call log_allocate(LU_levs,LU_iwk) 
!  !End Allocation
!    
!write(*,*)  "Call LU0 preconditioning for pgmres solver"
!  call  ziluk_st(A_csr, 0, LU_msr, LU_ju, LU_levs, LU_iwk) 
!  
!  call log_deallocate(LU_levs)
!  
!  nzval_address=1
!  INV_csc%rowind=M_csc%rowind
!  INV_csc%colpnt=M_csc%colpnt
!
!write(*,*) "Allocations and initializations for PGMRES solver"
!write(*,*) im*A_csr%nrow+1
!  iout=0
!
!  call log_allocate(vv,im*A_csr%nrow+1) 
!  call log_allocate(x,A_csr%nrow)
!  call log_allocate(y,A_csr%nrow)  
!  
!  x(:)=(0.d0, 0.d0) !Solution initial guess
!
!  do i=1,A_csr%nrow
!
!     y(:)=(0.d0, 0.d0)
!     y(i)=(1.d0, 0.d0)
!
!     !Call iterative approssimative solver PGMRES
!     call zpgmres(A_csr%nrow, im, y, x, vv, eps, maxits, iout, A_csr%nzval, A_csr%colind, &
!          A_csr%rowpnt, LU_msr%nzval, LU_msr%index, LU_ju, ierr)
!
!     if (ierr.ne.0) then
!        write(*,*)'Error number ',ierr,' in pgmres linear solver'
!     endif
!
!     index_start=M_csc%colpnt(i)
!     index_end=M_csc%colpnt(i+1)
!
!     do j=index_start,index_end-1
!
!         row_index=M_csc%rowind(j)
!        INV_csc%nzval(nzval_address)=x(row_index)
!        nzval_address=nzval_address+1
!     enddo

!  enddo
  
! write(*,*) "Work arrays and structures deallocation"
!  call destroy(LU_msr)
!  call log_deallocate(vv)
!  call log_deallocate(LU_ju)
!  call log_deallocate(x)
!  call log_deallocate(y)
!
!
!END subroutine ZMSKINVP_MA


!**************************************************************
!                                                             |
!  PGMRES based inversion without masking and more arguments  |
!                                                             |
!**************************************************************
!
!subroutine zINVP_MA(A_csr, INV, nrow, im, maxits, eps)
!
!!*************************************************************************
!!Calculate the inverted of A_csr matrix in CSR format                    |
!!Result is INV in dense format.                                          |
!!                                                                        |
!!Input:                                                                  |
!!A_csr: Matrix to be inverted in Real CSR format                         |                                       
!!im: krylov subspace dimension for PGMRES solver                         |
!!maxits: maximum iterations number for PGMRES solver                     |
!!eps: error allowed in PGMRES solver                                     |
!!nrow: INV matrix number of rows                                         |
!!Output:                                                                 |
!!INV: Inverted matrix                                                    |
 !*************************************************************************
!
!integer :: i
!integer :: nrow, ierr, LU_iwk, iout
!type(z_CSR) :: A_csr
!complex(kind=dp), DIMENSION(nrow,nrow) :: INV
!complex(kind=dp), DIMENSION(:), ALLOCATABLE :: y 
!complex(kind=dp), DIMENSION(:), ALLOCATABLE :: x 
!type(z_MSR) :: LU_msr
!integer, DIMENSION(:), ALLOCATABLE :: LU_levs 
!integer, DIMENSION(:), ALLOCATABLE :: LU_ju 
!complex(kind=dp), DIMENSION(:), ALLOCATABLE :: vv 
!
!integer :: im, maxits
!real(kind=dp) :: eps
!
!if (nrow.ne.A_csr%nrow) then
!  stop 'Error in INVP_MA: INV must have same number of rows of A_csr'
!endif
!  
!write(*,*) "Allocation of LU_msr, LU_ju e LU_levs for preconditioning"
!  
!  LU_iwk=A_csr%nnz+1
!  call create(LU_msr,A_csr%nrow,A_csr%ncol,LU_iwk-1)
!  call log_allocate(LU_ju,A_csr%nrow)
!  call log_allocate(LU_levs,LU_iwk) 
!  !End Allocation
!    
!write(*,*)  "Call LU0 preconditioning for pgmres solver"
!  call  ziluk_st(A_csr, 0, LU_msr, LU_ju, LU_levs, LU_iwk) 
!  
!  call log_deallocate(LU_levs)
!  
!!Debug
!write(*,*) "Allocations and initializations for PGMRES solver"
!!EndDebug
!  iout=0
!  
!  call log_allocate(vv,im*A_csr%nrow+1) 
!  call log_allocate(x,A_csr%nrow)
!  call log_allocate(y,A_csr%nrow)  
!  
!  x(:)=(0.d0, 0.d0) !Solution initial guess
!
!  do i=1,A_csr%nrow
!
!     y(:)=(0.d0, 0.d0)
!     y(i)=(1.d0, 0.d0)
!
!     !Call iterative approssimative solver PGMRES
!     call zpgmres(A_csr%nrow, im, y, x, vv, eps, maxits, iout, A_csr%nzval, A_csr%colind, &
!          A_csr%rowpnt, LU_msr%nzval, LU_msr%index, LU_ju, ierr)
!
!     if (ierr.ne.0) then
!        write(*,*)'Error number ',ierr,' in pgmres linear solver'
!        if (ierr.eq.1) then
!          write(*,*) 'Convergence not achieved in the allowed number of iterations'
!        endif
!        if (ierr.eq.(-1)) then
!          write(*,*) 'The initial guess (x=0) seems to be the exact solution'
!        endif
!     endif
!
!     INV(:,i)=x(:)  
!
!  enddo
!  
!  call destroy(LU_msr)
!  call log_deallocate(LU_ju)
!  call log_deallocate(x)
!  call log_deallocate(y)
!  call log_deallocate(vv)
   !Debug
!  write(*,*)'INVP_MA Done'
   !EndDebug
!
!end subroutine zINVP_MA
!

#:if defined("__SUPERLU")
!*********************************************
!                                            |
!  SuperLU based inversion without masking   |
!                                            | 
!*********************************************
subroutine zINV_LU(A_csr, INV)

!***********************************************************
!Calculate the inverse matrix                              |
!                                                          |
!Input:                                                    |
!A_csr: matrix to be inverted in CSR format                | 
!       (must be square)                                   |
!                                                          |
!Output:                                                   |
!INV: inverse matrix in dense format                       |
!                                                          |
!***********************************************************

  type(z_CSR) :: A_csr
  type(z_CSC) :: A_csc
  integer :: factors, iopt, ldb, nrhs, info

  complex(kind=dp), DIMENSION(:,:) :: INV
  integer :: i, mem
 
  if (timing) call SYSTEM_CLOCK(t1_i,cr_i,cm_i)
  if (timing) call SYSTEM_CLOCK(t1_ii,cr_ii,cm_ii)

  !Trasposizione da csr a csc 
  !(nota: A_csr viene lasciata allocata)
  call create(A_csc,A_csr%nrow,A_csr%ncol,A_csr%nnz)
  call csr2csc(A_csr,A_csc) 

  !Parameter and rhs inizialization
  nrhs=A_csc%nrow
  ldb=A_csc%nrow
  
  INV=(0.d0, 0.d0)
  
  FORALL (i=1:nrhs)  INV(i,i)=(1.D0,0.D0)
  
  if (timing) call SYSTEM_CLOCK(t2_ii,cr_ii,cm_i)
  IF (timing) WRITE(*,*) 'SuperLU pre-inversion operations',&
       (t2_ii-t1_ii)*1.0/cr_ii,'sec'
  
  !Chiamata Superlu 2.0
  call c_bridge_zgssv(nrhs, A_csc%nnz, nrhs, &
       A_csc%nzval, A_csc%rowind, A_csc%colpnt, &
       INV, ldb, info, mem)      

  !alloc_mem=alloc_mem+mem
  !if (alloc_mem.gt.peak_mem) peak_mem = alloc_mem
  !alloc_mem=alloc_mem-mem         
  
  !Deallocazione della matrice CSC
  call destroy(A_csc)
  
  if (timing) call SYSTEM_CLOCK(t2_i,cr_i,cm_i)
  IF (timing) WRITE(*,*) 'SuperLU inversions done in ',&
       (t2_i-t1_i)*1.0/cr_i,'sec'

end subroutine zINV_LU

#:endif

!!$!************************************************
!!$!                                               |
!!$!  SuperLU3.0 based inversion without masking   |
!!$!                                               | 
!!$!************************************************
!!$subroutine zINV_LU(A_csr, INV)
!!$  
!!$!***********************************************************
!!$!Calculate the inverse matrix                              |
!!$!                                                          |
!!$!Input:                                                    |
!!$!A_csr: matrix to be inverted in CSR format                | 
!!$!       (must be square)                                   |
!!$!                                                          |
!!$!Output:                                                   |
!!$!INV: inverse matrix in dense format                       |
!!$!                                                          |
!!$!***********************************************************
!!$
!!$  type(z_CSR) :: A_csr
!!$  type(z_CSC) :: A_csc
!!$  integer :: factors, iopt, ldb, nrhs, info
!!$
!!$  complex(kind=dp), DIMENSION(:,:) :: INV
!!$  integer :: i, mem
!!$ 
!!$if (timing) call SYSTEM_CLOCK(t1_i,cr_i,cm_i)
!!$if (timing) call SYSTEM_CLOCK(t1_ii,cr_ii,cm_ii)
!!$
!!$!Trasposizione da csr a csc (nota: A_csr viene lasciata allocata)
!!$call create(A_csc,A_csr%nrow,A_csr%ncol,A_csr%nnz)
!!$call zcsrcsc_st(A_csr,A_csc) 
!!$
!!$!Parameter and rhs inizialization
!!$nrhs=A_csc%nrow
!!$ldb=A_csc%nrow
!!$
!!$INV=(0.d0, 0.d0)
!!$FORALL (i=1:nrhs) INV(i,i)=(1.D0,0.D0)
!!$  
!!$IF (timing) call SYSTEM_CLOCK(t2_ii,cr_ii,cm_i)  
!!$IF (timing) WRITE(*,*) 'SuperLU pre-inversion operations',(t2_ii-t1_ii)*1.0/cr_ii,'sec'
!!$
!!$!First call, LU factorization
!!$iopt = 1
!!$call c_fortran_zgssv( iopt, A_csc%nrow, A_csc%nnz, nrhs, A_csc%nzval, A_csc%rowind, A_csc%colpnt, & 
!!$    INV, ldb, factors, info, mem)
!!$
!!$!Second call, solve system
!!$iopt=2
!!$call c_fortran_zgssv( iopt, A_csc%nrow, A_csc%nnz, nrhs, A_csc%nzval, A_csc%rowind, A_csc%colpnt, & 
!!$    INV, ldb, factors, info, mem)
!!$
!!$!Third call, free internal allocated memory
!!$iopt=3
!!$call c_fortran_zgssv( iopt, A_csc%nrow, A_csc%nnz, nrhs, A_csc%nzval, A_csc%rowind, A_csc%colpnt, & 
!!$    INV, ldb, factors, info, mem)
!!$
!!$alloc_mem=alloc_mem+mem
!!$ if (alloc_mem.gt.peak_mem) peak_mem = alloc_mem
!!$alloc_mem=alloc_mem-mem         
!!$
!!$!Deallocazione della matrice CSC
!!$call zdestroy_CSC(A_csc)
!!$
!!$if (timing) call SYSTEM_CLOCK(t2_i,cr_i,cm_i)
!!$IF (timing) WRITE(*,*) 'SuperLU inversions done in ',(t2_i-t1_i)*1.0/cr_i,'sec'
!!$
!!$end subroutine zINV_LU

#:if defined("__PARDISO") 
!***********************************************************
!
!  PARDISO Direct inversion
!
!***********************************************************
SUBROUTINE zINV_PARDISO(A_csr, ndim, INV)

  !***********************************************************
  !Calculate the inverse matrix                              |
  !                                                          |
  !Input:                                                    |
  !A_csr: matrix to be inverted in CSR format                | 
  !       (must be square)                                   |
  !                                                          |
  !Output:                                                   |
  !INV: inverse matrix in dense format                       |
  !                                                          |
  !***********************************************************

  type(z_CSR) :: A_csr
  integer :: ndim
  type(z_CSR) :: INV

  COMPLEX(dp), ALLOCATABLE, DIMENSION(:,:) :: X,B
  TYPE(z_vec), ALLOCATABLE, DIMENSION(:) :: RM
  INTEGER :: i1,i2

  integer(4) :: MTYPE,NRHS,cnt,indx,nnz
  INTEGER(4) :: info,MSGLVL,PHASE,MAXFCT,MNUM
  INTEGER(4) :: PT(64)
  INTEGER(4) :: IPARM(64)

  INTEGER, DIMENSION(:), ALLOCATABLE :: PERM, iwork, wind
  INTEGER :: iwork_lenght

  iwork_lenght=MAX(A_csr%nrow+1,2*A_csr%nnz)
  CALL log_allocate(iwork,iwork_lenght)
  CALL zcsort(A_csr%nrow,A_csr%nzval,A_csr%colind,A_csr%rowpnt,iwork,.TRUE.) 
  CALL log_deallocate(iwork) 

  NRHS = 1  ! Number or RHS
 
  CALL log_allocate(PERM,ndim)
  CALL log_allocate(B,ndim,NRHS)
  CALL log_allocate(X,ndim,NRHS)
  CALL log_allocate(wind,ndim)

  MTYPE=13 ! 6: Complex, symmetric; 13: Complex, non symmetric 

  IPARM(:)=0
  PT(:)=0
  do i1=1,ndim
     PERM(i1)=i1
  enddo

  !Setting IPARM
  IPARM(1)=1;  ! 0 sets all default values
  IPARM(2)=2;  ! 0 MDA, 2 ND (metis)
  IPARM(3)=1;  ! is number of processors - must be set
  IPARM(4)=0;  ! 10*L+K; K=0: LU; K=1: CGS; K=2: CG  acc (1E-L) 
  IPARM(5)=0;  ! 1 if user supplies permutation (PERM)
  IPARM(6)=0;  ! 1 => the solution in on B
  IPARM(8)=0;  ! input: max num of iterative refinement steps
  IPARM(7)=0;  ! output: num of iterative steps
  IPARM(10)=13;! |A|*1.0E-IPARM(10) pivot threshold
  IPARM(11)=1; ! Scaling vectors (For unsym or indef sym mat) 
  IPARM(13)=1; ! move large close to diag (use for highly indef)
  IPARM(18)=0; ! if <0 the solver will report nnz in L U
  IPARM(19)=0; ! if <0 the solver will report factoring MFlops
  IPARM(21)=0; ! (for sym indef) 0: 1x1 diagonal; 1: 2x2 pivoting
  IPARM(27)=0; ! 1 checks sorting within the rows
  IPARM(28)=0; ! 1 switch to single precision
  IPARM(60)=0; ! 0: in-core; 1|2: out-of-core (disk)


  MSGLVL = 0   ! No output prints

  MAXFCT = 1   ! Maximum number of factorizations in memory

  MNUM = 1     ! Actual matrix to factorize: 0<MNUM<=MAXFCT

  B = (0.d0, 0.d0)
  B(1,1) = (1.D0, 0.D0)

  if (timing) call SYSTEM_CLOCK(t2_i,cr_i,cm_i)
  ! ----------------------------------------------------------------------------------------
  PHASE = 12 ! Analysis, numerical factorization
  call PARDISO(PT, MAXFCT, MNUM, MTYPE, PHASE, ndim, A_csr%nzval,A_csr%rowpnt,A_csr%colind, &
       PERM, NRHS, IPARM, MSGLVL, B, X, info)

  IF(info.NE.0) THEN 
     WRITE(*,*) 'PARDISO Inversion error number ',info, ' in phase ',PHASE
     SELECT CASE (info)
     CASE (-1) 
        WRITE(*,*) 'Input inconsistent'
     CASE (-2)
        WRITE(*,*) 'No enough memory'
     CASE (-3)
        WRITE(*,*) 'Reordering problem'
     CASE (-4) 
        WRITE(*,*) 'Zero pivot'
     CASE (-5)
        WRITE(*,*) 'Unclassified internal error'
     CASE (-6)
        WRITE(*,*) 'Preordering failed'
     CASE (-7)
        WRITE(*,*) 'Diagonal matrix problem'
     CASE default
        WRITE(*,*) 'Invalid info value returned'
     END SELECT
  ENDIF

  if (timing) call SYSTEM_CLOCK(t2_i,cr_i,cm_i)
  IF (timing) WRITE(*,*) 'Pardiso LU done in ',(t2_i-t1_i)*1.0/cr_i,'sec'
  if (timing) call SYSTEM_CLOCK(t2_ii,cr_ii,cm_ii)

  ! --------------------------------------------------------------------------------------
  ! SOLVE ndim linear systems and compress the results in a CSR matrix
  allocate(RM(ndim))  
  PHASE = 33 ! solve
  nnz = 0

  do i1 = 1, ndim

     B(i1,1)=(1.d0,0.d0)  ! set rhs column
        
     call PARDISO(PT, MAXFCT, MNUM, MTYPE, PHASE, ndim, A_csr%nzval,A_csr%rowpnt, &
          A_csr%colind, PERM, NRHS, IPARM, MSGLVL, B, X, info)

     IF(info.NE.0) THEN
         WRITE(*,*) 'Error in solve'
         STOP
     ENDIF  
     ! counts all non-zero elements in this column 
     cnt = 0
     do i2=1,ndim
        if(ABS(X(i2,1)).gt.EPS10) then
            cnt = cnt + 1
            wind(cnt)=i2
         end if
     end do
     ! store the column in the ragged matrix RM
     call log_allocate(RM(i1)%val,cnt)
     call log_allocate(RM(i1)%ind,cnt)     
     RM(i1)%len=cnt 
     do i2=1,cnt
         nnz = nnz + 1
         RM(i1)%ind(i2) = wind(i2)
         RM(i1)%val(i2) = X(wind(i2),1)
     end do            

     B(i1,1)=(0.d0,0.d0)  ! unset rhs column

  end do

  ! unroll the ragged matrix into the CSR matrix.
  call create(INV,ndim,ndim,nnz)
  INV%rowpnt(1)=1
  indx = 0
  do i1 = 1, ndim
     do i2 = 1,RM(i1)%len   
        indx = indx + 1
        INV%colind(indx) = RM(i1)%ind(i2)
        INV%nzval(indx) = RM(i1)%val(i2)
     end do
     INV%rowpnt(i1+1) = indx+1
     call log_deallocate(RM(i1)%ind)
     call log_deallocate(RM(i1)%val)
  enddo
  deallocate(RM)

  ! -----------------------------------------------------------------------------------------
  !alloc_mem=alloc_mem+IPARM(15)*1024
  ! if (alloc_mem.gt.peak_mem) peak_mem = alloc_mem
  !alloc_mem=alloc_mem-IPARM(15)*1024
  ! -----------------------------------------------------------------------------------------
  PHASE = -1 !release all internal memory

  call PARDISO(PT, MAXFCT, MNUM, MTYPE, PHASE, ndim, A_csr%nzval,A_csr%rowpnt,A_csr%colind, &
       PERM, NRHS, IPARM, MSGLVL, B, INV, info)

  CALL log_deallocate(B)
  CALL log_deallocate(X)
  CALL log_deallocate(PERM)
  CALL log_deallocate(wind)

  if (timing) call SYSTEM_CLOCK(t2_ii,cr_ii,cm_ii)
  IF (timing) WRITE(*,*) 'Pardiso inversions done in ',(t2_ii-t1_ii)*1.0/cr_ii,'sec'

END SUBROUTINE zINV_PARDISO

#:endif

!---------- INTERFACE FOR LAPACK INVERSION (Complex MATRICES) -------------
subroutine zINV_LAPACK(A_csr, INV)

  type(z_CSR) :: A_csr
  integer :: ndim
  complex(kind=dp), DIMENSION(:,:) :: INV
  type(z_DNS) :: A_dns

  ! transform sparse into dense matrix
  call create(A_dns,A_csr%nrow,A_csr%ncol)
  call csr2dns(A_csr,A_dns)

  if(A_dns%nrow.ne.A_dns%ncol) then
     write(*,*) 'ERROR: nrow != ncol'
     stop
  endif

  ndim=A_dns%nrow

  call zINV(INV,A_dns%val,ndim)

  call destroy(A_dns)

end subroutine zINV_LAPACK


! --------------------------------------------------------------------------
! Single precision inversion routine - modified as zinv
! 
subroutine cinv(inA,A,n)
  implicit none
  complex(sp), dimension(:,:) :: inA, A
  integer :: n

  INTEGER :: ipiv(n),info,i 
  complex(sp), dimension(:,:), allocatable  :: LU
 
  call log_allocate(LU,n,n)
  
  LU=A
  call cgetrf( n, n, inA, n, ipiv, info )

  if (info.ne.0)  then
     write(*,*)
     write(*,*) 'ERROR in LU factorization (cgetrf)',info
     stop
  end if
  
  inA = (0.0_sp, 0.0_sp)
  do i = 1, n
    inA(i,i) = (1.0_sp, 0.0_sp)
  end do

  call cgetrs( 'N', n, n, LU, n, ipiv, inA, n, info )

  if (info.ne.0)  then
     write(*,*)
     write(*,*) 'ERROR in INVERSION (cgetrs)',info
     stop
  end if

  call log_deallocate(LU)
  
end subroutine cinv  


! --------------------------------------------------------------------------
subroutine zinv(inA,A,n)
   
  Implicit none
  complex(kind=dp), dimension(:,:) :: inA, A
  integer :: n

  INTEGER :: ipiv(n), info, i
  complex(dp), dimension(:,:), allocatable  :: LU 

  call log_allocate(LU,n,n)

  LU=A
  call zgetrf( n, n, LU, n, ipiv, info )

  if (info.ne.0)  then
     write(*,*)
     write(*,*) 'ERROR in LU factorization (zgetrf)',info
     stop
  end if

  inA = (0.0_dp, 0.0_dp)
  do i = 1, n
    inA(i,i) = (1.0_dp, 0.0_dp)
  end do

  call zgetrs( 'N', n, n, LU, n, ipiv, inA, n, info )
  if (info.ne.0)  then
     write(*,*)
     write(*,*) 'ERROR in INVERSION (zgetrs)',info
     stop
  end if

  call log_deallocate(LU)
  
end subroutine zinv

!---------------------------------------------------------------------------
! Modified like zinv
!---------------------------------------------------------------------------
subroutine rinv(inA,A,n)
  Implicit none
  real(kind=dp), dimension(:,:) :: inA, A
  integer :: n

  INTEGER :: ipiv(n), info, i
  complex(dp), dimension(:,:), allocatable  :: LU 
 
  call log_allocate(LU,n,n)
  
  LU=A
  call  dgetrf(n, n, LU, n, ipiv, info )

  if (info.ne.0)  then
     write(*,*)
     write(*,*) 'ERROR in LU factorization (dgetrf)',info
     stop
  end if
  
  inA = 0.0_dp
  do i = 1, n
    inA(i,i) = 1.0_dp
  end do

  call dgetrs( 'N', n, n, LU, n, ipiv, inA, n, info )
  if (info.ne.0)  then
     write(*,*)
     write(*,*) 'ERROR in INVERSION (dgetrs)',info
     stop
  end if

  call log_deallocate(LU)

end subroutine rinv


  !--------------------------------------------------------------------------
  !
  !  (G11   )|(      )
  !  (   G22)|(      )
  !  ----------------- 
  !  (      )|(G11   )
  !  (      )|(   G22)
  !
  recursive subroutine block2Green(G,A,n,iter)

    complex(dp), dimension(:,:), intent(out) :: G
    complex(dp), dimension(:,:), intent(in) :: A
    integer, intent(in) :: n
    integer, intent(inout) :: iter
    ! Locals
    complex(dp), dimension(:,:), allocatable :: B
    Type(z_DNS) :: G11,G12,G21,G22
    Type(z_dns) :: h11,gr22
    Type(z_DNS) :: work1,work2,work4
    Integer :: bl1, b22, n2
    Integer, parameter :: maxiter=1

    if (n.le.10) then
       call inverse(G,A,n)
       return
    endif

    bl1 = n/2
    b22 = n - bl1
    n2 = bl1 + 1 

    call create(gr22,b22,b22)

    if (iter.lt.maxiter) then
       iter = iter + 1
       allocate(B(n2:n,n2:n))
       B=A(n2:n,n2:n)
       call block2Green(gr22%val,B,b22,iter)
       deallocate(B)     
       iter = iter - 1
    else
       allocate(B(n2:n,n2:n))
       B=A(n2:n,n2:n)
       call inverse(gr22%val,B,b22)    !g22 in uscita
       deallocate(B)     
    endif

    allocate(B(1:bl1,n2:n))
    B=A(1:bl1,n2:n)
    call prealloc_mult(B,gr22%val,work1)
    deallocate(B)

    allocate(B(n2:n,1:bl1))
    B=A(n2:n,1:bl1)
    call prealloc_mult(work1%val,B,(-1.d0,0.d0),work2)
    deallocate(B)

    !keep work1 = T12 * g22

    call create(h11,bl1,bl1)

    h11%val = A(1:bl1,1:bl1) + work2%val

    call destroy(work2)

    call create(G11,bl1,bl1)

    if (iter.lt.maxiter) then
       iter = iter + 1
       call block2Green(G11%val,h11%val,bl1,iter)     
       iter = iter - 1
    else
       call inverse(G11%val,h11%val,bl1)        !blocco G11
    endif

    call destroy(h11)

    call prealloc_mult(G11%val,work1%val,(-1.d0,0.d0),G12)

    allocate(B(n2:n,1:bl1))
    B=A(n2:n,1:bl1)
    call prealloc_mult(gr22%val,B,work2)
    deallocate(B)

    call prealloc_mult(work2,G11,G21)    !blocco G21 

    call destroy(work2)  

    call prealloc_mult(G21%val,work1%val,work4)

    call destroy(work1)

    call prealloc_sum(gr22,work4,G22)          !blocco G22

    call destroy(work4)  

    call destroy(gr22)


    G(1:bl1,1:bl1) = G11%val
    G(1:bl1,n2:n) = G12%val
    G(n2:n,1:bl1) = -G21%val
    G(n2:n,n2:n) = G22%val


    call destroy(G11)
    call destroy(G12)
    call destroy(G21)
    call destroy(G22)

  end subroutine block2Green

  !--------------------------------------------------------------------------

  subroutine block3Green(G,A,n)

    complex(dp), dimension(:,:) :: A, G
    Type(z_DNS) :: G11,G12,G21,G22
    Type(z_DNS) :: G33,G13,G23,G32,G31
    Type(z_dns) :: h22,gr11,gr33
    Type(z_DNS) :: work1,work2,work3,work4
    Integer :: n, bl1, bl2, bl3, n2

    !bl1 = n/3
    !bl2 = n/3
    !bl3 = n - (n/3)*2
    bl1 = 428
    bl2 = 424
    bl3 = 428
    n2 = bl1 + 1
    n3 = bl1 + bl2 + 1

    call create(gr33,bl3,bl3)
    call create(gr11,bl1,bl1)
    call inverse(gr33%val,A(n3:n,n3:n),bl3)       !Blocco g33 in uscita
    call inverse(gr11%val,A(1:bl1,1:bl1),bl1)     !Blocco g11 in uscita

    call prealloc_mult(A(n2:n3-1,n3:n),gr33%val,work1)                          !Moltiplicatione
    call prealloc_mult(work1%val,A(n3:n,n2:n3-1),(-1.d0,0.d0),work2)            !-H23g33H32  
    call destroy(work1)
    call prealloc_mult(A(n2:n3-1,1:bl1),gr11%val,work1)                         !Moltiplicazione
    call prealloc_mult(work1%val,A(1:bl1,n2:n3-1),(-1.d0,0.d0),work3)           !-H21g11H12
    call destroy(work1)
    call create(h22,bl2,bl2)
    h22%val = A(n2:n3-1,n2:n3-1) + work2%val + work3%val                         !h22=H22-H23g33H32-H21g11H12
    call destroy(work2)
    call destroy(work3)
    call create(G22,bl2,bl2)
    call inverse(G22%val,h22%val,bl2)            !Blocco G22 in uscita, G22=(h22)^-1
    call destroy(h22)
    
    call prealloc_mult(gr33%val,A(n3:n,n2:n3-1),work1)                         
    call prealloc_mult(work1,G22,work2)
    call destroy(work1)
    call prealloc_mult(work2%val,A(n2:n3-1,n3:n),work3)
    call destroy(work2)
    call prealloc_mult(work3,gr33,work4)                      !Moltiplicazione
    call destroy(work3)                                       !g33*H32*G22*H23*g33          
    call prealloc_sum(gr33,work4,G33)                          !Blocco G33 in uscita  
    call destroy(work4)                               !G33=g33+g33H32G22H23g33

    call prealloc_mult(gr11%val,A(1:bl1,n2:n3-1),work1)
    call prealloc_mult(work1,G22,work2)
    call destroy(work1)
    call prealloc_mult(work2%val,A(n2:n3-1,1:bl1),work3)
    call destroy(work2)
    call prealloc_mult(work3,gr11,work4)                    !Moltiplicazione
    call destroy(work3)                                     !g11*H12*G22*H21*g11
    call prealloc_sum(gr11,work4,G11)                       !Blocco G11 in uscita
    call destroy(work4)                                 !G11=g11+g11H12G22H21g11

    call prealloc_mult(gr11%val,A(1:bl1,n2:n3-1),work1)
    call prealloc_mult(work1,G22,(-1.d0,0.d0),G12)          !Blocco G12 in uscita, G12=-g11*H12*G22  
    call destroy(work1) 
    
    call prealloc_mult(G22%val,A(n2:n3-1,1:bl1),work1)
    call prealloc_mult(work1,gr11,(-1.d0,0.d0),G21)         !Blocco G21 in uscita, G21=-G22*H21*g11
    call destroy(work1)

    call prealloc_mult(G22%val,A(n2:n3-1,n3:n),work1)
    call prealloc_mult(work1,gr33,(-1.d0,0.d0),G23)         !Blocco G23 in uscita, G23=-G22*H23*g33
    call destroy(work1)

    call prealloc_mult(gr33%val,A(n3:n,n2:n3-1),work1)
    call prealloc_mult(work1,G22,(-1.d0,0.d0),G32)          !Blocco G32 in uscita, G32=-g33*H32*G22
    call destroy(work1)  

    call prealloc_mult(gr11%val,A(1:bl1,n2:n3-1),work1)
    call prealloc_mult(work1,G23,(-1.d0,0.d0),G13)          !Blocco G13 in uscita, G13=-g11*H12*G23
    call destroy(work1)

    call prealloc_mult(gr33%val,A(n3:n,n2:n3-1),work1)
    call prealloc_mult(work1,G21,(-1.d0,0.d0),G31)          !Blocco G31 in uscita, G31=-g33*H32*G21
    call destroy(work1) 

    call destroy(gr11)
    call destroy(gr33)

    G(1:bl1,1:bl1) = G11%val
    G(n2:n3-1,n2:n3-1) = G22%val
    G(n3:n,n3:n) = G33%val
    G(1:bl1,n2:n3-1) = G12%val
    G(n2:n3-1,1:bl1) = G21%val
    G(n2:n3-1,n3:n) = G23%val
    G(n3:n,n2:n3-1) = G32%val
    G(1:bl1,n3:n) = G13%val
    G(n3:n,1:bl1) = G31%val

    call destroy(G11)
    call destroy(G22)
    call destroy(G33)
    call destroy(G12)
    call destroy(G21)
    call destroy(G23)
    call destroy(G32)
    call destroy(G13)
    call destroy(G31)

  end subroutine block3Green

  !---------------------------------------------------
  ! SOLVE linear system with matrix as RHS
  !
  ! A U = T
  !
  ! A is an n x n complex matrix
  ! 
  subroutine zlin(A,T,gT,n)
    complex(dp), dimension(:,:) :: A, T, gT
    integer :: n
  
    INTEGER :: ipiv(n),info, lwork, NB, nrhs
    integer, external :: ilaenv
    complex(dp), dimension(:), allocatable  :: work
    complex(dp), dimension(:,:), allocatable :: LU

    call log_allocate(LU,n,n)
    LU=A

    ! LU factorization
    call zgetrf( n, n, LU, n, ipiv, info )
  
    if (info.ne.0)  then
       write(*,*)
       write(*,*) 'ERROR in INVERSION part 1',info
       stop
    end if
  
    nrhs = size(T,2)
    lwork = n*ilaenv(1,'zgetrs','',N,-1,-1,-1)
    call log_allocate(work,lwork)

    gT = T
    call zgetrs( 'N', n, nrhs, LU, n, ipiv, gT, n, info )
    if (info.ne.0)  then
       write(*,*)
       write(*,*) 'ERROR in INVERSION part 2',info
       stop
    end if
  
    call log_deallocate(work)
    call log_deallocate(LU)

  end subroutine zlin

END MODULE
