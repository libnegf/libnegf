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

public :: ZMSKINVP_LA,  ZMSKINVP_MA, zINVP_MA
public :: zINV_LAPACK
#ifdef __SUPERLU
public :: zINV_LU
#endif
#ifdef __PARDISO
public :: zINV_PARDISO
#endif

public :: inverse

interface inverse
   module procedure zinv, rinv
end interface

contains

!DEBUGGED AND COMPLETE***************************************************

!******************************************************
!                                                     | 
!  PGMRES based inversion subroutine with masking     |
!                                                     |
!******************************************************

subroutine ZMSKINVP_LA(A_csr, M_csc, INV_csc)

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

integer :: i, j, INV_address, index_start, index_end, row_index, nzval_address  
integer :: ierr, LU_iwk, iout
type(z_CSR) :: A_csr
complex(kind=dp), DIMENSION(:), ALLOCATABLE :: y 
complex(kind=dp), DIMENSION(:), ALLOCATABLE :: x 
type(z_MSR) :: LU_msr
integer, DIMENSION(:), ALLOCATABLE :: LU_levs 
integer, DIMENSION(:), ALLOCATABLE :: LU_ju 
complex(kind=dp), DIMENSION(:), ALLOCATABLE :: vv 
type(z_CSC) :: M_csc, INV_csc

!Define parameters for PGMRES and FILTER calls
!im = krylov subspace dimension for PGMRES solver
!maxits = maximum iterations number allowed in PGMRES solver
!eps = maximum allowed error in PGMRES solver

integer :: im=15, maxits=35
real(kind=dp), PARAMETER :: eps=1e-16


write(*,*) "Allocation of LU_msr, LU_ju e LU_levs for preconditioning"
  
  LU_iwk=A_csr%nnz+1
  call create(LU_msr,A_csr%nrow,A_csr%ncol,LU_iwk)
  call log_allocate(LU_ju,A_csr%nrow)
  call log_allocate(LU_levs,LU_iwk) 
  !End Allocation
    
write(*,*)  "Call LU0 preconditioning for pgmres solver"
  call  ziluk_st(A_csr, 0, LU_msr, LU_ju, LU_levs, LU_iwk) 
  
  call log_deallocate(LU_levs)
  
  nzval_address=1
  INV_csc%rowind=M_csc%rowind
  INV_csc%colpnt=M_csc%colpnt

write(*,*) "Allocations and initializations for PGMRES solver"
write(*,*) im*A_csr%nrow+1
  iout=0

  call log_allocate(vv,im*A_csr%nrow+1) 
  call log_allocate(x,A_csr%nrow)
  call log_allocate(y,A_csr%nrow)  
  
  x(:)=(0.d0,0.d0) !Solution initial guess

  do i=1,A_csr%nrow

     y(:)=(0.d0,0.d0)
     y(i)=(1.d0,0.d0)

     !Call iterative approssimative solver PGMRES
     call zpgmres(A_csr%nrow, im, y, x, vv, eps, maxits, iout, A_csr%nzval, A_csr%colind, &
          A_csr%rowpnt, LU_msr%nzval, LU_msr%index, LU_ju, ierr)

     if (ierr.ne.0) then
        write(*,*)'Error number ',ierr,' in pgmres linear solver'
     endif

     index_start=M_csc%colpnt(i)
     index_end=M_csc%colpnt(i+1)

     do j=index_start,index_end-1

        row_index=M_csc%rowind(j)
        INV_csc%nzval(nzval_address)=x(row_index)
        nzval_address=nzval_address+1
     enddo

  enddo
  
 write(*,*) "Work arrays and structures deallocation"
  call destroy(LU_msr)
  call log_deallocate(vv)
  call log_deallocate(LU_ju)
  call log_deallocate(x)
  call log_deallocate(y)


END subroutine ZMSKINVP_LA


!*************************************************************************
!                                                                        | 
!  PGMRES based inversion subroutine with masking and more arguments     |
!                                                                        |
!*************************************************************************

subroutine ZMSKINVP_MA(A_csr, M_csc, INV_csc, im, maxits, eps)


!*************************************************************************
!Calculate the inverted of A_csr matrix in CSR format masked             |
!with M_csc matrix in CSC format. Result is INV_csc in CSC format.       |
!                                                                        |
!Input:                                                                  |
!A_csr: Matrix to be inverted in Real CSR format                         |
!M_csc: Mask matrix in csc format                                        |
!im: krylov subspace dimension for PGMRES solver                         |
!maxits: maximum iterations number for PGMRES solver                     |
!eps: error allowed in PGMRES solver                                     |
!Output:                                                                 |
!INV_csc: Inverted masked matrix                                         |
!*************************************************************************

integer :: i, j, INV_address, index_start, index_end, row_index, nzval_address  
integer :: ierr, LU_iwk, iout
type(z_CSR) :: A_csr
complex(kind=dp), DIMENSION(:), ALLOCATABLE :: y 
complex(kind=dp), DIMENSION(:), ALLOCATABLE :: x 
type(z_MSR) :: LU_msr
integer, DIMENSION(:), ALLOCATABLE :: LU_levs 
integer, DIMENSION(:), ALLOCATABLE :: LU_ju 
complex(kind=dp), DIMENSION(:), ALLOCATABLE :: vv 
type(z_CSC) :: M_csc, INV_csc

integer :: im, maxits
real(kind=dp) :: eps


write(*,*) "Allocation of LU_msr, LU_ju e LU_levs for preconditioning"
  
  LU_iwk=A_csr%nnz+1
  call create(LU_msr,A_csr%nrow,A_csr%ncol,LU_iwk)
  call log_allocate(LU_ju,A_csr%nrow)
  call log_allocate(LU_levs,LU_iwk) 
  !End Allocation
    
write(*,*)  "Call LU0 preconditioning for pgmres solver"
  call  ziluk_st(A_csr, 0, LU_msr, LU_ju, LU_levs, LU_iwk) 
  
  call log_deallocate(LU_levs)
  
  nzval_address=1
  INV_csc%rowind=M_csc%rowind
  INV_csc%colpnt=M_csc%colpnt

write(*,*) "Allocations and initializations for PGMRES solver"
write(*,*) im*A_csr%nrow+1
  iout=0

  call log_allocate(vv,im*A_csr%nrow+1) 
  call log_allocate(x,A_csr%nrow)
  call log_allocate(y,A_csr%nrow)  
  
  x(:)=(0.d0, 0.d0) !Solution initial guess

  do i=1,A_csr%nrow

     y(:)=(0.d0, 0.d0)
     y(i)=(1.d0, 0.d0)

     !Call iterative approssimative solver PGMRES
     call zpgmres(A_csr%nrow, im, y, x, vv, eps, maxits, iout, A_csr%nzval, A_csr%colind, &
          A_csr%rowpnt, LU_msr%nzval, LU_msr%index, LU_ju, ierr)

     if (ierr.ne.0) then
        write(*,*)'Error number ',ierr,' in pgmres linear solver'
     endif

     index_start=M_csc%colpnt(i)
     index_end=M_csc%colpnt(i+1)

     do j=index_start,index_end-1

        row_index=M_csc%rowind(j)
        INV_csc%nzval(nzval_address)=x(row_index)
        nzval_address=nzval_address+1
     enddo

  enddo
  
 write(*,*) "Work arrays and structures deallocation"
  call destroy(LU_msr)
  call log_deallocate(vv)
  call log_deallocate(LU_ju)
  call log_deallocate(x)
  call log_deallocate(y)


END subroutine ZMSKINVP_MA


!**************************************************************
!                                                             |
!  PGMRES based inversion without masking and more arguments  |
!                                                             |
!**************************************************************

subroutine zINVP_MA(A_csr, INV, nrow, im, maxits, eps)

!*************************************************************************
!Calculate the inverted of A_csr matrix in CSR format                    |
!Result is INV in dense format.                                          |
!                                                                        |
!Input:                                                                  |
!A_csr: Matrix to be inverted in Real CSR format                         |                                       
!im: krylov subspace dimension for PGMRES solver                         |
!maxits: maximum iterations number for PGMRES solver                     |
!eps: error allowed in PGMRES solver                                     |
!nrow: INV matrix number of rows                                         |
!Output:                                                                 |
!INV: Inverted matrix                                                    |
!*************************************************************************

integer :: i, j 
integer :: nrow, ierr, LU_iwk, iout
type(z_CSR) :: A_csr
complex(kind=dp), DIMENSION(nrow,nrow) :: INV
complex(kind=dp), DIMENSION(:), ALLOCATABLE :: y 
complex(kind=dp), DIMENSION(:), ALLOCATABLE :: x 
type(z_MSR) :: LU_msr
integer, DIMENSION(:), ALLOCATABLE :: LU_levs 
integer, DIMENSION(:), ALLOCATABLE :: LU_ju 
complex(kind=dp), DIMENSION(:), ALLOCATABLE :: vv 

integer :: im, maxits
real(kind=dp) :: eps

if (nrow.ne.A_csr%nrow) then
  stop 'Error in INVP_MA: INV must have same number of rows of A_csr'
endif
  
write(*,*) "Allocation of LU_msr, LU_ju e LU_levs for preconditioning"
  
  LU_iwk=A_csr%nnz+1
  call create(LU_msr,A_csr%nrow,A_csr%ncol,LU_iwk-1)
  call log_allocate(LU_ju,A_csr%nrow)
  call log_allocate(LU_levs,LU_iwk) 
  !End Allocation
    
write(*,*)  "Call LU0 preconditioning for pgmres solver"
  call  ziluk_st(A_csr, 0, LU_msr, LU_ju, LU_levs, LU_iwk) 
  
  call log_deallocate(LU_levs)
  
!Debug
write(*,*) "Allocations and initializations for PGMRES solver"
!EndDebug
  iout=0
  
  call log_allocate(vv,im*A_csr%nrow+1) 
  call log_allocate(x,A_csr%nrow)
  call log_allocate(y,A_csr%nrow)  
  
  x(:)=(0.d0, 0.d0) !Solution initial guess

  do i=1,A_csr%nrow

     y(:)=(0.d0, 0.d0)
     y(i)=(1.d0, 0.d0)

     !Call iterative approssimative solver PGMRES
     call zpgmres(A_csr%nrow, im, y, x, vv, eps, maxits, iout, A_csr%nzval, A_csr%colind, &
          A_csr%rowpnt, LU_msr%nzval, LU_msr%index, LU_ju, ierr)

     if (ierr.ne.0) then
        write(*,*)'Error number ',ierr,' in pgmres linear solver'
        if (ierr.eq.1) then
          write(*,*) 'Convergence not achieved in the allowed number of iterations'
        endif
        if (ierr.eq.(-1)) then
          write(*,*) 'The initial guess (x=0) seems to be the exact solution'
        endif
     endif

     INV(:,i)=x(:)  

  enddo
  
  call destroy(LU_msr)
  call log_deallocate(LU_ju)
  call log_deallocate(x)
  call log_deallocate(y)
  call log_deallocate(vv)
  !Debug
  write(*,*)'INVP_MA Done'
  !EndDebug

end subroutine zINVP_MA

#ifdef __SUPERLU
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

  alloc_mem=alloc_mem+mem
  if (alloc_mem.gt.peak_mem) peak_mem = alloc_mem
  alloc_mem=alloc_mem-mem         
  
  !Deallocazione della matrice CSC
  call destroy(A_csc)
  
  if (timing) call SYSTEM_CLOCK(t2_i,cr_i,cm_i)
  IF (timing) WRITE(*,*) 'SuperLU inversions done in ',&
       (t2_i-t1_i)*1.0/cr_i,'sec'

end subroutine zINV_LU

#endif

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

#ifdef __PARDISO 
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
  INTEGER :: i1,i2, mem

  integer(4) :: MTYPE,NRHS,cnt,indx,nnz
  INTEGER(4) :: info,MSGLVL,PHASE,MAXFCT,MNUM
  INTEGER(4) :: PT(64)
  INTEGER(4) :: IPARM(64)

  INTEGER, DIMENSION(:), ALLOCATABLE :: PERM, iwork, wind
  INTEGER :: iwork_lenght
  real(dp), parameter :: EPS10 = 1.d-10

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
  alloc_mem=alloc_mem+IPARM(15)*1024
   if (alloc_mem.gt.peak_mem) peak_mem = alloc_mem
   alloc_mem=alloc_mem-IPARM(15)*1024
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

#endif

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
subroutine zinv(inA,A,n)
   
  Implicit none
  complex(kind=dp), dimension(:,:) :: inA, A
  integer :: n

  !integer :: lwork=n
  INTEGER :: ipiv(n),info
  complex(kind=dp) :: work(n)
  
  inA=A
  call zgetrf( n, n, inA, n, ipiv, info )

  if (info.ne.0)  then
     write(*,*)
     write(*,*) 'ERROR in INVERSION part 1',info
     stop
  end if

  call zgetri( n, inA, n, ipiv, work, n, info )
  if (info.ne.0)  then
     write(*,*)
     write(*,*) 'ERROR in INVERSION part 2',info
     stop
  end if
  
  return
  
end subroutine zinv
!---------------------------------------------------------------------------
!---------- INTERFACE FOR LAPACK INVERSION (SYMMETRIC MATRICES) -------------
subroutine rinv(inA,A,n)
  Implicit none
  real(kind=dp), dimension(:,:) :: inA, A
  integer :: n

  INTEGER :: ipiv(n),info
  real(dp) :: work(n)
  
  inA=A
  call  dgetrf(n, n, inA, n, ipiv, info )
  if (info.ne.0)  then
     write(*,*) 'ERROR in INVERSION part 1',info
     stop
  end if
  call dgetri( n, inA, n, ipiv, work, n, info )
  if (info.ne.0)  then
     write(*,*) 'ERROR in INVERSION part 2',info
     stop
  end if
  
  return
  
end subroutine rinv



END MODULE
