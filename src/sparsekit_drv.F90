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



MODULE sparsekit_drv
  USE ln_precision
  USE ln_allocation
  USE mat_def
  use skit_module

  private

  public :: dns2csr, csr2dns, coo2csr, csr2coo, csr2csc, csc2dns, dns2csc
  public :: clone, extract, concat, msort, mask
  public :: prealloc_sum, prealloc_mult
  public :: check_nnz, nzdrop, getdiag, trace, getelement
  public :: check_if_hermitian
  public :: zcsr2blk_sod

  public :: zsumcsr1
  public :: ramub_st
  public :: zprint_csrdns, zprint_csrcoo, ziluk_st, ztransp_st, ztransp2_st
  public :: zsubmat_st, zcopymat_st, zamub_st, zaplb_st, zcplsamub_st
  public :: zdagger, zspectral
  public :: zmask_realloc

  public :: nnz

  private :: zrconcatm_csr, zconcat_csr, zconcatm_csr
  private :: rcoocsr_st, rcsrcoo_st, rcsrdns_st, zdnscsr_st, zdnscsc_st, zcooxcsr_st
  private :: zcsrdns_st, zcscdns_st, zcoocsr_st, zcsrcoo_st, zcsrcsc_st, zcsccsr_st
  private :: rclone_st, zclone_st
  private :: rsumcsr, rsumcsrs, zsumcsr, zsumcsrs, zsumcsrs1s2, zsumdns, zsumdnss
  private :: zmultcsr, zmultcsrs
  private :: nnz_sum, nnz_sum1, nnz_mult

  !*************************************************************************
  !                                                                        |
  !Subroutines di chiamata compatta di Sparsekit e altre utilities         |
  !per la gestione delle matrici sparse                                    |
  !                                                                        |
  !*************************************************************************
  ! REAL ROUTINES:
  ! -----------------------------------------------------------------------
! !rcoocsr_st(coo,sp)           ::  converts coo to csr
! !rcsrcoo_st(sp,coo)           ::  converts csr to coo
! !rcsrdns_st(sp,dense)         ::  converts csr to dense
  !ramub_st(A,B,C)              ::  C = A*B  (C already allocated)
! !rsumcsr(A,B,C)               ::  preallocates C and performs C = A + B
! !rsumcsrs(A,B,s,C)            ::  preallocates C and performs C = A + s B
! !rzmask(A,B,C)                ::  mask A with B and output on C
! !rclone_st(sp_in,sp_out)      ::  creates and clone sp_in into sp_out
! !rcsort(A)                    ::  sort real matrix A in CSR format   !
  !
  ! COMPLEX ROUTINES:
  ! -----------------------------------------------------------------------
! !zclone_st(sp_in,sp_out)      ::  creates and clone the matrix sp_in
! !zdnscsr_st(dense,sp)         ::  converts dense to csr
! !zdnscsc_st(dense,sp)         ::  converts dense to csc
! !zcsrdns_st(sp,dense)         ::  converts csr to dense
! !zcscdns_st(sp,dense)         ::  converts csc to dense
! !zcoocsr_st(coo,sp)           ::  converts coo to csr
! !zcsrcoo_st(csr,coo)          ::  converts csr to coo
! !zcsrcsc_st(csr,csc)          ::  converts csr into csc
! !zcsccsr_st(csc,csr)          ::  converts csc into csr
! !zmultcsr(A,B,C)               :: preallocates C and performs C = A * B
! !zmultcsrs(A,B,s,C)            :: preallocates C and performs C = s A * B
! !zsumcsr(A,B,C)                :: preallocates C and performs C = A + B
! !zsumcsr1(A,B,C)               :: versione sorted di zsumcsr
! !zsumcsrs(A,B,s,C)             :: preallocates C and performs C = A + s B
! !zsumcsrs1s2(A,B,s1,s2,C)      :: preallocates C and performs C = s1 A + s2 B
! !zconcat_csr(A,B,row,col)      :: concat matrix B to A at coord (row,col)
! !zrconcatm_csr(A,s,B,'T',row,col)  :: concat in place Re(s B) or Im(s B) to real A
! !zconcatm_csr(A,s,B,row,col)   :: concat in place at coord (row,col)
! !zcsort_st(A)                  :: sort complex matrix A in CSR format
! !zmask(A,B,C)                  :: mask A with B and output on C
  !nnz=nnz(A,'+*',B)
  !nnz=nnz_sum(A,B)              :: computes nnz of C = A * B
  !nnz=nnz_mult(A,B)             :: computes nnz of C = A + B
  !nnz=nnz_sum1(A,B)             :: sorted version of nnz_sum
  !nnz=zcheck_nnz(A,r1,r2,c1,c2,sor):: checks nnz values of a sub-matrix in a csr
  !nnz=zchkdrp(A,drop)           :: drops values < drop from a dense matrix

  !ztransp_st(csr)              ::  compute transpose in place
  !zdagacsr(A,B)                 :: preallocates B = Hermitian conjugate of A
  !ztransp2_st(A_csr,B_csr)     ::  B = A^T

  !getdiag(A,D)                :: A CSR matrix, D vector of size nrow

  !zextract(A,r1,r2,c1,c2,S)     :: preallocates sub-block and extracts

  !zsubmat_st(csr,r1,r2,c1,c2,sub) :: extracts a subblock from csr
  !zcopymat_st(origin,dest)     ::  copy a csr into another
  !zprint_csrdns(io,csr)        ::  converts a csr in dense and prints it
  !zprint_csrcoo(io,csr)        ::  print a csr in coo format
  !zamub_st(A,B,C)              ::  C = A*B  (C already allocated)
  !zaplb_st(A,B,C)              ::  C = A+B  (C already allocated)
  !zcplsamub_st(A,B,s,C)        ::  C = C + s * A * B  (C already allocated)
  !ziluk_st()                   ::  incomplete LU factorization
  !zspectral(G1,G2,fl,A)         :: computes A = j(Gr - Ga)
  !zmask_realloc(G,M)            :: mask routine in place with reallocation
  !                                 (G masked by M)

  interface dns2csr
     module procedure rdnscsr_st
     module procedure zdnscsr_st
  end interface
  interface csr2dns
     module procedure rcsrdns_st
     module procedure zcsrdns_st
  end interface
  interface coo2csr
     module procedure rcoocsr_st
     module procedure zcoocsr_st
     module procedure zcooxcsr_st
  end interface
  interface csr2coo
     module procedure rcsrcoo_st
     module procedure zcsrcoo_st
  end interface
  interface csr2csc
     module procedure zcsrcsc_st
  end interface
  interface csc2csr
     module procedure zcsccsr_st
  end interface
  interface csc2dns
     module procedure zcscdns_st
  end interface
  interface dns2csc
     module procedure zdnscsc_st
  end interface
  interface clone
     module procedure rclone_st
     module procedure zclone_st
  end interface
  interface msort
     module procedure rcsort_st
     module procedure zcsort_st
  end interface
  interface mask
     !module procedure rmask
     module procedure zmask
     module procedure rzmask
     module procedure zmask_dns
  end interface
  interface prealloc_sum
     module procedure rsumcsr
     module procedure rsumcsrs
     module procedure zsumcsr
     module procedure zsumcsrs
     module procedure zsumcsrs1s2
     module procedure zsumdns
     module procedure zsumdnss
     module procedure zsumdnss1s2
  end interface
  interface prealloc_mult
     module procedure zmultcsr
     module procedure zmultcsrs
     module procedure zmultdns
     module procedure zmultdnss
     module procedure zmatmul
     module procedure zmatmuls
  end interface
  interface extract
     module procedure zextract_csr
     module procedure zextract_dns
  end interface
  interface concat
     module procedure  zconcat_csr
     module procedure  zrconcatm_csr
     module procedure  zconcatm_csr
  end interface
  interface nzdrop
     module procedure zchkdrp
  end interface
  interface check_nnz
     module procedure zcheck_nnz
  end interface
  interface getdiag
     module procedure getdiag_csr
     module procedure rgetdiag_csr
  end interface
  interface trace
     module procedure ztrace_csr
     module procedure ztrace_dns
     module procedure ztrace_arr
  end interface
  interface getelement
     module procedure rgetelm_csr
     module procedure zgetelm_csr
  end interface
  interface zdagger
     module procedure zdagacsr
     module procedure zdagadns
  end interface
  interface zspectral
     module procedure zspectral_csr
     module procedure zspectral_dns
  end interface
  interface check_if_hermitian
     module procedure check_hermitian_csr
     module procedure check_hermitian_dns
  end interface

  integer, parameter :: MISMATCH = 1
  integer, parameter :: CONVERR = 2
  integer, parameter :: BADINDEX = 3
  integer, parameter :: OUTOFBOUND = 4
  integer, parameter :: ALLOCERR = 5



CONTAINS
  ! -------------------------------------------------------------------------
  !   REAL ROUTINES:
  ! -------------------------------------------------------------------------
  !
  SUBROUTINE rcoocsr_st(coo,sp)

    !Short call for dsncsr:
    !convert matrix from COO to CSR format
    !Input:
    !coo: coordinate matrix to be converted
    !Output
    !sp: csr sparse matrix
    !

    IMPLICIT NONE

    TYPE(r_COO) :: coo
    TYPE(r_CSR) :: sp

    IF ((coo%nrow.NE.sp%nrow).OR.(coo%ncol.NE.sp%ncol)) then
        call error_msg('(rcoocsr_st)',MISMATCH)
    ENDIF

    IF (coo%nnz.EQ.0) THEN
       sp%rowpnt=1
    ELSE

       CALL coocsr(coo%nrow,sp%nnz,coo%nzval,coo%index_i,coo%index_j,sp%nzval, &
            sp%colind,sp%rowpnt)

    ENDIF

  END SUBROUTINE rcoocsr_st
  ! -------------------------------------------------------------------------
  SUBROUTINE rcsrcoo_st(sp,coo)

    !Short call for dsncsr:
    !convert matrix from CSR to COO format
    !Input:
    !sp: csr sparse matrix
    !Output
    !coo: coordinate matrix to be converted
    !

    IMPLICIT NONE

    TYPE(r_COO) :: coo
    TYPE(r_CSR) :: sp
    INTEGER :: ierr

    IF ((coo%nrow.NE.sp%nrow).OR.(coo%ncol.NE.sp%ncol)) THEN
       call error_msg('(rcsrcoo_st)',MISMATCH)
    ENDIF

    IF (sp%nnz.NE.0) THEN

       CALL csrcoo(coo%nrow,3,sp%nnz,sp%nzval,sp%colind,sp%rowpnt,coo%nnz,coo%nzval, &
             & coo%index_i,coo%index_j,ierr)
    ENDIF

  END SUBROUTINE rcsrcoo_st
  ! -------------------------------------------------------------------------
  SUBROUTINE rcsrdns_st(sp,dense)

    !Short call for csrdns
    !Input:
    !sp: sparse matrix to be converted
    !Output:
    !dense: dense matrix
    !NOTE: for squared matrix
    IMPLICIT NONE

    type(r_CSR) :: sp
    type(r_DNS) :: dense

    integer :: ierr

    if ((dense%nrow.ne.sp%nrow).or.(dense%ncol.ne.sp%ncol)) then
        call error_msg('(rcsrdns_st)',MISMATCH)
    endif

    dense%val(:,:)=0.0_dp

    IF (sp%nnz.NE.0) THEN

       CALL csrdns(sp%nrow,sp%ncol,sp%nzval,sp%colind,sp%rowpnt,dense%val,ierr)

    ENDIF

    if (ierr.ne.0) call error_msg('(rcsrdns_st)',CONVERR)

  END SUBROUTINE rcsrdns_st
  ! -------------------------------------------------------------------------
  SUBROUTINE rdnscsr_st(dense,sp)

    !Short call for csrdns
    !Input:
    !sp: sparse matrix to be converted
    !Output:
    !dense: dense matrix
    !NOTE: for squared matrix
    IMPLICIT NONE

    type(r_CSR) :: sp
    type(r_DNS) :: dense

    integer :: ierr

    if ((dense%nrow.ne.sp%nrow).or.(dense%ncol.ne.sp%ncol)) then
        call error_msg('(rdnscsr_st)',MISMATCH)
    endif

    IF (sp%nnz.NE.0) THEN
       call dnscsr(dense%nrow,dense%ncol,sp%nnz,dense%val,&
            sp%nzval,sp%colind,sp%rowpnt,ierr)

       sp%nnz=sp%rowpnt(sp%nrow+1)-1

    ELSE

       sp%rowpnt=1

    ENDIF

    if (ierr.ne.0) call error_msg('(rdnscsr_st)',CONVERR)

  END SUBROUTINE rdnscsr_st
  ! -------------------------------------------------------------------------

  SUBROUTINE ramub_st(A_csr,B_csr,C_csr)

    !Short call for amub:
    !A-csr+B_csr without preallocation
    !Input:
    !A_csr: CSR sparse matrix
    !B_csr: CSR sparse matrix
    !Output:
    !C_csr: CSR sparse matrix (allocated externally)

    implicit none

    type(r_CSR) :: A_csr,B_csr,C_csr
    integer :: ierr,B_ncol
    integer, DIMENSION(:), ALLOCATABLE :: iw

    if(C_csr%nrow.ne.A_csr%nrow .or. C_csr%ncol.ne.B_csr%ncol) THEN
       call error_msg('(ramub_st)',MISMATCH)
    endif

    B_ncol=B_csr%ncol
    C_csr%ncol=B_ncol
    !Allocazione work array iw
    call log_allocate(iw,B_ncol)

    call amub(A_csr%nrow,B_ncol,1,A_csr%nzval,A_csr%colind,A_csr%rowpnt,&
              B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind,&
              C_csr%rowpnt,C_csr%nnz,iw,ierr)

    if (ierr.ne.0) call error_msg('(rmaub_st)',OUTOFBOUND)

    call log_deallocate(iw)

  end subroutine ramub_st

  !***********************************************************************
  !
  !  Subroutine di somma sparsa compatta reale (CSR unsorted)
  !
  !***********************************************************************

  SUBROUTINE rsumcsr(A_csr,B_csr,C_csr)

    IMPLICIT NONE

    TYPE(r_CSR) :: A_csr,B_csr,C_csr
    INTEGER, DIMENSION(:), ALLOCATABLE :: iw
    INTEGER :: ierr,A_ncol

    IF(A_csr%nrow.NE.B_csr%nrow .or. A_csr%ncol.NE.B_csr%ncol) THEN
         call error_msg('(rsumcsr)',MISMATCH)
    endif

    IF ((A_csr%nnz.EQ.0).AND.(B_csr%nnz.EQ.0)) THEN
       CALL  create(C_csr,A_csr%nrow,A_csr%ncol,0)
       C_csr%rowpnt=1
    ELSEIF (A_csr%nnz.EQ.0) THEN
       CALL create(C_csr,B_csr%nrow,B_csr%ncol,B_csr%nnz)
       C_csr%nzval=B_csr%nzval; C_csr%colind=B_csr%colind; C_csr%rowpnt=B_csr%rowpnt;
    ELSEIF (B_csr%nnz.EQ.0) THEN
       CALL create(C_csr,A_csr%nrow,A_csr%ncol,A_csr%nnz)
       C_csr%nzval=A_csr%nzval; C_csr%colind=A_csr%colind; C_csr%rowpnt=A_csr%rowpnt;

    ELSE

       A_ncol = A_csr%ncol

       !Alloca le parti di C_csr di interesse
       C_csr%nrow=A_csr%nrow
       C_csr%ncol=A_csr%ncol
       C_csr%nnz=(A_csr%nnz+B_csr%nnz)

       call log_allocate(C_csr%nzval,MISMATCH)
       call log_allocate(C_csr%colind,C_csr%nnz)
       call log_allocate(C_csr%rowpnt,(C_csr%nrow+1))

       !Allocazione work array iw
       call log_allocate(iw,A_ncol)

       call aplb (A_csr%nrow,A_ncol,0,A_csr%nzval,A_csr%colind,A_csr%rowpnt,&
            B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind,&
            C_csr%rowpnt,C_csr%nnz,iw,ierr)

       if (ierr.ne.0) call error_msg('(rsumcsr)',OUTOFBOUND)

       C_csr%nnz=C_csr%rowpnt(A_csr%nrow+1)-1

       call log_deallocate(C_csr%nzval)
       call log_deallocate(C_csr%colind)

       call log_allocate(C_csr%nzval,C_csr%nnz)
       call log_allocate(C_csr%colind,C_csr%nnz)

       call aplb (A_csr%nrow,A_ncol,1,A_csr%nzval,A_csr%colind,A_csr%rowpnt,&
            B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind,&
            C_csr%rowpnt,C_csr%nnz,iw,ierr)

       if (ierr.ne.0) call error_msg('(rsumcsr)',OUTOFBOUND)
       call log_deallocate(iw)

    ENDIF

  end subroutine rsumcsr

  !---------------------------------------------------------------------

  SUBROUTINE rclone_st(sp_in,sp_out)

    IMPLICIT NONE

    TYPE(r_CSR) :: sp_in,sp_out
    INTEGER :: i

    CALL create(sp_out,sp_in%nrow,sp_in%ncol,sp_in%nnz)

    IF (sp_in%nnz.NE.0) THEN

       DO i=1,sp_in%nnz
          sp_out%nzval(i) = sp_in%nzval(i)
          sp_out%colind(i) = sp_in%colind(i)
       ENDDO

    ENDIF

    DO i=1,sp_in%nrow+1
       sp_out%rowpnt(i) = sp_in%rowpnt(i)
    ENDDO

  END SUBROUTINE rclone_st

  !***********************************************************************
  !
  !  Subroutine di somma sparsa compatta con prodotto per scalare reale
  !
  !***********************************************************************

  subroutine rsumcsrs(A_csr,B_csr,s,C_csr)


    implicit none

    type(r_CSR) :: A_csr,B_csr,C_csr
    real(kind=dp) :: s
    integer, DIMENSION(:), ALLOCATABLE :: iw
    integer :: ierr,A_ncol

    ierr=0
    if(A_csr%nrow.ne.B_csr%nrow .or. A_csr%ncol.ne.B_csr%ncol) then
            call error_msg('(rsumcsrs)',MISMATCH)
    endif

    IF ((A_csr%nnz.EQ.0).AND.(B_csr%nnz.EQ.0)) THEN
       CALL  create(C_csr,A_csr%nrow,A_csr%ncol,0)
       C_csr%rowpnt=1
    ELSEIF (A_csr%nnz.EQ.0) THEN
       CALL create(C_csr,B_csr%nrow,B_csr%ncol,B_csr%nnz)
       C_csr%nzval=B_csr%nzval; C_csr%colind=B_csr%colind; C_csr%rowpnt=B_csr%rowpnt;
    ELSEIF (B_csr%nnz.EQ.0) THEN
       CALL create(C_csr,A_csr%nrow,A_csr%ncol,A_csr%nnz)
       C_csr%nzval=A_csr%nzval; C_csr%colind=A_csr%colind; C_csr%rowpnt=A_csr%rowpnt;

    ELSE

    C_csr%nrow=A_csr%nrow
    C_csr%ncol=A_csr%ncol

    A_ncol=A_csr%ncol
    !Alloca le parti di C_csr di interesse
    C_csr%nnz=(A_csr%nnz+B_csr%nnz)

    call log_allocate(C_csr%nzval,MISMATCH)
    call log_allocate(C_csr%colind,C_csr%nnz)
    call log_allocate(C_csr%rowpnt,(C_csr%nrow+1))

    call log_allocate(iw,A_ncol)

    call aplb(A_csr%nrow,A_ncol,0,A_csr%nzval,A_csr%colind,A_csr%rowpnt,&
              B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind,&
              C_csr%rowpnt,C_csr%nnz,iw,ierr)

    if (ierr.ne.0) call error_msg('(rsumcsrs)',OUTOFBOUND)

    C_csr%nnz=C_csr%rowpnt(A_csr%nrow+1)-1

    call log_deallocate(C_csr%nzval)
    call log_deallocate(C_csr%colind)

    call log_allocate(C_csr%nzval,C_csr%nnz)
    call log_allocate(C_csr%colind,C_csr%nnz)

    call aplsb(A_csr%nrow,A_ncol,1,A_csr%nzval,A_csr%colind,A_csr%rowpnt,s,&
               B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind,&
               C_csr%rowpnt,C_csr%nnz,iw,ierr)

    if (ierr.ne.0) call error_msg('(rsumcsrs)',OUTOFBOUND)

        call log_deallocate(iw)

    ENDIF

  end subroutine rsumcsrs

  !*************************************************************************
  !
  !  Subroutine per il sorting delle matrici csr
  !
  !*************************************************************************

  subroutine rcsort_st(A)

  implicit none

  TYPE(r_CSR) :: A
  INTEGER, DIMENSION(:), ALLOCATABLE :: iwork

   if (.not.A%sorted) then
     CALL log_allocate(iwork,MAX(A%nrow+1,2*A%nnz))
     CALL csort(A%nrow,A%nzval,A%colind,A%rowpnt,iwork,.true.)
     call log_deallocate(iwork)
     A%sorted=.true.
   endif

  end subroutine rcsort_st

  !*************************************************************************
  !
  !  Subroutine per il masking delle matrici csr
  !  (reale mascherata da complessa)
  !  (utilizzabile anche in place)
  !
  !*************************************************************************

  SUBROUTINE rzmask(A,M,C)

  implicit none

  TYPE(z_CSR) :: M
  TYPE(r_CSR) :: A,C
  INTEGER :: ierr
  logical, ALLOCATABLE, DIMENSION(:) :: iw

  allocate(iw(A%ncol),stat=ierr)
  IF (ierr.ne.0) call error_msg('(rzmask)',ALLOCERR)

  CALL amask(A%nrow,A%ncol,A%nzval,A%colind,A%rowpnt,M%colind,&
              M%rowpnt,C%nzval,C%colind,C%rowpnt,iw,C%nnz,ierr)

  IF (ierr.GT.1) call error_msg('(rzmask)',CONVERR)

  deallocate(iw)

  END SUBROUTINE rzmask

  ! -------------------------------------------------------------------------
  !   COMPLEX ROUTINES:
  ! -------------------------------------------------------------------------
  !


  !*************************************************************************
  !
  !  Subroutine per il masking delle matrici csr
  !  (reale mascherata da complessa)
  !  (utilizzabile anche in place)
  !
  !*************************************************************************

  SUBROUTINE zmask(A,M,C)

  implicit none

  TYPE(z_CSR) :: M
  TYPE(z_CSR) :: A,C
  INTEGER :: ierr
  logical, ALLOCATABLE, DIMENSION(:) :: iw

  allocate(iw(A%ncol),stat=ierr)
  IF (ierr.ne.0) call error_msg('(zmask)',ALLOCERR)

  CALL amask(A%nrow,A%ncol,A%nzval,A%colind,A%rowpnt,M%colind,&
              M%rowpnt,C%nzval,C%colind,C%rowpnt,iw,C%nnz,ierr)

  IF (ierr.GT.1) call error_msg('(zmask)',CONVERR)

  deallocate(iw)

  END SUBROUTINE zmask

  !*************************************************************************
  ! Mask of dense matrix (set elements to 0.0)
  !*************************************************************************
  SUBROUTINE zmask_dns(A,M)

  implicit none

  TYPE(z_DNS) :: M
  TYPE(z_DNS) :: A

  INTEGER :: i, k

   if(A%nrow .ne. M%nrow .and. A%ncol.ne.M%ncol) then
       call error_msg('(zmask_dns)',MISMATCH)
   endif

   do i = 1, A%ncol
     do k = 1, A%nrow
       if(abs(M%val(k,i)).eq.0.0_dp) A%val(k,i)=(0.0_dp,0.0_dp)
     enddo
   enddo

  END SUBROUTINE zmask_dns

  !**********************************************************************
  !
  !  Masking routine with reallocation
  !
  !**********************************************************************

  SUBROUTINE zmask_realloc(G,M)

    !**********************************************************************
    !Input/Output:
    !G:matrix to be masked and reallocated
    !M: masking matrix
    !
    ! It differs by zamask because it acts as a masking routine in
    ! place, but reallocating G with the right size
    ! (avoiding memory waste)
    !**********************************************************************

    IMPLICIT NONE

    TYPE(z_CSR) :: G,M
    TYPE(z_CSR) :: work

    CALL create(work, M%nrow, M%ncol, M%nnz)
    CALL zmask(G,M,work)
    CALL destroy(G)
    CALL zclone_st(work, G)
    CALL destroy(work)

  END SUBROUTINE zmask_realloc


!-----------------------------------------------------------------------------------

  SUBROUTINE zclone_st(sp_in,sp_out)

    IMPLICIT NONE

    TYPE(z_CSR) :: sp_in,sp_out
    INTEGER :: i

    CALL create(sp_out,sp_in%nrow,sp_in%ncol,sp_in%nnz)

    IF (sp_in%nnz.NE.0) THEN

       DO i=1,sp_in%nnz
          sp_out%nzval(i) = sp_in%nzval(i)
          sp_out%colind(i) = sp_in%colind(i)
       ENDDO

    ENDIF

    DO i=1,sp_in%nrow+1
       sp_out%rowpnt(i) = sp_in%rowpnt(i)
    ENDDO

  END SUBROUTINE zclone_st


  !----------------------------------------------------------------------------
  SUBROUTINE zdnscsr_st(dense,sp)

    !Short call for dsncsr:
    !convert matrix from dense to CSR format
    !Input:
    !dense: dense matrix to be converted
    !nrow: row number of dense matrix (input)
    !Output
    !sp: csr sparse matrix
    !
    !NOTE: works on square matrix

    IMPLICIT NONE

    TYPE(z_DNS) :: dense
    TYPE(z_CSR) :: sp
    INTEGER :: ierr

    IF(dense%nrow.ne.sp%nrow) THEN
        call error_msg('(zdnscsr_st)',MISMATCH)
    ENDIF

    IF (sp%nnz.NE.0) THEN

       call dnscsr(dense%nrow,dense%ncol,sp%nnz,dense%val,sp%nzval,sp%colind,sp%rowpnt,ierr)

       sp%nnz=sp%rowpnt(sp%nrow+1)-1

    ELSE

       sp%rowpnt=1

    ENDIF

  END SUBROUTINE zdnscsr_st

  ! -------------------------------------------------------------------------

  SUBROUTINE zdnscsc_st(dense,sp)

    !convert matrix from dense to CSR format
    !Input:
    !dense: dense matrix to be converted
    !nrow: row number of dense matrix (input)
    !Output
    !sp: csr sparse matrix
    !
    !NOTE: works on square matrix

    IMPLICIT NONE

    TYPE(z_DNS) :: dense
    TYPE(z_CSC) :: sp
    TYPE(z_CSR) :: sp_csr
    INTEGER :: ierr

    IF (sp%nnz.NE.0) THEN

       call create(sp_csr,dense%nrow,dense%ncol,sp%nnz)
       call dnscsr(dense%nrow,dense%ncol,sp_csr%nnz,dense%val,sp_csr%nzval,sp_csr%colind,sp_csr%rowpnt,ierr)
       call zcsrcsc_st(sp_csr,sp)
       call destroy(sp_csr)

    ELSE

       sp%colpnt=1

    ENDIF

  END SUBROUTINE zdnscsc_st

  !---------------------------------------------------------------------------

  SUBROUTINE zcsrdns_st(sp,dense)

    !Short call for csrdns
    !Input:
    !sp: sparse matrix to be converted
    !Output:
    !dense: dense matrix
    !NOTE: for squared matrix
    IMPLICIT NONE

    type(z_CSR) :: sp
    type(z_DNS) :: dense

    integer :: ierr

    if ((dense%nrow.ne.sp%nrow).or.(dense%ncol.ne.sp%ncol)) then
        call error_msg('(rcsrdns_st)',MISMATCH)
    endif

    dense%val(:,:)=(0.0_dp, 0.0_dp)
    ierr=0
    IF (sp%nnz.NE.0) THEN

       CALL csrdns(sp%nrow,sp%ncol,sp%nzval,sp%colind,sp%rowpnt,dense%val,ierr)

    ENDIF

    if (ierr.ne.0) call error_msg('(zcsrdns_st)',CONVERR)

  END SUBROUTINE zcsrdns_st

  !--------------------------------------------------------------------------

  subroutine zcscdns_st(sp,dense)

    !Short call for csrdns
    !Input:
    !sp: sparse matrix to be converted
    !Output:
    !dense: dense matrix
    !NOTE: for squared matrix
    implicit none

    type(z_CSC) :: sp
    type(z_CSR) :: sp_csr
    type(z_DNS) :: dense

    integer :: ierr

    if ((dense%nrow.ne.sp%nrow).or.(dense%ncol.ne.sp%ncol)) then
        call error_msg('(rcoocsr_st)',MISMATCH)
    endif

    IF (sp%nnz.NE.0) THEN

       call create(sp_csr,sp%nrow,sp%ncol,sp%nnz)

       call zcsccsr_st(sp, sp_csr)
       call csrdns(sp_csr%nrow,sp_csr%ncol,sp_csr%nzval,sp_csr%colind,sp_csr%rowpnt,dense%val,ierr)
       call destroy(sp_csr)

    else

       dense%val(:,:)=(0.0_dp, 0.0_dp)

    endif

    if (ierr.ne.0) call error_msg('(zcscdns_st)',CONVERR)

  end subroutine zcscdns_st

  !--------------------------------------------------------------------------
  subroutine zcoocsr_st(coo,sp)

    !Short call for dsncsr:
    !convert matrix from COO to CSR format
    !Input:
    !coo: coordinate matrix to be converted
    !Output
    !sp: csr sparse matrix
    !

    implicit none

    type(z_COO) :: coo
    type(z_CSR) :: sp
    integer, ALLOCATABLE, DIMENSION(:) :: iwork
    logical :: values
    integer :: ierr

    if ((coo%nrow.ne.sp%nrow).or.(coo%ncol.ne.sp%ncol)) then
        call error_msg('(rcoocsr_st)',MISMATCH)
    endif

    if (coo%nnz.ne.0) then

       call coocsr(coo%nrow,sp%nnz,coo%nzval,coo%index_i,coo%index_j,sp%nzval,sp%colind,sp%rowpnt)

       allocate(iwork(2*sp%nnz), stat=ierr)
       if (ierr.ne.0) call error_msg('(zcoocsr_st)',ALLOCERR)
       values= .true.

       call csort(sp%nrow,sp%nzval,sp%colind,sp%rowpnt,iwork,values)

       sp%sorted = .true.
       deallocate(iwork)

    else

       sp%rowpnt=1

    endif


  end subroutine zcoocsr_st
  ! -------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  subroutine zcsrcoo_st(csr,coo)

    !Short call for zcsrcoo:
    !convert matrix from CSR to COO format

    implicit none

    type(z_COO) :: coo
    type(z_CSR) :: csr
    integer :: ierr

    if ((coo%nrow.ne.csr%nrow).or.(coo%ncol.ne.csr%ncol)) then
        call error_msg('(zcsrcoo_st)',MISMATCH)
    endif

    if (csr%nnz.ne.0) then

       call csrcoo(csr%nrow,3,coo%nnz,csr%nzval,csr%colind,csr%rowpnt,&
            csr%nnz,coo%nzval,coo%index_i,coo%index_j,ierr)

       if(ierr.ne.0) call error_msg('(zcsrcoo_st)',CONVERR)

    endif

  end subroutine zcsrcoo_st

  !---------------------------------------------------------------------------
  subroutine zcsrcsc_st(A_csr, A_csc)

    !Short call for csrcsc (conversion CSR-CSC, transposition)
    !Input:
    !A_csr: matrix in CSR format
    !A_csc: matrix in CSC format

    implicit none

    type(z_CSR) :: A_csr
    type(z_CSC) :: A_csc
    integer :: job, ipos

    if ((A_csr%nrow.ne.A_csc%nrow).or.(A_csr%ncol.ne.A_csc%ncol)) then
        call error_msg('(rcsrcsc_st)',MISMATCH)
    endif

    IF (A_csr%nnz.NE.0) THEN

       job=1
       ipos=1
       call csrcsc2(A_csr%nrow, A_csr%ncol, job, ipos, A_csr%nzval, A_csr%colind, &
            A_csr%rowpnt, A_csc%nzval, A_csc%rowind, A_csc%colpnt)

    ELSE

       A_csc%colpnt=1

    ENDIF

  end subroutine zcsrcsc_st

  !---------------------------------------------------------------------------

  subroutine zcsccsr_st(A_csc, A_csr)

    !Short call for csrcsc (conversion CSC-CSR, transposition)
    !Input:
    !A_csr: matrix in CSR format
    !Output:
    !A_csc: matrix in CSC format
    !Works on square matrixes

    implicit none

    type(z_CSR) :: A_csr
    type(z_CSC) :: A_csc
    integer :: job, ipos

    if ((A_csr%nrow.ne.A_csc%nrow).or.(A_csr%ncol.ne.A_csc%ncol)) then
        call error_msg('(zcsrcsc_st)',MISMATCH)
    endif

    IF (A_csr%nnz.NE.0) THEN

       job=1
       ipos=1
       call csrcsc(A_csc%nrow, job, ipos, A_csc%nzval, A_csc%rowind, A_csc%colpnt, &
            A_csr%nzval, A_csr%colind, A_csr%rowpnt)

    ELSE

       A_csr%rowpnt=1

    ENDIF

  end subroutine zcsccsr_st
  !---------------------------------------------------------------------------

  subroutine zprint_csrdns(iofile,A_csr,flag)

    type(z_CSR) :: A_csr
    integer :: iofile
    character(1) :: flag

    type(z_DNS) :: A
    integer :: i,j

    call create(A,A_csr%nrow,A_csr%ncol)

    call zcsrdns_st(A_csr,A)

    do j=1,A%ncol
       do i=1,A%nrow-1
             if (flag.eq.'r') then
                write(iofile,'(f10.5)', advance='NO') real( A%val(i,j) )
             elseif (flag.eq.'i') then
                write(iofile,'(f10.5)', advance='NO') aimag( A%val(i,j) )
             elseif (flag.eq.'c') then
                write(iofile,'(2E23.15)', advance='NO')      ( A%val(i,j) )
             endif
       enddo
             if (flag.eq.'r') then
                write(iofile,'(f10.5)', advance='YES') real( A%val(i,A%ncol) )
             elseif (flag.eq.'i') then
                write(iofile,'(f10.5)', advance='YES') aimag( A%val(i,j) )
             elseif (flag.eq.'c') then
                write(iofile,'(2E23.15)', advance='YES')      ( A%val(i,j) )
             endif
    enddo

    call destroy(A)

  end subroutine zprint_CSRDNS


  !subroutine rprint_csrdns(iofile,A_csr,flag)
  !
  !  type(r_CSR) :: A_csr
  !  integer :: i,j,A_ncol,ierr,iofile
  !  character(1) :: flag
  !
  !  type(r_DNS) :: A
  !
  !  call create(A,A_csr%nrow,A_csr%ncol)
  !
  !  call rcsrdns_st(A_csr,A)
  !
  !  do i=1,A%nrow
  !     do j=1,A%ncol-1
  !           if (flag.eq.'r') then
  !              write(iofile,'(f10.5)', advance='NO') dreal( A%val(i,j) )
  !           elseif (flag.eq.'i') then
  !              write(iofile,'(f10.5)', advance='NO') dimag( A%val(i,j) )
  !           elseif (flag.eq.'c') then
  !              write(iofile,'(E23.15)', advance='NO')      ( A%val(i,j) )
  !           endif
  !     enddo
  !           if (flag.eq.'r') then
  !              write(iofile,'(f10.5)', advance='YES') real( A%val(i,A%ncol) )
  !           elseif (flag.eq.'i') then
  !              write(iofile,'(f10.5)', advance='YES') dimag( A%val(i,j) )
  !           elseif (flag.eq.'c') then
  !              write(iofile,'(E23.15)', advance='YES')      ( A%val(i,j) )
  !           endif
  !  enddo
  !
  !  call destroy(A)
  !
  !end subroutine rprint_csrdns




  !---------------------------------------------------------------------------
  subroutine zprint_csrcoo(iofile,A_csr,flag)

    type(z_CSR) :: A_csr
    integer :: iofile
    character(1) :: flag

    type(z_COO) :: A
    integer :: k

    call create(A,A_csr%nrow,A_csr%ncol,A_csr%nnz)

    A%nzval=(0.0_dp,0.0_dp)

    call csr2coo(A_csr,A)

    !write (iofile,*) '# matrix dimension = ', A%nrow, ' x ', A%ncol
    !write (iofile,*) '#'

    do k = 1, A%nnz

       if (flag.eq.'r') then
          write(iofile,'(2i8,f20.10)') A%index_i(k), A%index_j(k), real(A%nzval(k))
       elseif (flag.eq.'i') then
          write(iofile,'(2i8,f20.10)') A%index_i(k), A%index_j(k), aimag(A%nzval(k))
       elseif (flag.eq.'c') then
          write(iofile,'(2i8,(f20.10,f20.10))')  A%index_i(k), A%index_j(k), A%nzval(k)
       endif

    enddo

    call destroy(A)

  end subroutine zprint_csrcoo

  !---------------------------------------------------------------------------

  subroutine ziluk_st(A_csr, lfil, LU_msr, LU_ju, LU_levs, LU_iwk)

    !Short call for iluk:
    !LU preconditioning for PGMRES solver
    !Input:
    !A_csr: matrix to be preconditioned
    !lfil: order of factorization (min 0)
    !Output:
    !LU_msr: LU matrix stored in MSR format
    !LU_ju: integer array of lenght n containing pointers to the begining
    !of each row of U in LU MSR matrix
    !Work:
    !LU_levs: integer array (of minimum lenght of LU_msr%nzval) containing
    !the level of each element in LU_msr
    !LU_iwk: integer. Minimum lenght of LU_msr
    !NOTE: work arrays must be external as depend on lfil

    implicit none

    integer :: LU_iwk, lfil
    type(z_CSR) :: A_csr
    type(z_MSR) :: LU_msr
    integer, DIMENSION(LU_iwk) :: LU_levs
    integer, DIMENSION(A_csr%nrow) :: LU_ju
    integer, DIMENSION(:), ALLOCATABLE :: jw
    complex(kind=dp), DIMENSION(:), ALLOCATABLE :: w
    integer :: ierr

    call log_allocate(jw,3*A_csr%nrow)
    call log_allocate(w, A_csr%nrow)

    !call ziluk(A_csr%nrow, A_csr%nzval, A_csr%colind, A_csr%rowpnt, lfil, &
    !     LU_msr%nzval, LU_msr%index, LU_ju, LU_levs, LU_iwk, w, jw, ierr)

    if (ierr.ne.0) then
       write(*,*)'ERROR: LU FACTORIZATION FAILED'
       write(*,*)'ERROR NUMBER ',ierr,' IN ILUK SUBROUTINE'
       if (ierr.gt.0) then
         write(*,*) 'Zero pivot encountered at step number ',ierr
       endif
       if (ierr.eq.(-1)) then
         write(*,*) 'Input matrix may be wrong'
       endif
       if (ierr.eq.(-2)) then
         write(*,*) 'Matrix L overflows array al'
       endif
       if (ierr.eq.(-3)) then
         write(*,*) 'Matrix U overflows array alu'
       endif
       if (ierr.eq.(-4)) then
         write(*,*) 'Illegal value for lfil'
       endif
       if (ierr.eq.(-5)) then
         write(*,*) 'zero row encounterd in A or U'
       endif
       STOP
    endif

  call log_deallocate(jw)
  call log_deallocate(w)

  end subroutine ziluk_st

  !----------------------------------------------------------------------------

  subroutine ztransp_st(A_csr)

    !Short call for transp (matrix in-place transposition)
    !Input:
    !A_csr: matrix in CSR format
    !nrow: number of rows in A_csr
    !ncol: number of columns in A_csr

    implicit none

    type(z_CSR) :: A_csr
    integer :: nrow,ncol,ierr
    integer, DIMENSION(:), ALLOCATABLE :: iwk

    nrow=A_csr%nrow
    ncol=A_csr%ncol

    if (A_csr%nnz.eq.0) then
       A_csr%nrow=ncol
       A_csr%ncol=nrow
    else

       call log_allocate(iwk,A_csr%nnz)
       call transp(nrow,ncol,A_csr%nzval,A_csr%colind,A_csr%rowpnt,iwk,ierr)

       call log_deallocate(iwk)
       if (ierr.ne.0) call error_msg('(ztransp_st)',CONVERR)

    endif

  end subroutine ztransp_st
!------------------------------------------------------------------------

  subroutine ztransp2_st(A_csr,B_csr)

    !Short call for transp (matrix transposition)
    !Input:
    !A_csr: matrix in CSR format
    !ncol: number of columns in A_csr
    !B_csr: output matrix, transposed of A

    implicit none

    type(z_CSR) :: A_csr, B_csr
    integer :: ncol

    ncol=A_csr%ncol

    call create(B_csr,ncol,A_csr%nrow,A_csr%nnz)

    IF (B_csr%nnz.EQ.0) THEN
       B_csr%nrow=A_csr%ncol
       B_csr%ncol=A_csr%nrow
       B_csr%rowpnt=1
    ELSE

    call csrcsc2(A_csr%nrow,ncol,1,1,A_csr%nzval,A_csr%colind,A_csr%rowpnt,&
                     B_csr%nzval,B_csr%colind,B_csr%rowpnt)

    ENDIF

  end subroutine ztransp2_st
  !------------------------------------------------------------------------------

  subroutine zsubmat_st(A_csr,i1,i2,j1,j2,Asub_csr)

    !Short call for submat (extraction of submatrix in csr format from matrix in csr format)
    !Input:
    !A_csr: complete matrix
    !i1,i2: row to be extracted
    !j1,j2: columns to be extracted
    !Output:
    !Asub_csr: extracted submatrix

    implicit none

    type(z_CSR) :: A_csr,Asub_csr
    integer :: i1,i2,j1,j2,nr,nc
    integer, PARAMETER :: job=1

    nr=i2-i1+1
    nc=j2-j1+1

    IF ((Asub_csr%nrow.NE.nr).OR.(Asub_csr%ncol.NE.nc)) THEN
       call error_msg('(zsubmat_st)',MISMATCH)
    ENDIF

    IF (A_csr%nnz.NE.0) THEN

       call submat(A_csr%nrow,job,i1,i2,j1,j2,A_csr%nzval,A_csr%colind,A_csr%rowpnt,nr,nc,&
            Asub_csr%nzval,Asub_csr%colind,Asub_csr%rowpnt)

    ELSE

       Asub_csr%rowpnt=1

    ENDIF

  end subroutine zsubmat_st

  !----------------------------------------------------------------------------

  subroutine zcopymat_st(A_csr,B_csr)

    !Short call for copmat: copy matrix A_csr in matrix B_csr

    implicit none

    type(z_CSR) :: A_csr,B_csr

    IF ((A_csr%nrow.NE.B_csr%nrow).OR.(A_csr%ncol.NE.B_csr%ncol).OR.(A_csr%nnz.NE.B_csr%nnz)) THEN
        call error_msg('(zcopmat_st)',MISMATCH)
    ENDIF

    IF (A_csr%nnz.NE.0) THEN

       call copymat(A_csr%nrow,A_csr%nzval,A_csr%colind,A_csr%rowpnt,B_csr%nzval,&
            B_csr%colind,B_csr%rowpnt,1,MISMATCH)

    ELSE

       B_csr%rowpnt=1

    ENDIF

  end subroutine zcopymat_st


  !****************************************************************************
  !Utilities varie
  !****************************************************************************
  subroutine zamub_st(A_csr,B_csr,C_csr)

    !*****************************************************************
    !
    !Input:
    !A_csr: primo fattore in formato CSR
    !B_csr: secondo fattore in formato  CSR
    !C_csr: risultato in formato CSR (l'allocazione esatta viene eseguita
    !nella subroutine
    !
    !****************************************************************
    implicit none

    type(z_CSR) :: A_csr,B_csr,C_csr
    integer :: ierr,B_ncol
    integer, DIMENSION(:), ALLOCATABLE :: iw


    if(C_csr%nrow.ne.A_csr%nrow .or. C_csr%ncol.ne.B_csr%ncol) THEN
            call error_msg('(zamub_st)',MISMATCH)
    endif

    !Allocazione work array iw

    B_ncol=B_csr%ncol
    call log_allocate(iw,B_ncol)

    call amub(A_csr%nrow,B_ncol,1,A_csr%nzval,A_csr%colind,A_csr%rowpnt,&
         B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind,C_csr%rowpnt,&
         C_csr%nnz,iw,ierr)

    if (ierr.ne.0) call error_msg('(zamub_st)',OUTOFBOUND)

    call log_deallocate(iw)

  end subroutine zamub_st
  !****************************************************************************
  subroutine zaplb_st(A_csr,B_csr,C_csr)

    !*****************************************************************
    !
    !Input:
    !A_csr: primo fattore in formato CSR
    !B_csr: secondo fattore in formato  CSR
    !B_ncol: numero di colonne in A_csr e B_csr
    !C_csr: risultato in formato CSR (l'allocazione esatta viene eseguita
    !nella subroutine
    !
    !****************************************************************
    implicit none

    type(z_CSR) :: A_csr,B_csr,C_csr
    integer :: ierr,A_ncol
    integer, DIMENSION(:), ALLOCATABLE :: iw


    if(A_csr%nrow.ne.B_csr%nrow .or. A_csr%ncol.ne.B_csr%ncol) then
            call error_msg('(zaplb_st)',MISMATCH)
    endif

    !Allocazione work array iw
    A_ncol=A_csr%ncol
    C_csr%nrow=A_csr%nrow
    C_csr%ncol=A_csr%ncol

    call log_allocate(iw,A_ncol)

    call aplb (A_csr%nrow,A_ncol,1,A_csr%nzval,A_csr%colind,A_csr%rowpnt,&
                B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind,&
                C_csr%rowpnt,C_csr%nnz,iw,ierr)

    if (ierr.ne.0) call error_msg('(zaplb_st)',OUTOFBOUND)

    call log_deallocate(iw)

  end subroutine zaplb_st
  !****************************************************************************

  !-------------------------------------------------------------------
  ! performs the matrix by matrix product C = C + s A B
  !-------------------------------------------------------------------
  subroutine zcplsamub_st(A_csr,B_csr,s,C_csr)

    !*****************************************************************
    !
    !Input:
    !A_csr: primo fattore in formato CSR
    !B_csr: secondo fattore in formato  CSR
    !  s  : scalare
    !C_csr: risultato in formato CSR (l'allocazione esatta viene eseguita
    !nella subroutine
    !
    !****************************************************************

    implicit none

    type(z_CSR) :: A_csr,B_csr,C_csr
    complex(kind=dp) :: s
    integer :: ierr,B_ncol
    integer, DIMENSION(:), ALLOCATABLE :: iw

    !Allocazione work array iw

    B_ncol=B_csr%ncol
    call log_allocate(iw,B_ncol)

     ! zcplsamub(nrow,ncol,job,a,ja,ia,s,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
     ! complex*16 a(*), b(*), c(*), s
     ! integer ja(*),jb(*),jc(*),ia(nrow+1),ib(*),ic(*),iw(ncol)
     ! integer job : 0 = no comp , 1 = comp

    call cplsamub(A_csr%nrow,B_ncol,1,A_csr%nzval,A_csr%colind,A_csr%rowpnt,s, &
         B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind,C_csr%rowpnt,&
         C_csr%nnz,iw,ierr)

    if (ierr.ne.0) call error_msg('(zcplamub_st)',OUTOFBOUND)

    call log_deallocate(iw)

  end subroutine zcplsamub_st

  function nnz(A_csr,op,B_csr) result(C_nnz)
    type(z_CSR) :: A_csr, B_csr
    character(1) :: op
    integer :: C_nnz


    select case(op)
    case('+')
      if (A_csr%sorted .and. B_csr%sorted) then
        C_nnz = nnz_sum1(A_csr, B_csr)
      else
        C_nnz = nnz_sum(A_csr, B_csr)
      end if
    case('*')
      C_nnz = nnz_mult(A_csr, B_csr)
    end select

  end function nnz

  !***************************************************************************************
  !Subroutine per l'allocazione del prodotto sparso A x B
  !Verifica quanti elementi non nulli si avranno dal prodotto di due matrici sparse CSR
  !***************************************************************************************
  function nnz_mult(A_csr,B_csr) result (C_nnz)

    !************************************************************************************
    !
    !Input:
    !A_csr: A matrix in CSR format
    !B_csr: B_matrix in csr format
    !A_ncol: number of columns of A matrix
    !Output:
    !C_nnz: number of nonzero values needed in C_nnz
    !
    !************************************************************************************

    implicit none

    !Input/Output arguments
    type(z_CSR) :: A_csr, B_csr, C_csr
    integer :: A_ncol, B_ncol, C_nnz

    !Work Arguments
    integer, DIMENSION(:), ALLOCATABLE :: iw
    integer :: ierr

    !Alloca le parti di C_csr di interesse
    A_ncol=A_csr%ncol
    B_ncol=B_csr%ncol

    call log_allocate(C_csr%nzval,MISMATCH)
    call log_allocate(C_csr%colind,A_csr%nrow*B_ncol)
    call log_allocate(C_csr%rowpnt,A_csr%nrow+1)
    !Allocazione work array iw
    call log_allocate(iw,A_ncol)

    call amub(A_csr%nrow,B_ncol,0,A_csr%nzval,A_csr%colind,A_csr%rowpnt,&
               B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind, &
               C_csr%rowpnt,A_csr%nrow*B_ncol,iw,ierr)

    if (ierr.ne.0) call error_msg('(nnz_mult)',OUTOFBOUND)

    call log_deallocate(C_csr%nzval)
    call log_deallocate(C_csr%colind)
    call log_deallocate(iw)

    C_nnz=C_csr%rowpnt(A_csr%nrow+1)-1

    call log_deallocate(C_csr%rowpnt)

  end function nnz_mult

  !***************************************************************************************
  !Subroutine per l'allocazione della somma sparsa A+B
  !Verifica quanti elementi non nulli si avranno dalla somma di due matrici sparse CSR
  !***************************************************************************************
  function nnz_sum(A_csr,B_csr) result(C_nnz)

    !************************************************************************************
    !
    !Input:
    !A_csr: A matrix in CSR format
    !B_csr: B_matrix in csr format
    !A_ncol: number of columns of A matrix
    !B_ncol: number of columns of B matrix
    !Output:
    !C_nnz: number of nonzero values needed in C_nnz
    !************************************************************************************


    implicit none

    !Input/Output arguments
    type(z_CSR) :: A_csr, B_csr, C_csr
    integer :: A_ncol, C_nnz

    !Work Arguments
    integer, DIMENSION(:), ALLOCATABLE :: iw
    integer :: ierr

    !Alloca le parti di C_csr di interesse
    A_ncol=A_csr%ncol

    call log_allocate(C_csr%nzval,MISMATCH)
    call log_allocate(C_csr%colind,A_csr%nrow*A_ncol)
    call log_allocate(C_csr%rowpnt,A_csr%nrow+1)
    !Allocazione work array iw
    call log_allocate(iw,A_ncol)

    call aplb (A_csr%nrow,A_ncol,0,A_csr%nzval,A_csr%colind,A_csr%rowpnt,&
                B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind,&
                C_csr%rowpnt,A_csr%nrow*A_ncol,iw,ierr)

    if (ierr.ne.0) call error_msg('(nnz_sum)',OUTOFBOUND)

    C_nnz = C_csr%rowpnt(A_csr%nrow+1)-1

    call log_deallocate(C_csr%nzval)
    call log_deallocate(C_csr%colind)
    call log_deallocate(iw)

    call log_deallocate(C_csr%rowpnt)

  end function nnz_sum

  !************************************************************************************
  !Subroutine per l'allocazione della somma sparsa A+B (VERSIONE SORTED)
  !Verifica quanti elementi non nulli si avranno dalla somma di due matrici sparse CSR
  !************************************************************************************
  function nnz_sum1(A_csr,B_csr) result(C_nnz)

    !***********************************************************************************
    !
    !Input:
    !A_csr: A matrix in CSR format
    !B_csr: B_matrix in csr format
    !A_ncol: number of columns of A matrix
    !Output:
    !C_nnz: number of nonzero values needed in C_nnz
    !
    !***********************************************************************************


    implicit none

    !Input/Output arguments
    type(z_CSR) :: A_csr, B_csr
    complex(kind=dp), DIMENSION(:), ALLOCATABLE :: nzval
    integer, DIMENSION(:), ALLOCATABLE :: colind
    integer, DIMENSION(:), ALLOCATABLE :: rowpnt
    integer :: A_ncol, C_nnz

    !Work Arguments
    integer :: ierr


    !Alloca le parti di C_csr di interesse
    A_ncol=A_csr%ncol

    call log_allocate(nzval,MISMATCH)
    call log_allocate(colind,A_csr%nrow*A_ncol)
    call log_allocate(rowpnt,A_csr%nrow+1)

    call aplb1(A_csr%nrow,A_ncol,0,A_csr%nzval,A_csr%colind,A_csr%rowpnt,&
                 B_csr%nzval,B_csr%colind,B_csr%rowpnt,nzval,colind,&
                 rowpnt,A_csr%nrow*A_ncol,ierr)

    if (ierr.ne.0) call error_msg('(nnz_sum1)',OUTOFBOUND)

    C_nnz=rowpnt(A_csr%nrow+1)-1

    call log_deallocate(nzval)
    call log_deallocate(colind)

    call log_deallocate(rowpnt)

  end function nnz_sum1

  !*********************************************************************************
  !                                                                                |
  !  Function for checking non zero values in a submatrix of a sparse CSR matrix   |
  !                                                                                |
  !                                                                                |
  !*********************************************************************************

  function zcheck_nnz(A_csr, i1, i2, j1, j2) result(nnz)

    !*********************************************************************************
    !                                                                                |
    !Input:                                                                          |
    !A_csr: sparse CSR matrix                                                        |
    !i1: starting row                                                                |
    !i2: ending row                                                                  |
    !j1: starting column                                                             |
    !j2: ending column                                                               |
    !                                                                                |
    !Output:                                                                         |
    !nzval: non zero values found in submatrix specified by i1,i2,j1,j2              |
    !                                                                                |
    !*********************************************************************************

    type(z_CSR) :: A_csr
    integer :: i1,i2,j1,j2,sorted,nnz
    integer :: i,j

    !Check i1, i2, j1, j2 validity
    if ((i1.lt.1).or.(i2.gt.A_csr%nrow).or.(j2.lt.j1) &
        & .or.(i2.lt.i1).or.(j1.lt.1).or.(j2.gt.A_csr%ncol)) then
       call error_msg('(zcheck_nnz)',BADINDEX)
    endif

    nnz=0

    IF (A_csr%nnz.NE.0) THEN

       do i=i1,i2

          if (A_csr%rowpnt(i).eq.A_csr%rowpnt(i+1)) then
             cycle
          else

             do j=A_csr%rowpnt(i),A_csr%rowpnt(i+1)-1

                if (A_csr%colind(j).ge.j1) then
                   if (A_csr%colind(j).le.j2) then
                      nnz=nnz+1
                   else
                      if (A_csr%sorted) exit
                   endif
                else
                   cycle
                endif

             enddo

          endif

       enddo

    ELSE

       nnz=0

    ENDIF

  end function zcheck_nnz


  !********************************************************************************
  !
  !  Function for dropping quasi-zero values and checking number of non-zero
  !  values in a dense matrix
  !
  !********************************************************************************
  !
  !Nota: la funzione e' implementata col valore assoluto, pertanto dovrebbe andare
  !bene pure coi complessi

  function zchkdrp(A,drop) result(nnz)

    implicit none
    integer :: nnz,i,j

    Type(z_DNS) :: A
    !Questo deve rimanere reale
    real(kind=dp) :: drop

    nnz=0
    do j=1,A%ncol
       do i=1,A%nrow
          if (abs(A%val(i,j)).lt.drop) then
             A%val(i,j)=(0.0_dp, 0.0_dp)
          else
             nnz=nnz+1
          endif
       enddo
    enddo


  end function zchkdrp
  !--------------------------------------------------------------------------------

  subroutine zmultcsr(A_csr,B_csr,C_csr)

    !*****************************************************************
    !
    !Input:
    !A_csr: primo fattore in formato CSR
    !B_csr: secondo fattore in formato  CSR
    !C_csr: risultato in formato CSR (l'allocazione esatta viene eseguita
    !nella subroutine
    !
    !****************************************************************

    implicit none

    type(z_CSR) :: A_csr,B_csr,C_csr
    integer, DIMENSION(:), ALLOCATABLE :: iw
    integer :: ierr,nnz
    integer :: a_bw, b_bw, ml, mu, iband
    real(dp) :: bndav

    IF (A_csr%ncol.NE.B_csr%nrow) THEN
       call error_msg('(zmult_csr)',MISMATCH)
    ENDIF

    IF ((A_csr%nnz.EQ.0).OR.(B_csr%nnz.EQ.0)) THEN

       CALL create(C_csr,A_csr%nrow,B_csr%ncol,0)
       C_csr%rowpnt=1

    ELSE

       ! G. P. preallocation for exact calculation of nonzero values.
       ! The first guess is built depending on A and B bandwidth
       ! I consider the maximum bandwidth
       call bandwidth(A_csr%ncol, A_csr%colind, A_csr%rowpnt, ml, mu, a_bw, bndav)
       a_bw = max(ml + 1, mu + 1) * 2
       call bandwidth(B_csr%ncol, B_csr%colind, B_csr%rowpnt, ml, mu, b_bw, bndav)
       b_bw = max(ml + 1, mu + 1) * 2

       ! Bandwidth of products is two time the bandwith of factors
       nnz = max(a_bw, b_bw) * max(A_csr%nrow, B_csr%ncol) * 2

       !Preliminar product on indexes only. This is used to determine the exact amount
       !of memory which needs to be allocate, avoiding a temporary unused heavy
       !complex array
       call log_allocate(C_csr%nzval,MISMATCH)
       call log_allocate(C_csr%colind,nnz)
       call log_allocate(C_csr%rowpnt,(A_csr%nrow+1))
       C_csr%rowpnt=0

       !Allocazione work array iw
       call log_allocate(iw,B_csr%ncol)
       call amub(A_csr%nrow,B_csr%ncol,0,A_csr%nzval,A_csr%colind,A_csr%rowpnt,&
            B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind,C_csr%rowpnt,&
            nnz,iw,ierr)


       if (ierr.ne.0) call error_msg('(zamub)',OUTOFBOUND)

       !Riallocazioni esatte di C_csr
       call log_deallocate(C_csr%colind)
       call log_deallocate(C_csr%nzval)

       nnz=C_csr%rowpnt(A_csr%nrow+1)-1

       call log_deallocate(C_csr%rowpnt)

       if(nnz.ne.0) then

          call create(C_csr,A_csr%nrow,B_csr%ncol,nnz)

          !Prodotto C_csr=A_csr*B_csr
          call amub(A_csr%nrow,B_csr%ncol,1,A_csr%nzval,A_csr%colind,A_csr%rowpnt,&
               B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind, &
               C_csr%rowpnt,C_csr%nnz,iw,ierr)

          if (ierr.ne.0) call error_msg('(zamub)',OUTOFBOUND)

       else

          call create(C_csr,A_csr%nrow,B_csr%ncol,MISMATCH)
          C_csr%nnz=0

       endif

       call log_deallocate(iw)

    ENDIF

  end subroutine zmultcsr

  !*****************************************************************
  !
  !  Subroutine di moltiplicazione sparsa compatta (CSR unsorted)
  !  con prodotto per scalare  C=s*AxB
  !
  !*****************************************************************

  subroutine zmultcsrs(A_csr,B_csr,s,C_csr)

    !*****************************************************************
    !
    !Input:
    !A_csr: primo fattore in formato CSR
    !B_csr: secondo fattore in formato  CSR
    !C_csr: risultato in formato CSR (l'allocazione esatta viene eseguita
    !nella subroutine
    !
    !****************************************************************

    implicit none

    type(z_CSR) :: A_csr,B_csr,C_csr
    complex(kind=dp) :: s
    integer :: ierr,nnz
    integer :: a_bw, b_bw, ml, mu
    real(dp) :: bndav
    integer, DIMENSION(:), ALLOCATABLE :: iw
    IF (A_csr%ncol.NE.B_csr%nrow) THEN
        call error_msg('(zmultccsr)',MISMATCH)
    ENDIF
        
    !Comment: routine "bandwidth()" does not longer accept argument bndav as integer, so it has been defined as real(dp)

    IF ((A_csr%nnz.EQ.0).OR.(B_csr%nnz.EQ.0).OR.(ABS(s).EQ.0)) THEN

       CALL create(C_csr,A_csr%nrow,B_csr%ncol,0)
       C_csr%rowpnt=1

    else


       ! G. P. preallocation for exact calculation of nonzero values.
       ! The first guess is built depending on A and B bandwidth
       ! I consider the maximum bandwidth
       call bandwidth(A_csr%ncol, A_csr%colind, A_csr%rowpnt, ml, mu, a_bw, bndav)
       a_bw = max(ml + 1, mu + 1) * 2
       call bandwidth(B_csr%ncol, B_csr%colind, B_csr%rowpnt, ml, mu, b_bw, bndav)
       b_bw = max(ml + 1, mu + 1) * 2

       ! Bandwidth of products is two time the bandwith of factors
       nnz = max(a_bw, b_bw) * max(A_csr%nrow, B_csr%ncol) * 2

       !Preliminar product on indexes only. This is used to determine the exact amount
       !of memory which needs to be allocate, avoiding a temporary unused heavy
       !complex array
       call log_allocate(C_csr%nzval,1)
       call log_allocate(C_csr%colind,nnz)
       call log_allocate(C_csr%rowpnt,(A_csr%nrow+1))
       C_csr%rowpnt=0

       !Allocazione work array iw
       call log_allocate(iw,B_csr%ncol)

       call amub(A_csr%nrow,B_csr%ncol,0,A_csr%nzval,A_csr%colind,A_csr%rowpnt,&
            B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind,&
            C_csr%rowpnt,nnz,iw,ierr)

       if (ierr.ne.0) call error_msg('(zamub)',OUTOFBOUND)

       !Riallocazioni esatte di C_csr
       nnz=C_csr%rowpnt(A_csr%nrow+1)-1

       call log_deallocate(C_csr%colind)
       call log_deallocate(C_csr%nzval)
       call log_deallocate(C_csr%rowpnt)

       if(nnz.ne.0) then

          call create(C_csr,A_csr%nrow,B_csr%ncol,nnz)

          !Prodotto C_csr=A_csr*B_csr
          call amubs(A_csr%nrow,B_csr%ncol,1,A_csr%nzval,A_csr%colind,A_csr%rowpnt,s, &
               B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind, &
               C_csr%rowpnt,C_csr%nnz,iw,ierr)

          if (ierr.ne.0) call error_msg('(zamubs)',OUTOFBOUND)

       else

          call create(C_csr,A_csr%nrow,B_csr%ncol,MISMATCH)
          C_csr%nnz=0

       endif

       call log_deallocate(iw)

    ENDIF

  end subroutine zmultcsrs

  !*****************************************************************
  !
  !  Subroutine di moltiplicazione densa compatta (DNS)
  !  con prodotto per scalare  C=s*A*B
  !
  !*****************************************************************
  subroutine zmultdns(A_dns,B_dns,C_dns)

    !*****************************************************************
    !
    !Input:
    !A_dns: primo fattore in formato DNS
    !B_dns: secondo fattore in formato  DNS
    !C_dns: risultato in formato DNS (l'allocazione esatta viene eseguita
    !nella subroutine
    !
    !****************************************************************

    implicit none

    type(z_DNS) :: A_dns,B_dns,C_dns
    integer :: M,N,K
    complex(dp), parameter :: s = (1.0_dp, 0.0_dp)
    complex(dp) :: beta

    IF (A_dns%ncol.NE.B_dns%nrow) THEN
       call error_msg('(zmultdns) A B',MISMATCH)
    ENDIF

    M = A_dns%nrow
    N = B_dns%ncol
    K = A_dns%ncol

    IF (allocated(C_dns%val)) THEN
      IF(C_dns%nrow .ne. M .or. C_dns%ncol .ne. N) THEN
         call error_msg('(zmultdns) C',MISMATCH)
      ENDIF
      beta = (1.0_dp,0.0_dp)
    ELSE
      CALL create(C_dns,M,N)
      beta = (0.0_dp,0.0_dp)
    ENDIF

    CALL ZGEMM('N','N', M, N, K, s, A_dns%val, M, &
            B_dns%val, B_dns%nrow, beta, C_dns%val, M)

  end subroutine zmultdns

  !***********************************************************************

  subroutine zmatmul(A,B,C_dns)

    !*****************************************************************
    !
    !Input:
    !A_dns: primo fattore in formato DNS
    !B_dns: secondo fattore in formato  DNS
    !C_dns: risultato in formato DNS (l'allocazione esatta viene eseguita
    !nella subroutine
    !
    !****************************************************************

    implicit none

    complex(dp), Dimension(:,:) :: A,B
    type(z_DNS) :: C_dns
    integer :: M,N,K
    complex(dp), parameter :: s = (1.0_dp, 0.0_dp)
    complex(dp) :: beta

    IF (size(A,2).NE.size(B,MISMATCH)) THEN
       call error_msg('(zmatmul) A B',MISMATCH)
    ENDIF

    M = size(A,MISMATCH)
    N = size(B,2)
    K = size(A,2)

    IF (allocated(C_dns%val)) THEN
      IF(C_dns%nrow .ne. M .or. C_dns%ncol .ne. N) THEN
         call error_msg('(zmultdnss) C',MISMATCH)
      ENDIF
      beta = (1.0_dp,0.0_dp)
    ELSE
      CALL create(C_dns,M,N)
      beta = (0.0_dp,0.0_dp)
    ENDIF

    CALL ZGEMM('N','N', M, N, K, s, A, M, &
            B, K, beta, C_dns%val, M)

  end subroutine zmatmul

  !***********************************************************************

  subroutine zmatmuls(A,B,s,C_dns)

    !*****************************************************************
    !
    !Input:
    !A_dns: primo fattore in formato DNS
    !B_dns: secondo fattore in formato  DNS
    !C_dns: risultato in formato DNS (l'allocazione esatta viene eseguita
    !nella subroutine
    !
    !****************************************************************

    implicit none

    complex(dp), Dimension(:,:) :: A,B
    type(z_DNS) :: C_dns
    complex(dp) :: s

    integer :: M,N,K
    complex(dp) ::  beta

    IF (size(A,2).NE.size(B,MISMATCH)) THEN
       call error_msg('(zmatmuls) A B',MISMATCH)
    ENDIF

    M = size(A,MISMATCH)
    N = size(B,2)
    K = size(A,2)

    !L = size(B,MISMATCH) = K)
    IF (allocated(C_dns%val)) THEN
      IF(C_dns%nrow .ne. M .or. C_dns%ncol .ne. N) THEN
       call error_msg('(zmatmuls) C',MISMATCH)
      ENDIF
      beta = (1.0_dp,0.0_dp)
    ELSE
      CALL create(C_dns,M,N)
      beta = (0.0_dp,0.0_dp)
    ENDIF

    CALL ZGEMM('N','N', M, N, K, s, A, M, &
            B, K, beta, C_dns%val, M)

  end subroutine zmatmuls

  !***********************************************************************



  !*****************************************************************
  !
  !  Subroutine di moltiplicazione densa compatta (DNS)
  !  con prodotto per scalare  C=C+s*A*B
  !
  !*****************************************************************
  subroutine zmultdnss(A_dns,B_dns,s,C_dns)

    !*****************************************************************
    !
    !Input:
    !A_dns: primo fattore in formato DNS
    !B_dns: secondo fattore in formato  DNS
    !C_dns: risultato in formato DNS (l'allocazione esatta viene eseguita
    !nella subroutine
    !
    !****************************************************************

    implicit none

    type(z_DNS) :: A_dns,B_dns,C_dns
    complex(dp) :: s
    integer :: M,N,K
    complex(dp) :: beta

    IF (A_dns%ncol.NE.B_dns%nrow) THEN
       call error_msg('(zmultdnss) A B',MISMATCH)
    ENDIF

    M = A_dns%nrow
    N = B_dns%ncol
    K = A_dns%ncol

    IF (allocated(C_dns%val)) THEN
      IF(C_dns%nrow .ne. M .or. C_dns%ncol .ne. N) THEN
        call error_msg('(zmultdnss) C',MISMATCH)
      ENDIF
      beta = (1.0_dp,0.0_dp)
    ELSE
      CALL create(C_dns,M,N)
      beta = (0.0_dp,0.0_dp)
    ENDIF

    ! C = beta C + s A * B
    CALL ZGEMM('N','N', M, N, K, s, A_dns%val, M, &
            B_dns%val, B_dns%nrow, beta, C_dns%val, M)

  end subroutine zmultdnss

  !***********************************************************************

  subroutine zsumcsr(A_csr,B_csr,C_csr)

    implicit none

    type(z_CSR) :: A_csr,B_csr,C_csr
    integer, DIMENSION(:), ALLOCATABLE :: iw
    integer :: ierr,A_ncol

    if(A_csr%nrow.ne.B_csr%nrow .or. A_csr%ncol.ne.B_csr%ncol) THEN
       call error_msg('(zsumcsr) ',MISMATCH)
    endif

    IF ((A_csr%nnz.EQ.0).AND.(B_csr%nnz.EQ.0)) THEN
       CALL  create(C_csr,A_csr%nrow,A_csr%ncol,0)
       C_csr%rowpnt=1
    ELSEIF (A_csr%nnz.EQ.0) THEN
       CALL create(C_csr,B_csr%nrow,B_csr%ncol,B_csr%nnz)
       C_csr%nzval=B_csr%nzval; C_csr%colind=B_csr%colind; C_csr%rowpnt=B_csr%rowpnt;
    ELSEIF (B_csr%nnz.EQ.0) THEN
       CALL create(C_csr,A_csr%nrow,A_csr%ncol,A_csr%nnz)
       C_csr%nzval=A_csr%nzval; C_csr%colind=A_csr%colind; C_csr%rowpnt=A_csr%rowpnt;

    ELSE

       C_csr%nrow=A_csr%nrow
       C_csr%ncol=A_csr%ncol

       !Alloca le parti di C_csr di interesse
       A_ncol=A_csr%ncol
       C_csr%nnz=(A_csr%nnz+B_csr%nnz)

       call log_allocate(C_csr%nzval,MISMATCH)
       call log_allocate(C_csr%colind,C_csr%nnz)
       call log_allocate(C_csr%rowpnt,(C_csr%nrow+1))

       !Allocazione work array iw
       call log_allocate(iw,A_ncol)

       call aplb(A_csr%nrow,A_ncol,0,A_csr%nzval,A_csr%colind,A_csr%rowpnt,&
            B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind,&
            C_csr%rowpnt,C_csr%nnz,iw,ierr)

       if (ierr.ne.0) call error_msg('()',3)

       C_csr%nnz=C_csr%rowpnt(A_csr%nrow+1)-1

       call log_deallocate(C_csr%nzval)
       call log_deallocate(C_csr%colind)

       call log_allocate(C_csr%nzval,C_csr%nnz)
       call log_allocate(C_csr%colind,C_csr%nnz)

       call aplb(A_csr%nrow,A_ncol,1,A_csr%nzval,A_csr%colind,A_csr%rowpnt,&
            B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind,&
            C_csr%rowpnt,C_csr%nnz,iw,ierr)

       if (ierr.ne.0) call error_msg('(zsmcsr)',OUTOFBOUND)

       call log_deallocate(iw)

    ENDIF

  end subroutine zsumcsr


  !***********************************************************************
  !
  !  Subroutine di somma sparsa compatta (CSR sorted)
  !
  !***********************************************************************

  subroutine zsumcsr1(A_csr,B_csr,C_csr)


    implicit none

    type(z_CSR) :: A_csr,B_csr,C_csr
    integer :: ierr,A_ncol

    if(A_csr%nrow.ne.B_csr%nrow) STOP 'Error in aplb subroutine: nrow differ'
    if(A_csr%ncol.ne.B_csr%ncol) STOP 'Error in aplb subroutine: ncol differ'

    IF ((A_csr%nnz.EQ.0).AND.(B_csr%nnz.EQ.0)) THEN
       CALL  create(C_csr,A_csr%nrow,A_csr%ncol,0)
       C_csr%rowpnt=1
    ELSEIF (A_csr%nnz.EQ.0) THEN
       CALL create(C_csr,B_csr%nrow,B_csr%ncol,B_csr%nnz)
       C_csr%nzval=B_csr%nzval; C_csr%colind=B_csr%colind; C_csr%rowpnt=B_csr%rowpnt;
    ELSEIF (B_csr%nnz.EQ.0) THEN
       CALL create(C_csr,A_csr%nrow,A_csr%ncol,A_csr%nnz)
       C_csr%nzval=A_csr%nzval; C_csr%colind=A_csr%colind; C_csr%rowpnt=A_csr%rowpnt;

    ELSE

       C_csr%nrow=A_csr%nrow
       C_csr%ncol=A_csr%ncol

       !Alloca le parti di C_csr di interesse
       A_ncol=A_csr%ncol
       C_csr%nnz =(A_csr%nnz+B_csr%nnz)
       call log_allocate(C_csr%nzval,MISMATCH)
       call log_allocate(C_csr%colind,C_csr%nnz)
       call log_allocate(C_csr%rowpnt,(C_csr%nrow+1))

       call aplb1(A_csr%nrow,A_ncol,0,A_csr%nzval,A_csr%colind,A_csr%rowpnt,&
            B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind,&
            C_csr%rowpnt,C_csr%nnz,ierr)

       if (ierr.ne.0) call error_msg('(aplb1)',OUTOFBOUND)

       C_csr%nnz=C_csr%rowpnt(A_csr%nrow+1)-1

       call log_deallocate(C_csr%nzval)
       call log_deallocate(C_csr%colind)

       call log_allocate(C_csr%nzval,C_csr%nnz)
       call log_allocate(C_csr%colind,C_csr%nnz)

       call aplb1(A_csr%nrow,A_ncol,1,A_csr%nzval,A_csr%colind,A_csr%rowpnt,&
            B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind, &
            C_csr%rowpnt,C_csr%nnz,ierr)

       if (ierr.ne.0) call error_msg('(aplb1)',OUTOFBOUND)

    ENDIF

  end subroutine zsumcsr1

  !***********************************************************************
  !
  !  Subroutine di somma sparsa compatta con prodotto per scalare complessa
  !
  !***********************************************************************

  subroutine zsumcsrs(A_csr,B_csr,s,C_csr)


    implicit none

    type(z_CSR) :: A_csr,B_csr,C_csr
    complex(kind=dp) :: s
    integer, DIMENSION(:), ALLOCATABLE :: iw
    integer :: ierr,A_ncol

    if(A_csr%nrow.ne.B_csr%nrow) STOP 'Error in aplb subroutine: nrow differ'
    if(A_csr%ncol.ne.B_csr%ncol) STOP 'Error in aplb subroutine: ncol differ'

    IF ((A_csr%nnz.EQ.0).AND.(B_csr%nnz.EQ.0)) THEN
       CALL  create(C_csr,A_csr%nrow,A_csr%ncol,0)
       C_csr%rowpnt=1
    ELSEIF (A_csr%nnz.EQ.0) THEN
       CALL create(C_csr,B_csr%nrow,B_csr%ncol,B_csr%nnz)
       C_csr%nzval=B_csr%nzval; C_csr%colind=B_csr%colind; C_csr%rowpnt=B_csr%rowpnt;
    ELSEIF (B_csr%nnz.EQ.0) THEN
       CALL create(C_csr,A_csr%nrow,A_csr%ncol,A_csr%nnz)
       C_csr%nzval=A_csr%nzval; C_csr%colind=A_csr%colind; C_csr%rowpnt=A_csr%rowpnt;

    ELSE

       C_csr%nrow=A_csr%nrow
       C_csr%ncol=A_csr%ncol

       !Alloca le parti di C_csr di interesse
       A_ncol=A_csr%ncol
       C_csr%nnz=(A_csr%nnz+B_csr%nnz)

       call log_allocate(C_csr%nzval,MISMATCH)
       call log_allocate(C_csr%colind,C_csr%nnz)
       call log_allocate(C_csr%rowpnt,(C_csr%nrow+1))

       call log_allocate(iw,A_ncol)

       call aplb(A_csr%nrow,A_ncol,0,A_csr%nzval,A_csr%colind,A_csr%rowpnt,&
            B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind,&
            C_csr%rowpnt,C_csr%nnz,iw,ierr)

       if (ierr.ne.0) call error_msg('(zaplb)',OUTOFBOUND)

       C_csr%nnz=C_csr%rowpnt(A_csr%nrow+1)-1

       if( C_csr%nnz.gt.(A_csr%nnz+B_csr%nnz) ) write(*,*) C_csr%nnz,'!!!!'

       call log_deallocate(C_csr%nzval)
       call log_deallocate(C_csr%colind)

       call log_allocate(C_csr%nzval,C_csr%nnz)
       call log_allocate(C_csr%colind,C_csr%nnz)

       call aplsb(A_csr%nrow,A_ncol,1,A_csr%nzval,A_csr%colind,A_csr%rowpnt,s,&
            B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind,&
            C_csr%rowpnt,C_csr%nnz,iw,ierr)

       if (ierr.ne.0) call error_msg('(zaplb)',OUTOFBOUND)

       call log_deallocate(iw)

    ENDIF

  end subroutine zsumcsrs



  !***********************************************************************
  !
  !  Subroutine di somma sparsa compatta con prodotto per scalare
  !
  !***********************************************************************

  subroutine zsumcsrs1s2(A_csr,B_csr,s1,s2,C_csr)

    implicit none

    type(z_CSR) :: A_csr,B_csr,C_csr
    complex(kind=dp) :: s1,s2
    integer, DIMENSION(:), ALLOCATABLE :: iw
    integer :: ierr,A_ncol


    if(A_csr%nrow.ne.B_csr%nrow) STOP 'Error in aplb subroutine: nrow differ'
    if(A_csr%ncol.ne.B_csr%ncol) STOP 'Error in aplb subroutine: ncol differ'

    IF ((A_csr%nnz.EQ.0).AND.(B_csr%nnz.EQ.0)) THEN
       CALL  create(C_csr,A_csr%nrow,A_csr%ncol,0)
       C_csr%rowpnt=1
    ELSEIF (A_csr%nnz.EQ.0) THEN
       CALL create(C_csr,B_csr%nrow,B_csr%ncol,B_csr%nnz)
       C_csr%nzval=B_csr%nzval; C_csr%colind=B_csr%colind; C_csr%rowpnt=B_csr%rowpnt;
    ELSEIF (B_csr%nnz.EQ.0) THEN
       CALL create(C_csr,A_csr%nrow,A_csr%ncol,A_csr%nnz)
       C_csr%nzval=A_csr%nzval; C_csr%colind=A_csr%colind; C_csr%rowpnt=A_csr%rowpnt;

    ELSE

       !Alloca le parti di C_csr di interesse
       C_csr%nrow=A_csr%nrow
       C_csr%ncol=A_csr%ncol
       A_ncol=A_csr%ncol
       C_csr%nnz=A_csr%nnz+B_csr%nnz
       call log_allocate(C_csr%nzval,MISMATCH)
       call log_allocate(C_csr%colind, C_csr%nnz)
       call log_allocate(C_csr%rowpnt,(C_csr%nrow+1))
       call log_allocate(iw,A_ncol)

       ierr=0
       !call zas1pls2b(A_csr%nrow,A_ncol,0,A_csr%nzval,A_csr%colind,A_csr%rowpnt,s1,s2,&
       !     B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind,C_csr%rowpnt,&
       !     C_csr%nnz,iw,ierr)
       call aplb(A_csr%nrow,A_ncol,0,A_csr%nzval,A_csr%colind,A_csr%rowpnt,&
            B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind,&
            C_csr%rowpnt,C_csr%nnz,iw,ierr)

       if (ierr.ne.0) call error_msg('(zsumcsrs1s2)',OUTOFBOUND)

       C_csr%nnz=C_csr%rowpnt(A_csr%nrow+1)-1

       call log_deallocate(C_csr%nzval)
       call log_deallocate(C_csr%colind)

       call log_allocate(C_csr%nzval,C_csr%nnz)
       call log_allocate(C_csr%colind,C_csr%nnz)

       call as1pls2b(A_csr%nrow,A_ncol,1,A_csr%nzval,A_csr%colind,A_csr%rowpnt,s1,s2,&
            B_csr%nzval,B_csr%colind,B_csr%rowpnt,C_csr%nzval,C_csr%colind,C_csr%rowpnt,&
            C_csr%nnz,iw,ierr)

       if (ierr.ne.0) call error_msg('(zsumcsrs1s2)',OUTOFBOUND)

       call log_deallocate(iw)

    ENDIF

  end subroutine zsumcsrs1s2

    !***********************************************************************
    !
    !  Subroutine di somma densa
    !
    !***********************************************************************

  subroutine zsumdns(A_dns,B_dns,C_dns)

    !*****************************************************************
    !
    !Input:
    !A_dns: primo fattore in formato DNS
    !B_dns: secondo fattore in formato  DNS
    !C_dns: risultato in formato DNS (l'allocazione esatta viene eseguita
    !nella subroutine
    !
    !****************************************************************

    implicit none

    type(z_DNS) :: A_dns,B_dns,C_dns

    IF (A_dns%ncol.NE.B_dns%ncol .AND. A_dns%nrow.NE.B_dns%nrow) THEN
        call error_msg('(zsumdns)',MISMATCH)
    ENDIF


    CALL create(C_dns,A_dns%nrow,B_dns%ncol)

    C_dns%val = A_dns%val + B_dns%val

  end subroutine zsumdns

  !***********************************************************************
  subroutine zsumdnss(A_dns,B_dns,s,C_dns)

    implicit none

    type(z_DNS) :: A_dns,B_dns,C_dns
    complex(dp) :: s

    IF (A_dns%ncol.NE.B_dns%ncol .AND. A_dns%nrow.NE.B_dns%nrow) THEN
       call error_msg('(zsumdnss)',MISMATCH)
    ENDIF


    CALL create(C_dns,A_dns%nrow,B_dns%ncol)

    C_dns%val = A_dns%val + s * B_dns%val

  end subroutine zsumdnss

  !***********************************************************************
  subroutine zsumdnss1s2(A_dns,B_dns,s1,s2,C_dns)

    implicit none

    type(z_DNS) :: A_dns,B_dns,C_dns
    complex(dp) :: s1,s2

    IF (A_dns%ncol.NE.B_dns%ncol .AND. A_dns%nrow.NE.B_dns%nrow) THEN
        call error_msg('(zsumdnss1s2)',MISMATCH)
    ENDIF


    CALL create(C_dns,A_dns%nrow,B_dns%ncol)

    C_dns%val = s1 * A_dns%val + s2 * B_dns%val

  end subroutine zsumdnss1s2


  !***********************************************************************

  !**************************************************************************
  !
  !  Subroutine di estrazione di sottomatrice, allocazione in loco
  !  (A_csr unsorted)
  !
  !**************************************************************************

  subroutine zextract_csr(A_csr,i1,i2,j1,j2,A_sub)

    !********************************************************************
    !
    !Input:
    !A_csr: matrice da cui si vuole estrarre la sottomatrice (unsorted)
    !i1,i2: righe iniziale e finale (incluse)
    !j1,j2: colonne iniziale e finale (incluse)
    !Output:
    !A_sub: matrice estratta (allocata esatta nella subroutine)
    !
    !********************************************************************

    implicit none

    type(z_CSR) :: A_csr,A_sub
    integer :: i1,i2,j1,j2
    integer :: nnz

    IF ((i1.GT.i2).OR.(j1.GT.j2).OR.(i2.GT.A_csr%nrow).OR.(j2.GT.A_csr%ncol)) THEN
       print*, 'ERROR (zextract_csr): bad indeces specification';
       print*, 'Trying to extract block from matrix',A_csr%nrow,'x',A_csr%ncol
       print*, 'Indices Rows',i1,i2,'Cols',j1,j2
       STOP
    ENDIF

    nnz = zcheck_nnz(A_csr, i1, i2, j1, j2);

    IF (nnz.NE.0) THEN
       CALL create(A_sub,(i2-i1+1),(j2-j1+1),nnz)
       call zsubmat_st(A_csr,i1,i2,j1,j2,A_sub)

    ELSE
       CALL create(A_sub,(i2-i1+1),(j2-j1+1),MISMATCH)
       A_sub%rowpnt=1

    ENDIF

  end subroutine zextract_csr

  !------------------------------------------------------------------------------
  subroutine zextract_dns(A_csr,i1,i2,j1,j2,A_dns)
    type(z_CSR) :: A_csr
    type(z_DNS) :: A_dns
    integer :: i1,i2,j1,j2

    integer :: i,k

    IF ((i1.GT.i2).OR.(j1.GT.j2).OR.(i2.GT.A_csr%nrow).OR.(j2.GT.A_csr%ncol)) THEN
       print*, 'ERROR (zextract_dns): bad indeces specification';
       print*, 'Trying to extract block from matrix',A_csr%nrow,'x',A_csr%ncol
       print*, 'Indices Rows',i1,i2,'Cols',j1,j2
       STOP
    ENDIF

    call create(A_dns,(i2-i1+1),(j2-j1+1))
    A_dns%val=(0.0_dp,0.0_dp)

    do i = i1, i2
       do k = A_csr%rowpnt(i), A_csr%rowpnt(i+1)-1
          if(A_csr%colind(k).ge.j1 .and. A_csr%colind(k).le.j2) then
             A_dns%val(i-i1+1,A_csr%colind(k)-j1+1) = A_csr%nzval(k)
          endif
       enddo
    enddo

  end subroutine zextract_dns
  !------------------------------------------------------------------------------

  !**************************************************************************
  !
  !  SUBROUTINE per la concatenazione di matrici sparse (concatena blocchi in place)
  !  o per la somma di un sottoblocco all'interno di una matrice pi grande
  !
  !*************************************************************************
  SUBROUTINE zconcat_csr(A_csr,B_csr,i1,j1)

    !*************************************************************************
    !Nota: per A quadrata
    !Input:
    !A_csr: matrice di partenza (contiene i blocchi precedenti)
    !B_csr: matrice da concatenare
    !i1,j1: indice di riga e di colonna in cui viene inserita B_csr
    !
    !A_csr viene opportunamente riallocata per contenere anche i nuovi valori
    !
    !Nota: non viene effettuato nessun controllo su A_csr, se A_csr e' non nulla
    !nella regione di concatenazione la SUBROUTINE effettua la somma dei valori
    !gia' presenti con quelli di B_csr
    !*************************************************************************

    IMPLICIT NONE

    TYPE(z_CSR) :: A_csr, B_csr,D_csr
    TYPE(z_CSR) :: C_csr
    INTEGER :: i1,j1

    IF (A_csr%nrow.lt.(B_csr%nrow+i1-1)) THEN
        call error_msg('(zconcat)',BADINDEX)
    ENDIF

    IF (B_csr%nnz.eq.0) RETURN

    !Crea una matrice C sparsa di zeri che contiene B come sottomatrice
    !nella regione di interesse
    CALL create(C_csr,A_csr%nrow,max(B_csr%ncol+j1-1,A_csr%ncol),B_csr%nnz)

    !Inizializza C_csr%rowpnt
    C_csr%rowpnt=1
    C_csr%colind=0
    !Assign B_csr elements and indexes to C_csr
    C_csr%nzval=B_csr%nzval
    C_csr%colind=B_csr%colind+j1-1
    C_csr%rowpnt(i1:B_csr%nrow+i1)=B_csr%rowpnt


    IF ( (B_csr%nrow+i1).lt.C_csr%nrow ) THEN
       ! repeat last rowpnt for the remaining rows.
       C_csr%rowpnt((B_csr%nrow+i1+1):(C_csr%nrow+1))=C_csr%rowpnt(B_csr%nrow+i1)
    ENDIF

    call zsumcsr(A_csr,C_csr,D_csr)

    !Assegna in A_csr i valori di D_csr e distrugge D_csr
    CALL destroy(C_csr)
    CALL destroy(A_csr)
    CALL create(A_csr,D_csr%nrow,D_csr%ncol,D_csr%nnz)

    A_csr%nzval=D_csr%nzval
    A_csr%colind=D_csr%colind
    A_csr%rowpnt=D_csr%rowpnt

    CALL destroy(D_csr)

  END SUBROUTINE zconcat_csr

  !**************************************************************************
  !
  !  SUBROUTINE per la concatenazione di matrici sparse (concatena blocchi in place)
  !  o per somma con prodotto per scalare di un sottoblocco all'interno di una matrice pi grande
  !  A+sB
  !*************************************************************************
  SUBROUTINE zconcatm_csr(A_csr,s,B_csr,i1,j1)

    !*************************************************************************
    !Nota: per A quadrata
    !Input:
    !A_csr: matrice di partenza (contiene i blocchi precedenti)
    !B_csr: matrice da concatenare
    !i1,j1: indice di riga e di colonna in cui viene inserita B_csr
    !
    !A_csr viene opportunamente riallocata per contenere anche i nuovi valori
    !
    !Nota: non viene effettuato nessun controllo su A_csr, se A_csr e' non nulla
    !nella regione di concatenazione la SUBROUTINE effettua la somma dei valori
    !gia' presenti con quelli di B_csr
    !*************************************************************************

    IMPLICIT NONE

    TYPE(z_CSR) :: A_csr, B_csr,D_csr
    TYPE(z_CSR) :: C_csr
    INTEGER :: i1,i2,j1
    COMPLEX(kind=dp) :: s

    IF (A_csr%nrow.lt.(B_csr%nrow+i1-1)) THEN
       call error_msg('(zconcatm_csr)',BADINDEX)
    ENDIF

    IF (B_csr%nnz.EQ.0) RETURN

    !Crea una matrice C sparsa di zeri che contiene B come sottomatrice
    !nella regione di interesse
    CALL create(C_csr,A_csr%nrow,max(B_csr%ncol+j1-1,A_csr%ncol),B_csr%nnz)

    !Inizializza C_csr%rowpnt
    C_csr%rowpnt(1:C_csr%nrow+1)=1
    C_csr%colind(1:C_csr%nnz)=0

    !Assign B_csr elements and indexes to C_csr
    do i2=1,C_csr%nnz
       C_csr%nzval(i2)=B_csr%nzval(i2)
       C_csr%colind(i2)=B_csr%colind(i2)+j1-1
    enddo
    do i2=1,B_csr%nrow+1
       C_csr%rowpnt(i2+i1-1)=B_csr%rowpnt(i2)
    enddo

    IF ((i1+B_csr%nrow).lt.C_csr%nrow) THEN
       C_csr%rowpnt((B_csr%nrow+i1+1):(C_csr%nrow+1))=C_csr%rowpnt(B_csr%nrow+i1)
    ENDIF

    CALL zsumcsrs(A_csr,C_csr,s,D_csr)

    !Assegna in A_csr i valori di D_csr e distrugge D_csr
    CALL destroy(C_csr)
    CALL destroy(A_csr)
    CALL create(A_csr,D_csr%nrow,D_csr%ncol,D_csr%nnz)

    A_csr%nzval=D_csr%nzval
    A_csr%colind=D_csr%colind
    A_csr%rowpnt=D_csr%rowpnt

    CALL destroy(D_csr)

  END SUBROUTINE zconcatm_csr


  !**************************************************************************
  !
  !  SUBROUTINE per la concatenazione di matrici sparse (concatena blocchi in
  !  place) o per somma con prodotto per scalare di un sottoblocco all'interno
  !  di una matrice A+sB.
  !*************************************************************************
  subroutine zrconcatm_csr(A_csr,s,B_csr,ty,i1,j1)

    !*************************************************************************
    !Nota: per A quadrata
    !Input:
    !A_csr: matrice reale di partenza (contiene i blocchi precedenti)
    !B_csr: matrice complessa da concatenare
    !i1,j1: indice di riga e di colonna in cui viene inserita B_csr
    !
    !A_csr viene opportunamente riallocata per contenere anche i nuovi valori
    !
    !Nota: non viene effettuato nessun controllo su A_csr, se A_csr e` non nulla
    !nella regione di concatenazione la SUBROUTINE effettua la somma dei valori
    !gia` presenti con quelli di B_csr
    !*************************************************************************

    IMPLICIT NONE

    TYPE(r_CSR) :: A_csr, D_csr, C_csr
    TYPE(z_CSR) :: B_csr
    CHARACTER(1) :: Ty
    INTEGER :: i1,i2,j1
    COMPLEX(kind=dp) :: s

    IF (A_csr%nrow.lt.(B_csr%nrow+i1-1)) THEN
       call error_msg('(zrconcatm_csr)',BADINDEX)
    ENDIF

    IF (B_csr%nnz.EQ.0) RETURN

    !Crea una matrice C sparsa di zeri che contiene B come sottomatrice
    !nella regione di interesse
    CALL create(C_csr,A_csr%nrow,max(B_csr%ncol+j1-1,A_csr%ncol),B_csr%nnz)
    !Inizializza C_csr%rowpnt
    C_csr%rowpnt(1:C_csr%nrow+1)=1
    C_csr%colind(1:C_csr%nnz)=0

    !Assign B_csr elements and indexes to C_csr
    if (Ty.eq.'I') then
       do i2=1,C_csr%nnz
          C_csr%nzval(i2)=aimag(s*B_csr%nzval(i2))
          C_csr%colind(i2)=B_csr%colind(i2)+j1-1
       enddo
    else
       do i2=1,C_csr%nnz
          C_csr%nzval(i2)=real(s*B_csr%nzval(i2))
          C_csr%colind(i2)=B_csr%colind(i2)+j1-1
       enddo
    endif

    do i2=1,B_csr%nrow+1
       C_csr%rowpnt(i2+i1-1)=B_csr%rowpnt(i2)
    enddo

    IF ((i1+B_csr%nrow).lt.C_csr%nrow) THEN
       C_csr%rowpnt((B_csr%nrow+i1+1):(C_csr%nrow+1))=C_csr%rowpnt(B_csr%nrow+i1)
    ENDIF

    call rsumcsr(A_csr,C_csr,D_csr)

    !Assegna in A_csr i valori di D_csr e distrugge D_csr
    CALL destroy(C_csr)
    CALL destroy(A_csr)
    CALL create(A_csr,D_csr%nrow,D_csr%ncol,D_csr%nnz)

    A_csr%nzval=D_csr%nzval
    A_csr%colind=D_csr%colind
    A_csr%rowpnt=D_csr%rowpnt

    CALL destroy(D_csr)

  end subroutine zrconcatm_csr


  !******************************************************************
  !
  !  Subroutine per la conversione di un sottoblocco della Green
  !  da retarded ad advanced (Hermitiano) con allocazione in loco
  !
  !******************************************************************

  subroutine zdagacsr(A_csr,B_csr)

    !************************************************************************
    !
    !Input:
    !A_csr: matrice di partenza in forma sparsa csr
    !A_ncol: numero di colonne di A_csr
    !B_csr: matrice ottenuta con l'operazione di Hermitiano (allocata in loco)
    !
    !****************************************************************************

    implicit none

    type(z_CSR) :: A_csr,B_csr
    integer :: A_ncol,i
    complex(kind=dp) :: temp


    A_ncol=A_csr%ncol

    call create(B_csr,A_ncol,A_csr%nrow,A_csr%nnz)

    if (A_csr%nnz.eq.0) then

       B_csr%rowpnt=1

    else

       call csrcsc2(A_csr%nrow,A_ncol,1,1,A_csr%nzval,A_csr%colind,A_csr%rowpnt,&
            B_csr%nzval,B_csr%colind,B_csr%rowpnt)

       !Operazione di coniugio su B_csr
       do i=1,B_csr%nnz
          temp=conjg(B_csr%nzval(i))
          B_csr%nzval(i)=temp
       enddo

    endif

  end subroutine zdagacsr

  !******************************************************************
  !
  !  Subroutine per la conversione di un sottoblocco della Green
  !  da retarded ad advanced (Hermitiano) con allocazione in loco
  !  per matrici dense
  !
  !******************************************************************

  subroutine zdagadns(A_dns,B_dns)

    implicit none

    Type(z_dns) :: A_dns, B_dns

    call create(B_dns,A_dns%ncol,A_dns%nrow)

    B_dns%val = conjg(transpose(A_dns%val))


  end subroutine zdagadns

  !*************************************************************************
  !
  !  Subroutine per il calcolo della spectral density associata a un blocco
  !  della funzione di Green data solo le retarded
  !
  !*************************************************************************

  subroutine zspectral_csr(GreenR1,GreenR2,flagR,A)

    !*************************************************************************
    !
    !Input:
    !GreenR1: blocco della Green Retarded
    !GreenR2: blocco della Green Retarded corrispondente alla Green Advanced
    !         che viene calcolata all'interno della subroutine
    !GreenR2_ncol: numero di colonne di GreenR2
    !flag: Se il flag =  1 viene deallocata GreenR2 internamente appena e`
    !      disponibile l'advanced corrispondente, altrimenti no
    !
    !A: Spectral Density corrispondente
    !La funzione calcola A=j*(GreenR1-dagacsr(GreenR2))
    !
    !**************************************************************************

    implicit none

    type(z_CSR) :: GreenR1, GreenR2, GreenA, A
    integer :: flagR

    integer :: i

    !Hermitiano di GreenR2 passato in GreenA
    call zdagacsr(GreenR2,GreenA)

    !Cambio segno a GreenA -> -GreenA
    do i=1,GreenA%nnz
      GreenA%nzval(i)=(-1.0_dp, 0.0_dp)*GreenA%nzval(i)
    enddo

    !Esegue A=GreenR1-GreenA
    call zsumcsr(GreenR1,GreenA,A)

    !A=j(GreenR-GreenA)
    do i=1,A%nnz
      A%nzval(i)=(0.0_dp, 1.0_dp)*(A%nzval(i))
    enddo

    call destroy(GreenA)

    if (flagR.eq.1) then
       call destroy(GreenR2)
    endif

  end subroutine zspectral_csr

  ! DNS version  **********************************************************

  subroutine zspectral_dns(GreenR1,GreenR2,flagR,A)
    implicit none

    type(z_DNS) :: GreenR1, GreenR2, GreenA, A
    integer :: flagR


    !Hermitiano di GreenR2 passato in GreenA
    call zdagger(GreenR2,GreenA)

    call create(A,GreenA%nrow,GreenA%ncol)

    A%val=(0.0_dp,1.0_dp)*(GreenR1%val - GreenA%val)

    call destroy(GreenA)

    if (flagR.eq.1) then
       call destroy(GreenR2)
    endif

  end subroutine zspectral_dns
  !*************************************************************************
  !
  !  Subroutine per il sorting delle matrici csr
  !
  !*************************************************************************

  subroutine zcsort_st(A)

  implicit none

  TYPE(z_cSR) :: A
  INTEGER, DIMENSION(:), ALLOCATABLE :: iwork

   if (.not.A%sorted) then
     CALL log_allocate(iwork,MAX(A%nrow+1,2*A%nnz))
     CALL csort(A%nrow,A%nzval,A%colind,A%rowpnt,iwork,.true.)
     call log_deallocate(iwork)
     A%sorted = .true.
   endif

  end subroutine zcsort_st

  !------------------------------------------------------------------------

  subroutine zcooxcsr_st(coo,sp)

    implicit none

    type(z_EXT_COO) :: coo
    type(z_CSR) :: sp

    integer :: n,nnz,nnz_p

    if ((coo%nrow.ne.sp%nrow).or.(coo%ncol.ne.sp%ncol)) then
        call error_msg('(zcooxcsr_st)',MISMATCH)
    endif

    nnz=coo%nnz
    ! First count primitive non-zero values
    nnz_p=0
    do n=1,nnz
      if(coo%first(n)) nnz_p=nnz_p+1
    enddo

    !Allocate matrix B
    call create(sp,coo%nrow,coo%ncol,nnz_p)

    call zcooxcsr(coo%nrow,coo%nnz,coo%nzval,coo%index_i,coo%index_j,coo%first, &
        sp%nzval,sp%colind,sp%rowpnt)

  end subroutine zcooxcsr_st

  !------------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !  Extended Coordinate   to   Compressed Sparse Row
  !-----------------------------------------------------------------------
  ! converts a matrix that is stored in coordinate format
  !  a, ir, jc into a row general sparse ao, jao, iao format.
  !
  ! on entry:
  !---------
  ! nrow	= dimension of the matrix
  ! nnz	= number of nonzero elements in matrix
  ! a,
  ! ir,
  ! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
  !         nonzero elements of the matrix with a(k) = actual real value of
  ! 	  the elements, ir(k) = its row number and jc(k) = its column
  !	  number. The order of the elements is arbitrary.
  ! lp    = logical flag indicating primitives.
  !
  ! on return:
  !-----------
  ! ir 	is destroyed!
  !
  ! ao, jao, iao = matrix in general sparse matrix format with ao
  ! 	continung the real values, jao containing the column indices,
  !	and iao being the pointer to the beginning of the row,
  !	in arrays ao, jao.
  !
  !Notes:
  !------ This routine is NOT in place.  See coicsr
  !------------------------------------------------------------------------

  subroutine zcooxcsr(nrow,nnz,a,ir,jc,lp,ao,jao,iao)
    integer :: nrow, nnz
    complex(dp)  :: a(*),ao(*),x
    integer     :: ir(*),jc(*),jao(*),iao(*)
    logical     :: lp(*)

    integer :: l, k, j, k0, i, iad

    do k=1,nrow+1
      iao(k) = 0
    enddo
    ! determine row-lengths.
    do k=1, nnz
      if (lp(k)) iao(ir(k)) = iao(ir(k))+1
    enddo
    !starting position of each row..
    k = 1
    do j=1,nrow+1
      k0 = iao(j)
      iao(j) = k
      k = k+k0
    enddo
    ! go through the structure  once more. Fill in primitives.
    do k=1, nnz
      if (lp(k)) then
        i = ir(k)
        j = jc(k)
        x = a(k)

        iad = iao(i)
        ao(iad) = x
        jao(iad) = j
        iao(i) = iad+1
      endif
    enddo
    ! shift back iao
    do j=nrow,1,-1
      iao(j+1) = iao(j)
    enddo
    iao(1) = 1
    ! go through the structure  once more. to add non primitives
    do k=1, nnz
      if (.not.lp(k)) then
        i = ir(k)
        j = jc(k)
        x = a(k)

        do l=iao(i),iao(i+1)-1
          if (jao(l).eq.j) then
            iad=l
            exit
          endif
        enddo

        ao(iad) = ao(iad) + x

      endif
    enddo

    return

  end subroutine zcooxcsr

  !------------- end of coocsr -------------------------------------------
  !-----------------------------------------------------------------------
  subroutine getdiag_csr(A_csr,D_vec)

    !************************************************************************
    !
    !Input:
    !A_csr: matrice di partenza in forma sparsa csr
    !D_vec: vector of dimension nrow
    !
    !****************************************************************************

    implicit none

    type(z_CSR) :: A_csr
    complex(kind=dp), dimension(:) :: D_vec

    integer :: len
    integer, allocatable, dimension(:) :: idiag

    IF(SIZE(D_vec).lt.A_csr%nrow) THEN
       STOP 'Error in getdiag. D_vec has dimension lower than nrow'
    ENDIF

    call log_allocate(idiag,A_csr%nrow)

    call getdia(A_csr%nrow,A_csr%ncol,0,A_csr%nzval,A_csr%colind, &
                 A_csr%rowpnt,len,D_vec,idiag,0)

    call log_deallocate(idiag)

  end subroutine getdiag_csr

  !-----------------------------------------------------------------------
  subroutine rgetdiag_csr(A_csr,D_vec)

    !************************************************************************
    !
    !Input:
    !A_csr: matrice di partenza in forma sparsa csr
    !D_vec: vector of dimension nrow
    !
    !****************************************************************************

    implicit none

    type(r_CSR) :: A_csr
    real(kind=dp), dimension(:) :: D_vec

    integer :: len
    integer, allocatable, dimension(:) :: idiag

    IF(SIZE(D_vec).lt.A_csr%nrow) THEN
       STOP 'Error in getdiag. D_vec has dimension lower than nrow'
    ENDIF

    call log_allocate(idiag,A_csr%nrow)

    call getdia(A_csr%nrow,A_csr%ncol,0,A_csr%nzval,A_csr%colind, &
                 A_csr%rowpnt,len,D_vec,idiag,0)

    call log_deallocate(idiag)

  end subroutine rgetdiag_csr
  !---------------------------------------------

  function ztrace_csr(mat, mask) result(trace)
    type(z_CSR), intent(in) :: mat
    logical, intent(in), optional :: mask(:)
    complex(dp) :: trace

    complex(kind=dp), dimension(:), allocatable :: D_vec

    if (present(mask)) then
       if (size(mask) /= mat%nrow) then
          stop 'Error in ztrace_csr: size(mask) /= nrow'
       end if
    end if

    call log_allocate(D_vec,mat%nrow)

    call getdiag(mat,D_vec)

    if (present(mask)) then
      trace = sum(D_vec, mask)
    else
      trace = sum(D_vec)
    end if

    call log_deallocate(D_vec)

  end function ztrace_csr

!---------------------------------------------

  function ztrace_dns(mat, mask) result(trace)
    type(z_DNS), intent(in) :: mat
    logical, intent(in), optional :: mask(:)
    complex(dp) :: trace

    integer :: i

    trace = (0.0_dp,0.0_dp)
    if (present(mask)) then
       if (size(mask) /= mat%nrow) then
          stop 'Error in ztrace_csr: size(mask) /= nrow'
       end if
       do i = 1,mat%nrow
          if (mask(i)) then
             trace = trace + mat%val(i,i)
          end if
       end do
    else
       do i = 1,mat%nrow
         trace = trace + mat%val(i,i)
       end do
    end if

  end function ztrace_dns

!---------------------------------------------

  function ztrace_arr(mat, mask) result (tr)
    complex(dp), intent(in) :: mat(:,:)
    logical, intent(in), optional :: mask(:)
    complex(dp) :: tr

    integer :: ii
    tr = (0.0_dp, 0.0_dp)
    if (present(mask)) then
       if (size(mask) /= size(mat,1)) then
          stop 'Error in ztrace_csr: size(mask) /= nrow'
       end if
       do ii = 1, size(mat,1)
         if (mask(ii)) then
           tr = tr + mat(ii, ii)
         end if
       end do
    else
       do ii = 1, size(mat,1)
           tr = tr + mat(ii, ii)
       end do
    end if

  end function ztrace_arr

!---------------------------------------------

  function rgetelm_csr(i1,i2,Grm) result(getelm_out)
    type(r_CSR) :: Grm
    integer :: i1, i2, iadd
    real(dp) :: getelm_out

    if (Grm%sorted) then
      getelm_out = getelm(i1,i2,Grm%nzval,Grm%colind,Grm%rowpnt,iadd,.true.)
    else
      getelm_out = getelm(i1,i2,Grm%nzval,Grm%colind,Grm%rowpnt,iadd,.false.)
    endif
  end function rgetelm_csr

  function zgetelm_csr(i1,i2,Grm) result(getelm_out)
    type(z_CSR) :: Grm
    integer :: i1, i2, iadd
    complex(dp) :: getelm_out

    if (Grm%sorted) then
      getelm_out = getelm(i1,i2,Grm%nzval,Grm%colind,Grm%rowpnt,iadd,.true.)
    else
      getelm_out = getelm(i1,i2,Grm%nzval,Grm%colind,Grm%rowpnt,iadd,.false.)
    endif

  end function zgetelm_csr
!---------------------------------------------

  subroutine check_hermitian_dns(M)
    TYPE(z_DNS) :: M

    integer :: ii, jj, file_num

    open(newunit=file_num,file='herm_check.dat')
    do jj = 1, size(M%val,2)
      do ii = jj+1, size(M%val,1)
         if (abs(M%val(ii,jj)-conjg(M%val(jj,ii))) > 1e-10 ) then
            write(file_num,*) 'elements',ii,jj,'non-hermitian:'
            write(file_num,*) M%val(ii,jj)
            write(file_num,*) M%val(jj,ii)
         end if
       end do
    end do
    close(file_num)
  end subroutine check_hermitian_dns

  SUBROUTINE check_hermitian_csr(ham)

    TYPE(z_CSR) :: ham

    ! Local variables:

    INTEGER :: row, p, col, file_num, count
    COMPLEX( dp ) :: matel1, matel2

    call msort(ham)

    open(newunit=file_num,file='herm_check.dat')

    count = 0
    do row = 1, ham%nrow

       do p = ham%rowpnt(row), ham%rowpnt(row+1) - 1 !row+1,n_ham

          col = ham%colind(p)

          matel1 = sprs_element(ham%nzval, ham%colind, ham%rowpnt, row, col)
          matel2 = sprs_element(ham%nzval, ham%colind, ham%rowpnt, col, row)

          if( abs(matel1-conjg(matel2)).gt. 1.d-10 ) then

             count = count + 1
             write(file_num,*) row,col,matel1
             write(file_num,*) col,row,matel2
             write(file_num,*)

          end if

       enddo
    enddo
    if (count .ne. 0) then
       write(*,*) 'Found',count,'wrong elements'
    end if

    close(file_num)

  END SUBROUTINE check_hermitian_csr

  !--------------------------------------------------------------
  ! Returns value of element M(row, col) from sparse matrix
  !
  COMPLEX (dp) FUNCTION sprs_element(M, colind, rowpnt, row, col)

    !--------------------------------------------------
    ! IN data
    INTEGER, INTENT( IN ) :: row, col
    !-------------------------------------------------
    ! OUT data
    INTEGER, DIMENSION(:) :: colind
    INTEGER, DIMENSION(:) :: rowpnt
    COMPLEX (dp), DIMENSION(:) :: M
    !-------------------------------------------------
    INTEGER :: index, ibeg, iend, imid


    index = 0
    ibeg = rowpnt( row )
    iend = rowpnt( row + 1 ) - 1

    DO WHILE (iend.GE.ibeg)

       imid = (ibeg + iend) / 2

       IF ( colind(imid) .eq. col ) THEN
          index = imid
          exit
       END IF

       IF ( colind(imid) .GT. col) then
          iend = imid - 1
       ELSE
          ibeg = imid + 1
       END IF

    END DO

    IF ( index .EQ. 0 ) THEN
       sprs_element = 0.0_dp
    ELSE
       sprs_element = M(index)
    END IF

  END FUNCTION sprs_element
  !-------------------------------------------------

  subroutine error_msg(string,err)
    integer, intent(in) :: err
    character(*), intent(in) :: string

    write(*,*)

    select case(err)
    case(MISMATCH)
       WRITE(*,*) 'ERROR '//trim(string)//' matrices don''t match'
    case(OUTOFBOUND)
       WRITE(*,*) 'ERROR '//trim(string)//' exceeding nnz of destination matrix'
    case(BADINDEX)
       WRITE(*,*) 'ERROR '//trim(string)//' bad indeces'
    case(CONVERR)
       WRITE(*,*) 'ERROR '//trim(string)//' conversion error'
    end select

  end subroutine error_msg

  !------------------------------------------------------------
  !> Convert a CSR in a block dense matrix, only taking the
  !  diagonal and sub/over diag
  subroutine zcsr2blk_sod(Acsr,Ablk,indblk)

        IMPLICIT NONE

    INTEGER :: i
    TYPE(z_CSR) :: Acsr
    INTEGER :: nbl
    TYPE(z_DNS), DIMENSION(:,:) :: Ablk
    INTEGER, DIMENSION(:) :: indblk

    nbl = size(Ablk,1)

    DO i=1,nbl

      CALL extract(Acsr,indblk(i),indblk(i+1)-1,indblk(i),indblk(i+1)-1,Ablk(i,i))

    END DO

    DO i=2,nbl

      CALL extract(Acsr,indblk(i-1),indblk(i)-1,indblk(i),indblk(i+1)-1,Ablk(i-1,i))
      CALL extract(Acsr,indblk(i),indblk(i+1)-1,indblk(i-1),indblk(i)-1,Ablk(i,i-1))

    END DO


  end subroutine zcsr2blk_sod


end module sparsekit_drv
