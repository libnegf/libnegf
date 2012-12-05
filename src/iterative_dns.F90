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


#ifdef __PARDISO
#  undef __LAPACK
#  undef __SUPERLU
#endif
#ifdef __LAPACK
#  undef __PARDISO
#  undef __SUPERLU
#endif
#ifdef __SUPERLU
#  undef __PARDISO
#  undef __LAPACK
#endif




MODULE iterative_dns

  USE ln_precision
  USE ln_constants, only : pi
  USE ln_allocation
  USE mat_def
  USE sparsekit_drv
  USE inversions
  USE ln_structure, only : TStruct_Info
  USE lib_param, only : MAXNCONT, Tnegf
  USE outmatrix, only : outmat_c, inmat_c, direct_out_c, direct_in_c 
  USE clock
  USE transform

  IMPLICIT NONE
  private

  public :: allocate_gsmr_dns
  public :: allocate_gsml_dns
  public :: allocate_Gr_dns  

  public :: deallocate_gsmr_dns
  public :: deallocate_gsml_dns
  public :: deallocate_Gr_dns

  public :: tunneling_dns
  public :: tun_and_dos

  public :: compGreen
 
  public :: calls_eq_mem_dns
  public :: calls_neq_mem_dns
  public :: calls_neq_ph

  public :: sigma_ph_n
  public :: sigma_ph_p
  public :: sigma_ph_r
  public :: sigma_ph_r_z
  public :: check_sigma_ph_r

  public :: sub_ESH_dns
  public :: rebuild_dns
  public :: Make_gsmr_mem_dns
  public :: Make_gsml_mem_dns
  public :: Make_Grdiag_mem_dns
  public :: Make_Grcol_mem_dns
  !public :: Make_Spectral_mem_dns
  public :: Make_Gr_mem_dns  ! New. Best.
  public :: Make_Gn_mem_dns
  !public :: Outer_A_mem_dns
  public :: Outer_Gr_mem_dns
  public :: Outer_Gn_mem_dns

  public :: complete_sigma_ph_r

  LOGICAL, PARAMETER :: debug=.false. 
  !Dropout value
  REAL(dp), PARAMETER :: drop=1e-20

  TYPE(z_DNS), DIMENSION(:), ALLOCATABLE :: gsmr
  TYPE(z_DNS), DIMENSION(:), ALLOCATABLE :: gsml
  TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Gr

CONTAINS

  !****************************************************************************
  ! 
  ! Driver for computing Equilibrium Green Retarded (Outer blocks included)
  ! writing on memory
  !
  !****************************************************************************

  SUBROUTINE calls_eq_mem_dns(pnegf,E,SelfEneR,Tlc,Tcl,gsurfR,A,struct,outer)

    !****************************************************************************
    !
    !Input
    !H: sparse matrix contaning Device Hamiltonian
    !S: sparse matrix containing Device Overlap
    !SelfEneR: sparse matrices array containing contacts Self Energy
    !Tlc: sparse matrices array containing contacts-device interaction blocks (ES-H)
    !Tcl: sparse matrices array containing device-contacts interaction blocks (ES-H)
    !gsurfR: sparse matrices array containing contacts surface green
    !ref: reference contact 
    !lower_outer: optional parameter. If defined (and true), Aout contains also 
    !the lower outer parts (needed for K-points calculations)
    !
    !Output:
    !A: Spectral function (Device + Contacts overlap regions -> effective conductor)
    !
    !*****************************************************************************

    IMPLICIT NONE

    !In/Out
    TYPE(Tnegf), pointer :: pnegf
    COMPLEX(dp) :: E
    TYPE(z_DNS), DIMENSION(:) :: SelfEneR
    TYPE(z_DNS), DIMENSION(:) :: Tlc, Tcl, gsurfR
    TYPE(z_CSR) :: A
    TYPE(Tstruct_info) :: struct
    INTEGER :: outer

    !Work
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: ESH
    TYPE(z_CSR) :: ESH_tot, Ain
    INTEGER :: i,ierr, nbl, ncont
    INTEGER, DIMENSION(:), POINTER :: cblk, indblk

    nbl = struct%num_PLs
    ncont = struct%num_conts
    cblk => struct%cblk
    indblk => struct%mat_PL_start

    !Costruiamo matrice densa ESH_tot
    CALL prealloc_sum(pnegf%HM,pnegf%SM,(-1.d0, 0.d0),E,ESH_tot)

    !Costruiamo l'array di dense ESH
    !Allocazione dell'array di dense ESH
    ALLOCATE(ESH(nbl,nbl),stat=ierr)
    IF (ierr.EQ.1) STOP 'Error in ESH allocation'

    CALL sub_ESH_dns(ESH_tot,ESH,indblk)

    CALL destroy(ESH_tot)

    DO i=1,ncont
       !call create(SelfEner_d, SelfEner(i)%nrow, SelfEner(i)%ncol)
       !call csr2dns(SelfEneR(i),SelfEneR_d)
       ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
       !call destroy(SelfEneR_d)
    ENDDO

    if (pnegf%elph%numselmodes.gt.0) then
       ! Reload and add Sigma_ph_r to ESH
       call add_sigma_ph_r(pnegf, ESH, pnegf%elph%scba_iter)
    endif

    !Allocazione delle gsmr
    ALLOCATE(gsmr(nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate gsmr'

    !Chiamata di Make_gsmr_mem

    CALL Make_gsmr_mem_dns(ESH,nbl,2)


    !Allocazione delle Gr
    ALLOCATE(Gr(nbl,nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Gr'

    CALL Make_Grdiag_mem_dns(ESH,indblk)


    !Distruzione delle gsmall
    do i=2,nbl 
       call destroy(gsmr(i))
    enddo
    DEALLOCATE(gsmr)

    SELECT CASE (outer)
    CASE(0)
    CASE(1)
       CALL Outer_Gr_mem_dns(Tlc,Tcl,gsurfR,struct,.FALSE.,A)   
    CASE(2)
       CALL Outer_Gr_mem_dns(Tlc,Tcl,gsurfR,struct,.TRUE.,A) 
    END SELECT

    !CALL Make_Spectral_mem(struct,0,A)
    !if (debug) write(*,*) 'Compute GreenR'   
    !Nota: Il flag 0 indica che non teniamo le Gr calcolate 
    !E' necessario metterlo a 1 se dopo il calcolo della spectral 
    !effettuiamo il calcolo dei contributi aggiunti della Gn=-iG<
    call Make_Gr_mem_dns(pnegf%SM,nbl,indblk,Ain)

    if (pnegf%elph%numselmodes.gt.0) then
      ! save diagonal blocks of Gn = -i G<
      DO i = 1, nbl
         call write_blkmat(Gr(i,i),pnegf%scratch_path,'G_r_',i,i,pnegf%iE)
      ENDDO
    endif   

    !Distruzione dell'array Gr
    CALL destroy_blk(Gr)
    DEALLOCATE(Gr)

    CALL destroy_ESH(ESH)
    DEALLOCATE(ESH)

    !Concatenazioone di A in Aout

    SELECT CASE(outer)
    CASE(0)
       call clone(Ain,A)
    CASE(1:2)
       call concat(A,Ain,1,1)
    END SELECT

    call destroy(Ain)


  END SUBROUTINE calls_eq_mem_dns




  !****************************************************************************
  !
  ! Driver for computing G_n contributions due to all contacts but reference:
  !
  !   Sum   [f_j(E)-f_r(E)] Gr Gam_j Ga
  !   j!=r
  !
  ! NOTE: The subroutine assumes that 
  !
  !****************************************************************************

  SUBROUTINE calls_neq_mem_dns(pnegf,E,SelfEneR,Tlc,Tcl,gsurfR,struct,frm,Glout,out)

    !****************************************************************************
    !
    !Input
    !H: sparse matrix contaning Device Hamiltonian
    !S: sparse matrix containing Device Overlap
    !SelfEneR: sparse matrices array containing contacts Self Energy
    !Tlc: sparse matrices array containing contacts-device interaction blocks (ES-H)
    !Tcl: sparse matrices array containing device-contacts interaction blocks (ES-H)
    !gsurfR: sparse matrices array containing contacts surface green
    !frm: array containing Fermi distribution values for all contacts
    !ref: reference contact excluded from summation
    !
    !Output:
    !Aout: G_n contributions (Device + Contacts overlap regions -> effective conductor)
    !
    !*****************************************************************************


    IMPLICIT NONE

    !In/Out
    TYPE(Tnegf), pointer :: pnegf
    TYPE(z_CSR) :: Glout
    TYPE(z_DNS), DIMENSION(:) :: SelfEneR, gsurfR, Tlc, Tcl
    !TYPE(z_DNS) :: SelfEner_d 
    REAL(dp) :: E
    TYPE(Tstruct_info) :: struct
    REAL(dp), DIMENSION(:) :: frm
    INTEGER :: out

    !Work
    INTEGER :: ref
    COMPLEX(dp) :: Ec
    INTEGER :: i,ierr,i1,ncont,nbl, lbl
    INTEGER, DIMENSION(:), POINTER :: cblk, indblk
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: ESH
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Gn
    TYPE(z_CSR) :: ESH_tot, Gl
    LOGICAL :: destr, mask(MAXNCONT)

    nbl = struct%num_PLs
    ncont = struct%num_conts
    indblk => struct%mat_PL_start
    cblk => struct%cblk
    ref = pnegf%refcont

    !if (debug) write(*,*) '----------------------------------------------------'
    !if (debug) call writeMemInfo(6)
    !if (debug) call writePeakInfo(6)
    !if (debug) write(*,*) '----------------------------------------------------'
    Ec=cmplx(E,0.d0,dp)

    !Costruiamo la matrice sparsa ESH
    CALL prealloc_sum(pnegf%HM,pnegf%SM,(-1.d0, 0.d0),Ec,ESH_tot)

    !Costruiamo l'array di sparse ESH
    !Allocazione dell'array di sparse ESH
    ALLOCATE(ESH(nbl,nbl),stat=ierr)
    IF (ierr.EQ.1) STOP 'Error in ESH allocation'

    CALL sub_ESH_dns(ESH_tot,ESH,indblk)

    CALL destroy(ESH_tot)

    DO i=1,ncont
       !call create(SelfEner_d, SelfEner(i)%nrow, SelfEner(i)%ncol)
       !call csr2dns(SelfEneR(i),SelfEneR_d)
       ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
       !call destroy(SelfEneR_d)
    ENDDO

    !if (debug) write(*,*) '----------------------------------------------------'
    !if (debug) call writeMemInfo(6)
    !if (debug) call writePeakInfo(6)
    !if (debug) write(*,*) '----------------------------------------------------'

    !Allocazione delle gsmr
    ALLOCATE(gsmr(nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate gsmr'
    ALLOCATE(gsml(nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate gsml'

    !Chiamata di Make_gsmr_mem
    CALL Make_gsmr_mem_dns(ESH,nbl,2)

    !Chiamata di Make_gsml_mem solo per i blocchi 1..lbl dove  
    ! lbl = maxval(cblk,mask) - 2 
    mask = .true.
    mask(ref) = .false. 
    lbl = maxval(cblk(1:ncont),mask(1:ncont)) - 2

    if( ncont.gt.1 ) then
       CALL Make_gsml_mem_dns(ESH,1,lbl)    
    endif
    !if (debug) write(*,*) '----------------------------------------------------'
    !if (debug) call writeMemInfo(6)
    !if (debug) call writePeakInfo(6)
    !if (debug) write(*,*) '----------------------------------------------------'

    !Allocazione delle Gr
    ALLOCATE(Gr(nbl,nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Gr'

    CALL Make_Grdiag_mem_dns(ESH,indblk)

    !if (debug) write(*,*) '----------------------------------------------------'
    !if (debug) call writeMemInfo(6)
    !if (debug) call writePeakInfo(6)
    !if (debug) write(*,*) '----------------------------------------------------'


    !Chiamata di Make_Grcol_mem per i contatti necessari 
    DO i=1,ncont
       IF (i.NE.ref) THEN
          CALL Make_Grcol_mem_dns(ESH,cblk(i),indblk)
       ENDIF
    ENDDO

    !Distruzione delle gsmall
    do i=2,nbl 
       call destroy(gsmr(i))
    enddo
    DEALLOCATE(gsmr)
    DO i=1,nbl-1 
       IF(gsml(i)%nrow.NE.0) CALL destroy(gsml(i))
    ENDDO
    DEALLOCATE(gsml)

    !if (debug) write(*,*) '----------------------------------------------------'
    !if (debug) call writeMemInfo(6)
    !if (debug) call writePeakInfo(6)
    !if (debug) write(*,*) '----------------------------------------------------'

    !Calcolo degli outer blocks 
    !if (debug) write(*,*) 'Compute Outer_G_n' 
    SELECT CASE (out)
    CASE(0)
    CASE(1)
       CALL Outer_Gn_mem_dns(Tlc,gsurfR,SelfEneR,struct,frm,ref,.false.,Glout)
    CASE(2)
       CALL Outer_Gn_mem_dns(Tlc,gsurfR,SelfEneR,struct,frm,ref,.true.,Glout)
    END SELECT

    !if (debug) write(*,*) '----------------------------------------------------'
    !if (debug) call writePeakInfo(6)
    !if (debug) write(*,*) '----------------------------------------------------'

    !Calcolo della G_n nel device
    !if (debug) write(*,*) 'Compute G_n' 
    !Allocazione degli array di sparse
    ALLOCATE(Gn(nbl,nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Gn(nbl,nbl)'
 
    call init_blkmat(Gn,ESH)

    CALL Make_Gn_mem_dns(ESH,SelfEneR,frm,ref,struct,Gn)
    
    call blk2csr(Gn,struct,pnegf%SM,Gl)

    DEALLOCATE(Gn)

    !Distruzione dell'array Gr
    CALL destroy_blk(Gr)
    DEALLOCATE(Gr)

    CALL destroy_ESH(ESH)
    DEALLOCATE(ESH)

    !Concatenazione di Gl in Glout
    SELECT CASE(out)
    CASE(0)
       call clone(Gl,Glout)
    CASE(1:2)
       call concat(Glout,Gl,1,1)
    END SELECT

    call destroy(Gl)
    !if (debug) write(*,*) '----------------------------------------------------'
    !if (debug) call writeMemInfo(6)
    !if (debug) call writePeakInfo(6)
    !if (debug)  write(*,*) '----------------------------------------------------'

  END SUBROUTINE calls_neq_mem_dns

!****************************************************************************
!
! Driver for computing G_n = -iG< contributions including el-ph interactions
!
!   Sum   f_j(E) Gr Gam_j Ga +   Gr Sigma_ph< Ga
!    j
!
! NOTE: The subroutine assumes that
!
!****************************************************************************

SUBROUTINE calls_neq_ph(pnegf,E,SelfEneR,Tlc,Tcl,gsurfR,frm,Glout,out)

!****************************************************************************
!
!Input
!H: sparse matrix contaning Device Hamiltonian
!S: sparse matrix containing Device Overlap
!SelfEneR: sparse matrices array containing contacts Self Energy
!Tlc: sparse matrices array containing contacts-device interaction blocks (ES-H)
!Tcl: sparse matrices array containing device-contacts interaction blocks (ES-H)
!gsurfR: sparse matrices array containing contacts surface green
!frm: array containing Fermi distribution values for all contacts
!ref: reference contact excluded from summation
!
!Output:
!Aout: G_n contributions (Device + Contacts overlap regions -> effective conductor)
!
!*****************************************************************************
IMPLICIT NONE

!In/Out
TYPE(Tnegf), pointer :: pnegf
TYPE(z_CSR) :: Glout
TYPE(z_DNS), DIMENSION(:) :: SelfEneR, gsurfR, Tlc, Tcl
REAL(dp) :: E
REAL(dp), DIMENSION(:) :: frm
INTEGER :: out

!Work
COMPLEX(dp) :: Ec
INTEGER :: i,ierr,i1,ncont,nbl, lbl
INTEGER :: ref, iter
INTEGER, DIMENSION(:), POINTER :: cblk, indblk
TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: ESH, Gn, Gp
TYPE(z_CSR) :: ESH_tot, Gl
LOGICAL :: destr, mask(MAXNCONT)
TYPE(Tstruct_info) :: struct
REAL(dp) :: Iloc
REAL(dp), DIMENSION(:), allocatable :: cfrm

struct = pnegf%str
nbl = struct%num_PLs
ncont = struct%num_conts
indblk => struct%mat_PL_start
cblk => struct%cblk
ref = pnegf%refcont
iter = pnegf%elph%scba_iter

Ec=cmplx(E,0.d0,dp)

!if (debug) write(*,*) '----------------------------------------------------'
!if (debug) call writeMemInfo(6)
!if (debug) call writePeakInfo(6)
!if (debug) write(*,*) '----------------------------------------------------'

 !Costruiamo la matrice sparsa ESH
 CALL prealloc_sum(pnegf%HM,pnegf%SM,(-1.d0, 0.d0),Ec,ESH_tot)

 !Costruiamo l'array di sparse ESH
 !Allocazione dell'array di sparse ESH
 ALLOCATE(ESH(nbl,nbl),stat=ierr)
 IF (ierr.EQ.1) STOP 'Error in ESH allocation'

 CALL sub_ESH_dns(ESH_tot,ESH,indblk)

 call destroy(ESH_tot)

 DO i=1,ncont
       ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
 ENDDO

 ! Reload and add Sigma_ph_r to ESH
 call add_sigma_ph_r(pnegf, ESH, iter)

 !Allocazione delle gsmr
 ALLOCATE(gsmr(nbl),stat=ierr)
 IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate gsmr'
 ALLOCATE(gsml(nbl),stat=ierr)
 IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate gsml'

 !Chiamata di Make_gsmr_mem
 CALL Make_gsmr_mem_dns(ESH,nbl,2)

 !Chiamata di Make_gsml_mem solo per i blocchi 1..lbl dove  
 ! lbl = maxval(cblk,mask) - 2 
 mask = .true.
 mask(ref) = .false. 
 lbl = maxval(cblk(1:ncont),mask(1:ncont)) - 2

 if( ncont.gt.1 ) then
  CALL Make_gsml_mem_dns(ESH,1,lbl)    
 endif
 !if (debug) write(*,*) '----------------------------------------------------'
 !if (debug) call writeMemInfo(6)
 !if (debug) call writePeakInfo(6)
 !if (debug) write(*,*) '----------------------------------------------------'

 !Allocazione delle Gr
 ALLOCATE(Gr(nbl,nbl),stat=ierr)
 IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Gr'

 CALL Make_Grdiag_mem_dns(ESH,indblk)

!if (debug) write(*,*) '----------------------------------------------------'
!if (debug) call writeMemInfo(6)
!if (debug) call writePeakInfo(6)
!if (debug) write(*,*) '----------------------------------------------------'


 !Chiamata di Make_Grcol_mem per i contatti necessari 
 !DO i=1,ncont 
 !  IF (i.NE.ref) THEN
 !    CALL Make_Grcol_mem_dns(ESH,cblk(i),indblk)
 !  ENDIF
 !ENDDO

 !print*
 !print*,'Column-blocks'
 DO i=1,nbl
   CALL Make_Grcol_mem_dns(ESH,i,indblk)
 ENDDO  


 !Distruzione delle gsmall
 do i=2,nbl 
    call destroy(gsmr(i))
 enddo
 DEALLOCATE(gsmr)
 DO i=1,nbl-1 
   IF(gsml(i)%nrow.NE.0) CALL destroy(gsml(i))
 ENDDO
 DEALLOCATE(gsml)

!if (debug) write(*,*) '----------------------------------------------------'
!if (debug) call writeMemInfo(6)
!if (debug) call writePeakInfo(6)
!if (debug) write(*,*) '----------------------------------------------------'

 !Calcolo degli outer blocks 
 !if (debug) write(*,*) 'Compute Outer_G_n' 
 SELECT CASE (out)
 CASE(0)
 CASE(1)
   CALL Outer_Gn_mem_dns(Tlc,gsurfR,SelfEneR,pnegf%str,frm,ref,.false.,Glout)
 CASE(2)
   CALL Outer_Gn_mem_dns(Tlc,gsurfR,SelfEneR,pnegf%str,frm,ref,.true.,Glout)
 END SELECT

!if (debug) write(*,*) '----------------------------------------------------'
!if (debug) call writePeakInfo(6)
!if (debug) write(*,*) '----------------------------------------------------'

 !Calcolo della G_n nel device
 !if (debug) write(*,*) 'Compute G_n' 
 !Allocazione degli array di sparse
 ALLOCATE(Gn(nbl,nbl),stat=ierr)
 IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Gn(nbl,nbl)'

 call init_blkmat(Gn,ESH)

 CALL Make_Gn_mem_dns(ESH,SelfEneR,frm,ref,struct,Gn)

 call Make_Gn_ph(pnegf,ESH,iter,Gn)

 !print*
 ! save diagonal blocks of G_n
 DO i = 1, nbl
   !print*,'(G_n) G_n',minval(abs(Gn(i,i)%val)), maxval(abs(Gn(i,i)%val))
   call write_blkmat(Gn(i,i),pnegf%scratch_path,'G_n_',i,i,pnegf%iE)
   !print*,'(G_r) G_r',minval(abs(Gn(i,i)%val)),maxval(abs(Gr(i,i)%val))
   call write_blkmat(Gr(i,i),pnegf%scratch_path,'G_r_',i,i,pnegf%iE)
 ENDDO


 !open(5001,file=trim(pnegf%out_path)//'Lcurr.dat',position='APPEND')
 !write(5001,'(ES15.5)',advance='NO') real(Ec)
 !do i = 1, nbl
 !  do i1 = 1, ESH(i,i)%nrow-1
 !     Iloc = 2*aimag(Glsub(i,i)%val(i1,i1+1))*ESH(i,i)%val(i1+1,i1)
 !     write(5001,'(ES15.5)',advance='NO') Iloc
 !  enddo
 !enddo    
 !write(5001,*)
 !close(5001)

 ! the following destroys Gn
 call blk2csr(Gn,pnegf%str,pnegf%SM,Gl)

 DEALLOCATE(Gn)
 
 ! Temporary -----------------------------------------------------------
 ! Computation of G_p for Sigma_p
 !
    allocate(cfrm(10))
    cfrm = 1.0_dp - frm
    cfrm(ref) = 0.0_dp
 !
    ALLOCATE(Gp(nbl,nbl),stat=ierr)
 !
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Gn(nbl,nbl)'
    call init_blkmat(Gp,ESH)
 !   
    CALL Make_Gn_mem_dns(ESH,SelfEneR,cfrm,ref,struct,Gp)
    deallocate(cfrm)
 !   
    call Make_Gp_ph(pnegf,ESH,iter,Gp)
 !   
    ! save diagonal blocks of G_p
    DO i = 1, nbl
      !print*
      !print*,'(G_p) G_p',minval(abs(Gn(i,i)%val)),maxval(abs(Gn(i,i)%val))
      call write_blkmat(Gp(i,i),pnegf%scratch_path,'G_p_',i,i,pnegf%iE)
    ENDDO
 !   
    call destroy_blk(Gp)
    DEALLOCATE(Gp)
 ! ------------------------------------------------------------------------

 !Distruzione dell'array Gr
 CALL destroy_blk(Gr)
 DEALLOCATE(Gr)

 CALL destroy_ESH(ESH)
 DEALLOCATE(ESH)

 !Concatenazione di Gl in Glout
 SELECT CASE(out)
 CASE(0)
    call clone(Gl,Glout)
 CASE(1:2)
    call concat(Glout,Gl,1,1)
 END SELECT

 call destroy(Gl)
 !if (debug) write(*,*) '----------------------------------------------------'
 !if (debug) call writeMemInfo(6)
 !if (debug) call writePeakInfo(6)
 !if (debug)  write(*,*) '----------------------------------------------------'

  END SUBROUTINE calls_neq_ph

  !**********************************************************************
  SUBROUTINE init_blkmat(Matrix,S)
    type(z_DNS), dimension(:,:) :: Matrix,S

    integer :: nbl, j

    nbl = SIZE(Matrix,1)

    call create(Matrix(1,1),S(1,1)%nrow,S(1,1)%ncol)      
    Matrix(1,1)%val=(0.d0,0.d0)
    DO j=2,nbl-1
       call create(Matrix(j-1,j),S(j-1,j)%nrow,S(j-1,j)%ncol)      
       Matrix(j-1,j)%val=(0.d0,0.d0)
       call create(Matrix(j,j),S(j,j)%nrow,S(j,j)%ncol)      
       Matrix(j,j)%val=(0.d0,0.d0)
       call create(Matrix(j,j-1),S(j,j-1)%nrow,S(j,j-1)%ncol)      
       Matrix(j,j-1)%val=(0.d0,0.d0)
    ENDDO
    IF (nbl.gt.1) then
       call create(Matrix(nbl,nbl),S(1,1)%nrow,S(1,1)%ncol)      
       Matrix(nbl,nbl)%val=(0.d0,0.d0)
    ENDIF

  END SUBROUTINE init_blkmat  
  !**********************************************************************
  SUBROUTINE destroy_blk(M)
    type(z_DNS), DIMENSION(:,:) :: M
    integer :: i, i1, nbl

    nbl=size(M,1)

    DO i=1,nbl
       DO i1=1,nbl
          IF (ALLOCATED(M(i,i1)%val)) THEN
             CALL destroy(M(i,i1))
             !if (debug) WRITE(*,*) 'Deallocating Gr out of Make_Gn_mem'
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE destroy_blk

  !**********************************************************************
  SUBROUTINE destroy_ESH(ESH)

    integer :: i, nbl
    type(z_DNS), dimension(:,:) :: ESH

    nbl=size(ESH,1)

    DO i=1,nbl
       CALL destroy(ESH(i,i))
    ENDDO
    DO i=2,nbl
       CALL destroy(ESH(i-1,i))
       CALL destroy(ESH(i,i-1))
    ENDDO

  END SUBROUTINE destroy_ESH

  !**********************************************************************
  !
  !  Divides sparse matrix ES-H in the device region in a sparse matrices
  !  array ESH(nbl,nbl)  (needs global variable indblk)
  !
  !**********************************************************************

  SUBROUTINE sub_ESH_dns(ESH_tot,ESH,indblk)

    !**********************************************************************
    !Input:
    !ESH_tot: sparse matrix ES-H related to device
    !
    !Output:
    !ESH(nbl,nbl): dense matrix array -> single matrices allocated 
    !              internally, array ESH(nbl,nbl) allocated externally
    !**********************************************************************

    IMPLICIT NONE 

    INTEGER :: i
    TYPE(z_CSR) :: ESH_tot
    INTEGER :: nbl
    TYPE(z_DNS), DIMENSION(:,:) :: ESH
    INTEGER, DIMENSION(:) :: indblk

    nbl = size(ESH,1)

    DO i=1,nbl

       CALL extract(ESH_tot,indblk(i),indblk(i+1)-1,indblk(i),indblk(i+1)-1,ESH(i,i))

    ENDDO

    DO i=2,nbl

       CALL extract(ESH_tot,indblk(i-1),indblk(i)-1,indblk(i),indblk(i+1)-1,ESH(i-1,i))
       CALL extract(ESH_tot,indblk(i),indblk(i+1)-1,indblk(i-1),indblk(i)-1,ESH(i,i-1))

    ENDDO

  END SUBROUTINE sub_ESH_dns

  !***********************************************************************
  !
  !  Reloads the retarded elph self-energy and add it to ES-H
  !
  !***********************************************************************
  SUBROUTINE add_sigma_ph_r(pnegf, ESH, iter)

     TYPE(Tnegf), pointer :: pnegf
     TYPE(z_DNS), DIMENSION(:,:) :: ESH
     INTEGER :: iter

     TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_ph_r
     INTEGER :: n, nbl, nrow, ierr

     nbl = pnegf%str%num_PLs
    
     ALLOCATE(Sigma_ph_r(nbl,nbl),stat=ierr)
     IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_r'

     if (pnegf%elph%diagonal) THEN

        DO n = 1, nbl

           nrow = ESH(n,n)%nrow

           call create(Sigma_ph_r(n,n), nrow, nrow)

           Sigma_ph_r(n,n)%val = (0.0_dp, 0.0_dp)
           if (iter .gt. 0) then
              call read_blkmat(Sigma_ph_r(n,n),pnegf%scratch_path,'Sigma_ph_r_',n,n,pnegf%iE)
           else
              call write_blkmat(Sigma_ph_r(n,n),pnegf%scratch_path,'Sigma_ph_r_',n,n,pnegf%iE)
           endif

           ESH(n,n)%val = ESH(n,n)%val - Sigma_ph_r(n,n)%val

           call destroy(Sigma_ph_r(n,n))

         END DO

     ELSE

     ENDIF

     DEALLOCATE(Sigma_ph_r)

  END SUBROUTINE add_sigma_ph_r 

  !***********************************************************************
  !
  !  Construct a sparse matrix starting from the sparse matrices array
  !  related to the blocks  
  !
  !***********************************************************************

  SUBROUTINE rebuild_dns(Atot,A,n,indb)

    !***********************************************************************
    !Input:
    !A: sparse matrices array
    !n: A dimension (A(n,n))
    !indb: blocks indeces array (contains position of first elements of diagonal
    !      blocks, included an additional ending address "last row+1")
    !
    !Output:
    !Atot: total sparse matrix (allocated internally)
    !***********************************************************************  

    IMPLICIT NONE

    !In/Out
    TYPE(z_CSR) :: Atot  !Allocato internamente
    INTEGER :: n
    INTEGER, DIMENSION(n+1) :: indb
    TYPE(z_cSR), DIMENSION(n,n) :: A

    !Work
    INTEGER :: i,j,Atot_nrow,i1,j1

    Atot_nrow=indb(n+1)-1

    CALL create(Atot,Atot_nrow,Atot_nrow,0)
    Atot%rowpnt(:)=1

    DO i=1,n
       DO j=1,n

          IF (A(i,j)%nrow.GT.0) THEN
             i1=indb(i)
             j1=indb(j)
             CALL concat(Atot,A(i,j),i1,j1)
          ENDIF

       ENDDO
    ENDDO

  END SUBROUTINE rebuild_dns

  !***********************************************************************
  !
  !  g_small right (gsmr) calculation - write on memory
  !
  !***********************************************************************

  SUBROUTINE Make_gsmr_mem_dns(ESH,sbl,ebl)

    !***********************************************************************
    !Input:
    !ESH: sparse matrices array ESH(nbl,nbl)
    !
    !nbl (number of layers)  global needed   (indblk not anymore)
    !
    !Output:
    !sparse matrices array global variable gsmr(nbl) is available in memory
    !single blocks are allocated internally, array Gr(nbl,nbl) 
    !must be allocated externally
    !***********************************************************************

    IMPLICIT NONE

    !In/Out
    TYPE(z_DNS), DIMENSION(:,:) :: ESH
    integer :: sbl,ebl                            ! start block, end block

    !Work
    !TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: INV
    TYPE(z_DNS) :: work1, work2
    INTEGER :: nrow
    INTEGER :: i, nbl


    nbl = size(ESH,1)

    if (sbl.lt.ebl) return

    !***
    !gsmr(sbl)
    !***
    
    if (nbl.gt.1) then

       nrow=ESH(sbl,sbl)%nrow    !indblk(nbl+1)-indblk(nbl)

       call create(gsmr(sbl),nrow,nrow)

       call compGreen(gsmr(sbl),ESH(sbl,sbl),nrow)

    endif

    !***
    !gsmr(sbl-1):gsmr(ebl)
    !***

    DO i=sbl-1,ebl,(-1)

       nrow=ESH(i,i)%nrow          !indblk(i+1)-indblk(i)
       CALL prealloc_mult(ESH(i,i+1),gsmr(i+1),(-1.d0, 0.d0),work1)
       CALL prealloc_mult(work1,ESH(i+1,i),work2)

       CALL destroy(work1)

       CALL prealloc_sum(ESH(i,i),work2,work1)

       CALL destroy(work2)

       CALL create(gsmr(i),nrow,nrow)

       call compGreen(gsmr(i),work1,nrow)

       CALL destroy(work1)

    ENDDO

    if (debug) then
       WRITE(*,*)
       WRITE(*,*) '********************'
       WRITE(*,*) 'Make_gsmr_mem done'
       WRITE(*,*) '********************'
    endif
 

  END SUBROUTINE Make_gsmr_mem_dns

 


  !***********************************************************************
  !
  !  g_small left (gsml) calculation - write on memory
  !
  !***********************************************************************

  SUBROUTINE Make_gsml_mem_dns(ESH,sbl,ebl)

    !***********************************************************************
    !Input:
    !ESH: sparse matrices array ESH(nbl,nbl)
    !
    !nbl (number of layers), indblk(nbl+1) global variables needed
    !
    !Output:
    !sparse matrices array global variable gsml(nbl) is available in memory
    !single blocks are allocated internally, array Gr(nbl,nbl) 
    !must be allocated externally
    !***********************************************************************


    IMPLICIT NONE

    !In/Out
    TYPE(z_DNS), DIMENSION(:,:) :: ESH
    integer :: sbl,ebl                       ! start block, end block

    !Work
    TYPE(z_DNS) :: work1, work2
    INTEGER :: nrow
    INTEGER :: i, nbl
   !TYPE(z_DNS) :: INV(sbl,sbl)

    if (sbl.gt.ebl) return

    !***
    !gsml(sbl)
    !***
    nbl = size(ESH,1)
    nrow=ESH(sbl,sbl)%nrow  !indblk(2)-indblk(1)

    CALL create(gsml(sbl),nrow,nrow)

    call compGreen(gsml(sbl),ESH(sbl,sbl),nrow)


    DO i=sbl+1,ebl

       nrow=ESH(i,i)%nrow   !indblk(i+1)-indblk(i)

       CALL prealloc_mult(ESH(i,i-1),gsml(i-1),(-1.d0, 0.d0),work1)

       CALL prealloc_mult(work1,ESH(i-1,i),work2)

       CALL destroy(work1)

       CALL prealloc_sum(ESH(i,i),work2,work1)

       CALL destroy(work2)

       call create(gsml(i),nrow,nrow)

       call compGreen(gsml(i),work1,nrow)

       CALL destroy(work1)

    ENDDO

    if (debug) then
       WRITE(*,*) '********************'
       WRITE(*,*) 'Make_gsml_mem done'
       WRITE(*,*) '********************'
    endif

  END SUBROUTINE Make_gsml_mem_dns





  !***********************************************************************
  !
  !  Diagonal, Subdiagonal, Superdiagonal blocks of Green Retarded 
  !  Gr(nbl,nbl) - writing on memory
  !
  !*********************************************************************** 

  SUBROUTINE Make_Grdiag_mem_dns(ESH,indblk,mybls)

    !***********************************************************************
    !Input:
    !ESH: dense matrices array ESH(nbl,nbl)
    !
    !global variables needed: nbl (number of layers), indblk(nbl+1),
    !gsmr(:) 
    !
    !Output:
    !sparse matrices array global variable Gr(nbl,nbl) is available in 
    !memory - single blocks are allocated internally, array Gr(nbl,nbl) 
    !must be allocated externally
    !***********************************************************************

    IMPLICIT NONE 

    !In/Out
    TYPE(z_DNS), DIMENSION(:,:) :: ESH
    INTEGER, DIMENSION(:), POINTER :: indblk
    Integer, optional :: mybls

    !Work
    INTEGER :: i,nrow,nbl,mybl2
    TYPE(z_DNS) :: work1, work2, work3

    !***
    !Gr(1,1)
    !***
 
    nbl = size(ESH,1)
    nrow=indblk(2)-indblk(1)  
    if(.not.present(mybls)) then
       mybl2 = nbl
    else
       mybl2 = mybls
    endif

    if(nbl.gt.1) then
       CALL prealloc_mult(ESH(1,2),gsmr(2),work1)
       CALL prealloc_mult(work1,ESH(2,1),work2)

       CALL destroy(work1)

       CALL prealloc_sum(ESH(1,1),work2,(-1.d0, 0.d0),work1)
       
       CALL destroy(work2)

    else

       call create(work1,ESH(1,1)%nrow,ESH(1,1)%ncol)
       work1%val = ESH(1,1)%val

    endif
     
    CALL create(Gr(1,1),nrow,nrow)

    call compGreen(Gr(1,1),work1,nrow)

    CALL destroy(work1)

    !***
    !Diagonal, Subdiagonal and Superdiagonal blocks
    !***

    DO i=2,mybl2

       nrow=indblk(i+1)-indblk(i)

       !***
       !Diagonal blocks
       !***
       CALL prealloc_mult(gsmr(i),ESH(i,i-1),work1)
       CALL prealloc_mult(work1,Gr(i-1,i-1),work2)

       CALL destroy(work1)
       CALL prealloc_mult(ESH(i-1,i),gsmr(i),work3)
       CALL prealloc_mult(work2,work3,work1)

       CALL destroy(work2)

       CALL prealloc_sum(gsmr(i),work1,Gr(i,i))

       CALL destroy(work1) 

       !***
       !Superdiagonal blocks
       !***
       CALL prealloc_mult(Gr(i-1,i-1),work3,(-1.d0, 0.d0),Gr(i-1,i))

       CALL destroy(work3)

       !***
       !Subdiagonal blocks
       !***
       CALL prealloc_mult(gsmr(i),ESH(i,i-1),(-1.d0, 0.d0),work1)
       CALL prealloc_mult(work1,Gr(i-1,i-1),Gr(i,i-1))

       CALL destroy(work1)

    ENDDO

    if (debug) then
       WRITE(*,*) '********************'
       WRITE(*,*) 'Make_Grdiag_mem done'
       WRITE(*,*) '********************'
    endif

  END SUBROUTINE Make_Grdiag_mem_dns





  !**************************************************************************
  !
  !  Calculate Green Retarded column "n" - writing on memory 
  !
  !**************************************************************************

  SUBROUTINE Make_Grcol_mem_dns(ESH,n,indblk)

    !***********************************************************************
    !Input:
    !ESH: sparse matrices array ESH(nbl,nbl)
    !n: n umber of column to be calculated
    !
    !global variables needed: nbl (number of layers), indblk(nbl+1), 
    !Gr diagonal, subadiagonal and superdiagonal, gsmr(:) for 
    !downgoing and gsml(:) for upgoing 
    !
    !Output:
    !sparse matrices array global variable Gr(:,n) is available in 
    !memory - single blocks are allocated internally, array Gr(nbl,nbl) 
    !must be allocated externally
    !*********************************************************************** 

    IMPLICIT NONE

    !In/Out
    TYPE(z_DNS), DIMENSION(:,:) :: ESH
    INTEGER :: n
    INTEGER, DIMENSION(:), POINTER :: indblk 

    !Work
    INTEGER :: i,nrow,ncol,nbl
    TYPE(z_DNS) :: work1
    REAL(dp) :: max

    nbl = size(ESH,1)

    IF (n.GT.nbl) THEN
       STOP 'Error in Make_Grcol_mem : n is greater than nbl'
    ENDIF

    !***************************************
    !Downgoing (j>=n+2 && n<nbl-1)
    !
    !   G_j,n = -gR_jj T_j,j-1 G_j-1,n
    !
    !***************************************
    ncol=indblk(n+1)-indblk(n)

    IF (n.LT.(nbl-1)) THEN

       DO i=n+2,nbl

          nrow=indblk(i+1)-indblk(i)

          max=MAXVAL(ABS(Gr(i-1,n)%val))

          IF (max.GT.drop) THEN

             CALL prealloc_mult(gsmr(i),ESH(i,i-1),(-1.d0, 0.d0),work1)
             CALL prealloc_mult(work1,Gr(i-1,n),Gr(i,n))
             CALL destroy(work1)

          ELSE

             CALL create(Gr(i,n),nrow,ncol)
             Gr(i,n)%val(:,:)=(0.d0,0.d0)

          ENDIF

       ENDDO

    ENDIF

    !*************************************
    !Downgoing (j<=n-2 && n>2)
    !
    !   G_j,n = -gL_jj T_j,j+1 G_j+1,n
    !
    !*************************************

    IF (n.GT.2) THEN

       DO i=n-2,1,(-1)

          nrow=indblk(i+1)-indblk(i)

          max=MAXVAL(ABS(Gr(i+1,n)%val))

          IF (max.GT.drop) THEN

             CALL prealloc_mult(gsml(i),ESH(i,i+1),(-1.d0, 0.d0),work1)
             CALL prealloc_mult(work1,Gr(i+1,n),Gr(i,n))
             CALL destroy(work1)

          ELSE

             CALL create(Gr(i,n),nrow,ncol)
             Gr(i,n)%val(:,:)=(0.d0,0.d0)

          ENDIF

       ENDDO

    ENDIF

    if (debug) then
       WRITE(*,*) '********************'
       WRITE(*,*) 'Make_Grcol_mem done column',n
       WRITE(*,*) '********************'
    endif

  END SUBROUTINE Make_Grcol_mem_dns




  !****************************************************************************
  !
  !  Calculate Spectral Density - writing on memory
  !
  !****************************************************************************

!  SUBROUTINE Make_Spectral_mem(struct,A)

    !****************************************************************************
    !Input:
    !nrow_tot: total device matrix rows
    !
    !global variable needed: nbl, indblk(nbl+1), Gr(:,:) 
    !
    !Output:
    !A: sparse matrix containing spectral density of device (allocated internally)
    !****************************************************************************

!    IMPLICIT NONE

    !In/Out
!    TYPE(z_CSR) :: A
!    TYPE(Tstruct_info) :: struct

    !Work
!    INTEGER :: i,nrow,nrow_prev,i1,j1,ierr,nrow_tot,nbl
!    TYPE(z_CSR), DIMENSION(:,:), ALLOCATABLE :: Asub
!    INTEGER, DIMENSION(:), POINTER :: indblk

!    nbl = struct%num_PLs
!    indblk => struct%mat_PL_start
!    nrow_tot = struct%total_dim


    !Allocazione dell'array di sparse Asub e di A
!    ALLOCATE(Asub(nbl,nbl),stat=ierr)
!    IF (ierr.NE.0) THEN
!       STOP 'ALLOCATION ERROR: could not allocate Asub(nbl,nbl)'
!    ENDIF

!    CALL create(A,nrow_tot,nrow_tot,0)
!    A%rowpnt(:)=1

    !***
    !A(1,1) 
    !***
!    nrow=indblk(2)-indblk(1)

!    CALL zspectral(Gr(1,1),Gr(1,1),0,Asub(1,1))

!    CALL concat(A,Asub(1,1),1,1)
!    CALL destroy(Asub(1,1))

    !***
    !Diagonal, Subdiagonal and Super diagonal blocks
    !***

!    DO i=2,nbl

!       nrow=indblk(i+1)-indblk(i)

!       i1=indblk(i)
!       j1=indblk(i)

!       CALL zspectral(Gr(i,i),Gr(i,i),0,Asub(i,i))

!       CALL concat(A,Asub(i,i),i1,j1)
!       CALL destroy(Asub(i,i))

!       i1=indblk(i-1)
!       j1=indblk(i)

!       CALL zspectral(Gr(i-1,i),Gr(i,i-1),0,Asub(i-1,i))
!       CALL concat(A,Asub(i-1,i),i1,j1)
!       CALL destroy(Asub(i-1,i))  

!       i1=indblk(i)
!       j1=indblk(i-1)

!       CALL zspectral(Gr(i,i-1),Gr(i-1,i),0,Asub(i,i-1))

       !Blocks concatenation
!       CALL concat(A,Asub(i,i-1),i1,j1)
!       CALL destroy(Asub(i,i-1))

!    ENDDO

!    DEALLOCATE(Asub)

!    if (debug) then    
!       WRITE(*,*) '**********************'
!       WRITE(*,*) 'Make_Spectral_mem done'
!       WRITE(*,*) '**********************'
!    endif

!  END SUBROUTINE Make_Spectral_mem





  !****************************************************************************
  !
  !  Calculate Green Retarded - writing on memory
  !
  !****************************************************************************

!!$  SUBROUTINE Make_GreenR_mem(nrow_tot,A)
!!$
!!$    !****************************************************************************
!!$    !Input:
!!$    !nrow_tot: total device matrix rows
!!$    !
!!$    !global variable needed: nbl, indblk(nbl+1), Gr(:,:) 
!!$    !
!!$    !Output:
!!$    !A: sparse matrix containing Green Retarded of device (allocated internally)
!!$    !****************************************************************************
!!$
!!$    IMPLICIT NONE
!!$
!!$    !In/Out
!!$    INTEGER :: nrow_tot
!!$    TYPE(z_CSR) :: A
!!$
!!$    !TYPE(z_CSR), DIMENSION(nbl,nbl) :: ESH
!!$
!!$    !Work
!!$    INTEGER :: i,nrow,i1,j1,ierr
!!$
!!$    CALL create(A,nrow_tot,nrow_tot,0)
!!$    A%rowpnt(:)=1
!!$
!!$    !write(*,*) 'A created'
!!$
!!$    !***
!!$    !A(1,1)
!!$    !***
!!$    nrow=indblk(2)-indblk(1)
!!$
!!$    CALL concat(A,Gr(1,1),1,1)
!!$
!!$    !***
!!$    !Diagonal, Subdiagonal and Superdiagonal blocks
!!$    !***
!!$    DO i=2,nbl
!!$
!!$       nrow=indblk(i+1)-indblk(i)
!!$
!!$       !       write(*,*) 'nrow=',nrow
!!$
!!$       i1=indblk(i)
!!$       j1=indblk(i)
!!$
!!$       CALL concat(A,Gr(i,i),i1,j1)
!!$
!!$       !       write(*,*) 'Gr',i,i,'concat'
!!$
!!$       i1=indblk(i-1)
!!$       j1=indblk(i)
!!$
!!$       CALL concat(A,Gr(i-1,i),i1,j1)
!!$
!!$       !       write(*,*) 'Gr',i-1,i,'concat'
!!$
!!$       i1=indblk(i)
!!$       j1=indblk(i-1)
!!$
!!$       CALL concat(A,Gr(i,i-1),i1,j1)
!!$
!!$       !       write(*,*) 'Gr',i,i-1,'concat'
!!$
!!$
!!$       !       write(*,*) 'Gr dealloc'
!!$
!!$    ENDDO
!!$
!!$    !if (debug) call writePeakInfo(6)    
!!$    if (debug) then
!!$       WRITE(*,*) '**********************'
!!$       WRITE(*,*) 'Make_GreenR_mem done'
!!$       WRITE(*,*) '**********************'
!!$    endif
!!$
!!$  END SUBROUTINE Make_GreenR_mem


  !****************************************************************************
  !
  !  Calculate Green Retarded - writing on memory (optimized on mask)
  !
  !****************************************************************************

  SUBROUTINE Make_Gr_mem_dns(P,nbl,indblk,A)

    !****************************************************************************
    !Input:
    !P: CSR matrix containing masking pattern
    !
    !global variable needed: nbl, indblk(nbl+1), Gr(:,:) 
    !
    !Output:
    !A: sparse matrix containing Green Retarded of device (allocated internally)
    !****************************************************************************

    IMPLICIT NONE

    !In/Out
    INTEGER :: nbl
    INTEGER, DIMENSION(:), POINTER :: indblk
    TYPE(z_CSR) :: A, P, GrCsr

    !Work
    INTEGER :: i, j, i1, ix, iy, x, y, col, oldx

    !create A with same pattern of P
    CALL create(A,P%nrow,P%ncol,P%nnz)
    A%rowpnt(:)=P%rowpnt(:)
    A%colind(:)=P%colind(:)


    !If only one block is present, concatenation is not needed 
    !and it's implemented in a more trivial way
    IF (nbl.EQ.1) THEN

       call create(GrCsr,Gr(1,1)%nrow,Gr(1,1)%ncol,Gr(1,1)%nrow*Gr(1,1)%ncol)
       call dns2csr(Gr(1,1),GrCsr)

       call mask(GrCsr,P,A)
       call destroy(GrCsr)

    ELSE  

       !Cycle upon all rows
       x = 1    
       DO i = 1, A%nrow
          !Choose which block (row) we're dealing with
          oldx = x

          !Check if row is in same block of previous or in next block. Not needed 
          !(and not allowed not to exceed indblk index boundaries) if we're in the last block
          IF (oldx.EQ.nbl) THEN 
             x = oldx
          ELSE
             DO ix = oldx, oldx+1
                IF ( (i.GE.indblk(ix)).AND.(i.LT.indblk(ix+1)) ) x = ix
             ENDDO
          ENDIF

          !Offset: i1 is the index for separate blocks
          i1 = i - indblk(x) + 1
          !Cycle upon columns 
          DO j = A%rowpnt(i), A%rowpnt(i+1) -1
             !Choose which block column we're dealing with
             y = 0
             IF (x.EQ.1) THEN
                IF ( (A%colind(j).GE.indblk(x)).AND.(A%colind(j).LT.indblk(x + 1)) ) then 
                   y = 1
                ELSEIF ( (A%colind(j).GE.indblk(x + 1)).AND.(A%colind(j).LT.indblk(x + 2)) ) then 
                   y = 2
                ENDIF
             elseif (x.eq.nbl) then
                IF ( (A%colind(j).GE.indblk(x)).AND.(A%colind(j).LT.indblk(x + 1)) ) then 
                   y = nbl
                ELSEIF ( (A%colind(j).GE.indblk(x - 1)).AND.(A%colind(j).LT.indblk(x)) ) then 
                   y = nbl - 1
                ENDIF
             ELSE
                DO iy = x-1, x+1
                   if ( (A%colind(j).GE.indblk(iy)).AND.(A%colind(j).LT.indblk(iy + 1)) ) y = iy
                ENDDO
             ENDIF
             IF (y.eq.0) then
                write(*,*) 'ERROR in Make_GreenR_mem2: probably wrong PL size',x
                write(*,*) 'row',i,A%colind(j)
                write(*,*) 'block indeces:',indblk(1:nbl)
                stop
             ENDIF

             col = A%colind(j) - indblk(y) + 1

             A%nzval(j) = (0.0, 0.0)

             A%nzval(j) = Gr(x,y)%val(i1,col) 

          ENDDO

       ENDDO

    ENDIF

    !if (debug) call writePeakInfo(6)    
    if (debug) then
       WRITE(*,*) '**********************'
       WRITE(*,*) 'Make_GreenR_mem2 done'
       WRITE(*,*) '**********************'
    endif

  END SUBROUTINE Make_Gr_mem_dns





  !****************************************************************************
  !
  ! Calculate G_n contributions for all contacts (except reference)
  ! Writing on memory
  !
  !****************************************************************************

  SUBROUTINE Make_Gn_mem_dns(ESH,SelfEneR,frm,ref,struct,Gn)

    !******************************************************************************
    !Input:  
    !ESH(nbl,nbl): sparse matrices array ES-H
    !SelfEneR(ncont): sparse matrices array with contacts Self Energy 
    !frm(ncont): Fermi diistribution value for each contact 
    !ref:  reference contact 
    !
    !global variables needed: nbl, indblk(nbl+1), ncont, cblk(ncont), Gr(:,:) 
    !diagonal, subdiagonal, overdiagonal and in colums Gr(:,cb) where cb are
    !layers interacting with all contacts but collector
    !
    !Output:
    !Gl: sparse matrix containing G  contacts contribution
    !*******************************************************************************  

    IMPLICIT NONE

    !In/Out
    TYPE(z_DNS), DIMENSION(:,:) :: ESH, Gn
    TYPE(z_DNS), DIMENSION(:) :: SelfEneR
    TYPE(Tstruct_info), intent(in) :: struct
    REAL(dp), DIMENSION(:) :: frm
    INTEGER :: ref 

    !Work
    Type(z_DNS) :: Gam
    TYPE(z_DNS) :: work1,Ga
    INTEGER :: ierr,i,j,cb
    INTEGER :: ncont, nbl
    INTEGER :: oldx, row, col, iy, ix, x, y, ii, jj
    INTEGER, DIMENSION(:), POINTER :: cblk
    COMPLEX(dp) :: frmdiff

    ncont = struct%num_conts    
    nbl = struct%num_PLs
    cblk => struct%cblk

    !*******************************************
    ! Contact Iteration
    !*******************************************
    DO j=1,ncont
   
       IF (j.NE.ref .AND. ABS(frm(j)-frm(ref)).GT.drop) THEN
 
          cb=cblk(j) ! block corresponding to contact j

          CALL zspectral(SelfEneR(j),SelfEneR(j),0,Gam)

          frmdiff = cmplx(frm(j)-frm(ref),0.d0,dp)

          ! Computation of Gl(1,1) 
          CALL prealloc_mult(Gr(1,cb),Gam,work1)    
          CALL zdagger(Gr(1,cb),Ga)
          CALL prealloc_mult(work1,Ga,frmdiff,Gn(1,1))
          CALL destroy(work1, Ga)

          ! Computation of all tridiagonal blocks
          DO i=2,nbl

             CALL prealloc_mult(Gr(i-1,cb),Gam,work1)
             CALL zdagger(Gr(i,cb),Ga)
             CALL prealloc_mult(work1,Ga,frmdiff,Gn(i-1,i))
             CALL destroy(work1)
             
             CALL prealloc_mult(Gr(i,cb),Gam,work1)
             CALL prealloc_mult(work1,Ga,frmdiff,Gn(i,i))

             CALL destroy(work1, Ga)

             CALL prealloc_mult(Gr(i,cb),Gam,work1)
             CALL zdagger(Gr(i-1,cb),Ga)
             CALL prealloc_mult(work1,Ga,frmdiff,Gn(i,i-1))

             CALL destroy(work1, Ga)

          ENDDO

          call destroy(Gam)

      ENDIF 

    ENDDO

  END SUBROUTINE Make_Gn_mem_dns


  !****************************************************************************
  !
  ! Calculate G_n contributions due to el-ph
  ! Writing on memory
  !
  !****************************************************************************
  SUBROUTINE Make_Gn_ph(pnegf,ESH,iter,Gn)

    TYPE(Tnegf), pointer :: pnegf
    TYPE(z_DNS), DIMENSION(:,:) :: ESH, Gn
    INTEGER :: iter

    Type(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_ph_n
    Type(z_DNS) :: Ga, work1, work2
    INTEGER :: n, k, nbl, nrow, ierr

    nbl = pnegf%str%num_PLs
    ALLOCATE(Sigma_ph_n(nbl,nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_n'

    DO n = 1, nbl
      nrow = ESH(n,n)%nrow
      call create(Sigma_ph_n(n,n), nrow, nrow)
    
      Sigma_ph_n(n,n)%val = (0.0_dp, 0.0_dp)
    
      if (iter .gt. 0) then
         call read_blkmat(Sigma_ph_n(n,n),pnegf%scratch_path,'Sigma_ph_n_',n,n,pnegf%iE)
      else 
         call write_blkmat(Sigma_ph_n(n,n),pnegf%scratch_path,'Sigma_ph_n_',n,n,pnegf%iE)
      endif
    END DO
    
    ! Computing diagonal blocks of Gn(n,n)

    DO n = 1, nbl

      DO k = 1, nbl

         CALL zdagger(Gr(n,k),Ga)
         CALL prealloc_mult(Gr(n,k), Sigma_ph_n(k,k), work1)
         CALL prealloc_mult(work1, Ga, work2)
         Gn(n,n)%val = Gn(n,n)%val + work2%val
         call destroy(work1, work2, Ga)

      END DO    
    
    END DO
    
    DO n = 1, nbl
      CALL destroy(Sigma_ph_n(n,n))
    END DO

    DEALLOCATE(Sigma_ph_n)

  END SUBROUTINE Make_Gn_ph
  
  !****************************************************************************
  !
  ! Calculate G_p=iG> contributions due to el-ph
  ! Writing on memory
  !
  !****************************************************************************
  SUBROUTINE Make_Gp_ph(pnegf,ESH,iter,Gp)

    TYPE(Tnegf), pointer :: pnegf
    TYPE(z_DNS), DIMENSION(:,:) :: ESH, Gp
    INTEGER :: iter

    Type(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_ph_p
    Type(z_DNS) :: Ga, work1, work2 
    INTEGER :: n, k, nbl, nrow, ierr

    nbl = pnegf%str%num_PLs
    ALLOCATE(Sigma_ph_p(nbl,nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_n'


    DO n = 1, nbl

      nrow = ESH(n,n)%nrow

      call create(Sigma_ph_p(n,n), nrow, nrow)

      Sigma_ph_p(n,n)%val = (0.0_dp, 0.0_dp)
      if (iter .gt. 0) then
         call read_blkmat(Sigma_ph_p(n,n),pnegf%scratch_path,'Sigma_ph_p_',n,n,pnegf%iE)
      else 
         call write_blkmat(Sigma_ph_p(n,n),pnegf%scratch_path,'Sigma_ph_p_',n,n,pnegf%iE)
      endif
    
    END DO

    DO n = 1, nbl

      DO k = 1, nbl

         CALL zdagger(Gr(n,k),Ga)
         CALL prealloc_mult(Gr(n,k), Sigma_ph_p(k,k), work1)
         CALL prealloc_mult(work1, Ga, work2)
         Gp(n,n)%val = Gp(n,n)%val + work2%val
         call destroy(work1, work2, Ga)

      END DO    
    
    END DO
    
    DO n = 1, nbl
      CALL destroy(Sigma_ph_p(n,n))
    END DO

    DEALLOCATE(Sigma_ph_p)

  END SUBROUTINE Make_Gp_ph

  !Concatenation for every contact in G_n. Performs a sum on elements, not a replacement
  !Similar to Make_GreenR_mem2, except for sum of elements
  !Note: to backup old version zconcat calls (and Glsub deallocations) must be
  !      uncommented and all this part removed
  !If only one block is present, concatenation is not needed and it's implemented in a
  !more trivial way
  SUBROUTINE blk2csr(Gn,struct,P,Gl)

    TYPE(z_DNS), DIMENSION(:,:) :: Gn
    TYPE(Tstruct_info), intent(in) :: struct
    TYPE(z_CSR) :: Gl
    TYPE(z_CSR) :: P, Gl_sp

    INTEGER, DIMENSION(:), POINTER :: indblk, cblk
    INTEGER :: nbl, oldx, row, col, iy, ix, x, y, ii, jj

    nbl = struct%num_PLs
    cblk => struct%cblk
    indblk => struct%mat_PL_start

    !create Gl with same pattern of P

    CALL create(Gl,P%nrow,P%ncol,P%nnz)
    Gl%rowpnt = P%rowpnt
    Gl%colind = P%colind
    Gl%nzval = (0.d0, 0.d0)

    IF (nbl.EQ.1) THEN
        IF (ALLOCATED(Gn(1,1)%val)) THEN
            ii = nzdrop(Gn(1,1),drop) 
        ELSE   
            ii = 0 
        ENDIF    

        IF (ii.gt.0) then
           call create(Gl_sp, Gn(1,1)%nrow, Gn(1,1)%ncol, ii)
           call dns2csr(Gn(1,1),Gl_sp)
           call zmask_realloc(Gl_sp,P)
           call concat(Gl,Gl_sp,1,1)
           call destroy(Gl_sp)
        ENDIF  

        call destroy(Gn(1,1))

    ELSE

        !Cycle upon all rows
        x = 1
        DO ii = 1, Gl%nrow
                !Choose which block (row) we're dealing with
                oldx = x
                IF (oldx.EQ.nbl) THEN 
                   x = oldx
                ELSE
                   DO ix = oldx, oldx+1
                      IF ( (ii.GE.indblk(ix)).AND.(ii.LT.indblk(ix+1)) ) x = ix
                   ENDDO
                ENDIF

                IF (x.EQ.0) THEN
                   WRITE(*,*) 'Error in Make_GreenR_mem2: &
                               &could not find row block of index ',ii,ix
                   STOP
                ENDIF

                !Offset: row is the index for separate blocks
                row = ii - indblk(x) + 1

                !Cycle upon columns of Gl (which has been ALREADY MASKED by ESH)
                DO jj = Gl%rowpnt(ii), Gl%rowpnt(ii+1) -1
                   !Choose which block column we're dealing with
                   y = 0
                   if (x.eq.1) then
                      IF ( (Gl%colind(jj).GE.indblk(x)).AND.(Gl%colind(jj).LT.indblk(x + 1)) ) then 
                         y = 1
                      ELSEIF ( (Gl%colind(jj).GE.indblk(x + 1)).AND.(Gl%colind(jj).LT.indblk(x + 2)) ) then 
                         y = 2
                      ENDIF
                   elseif (x.eq.nbl) then
                      IF ( (Gl%colind(jj).GE.indblk(x)).AND.(Gl%colind(jj).LT.indblk(x + 1)) ) then 
                         y = nbl
                      ELSEIF ( (Gl%colind(jj).GE.indblk(x - 1)).AND.(Gl%colind(jj).LT.indblk(x)) ) then 
                         y = nbl - 1
                      ENDIF
                   else
                      DO iy = x-1, x+1
                         if ( (Gl%colind(jj).GE.indblk(iy)).AND.(Gl%colind(jj).LT.indblk(iy + 1)) ) y = iy
                      ENDDO
                   ENDIF
                   IF (y.EQ.0) THEN
                      WRITE(*,*) 'ERROR in Make_Gn_mem: probably wrong PL size', x
                      write(*,*) 'row',ii,Gl%colind(jj)
                      write(*,*) 'block indeces:',indblk(1:nbl)
                      STOP
                   ENDIF
                   
                   col = Gl%colind(jj) - indblk(y) + 1

                   IF (allocated(Gn(x,y)%val)) THEN
                       Gl%nzval(jj) = Gl%nzval(jj) + Gn(x,y)%val(row,col)
                   ENDIF

                ENDDO

                IF(oldx.NE.x) THEN

                   IF (x.eq.2)  THEN 
                      CALL destroy(Gn(x-1,x-1))
                      CALL destroy(Gn(x-1,x))

                   ELSEIF  (x.GT.2)  THEN 
                      CALL destroy(Gn(x-1,x-1))
                      CALL destroy(Gn(x-1,x-2))
                      CALL destroy(Gn(x-1,x))
                   endif

                ENDIF

        ENDDO

        CALL destroy(Gn(nbl, nbl))
        CALL destroy(Gn(nbl, nbl - 1))

    ENDIF

  END SUBROUTINE blk2csr

  ! ******************************************************************************
  ! Computes Sigma_ph_n and save it file
  ! ******************************************************************************
  SUBROUTINE Sigma_ph_n(pnegf,Epnt)

    TYPE(Tnegf),POINTER :: pnegf
    REAL(dp), DIMENSION(:),allocatable :: Epnt

    !Local variables
    TYPE(z_DNS) :: work1,work2
    TYPE(z_DNS) :: Mq_mat
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_n, G_n_interP, G_n_interN

    REAL(dp), DIMENSION(:), POINTER :: Wq
    REAL(dp), DIMENSION(:), POINTER :: Nq
    REAL(dp), DIMENSION(:), POINTER :: Mq
    LOGICAL, DIMENSION(:), POINTER :: selmodes

    INTEGER :: i, m, iE, i1,i2
    INTEGER :: nummodes, numselmodes, nbl
    INTEGER, DIMENSION(:), pointer :: indblk
    REAL(dp) :: E1, E2, En

    nbl = pnegf%str%num_PLs
    indblk => pnegf%str%mat_PL_start

    selmodes => pnegf%elph%selmodes
    Wq => pnegf%elph%Wq
    Mq => pnegf%elph%Mq
    Nq => pnegf%elph%Nq
    numselmodes = pnegf%elph%numselmodes
    nummodes = pnegf%elph%nummodes

    iE = pnegf%iE

    allocate(G_n_interP(nbl,nbl))
    allocate(G_n_interN(nbl,nbl))
    allocate(Sigma_n(nbl,nbl))

    do i = 1, nbl
      m = indblk(i+1)-indblk(i)
      call create(Sigma_n(i,i), m, m)
      Sigma_n(i,i)%val = (0.d0, 0.d0)
      call create(G_n_interP(i,i), m, m)
      call create(G_n_interN(i,i), m, m)
    enddo
    
    do m = 1 , nummodes

      if (.not.selmodes(m)) cycle
      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(pnegf,Wq(m),Epnt,i1,i2,E1,E2)
      
      En = real(pnegf%Epnt)+Wq(m)
      
      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
                             'G_n_', G_n_interP)

      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(pnegf,-Wq(m),Epnt,i1,i2,E1,E2)
      
      En = real(pnegf%Epnt)-Wq(m)
      
      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
                              'G_n_', G_n_interN)


      do i = 1, nbl

        i1 =  G_n_interN(i,i)%nrow

        call create(work1, i1, i1)       

        work1%val = (Nq(m)+1.0_dp)*G_n_interP(i,i)%val + Nq(m)* G_n_interN(i,i)%val

        call create_id(Mq_mat,i1,Mq(m))

        call prealloc_mult( Mq_mat, work1, work2)

        call destroy(work1)

        call prealloc_mult( work2, Mq_mat, work1)

        call destroy(work2)

        Sigma_n(i,i)%val = Sigma_n(i,i)%val + work1%val 

        call destroy(work1, Mq_mat)

     end do

   end do

   !Throw away all non diagonal parts
   !do i = 1, nbl
   !   i1 = Sigma_n(i,i)%nrow
   !   do m = 1, i1 
   !      Sigma_n(i,i)%val(1:m-1,m) = (0.0_dp,0.0_dp) 
   !      Sigma_n(i,i)%val(m+1:i1,m) = (0.0_dp,0.0_dp) 
   !   enddo
   !enddo
   !

   do i = 1, nbl
     call write_blkmat(Sigma_n(i,i),pnegf%scratch_path,'Sigma_ph_n_',i,i,iE)
     call destroy(Sigma_n(i,i))
     call destroy(G_n_interN(i,i))
     call destroy(G_n_interP(i,i))
   enddo

   deallocate(Sigma_n,G_n_interN, G_n_interP)

 END SUBROUTINE Sigma_ph_n


  ! ******************************************************************************
  ! Computes Sigma_ph> and save it file
  ! ******************************************************************************
  SUBROUTINE Sigma_ph_p(pnegf,Epnt)

    TYPE(Tnegf),POINTER :: pnegf
    REAL(dp), DIMENSION(:),allocatable :: Epnt

    !Local variables
    TYPE(z_DNS) :: work1,work2
    TYPE(z_DNS) :: Mq_mat
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_p, G_p_interP, G_p_interN

    REAL(dp), DIMENSION(:), POINTER :: Wq
    REAL(dp), DIMENSION(:), POINTER :: Nq
    REAL(dp), DIMENSION(:), POINTER :: Mq
    LOGICAL, DIMENSION(:), POINTER :: selmodes

    INTEGER :: i, m, iE, i1,i2
    INTEGER :: nummodes, numselmodes, nbl
    INTEGER, DIMENSION(:), pointer :: indblk
    REAL(dp) :: E1, E2, En

    nbl = pnegf%str%num_PLs
    indblk => pnegf%str%mat_PL_start

    selmodes => pnegf%elph%selmodes
    Wq => pnegf%elph%Wq
    Mq => pnegf%elph%Mq
    Nq => pnegf%elph%Nq
    numselmodes = pnegf%elph%numselmodes
    nummodes = pnegf%elph%nummodes

    iE = pnegf%iE

    allocate(G_p_interP(nbl,nbl))
    allocate(G_p_interN(nbl,nbl))
    allocate(Sigma_p(nbl,nbl))

    do i = 1, nbl
      m = indblk(i+1)-indblk(i)
      call create(Sigma_p(i,i), m, m)
      Sigma_p(i,i)%val = (0.d0, 0.d0)
      call create(G_p_interP(i,i), m, m)
      call create(G_p_interN(i,i), m, m)
    enddo
    
    do m = 1 , nummodes

      if (.not.selmodes(m)) cycle
      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(pnegf,Wq(m),Epnt,i1,i2,E1,E2)
    
      En = real(pnegf%Epnt)+Wq(m)
      
      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
                             'G_p_', G_p_interP)
      
      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(pnegf,-Wq(m),Epnt,i1,i2,E1,E2)
      En = real(pnegf%Epnt)-Wq(m)
      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
                              'G_p_', G_p_interN)

      do i = 1, nbl

        !print*
        !print*,'(sigma_r) G_p+',maxval(abs(G_p_interP(i,i)%val))
        !print*,'(sigma_r) G_p-',maxval(abs(G_p_interN(i,i)%val))

        i1 =  G_p_interN(i,i)%nrow

        call create(work1, i1, i1)       

        work1%val = (Nq(m)+1.0_dp)*G_p_interN(i,i)%val + Nq(m)* G_p_interP(i,i)%val

        call create_id(Mq_mat,i1,Mq(m))

        call prealloc_mult( Mq_mat, work1, work2)

        call destroy(work1)

        call prealloc_mult ( work2, Mq_mat, work1)

        call destroy(work2)

        Sigma_p(i,i)%val = Sigma_p(i,i)%val + work1%val 

        call destroy(work1, Mq_mat)

     end do

   end do
   
   ! throw away all non diagonal parts
   !do i = 1, nbl
   !   i1 = Sigma_p(i,i)%nrow
   !   do m = 1, i1 
   !      Sigma_p(i,i)%val(1:m-1,m) = (0.0_dp,0.0_dp) 
   !      Sigma_p(i,i)%val(m+1:i1,m) = (0.0_dp,0.0_dp) 
   !   enddo
   !enddo
   !

   do i = 1, nbl
     call write_blkmat(Sigma_p(i,i),pnegf%scratch_path,'Sigma_ph_p_',i,i,iE)
     call destroy(Sigma_p(i,i))
     call destroy(G_p_interN(i,i))
     call destroy(G_p_interP(i,i))
   enddo

  deallocate(Sigma_p,G_p_interN, G_p_interP)

END SUBROUTINE Sigma_ph_p

!----------------------------------------------------------------------------------
SUBROUTINE Sigma_ph_r(pnegf,Epnt)

  TYPE(Tnegf),POINTER :: pnegf
  REAL(dp), DIMENSION(:),allocatable :: Epnt

  !Local variables
  TYPE(z_DNS) :: work1,work2,work3
  TYPE(z_DNS) :: Mq_mat
  TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_r, G_r_interP, G_r_interN
  TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: G_n_interP, G_n_interN, G_r

  REAL(dp), DIMENSION(:), POINTER :: Wq
  REAL(dp), DIMENSION(:), POINTER :: Nq
  REAL(dp), DIMENSION(:), POINTER :: Mq
  LOGICAL, DIMENSION(:), POINTER :: selmodes
  INTEGER :: i, m, iE, i1, i2
  INTEGER :: nummodes, numselmodes, nbl
  INTEGER, DIMENSION(:), pointer :: indblk
  REAL(dp) :: E1, E2, En

  nbl = pnegf%str%num_PLs
  indblk => pnegf%str%mat_PL_start


  selmodes => pnegf%elph%selmodes
  Mq => pnegf%elph%Mq
  Wq => pnegf%elph%Wq
  Nq => pnegf%elph%Nq
  nummodes = pnegf%elph%nummodes
  numselmodes = pnegf%elph%numselmodes

  iE=pnegf%iE

  allocate(G_r_interP(nbl,nbl))
  allocate(G_r_interN(nbl,nbl))
  allocate(G_n_interP(nbl,nbl))
  allocate(G_n_interN(nbl,nbl))
  allocate(Sigma_r(nbl,nbl))
  !allocate(G_r(nbl,nbl))

    
      
  do i = 1, nbl
      m = indblk(i+1)-indblk(i)
      call create(Sigma_r(i,i), m, m)
      Sigma_r(i,i)%val = (0.d0, 0.d0)
      call create(G_r_interP(i,i), m, m)
      call create(G_r_interN(i,i), m, m)
      call create(G_n_interP(i,i), m, m)
      call create(G_n_interN(i,i), m, m)
      !call create(G_r(i,i), m, m)
      !call read_blkmat(G_r(i,i),pnegf%scratch_path,'Gr_',i,i,iE)
  enddo

  do m = 1 , nummodes

      if (.not.selmodes(m)) cycle

      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(pnegf,Wq(m),Epnt,i1,i2,E1,E2)
      En = real(pnegf%Epnt)+Wq(m)
      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
                             'G_n_', G_n_interP)
      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
                             'G_r_', G_r_interP)

      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(pnegf,-Wq(m),Epnt,i1,i2,E1,E2)
      En = real(pnegf%Epnt)-Wq(m)
      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
                             'G_n_', G_n_interN)
      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
                             'G_r_', G_r_interN)

      do i = 1, nbl
         !print*
         !print*,'(sigma_r) G_r+',maxval(abs(G_r_interP(i,i)%val))
         !print*,'(sigma_r) G_r-',maxval(abs(G_r_interN(i,i)%val))
         
         !print*,'(sigma_r) G_n+',maxval(abs(G_n_interP(i,i)%val))
         !print*,'(sigma_r) G_n-',maxval(abs(G_n_interN(i,i)%val))
                         
         i1 = Sigma_r(i,i)%nrow

         call create(work1,i1,i1)
 
         work1%val = (0.d0,0.d0)

         if (pnegf%elph%selfene_gr) then
             !print*,'SelfEneR_Gr'    
             ! Via I: should be exact for E >> Ef_max
             work1%val = (Nq(m)+1.0_dp)*G_r_interN(i,i)%val + Nq(m)*G_r_interP(i,i)%val
             ! Via II: should be exact for E << Ef_min
             !work1%val = (Nq(m)+1.0_dp)*G_r_interP(i,i)%val + Nq(m)*G_r_interN(i,i)%val
             ! Via III: should work as a compromise
             !work1%val = Nq(m)*(G_r_interP(i,i)%val + G_r_interN(i,i)%val) + G_r(i,i)%val
         endif   
            
         if (pnegf%elph%selfene_gless) then
             !print*,'SelfEneR_G<'    
             ! should be 1/2 [G<- - G<+] == i/2 [Gn- - Gn+]    
             work1%val = work1%val + &
                    (0.0_dp,0.5_dp) * (G_n_interN(i,i)%val - G_n_interP(i,i)%val)
         endif

         if (pnegf%elph%selfene_gless .or. pnegf%elph%selfene_gr) then

             call create_id(Mq_mat,i1,Mq(m))
                      
             call prealloc_mult(Mq_mat, work1, work2)

             call destroy(work1)

             call prealloc_mult(work2, Mq_mat, work1)

             call destroy(work2)

             Sigma_r(i,i)%val = Sigma_r(i,i)%val + work1%val 

             call destroy(work1, Mq_mat)

          else

             Sigma_r(i,i)%val = (0.d0,0.d0) 
             call destroy(work1)

          endif   

      enddo

  enddo
   
  ! throw away all non diagonal parts
  !do i = 1, nbl
  !    i1 = Sigma_r(i,i)%nrow
  !    do m = 1, i1 
  !       Sigma_r(i,i)%val(1:m-1,m) = (0.0_dp,0.0_dp) 
  !       Sigma_r(i,i)%val(m+1:i1,m) = (0.0_dp,0.0_dp) 
  !    enddo
  ! enddo
  !

  do i = 1, nbl
   !print*          
   !print*,'(sigma_r) Sigma_ph_r',maxval(abs(Sigma_r(i,i)%val)), iE
    call write_blkmat(Sigma_r(i,i),pnegf%scratch_path,'Sigma_ph_r_',i,i,iE)
    call destroy(Sigma_r(i,i))
    call destroy(G_r_interP(i,i))
    call destroy(G_r_interN(i,i))
    call destroy(G_n_interP(i,i))
    call destroy(G_n_interN(i,i))
    !call destroy(G_r(i,i))
  enddo

  deallocate(Sigma_r,G_r_interP,G_r_interN,G_n_interP,G_n_interN)
  !deallocate(G_r)

END SUBROUTINE Sigma_ph_r

SUBROUTINE Sigma_ph_r_z(pnegf,z)

  TYPE(Tnegf),POINTER :: pnegf
  COMPLEX(dp) :: z
  
  !Local variables
  TYPE(z_DNS) :: work1,work2,work3
  TYPE(z_DNS) :: Mq_mat
  TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_r, G_r
  REAL(dp), DIMENSION(:), POINTER :: Wq
  REAL(dp), DIMENSION(:), POINTER :: Nq
  REAL(dp), DIMENSION(:), POINTER :: Mq
  LOGICAL, DIMENSION(:), POINTER :: selmodes
  INTEGER :: i, m, iE,i1
  INTEGER :: nummodes, numselmodes, nbl
  INTEGER, DIMENSION(:), pointer :: indblk

  nbl = pnegf%str%num_PLs
  indblk => pnegf%str%mat_PL_start


  selmodes => pnegf%elph%selmodes
  Mq => pnegf%elph%Mq
  Wq => pnegf%elph%Wq
  Nq => pnegf%elph%Nq
  nummodes = pnegf%elph%nummodes
  numselmodes = pnegf%elph%numselmodes

  iE=pnegf%iE

  allocate(Sigma_r(nbl,nbl))
  allocate(G_r(nbl,nbl))

  do i = 1, nbl
      m = indblk(i+1)-indblk(i)
      call create(Sigma_r(i,i), m, m)
      Sigma_r(i,i)%val = (0.d0, 0.d0)
      call create(G_r(i,i), m, m)
  enddo 

  do m = 1 , nummodes

    if (.not.selmodes(m)) cycle
    
    do i = 1, nbl
    
      call read_blkmat(G_r(i,i),pnegf%scratch_path,'Gr_',i,i,iE)

      i1 = Sigma_r(i,i)%nrow

      call create(work1,i1,i1)
 
      work1%val = (0.d0,0.d0)

      if (pnegf%elph%selfene_gr) then
          !print*,'SelfEneR_Gr'    
          work1%val = (2*Nq(m)+1.0_dp)*G_r(i,i)%val

          call create_id(Mq_mat,i1,Mq(m))
                   
          call prealloc_mult(Mq_mat, work1, work2)

          call destroy(work1)

          call prealloc_mult(work2, Mq_mat, work1)

          call destroy(work2)

          Sigma_r(i,i)%val = Sigma_r(i,i)%val +  work1%val 

          call destroy(work1, Mq_mat)
       
       else
          Sigma_r(i,i)%val = (0.d0,0.d0) 
       endif
    
     enddo  
     
   enddo  

   do i = 1, nbl
    !print*          
    !print*,'(sigma_r) Sigma_ph_r',maxval(abs(Sigma_r(i,i)%val)), pnegf%iE
    call write_blkmat(Sigma_r(i,i),pnegf%scratch_path,'Sigma_ph_r_',i,i,iE)
    call destroy(Sigma_r(i,i))
    call destroy(G_r(i,i))
  enddo

END SUBROUTINE Sigma_ph_r_z

! ----------------------------------------------------------
SUBROUTINE check_sigma_ph_r(pnegf)
  TYPE(Tnegf),POINTER :: pnegf
  REAL(dp), DIMENSION(:),allocatable :: Epnt
  integer :: ioffset
  
  TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_r, Sigma_p, Sigma_n
  TYPE(z_DNS) :: Gam, T 
  integer :: nbl, n, sizebl, i_start, i_stop, psize, maxpos(2)
  real(dp) :: Wmax, maxdev, tmp, maxG
  INTEGER, DIMENSION(:), pointer :: indblk

  nbl = pnegf%str%num_PLs
  indblk => pnegf%str%mat_PL_start

  allocate(Sigma_n(nbl,nbl))
  allocate(Sigma_p(nbl,nbl))
  allocate(Sigma_r(nbl,nbl))

  maxdev = 0.0_dp
  maxG = 0.0_dp
  psize = 0

  do n = 1, nbl
    sizebl = indblk(n+1)-indblk(n)
    call create(Sigma_r(n,n),sizebl,sizebl)
    call create(Sigma_n(n,n),sizebl,sizebl)
    call create(Sigma_p(n,n),sizebl,sizebl)
    call read_blkmat(Sigma_r(n,n), pnegf%scratch_path, 'Sigma_ph_r_',n,n, pnegf%iE)
    call read_blkmat(Sigma_n(n,n),pnegf%scratch_path,'Sigma_ph_n_',n,n,pnegf%iE)
    call read_blkmat(Sigma_p(n,n),pnegf%scratch_path,'Sigma_ph_p_',n,n,pnegf%iE)

    call zspectral(Sigma_r(n,n),Sigma_r(n,n),0,Gam)

    call create(T,sizebl,sizebl)

    !print*,'CHECK Sigma_ph_n+Sigma_ph_p',maxval(abs(Sigma_n(n,n)%val+Sigma_p(n,n)%val))
  
    T%val = Sigma_n(n,n)%val + Sigma_p(n,n)%val - Gam%val 
    
    tmp = maxval(abs(Gam%val))
    if (tmp .gt. maxG) maxG=tmp

    tmp = maxval(abs(T%val))/maxval(abs(Gam%val)) 
    if (tmp .gt. maxdev) then
        maxdev = tmp
        maxpos = maxloc(abs(T%val)) + psize
    endif     

    !print*
  
    psize = psize + sizebl

    call destroy(Sigma_r(n,n))
    call destroy(Sigma_n(n,n))
    call destroy(Sigma_p(n,n))
    call destroy(Gam)
    call destroy(T)

  enddo

  print*,'CHECK Sigma_ph_r',pnegf%iE, maxG, maxdev

  deallocate(Sigma_n, Sigma_p, Sigma_r)


END SUBROUTINE check_sigma_ph_r


! ----------------------------------------------------------
! Search points for interpolations. 
! Wq has a sign (+/- Wq)
!
SUBROUTINE search_points(pnegf, Wq, Epnt, i1, i2, E1, E2)
 use energy_mesh, only : elem 
 type(TNegf), pointer :: pnegf
 real(dp) :: Wq                
 real(dp), dimension(:), allocatable :: Epnt
 integer, intent(out) :: i1,i2 
 real(dp), intent(out) :: E1, E2

 integer :: iE, iel, istart, iend, ip
 Type(elem), pointer :: pel
 real(dp) :: En

 if (Wq.eq.0) then
    i1 = pnegf%iE
    i2 = pnegf%iE
    E1 = pnegf%Epnt
    E2 = pnegf%Epnt
    return
 endif

 if (allocated(Epnt)) then
   ! Remove offset such that search can work on Epnt(1..N)      
   iE = pnegf%iE - pnegf%Np_n(1) - pnegf%Np_n(2) - pnegf%n_poles
   En = real(pnegf%Epnt) + Wq !Wq carry the right sign      
   !print*
   !print*,'iE', iE, real(pnegf%Epnt) + Wq 

   if (sign(1.0_dp,Wq) .gt. 0) then
     i2 = iE + 1
     if (i2.gt.size(Epnt)) then
        i2 = size(Epnt)
     else    
        do while (Epnt(i2) .lt. En)
           i2 = i2 + 1
        end do
     endif
     i1 = i2 - 1
   else
     i1 = iE - 1
     if (i1.lt.1) then
        i1 = 1 
     else
        do while (Epnt(i1) .gt. En)
           i1 = i1 - 1
        end do
     endif   
     i2 = i1 + 1
   endif
   E1 = Epnt(i1)
   E2 = Epnt(i2)
   ! add back the offset to the point
   i1 = i1 + pnegf%Np_n(1) + pnegf%Np_n(2) + pnegf%n_poles
   i2 = i2 + pnegf%Np_n(1) + pnegf%Np_n(2) + pnegf%n_poles
 
 else  
 
   !if (.not.allocated(pnegf%emesh%pactive)) STOP 'emesh not initialized'

   En = real(pnegf%Epnt) + Wq       
         
   if (sign(1.0_dp,Wq) .gt. 0) then
     istart = pnegf%emesh%iactive
     iend = pnegf%emesh%maxind

     elloop1: do iel = istart, iend
        pel => pnegf%emesh%pactive(iel)%pelem
        do ip = 1, 3
           if (pel%pnt(ip) .gt. En) then
              exit elloop1 
           endif     
        end do
     end do elloop1   
     i1 = pel%map(ip-1)
     i2 = pel%map(ip)
     E1 = pel%pnt(ip-1)
     E2 = pel%pnt(ip)
   else
     istart = pnegf%emesh%iactive
     iend = 1

     elloop2: do iel = istart, iend, -1
        pel => pnegf%emesh%pactive(iel)%pelem
        do ip = 3, 1, -1
           if (pel%pnt(ip) .lt. En) then
              exit elloop2 
           endif     
        end do
     end do elloop2  
     i1 = pel%map(ip)
     i2 = pel%map(ip+1)
     E1 = pel%pnt(ip)
     E2 = pel%pnt(ip+1)
   end if
 end if
 !print*
 !print*,E1,En,E2
 !print*,'interpolate between:',i1,i2

END SUBROUTINE search_points 


SUBROUTINE interpolation(i1,i2, E1, E2, E, path, name, G_interp)
     INTEGER, intent(in) :: i1, i2
     REAL(dp) :: E1, E2, E
     TYPE(z_DNS), DIMENSION(:,:) :: G_interp
     CHARACTER(*) :: path
     CHARACTER(*) :: name

     !local variables
     TYPE(z_DNS) :: work1,work2
     INTEGER :: i

     do i = 1, size(G_interp,1)
        
        call create(work1,G_interp(i,i)%nrow,G_interp(i,i)%ncol)
        work1%val = (0.d0,0.d0)
        call read_blkmat(work1, path, name, i, i, i1)

        if (E1.ne.E2 .and. i1.ne.i2) then
         
          call create(work2,G_interp(i,i)%nrow,G_interp(i,i)%ncol)
          work2%val = (0.d0,0.d0)
          call read_blkmat(work2, path, name, i, i, i2)

          G_interp(i,i)%val = ((E-E1)*work2%val + (E2-E)*work1%val)/(E2-E1)

          call destroy(work2)
      
        else

          G_interp(i,i)%val = work1%val
        
        endif
        
        call destroy(work1)

     end do


END SUBROUTINE interpolation


!*********************************************************************
!ADD Hilbert-transform part to Sigma_ph_r  
!
!********************************************************************
SUBROUTINE complete_sigma_ph_r(pnegf, Epnt, ioffset)

  TYPE(Tnegf),POINTER :: pnegf
  REAL(dp), DIMENSION(:) :: Epnt
  INTEGER :: ioffset

  ! Locals
  INTEGER :: i,j,k,n,m, iE, bl, sizebl,i_start,i_stop, ierr, nummodes
  REAL(dp), DIMENSION(:), POINTER :: Wq
  REAL(dp), DIMENSION(:), POINTER :: Mq
  REAL(dp) :: Wmax, dE, tmp
  
  CHARACTER(64) :: filename
  character(4) :: ofblki, ofblkj
  character(10) :: ofpnt
  complex(dp), DIMENSION(:), allocatable :: temp1 
  type(z_DNS), DIMENSION(:), allocatable :: Gn_E, Sigma_r_E
  type(z_DNS) :: work1, work2, work3, Mq_mat 
  complex(dp) :: temp2 
  LOGICAL, DIMENSION(:), POINTER :: selmodes

  Mq => pnegf%elph%Mq
  Wq => pnegf%elph%Wq
  selmodes => pnegf%elph%selmodes

  n = size(Epnt)
  nummodes = size(Wq)

  Wmax = maxval(Wq)
  dE = Epnt(2)-Epnt(1)
  ! ENERGY-INTERVAL FOR SELF-ENERGIES:
  i_start = aint(Wmax/dE) + 1
  i_stop =  (n - i_start) 
  i_start = i_start + 1
  
  ! PROBLEM: 
  ! The Hilb-Transf should be resticted on the sub interval
  ! Emin+m*Wmax ... Emax-m*Wmax
  ! However the number of points should be always 2^p
  ! This condition is difficult to fulfill.
  ! Possible solution: interpolation between grids.

  ! CREATE THE ARRAY OF MATRICES (one for each E-point)
  allocate(Gn_E(n),stat=ierr)
  allocate(Sigma_r_E(n),stat=ierr)
  if(ierr.ne.0) STOP 'ERROR in allocation of Gn_E'
  call log_allocate(temp1,n)

  print*
  print*,'HILBERT TRANSFORM memory:',sizebl*sizebl*n*16
  print*,'HILBERT TRANSFORM interval:',i_start,i_stop
 
  ! LOOP ON blocks
  do bl = 1, pnegf%str%num_PLs

    sizebl =  pnegf%str%mat_PL_start(bl+1) - pnegf%str%mat_PL_start(bl)
    if (bl.le.9999) write(ofblki,'(i4.4)') bl 
    if (bl.gt.9999) stop 'ERROR: too many blks (> 9999)'
 
    ! LOAD ALL G< FROM FILES and store in Gn_E(iE)
    do iE = 1, n
      call create(Gn_E(iE),sizebl,sizebl)
      call read_blkmat(Gn_E(iE), pnegf%scratch_path, 'G_n_', bl, bl,iE+ioffset)
    enddo
    
    ! LOAD ALL Sigma_r in the right interval
    do iE = i_start, i_stop
      call create(Sigma_r_E(iE),sizebl,sizebl)
      call read_blkmat(Sigma_r_E(iE), pnegf%scratch_path, 'Sigma_ph_r_', bl, bl, iE+ioffset)
    enddo

    do m = 1, nummodes 
      if (.not.selmodes(m)) cycle

      ! COMPUTE   Mq G< Mq   (assume diagonal now) 
      do iE = 1, n
        !call create_id(Mq_mat,sizebl,Mq(m))
        !call prealloc_mult(Mq_mat, Gn_E(iE), work1)
        !call prealloc_mult(work1, Mq_mat, work2)
        !call destroy(work1, Mq_mat)
        !Gn_E(iE)%val = work2%val
        !call destroy(work2)
        tmp = Mq(m)*Mq(m)
        Gn_E(iE)%val = tmp*Gn_E(iE)%val
      enddo

      ! PERFORM Hilbert in Energy for each (i,j)     
      do j = 1, sizebl  
        do i = 1, sizebl

!          k = (j-1)*sizebl+i
!          do iE = 1, n
!             if (iE.le.99999) write(ofpnt,'(i5.5)') iE+ioffset
!             filename = 'G_n_'//trim(ofblki)//'_'//trim(ofblki)//'_'//trim(ofpnt)//'.dat'
!             open(100,file=trim(pnegf%scratch_path)//trim(filename), access='DIRECT', recl=4)
!             READ(100, rec = k)  temp1(iE)
!             close(100)
!          enddo

           ! SETUP a vector out of all G<_ij(E)
           ! Here we could perform an efficient interpolation on a regular grid 2^p
           do iE = 1, n
             temp1(iE) = Gn_E(iE)%val(i,j)
           enddo
 
           Wq(m) = Wq(m)*2.0_dp*pi/((n-1)*dE)
 
           call Hilbert_shift(temp1, Wq(m))
         
           Wq(m) = Wq(m)*((n-1)*dE)/(2.0_dp*pi)
 
           ! UPDATE the self-energies with the Hilbert part. 
           do iE = i_start, i_stop
             ! Should be  -i/2 H[ G<+ - G<-  ] = 1/2 H[ Gn+ - Gn- ]
             Sigma_r_E(iE)%val(i,j) =  Sigma_r_E(iE)%val(i,j) + (0.5_dp, 0.0)* temp1(iE) 
           enddo


!          do iE = i_start+1, i_stop
!             if (iE.le.99999) write(ofpnt,'(i5.5)') iE+ioffset
!             filename = 'Sigma_ph_r_'//trim(ofblki)//'_'//trim(ofblki)//'_'//trim(ofpnt)//'.dat'
!             open(100,file=trim(pnegf%scratch_path)//trim(filename), access='DIRECT', recl=4)
!             READ (100, rec = k) temp2 
!             temp2 = temp2 - (0.0_dp, 0.5_dp)* temp1(iE)   
!             WRITE (100, rec = k) temp2
!             close(100)
!          enddo

        enddo ! Loop on block size   
      enddo ! Loop on block size 

    enddo !Loop on modes 

    do iE = i_start, i_stop
      call write_blkmat(Sigma_r_E(iE), pnegf%scratch_path, 'Sigma_ph_r_', bl, bl, iE+ioffset)
    enddo

    do iE = 1, n
      call destroy(Gn_E(iE))
    enddo
    do iE = i_start, i_stop
      call destroy(Sigma_r_E(iE))
    enddo

   enddo !Loop on blocks 
  
   call log_deallocate(temp1)
   deallocate(Gn_E)
   deallocate(Sigma_r_E)

END SUBROUTINE complete_sigma_ph_r 

  
! READ Matrices
SUBROUTINE read_blkmat(Matrix, path, name, i, j, iE)

    TYPE(z_DNS) :: Matrix
    CHARACTER(*) :: path
    CHARACTER(*) :: name
    INTEGER :: i, j, iE

    CHARACTER(64) :: filename
    character(4) :: ofblki, ofblkj
    character(10) :: ofpnt
    logical :: lex
    integer :: i1,i2
    complex(dp) :: mat_el

    Matrix%val = (0.d0,0.d0)

    if (i.le.9999) write(ofblki,'(i4.4)') i
    if (i.gt.9999) stop 'ERROR: too many blks (> 9999)'
    if (j.le.9999) write(ofblkj,'(i4.4)') j 
    if (j.gt.9999) stop 'ERROR: too many blks (> 9999)'

    if (iE.le.99999) write(ofpnt,'(i5.5)') iE

    filename = trim(name)//trim(ofblki)//'_'//trim(ofblkj)//'_'//trim(ofpnt)//'.dat'
    
    inquire(file=trim(path)//trim(filename),EXIST=lex)
    if (.not.lex) then
       RETURN
       !WRITE(*,*) 'ERROR: FILE '//trim(filename)//' DOES NOT EXIST'
       !STOP
    endif   

    open(9091,file=trim(path)//trim(filename), form='UNFORMATTED')
    
    call inmat_c(9091,.false.,Matrix%val,Matrix%nrow,Matrix%ncol)

    !open(9001,file=trim(path)//trim(filename), access='DIRECT', recl=4)

    !call direct_in_c(9001,Matrix%val,Matrix%nrow)

    close(9091)

  END SUBROUTINE read_blkmat

  ! WRITE Matrices
  SUBROUTINE write_blkmat(Matrix, path, name, i, j, iE)

    TYPE(z_DNS) :: Matrix
    CHARACTER(*) :: path
    CHARACTER(*) :: name
    INTEGER :: i, j, iE

    CHARACTER(64) :: filename
    character(4) :: ofblki, ofblkj
    character(10) :: ofpnt

    if (i.le.9999) write(ofblki,'(i4.4)') i
    if (i.gt.9999) stop 'ERROR: too many blks (> 9999)'
    if (j.le.9999) write(ofblkj,'(i4.4)') j
    if (j.gt.9999) stop 'ERROR: too many blks (> 9999)'

    if (iE.le.99999) write(ofpnt,'(i5.5)') iE

    filename = trim(name)//trim(ofblki)//'_'//trim(ofblkj)//'_'//trim(ofpnt)//'.dat'

    open(9001,file=trim(path)//trim(filename), form='UNFORMATTED', status='REPLACE')

    call outmat_c(9001,.false.,Matrix%val,Matrix%nrow,Matrix%ncol) !,1.0d-36)
    
    !open(9001,file=trim(path)//trim(filename), status='REPLACE', access='DIRECT', recl=4)
   
    !call direct_out_c(9001,Matrix%val,Matrix%nrow)

    close(9001)

  END SUBROUTINE write_blkmat

  !****************************************************************************
  !
  !  Calculate Gn contributions for all contacts but collector, in the 
  !  contacts regions, where overlap with device orbitals is non-zero
  !  writing on memory
  !
  !****************************************************************************

  !SUBROUTINE Outer_A_mem_dns(Tlc,Tcl,gsurfR,struct,Aout)

    !****************************************************************************
    !Input:
    !Tlc: sparse arraymatrices array containing contact-device interacting blocks 
    !     ESH-H
    !Tlc: sparse arraymatrices array containing device-contact interacting blocks 
    !     ESH-H
    !gsurfR: sparse matrices array containing contacts surface green
    !
    !global variables needed: nbl, indblk(nbl+1), ncont, cblk(ncont), 
    !cindblk(ncont), Gr in the diagonal block related to the interaction layer
    ! 
    !Output:
    !Aout: sparse matrix containing density matrix in the region 
    !      corresponding to non-zero overlap
    !
    !****************************************************************************

!    IMPLICIT NONE 

    !In/Out
!    TYPE(z_CSR), DIMENSION(:) :: Tlc,Tcl,gsurfR
!    TYPE(Tstruct_Info), intent(in) :: struct
!    TYPE(z_CSR), intent(out) :: Aout

    !Work
!    TYPE(z_CSR) :: work1,GrCSR
!    TYPE(z_CSR) :: Asub, Grlc, Grcl
!    INTEGER :: i,cb,i1,j1
!    INTEGER :: ncont,nrow_tot,nbl
!    INTEGER, DIMENSION(:), POINTER :: indblk, cblk

!    ncont = struct%num_conts
!    nbl = struct%num_PLs
!    nrow_tot = struct%total_dim
!    indblk => struct%mat_PL_start
!    cblk => struct%cblk

    !Righe totali del conduttore effettivo nrow_tot 
    !nrow_tot=indblk(nbl+1)-1
    !DO i=1,ncont
    !   nrow_tot=nrow_tot+ncdim(i)   !gsurfR(i)%nrow
    !ENDDO
!    CALL create(Aout,nrow_tot,nrow_tot,0)
!    Aout%rowpnt(:)=1



!    DO i=1,ncont

       !Numero di blocco del contatto
!       cb=cblk(i)

       !converto Gr(cb,cb) da denso a sparso
!       call create(GrCSR,Gr(cb,cb)%nrow,Gr(cb,cb)%ncol,Gr(cb,cb)%nrow*Gr(cb,cb)%ncol)
!       call dns2csr(Gr(cb,cb),GrCSR)

       !Calcolo di Grlc
!       CALL prealloc_mult(GrCSR,Tlc(i),(-1.d0, 0.d0),work1)
       !Nota: numero colonne di Tlc = numero di righe di gsurf(i)

!       CALL prealloc_mult(work1,gsurfR(i),Grlc)
!       CALL destroy(work1)

       !Calcolo di Grcl
!       CALL prealloc_mult(gsurfR(i),Tcl(i),(-1.d0, 0.d0),work1)

!       CALL prealloc_mult(work1,GrCSR,Grcl)
!       CALL destroy(work1)

       !Calcolo della spectral density
!       CALL zspectral(Grlc,Grcl,0,Asub)

       !Dellocazione delle Green Retarded corrispondenti
!       CALL destroy(Grlc)
!       CALL destroy(Grcl)

       !Concatenazione di Asub nella posizione corrispondente
!       i1=indblk(cb)
!       j1=struct%mat_B_start(i)-struct%central_dim+indblk(nbl+1)-1
!       CALL concat(Aout,Asub,i1,j1)
!       CALL destroy(Asub)  
       
!       call destroy(GrCSR)

!    ENDDO



!    if (debug) then
!       WRITE(*,*) '********************'
!       WRITE(*,*) 'Outer_A_mem done'
!       WRITE(*,*) '********************'
!    endif

!  END SUBROUTINE Outer_A_mem_dns


  !****************************************************************************
  !
  !  Calculate Green Retarded in the 
  !  contacts regions, where overlap with device orbitals is non-zero
  !  (only the upper part, needed in both K=0 and K points calculations,
  !   or both upper and lower parts)
  !  writing on memory
  !
  !****************************************************************************

  SUBROUTINE Outer_Gr_mem_dns(Tlc,Tcl,gsurfR,struct,lower,Aout)

    !****************************************************************************
    !Input:
    !Tlc: sparse arraymatrices array containing contact-device interacting blocks 
    !     ESH-H
    !Tlc: sparse arraymatrices array containing device-contact interacting blocks 
    !     ESH-H
    !gsurfR: sparse matrices array containing contacts surface green
    !lower: if .true., also lower parts are calculated and concatenated
    !
    !global variables needed: nbl, indblk(nbl+1), ncont, cblk(ncont), 
    !cindblk(ncont), Gr in the diagonal block related to the interaction layer
    ! 
    !Output:
    !Aout: sparse matrix containing density matrix in the region 
    !      corresponding to non-zero overlap
    !
    !****************************************************************************

    IMPLICIT NONE 

    !In/Out
    TYPE(z_DNS), DIMENSION(:) :: Tlc,Tcl,gsurfR
    LOGICAL :: lower
    TYPE(Tstruct_info), intent(in) :: struct
    TYPE(z_CSR), intent(out) :: Aout


    !Work
    TYPE(z_DNS) :: work1, Grcl, Grlc 
    TYPE(z_CSR) :: GrCSR, TCSR
    INTEGER :: i,cb,nrow_tot,i1,j1
    INTEGER :: ncont, nbl
    INTEGER, DIMENSION(:), POINTER :: indblk, cblk

    ncont = struct%num_conts
    nbl = struct%num_PLs
    nrow_tot = struct%total_dim
    indblk => struct%mat_PL_start
    cblk => struct%cblk

    !Righe totali del conduttore effettivo nrow_tot 
    !nrow_tot=indblk(nbl+1)-1
    !DO i=1,ncont
    !   nrow_tot=nrow_tot+ncdim(i) !gsurfR(i)%nrow
    !ENDDO

    CALL create(Aout,nrow_tot,nrow_tot,0)

    Aout%rowpnt(:)=1

    DO i=1,ncont

       !Numero di blocco del contatto
       cb=cblk(i)
       !converto Gr(cb,cb) da denso a sparso
      !call create(GrCSR,Gr(cb,cb)%nrow,Gr(cb,cb)%ncol,Gr(cb,cb)%nrow*Gr(cb,cb)%ncol)
      !call dns2csr(Gr(cb,cb),GrCSR)
       !Calcolo di Grlc
      !CALL prealloc_mult(GrCSR,Tlc(i),(-1.d0, 0.d0),work1)
       CALL prealloc_mult(Gr(cb,cb),Tlc(i),(-1.d0, 0.d0),work1)
       CALL prealloc_mult(work1,gsurfR(i),Grlc)

       CALL destroy(work1)

       j1=nzdrop(Grlc,EPS) 
       CALL create(GrCSR,Grlc%nrow,Grlc%ncol,j1)
       CALL dns2csr(Grlc,GrCSR)
       CALL destroy(Grlc)
       j1=nzdrop(Tlc(i),EPS)
       CALL create(TCSR,Tlc(i)%nrow,Tlc(i)%ncol,j1)
       CALL dns2csr(Tlc(i),TCSR)
       CALL zmask_realloc(GrCSR,TCSR)
       CALL destroy(TCSR)
      !CALL zmask_realloc(Grlc, Tlc(i))

       !Concatenazione di Asub nella posizione corrispondente
       i1=indblk(cb)
       j1=struct%mat_B_start(i)-struct%central_dim+indblk(nbl+1)-1

      ! !CALL concat(Aout,Grlc,i1,j1)
       CALL concat(Aout,GrCSR,i1,j1)   

       CALL destroy(GrCSR)

       IF (lower) THEN

          CALL prealloc_mult(gsurfR(i),Tcl(i),(-1.d0, 0.d0), work1)
          CALL prealloc_mult(work1, Gr(cb,cb), Grcl)
        !  !CALL prealloc_mult(work1, GrCSR, Grcl)
        !  !CALL zmask_realloc(Grcl, Tcl(i))

          CALL destroy(work1)

          j1=nzdrop(Grcl,EPS) 
          CALL create(GrCSR,Grcl%nrow,Grcl%ncol,j1)
          CALL dns2csr(Grcl,GrCSR)
          CALL destroy(Grcl)
          j1=nzdrop(Tcl(i),EPS)
          CALL create(TCSR,Tcl(i)%nrow,Tcl(i)%ncol,j1)
          CALL dns2csr(Tcl(i),TCSR)
          CALL zmask_realloc(GrCSR,TCSR)
          CALL destroy(TCSR)
        
          i1 = struct%mat_B_start(i)-struct%central_dim+indblk(nbl+1)-1
          j1 = indblk(cb)

         ! !CALL concat(Aout,Grcl,i1,j1)
          CALL concat(Aout,GrCSR,i1,j1)

          CALL destroy(GrCSR)

       ENDIF

       ! !CALL destroy(GrCSR)
      
    ENDDO

    if (debug) then
       WRITE(*,*) '********************'
       WRITE(*,*) 'Outer_GreenR_mem done'
       WRITE(*,*) '********************'
    endif

  END SUBROUTINE Outer_Gr_mem_dns



  !****************************************************************************
  !
  !  Calculate Gn contributions for all contacts except collector, in the
  !  outer region where contacts-device overlap is non-zero - writing on 
  !  memory
  !
  !****************************************************************************

  SUBROUTINE Outer_Gn_mem_dns(Tlc,gsurfR,SelfEneR,struct,frm,ref,lower,Glout)

    !****************************************************************************
    !Input:
    !Tlc: sparse matrices array containing contacts-device interaction blocks (ES-H)
    !gsurfR: sparse matrices array containing contacts surface green
    !SelfEneR: sparse matrices array containing contacts Self Energy
    !frm: array containing Fermi distribution values for all contacts
    !ref: reference contact
    !
    !global variables needed: nbl, indblk(nbl+1), cindblk(ncont), ncont, 
    !cblk(ncont), Gr(:,:), diagonal, subdiagonal, overdiagonal and 
    !in colums Gr(:,cb) where cb are layers interacting with all contacts 
    !but collector 
    !
    !Output:
    !Glout: sparse matrix containing G_n  contributions in the region 
    !       corresponding to non-zero overlap
    !
    !****************************************************************************

    IMPLICIT NONE

    !In/Out
    TYPE(z_DNS), DIMENSION(:) :: Tlc, gsurfR, SelfEneR
    REAL(dp), DIMENSION(:) :: frm
    TYPE(Tstruct_info), intent(in) :: struct
    INTEGER :: ref 
    LOGICAL :: lower
    TYPE(z_CSR) :: Glout

    !Work
    TYPE(z_DNS) :: Gam, gsurfA, Ga, work1, work2, work3, Glsub
    TYPE(z_CSR) :: GlCSR, TCSR
    INTEGER :: j,k,cb,cbj,i1,j1,nrow_tot
    INTEGER :: ncont, nbl
    INTEGER, DIMENSION(:), POINTER :: indblk, cblk
    COMPLEX(dp) :: frmdiff
    TYPE(z_CSR) :: Gr_sparse 


    ncont = struct%num_conts
    nbl = struct%num_PLs
    nrow_tot = struct%total_dim
    indblk => struct%mat_PL_start
    cblk => struct%cblk
   
    !Allocazione della Glout
    !Righe totali del conduttore effettivo nrow_tot 
    !nrow_tot=indblk(nbl+1)-1
    !DO i=1,ncont
    !   nrow_tot=nrow_tot+ncdim(i) !gsurfR(i)%nrow
    !ENDDO
    CALL create(Glout,nrow_tot,nrow_tot,0)
    Glout%rowpnt(:)=1

    !***
    !Iterazione su tutti i contatti "k" 
    !***
    DO k=1,ncont

       !Esegue le operazioni relative al contatto solo se e` valida la condizione
       !sulle distribuzioni di Fermi e se non si tratta del contatto iniettante (ref)
       IF ((ABS(frm(k)-frm(ref)).GT.drop).AND.(k.NE.ref)) THEN

          !Calcolo della Gamma corrispondente
          cb=cblk(k)
          !nota:cb indica l'indice di blocco corrispondente al contatto j-esimo
          !nrow_cont=ncdim(k) !gsurfR(k)%nrow          
          !nrow_cb=indblk(cb+1)-indblk(cb)

          !Conversione densa-sparsa della green
         ! !call create(Gr_sparse,Gr(cb,cb)%nrow,Gr(cb,cb)%ncol,Gr(cb,cb)%nrow*Gr(cb,cb)%ncol)
         ! !call dns2csr(Gr(cb,cb),Gr_sparse)

          CALL zspectral(SelfEneR(k),SelfEneR(k),0,Gam)

          !***
          !Calcolo del contributo sulla proria regione
          !***          
          !Primo addendo

          frmdiff=cmplx(frm(ref)-frm(k))

          CALL zspectral(gsurfR(k),gsurfR(k),0,work1)
        
          !print*, 'work2=Tlc*work1=Tlc*j(gsurfR-gsurfA)'
          CALL prealloc_mult(Tlc(k),work1,work2)
          CALL destroy(work1)

          !print *, 'work1=-Gr(cb,cb)*work2=-Gr(cb,cb)*Tlc*j(gsurfR-gsurfA)'
         ! !CALL prealloc_mult(Gr_sparse,work2,frmdiff,work1)
          CALL prealloc_mult(Gr(cb,cb),work2,frmdiff,work1)
          CALL destroy(work2)

          !if (debug) write(*,*) 'cont',k,'primo addendo ok'
          !Secondo addendo

          !print*,'work2=Tlc*gsurfA'
          CALL zdagger(gsurfR(k),gsurfA)
          CALL prealloc_mult(Tlc(k),gsurfA,work2)
          CALL destroy(gsurfA)

          !print*,'work3=Ga*work2=Ga*Tlc*gsurfA'           
         ! !CALL zdagger(Gr_sparse,Ga)
          CALL zdagger(Gr(cb,cb),Ga)
          CALL prealloc_mult(Ga,work2,work3)

          CALL destroy(Ga)
          CALL destroy(work2)

          !print*,'work2=Gam*work3=Gam*Ga*Tlc*gsurfA'          
          CALL prealloc_mult(Gam,work3,work2)
          CALL destroy(work3)

          !print *,'work3=-Gr*work2=-Gr*Gam*Ga*Tlc*gsurfA'        
          CALL prealloc_mult(Gr(cb,cb),work2,frmdiff,work3)
          CALL destroy(work2)

          !if (debug) write(*,*) 'cont',k,'secondo addendo ok'

          !Contributo totale sulla propria regione
          CALL prealloc_sum(work3,work1,Glsub)
          CALL destroy(work1)
          CALL destroy(work3)

          j1=nzdrop(Glsub,EPS) 
          CALL create(GlCSR,Glsub%nrow,Glsub%ncol,j1)
          CALL dns2csr(Glsub,GlCSR)
          CALL destroy(Glsub)
          j1=nzdrop(Tlc(k),EPS)
          CALL create(TCSR,Tlc(k)%nrow,Tlc(k)%ncol,j1)
          CALL dns2csr(Tlc(k),TCSR)
          CALL zmask_realloc(GlCSR,TCSR)
          CALL destroy(TCSR)
          
         ! !CALL zmask_realloc(Glsub, Tlc(k)) 

          !Concatenazione di Glsub nella matrice globale Glout
          i1=indblk(cb)
          j1=struct%mat_B_start(k)-struct%central_dim+indblk(nbl+1)-1
          CALL concat(Glout,GlCSR,i1,j1)

          ! compute lower outer part using (iG<)+ = iG<    
          IF (lower) THEN
             call zdagger(GlCSR,TCSR)
             call concat(Glout,TCSR,j1,i1)
             call destroy(TCSR)
          ENDIF

          CALL destroy(GlCSR)
        ! !CALL destry(Gr_sparse)

          !if (debug) write(*,*) 'cont',k,'concat ok'
          !***
          !Ciclo per il calcolo dei contributi sulle regioni degli altri contatti
          !***

          DO j=1,ncont

             cbj=cblk(j)
             !nrow_cbj=indblk(cbj+1)-indblk(cbj)
             !nrow_contj=ncdim(j) !gsurfR(j)%nrow
             !Esegue le operazioni del ciclo solo se il j.ne.k o se
             !il blocco colonna di Gr e` non nullo (altrimenti il contributo e` nullo) 

             IF ((j.NE.k).AND.(Gr(cbj,cb)%nrow.NE.0 .AND. (Gr(cbj,cb)%ncol.NE.0))) THEN

               ! !call create(Gr_sparse,Gr(cbj,cb)%nrow,Gr(cbj,cb)%ncol, &
               ! !            Gr(cbj,cb)%nrow*Gr(cbj,cb)%ncol)
               ! !call dns2csr(Gr(cbj,cb),Gr_sparse)

                !print*,'work1=Tlc*gsurfA'  
                CALL zdagger(gsurfR(j),gsurfA)
                CALL prealloc_mult(Tlc(j),gsurfA,work1)
                CALL destroy(gsurfA)

                !if (debug) write(*,*) 'cont',j,'T*g'

                !print*,'work2=Ga*work1=Ga*Tlc*gsurfA'  
               ! !CALL zdagger(Gr_sparse,Ga)
                CALL zdagger(Gr(cbj,cb),Ga)
                CALL prealloc_mult(Ga,work1,work2)

                !if (debug) write(*,*) 'cont',j,'Ga*T*g'

                CALL destroy(Ga)
                CALL destroy(work1)

                !if (debug) write(*,*) 'Gam',Gam%nrow,Gam%ncol,Gam%nnz
                !if (debug) write(*,*) 'work2',work2%nrow,work2%ncol,work2%nnz
                !print*,'work1=Gam*work2=Gam*Ga*Tlc*gsurfA'  
                CALL prealloc_mult(Gam,work2,work1)
                CALL destroy(work2)

                !if (debug) write(*,*) 'cont',j,'Gam*Ga*T*g'

                !if (debug) write(*,*) 'work1',work1%nrow,work1%ncol,work1%nnz
                !if (debug) write(*,*) 'Gr',Gr(cbj,cb)%nrow,Gr(cbj,cb)%ncol,Gr(cbj,cb)%nnz

                !print*,'Glsub=-Gr*work1=-Gr*Gam*Ga*Tlc*gsurfA'  
                CALL prealloc_mult(Gr(cbj,cb),work1,frmdiff,Glsub)
                CALL destroy(work1)

                j1=nzdrop(Glsub,EPS) 
                CALL create(GlCSR,Glsub%nrow,Glsub%ncol,j1)
                CALL dns2csr(Glsub,GlCSR)
                CALL destroy(Glsub)
                j1=nzdrop(Tlc(j),EPS)
                CALL create(TCSR,Tlc(j)%nrow,Tlc(j)%ncol,j1)
                CALL dns2csr(Tlc(j),TCSR)
                CALL zmask_realloc(GlCSR,TCSR)
                CALL destroy(TCSR)
             
                ! !CALL zmask_realloc(Glsub, Tlc(j))
 
                !if (debug) write(*,*) 'cont',j,'Gr*Gam*Ga*T*g'

                !Concatenazione di Glsub nella posizione corrispondente al contatto "j"
                i1=indblk(cbj)
                j1=struct%mat_B_start(j)-struct%central_dim+indblk(nbl+1)-1
                CALL concat(Glout,GlCSR,i1,j1)

                ! compute lower outer part using (iG<)+ = iG<    
                IF (lower) THEN
                   call zdagger(GlCSR,TCSR)
                   call concat(Glout,TCSR,j1,i1)
                   call destroy(TCSR)
                ENDIF
               
                CALL destroy(GlCSR)
                !CALL destroy(Gr_sparse)
                !if (debug) write(*,*) 'cont',k,'concat ok'

             ENDIF

          ENDDO

          CALL destroy(Gam)

       ENDIF

    ENDDO

    if (debug) then
       WRITE(*,*) '********************'
       WRITE(*,*) 'Outer_Gl_mem done'
       WRITE(*,*) '********************'
    endif

  END SUBROUTINE Outer_Gn_mem_dns


  !---------------------------------------------------

  SUBROUTINE tunneling_dns(HM,SM,Ec,SelfEneR,ni,nf,size_ni,str,TUN_MAT)

    implicit none

    Type(z_CSR) :: HM
    Type(z_CSR) :: SM           
    Complex(dp) :: Ec
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR
    !Type(z_DNS) :: SelfEner_d 
    Real(dp), Dimension(:) :: TUN_MAT
    Type(z_CSR) :: ESH_tot
    Type(TStruct_Info) :: str

    ! Local variables
    Type(z_DNS), Dimension(:,:), allocatable :: ESH
    Real(dp) :: tun
    Integer :: ni(MAXNCONT)
    Integer :: nf(MAXNCONT)
    Integer :: nbl,ncont,size_ni
    Integer :: i, ierr, icpl, nit, nft, nt, nt1

    nbl = str%num_PLs
    ncont = str%num_conts

    !Calculation of ES-H and brak into blocks
    call prealloc_sum(HM,SM,(-1.d0, 0.d0),Ec,ESH_tot)    

    allocate(ESH(nbl,nbl),stat=ierr)
    if (ierr.eq.1) stop 'Error in ESH allocation'
    call sub_ESH_dns(ESH_tot,ESH,str%mat_PL_start)
    call destroy(ESH_tot)

    !Inclusion of the contact Self-Energies to the relevant blocks
    do i=1,ncont
       !call create(SelfEner_d, SelfEner(i)%nrow, SelfEner(i)%ncol)
       !call csr2dns(SelfEneR(i),SelfEneR_d)
       ESH(str%cblk(i),str%cblk(i))%val = ESH(str%cblk(i),str%cblk(i))%val-SelfEneR(i)%val
       !call destroy(SelfEneR_d)
    enddo
    call allocate_gsmr_dns(nbl)

    !Iterative calculation up with gsmr 
    call Make_gsmr_mem_dns(ESH,nbl,2)

    !Iterative calculation down for Gr
    call allocate_Gr_dns(nbl)

    !Computation of transmission(s) between contacts ni(:) -> nf(:)  
    nit=ni(1)
    nft=nf(1)
  
    !Arrange contacts in a way that the order between first and second is always the
    !same (always ct1 < ct2)
    
    if (str%cblk(nit).lt.str%cblk(nft)) then
       nt1 = str%cblk(nit)
    else
       nt1 = str%cblk(nft)
    endif
   
    ! Build Gr up to the lowest contact  block
    call Make_Grdiag_mem_dns(ESH,str%mat_PL_start,nt1)
    
    if (size_ni.eq.1) then
       call trasmission_dns(nit,nft,ESH,SelfEneR,nbl,str%cblk,str%mat_PL_start,tun) 
    else
       call trasmission_old(nit,nft,ESH,SelfEneR,nbl,str%cblk,str%mat_PL_start,tun) 
    endif

    TUN_MAT(1) = tun 

    ! When more contacts are present sometimes we can re-use previous Gr 
    do icpl = 2, size_ni
       
       nit=ni(icpl)
       nft=nf(icpl)

       if (str%cblk(nit).lt.str%cblk(nft)) then
          nt = str%cblk(nit)
       else
          nt = str%cblk(nft)
       endif

       if (nt .ne. nt1) then
          !Distruzione dei blocchi fuori-diagonale
          do i=2,nt1
             call destroy(Gr(i,i))
             call destroy(Gr(i-1,i))
             call destroy(Gr(i,i-1))
          enddo
          call destroy(Gr(1,1))

          call Make_Grdiag_mem_dns(ESH,str%mat_PL_start,nt)
       endif

       if (size_ni.eq.1) then
          call trasmission_dns(nit,nft,ESH,SelfEneR,nbl,str%cblk,str%mat_PL_start,tun) 
       else
          call trasmission_old(nit,nft,ESH,SelfEneR,nbl,str%cblk,str%mat_PL_start,tun) 
       endif
  
       TUN_MAT(icpl) = tun 
    
       nt1 = nt

    enddo

    !Distruzione delle Green
    do i=2,nt1
       call destroy(Gr(i,i))
       call destroy(Gr(i-1,i))
       call destroy(Gr(i,i-1))
    enddo
    call destroy(Gr(1,1))

    call deallocate_Gr_dns

    !Deallocate energy-dependent matrices
    do i=2,nbl 
       call destroy(gsmr(i))
    enddo

    call deallocate_gsmr_dns

    !Deallocate matrices

    call destroy(ESH(1,1))
    do i=2,nbl
       call destroy(ESH(i,i))
       call destroy(ESH(i-1,i))
       call destroy(ESH(i,i-1))
    enddo

    deallocate(ESH,stat=ierr)
    if (ierr.ne.0) stop 'DEALLOCATION ERROR: could not deallocate ESH'

  end SUBROUTINE tunneling_dns

  !************************************************************************
  !
  ! Subroutine for transmission calculation
  !
  !************************************************************************
  ! NOTE:
  !
  !  This subroutine was hacked quickly to obain effiecient tunneling calcs
  !  Useful only when there are 2 contacts
  !                ===================
  !************************************************************************


  subroutine trasmission_dns(ni,nf,ESH,SelfEneR,nbl,cblk,indblk,TUN)

    implicit none
    
    !In/Out
    Integer :: ni,nf, nbl
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR
    Type(z_DNS), Dimension(:,:) :: ESH
    Integer, Dimension(:), pointer :: cblk, indblk
    Real(dp) :: TUN
    
    !Work variables
    Integer :: ct1, nt1, i
    Type(z_DNS) :: work1, work2, GAM1_dns, GA, TRS, AA
    Type(z_CSR) :: GAM1
    Complex(dp), PARAMETER ::    j = (0.d0,1.d0)  ! CMPX unity
   
    if (size(cblk).gt.2) then
       write(*,*) "ERROR: transmission_dns is valid only for 2 contacts"
       return
    endif

    !Arrange contacts in way that order between first and second is always the
    !same (always ct1 < ct2)

    if (cblk(ni).lt.cblk(nf)) then
       ct1=ni;
    else
       ct1=nf;
    endif

    nt1=cblk(ct1); 

    call zdagger(Gr(nt1,nt1),GA)

    ! Computes the Gamma matrices
    call zspectral(SelfEneR(ct1),SelfEneR(ct1),0,GAM1_dns)

    !call create(GAM1_dns,GAM1%nrow,GAM1%ncol)

    !call csr2dns(GAM1,GAM1_dns)

    !call destroy(GAM1)
    
    ! Work to compute transmission matrix (Gamma G Gamma G)
    call prealloc_mult(GAM1_dns,Gr(nt1,nt1),work1)
    
    call prealloc_mult(work1,GAM1_dns,work2)
    
    call destroy(work1)

    call prealloc_mult(work2,GA,work1)
 
    call destroy(work2)

    call create(AA,GA%nrow,GA%ncol)

    AA%val = j * (Gr(nt1,nt1)%val-GA%val)

    call destroy(GA) 

    call prealloc_mult(GAM1_dns,AA,work2) 

    call destroy(GAM1_dns,AA)

    call create(TRS,work1%nrow,work1%ncol)

    TRS%val = work2%val - work1%val

    TUN = abs( real(trace(TRS)) )  

    call destroy(TRS,work1,work2)
   
  end subroutine trasmission_dns

  !************************************************************************
  !
  ! Subroutine for transmission calculation (GENERIC FOR N CONTACTS)
  !
  !************************************************************************ 
  subroutine trasmission_old(ni,nf,ESH,SelfEneR,nbl,cblk,indblk,TUN)
    
    !In/Out
    Integer :: ni,nf
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR
    Type(z_DNS), Dimension(:,:) :: ESH
    Integer, Dimension(:), pointer :: cblk, indblk
    Real(kind=dp) :: TUN
    
    !Work variables
    Integer :: ct1, ct2, nt1, nt2, i, nrow, ncol, nbl
    Type(z_DNS) :: work1, work2, GAM1_dns, GAM2_dns, GA, TRS, AA
    Type(z_CSR) :: GAM1, GAM2
    Real(kind=dp) :: max
   
    !Arrange contacts in way that order between first and second is always the
    !same (always ct1 < ct2)

    if (cblk(ni).lt.cblk(nf)) then
       ct1=ni;ct2=nf;
    else
       ct1=nf;ct2=ni;
    endif
    
    nt1=cblk(ct1); nt2=cblk(ct2);

    ! in this way nt1 < nt2 by construction

    ncol=indblk(nt1+1)-indblk(nt1)
    
    if ( nbl.gt.1 .and. (nt2-nt1).gt.1) then
       
       !Calcolo dei blocchi colonna (nt1) fino al desiderato (riga nt2)

       do i = nt1+1, nbl

          !Checks whether previous block is non null. 
          !If so next block is also null => TUN = 0       
          max=maxval(abs(Gr(i-1,nt1)%val))

          if (max.lt.drop) then
             TUN = 0.d0
             !Destroy also the block adjecent to diagonal since 
             !this is not deallocated anymore in calling subroutine
             if (i.gt.(nt1+1)) call destroy(Gr(i-1,nt1))
             return
          endif


          !Checks whether block has been created, if not do it  
          if (Gr(i,nt1)%nrow.eq.0 .or. Gr(i,nt1)%ncol.eq.0) then 

             call prealloc_mult(gsmr(i),ESH(i,i-1),(-1.d0, 0.d0),work1)
             
             call prealloc_mult(work1,Gr(i-1,nt1),Gr(i,nt1))
             
             call destroy(work1)

          endif
          
          !Destroy also the block adjecent to diagonal since 
          !this is not deallocated anymore in calling subroutine
          if (i.gt.(nt1+1)) call destroy(Gr(i-1,nt1))

          if (i.eq.nt2) exit

       enddo
       
    endif

    ! Computes the Gamma matrices
    call zspectral(SelfEneR(ct1),SelfEneR(ct1),0,GAM1_dns)
    call zspectral(SelfEneR(ct2),SelfEneR(ct2),0,GAM2_dns)
    
    !call create(GAM1_dns,GAM1%nrow,GAM1%ncol)

    !call csr2dns(GAM1,GAM1_dns)
    
    !call destroy(GAM1)

    !call create(GAM2_dns,GAM2%nrow,GAM2%ncol)

    !call csr2dns(GAM2,GAM2_dns)
    
    !call destroy(GAM2)

    ! Work to compute transmission matrix (Gamma2 Gr Gamma1 Ga)
    call prealloc_mult(GAM2_dns,Gr(nt2,nt1),work1)

    call destroy(GAM2_dns)

    call prealloc_mult(work1,GAM1_dns,work2)

    call destroy(work1)

    call destroy(GAM1_dns)

    call zdagger(Gr(nt2,nt1),GA)
    
    if (nt2.gt.nt1+1) call destroy( Gr(nt2,nt1) )

    call prealloc_mult(work2,GA,TRS)

    call destroy(work2)

    call destroy(GA) 
  
    TUN = real(trace(TRS))

    call destroy(TRS)
    
  end subroutine trasmission_old
 

  !---------------------------------------------------!
  !Subroutine for transmission and dos calculation    !
  !---------------------------------------------------!  

  subroutine tun_and_dos(HM,SM,Ec,SelfEneR,Gs,ni,nf,nLDOS,LDOS,size_ni,str,TUN_MAT,LEDOS)

    implicit none

    Type(z_CSR), intent(in) :: HM
    Type(z_CSR), intent(in) :: SM           
    Complex(dp), intent(in) :: Ec
    Type(z_DNS), Dimension(MAXNCONT), intent(in) :: SelfEneR, Gs
    Type(TStruct_Info), intent(in) :: str
    integer, intent(in)  :: nLdos, size_ni
    integer, Dimension(:,:), pointer :: LDOS      
    Real(dp), Dimension(:), intent(inout) :: TUN_MAT
    Real(dp), Dimension(:), intent(inout) :: LEDOS

    ! Local variables
    Type(z_CSR) :: ESH_tot, GrCSR
    Type(z_DNS) :: SelfEner_d 
    Type(z_DNS), Dimension(:,:), allocatable :: ESH
    Type(r_CSR) :: Grm                          ! Green Retarded nella molecola           
    Real(dp) :: tun
    Complex(dp) :: zc
    Integer :: ni(MAXNCONT)
    Integer :: nf(MAXNCONT)
    Integer :: nbl,ncont, ierr
    Integer :: nit, nft, icpl
    Integer :: iLDOS, i2, i
    Character(1) :: Im

    nbl = str%num_PLs
    ncont = str%num_conts
    Im = 'I'

    !Calculation of ES-H and brak into blocks
    call prealloc_sum(HM,SM,(-1.d0, 0.d0),Ec,ESH_tot)    

    allocate(ESH(nbl,nbl),stat=ierr)
    if (ierr.eq.1) stop 'Error in ESH allocation'
    call sub_ESH_dns(ESH_tot,ESH,str%mat_PL_start)
    call destroy(ESH_tot)

    !Inclusion of the contact Self-Energies to the relevant blocks
    do i=1,ncont
       !call create(SelfEner_d, SelfEner(i)%nrow, SelfEner(i)%ncol)
       !call csr2dns(SelfEneR(i),SelfEneR_d)
       ESH(str%cblk(i),str%cblk(i))%val = ESH(str%cblk(i),str%cblk(i))%val-SelfEneR(i)%val
       !call destroy(SelfEneR_d)
    enddo

    call allocate_gsmr_dns(nbl)

    !Iterative calculation up with gsmr 
    call Make_gsmr_mem_dns(ESH,nbl,2)

    !Iterative calculation down for Gr
    call allocate_Gr_dns(nbl)

    call Make_Grdiag_mem_dns(ESH,str%mat_PL_start)

    !Computation of transmission(s) between contacts ni(:) -> nf(:)
    do icpl=1,size_ni

       nit=ni(icpl)
       nft=nf(icpl)

       call trasmission_dns(nit,nft,ESH,SelfEneR,nbl,str%cblk,str%mat_PL_start,tun) 
  
       TUN_MAT(icpl) = tun 
    
    enddo
   
    !Deallocate energy-dependent matrices
    do i=2,nbl 
       call destroy(gsmr(i))
    enddo

    call deallocate_gsmr_dns
    !Distruzione dei blocchi fuori-diagonale
    do i=2,nbl
       call destroy(Gr(i-1,i))
       call destroy(Gr(i,i-1))
    enddo 

    call create(Grm,str%mat_PL_start(nbl+1)-1,str%mat_PL_start(nbl+1)-1,0)
    
    Grm%rowpnt(:)=1

    do i=1,nbl
       call create(GrCSR,Gr(i,i)%nrow,Gr(i,i)%ncol,Gr(i,i)%nrow*Gr(i,i)%ncol)
       call dns2csr(Gr(i,i),GrCSR)
       !Concatena direttamente la parte immaginaria per il calcolo della DOS
       zc=(-1.d0,0.d0)/3.14159265358979323844_dp

       call concat(Grm,zc,GrCSR,Im,str%mat_PL_start(i),str%mat_PL_start(i))
       call destroy(Gr(i,i))
       call destroy(GrCSR)
    enddo

    !Sort Grm for fast trace   
    call msort(Grm)

    !Compute LDOS on the specified intervals
    do iLDOS=1,nLDOS

       if( LDOS(2,iLDOS).le.str%central_dim ) then

          do i2 = LDOS(1,iLDOS), LDOS(2,iLDOS)

             LEDOS(iLDOS)=LEDOS(iLDOS) + getelment(i2,i2,Grm)

          enddo

       endif

    enddo

    call destroy(Grm)

    call deallocate_Gr_dns

    call destroy(ESH(1,1))
    do i=2,nbl
       call destroy(ESH(i,i))
       call destroy(ESH(i-1,i))
       call destroy(ESH(i,i-1))
    enddo

    deallocate(ESH,stat=ierr)
    if (ierr.ne.0) stop 'DEALLOCATION ERROR: could not deallocate ESH'

  end subroutine tun_and_dos


  !---------------------------------------------------


  subroutine allocate_gsmr_dns(nbl)
    integer :: nbl, ierr

    allocate(gsmr(nbl),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate gsmr'
    
  end subroutine allocate_gsmr_dns

  !---------------------------------------------------

  subroutine allocate_gsml_dns(nbl)
    integer :: nbl, ierr

    allocate(gsml(nbl),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate gsml'
    
  end subroutine allocate_gsml_dns

  !---------------------------------------------------

  subroutine allocate_Gr_dns(nbl)
    integer :: nbl, ierr

    allocate(Gr(nbl,nbl),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate Gr'

  end subroutine allocate_Gr_dns


  !---------------------------------------------------

  subroutine deallocate_gsmr_dns
    integer :: ierr

    deallocate(gsmr,stat=ierr)
    if (ierr.ne.0) stop 'DEALLOCATION ERROR: could not deallocate gsmr'

  end subroutine deallocate_gsmr_dns

  !---------------------------------------------------

  subroutine deallocate_gsml_dns
    integer :: ierr

    deallocate(gsml,stat=ierr)
    if (ierr.ne.0) stop 'DEALLOCATION ERROR: could not deallocate gsml'

  end subroutine deallocate_gsml_dns

  !---------------------------------------------------

  subroutine deallocate_Gr_dns
    integer :: ierr

    deallocate(Gr,stat=ierr)
    if (ierr.ne.0) stop 'DEALLOCATION ERROR: could not deallocate Gr'

  end subroutine deallocate_Gr_dns

  !---------------------------------------------------

  !--------------------------------------------------------------------------
  subroutine compGreen(G,A,n)
    Type(z_DNS) :: A, G
    Integer :: n

    Integer :: sel, iter


    sel = 2
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


  end subroutine compGreen





END MODULE iterative_dns
