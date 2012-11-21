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
  USE ln_allocation
  USE mat_def
  USE sparsekit_drv
  USE inversions
  USE ln_structure, only : TStruct_Info
  USE lib_param, only : MAXNCONT, Tnegf
  USE outmatrix, only : outmat_c
  USE clock

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

  public :: sub_ESH_dns
  public :: rebuild_dns
  public :: Make_gsmr_mem_dns
  public :: Make_gsml_mem_dns
  public :: Make_Grdiag_mem_dns
  public :: Make_Grcol_mem_dns
  public :: Make_Spectral_mem_dns
  public :: Make_GreenR_mem2_dns  ! New. Best.
  public :: Make_Gl_mem_dns
  !public :: Outer_A_mem_dns
  public :: Outer_GreenR_mem_dns
  public :: Outer_Gl_mem_dns

  LOGICAL, PARAMETER :: debug=.false. 
  !Dropout value
  REAL(dp), PARAMETER :: drop=1e-20

  TYPE(z_DNS), DIMENSION(:), ALLOCATABLE :: gsmr
  TYPE(z_DNS), DIMENSION(:), ALLOCATABLE :: gsml
  TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Gr
  TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Glsub

CONTAINS

  !****************************************************************************
  ! 
  ! Driver for computing Equilibrium Green Retarded (Outer blocks included)
  ! writing on memory
  !
  !****************************************************************************

  SUBROUTINE calls_eq_mem_dns(H,S,E,SelfEneR,Tlc,Tcl,gsurfR,A,struct,outer)

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
    TYPE(z_CSR) :: H, S, A
    TYPE(z_DNS), DIMENSION(:) :: SelfEneR
    TYPE(z_DNS), DIMENSION(:) :: Tlc, Tcl, gsurfR
    !TYPE(z_DNS) :: SelfEner_d 
    COMPLEX(dp) :: E
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
    CALL prealloc_sum(S,H,E,(-1.d0, 0.d0),ESH_tot)

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
       CALL Outer_GreenR_mem_dns(Tlc,Tcl,gsurfR,struct,.FALSE.,A)   
    CASE(2)
       CALL Outer_GreenR_mem_dns(Tlc,Tcl,gsurfR,struct,.TRUE.,A) 
    END SELECT

    !CALL Make_Spectral_mem(struct,0,A)
    !if (debug) write(*,*) 'Compute GreenR'   
    !Nota: Il flag 0 indica che non teniamo le Gr calcolate 
    !E' necessario metterlo a 1 se dopo il calcolo della spectral 
    !effettuiamo il calcolo dei contributi aggiunti della Gless
    !(omesso perchè è implementato un driver a parte)
    !CALL Make_GreenR_mem(H%nrow,0,A)
    call Make_GreenR_mem2_dns(S,nbl,indblk,Ain)

    !Distruzione dell'array Gr
    CALL destroy_Gr()
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
  ! Driver for computing Gless contributions due to all contacts but reference:
  !
  !   Sum   [f_j(E)-f_r(E)] Gr Gam_j Ga
  !   j!=r
  !
  ! NOTE: The subroutine assumes that 
  !
  !****************************************************************************

  SUBROUTINE calls_neq_mem_dns(H,S,E,SelfEneR,Tlc,Tcl,gsurfR,struct,frm,ref,Glout,out)

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
    !Aout: Gless contributions (Device + Contacts overlap regions -> effective conductor)
    !
    !*****************************************************************************


    IMPLICIT NONE

    !In/Out
    TYPE(z_CSR) :: H,S,Glout
    TYPE(z_DNS), DIMENSION(:) :: SelfEneR, gsurfR, Tlc, Tcl
    !TYPE(z_DNS) :: SelfEner_d 
    REAL(dp) :: E
    TYPE(Tstruct_info) :: struct
    REAL(dp), DIMENSION(:) :: frm
    INTEGER :: ref
    INTEGER :: out

    !Work
    COMPLEX(dp) :: Ec
    INTEGER :: i,ierr,i1,ncont,nbl, lbl
    INTEGER, DIMENSION(:), POINTER :: cblk, indblk
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: ESH
    TYPE(z_CSR) :: ESH_tot, Gl
    LOGICAL :: destr, mask(MAXNCONT)

    nbl = struct%num_PLs
    ncont = struct%num_conts
    indblk => struct%mat_PL_start
    cblk => struct%cblk
    Ec=cmplx(E,0.d0,dp)

    !if (debug) write(*,*) '----------------------------------------------------'
    !if (debug) call writeMemInfo(6)
    !if (debug) call writePeakInfo(6)
    !if (debug) write(*,*) '----------------------------------------------------'

    !Costruiamo la matrice sparsa ESH
    CALL prealloc_sum(H,S,(-1.d0, 0.d0),Ec,ESH_tot)

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
    !if (debug) write(*,*) 'Compute Outer_Gless' 
    SELECT CASE (out)
    CASE(0)
    CASE(1)
       CALL Outer_Gl_mem_dns(Tlc,gsurfR,SelfEneR,struct,frm,ref,.false.,Glout)
    CASE(2)
       CALL Outer_Gl_mem_dns(Tlc,gsurfR,SelfEneR,struct,frm,ref,.true.,Glout)
    END SELECT

    !if (debug) write(*,*) '----------------------------------------------------'
    !if (debug) call writePeakInfo(6)
    !if (debug) write(*,*) '----------------------------------------------------'

    !Calcolo della Gless nel device
    !if (debug) write(*,*) 'Compute Gless' 
    !Allocazione degli array di sparse
    ALLOCATE(Glsub(nbl,nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Glsub(nbl,nbl)'

    CALL Make_Gl_mem_dns(ESH,SelfEneR,frm,ref,struct,S,Gl)
    
    call blk2csr(struct,S,Gl)

    DEALLOCATE(Glsub)

    !Distruzione dell'array Gr
    CALL destroy_Gr()
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
! Driver for computing Gless contributions including el-ph interactions
!
!   Sum   f_j(E) Gr Gam_j Ga +   Gr Sigma_ph< Ga
!    j
!
! NOTE: The subroutine assumes that
!
!****************************************************************************

SUBROUTINE calls_neq_ph(pnegf,E,SelfEneR,Tlc,Tcl,gsurfR,frm,ref,Glout,out,iter)

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
!Aout: Gless contributions (Device + Contacts overlap regions -> effective conductor)
!
!*****************************************************************************
IMPLICIT NONE

!In/Out
TYPE(Tnegf), pointer :: pnegf
TYPE(z_CSR) :: Glout
TYPE(z_DNS), DIMENSION(:) :: SelfEneR, gsurfR, Tlc, Tcl
REAL(dp) :: E
REAL(dp), DIMENSION(:) :: frm
INTEGER :: ref
INTEGER :: out
INTEGER :: iter

!Work
COMPLEX(dp) :: Ec
INTEGER :: i,ierr,i1,ncont,nbl, lbl
INTEGER, DIMENSION(:), POINTER :: cblk, indblk
TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: ESH
TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_ph_r, Sigma_ph_less
TYPE(z_CSR) :: ESH_tot, Gl
LOGICAL :: destr, mask(MAXNCONT)
TYPE(Tstruct_info) :: struct

struct = pnegf%str
nbl = struct%num_PLs
ncont = struct%num_conts
indblk => struct%mat_PL_start
cblk => struct%cblk
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
       !call create(SelfEner_d, SelfEner(i)%nrow, SelfEner(i)%ncol)
       !call csr2dns(SelfEneR(i),SelfEneR_d)
       ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
       !call destroy(SelfEneR_d)
 ENDDO


 ALLOCATE(Sigma_ph_r(nbl,nbl),stat=ierr)
 IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_r'
 ALLOCATE(Sigma_ph_less(nbl,nbl),stat=ierr)
 IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_less'

 ! Reload and add Sigma_ph_r to ESH
 call reload_r(pnegf, ESH, Sigma_ph_r, iter)


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
 !if (debug) write(*,*) 'Compute Outer_Gless' 
 SELECT CASE (out)
 CASE(0)
 CASE(1)
   CALL Outer_Gl_mem_dns(Tlc,gsurfR,SelfEneR,pnegf%str,frm,ref,.false.,Glout)
 CASE(2)
   CALL Outer_Gl_mem_dns(Tlc,gsurfR,SelfEneR,pnegf%str,frm,ref,.true.,Glout)
 END SELECT

!if (debug) write(*,*) '----------------------------------------------------'
!if (debug) call writePeakInfo(6)
!if (debug) write(*,*) '----------------------------------------------------'

 !Calcolo della Gless nel device
 !if (debug) write(*,*) 'Compute Gless' 
 !Allocazione degli array di sparse
 ALLOCATE(Glsub(nbl,nbl),stat=ierr)
 IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Glsub(nbl,nbl)'

 CALL Make_Gl_mem_dns(ESH,SelfEneR,frm,ref,struct,pnegf%SM,Gl)

 call Make_Gl_ph(pnegf,ESH,Sigma_ph_less,iter)

 ! save diagonal blocks of Gless
 DO i = 1, nbl
   call write_blkmat(Glsub(i,i),pnegf%scratch_path,'G_less_',i,i,pnegf%Epnt)
   call write_blkmat(Gr(i,i),pnegf%scratch_path,'Gr_',i,i,pnegf%Epnt)
 ENDDO

 call blk2csr(pnegf%str,pnegf%SM,Gl)

 DEALLOCATE(Glsub)
 
 !Distruzione dell'array Gr
 CALL destroy_Gr()
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
  SUBROUTINE destroy_Gr()
    integer :: i, i1, nbl

    nbl=size(Gr,1)

    DO i=1,nbl
       DO i1=1,nbl
          IF (ALLOCATED(Gr(i,i1)%val)) THEN
             CALL destroy(Gr(i,i1))
             !if (debug) WRITE(*,*) 'Deallocating Gr out of Make_Gl_mem'
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE destroy_Gr

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
  SUBROUTINE reload_r(pnegf, ESH, Sigma_ph_r, iter)

  	TYPE(Tnegf), pointer :: pnegf
  	TYPE(z_DNS), DIMENSION(:,:) :: ESH, Sigma_ph_r

  	INTEGER :: n, nbl, nrow

  	nbl = pnegf%str%num_PLs

  	if (pnegf%elph%diagonal) THEN

  		DO n = 1, nbl

  		   nrow = ESH(n,n)%nrow

	       call create(Sigma_ph_r(n,n), nrow, nrow)

  		   if (iter .eq. 0) then
  		      Sigma_ph_r(n,n)%val = (0.0_dp, 0.0_dp)
  		   else
	          call read_blkmat(Sigma_ph_r(n,n),pnegf%scratch_path,'Sigma_ph_r_',n,n,pnegf%Epnt)
 		   endif

	       ESH(n,n)%val = ESH(n,n)%val - Sigma_ph_r(n,n)%val

	       call destroy(Sigma_ph_r(n,n))

        END DO

  	ELSE

  	ENDIF


  END SUBROUTINE reload_r

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

          max=MAXVAL(ABS(Gr(i-1,n)%val(:,:)))

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

          max=MAXVAL(ABS(Gr(i+1,n)%val(:,:)))

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

  SUBROUTINE Make_GreenR_mem2_dns(P,nbl,indblk,A)

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
       call destroy(Gr(1,1))
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

  END SUBROUTINE Make_GreenR_mem2_dns





  !****************************************************************************
  !
  ! Calculate Gless contributions for all contacts (except reference)
  ! Writing on memory
  !
  !****************************************************************************

  SUBROUTINE Make_Gl_mem_dns(ESH,SelfEneR,frm,ref,struct,P,Gl)

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
    !Gl: sparse matrix containing G less contacts contribution
    !*******************************************************************************  

    IMPLICIT NONE

    !In/Out
    TYPE(z_CSR) :: Gl
    TYPE(z_DNS), DIMENSION(:,:) :: ESH
    TYPE(z_DNS), DIMENSION(:) :: SelfEneR
    TYPE(Tstruct_info), intent(in) :: struct
    REAL(dp), DIMENSION(:) :: frm
    INTEGER :: ref

    !Work
    TYPE(z_CSR) :: P, Gam1, Gl_sp
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

    !***
    !Iterazione sui contatti
    !***

    DO j=1,ncont

       IF (j.NE.ref .AND. ABS(frm(j)-frm(ref)).GT.drop) THEN

          !Calcolo della Gamma corrispondente
          cb=cblk(j)
          !nota:cb indica l'indice di blocco corrispondente al contatto j-esimo

          CALL zspectral(SelfEneR(j),SelfEneR(j),0,Gam)

          !conversione in densa
         ! !call create(Gam,Gam1%nrow,Gam1%ncol)

         ! !call csr2dns(Gam1,Gam)

         ! !call destroy(Gam1)

          !***
          !Iterazione sui blocchi 
          !***          
          ! [frm(j)-frm(ref)]
          frmdiff = cmplx(frm(j)-frm(ref),0.d0,dp)

          !Calcolo del sottoblocco Gl(1,1) fuori iterazione

          CALL prealloc_mult(Gr(1,cb),Gam,work1)    
 
          CALL zdagger(Gr(1,cb),Ga)
  
          !G(1,j) va tenuta per la prima iterazione blocco sopradiagonale

          CALL prealloc_mult(work1,Ga,frmdiff,Glsub(1,1))
          CALL destroy(work1)
          CALL destroy(Ga)

          !Calcolo sottoblocchi diagonali, sopradiagonali e sottodiagonali
          DO i=2,nbl

             !Calcolo blocchi sopradiagonali
             CALL prealloc_mult(Gr(i-1,cb),Gam,work1)
   
             CALL zdagger(Gr(i,cb),Ga)
 
             CALL prealloc_mult(work1,Ga,frmdiff,Glsub(i-1,i))

             CALL destroy(work1)
             !Teniamo Ga per il calcolo del blocco diagonale

             !Calcolo blocchi diagonali
             CALL prealloc_mult(Gr(i,cb),Gam,work1)

             CALL prealloc_mult(work1,Ga,frmdiff,Glsub(i,i))

             CALL destroy(work1)
             call destroy(Ga)

             CALL prealloc_mult(Gr(i,cb),Gam,work1)
             CALL zdagger(Gr(i-1,cb),Ga)
             CALL prealloc_mult(work1,Ga,frmdiff,Glsub(i,i-1))

             CALL destroy(work1)
             call destroy(Ga)

          ENDDO

          call destroy(Gam)

       ENDIF

    ENDDO

  END SUBROUTINE Make_Gl_mem_dns


  !****************************************************************************
  !
  ! Calculate Gless contributions due to el-ph
  ! Writing on memory
  !
  !****************************************************************************
  SUBROUTINE Make_Gl_ph(pnegf,ESH,Sigma_ph_less,iter)

    TYPE(Tnegf), pointer :: pnegf
    TYPE(z_DNS), DIMENSION(:,:) :: ESH, Sigma_ph_less
    INTEGER :: iter


    Type(z_DNS) :: Ga, work1
    INTEGER :: n, nbl, nrow

    nbl = pnegf%str%num_PLs

    if (pnegf%elph%diagonal) THEN

      DO n = 1, nbl

  	nrow = ESH(n,n)%nrow

	call create(Sigma_ph_less(n,n), nrow, nrow)

  	if (iter .eq. 0) then
  	   Sigma_ph_less(n,n)%val = (0.0_dp, 0.0_dp)
  	else
	   call read_blkmat(Sigma_ph_less(n,n),pnegf%scratch_path,'Sigma_ph_less_',n,n,pnegf%Epnt)
  	endif

        CALL zdagadns(Gr(n,n),Ga)

        CALL prealloc_mult(Gr(n,n), Sigma_ph_less(n,n), work1)

        CALL prealloc_mult(work1, Ga, Glsub(n,n))

        CALL destroy(work1,Ga)

      END DO
 
    ELSE

    ENDIF


  END SUBROUTINE Make_Gl_ph

  ! Computes Sigma_ph< and save it file
  !
  SUBROUTINE Sigma_ph_less(pnegf,Epnt)

    TYPE(Tnegf),POINTER :: pnegf
    REAL(dp), DIMENSION(:) :: Epnt

    !Local variables
    TYPE(z_DNS) :: work1,work2,work3
    TYPE(z_DNS) :: Mq_mat
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_less, G_less_interP, G_less_interN

    REAL(dp), DIMENSION(:), POINTER :: Wq
    REAL(dp), DIMENSION(:), POINTER :: Nq
    REAL(dp), DIMENSION(:), POINTER :: Mq
    LOGICAL, DIMENSION(:), POINTER :: selmodes

    INTEGER :: i, m, iE
    INTEGER :: nummodes, numselmodes, nbl
    INTEGER, DIMENSION(:), pointer :: indblk


    nbl = pnegf%str%num_PLs
    indblk => pnegf%str%mat_PL_start

    selmodes => pnegf%elph%selmodes
    Wq => pnegf%elph%Wq
    Mq => pnegf%elph%Mq
    Nq => pnegf%elph%Nq
    numselmodes = pnegf%elph%numselmodes
    nummodes = pnegf%elph%nummodes

    iE = pnegf%Epnt

    allocate(G_less_interP(nbl,nbl))
    allocate(G_less_interN(nbl,nbl))
    allocate(Sigma_less(nbl,nbl))

    do i = 1, nbl
      m = indblk(i+1)-indblk(i)
      call create(Sigma_less(i,i), m, m)
      Sigma_less(i,i)%val = (0.d0, 0.d0)
      call create(G_less_interP(i,i), m, m)
      call create(G_less_interN(i,i), m, m)
    enddo

    !PROCEDURA PER CALCOLARE LA Sigma_ph_less
    do m = 1 , nummodes

      if (.not.selmodes(m)) cycle
      !Ora faccio il loop sui punti energetici.
      !Epnt e' l'energia del punto i-esimo
      i1 = iE + 1
      do while (Epnt(i1) .lt. Epnt(iE)+Wq(m))
          i1 = i1 + 1
      end do

      call interpolation(i1-1, i1, Epnt(iE)+Wq(m), Epnt, pnegf%scratch_path, &
                         'G_less_', G_less_interP)

      !Ora faccio il loop sui punti energetici.
      !Epnt e' l'energia del punto i-esimo
      i1 = iE - 1
      do while (Epnt(i1) .gt. Epnt(iE)-Wq(m))
          i1 = i1 - 1
      end do

      call interpolation(i1, i1+1, Epnt(iE)-Wq(m), Epnt, pnegf%scratch_path, &
                         'G_less_', G_less_interN)

      do i = 1, nbl

        call create_id(Mq_mat,G_less_interN(i,i)%nrow,Mq(m))

        call prealloc_mult( Mq_mat, G_less_interN(i,i), work1)

        call prealloc_mult ( work1, Mq_mat, work2)

        call destroy(work1)

        call prealloc_mult( Mq_mat, G_less_interP(i,i), work1)

        call prealloc_mult ( work1, Mq_mat, work3)

        call destroy(work1)

        Sigma_less(i,i)%val = Sigma_less(i,i)%val + (Nq(m)+1)*work2%val + Nq(m)*work3%val

        call destroy(work2, work3, Mq_mat)

     end do

  end do

  do i = 1, nbl
    call write_blkmat(Sigma_less(i,i),pnegf%scratch_path,'Sigma_ph_less_',i,i,iE)
    call destroy(Sigma_less(i,i))
    call destroy(G_less_interN(i,i))
    call destroy(G_less_interP(i,i))
  enddo

  deallocate(Sigma_less,G_less_interN, G_less_interP)

END SUBROUTINE Sigma_ph_less


SUBROUTINE Sigma_ph_r(pnegf,Epnt)

  TYPE(Tnegf),POINTER :: pnegf
  REAL(dp), DIMENSION(:) :: Epnt

  !Local variables
  TYPE(z_DNS) :: work1,work2,work3
  TYPE(z_DNS) :: Mq_mat
  TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_r

  REAL(dp), DIMENSION(:), POINTER :: Nq
  REAL(dp), DIMENSION(:), POINTER :: Mq
  LOGICAL, DIMENSION(:), POINTER :: selmodes
  INTEGER :: i, m, iE
  INTEGER :: nummodes, numselmodes, nbl
  INTEGER, DIMENSION(:), pointer :: indblk


  nbl = pnegf%str%num_PLs
  indblk => pnegf%str%mat_PL_start


  selmodes => pnegf%elph%selmodes
  Mq => pnegf%elph%Mq
  Nq => pnegf%elph%Nq
  nummodes = pnegf%elph%nummodes

  allocate(Sigma_r(nbl,nbl))
  allocate(Gr(nbl,nbl))

  do i = 1, nbl
     m = indblk(i+1)-indblk(i)
      call create(Sigma_r(i,i), m, m)
      Sigma_r(i,i)%val = (0.d0, 0.d0)
      call create(Gr(i,i),m,m)
  enddo

  call read_blkmat(Gr(i,i),pnegf%scratch_path,'Gr_',i,i,pnegf%Epnt)

  do m=1 , nummodes

      if (.not.selmodes(m)) cycle

      call create_id(Mq_mat,Sigma_r(i,i)%nrow,Mq(m))

      call prealloc_mult( Mq_mat, Gr(i,i), work1)

      call prealloc_mult ( work1, Mq_mat, work2)

      call destroy(work1)

      Sigma_r(i,i)%val = Sigma_r(i,i)%val + (2*Nq(m)+1)*work2%val

      call destroy(work1, Mq_mat)

  enddo

  do i = 1, nbl
    call write_blkmat(Sigma_r(i,i),pnegf%scratch_path,'Sigma_ph_r_',i,i,pnegf%Epnt)
    call destroy(Sigma_r(i,i))
    call destroy(Gr(i,i))
  enddo

  deallocate(Sigma_r,Gr)

END SUBROUTINE Sigma_ph_r


SUBROUTINE interpolation(i1,i2, E, pnts, path, name, G_less_interp)
     INTEGER, intent(in) :: i1, i2
     REAL(dp) :: E
     REAL(dp), DIMENSION(:) :: pnts
     TYPE(z_DNS), DIMENSION(:,:) :: G_less_interp
     CHARACTER(*) :: path
     CHARACTER(*) :: name

     !local variables
     TYPE(z_DNS) :: work1,work2
     INTEGER :: i

     do i = 1, size(G_less_interp,1)

        call create(work1,G_less_interp(i,i)%nrow,G_less_interp(i,i)%ncol)

        call read_blkmat(work1, path, name, i, i, i1)

        call create(work2,G_less_interp(i,i)%nrow,G_less_interp(i,i)%ncol)

        call read_blkmat(work1, path, name, i, i, i2)

        G_less_interp(i,i)%val = ((E-pnts(i1))*work1%val + (pnts(i2)-E)*work2%val)/(pnts(i2)-pnts(i1))

     end do

END SUBROUTINE interpolation

  ! READ Matrices
  SUBROUTINE read_blkmat(Matrix, path, name, i, j, iE)

  	TYPE(z_DNS) :: Matrix
  	CHARACTER(*) :: path
  	CHARACTER(*) :: name
  	INTEGER :: i, j, iE


  	CHARACTER(64) :: filename
    character(4) :: ofblki, ofblkj
    character(10) :: ofpnt

    if (n.le.9999) write(ofblki,'(i4.4)') i
	if (n.gt.9999) stop 'ERROR: too many blks (> 9999)'
    if (n.le.9999) write(ofblkj,'(i4.4)') j
	if (n.gt.9999) stop 'ERROR: too many blks (> 9999)'

	if (iE.le.99999) write(ofpnt,'(i5.5)') iE

	filename = trim(name)//ofblki//'_'//ofblkj//'_'//ofpnt//'.dat'

	open(1001,file=trim(path)//filename, form='UNFORMATTED')
    DO
        read (1001,end = 100) i1,i2,mat_el
        Matrix%val(i1,i2)=mat_el
    ENDDO
    100 close(1001)

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

    if (n.le.9999) write(ofblki,'(i4.4)') i
	if (n.gt.9999) stop 'ERROR: too many blks (> 9999)'
    if (n.le.9999) write(ofblkj,'(i4.4)') j
	if (n.gt.9999) stop 'ERROR: too many blks (> 9999)'

	if (iE.le.99999) write(ofpnt,'(i5.5)') iE

	filename = trim(name)//ofblki//'_'//ofblkj//'_'//ofpnt//'.dat'

	open(1001,file=trim(path)//filename, form='UNFORMATTED')

	call outmat_c(1001,.false.,Matrix%val,Matrix%nrow,Matrix%ncol)

	close(1001)

  END SUBROUTINE write_blkmat

  !****************************************************************************
  !
  !  Calculate Gless contributions for all contacts but collector, in the 
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

  !Concatenation for every contact in Gless. Performs a sum on elements, not a replacement
  !Similar to Make_GreenR_mem2, except for sum of elements
  !Note: to backup old version zconcat calls (and Glsub deallocations) must be
  !      uncommented and all this part removed
  !If only one block is present, concatenation is not needed and it's implemented in a
  !more trivial way
  SUBROUTINE blk2csr(struct,P,Gl)

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

        call create(Gl_sp, Glsub(1,1)%nrow, Glsub(1,1)%ncol, nzdrop(Glsub(1,1),drop) )
        call dns2csr(Glsub(1,1),Gl_sp)
        call destroy(Glsub(1,1))
        call zmask_realloc(Gl_sp,P)
        call concat(Gl,Gl_sp,1,1)
        call destroy(Gl_sp)

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
                      WRITE(*,*) 'ERROR in Make_Gl_mem: probably wrong PL size', x
                      write(*,*) 'row',ii,Gl%colind(jj)
                      write(*,*) 'block indeces:',indblk(1:nbl)
                      STOP
                   ENDIF
                   
                   col = Gl%colind(jj) - indblk(y) + 1

                   Gl%nzval(jj) = Gl%nzval(jj) + Glsub(x,y)%val(row,col)

                ENDDO

                IF(oldx.NE.x) THEN

                   IF (x.eq.2)  THEN 
                      CALL destroy(Glsub(x-1,x-1))
                      CALL destroy(Glsub(x-1,x))

                   ELSEIF  (x.GT.2)  THEN 
                      CALL destroy(Glsub(x-1,x-1))
                      CALL destroy(Glsub(x-1,x-2))
                      CALL destroy(Glsub(x-1,x))
                   endif

                ENDIF


        ENDDO

        CALL destroy(Glsub(nbl, nbl))
        CALL destroy(Glsub(nbl, nbl - 1))

    ENDIF

  END SUBROUTINE blk2csr

  !****************************************************************************
  !
  !  Calculate Green Retarded in the 
  !  contacts regions, where overlap with device orbitals is non-zero
  !  (only the upper part, needed in both K=0 and K points calculations,
  !   or both upper and lower parts)
  !  writing on memory
  !
  !****************************************************************************

  SUBROUTINE Outer_GreenR_mem_dns(Tlc,Tcl,gsurfR,struct,lower,Aout)

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

  END SUBROUTINE Outer_GreenR_mem_dns



  !****************************************************************************
  !
  !  Calculate Gless contributions for all contacts except collector, in the
  !  outer region where contacts-device overlap is non-zero - writing on 
  !  memory
  !
  !****************************************************************************

  SUBROUTINE Outer_Gl_mem_dns(Tlc,gsurfR,SelfEneR,struct,frm,ref,lower,Glout)

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
    !Glout: sparse matrix containing G less contributions in the region 
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

  END SUBROUTINE Outer_Gl_mem_dns


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
          max=maxval(abs(Gr(i-1,nt1)%val(:,:)))

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
