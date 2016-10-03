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
#  undef __PAR3DISO
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
  USE elph
  USE ln_structure, only : TStruct_Info
  USE lib_param, only : MAXNCONT, Tnegf, intarray
  USE outmatrix, only : outmat_c, inmat_c, direct_out_c, direct_in_c 
  USE clock
  !USE transform

  IMPLICIT NONE
  private

  public :: tunneling_dns
  public :: tun_and_dos

  public :: calls_eq_mem_dns
  public :: calls_neq_mem_dns
  public :: calls_neq_elph

  public :: sigma_ph_n
  public :: sigma_ph_p
  public :: sigma_ph_r
  public :: sigma_ph_r_z
  public :: check_sigma_ph_r
  public :: check_Gl_Gr

  public :: rebuild_dns
  public :: Make_gsmr_mem_dns
  public :: Make_gsml_mem_dns
  public :: Make_Gr_mem_dns
  public :: Make_Grcol_mem_dns
  public :: Make_Gn_mem_dns
  public :: Outer_Gr_mem_dns
  public :: Outer_Gn_mem_dns
  !public :: Outer_A_mem_dns
  public :: iterative_meir_wingreen

  public :: complete_sigma_ph_r

  public :: create_scratch
  public :: destroy_scratch

  LOGICAL, PARAMETER :: debug=.false. 

  TYPE(z_DNS), DIMENSION(:), ALLOCATABLE :: gsmr
  TYPE(z_DNS), DIMENSION(:), ALLOCATABLE :: gsml
  TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Gr

  TYPE(z_DNS), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: GGn
  TYPE(z_DNS), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: GGp
  TYPE(z_DNS), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: GGr
  TYPE(z_DNS), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: Sgn
  TYPE(z_DNS), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: Sgp
  TYPE(z_DNS), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: Sgr
  LOGICAL, PARAMETER :: memory = .true.

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
    !pnegf:    negf data container
    !E:        Energy point
    !SelfEneR: matrices array containing contacts Self Energy
    !Tlc:      matrices array containing contacts-device interaction blocks (ES-H)
    !Tcl:      matrices array containing device-contacts interaction blocks (ES-H)
    !gsurfR:   matrices array containing contacts surface green
    !outer:    optional parameter (0,1,2).  
    !
    !Output:
    !A: Spectral function (Device + Contacts overlap regions -> effective conductor)
    !   outer = 0  no outer parts are computed
    !   outer = 1  only D/C part is computed
    !   outer = 2  D/C and C/D parts are computed
    !              (needed for K-points calculations)
    !
    !*****************************************************************************

    IMPLICIT NONE

    !In/Out
    TYPE(Tnegf), intent(inout) :: pnegf
    COMPLEX(dp), intent(in) :: E
    TYPE(z_DNS), DIMENSION(:), intent(in) :: SelfEneR
    TYPE(z_DNS), DIMENSION(:), intent(in) :: Tlc, Tcl, gsurfR
    TYPE(z_CSR), intent(out) :: A
    TYPE(Tstruct_info), intent(in) :: struct
    INTEGER, intent(in) :: outer

    !Work
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: ESH
    TYPE(z_CSR) :: ESH_tot, Ain
    INTEGER :: i,ierr, nbl, ncont
    INTEGER, DIMENSION(:), POINTER :: cblk, indblk

    nbl = struct%num_PLs
    ncont = struct%num_conts
    cblk => struct%cblk
    indblk => struct%mat_PL_start

    ! Take CSR H,S and build ES-H in dense blocks
    CALL prealloc_sum(pnegf%H,pnegf%S,(-1.d0, 0.d0),E,ESH_tot)

    call allocate_blk_dns(ESH,nbl)

    CALL zcsr2blk_sod(ESH_tot,ESH,indblk)

    CALL destroy(ESH_tot)

    DO i=1,ncont
      ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
    ENDDO

    !! Add interaction self energy contribution, if any
    if (allocated(pnegf%inter)) call pnegf%inter%add_sigma_r(ESH)
    
    !----------------------------------

    call allocate_gsm_dns(gsmr,nbl)
    CALL Make_gsmr_mem_dns(ESH,nbl,2)

    call allocate_blk_dns(Gr,nbl)

    CALL Make_Gr_mem_dns(ESH,1)
    CALL Make_Gr_mem_dns(ESH,2,nbl)

    CALL destroy_ESH(ESH)
    CALL deallocate_blk_dns(ESH)

    call destroy_gsm(gsmr)
    call deallocate_gsm_dns(gsmr)


    ! SAVE ON FILES/MEMORY (for elph).........................
    if (pnegf%elph%numselmodes.gt.0 .and. pnegf%elph%model .eq. -1) then
      ! save diagonal blocks of Gn = -i G<
      DO i = 1, nbl
        !print*,'G_r ',minval(abs(Gr(i,i)%val)), maxval(abs(Gr(i,i)%val))
        call write_blkmat(Gr(i,i),pnegf%scratch_path,'G_r_',i,i,pnegf%iE)
      ENDDO
    endif
    !..........................................................
    !! Deliver Gr to interaction models if any
    if (allocated(pnegf%inter)) call pnegf%inter%set_Gr(Gr, pnegf%iE)
    !-----------------------------------------------------------

    call blk2csr(Gr,struct,pnegf%S,A)

    SELECT CASE (outer)
    CASE(0)
    CASE(1)
      CALL Outer_Gr_mem_dns(Tlc,Tcl,gsurfR,struct,.FALSE.,A)   
    CASE(2)
      CALL Outer_Gr_mem_dns(Tlc,Tcl,gsurfR,struct,.TRUE.,A) 
    END SELECT

    !Distruzione dell'array Gr
    CALL destroy_blk(Gr)
    DEALLOCATE(Gr)

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
    !pnegf:    negf data container
    !E:        Energy point
    !SelfEneR: matrices array containing contacts Self Energy
    !Tlc:      matrices array containing contacts-device interaction blocks (ES-H)
    !Tcl:      matrices array containing device-contacts interaction blocks (ES-H)
    !gsurfR:   matrices array containing contacts surface green
    !out:      optional parameter (0,1,2).  
    !
    !Output:
    !Gn: NE GF (Device + Contacts overlap regions -> effective conductor)
    !   out = 0  no outer parts are computed
    !   out = 1  only D/C part is computed
    !   out = 2  D/C and C/D parts are computed
    !              (needed for K-points calculations)
    !
    !*****************************************************************************

    IMPLICIT NONE

    !In/Out
    TYPE(Tnegf), intent(inout) :: pnegf
    TYPE(z_CSR), intent(inout)  :: Glout
    TYPE(z_DNS), DIMENSION(:), intent(in)  :: SelfEneR, gsurfR, Tlc, Tcl
    REAL(dp), intent(in)  :: E
    TYPE(Tstruct_info), intent(in)  :: struct
    REAL(dp), DIMENSION(:), intent(in)  :: frm
    INTEGER, intent(in)  :: out

    !Work
    INTEGER :: ref
    COMPLEX(dp) :: Ec
    INTEGER :: i,ierr,ncont,nbl, lbl, rbl
    INTEGER, DIMENSION(:), POINTER :: cblk, indblk
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: ESH
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Gn
    TYPE(z_CSR) :: ESH_tot, Gl
    LOGICAL :: mask(MAXNCONT)

    if (pnegf%elph%model .ne. 0) then
      !TODO: output logfile
      write(*,*) "Warning: calls_neq_mem is not compatible with el-ph models"
    endif

    nbl = struct%num_PLs
    ncont = struct%num_conts
    indblk => struct%mat_PL_start
    cblk => struct%cblk
    ref = pnegf%refcont

    Ec=cmplx(E,0.d0,dp)

    ! Take CSR H,S and build ES-H in dense blocks
    CALL prealloc_sum(pnegf%H,pnegf%S,(-1.d0, 0.d0),Ec,ESH_tot)

    call allocate_blk_dns(ESH,nbl)

    CALL zcsr2blk_sod(ESH_tot,ESH,indblk)

    CALL destroy(ESH_tot)

    DO i=1,ncont
      ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
    ENDDO

    call allocate_gsm_dns(gsmr,nbl)
    call allocate_gsm_dns(gsml,nbl)

    ! Compute blocks for gsmr and gsml
    mask = .true.
    mask(ref) = .false. 
    rbl = minval(cblk(1:ncont),mask(1:ncont)) + 1  
    lbl = maxval(cblk(1:ncont),mask(1:ncont)) - 1

    ! Fix to a bug when there are 2PLs
    ! later Make_Gr tries to compute Gr(1,1) but needs gsmr(2,2)
    ! 
    IF (nbl.eq.2) then
      CALL Make_gsmr_mem_dns(ESH,nbl,rbl-1)
      CALL Make_gsml_mem_dns(ESH,1,lbl+1)    
    ELSE
      CALL Make_gsmr_mem_dns(ESH,nbl,rbl)
      CALL Make_gsml_mem_dns(ESH,1,lbl)    
    ENDIF

    call allocate_blk_dns(Gr,nbl)

    ! -------------------------------------------------------------
    ! 1. rbl>lbl  => lbl+1=rbl-1 => compute first Gr(rbl-1,rbl-1)
    ! 2. rbl<lbl  => lbl=rbl-2 has been computed
    ! Make_Gr does not compute if sbl>nbl or sbl<1
    CALL Make_Gr_mem_dns(ESH,rbl-1)
    CALL Make_Gr_mem_dns(ESH,rbl,nbl)
    CALL Make_Gr_mem_dns(ESH,rbl-2,1)

    !Chiamata di Make_Grcol_mem per i contatti necessari 
    DO i=1,ncont
      IF (i.NE.ref) THEN
        CALL Make_Grcol_mem_dns(ESH,cblk(i),indblk)
      ENDIF
    ENDDO

    !Distruzione delle gsmall................................
    call destroy_gsm(gsmr)
    call deallocate_gsm_dns(gsmr)
    call destroy_gsm(gsml)
    call deallocate_gsm_dns(gsml)
    !.......................................................

    !Calcolo della G_n nel device
    call allocate_blk_dns(Gn,nbl)

    call init_blkmat(Gn,ESH)

    CALL Make_Gn_mem_dns(ESH,SelfEneR,frm,ref,struct,Gn)

    call blk2csr(Gn,struct,pnegf%S,Glout)

    !Calcolo degli outer blocks 
    SELECT CASE (out)
    CASE(0)
    CASE(1)
      CALL Outer_Gn_mem_dns(Tlc,gsurfR,SelfEneR,struct,frm,ref,.false.,Glout)
    CASE(2)
      CALL Outer_Gn_mem_dns(Tlc,gsurfR,SelfEneR,struct,frm,ref,.true.,Glout)
    END SELECT

    CALL destroy_blk(Gn)
    DEALLOCATE(Gn)

    CALL destroy_blk(Gr)
    DEALLOCATE(Gr)

    CALL destroy_ESH(ESH)
    DEALLOCATE(ESH)

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

  SUBROUTINE calls_neq_elph(pnegf,E,SelfEneR,Tlc,Tcl,gsurfR,frm,Glout,out)

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
    TYPE(Tnegf), intent(inout) :: pnegf
    TYPE(z_CSR) :: Glout
    TYPE(z_DNS), DIMENSION(:) :: SelfEneR, gsurfR, Tlc, Tcl
    REAL(dp) :: E
    REAL(dp), DIMENSION(:) :: frm
    INTEGER :: out

    !Work
    COMPLEX(dp) :: Ec
    INTEGER :: i,ierr,ncont,nbl, lbl
    INTEGER :: ref, iter
    INTEGER, DIMENSION(:), POINTER :: cblk, indblk
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: ESH, Gn, Gp
    TYPE(z_CSR) :: ESH_tot, Gl
    LOGICAL :: mask(MAXNCONT)
    TYPE(Tstruct_info) :: struct
    REAL(dp), DIMENSION(:), allocatable :: cfrm

    struct = pnegf%str
    nbl = struct%num_PLs
    ncont = struct%num_conts
    indblk => struct%mat_PL_start
    cblk => struct%cblk
    ref = pnegf%refcont
    iter = pnegf%inter%scba_iter

    Ec=cmplx(E,0.d0,dp)

    !Costruiamo la matrice sparsa ESH
    CALL prealloc_sum(pnegf%H,pnegf%S,(-1.d0, 0.d0),Ec,ESH_tot)

    call allocate_blk_dns(ESH,nbl)

    CALL zcsr2blk_sod(ESH_tot,ESH,indblk)

    call destroy(ESH_tot)
    DO i=1,ncont
      ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
    ENDDO

    !! Add interaction self energy if any
    if (allocated(pnegf%inter)) call pnegf%inter%add_sigma_r(ESH)

    !---------------------------------------------
    !Allocazione delle gsmr
    call allocate_gsm_dns(gsmr,nbl)
    call allocate_gsm_dns(gsml,nbl)
    !Chiamata di Make_gsmr_mem
    CALL Make_gsmr_mem_dns(ESH,nbl,2)

    !Chiamata di Make_gsml_mem solo per i blocchi 1..lbl dove  
    ! lbl = maxval(cblk,mask) - 2 
    lbl = nbl - 2    ! ALL BLOCKS are needed!!

    if( ncont.gt.1 ) then
      CALL Make_gsml_mem_dns(ESH,1,lbl)    
    endif
    !if (debug) write(*,*) '----------------------------------------------------'
    !if (debug) call writeMemInfo(6)
    !if (debug) call writePeakInfo(6)
    !if (debug) write(*,*) '----------------------------------------------------'

    call allocate_blk_dns(Gr,nbl)

    CALL Make_Gr_mem_dns(ESH,1)
    CALL Make_Gr_mem_dns(ESH,2,nbl)
    !! Update el-ph retarded self energy if any
    if (allocated(pnegf%inter)) call pnegf%inter%set_Gr(Gr, pnegf%iE)
    !--------------------------------------------------------
    !With el-ph we need all columns
    DO i=1,nbl
      CALL Make_Grcol_mem_dns(ESH,i,indblk)
    ENDDO

    !Distruzione delle gsmall
    call destroy_gsm(gsmr)
    call deallocate_gsm_dns(gsmr)
    call destroy_gsm(gsml)
    call deallocate_gsm_dns(gsml)

    !Calcolo degli outer blocks 
    !if (debug) write(*,*) 'Compute Outer_G_n' 
    SELECT CASE (out)
    CASE(0)
    CASE(1)
      CALL Outer_Gn_mem_dns(Tlc,gsurfR,SelfEneR,pnegf%str,frm,ref,.false.,Glout)
    CASE(2)
      CALL Outer_Gn_mem_dns(Tlc,gsurfR,SelfEneR,pnegf%str,frm,ref,.true.,Glout)
    END SELECT

    !Calcolo della G_n nel device
    !if (debug) write(*,*) 'Compute G_n' 
    !Allocazione degli array di sparse
    ALLOCATE(Gn(nbl,nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Gn(nbl,nbl)'

    call init_blkmat(Gn,ESH)

    CALL Make_Gn_mem_dns(ESH,SelfEneR,frm,ref,struct,Gn)
    call Make_Gn_ph(pnegf,ESH,iter,Gn)

    !! Pass Gr to interaction model
    if (allocated(pnegf%inter)) call pnegf%inter%set_Gn(Gn, pnegf%iE)
    !-----------------------------------------------------
    !! Skip this, old implementation
    !print*
    ! save diagonal blocks of G_n
    ! DO i = 1, nbl
    !   !print*,'(G_n) G_n',minval(abs(Gn(i,i)%val)), maxval(abs(Gn(i,i)%val))
    !   call write_blkmat(Gn(i,i),pnegf%scratch_path,'G_n_',i,i,pnegf%iE)
    !   !print*,'(G_r) G_r',minval(abs(Gn(i,i)%val)),maxval(abs(Gr(i,i)%val))
    !   call write_blkmat(Gr(i,i),pnegf%scratch_path,'G_r_',i,i,pnegf%iE)
    ! ENDDO

    ! the following destroys Gn
    call blk2csr(Gn,pnegf%str,pnegf%S,Gl)
    DEALLOCATE(Gn)

    !! Computation of G_p for Sigma_p
!!$ if (pnegf%elph%check) then
!!$    allocate(cfrm(10))
!!$    cfrm = 1.0_dp - frm
!!$    cfrm(ref) = 0.0_dp
!!$ 
!!$    ALLOCATE(Gp(nbl,nbl),stat=ierr)
!!$ 
!!$    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Gp(nbl,nbl)'
!!$    call init_blkmat(Gp,ESH)
!!$    
!!$    CALL Make_Gn_mem_dns(ESH,SelfEneR,cfrm,ref,struct,Gp)
!!$    deallocate(cfrm)
!!$    
!!$    call Make_Gp_ph(pnegf,ESH,iter,Gp)
!!$    
!!$   !  save diagonal blocks of G_p
!!$   DO i = 1, nbl
!!$     !print*
!!$     !print*,'(G_p) G_p',minval(abs(Gn(i,i)%val)),maxval(abs(Gn(i,i)%val))
!!$     call write_blkmat(Gp(i,i),pnegf%scratch_path,'G_p_',i,i,pnegf%iE)
!!$   ENDDO
!!$    
!!$   call destroy_blk(Gp)
!!$   DEALLOCATE(Gp)
!!$ endif
    !! ------------------------------------------------------------------------

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

  END SUBROUTINE calls_neq_elph

  !---------------------------------------------------------------------
  !>
  !  Iterative algorithm implementing Meir Wingreen formula for a given
  !  electrode
  !  Note: self consistent born approximation is not accounted for here
  !  It is assumed that the el-ph container already includes the 
  !  desired values. SCBA loop should be run outside
  !  Many operations from calls_neq_ph are repeated here, as it is
  !  assumed that A and Gn are not available at the time of the call
  !
  !  It implements the form without the reference electrode
  !  Iop = \Sigma_{i}^{n}A - \Gamma_{i}G^{n} =
  !      = \Gamma_{i}[f_{i}A - f_{ref}A - G^{n,l\neq ref} +
  !         - G^{n,l\neq ref}_{\phi}]
  !  where G^{n,l\neq ref} is the component including no el-ph
  !  If i=ref it reduces to
  !  Iop = \Gamma_{i}[-G^{n,l\neq ref}  - G^{n,l\neq ref}_{\phi}]
  !---------------------------------------------------------------------
  subroutine iterative_meir_wingreen(pnegf,E,SelfEneR,Tlc,Tcl,gsurfR,frm, ni,tun_mat)
    IMPLICIT NONE

    !In/Out
    TYPE(Tnegf), intent(inout) :: pnegf
    integer, dimension(:) :: ni
    TYPE(z_DNS), DIMENSION(:) :: SelfEneR, gsurfR, Tlc, Tcl
    Real(dp), Dimension(:), allocatable :: tun_mat
    REAL(dp) :: E
    REAL(dp), DIMENSION(:) :: frm
    INTEGER :: out, size_ni
    complex(dp) :: tmp

    !Work
    COMPLEX(dp) :: Ec
    INTEGER :: i,ierr,ncont,nbl,lbl
    INTEGER :: ref, iter, lead, lead_blk, ref_blk
    INTEGER, DIMENSION(:), POINTER :: cblk, indblk
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: ESH, Gn, Gp
    type(z_DNS) :: work1, work2, work3, Gamma, A
    TYPE(z_CSR) :: ESH_tot, Gl
    LOGICAL :: mask(MAXNCONT)
    TYPE(Tstruct_info) :: struct
    REAL(dp), DIMENSION(:), allocatable :: cfrm

    struct = pnegf%str
    nbl = struct%num_PLs
    ncont = struct%num_conts
    indblk => struct%mat_PL_start
    cblk => struct%cblk
    ref = pnegf%refcont
    ref_blk = pnegf%str%cblk(ref)
    iter = pnegf%inter%scba_iter

    Ec=cmplx(E,0.d0,dp)
    !Costruiamo la matrice sparsa ESH
    CALL prealloc_sum(pnegf%H,pnegf%S,(-1.d0, 0.d0),Ec,ESH_tot)

    call allocate_blk_dns(ESH,nbl)

    CALL zcsr2blk_sod(ESH_tot,ESH,indblk)

    call destroy(ESH_tot)

    DO i=1,ncont
      ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
    ENDDO

    !! Add el-ph self energy if any
    if (allocated(pnegf%inter)) call pnegf%inter%add_sigma_r(ESH)
    !------------------------------------------------
    !Allocazione delle gsmr
    call allocate_gsm_dns(gsmr,nbl)
    call allocate_gsm_dns(gsml,nbl)

    !Chiamata di Make_gsmr_mem
    CALL Make_gsmr_mem_dns(ESH,nbl,2)

    !Chiamata di Make_gsml_mem solo per i blocchi 1..lbl dove  
    ! lbl = maxval(cblk,mask) - 2 
    lbl = nbl - 2    ! ALL BLOCKS are needed!!

    if( ncont.gt.1 ) then
      CALL Make_gsml_mem_dns(ESH,1,lbl)    
    endif
    !if (debug) write(*,*) '----------------------------------------------------'
    !if (debug) call writeMemInfo(6)
    !if (debug) call writePeakInfo(6)
    !if (debug) write(*,*) '----------------------------------------------------'

    call allocate_blk_dns(Gr,nbl)

    CALL Make_Gr_mem_dns(ESH,1)
    CALL Make_Gr_mem_dns(ESH,2,nbl)

    !! Give Gr to interaction model if any
    if (allocated(pnegf%inter)) call pnegf%inter%set_Gr(Gr, pnegf%iE)
    !---------------------------------------------------
    !With el-ph we need all columns
    DO i=1,nbl
      CALL Make_Grcol_mem_dns(ESH,i,indblk)
    ENDDO

    !Distruzione delle gsmall
    call destroy_gsm(gsmr)
    call deallocate_gsm_dns(gsmr)
    call destroy_gsm(gsml)
    call deallocate_gsm_dns(gsml)

    !! Never calculate outer blocks

    !Calcolo della G_n nel device
    !if (debug) write(*,*) 'Compute G_n' 
    !Allocazione degli array di sparse
    ALLOCATE(Gn(nbl,nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Gn(nbl,nbl)'

    call init_blkmat(Gn,ESH)

    !! TEMPORARY AND INEFFICIENT: CALCULATE THE FULL GN WHEN A BLOCK 
    !! (THE CONTACT ONE) IS ENOUGH
    !! WE HAVE Gn WITHOUT REFERENCE CONTACT?? (usual neq contributions)
    CALL Make_Gn_mem_dns(ESH,SelfEneR,frm,ref,struct,Gn)

    call Make_Gn_ph(pnegf,ESH,iter,Gn)

    ! I probably need the next set_Gn in interaction for current conservation here

    do i=1,size(ni)
      if (ni(i) .eq. 0) then
        cycle
      endif
      lead = ni(i)
      lead_blk = pnegf%str%cblk(lead)
      call zspectral(SelfEneR(lead),SelfEneR(lead), 0, Gamma)
      if (lead.eq.ref) then
        call prealloc_mult(Gamma, Gn(lead_blk, lead_blk), (-1.d0, 0.d0), work1)
        tun_mat(i) = real(trace(work1))
        call destroy(work1)
      else
        call zspectral(Gr(ref, ref), Gr(ref, ref), 0, A)
        tmp = frm(lead)-frm(ref)
        call prealloc_sum(A, Gn(lead_blk, lead_blk), tmp, (-1.d0, 0.d0), work1)
        call destroy(A)
        call prealloc_mult(Gamma, work1, work2)
        call destroy(work1)
        tun_mat(i) = real(trace(work2))
        call destroy(work2)
      endif
      call destroy(Gamma)
    enddo

    !Convert to output CSR format.
    call blk2csr(Gn,pnegf%str,pnegf%S,Gl)
    DEALLOCATE(Gn)

    !Distruzione dell'array Gr
    CALL destroy_blk(Gr)
    DEALLOCATE(Gr)

    CALL destroy_ESH(ESH)
    DEALLOCATE(ESH)

  end subroutine iterative_meir_wingreen

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
      call create(Matrix(nbl,nbl),S(nbl,nbl)%nrow,S(nbl,nbl)%ncol)      
      Matrix(nbl,nbl)%val=(0.d0,0.d0)
    ENDIF

  END SUBROUTINE init_blkmat


  !**********************************************************************
  SUBROUTINE destroy_gsm(gsm)
    type(z_DNS), DIMENSION(:) :: gsm
    integer :: i, i1, nbl

    nbl=size(gsm,1)

    do i=1,nbl
      if (allocated(gsm(i)%val)) call destroy(gsm(i))
    end do

  END SUBROUTINE destroy_gsm

  !**********************************************************************
  SUBROUTINE destroy_blk(M)
    type(z_DNS), DIMENSION(:,:) :: M
    integer :: i, i1, nbl

    nbl=size(M,1)

    DO i=1,nbl
      DO i1=1,nbl
        IF (ALLOCATED(M(i1,i)%val)) THEN
          !print*,'kill Gr',i1,i
          CALL destroy(M(i1,i))
        END IF
      END DO
    END DO

  END SUBROUTINE destroy_blk

  !**********************************************************************
  SUBROUTINE destroy_ESH(ESH)

    integer :: i, nbl
    type(z_DNS), dimension(:,:) :: ESH

    nbl=size(ESH,1)

    DO i=1,nbl
      CALL destroy(ESH(i,i))
    END DO
    DO i=2,nbl
      CALL destroy(ESH(i-1,i))
      CALL destroy(ESH(i,i-1))
    END DO

  END SUBROUTINE destroy_ESH


  !***********************************************************************
  !> Update the value of el-ph Retarded  from block matrix
  !  Self Energy for the dephasing model
  !
  !***********************************************************************
  subroutine update_elph_r(negf, Gr)
    type(Tnegf), intent(inout) :: negf
    TYPE(z_DNS), dimension(:,:), intent(in) :: Gr

    TYPE(z_DNS) :: work1, work2, work3
    integer :: npl, nblk, n, ii, jj, norbs, indstart, indend, nnz

    npl = negf%str%num_PLs

    select case(negf%elph%model)
    case (0)
      return
    case (1)
      stop "Deprecated"
    case(2)
      stop "Deprecated"
    case(3)
      stop "Deprecated"
    case default
      write(*,*) 'Elph model not yet implemented'
      stop 0
    end select

  end subroutine update_elph_r

  !***********************************************************************
  !> Update the value of el-ph inscattering from block matrix
  !  Self Energy for the dephasing model
  !
  !***********************************************************************
  subroutine update_elph_n(negf, Gn)
    type(Tnegf), intent(inout) :: negf
    TYPE(z_DNS), dimension(:,:), intent(in) :: Gn
    
    type(z_DNS) :: work1, work2, work3
    integer :: npl, nblk, n, ii, jj, norbs, indstart, indend, nnz
    
    npl = negf%str%num_PLs

    select case(negf%elph%model)
    case (0)
      return
    case (1)
      stop "Deprecated"
    case(2)
      stop "Deprecated"
    case(3)
      stop "Deprecated"
    case default
      write(*,*) 'Elph model not yet implemented'
      stop 0
    end select
   
  end subroutine update_elph_n

  !***********************************************************************
  !
  !  Reloads the retarded elph self-energy and add it to ES-H
  !
  !***********************************************************************
  subroutine add_elph_sigma_r(pnegf, ESH, elph)
    type(Tnegf), intent(in) :: pnegf
    type(Telph), intent(in) :: elph
    type(z_DNS), dimension(:,:), intent(inout) :: ESH
    
    type(z_DNS), dimension(:,:), allocatable :: Sigma_ph_r, sigma_blk
    integer :: n, nbl, nrow, ierr, ii, jj, nblk, norbs, indstart, indend

    nbl = pnegf%str%num_PLs

    ! At first loop we should add no
    if (elph%scba_iter .eq. 0 .or. elph%model .eq. 0) then
      return
    end if
    
    select case(elph%model)
      
    case(0)
      return
    case(1)
      stop "Deprecated"
    case (2)
      stop "Deprecated"
    case(3)
      stop "Deprecated"
    ! Old Alex one
    case(-1)
      nbl = pnegf%str%num_PLs
      allocate(Sigma_ph_r(nbl,nbl),stat=ierr)
      IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_r'
      
      do n = 1, nbl
        nrow = ESH(n,n)%nrow
        call create(Sigma_ph_r(n,n), nrow, nrow)
        Sigma_ph_r(n,n)%val = (0.0_dp, 0.0_dp)
        if (elph%scba_iter .gt. 0) then
          call read_blkmat(Sigma_ph_r(n,n),pnegf%scratch_path,'Sigma_ph_r_',n,n,pnegf%iE)
        else
          call write_blkmat(Sigma_ph_r(n,n),pnegf%scratch_path,'Sigma_ph_r_',n,n,pnegf%iE)
        end if
        ESH(n,n)%val = ESH(n,n)%val - Sigma_ph_r(n,n)%val
        call destroy(Sigma_ph_r(n,n))
      end do
      deallocate(Sigma_ph_r)
    case default
      write(*,*) 'Elph model not yet implemented'
      stop 0
    end select
    
  end subroutine add_elph_sigma_r

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
    TYPE(z_DNS), DIMENSION(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl,ebl                 ! start block, end block

    !Work
    !TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: INV
    TYPE(z_DNS) :: work1, work2
    INTEGER :: nrow, M, N
    INTEGER :: i, nbl

    if (sbl.lt.ebl) return

    nbl = size(ESH,1)

    if (nbl.gt.1) then

      nrow=ESH(sbl,sbl)%nrow

      call create(gsmr(sbl),nrow,nrow)

      call compGreen(gsmr(sbl),ESH(sbl,sbl),nrow)

    endif

    DO i=sbl-1,ebl,-1

      CALL prealloc_mult(ESH(i,i+1),gsmr(i+1),(-1.d0, 0.d0),work1)

      CALL prealloc_mult(work1,ESH(i+1,i),work2)

      CALL destroy(work1)

      CALL prealloc_sum(ESH(i,i),work2,work1)

      CALL destroy(work2)

      CALL create(gsmr(i),work1%nrow,work1%nrow)

      CALL compGreen(gsmr(i),work1,work1%nrow)

      CALL destroy(work1)

    ENDDO


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
    TYPE(z_DNS), DIMENSION(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl,ebl                       ! start block, end block

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
    nrow=ESH(sbl,sbl)%nrow  

    CALL create(gsml(sbl),nrow,nrow)

    call compGreen(gsml(sbl),ESH(sbl,sbl),nrow)


    DO i=sbl+1,ebl

      nrow=ESH(i,i)%nrow   !indblk(i+1)-indblk(i)

      CALL prealloc_mult(ESH(i,i-1),gsml(i-1),(-1.d0, 0.d0),work1)

      CALL prealloc_mult(work1,ESH(i-1,i),work2)

      CALL destroy(work1)

      CALL prealloc_sum(ESH(i,i),work2,work1)

      CALL destroy(work2)

      call create(gsml(i),work1%nrow,work1%nrow)

      call compGreen(gsml(i),work1,work1%nrow)

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

  SUBROUTINE Make_Gr_mem_dns(ESH,sbl,ebl)

    !***********************************************************************
    !Input:
    !ESH: dense matrices array ESH(nbl,nbl)
    !sbl, ebl : block indexes
    ! If only sbl is specified, it calculates Gr(sbl, sbl)
    ! If sbl > ebl, it calculates Gr(ebl:sbl, ebl:sbl), Gr(ebl:sbl + 1, ebl:sbl),
    !    Gr(ebl:sbl, ebl:sbl + 1) (need gsml)          
    ! If sbl < ebl, it calculates Gr(sbl:ebl, sbl:ebl), Gr(sbl:ebl - 1, sbl:ebl),
    !    Gr(sbl:ebl, sbl:ebl - 1) (need gsmr)          
    !
    ! 
    !Output:
    !sparse matrices array global variable Gr(nbl,nbl) is available in 
    !memory - single blocks are allocated internally, array Gr(nbl,nbl) 
    !must be allocated externally
    !***********************************************************************

    IMPLICIT NONE 

    !In/Out
    TYPE(z_DNS), DIMENSION(:,:) :: ESH
    INTEGER :: sbl
    INTEGER, optional :: ebl

    !Work
    INTEGER :: i,nrow,nbl
    TYPE(z_DNS) :: work1, work2, work3

    nbl = size(ESH,1)

    if (sbl.gt.nbl) return 
    if (sbl.lt.1) return 

    if (.not.present(ebl)) then
      if (nbl.eq.1) then
        nrow = ESH(sbl,sbl)%nrow     
        CALL create(Gr(sbl,sbl),nrow,nrow)
        CALL compGreen(Gr(sbl,sbl),ESH(sbl,sbl),nrow)
      else
        nrow = ESH(sbl,sbl)%nrow     
        call create(work1,nrow,nrow)
        work1%val = ESH(sbl,sbl)%val
        if (sbl+1.le.nbl) then
          CALL prealloc_mult(ESH(sbl,sbl+1),gsmr(sbl+1),work2)
          CALL prealloc_mult(work2,ESH(sbl+1,sbl),work3)
          CALL destroy(work2)
          CALL prealloc_sum(work1,work3,(-1.d0, 0.d0),work2)
          CALL destroy(work3)
          work1%val = work2%val
          CALL destroy(work2)
        endif
        if (sbl-1.ge.1) then
          CALL prealloc_mult(ESH(sbl,sbl-1),gsml(sbl-1),work2)
          CALL prealloc_mult(work2,ESH(sbl-1,sbl),work3)
          CALL destroy(work2)
          CALL prealloc_sum(work1,work3,(-1.d0, 0.d0),work2)
          CALL destroy(work3)
          work1%val = work2%val
          CALL destroy(work2)
        endif

        CALL create(Gr(sbl,sbl),nrow,nrow)
        CALL compGreen(Gr(sbl,sbl),work1,nrow)
        CALL destroy(work1)
      endif
      return
    endif


    !***
    !Diagonal, Subdiagonal and Superdiagonal blocks
    !***
    IF ((ebl.ge.sbl).and.(ebl.gt.1).and.(sbl.gt.1)) THEN
      DO i=sbl,ebl,1
        CALL prealloc_mult(gsmr(i),ESH(i,i-1),work1)
        CALL prealloc_mult(work1,Gr(i-1,i-1),(-1.d0,0.d0),Gr(i,i-1))
        CALL destroy(work1)

        CALL prealloc_mult(ESH(i-1,i),gsmr(i),work2)
        CALL prealloc_mult(Gr(i-1,i-1),work2,(-1.d0, 0.d0),Gr(i-1,i))

        CALL prealloc_mult(Gr(i,i-1),work2,(-1.d0,0.d0),work1)
        CALL destroy(work2)

        CALL prealloc_sum(gsmr(i),work1,Gr(i,i))
        CALL destroy(work1) 
      ENDDO
    ELSE
      DO i=sbl,ebl,-1
        CALL prealloc_mult(gsml(i),ESH(i,i+1),work1)
        CALL prealloc_mult(work1,Gr(i+1,i+1),(-1.d0,0.d0),Gr(i,i+1))
        CALL destroy(work1)

        CALL prealloc_mult(ESH(i+1,i),gsml(i),work2)
        CALL prealloc_mult(Gr(i+1,i+1),work2,(-1.d0, 0.d0),Gr(i+1,i))

        CALL prealloc_mult(Gr(i,i+1),work2,(-1.d0,0.d0),work1)
        CALL destroy(work2)

        CALL prealloc_sum(gsml(i),work1,Gr(i,i))
        CALL destroy(work1) 
      ENDDO
    ENDIF

  END SUBROUTINE Make_Gr_mem_dns

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
    TYPE(z_DNS), DIMENSION(:,:), intent(in) :: ESH
    INTEGER, intent(in) :: n
    INTEGER, DIMENSION(:), intent(in) :: indblk 

    !Work
    INTEGER :: i,nrow,ncol,nbl
    TYPE(z_DNS) :: work1
    REAL(dp) :: max

    nbl = size(ESH,1)

    IF (n.GT.nbl) THEN
      STOP 'Error in Make_Grcol_mem : n is greater than nbl'
    ENDIF

    !***************************************
    !  Downgoing (j>=n+2 && n<nbl-1)
    !
    !   G_j,n = -gR_jj T_j,j-1 G_j-1,n
    !
    !***************************************
    ncol=indblk(n+1)-indblk(n)

    IF (n.LT.(nbl-1)) THEN

      DO i=n+2,nbl

        nrow=indblk(i+1)-indblk(i)

        max=MAXVAL(ABS(Gr(i-1,n)%val))

        IF (max.GT.EPS) THEN

          CALL prealloc_mult(gsmr(i),ESH(i,i-1),(-1.d0, 0.d0),work1)
          CALL prealloc_mult(work1,Gr(i-1,n),Gr(i,n))
          CALL destroy(work1)
        ENDIF

        ! NOW WHEN BLOCK IS SMALLER THAN EPS IT IS NOT CREATED
        !    CALL create(Gr(i,n),nrow,ncol)
        !    Gr(i,n)%val(:,:)=(0.d0,0.d0)
        ! ENDIF

      ENDDO

    ENDIF
    !*************************************
    !   Up-going (j<=n-2 && n>2)
    !
    !   G_j,n = -gL_jj T_j,j+1 G_j+1,n
    !
    !*************************************

    IF (n.GT.2) THEN

      DO i=n-2,1,(-1)

        nrow=indblk(i+1)-indblk(i)

        max=MAXVAL(ABS(Gr(i+1,n)%val))

        IF (max.GT.EPS) THEN

          CALL prealloc_mult(gsml(i),ESH(i,i+1),(-1.d0, 0.d0),work1)
          CALL prealloc_mult(work1,Gr(i+1,n),Gr(i,n))
          CALL destroy(work1)

        ENDIF

        ! NOW WHEN BLOCK IS SMALLER THAN EPS IT IS NOT CREATED
        !   CALL create(Gr(i,n),nrow,ncol)
        !   Gr(i,n)%val(:,:)=(0.d0,0.d0)
        !ENDIF

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
  !  Calculate Green Retarded - writing on memory (optimized on mask)
  !
  !****************************************************************************
  
  SUBROUTINE Gr_blk2csr(P,nbl,indblk,A)

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
    A%nzval = (0.d0,0.d0)

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
            write(*,*)     
            write(*,*) 'ERROR in blk2csr: probably wrong PL size',x
            write(*,*) 'row',i,A%colind(j)
            write(*,*) 'block indeces:',indblk(1:nbl)
            stop
          ENDIF

          col = A%colind(j) - indblk(y) + 1

          A%nzval(j) = Gr(x,y)%val(i1,col) 

        ENDDO

      ENDDO

    ENDIF

    !if (debug) call writePeakInfo(6)    
    if (debug) then
      WRITE(*,*) '**********************'
      WRITE(*,*) 'Make_GreenR_mem done'
      WRITE(*,*) '**********************'
    endif

  END SUBROUTINE Gr_blk2csr

  
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
    TYPE(z_DNS), DIMENSION(:,:), intent(in) :: ESH
    TYPE(z_DNS), DIMENSION(:,:), intent(inout) :: Gn
    TYPE(z_DNS), DIMENSION(:), intent(in) :: SelfEneR
    TYPE(Tstruct_info), intent(in) :: struct
    REAL(dp), DIMENSION(:), intent(in) :: frm
    INTEGER, intent(in) :: ref 

    !Work
    Type(z_DNS) :: Gam
    TYPE(z_DNS) :: work1,Ga
    INTEGER :: i,j,cb
    INTEGER :: ncont, nbl
    INTEGER, DIMENSION(:), POINTER :: cblk
    COMPLEX(dp) :: frmdiff

    ncont = struct%num_conts    
    nbl = struct%num_PLs
    cblk => struct%cblk

    !*******************************************
    ! Contact Iteration
    !*******************************************
    DO j=1,ncont

      ! NOTE: this soubroutine uses prealloc_mult that performs 
      ! C = C + A*B 
      IF (j.NE.ref .AND. ABS(frm(j)-frm(ref)).GT.EPS) THEN

        cb=cblk(j) ! block corresponding to contact j

        CALL zspectral(SelfEneR(j),SelfEneR(j),0,Gam)

        frmdiff = cmplx(frm(j)-frm(ref),0.d0,dp)
        ! Computation of Gl(1,1) = Gr(1,cb) Gam(cb) Ga(cb,1)
        if (allocated(Gr(1,cb)%val)) then
          CALL prealloc_mult(Gr(1,cb),Gam,work1)    
          CALL zdagger(Gr(1,cb),Ga)
          CALL prealloc_mult(work1,Ga,frmdiff,Gn(1,1))
          CALL destroy(work1, Ga)
        else
          Gn(1,1)%val=(0.d0,0.d0)
        endif

        ! Computation of all tridiagonal blocks
        DO i=2,nbl

          ! Computation of Gl(i,j) = Gr(i,cb) Gam(cb) Ga(cb,j)
          ! Both Gr(i,cb) and Gr(j,cb) must be non-zero
          if (Gr(i-1,cb)%nrow.gt.0 .and. Gr(i,cb)%nrow.gt.0) then
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
          else
            Gn(i-1,i)%val=(0.d0,0.d0)
            Gn(i,i-1)%val=(0.d0,0.d0)
          endif

        ENDDO

        call destroy(Gam)

      ENDIF

    ENDDO

  END SUBROUTINE Make_Gn_mem_dns

  !****************************************************************************
  !
  ! Calculate G_n contributions for all contacts (except reference)
  ! This version computes Grcol on the fly
  !
  !****************************************************************************
  
  SUBROUTINE Make_Gn_mem_dns2(ESH,SelfEneR,frm,ref,struct,Gn)

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
    TYPE(z_DNS), DIMENSION(:,:), intent(in) :: ESH
    TYPE(z_DNS), DIMENSION(:,:), intent(inout) :: Gn
    TYPE(z_DNS), DIMENSION(:), intent(in) :: SelfEneR
    TYPE(Tstruct_info), intent(in) :: struct
    REAL(dp), DIMENSION(:), intent(in) :: frm
    INTEGER, intent(in) :: ref 

    !Work
    Type(z_DNS) :: Gam
    TYPE(z_DNS) :: work1,Ga
    INTEGER :: i,j,cb
    INTEGER :: ncont, nbl
    INTEGER, DIMENSION(:), POINTER :: cblk
    COMPLEX(dp) :: frmdiff

    ncont = struct%num_conts    
    nbl = struct%num_PLs
    cblk => struct%cblk

    !*******************************************
    ! Contact Iteration
    !*******************************************
    DO j=1,ncont

      ! NOTE: this soubroutine uses prealloc_mult that performs 
      ! C = C + A*B 
      IF (j.NE.ref .AND. ABS(frm(j)-frm(ref)).GT.EPS) THEN

        cb=cblk(j) ! block corresponding to contact j

        CALL zspectral(SelfEneR(j),SelfEneR(j),0,Gam)

        frmdiff = cmplx(frm(j)-frm(ref),0.d0,dp)
        ! Computation of Gl(1,1) = Gr(1,cb) Gam(cb) Ga(cb,1)
        if (Gr(1,cb)%nrow.gt.0) then
          CALL prealloc_mult(Gr(1,cb),Gam,work1)    
          CALL zdagger(Gr(1,cb),Ga)
          CALL prealloc_mult(work1,Ga,frmdiff,Gn(1,1))
          CALL destroy(work1, Ga)
        else
          Gn(1,1)%val=(0.d0,0.d0)
        endif

        ! Computation of all tridiagonal blocks
        DO i=2,nbl

          ! Computation of Gr(i,cb) assuming Gr(i-1,cb) exists
          ! Assume downgoing: i > cb
          IF (Gr(i-1,cb)%nrow.GT.0) THEN
            CALL prealloc_mult(gsmr(i),ESH(i,i-1),(-1.d0, 0.d0),work1)
            CALL destroy(gsmr(i))
            CALL prealloc_mult(work1,Gr(i-1,cb),Gr(i,cb))
            CALL destroy(work1)
            IF (MAXVAL(ABS(Gr(i,cb)%val)).lt.EPS) call destroy(Gr(i,cb))
          ENDIF

          ! Computation of Gl(i,j) = Gr(i,cb) Gam(cb) Ga(cb,j)
          ! Both Gr(i,cb) and Gr(j,cb) must be non-zero
          IF (Gr(i-1,cb)%nrow.gt.0 .and. Gr(i,cb)%nrow.gt.0) THEN 
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
          ELSE
            Gn(i-1,i)%val=(0.d0,0.d0)
            Gn(i,i-1)%val=(0.d0,0.d0)
          ENDIF

          IF (Gr(i-1,cb)%nrow.gt.0) call destroy(Gr(i-1,cb))

        ENDDO

        call destroy(Gam)

      ENDIF

    ENDDO

  END SUBROUTINE Make_Gn_mem_dns2
  
  !****************************************************************************
  !>
  ! Calculate G_n contributions due to el-ph
  ! Writing on memory
  !
  !****************************************************************************
  SUBROUTINE Make_Gn_ph(pnegf,ESH,iter,Gn)

    TYPE(Tnegf), intent(in) :: pnegf
    TYPE(z_DNS), DIMENSION(:,:), intent(in) :: ESH
    TYPE(z_DNS), DIMENSION(:,:), intent(inout) :: Gn
    INTEGER, intent(in) :: iter

    INTEGER, DIMENSION(:), POINTER :: indblk
    Type(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_ph_n, sigma_blk
    Type(z_DNS) :: Ga, work1, work2, sigma_tmp
    INTEGER :: n, k, nbl, nrow, ierr, ii, jj, norbs, nblk, indstart, indend

    !! If this is the first scba cycle, there's nothing to do
    if (pnegf%inter%scba_iter .eq. 0) then
      return
    endif
    nbl = pnegf%str%num_PLs
    ALLOCATE(Sigma_ph_n(nbl,nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_n'

    ! The block sigma n is made available from el-ph model
    ! Note: the elph models could not keep a copy and calculate it 
    ! on the fly. You have to rely on the local copy
    if (allocated(pnegf%inter)) then
      call pnegf%inter%get_sigma_n(Sigma_ph_n, pnegf%ie)
    end if

    !! Make the sigma_ph_n available.
    !! The exact way depends on the el-ph model
    !! I could do it on the fly to save some memory
    select case(pnegf%elph%model)
    case (0)
      continue
    case (1)
      stop
    case (2)
      stop "Deprecated"
    case(3)
      stop "Deprecated"

    case default
      write(*,*) 'Not yet implemented'
      stop 0

      !! old alex implementation, not active
      if (iter .gt. 0) then
        call read_blkmat(Sigma_ph_n(n,n),pnegf%scratch_path,'Sigma_ph_n_',n,n,pnegf%iE)
      else 
        call write_blkmat(Sigma_ph_n(n,n),pnegf%scratch_path,'Sigma_ph_n_',n,n,pnegf%iE)
      endif
      !! old alex implementation, not active
    end select

    !! Calculate the diagonal and off diagonal (if needed) blocks of Gn
    !! in the assumption of diagonal self energy
    !! G(k,k) = Gr(k,i)Sigma_n(i,i)Ga(i,k)
    !! G(k,k+1) = Gr(k,i)Sigma_n(i,i)Ga(i,k+1)
    !! G(k,k-1) = Gr(k,i)Sigma_n(i,i)Ga(i,k-1)
    !! All the rows of Gr need to be available
    DO n = 1, nbl-1
      DO k = 1, nbl
        if (Gr(n,k)%nrow.gt.0) then
          CALL zdagger(Gr(n,k),Ga)
          CALL prealloc_mult(Gr(n,k), Sigma_ph_n(k,k), work1)
          CALL prealloc_mult(work1, Ga, work2)
          ! Computing diagonal blocks of Gn(n,n)
          Gn(n,n)%val = Gn(n,n)%val + work2%val
          call destroy(work2,Ga)
        endif
        ! Computing blocks of Gn(n,n+1)
        ! Only if S is not identity: Gn is initialized on ESH therefore
        ! we need to check the number of rows (or column)
        if (Gr(n+1,k)%nrow.gt.0 .and. Gn(n,n+1)%nrow .gt. 0) then
          CALL zdagger(Gr(n+1,k),Ga)
          CALL prealloc_mult(work1, Ga, work2)
          Gn(n,n+1)%val = Gn(n,n+1)%val + work2%val
          call destroy(work1,work2,Ga)         
          Gn(n+1,n)%val = conjg(transpose(Gn(n,n+1)%val))
        endif
      END DO
    END DO
    DO k = 1, nbl
      if (Gr(nbl,k)%nrow.gt.0) then
        CALL zdagger(Gr(nbl,k),Ga)
        CALL prealloc_mult(Gr(nbl,k), Sigma_ph_n(k,k), work1)
        CALL prealloc_mult(work1, Ga, work2)
        Gn(nbl,nbl)%val = Gn(nbl,nbl)%val + work2%val
        call destroy(work1,work2,Ga)
      endif
    END DO
    !! End Gn calculation


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
!!$  SUBROUTINE Make_Gp_ph(pnegf,ESH,iter,Gp)
!!$
!!$    TYPE(Tnegf) :: pnegf
!!$    TYPE(z_DNS), DIMENSION(:,:) :: ESH, Gp
!!$    INTEGER :: iter
!!$
!!$    Type(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_ph_p
!!$    Type(z_DNS) :: Ga, work1, work2 
!!$    INTEGER :: n, k, nbl, nrow, ierr
!!$
!!$    nbl = pnegf%str%num_PLs
!!$    ALLOCATE(Sigma_ph_p(nbl,nbl),stat=ierr)
!!$    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_n'
!!$
!!$
!!$    DO n = 1, nbl
!!$
!!$      nrow = ESH(n,n)%nrow
!!$
!!$      call create(Sigma_ph_p(n,n), nrow, nrow)
!!$
!!$      Sigma_ph_p(n,n)%val = (0.0_dp, 0.0_dp)
!!$      if (iter .gt. 0) then
!!$         call read_blkmat(Sigma_ph_p(n,n),pnegf%scratch_path,'Sigma_ph_p_',n,n,pnegf%iE)
!!$      else 
!!$         call write_blkmat(Sigma_ph_p(n,n),pnegf%scratch_path,'Sigma_ph_p_',n,n,pnegf%iE)
!!$      endif
!!$    
!!$    END DO
!!$
!!$    DO n = 1, nbl
!!$
!!$      DO k = 1, nbl
!!$
!!$         if (Gr(n,k)%nrow.gt.0) then
!!$            CALL zdagger(Gr(n,k),Ga)
!!$            CALL prealloc_mult(Gr(n,k), Sigma_ph_p(k,k), work1)
!!$            CALL prealloc_mult(work1, Ga, work2)
!!$            Gp(n,n)%val = Gp(n,n)%val + work2%val
!!$            call destroy(work1, work2, Ga)
!!$         endif
!!$
!!$      END DO    
!!$    
!!$    END DO
!!$    
!!$    DO n = 1, nbl
!!$      CALL destroy(Sigma_ph_p(n,n))
!!$    END DO
!!$
!!$    DEALLOCATE(Sigma_ph_p)
!!$
!!$  END SUBROUTINE Make_Gp_ph
  
  !Concatenation for every contact in G_n. Performs a sum on elements, not a replacement
  !Similar to Make_GreenR_mem2, except for sum of elements
  !Note: to backup old version zconcat calls (and Glsub deallocations) must be
  !      uncommented and all this part removed
  !If only one block is present, concatenation is not needed and it's implemented in a
  !more trivial way
  SUBROUTINE blk2csr(G,struct,P,Gcsr)

    TYPE(z_DNS), DIMENSION(:,:) :: G
    TYPE(Tstruct_info), intent(in) :: struct
    TYPE(z_CSR) :: Gcsr
    TYPE(z_CSR) :: P, G_sp

    INTEGER, DIMENSION(:), POINTER :: indblk, cblk
    INTEGER :: nbl, oldx, row, col, iy, ix, x, y, ii, jj, nrows

    nbl = struct%num_PLs
    cblk => struct%cblk
    indblk => struct%mat_PL_start
    nrows = struct%mat_PL_end(nbl)

    !create Gcsr with same pattern of P
    CALL create(Gcsr,P%nrow,P%ncol,P%nnz)
    Gcsr%rowpnt = P%rowpnt
    Gcsr%colind = P%colind
    Gcsr%nzval = (0.d0, 0.d0)

    !Cycle upon all rows
    x = 1
    DO ii = 1, nrows
      !Search block x containing row ii
      oldx = x
      IF (oldx.EQ.nbl) THEN 
        x = oldx
      ELSE
        DO ix = oldx, oldx+1
          IF ( (ii.GE.indblk(ix)).AND.(ii.LT.indblk(ix+1)) ) x = ix
        ENDDO
      ENDIF

      !Offset: row is the index for separate blocks
      row = ii - indblk(x) + 1

      !Cycle upon columns of Gcsr (which has been ALREADY MASKED by S)
      DO jj = Gcsr%rowpnt(ii), Gcsr%rowpnt(ii+1) -1
        IF (Gcsr%colind(jj).gt.nrows) CYCLE
        !Choose which block column we're dealing with
        y = 0
        IF (x.eq.1) then
          IF ( (Gcsr%colind(jj).GE.indblk(x)).AND.(Gcsr%colind(jj).LT.indblk(x + 1)) ) then 
            y = 1
          ELSEIF ( (Gcsr%colind(jj).GE.indblk(x + 1)).AND.(Gcsr%colind(jj).LT.indblk(x + 2)) ) then 
            y = 2
          ENDIF
        elseif (x.eq.nbl) then
          IF ( (Gcsr%colind(jj).GE.indblk(x)).AND.(Gcsr%colind(jj).LT.indblk(x + 1)) ) then 
            y = nbl
          ELSEIF ( (Gcsr%colind(jj).GE.indblk(x - 1)).AND.(Gcsr%colind(jj).LT.indblk(x)) ) then 
            y = nbl - 1
          ENDIF
        else
          DO iy = x-1, x+1
            if ( (Gcsr%colind(jj).GE.indblk(iy)).AND.(Gcsr%colind(jj).LT.indblk(iy + 1)) ) y = iy
          ENDDO
        ENDIF

        IF (y.EQ.0) THEN
          write(*,*)     
          write(*,*) 'ERROR in blk2csr: probably wrong PL size', x
          write(*,*) 'row',ii,nrows,Gcsr%colind(jj)
          write(*,*) 'block indeces:',indblk(1:nbl)
          STOP
        ENDIF

        col = Gcsr%colind(jj) - indblk(y) + 1

        IF (allocated(G(x,y)%val)) THEN
          Gcsr%nzval(jj) = Gcsr%nzval(jj) + G(x,y)%val(row,col)
        ENDIF

      ENDDO

    ENDDO


  END SUBROUTINE blk2csr

  ! ******************************************************************************
  ! Computes Sigma_ph_n and save it file
  ! ******************************************************************************
  SUBROUTINE Sigma_ph_n(pnegf,Epnt)

    TYPE(Tnegf) :: pnegf
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

    TYPE(Tnegf) :: pnegf
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

  !--------------------------------------------------------------------------------
  SUBROUTINE Sigma_ph_r(pnegf,Epnt)

    TYPE(Tnegf) :: pnegf
    REAL(dp), DIMENSION(:),allocatable :: Epnt

    !Local variables
    TYPE(z_DNS) :: work1,work2
    TYPE(z_DNS) :: Mq_mat
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_r, G_r_interP, G_r_interN
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: G_n_interP, G_n_interN

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
      !call read_blkmat(G_r(i,i),pnegf%scratch_path,'G_r_',i,i,iE)
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

    TYPE(Tnegf) :: pnegf
    COMPLEX(dp), intent(in) :: z

    !Local variables
    TYPE(z_DNS) :: work1,work2
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

        call read_blkmat(G_r(i,i),pnegf%scratch_path,'G_r_',i,i,iE)

        i1 = Sigma_r(i,i)%nrow

        if (pnegf%elph%selfene_gr) then

          call create(work1,i1,i1)

          work1%val = (0.d0,0.d0)
          !print*,'SelfEneR_Gr'    
          work1%val = (2.0_dp*Nq(m)+1.0_dp)*G_r(i,i)%val

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
      !print*,'(sigma_r) Sigma_ph_r', pnegf%iE, maxval(abs(Sigma_r(i,i)%val))
      call write_blkmat(Sigma_r(i,i),pnegf%scratch_path,'Sigma_ph_r_',i,i,iE)
      call destroy(Sigma_r(i,i))
      call destroy(G_r(i,i))
    enddo

  END SUBROUTINE Sigma_ph_r_z
  ! ----------------------------------------------------------
  SUBROUTINE check_Gl_Gr(pnegf)
    TYPE(Tnegf) :: pnegf
    REAL(dp), DIMENSION(:),allocatable :: Epnt
    integer :: ioffset

    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: G_r, G_p, G_n
    TYPE(z_DNS) :: A, T 
    integer :: nbl, n, sizebl, i_start, i_stop, psize, maxpos(2)
    real(dp) :: Wmax, maxdev, tmp, maxG
    INTEGER, DIMENSION(:), pointer :: indblk

    nbl = pnegf%str%num_PLs
    indblk => pnegf%str%mat_PL_start

    allocate(G_n(nbl,nbl))
    allocate(G_p(nbl,nbl))
    allocate(G_r(nbl,nbl))

    maxdev = 0.0_dp
    maxG = 0.0_dp
    psize = 0

    do n = 1, nbl
      sizebl = indblk(n+1)-indblk(n)
      call create(G_r(n,n),sizebl,sizebl)
      call create(G_n(n,n),sizebl,sizebl)
      call create(G_p(n,n),sizebl,sizebl)
      call read_blkmat(G_r(n,n), pnegf%scratch_path, 'G_r_',n,n, pnegf%iE)
      call read_blkmat(G_n(n,n),pnegf%scratch_path,'G_n_',n,n,pnegf%iE)
      call read_blkmat(G_p(n,n),pnegf%scratch_path,'G_p_',n,n,pnegf%iE)

      call zspectral(G_r(n,n),G_r(n,n),0,A)

      call create(T,sizebl,sizebl)

      !print*,'CHECK Sigma_ph_n+Sigma_ph_p',maxval(abs(Sigma_n(n,n)%val+Sigma_p(n,n)%val))

      T%val = G_n(n,n)%val + G_p(n,n)%val - A%val 

      tmp = maxval(abs(A%val))
      if (tmp .gt. maxG) maxG=tmp

      tmp = maxval(abs(T%val))/maxval(abs(A%val)) 
      if (tmp .gt. maxdev) then
        maxdev = tmp
        maxpos = maxloc(abs(T%val)) + psize
      endif

      !print*

      psize = psize + sizebl

      call destroy(G_r(n,n))
      call destroy(G_n(n,n))
      call destroy(G_p(n,n))
      call destroy(A)
      call destroy(T)

    enddo

    print*,'CHECK Gn+Gp=Gr-Ga',pnegf%iE, maxG, maxdev

    deallocate(G_n, G_p, G_r)

  END SUBROUTINE check_Gl_Gr

  ! ----------------------------------------------------------
  SUBROUTINE check_sigma_ph_r(pnegf)
    TYPE(Tnegf) :: pnegf

    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_r, Sigma_p, Sigma_n
    TYPE(z_DNS) :: Gam, T 
    integer :: nbl, n, sizebl, psize, maxpos(2)
    real(dp) :: maxdev, tmp, maxG
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
    type(TNegf) :: pnegf
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
  !Need to set back FFT transforms
  !********************************************************************
  SUBROUTINE complete_sigma_ph_r(pnegf, Epnt, ioffset)

    TYPE(Tnegf) :: pnegf
    REAL(dp), DIMENSION(:) :: Epnt
    INTEGER :: ioffset

    ! Locals
    INTEGER :: i,j,n,m, iE, bl, sizebl,i_start,i_stop, ierr, nummodes
    REAL(dp), DIMENSION(:), POINTER :: Wq
    REAL(dp), DIMENSION(:), POINTER :: Mq
    REAL(dp) :: Wmax, dE, tmp

    character(4) :: ofblki
    complex(dp), DIMENSION(:), allocatable :: temp1 
    type(z_DNS), DIMENSION(:), allocatable :: Gn_E, Sigma_r_E
    !type(z_DNS) :: work1, work2, work3, Mq_mat 
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

            !call Hilbert_shift(temp1, Wq(m))

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

  SUBROUTINE create_scratch(nbl, npoints)
    integer :: nbl, npoints

    integer :: i,j,k,err

    call destroy_scratch(nbl, npoints)

    ALLOCATE(GGn(nbl,nbl,npoints),stat=err)
    ALLOCATE(GGr(nbl,nbl,npoints),stat=err)
    ALLOCATE(GGp(nbl,nbl,npoints),stat=err)
    ALLOCATE(Sgn(nbl,nbl,npoints),stat=err)
    ALLOCATE(Sgr(nbl,nbl,npoints),stat=err)
    ALLOCATE(Sgp(nbl,nbl,npoints),stat=err)

    ! Initialize everything to 0 
    do k=1,npoints
      do j=1,nbl
        do i=1,nbl  
          GGn(i,j,k)%nrow=0
          GGn(i,j,k)%ncol=0
          GGr(i,j,k)%nrow=0
          GGr(i,j,k)%ncol=0
          GGp(i,j,k)%nrow=0
          GGp(i,j,k)%ncol=0
          Sgn(i,j,k)%nrow=0
          Sgn(i,j,k)%ncol=0
          Sgr(i,j,k)%nrow=0
          Sgr(i,j,k)%ncol=0
          Sgp(i,j,k)%nrow=0
          Sgp(i,j,k)%ncol=0
        enddo
      enddo
    enddo

    if(err.ne.0) then
      STOP 'ERROR: Cannot allocate GG'
    endif
    print*,'Created memory scratch',nbl,'x',nbl,'x',npoints

  END SUBROUTINE create_scratch

  SUBROUTINE destroy_scratch(nbl, npoints)
    integer :: nbl, npoints

    integer :: i, j, iE,  err

    err = 0


    do i = 1, nbl
      do j = 1, nbl
        do iE = 1, npoints

          if (allocated(GGn)) then
            if (allocated(GGn(i,j,iE)%val)) call destroy(GGn(i,j,iE))
          endif
          if (allocated(GGr)) then
            if (allocated(GGr(i,j,iE)%val)) call destroy(GGr(i,j,iE))
          endif
          if (allocated(GGp)) then
            if (allocated(GGp(i,j,iE)%val)) call destroy(GGp(i,j,iE))
          endif
          if (allocated(Sgn)) then
            if (allocated(Sgn(i,j,iE)%val)) call destroy(Sgn(i,j,iE))     
          endif
          if (allocated(Sgr)) then
            if (allocated(Sgr(i,j,iE)%val)) call destroy(Sgr(i,j,iE))
          endif
          if (allocated(Sgp)) then
            if (allocated(Sgp(i,j,iE)%val)) call destroy(Sgp(i,j,iE))     
          endif

        end do
      end do
    end do

    if (allocated(GGn)) DEALLOCATE(GGn,stat=err)
    if (allocated(GGr)) DEALLOCATE(GGr,stat=err)
    if (allocated(GGp)) DEALLOCATE(GGp,stat=err)
    if (allocated(Sgn)) DEALLOCATE(Sgn,stat=err)
    if (allocated(Sgr)) DEALLOCATE(Sgr,stat=err)
    if (allocated(Sgp)) DEALLOCATE(Sgp,stat=err)

    if(err.ne.0) then
      STOP 'ERROR: Cannot deallocate GG'
    endif

  END SUBROUTINE destroy_scratch

  ! READ Matrices
  SUBROUTINE read_blkmat(Matrix, path, name, i, j, iE)
    TYPE(z_DNS), intent(inout) :: Matrix
    CHARACTER(*), intent(in) :: path
    CHARACTER(*), intent(in) :: name
    INTEGER, intent(in) :: i, j, iE

    CHARACTER(64) :: filename
    character(4) :: ofblki, ofblkj
    character(10) :: ofpnt
    logical :: lex
    integer :: i1,i2
    complex(dp) :: mat_el

    if (memory) then
      select case(trim(name))
      case('G_n_')     
        Matrix%val = GGn(i,j,iE)%val 
      case('G_r_')     
        Matrix%val = GGr(i,j,iE)%val
      case('G_p_')     
        Matrix%val = GGp(i,j,iE)%val
      case('Sigma_ph_r_')     
        Matrix%val = Sgr(i,j,iE)%val 
      case('Sigma_ph_n_')     
        Matrix%val = Sgn(i,j,iE)%val
      case('Sigma_ph_p_')     
        Matrix%val = Sgp(i,j,iE)%val
      case default 
        stop 'internal error: read_blkmat does not correspond'
      end select
      return
    endif

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

    open(9091,file=trim(path)//trim(filename), access='STREAM')

    call inmat_c(9091,.false.,Matrix%val,Matrix%nrow,Matrix%ncol)

    !open(9001,file=trim(path)//trim(filename), access='DIRECT', recl=4)
    !call direct_in_c(9001,Matrix%val,Matrix%nrow)

    close(9091)

  END SUBROUTINE read_blkmat

  ! WRITE Matrices
  SUBROUTINE write_blkmat(Matrix, path, name, i, j, iE)
    TYPE(z_DNS), intent(in) :: Matrix
    CHARACTER(*), intent(in) :: path
    CHARACTER(*), intent(in) :: name
    INTEGER, intent(in) :: i, j, iE

    CHARACTER(64) :: filename
    character(4) :: ofblki, ofblkj
    character(10) :: ofpnt
    integer :: m,n

    if (memory) then

      select case(trim(name))
      case('G_n_')     
        if (.not.allocated(GGn(i,j,iE)%val)) then
          call create(GGn(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        GGn(i,j,iE)%val = Matrix%val
      case('G_r_')     
        if (.not.allocated(GGr(i,j,iE)%val)) then
          call create(GGr(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        GGr(i,j,iE)%val = Matrix%val
      case('G_p_')     
        if (.not.allocated(GGp(i,j,iE)%val)) then
          call create(GGp(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        GGp(i,j,iE)%val = Matrix%val 
      case('Sigma_ph_r_')     
        if (.not.allocated(Sgr(i,j,iE)%val)) then
          call create(Sgr(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        Sgr(i,j,iE)%val = Matrix%val
      case('Sigma_ph_n_')     
        if (.not.allocated(Sgn(i,j,iE)%val)) then
          call create(Sgn(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        Sgn(i,j,iE)%val = Matrix%val
      case('Sigma_ph_p_')     
        if (.not.allocated(Sgp(i,j,iE)%val)) then
          call create(Sgp(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        Sgp(i,j,iE)%val = Matrix%val
      case default 
        stop 'internal error: write_blkmat does not correspond'
      end select
      return
    endif

    if (i.le.9999) write(ofblki,'(i4.4)') i
    if (i.gt.9999) stop 'ERROR: too many blks (> 9999)'
    if (j.le.9999) write(ofblkj,'(i4.4)') j
    if (j.gt.9999) stop 'ERROR: too many blks (> 9999)'

    if (iE.le.99999) write(ofpnt,'(i5.5)') iE

    filename = trim(name)//trim(ofblki)//'_'//trim(ofblkj)//'_'//trim(ofpnt)//'.dat'

    open(9001,file=trim(path)//trim(filename), access='STREAM', status='REPLACE')

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
    TYPE(z_CSR) :: Aout


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

    IF (.not.allocated(Aout%nzval)) THEN
      CALL create(Aout,nrow_tot,nrow_tot,0)
      Aout%rowpnt(:)=1
    ENDIF

    DO i=1,ncont

      !Numero di blocco del contatto
      cb=cblk(i)
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

      !Concatenazione di Asub nella posizione corrispondente
      i1=indblk(cb)
      j1=struct%mat_B_start(i)

      CALL concat(Aout,GrCSR,i1,j1)   

      CALL destroy(GrCSR)

      IF (lower) THEN

        CALL prealloc_mult(gsurfR(i),Tcl(i),(-1.d0, 0.d0), work1)
        CALL prealloc_mult(work1, Gr(cb,cb), Grcl)

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

        CALL concat(Aout,GrCSR,i1,j1)

        CALL destroy(GrCSR)

      ENDIF

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
    IF (.not.allocated(Glout%nzval)) THEN
      CALL create(Glout,nrow_tot,nrow_tot,0)
      Glout%rowpnt(:)=1
    ENDIF
    !***
    !Iterazione su tutti i contatti "k" 
    !***
    DO k=1,ncont

      !Esegue le operazioni relative al contatto solo se e` valida la condizione
      !sulle distribuzioni di Fermi e se non si tratta del contatto iniettante (ref)
      IF ((ABS(frm(k)-frm(ref)).GT.EPS).AND.(k.NE.ref)) THEN

        !Calcolo della Gamma corrispondente
        cb=cblk(k)
        !nota:cb indica l'indice di blocco corrispondente al contatto j-esimo

        CALL zspectral(SelfEneR(k),SelfEneR(k),0,Gam)

        !***
        !Calcolo del contributo sulla proria regione
        !***          
        frmdiff=cmplx(frm(ref)-frm(k))

        CALL zspectral(gsurfR(k),gsurfR(k),0,work1)

        !print*, 'work2=Tlc*work1=Tlc*j(gsurfR-gsurfA)'
        CALL prealloc_mult(Tlc(k),work1,work2)
        CALL destroy(work1)

        !print *, 'work1=-Gr(cb,cb)*work2=-Gr(cb,cb)*Tlc*j(gsurfR-gsurfA)'
        CALL prealloc_mult(Gr(cb,cb),work2,frmdiff,work1)
        CALL destroy(work2)

        !print*,'work2=Tlc*gsurfA'
        CALL zdagger(gsurfR(k),gsurfA)
        CALL prealloc_mult(Tlc(k),gsurfA,work2)
        CALL destroy(gsurfA)

        !print*,'work3=Ga*work2=Ga*Tlc*gsurfA'           
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

        !Contributo totale sulla propria regione
        CALL prealloc_sum(work3,work1,Glsub)
        CALL destroy(work1)
        CALL destroy(work3)

        call mask(Glsub,Tlc(k)) 
        i1=nzdrop(Glsub,EPS) 

        IF (i1.gt.0) THEN
          CALL create(GlCSR,Glsub%nrow,Glsub%ncol,i1)
          CALL dns2csr(Glsub,GlCSR)

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
        END IF

        CALL destroy(Glsub)
        !***
        !Ciclo per il calcolo dei contributi sulle regioni degli altri contatti
        !***

        DO j=1,ncont

          cbj=cblk(j)
          !Esegue le operazioni del ciclo solo se il j.ne.k o se
          !il blocco colonna di Gr e` non nullo (altrimenti il contributo e` nullo) 

          IF ((j.NE.k).AND.(Gr(cbj,cb)%nrow.NE.0 .AND. (Gr(cbj,cb)%ncol.NE.0))) THEN

            !print*,'work1=Tlc*gsurfA'  
            CALL zdagger(gsurfR(j),gsurfA)
            CALL prealloc_mult(Tlc(j),gsurfA,work1)
            CALL destroy(gsurfA)

            !print*,'work2=Ga*work1=Ga*Tlc*gsurfA'  
            CALL zdagger(Gr(cbj,cb),Ga)
            CALL prealloc_mult(Ga,work1,work2)

            CALL destroy(Ga)
            CALL destroy(work1)

            !print*,'work1=Gam*work2=Gam*Ga*Tlc*gsurfA'  
            CALL prealloc_mult(Gam,work2,work1)
            CALL destroy(work2)

            !print*,'Glsub=-Gr*work1=-Gr*Gam*Ga*Tlc*gsurfA'  
            CALL prealloc_mult(Gr(cbj,cb),work1,frmdiff,Glsub)
            CALL destroy(work1)

            call mask(Glsub,Tlc(j)) 
            i1=nzdrop(Glsub,EPS) 

            IF (i1.gt.0) THEN
              CALL create(GlCSR,Glsub%nrow,Glsub%ncol,i1)
              CALL dns2csr(Glsub,GlCSR)

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
            ENDIF

            CALL destroy(Glsub)

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

  SUBROUTINE tunneling_dns(H,S,Ec,SelfEneR,ni,nf,size_ni,str,TUN_MAT)

    implicit none

    Type(z_CSR) :: H
    Type(z_CSR) :: S           
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
    call prealloc_sum(H,S,(-1.d0, 0.d0),Ec,ESH_tot)    

    call allocate_blk_dns(ESH,nbl)

    call zcsr2blk_sod(ESH_tot,ESH,str%mat_PL_start)
    call destroy(ESH_tot)

    !Inclusion of the contact Self-Energies to the relevant blocks
    do i=1,ncont
      ESH(str%cblk(i),str%cblk(i))%val = ESH(str%cblk(i),str%cblk(i))%val-SelfEneR(i)%val
    enddo
    call allocate_gsm_dns(gsmr,nbl)

    !Iterative calculation up with gsmr 
    call Make_gsmr_mem_dns(ESH,nbl,2)


    !Computation of transmission(s) between contacts ni(:) -> nf(:)  
    nit=ni(1)
    nft=nf(1)

    !Arrange contacts in a way that the order between first and second is always the
    !same (always ct1 < ct2), so nt is the largest contact block 

    if (str%cblk(nit).gt.str%cblk(nft)) then
      nt = str%cblk(nit)
    else
      nt = str%cblk(nft)
    endif

    ! Iterative calculation of Gr down to nt1 
    call allocate_blk_dns(Gr,nbl)
    call Make_Gr_mem_dns(ESH,1)
    if (nt.gt.1) call Make_Gr_mem_dns(ESH,2,nt)

    if (ncont.eq.2) then
      call trasmission_dns(nit,nft,ESH,SelfEneR,str%cblk,tun) 
    else
      call trasmission_old(nit,nft,ESH,SelfEneR,str%cblk,tun) 
    endif

    TUN_MAT(1) = tun 

    ! When more contacts are present sometimes we can re-use previous Gr 
    do icpl = 2, size_ni

      nit=ni(icpl)
      nft=nf(icpl)

      if (str%cblk(nit).gt.str%cblk(nft)) then
        nt1 = str%cblk(nit)
      else
        nt1 = str%cblk(nft)
      endif

      ! if nt1 > nt extend the Gr calculation
      if (nt1 .gt. nt) then
        call Make_Gr_mem_dns(ESH,nt+1,nt1)
        nt = nt1
      endif

      call trasmission_old(nit,nft,ESH,SelfEneR,str%cblk,tun) 

      TUN_MAT(icpl) = tun 

    enddo

    !Distruzione delle Green
    do i=2, nt 
      call destroy(Gr(i,i))
      call destroy(Gr(i-1,i))
      call destroy(Gr(i,i-1))
    enddo
    call destroy(Gr(1,1))

    call deallocate_blk_dns(Gr)

    !Deallocate energy-dependent matrices
    call destroy_gsm(gsmr)
    call deallocate_gsm_dns(gsmr)

    call destroy_ESH(ESH)
    call deallocate_blk_dns(ESH)

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


  subroutine trasmission_dns(ni,nf,ESH,SelfEneR,cblk,TUN)

    implicit none

    !In/Out
    Integer :: ni,nf
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR
    Type(z_DNS), Dimension(:,:) :: ESH
    Integer, Dimension(:), pointer :: cblk
    Real(dp) :: TUN

    !Work variables
    Integer :: ct1, bl1
    Type(z_DNS) :: work1, work2, GAM1_dns, GA, TRS, AA
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

    bl1=cblk(ct1); 

    call zdagger(Gr(bl1,bl1),GA)

    ! Computes the Gamma matrices
    call zspectral(SelfEneR(ct1),SelfEneR(ct1),0,GAM1_dns)

    ! Work to compute transmission matrix (Gamma G Gamma G)
    call prealloc_mult(GAM1_dns,Gr(bl1,bl1),work1)

    call prealloc_mult(work1,GAM1_dns,work2)

    call destroy(work1)

    call prealloc_mult(work2,GA,work1)

    call destroy(work2)

    call create(AA,GA%nrow,GA%ncol)

    AA%val = j * (Gr(bl1,bl1)%val-GA%val)

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
  subroutine trasmission_old(ni,nf,ESH,SelfEneR,cblk,TUN)

    !In/Out
    Integer :: ni,nf
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR
    Type(z_DNS), Dimension(:,:) :: ESH
    Integer, Dimension(:), pointer :: cblk
    Real(kind=dp) :: TUN

    !Work variables
    Integer :: ct1, ct2, bl1, bl2, i, nbl
    Type(z_DNS) :: work1, work2, GAM1_dns, GAM2_dns, GA, TRS
    Real(kind=dp) :: max

    !Arrange contacts in way that order between first and second is always the
    !same (always ct1 < ct2)

    if (cblk(ni).lt.cblk(nf)) then
      ct1=ni;ct2=nf;
    else
      ct1=nf;ct2=ni;
    endif

    bl1=cblk(ct1); bl2=cblk(ct2);
    nbl = size(cblk)
    ! in this way nt1 < nt2 by construction

    if ( nbl.gt.1 .and. (bl2-bl1).gt.1) then

      ! Compute column-blocks of Gr(i,bl1) up to i=bl2
      ! Gr(i,bl1) = -gr(i) T(i,i-1) Gr(i-1,bl1)
      do i = bl1+1, bl2
        !Checks whether previous block is non null. 
        !If so next block is also null => TUN = 0       
        max=maxval(abs(Gr(i-1,bl1)%val))

        if (max.lt.EPS) then
          TUN = EPS*EPS !for log plots 
          !Destroy also the block adjecent to diagonal since 
          !this is not deallocated anymore in calling subroutine
          if (i.gt.(bl1+1)) call destroy(Gr(i-1,bl1))
          return
        endif

        !Checks whether block has been created, if not do it 
        if (.not.allocated(Gr(i,bl1)%val)) then 

          call prealloc_mult(gsmr(i),ESH(i,i-1),(-1.d0, 0.d0),work1)

          call prealloc_mult(work1,Gr(i-1,bl1),Gr(i,bl1))

          call destroy(work1)

        endif

        ! avoid destroying blocks closer to diagonal
        if (i.gt.(bl1+2)) call destroy(Gr(i-1,bl1))

      enddo

    endif

    ! Computes the Gamma matrices
    call zspectral(SelfEneR(ct1),SelfEneR(ct1),0,GAM1_dns)
    call zspectral(SelfEneR(ct2),SelfEneR(ct2),0,GAM2_dns)

    ! Work to compute transmission matrix (Gamma2 Gr Gamma1 Ga)
    call prealloc_mult(GAM2_dns,Gr(bl2,bl1),work1)

    call destroy(GAM2_dns)

    call prealloc_mult(work1,GAM1_dns,work2)

    call destroy(work1)

    call destroy(GAM1_dns)

    call zdagger(Gr(bl2,bl1),GA)

    if (bl2.gt.bl1+1) call destroy( Gr(bl2,bl1) )

    call prealloc_mult(work2,GA,TRS)

    call destroy(work2)

    call destroy(GA) 

    TUN = real(trace(TRS))

    call destroy(TRS)

  end subroutine trasmission_old


  !---------------------------------------------------!
  !Subroutine for transmission and dos calculation    !
  !---------------------------------------------------!  

  subroutine tun_and_dos(H,S,Ec,SelfEneR,Gs,ni,nf,nLDOS,LDOS,size_ni,str,TUN_MAT,LEDOS)

    implicit none

    Type(z_CSR), intent(in) :: H
    Type(z_CSR), intent(in) :: S           
    Complex(dp), intent(in) :: Ec
    Type(z_DNS), Dimension(MAXNCONT), intent(in) :: SelfEneR, Gs
    Integer, intent(in) :: ni(MAXNCONT)
    Integer, intent(in) :: nf(MAXNCONT)
    Type(TStruct_Info), intent(in) :: str
    integer, intent(in)  :: nLdos, size_ni
    type(intarray), dimension(:), intent(in) :: LDOS      
    Real(dp), Dimension(:), intent(inout) :: TUN_MAT
    Real(dp), Dimension(:), intent(inout) :: LEDOS

    ! Local variables
    Type(z_CSR) :: ESH_tot, GrCSR
    Type(z_DNS), Dimension(:,:), allocatable :: ESH
    Type(r_CSR) :: Grm                          ! Green Retarded nella molecola           
    real(dp), dimension(:), allocatable :: diag
    Real(dp) :: tun
    Complex(dp) :: zc
    Integer :: nbl,ncont, ierr
    Integer :: nit, nft, icpl
    Integer :: iLDOS, i2, i
    Character(1) :: Im

    nbl = str%num_PLs
    ncont = str%num_conts
    Im = 'I'

    !Calculation of ES-H and brak into blocks
    call prealloc_sum(H,S,(-1.d0, 0.d0),Ec,ESH_tot)    

    call allocate_blk_dns(ESH,nbl)
    call zcsr2blk_sod(ESH_tot,ESH,str%mat_PL_start)
    call destroy(ESH_tot)

    !Inclusion of the contact Self-Energies to the relevant blocks
    do i=1,ncont
      ESH(str%cblk(i),str%cblk(i))%val = ESH(str%cblk(i),str%cblk(i))%val-SelfEneR(i)%val
    enddo

    call allocate_gsm_dns(gsmr,nbl)

    call Make_gsmr_mem_dns(ESH,nbl,2)

    call allocate_blk_dns(Gr,nbl)

    call Make_Gr_mem_dns(ESH,1)
    call Make_Gr_mem_dns(ESH,2,nbl)

    !Computation of transmission(s) between contacts ni(:) -> nf(:)
    do icpl=1,size_ni

      nit=ni(icpl)
      nft=nf(icpl)

      if (ncont.eq.2) then
        call trasmission_dns(nit,nft,ESH,SelfEneR,str%cblk,tun) 
      else
        call trasmission_old(nit,nft,ESH,SelfEneR,str%cblk,tun)
      endif

      TUN_MAT(icpl) = tun 

    enddo

    !Deallocate energy-dependent matrices
    call destroy_gsm(gsmr)
    call deallocate_gsm_dns(gsmr)
    !Distruzione dei blocchi fuori-diagonale
    do i=2,nbl
      call destroy(Gr(i-1,i))
      call destroy(Gr(i,i-1))
    enddo
    call destroy_ESH(ESH)
    call deallocate_blk_dns(ESH)

    call create(Grm,str%mat_PL_start(nbl+1)-1,str%mat_PL_start(nbl+1)-1,0)

    Grm%rowpnt(:)=1

    do i=1,nbl
      call create(GrCSR,Gr(i,i)%nrow,Gr(i,i)%ncol,Gr(i,i)%nrow*Gr(i,i)%ncol)
      call dns2csr(Gr(i,i),GrCSR)
      !Concatena direttamente la parte immaginaria per il calcolo della DOS
      zc=(-1.d0,0.d0)/pi

      call concat(Grm,zc,GrCSR,Im,str%mat_PL_start(i),str%mat_PL_start(i))
      call destroy(Gr(i,i))
      call destroy(GrCSR)
    enddo

    call deallocate_blk_dns(Gr)

    !Compute LDOS on the specified intervals
    if (nLDOS.gt.0) then
      call log_allocate(diag, Grm%nrow)
      call getdiag(Grm,diag)
      do iLDOS=1,nLDOS
        do i = 1, size(LDOS(iLDOS)%indexes)
          i2 = LDOS(iLDOS)%indexes(i)
          if (i2 .le. str%central_dim) then
            LEDOS(iLDOS) = LEDOS(iLDOS) + diag(i2)  
          end if
        end do
      enddo
      call log_deallocate(diag)
    endif

    call destroy(Grm)

  end subroutine tun_and_dos


  !---------------------------------------------------


  subroutine allocate_gsm_dns(gsm,nbl)
    type(z_DNS), dimension(:), allocatable :: gsm 
    integer :: nbl, ierr

    allocate(gsm(nbl),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate gsm'

  end subroutine allocate_gsm_dns

  !---------------------------------------------------

  subroutine allocate_blk_dns(blkM,nbl)
    type(z_DNS), dimension(:,:), allocatable :: blkM 
    integer :: nbl, ierr

    allocate(blkM(nbl,nbl),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate block-Matrix'

  end subroutine allocate_blk_dns

  !---------------------------------------------------

  subroutine deallocate_gsm_dns(gsm)
    type(z_DNS), dimension(:), allocatable :: gsm 
    integer :: ierr

    deallocate(gsm,stat=ierr)
    if (ierr.ne.0) stop 'DEALLOCATION ERROR: could not deallocate gsmr'

  end subroutine deallocate_gsm_dns

  !---------------------------------------------------

  subroutine deallocate_blk_dns(blkM)
    type(z_DNS), dimension(:,:), allocatable :: blkM 
    integer :: ierr

    deallocate(blkM,stat=ierr)
    if (ierr.ne.0) stop 'DEALLOCATION ERROR: could not deallocate block-Matrix'

  end subroutine deallocate_blk_dns

  !---------------------------------------------------






END MODULE iterative_dns




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
!!$    !Diagonal, Subdiagonal and Super953diagonal blocks
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

