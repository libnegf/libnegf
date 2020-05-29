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


 module negf_iterative_ph

  use negf_ln_precision
  use negf_ln_constants, only : pi
  use negf_ln_allocation
  use negf_mat_def
  use negf_sparsekit_drv
  USE inversions
  use negf_ln_structure, only : TStruct_Info
  use negf_lib_param, only : MAXNCONT, Tnegf, intarray
  use negf_outmatrix, only : outmat_c, inmat_c, direct_out_c, direct_in_c 
  use negf_clock
  !USE transform

  IMPLICIT NONE
  private

!  public :: calculate_transmissions
!  public :: calculate_transmissions_and_dos

  public :: calls_Dr_ph
  public :: calls_Dn_ph

!  public :: calculate_Gn_neq_components
!  public :: calls_neq_ph
!
!  public :: sigma_ph_n
!  public :: sigma_ph_p
!  public :: sigma_ph_r
!  public :: sigma_ph_r_z
!  public :: check_sigma_ph_r
!  public :: check_Gl_Gr
!
!  public :: csr2blkdns
!  public :: rebuild_dns
!  public :: calculate_gsmr_blocks
!  public :: calculate_gsml_blocks
!  public :: calculate_Gr_tridiag_blocks
!  public :: calculate_Gr_column_blocks
!  public :: calculate_Gn_tridiag_blocks
!  public :: calculate_Gr_outer
!  public :: calculate_Gn_outer
!  !public :: Outer_A_mem_dns
!
!  public :: complete_sigma_ph_r
!
  public :: create_scratch
  public :: destroy_scratch

  LOGICAL, PARAMETER :: debug=.false. 

  TYPE(z_DNS), DIMENSION(:), ALLOCATABLE :: gsmr
  TYPE(z_DNS), DIMENSION(:), ALLOCATABLE :: gsml
  TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Gr

  TYPE(z_DNS3), DIMENSION(:,:), ALLOCATABLE, SAVE :: DDn
  TYPE(z_DNS3), DIMENSION(:,:), ALLOCATABLE, SAVE :: DDp
  TYPE(z_DNS3), DIMENSION(:,:), ALLOCATABLE, SAVE :: DDr
  TYPE(z_DNS3), DIMENSION(:,:), ALLOCATABLE, SAVE :: Pin
  TYPE(z_DNS3), DIMENSION(:,:), ALLOCATABLE, SAVE :: Pip
  TYPE(z_DNS3), DIMENSION(:,:), ALLOCATABLE, SAVE :: Pir
  LOGICAL, PARAMETER :: memory = .true.

CONTAINS

  !****************************************************************************
  ! 
  ! Driver for computing Equilibrium Green Retarded + Sigma_c + Sigma_ph 
  !
  !****************************************************************************

  SUBROUTINE calls_Dr_ph(pnegf,E,SelfEneR,struct,comp_col,A)

    !****************************************************************************
    !
    !Input
    !pnegf:    negf data container
    !E:        Energy point
    !SelfEneR: matrices array containing contacts Self Energy
    !struct:   structure container to get nPls, nConts, indblk
    !comp_col: compute column GF
    !A:        optional, gives back the whole Gr masked by the overlap  
    !
    !*****************************************************************************

    IMPLICIT NONE
 
    !In/Out
    TYPE(Tnegf), intent(in) :: pnegf
    COMPLEX(dp), intent(in) :: E
    TYPE(z_DNS), DIMENSION(:), intent(in) :: SelfEneR
    TYPE(Tstruct_info), intent(in) :: struct
    logical, intent(in) :: comp_col 
    TYPE(z_CSR), intent(out), optional :: A

    !Work
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: ESH
    TYPE(z_CSR) :: ESH_csr 
    INTEGER :: i, m, cb, ierr, nbl, rbl, lbl, ncont, iE, np
    INTEGER, DIMENSION(:), POINTER :: cblk, indblk

    nbl = struct%num_PLs
    ncont = struct%num_conts
    cblk => struct%cblk
    indblk => struct%mat_PL_start
    iE = pnegf%iE
    np = pnegf%local_en_points

    ! Take CSR H,S and build ES-H in dense blocks
    CALL prealloc_sum(pnegf%H,pnegf%S,(-1.d0, 0.d0),E,ESH_csr)
    call allocate_blk_dns(ESH,nbl)
    CALL csr2blkdns(ESH_csr,ESH,indblk)
    CALL destroy(ESH_csr)
    
    ! ADDING CONTACT GF
    DO i=1,ncont
       cb = cblk(i)
       ESH(cb,cb)%val = ESH(cb,cb)%val-SelfEneR(i)%val
    ENDDO

    ! ADDING PH SELF ENERGY
    if (pnegf%phph%include_phph .and. pnegf%phph%scba_iter.gt.0) then 
      call add_sigma_ph_r(pnegf, ESH)
    end if
 
    ! BLOCK-ITERATIVE ALGORITHM (Build Gr tridiagonal blocks)
    call allocate_blk_dns(Gr,nbl)
    call allocate_gsm(gsmr,nbl)
    rbl = minval(cblk(1:ncont)) + 1  
    lbl = maxval(cblk(1:ncont)) - 1

    ! -------------------------------------------------------------
    ! MESSY PART TO COMPUTE THE COLUMNS of Gr AT THE CONTACTS
    ! Needed for:  Gn = Gr(i,cb) Gamma(cb,cb) Ga(cb,j) 
    ! NOTE the behaviour of Make_Gr_mem:
    !      (ESH,i)   computes block i,i by inversion
    !      (ESH,i,i) computes block i,i by iteration
    ! -------------------------------------------------------------
     
    IF (comp_col) THEN
      call allocate_gsm(gsml,nbl)

      ! Fix to a bug when there are 2PLs
      ! later Make_Gr tries to compute Gr(1,1) but needs gsmr(2,2)
      ! 
      IF (nbl.eq.2) then
        CALL calculate_gsmr_blocks(ESH,nbl,rbl-1)
        CALL calculate_gsml_blocks(ESH,1,lbl+1)    
      ELSE
        CALL calculate_gsmr_blocks(ESH,nbl,rbl)
        CALL calculate_gsml_blocks(ESH,1,lbl)    
      ENDIF

      ! 1. rbl>lbl  => lbl+1=rbl-1 => compute first Gr(rbl-1,rbl-1)
      ! 2. rbl<lbl  => lbl=rbl-2 has been computed
      ! Make_Gr does not compute if sbl>nbl or sbl<1
      CALL calculate_Gr_tridiag_blocks(ESH,rbl-1)
      CALL calculate_Gr_tridiag_blocks(ESH,rbl,nbl)
      CALL calculate_Gr_tridiag_blocks(ESH,rbl-2,1)
 
      ! build contact column green's for later use 
      DO i=1,ncont
         CALL calculate_Gr_column_blocks(ESH,cblk(i),indblk)
      ENDDO

      call destroy_gsm(gsml)
      call deallocate_gsm(gsml)

    ELSE
      
      CALL calculate_gsmr_blocks(ESH,nbl,2)
      CALL calculate_Gr_tridiag_blocks(ESH,1)
      CALL calculate_Gr_tridiag_blocks(ESH,2,nbl)

    ENDIF
 
    CALL destroy_ESH(ESH)
    CALL deallocate_blk_dns(ESH)
    call destroy_gsm(gsmr)
    call deallocate_gsm(gsmr)

    ! SAVE ON FILES/MEMORY (for phph).........................
    IF (pnegf%phph%include_phph) THEN
      DO i = 1, nbl
         !print*,'G_r ',minval(abs(Gr(i,i)%val)), maxval(abs(Gr(i,i)%val))
         call write_blkmat(Gr(i,i),pnegf%scratch_path,'G_r_',i,i,iE,np)
      ENDDO
      DO i = 1, nbl-1
         call write_blkmat(Gr(i,i+1),pnegf%scratch_path,'G_r_',i,i+1,iE,np)
         call write_blkmat(Gr(i+1,i),pnegf%scratch_path,'G_r_',i+1,i,iE,np)
      ENDDO
    END IF
   
    ! SAVES COLUMN Gr (needed for Gn =  Gr Gamma Ga) 
    IF (comp_col) THEN
      DO m=1,ncont
        cb = cblk(m)
        DO i = 1, cb-2
          call write_blkmat(Gr(i,cb),pnegf%scratch_path,'G_r_',i,cb,iE,np)
        ENDDO
        DO i = cb+2, nbl
          call write_blkmat(Gr(i,cb),pnegf%scratch_path,'G_r_',i,cb,iE,np)
        ENDDO
      ENDDO
    END IF
   
    !..........................................................
    if (present(A)) then
       call blk2csr(Gr,struct,pnegf%S,A)
    endif

    CALL destroy_blk(Gr)
    DEALLOCATE(Gr)

  END SUBROUTINE calls_Dr_ph

  !****************************************************************************
  !
  ! Driver for computing G_n = (+/-)iG< contributions including interactions
  !
  !    Sum   f_j(E) Gr Gam_j Ga +   Gr Sigma_ph< Ga
  !     j
  !
  ! NOTE: 
  !
  !****************************************************************************
  SUBROUTINE calls_Dn_ph(pnegf,E,SelfEneR,nB,struct,Gout)
    
    !****************************************************************************
    !
    !Input
    !pnegf:    negf data container
    !E:        Energy point
    !SelfEneR: matrices array containing contacts Self Energy
    !nB:       bose distribution of each contact
    !struct:   structure container to get nPls, nConts, indblk,
    !A:        optional, gives back the whole Gr masked by the overlap  
    !
    !*****************************************************************************

    IMPLICIT NONE

    !In/Out
    TYPE(Tnegf), intent(in) :: pnegf
    REAL(dp), intent(in)  :: E
    TYPE(z_DNS), DIMENSION(:), intent(in)  :: SelfEneR
    real(dp), DIMENSION(:), intent(in)  :: nB
    TYPE(Tstruct_info), intent(in)  :: struct
    TYPE(z_CSR), intent(inout), optional  :: Gout

    !Work
    COMPLEX(dp) :: Ec
    INTEGER :: i, k, ierr,ncont,nbl, lbl, rbl, ref, iE, np
    INTEGER, DIMENSION(:), POINTER :: cblk, indblk
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: ESH
    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Gn
    TYPE(z_CSR) :: ESH_csr, Gl
    LOGICAL :: mask(MAXNCONT)

    nbl = struct%num_PLs
    ncont = struct%num_conts
    indblk => struct%mat_PL_start
    cblk => struct%cblk
    ref = pnegf%refcont
    iE = pnegf%iE
    np = pnegf%local_en_points

    Ec=cmplx(E,0.d0,dp)


    DO i=1,ncont
       ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
    ENDDO

    call allocate_blk_dns(Gn,nbl)

    ! COMPUTE: Gn =  Sum_i [ Gr Gamma_i Ga nB_i ] 
    call init_blkmat(Gn,ESH)
    CALL calculate_Gn_tridiag_blocks(ESH,SelfEneR,nB,ref,struct,Gn)

    ! ADD THE SCATTERING SELF-ENERGY:  Gn = Gr Sigma_n Ga
    IF (pnegf%phph%include_phph .and. pnegf%phph%scba_iter.gt.0) THEN
       call calculate_Gn_tridiag_elph_contributions(pnegf,ESH,Gn)
    ENDIF

    ! SAVE tridiagonal blocks of G_n
    IF (pnegf%phph%include_phph) THEN
      DO i = 1, nbl
        !print*,'(G_n) G_n',minval(abs(Gn(i,i)%val)), maxval(abs(Gn(i,i)%val))
        call write_blkmat(Gn(i,i),pnegf%scratch_path,'G_n_',i,i,iE,np)
      ENDDO
      DO i = 1, nbl-1
        call write_blkmat(Gr(i,i+1),pnegf%scratch_path,'G_n_',i,i+1,iE,np)
        call write_blkmat(Gr(i+1,i),pnegf%scratch_path,'G_n_',i+1,i,iE,np)
      ENDDO
    ENDIF

    if (present(Gout)) then
      call blk2csr(Gn,struct,pnegf%S,Gout)
    endif
 
    CALL destroy_blk(Gn)
    DEALLOCATE(Gn)

    CALL destroy_blk(Gr)
    DEALLOCATE(Gr)

    CALL destroy_ESH(ESH)
    DEALLOCATE(ESH)

  END SUBROUTINE calls_Dn_ph

  !***********************************************************************
  !
  !  FEW utility subrutines to allocate/deallocate block-dense matrices
  ! 
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
    enddo

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
  !  Divides a sparse matrix A_csr into an array of dense matrices, 
  !  A(nbl,nbl) 
  !
  !**********************************************************************

  SUBROUTINE csr2blkdns(A_csr,A,indblk)

    !**********************************************************************
    !Input:
    !ESH_csr: sparse matrix ES-H related to device
    !
    !Output:
    !ESH(nbl,nbl): dense matrix array -> single matrices allocated 
    !              internally, array ESH(nbl,nbl) allocated externally
    !**********************************************************************

    IMPLICIT NONE 

    INTEGER :: i
    TYPE(z_CSR) :: A_csr
    INTEGER :: nbl
    TYPE(z_DNS), DIMENSION(:,:) :: A
    INTEGER, DIMENSION(:) :: indblk

    nbl = size(A,1)

    DO i=1,nbl
       CALL extract(A_csr,indblk(i),indblk(i+1)-1,indblk(i),indblk(i+1)-1,A(i,i))
    ENDDO

    DO i=2,nbl
       CALL extract(A_csr,indblk(i-1),indblk(i)-1,indblk(i),indblk(i+1)-1,A(i-1,i))
       CALL extract(A_csr,indblk(i),indblk(i+1)-1,indblk(i-1),indblk(i)-1,A(i,i-1))
    ENDDO

  END SUBROUTINE csr2blkdns

  !***********************************************************************
  !
  !  Reloads the retarded elph self-energy and add it to ES-H
  !
  !***********************************************************************
  SUBROUTINE add_sigma_ph_r(pnegf, ESH)
     TYPE(Tnegf), intent(in) :: pnegf
     TYPE(z_DNS), DIMENSION(:,:), intent(inout) :: ESH

     TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_ph_r
     INTEGER :: n, nbl, nrow, ierr

     nbl = pnegf%str%num_PLs
    
     ALLOCATE(Sigma_ph_r(nbl,nbl),stat=ierr)
     IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_r'


     DO n = 1, nbl
        call create(Sigma_ph_r(n,n), ESH(n,n)%nrow, ESH(n,n)%nrow)
        Sigma_ph_r(n,n)%val = (0.0_dp, 0.0_dp)
        call read_blkmat(Sigma_ph_r(n,n),pnegf%scratch_path,'Sigma_ph_r_',n,n,pnegf%iE)
        ESH(n,n)%val = ESH(n,n)%val - Sigma_ph_r(n,n)%val
        call destroy(Sigma_ph_r(n,n))
     END DO

     if (.not. pnegf%phph%diagonal) THEN
       DO n = 1, nbl-1       
          call create(Sigma_ph_r(n,n+1), ESH(n,n+1)%nrow, ESH(n,n+1)%ncol)
          Sigma_ph_r(n,n+1)%val = (0.0_dp, 0.0_dp)
          call read_blkmat(Sigma_ph_r(n,n+1),pnegf%scratch_path,'Sigma_ph_r_',n,n+1,pnegf%iE)
          ESH(n,n+1)%val = ESH(n,n+1)%val - Sigma_ph_r(n,n+1)%val
          call destroy(Sigma_ph_r(n,n+1))
          
          call create(Sigma_ph_r(n+1,n), ESH(n+1,n)%nrow, ESH(n+1,n)%ncol)
          Sigma_ph_r(n+1,n)%val = (0.0_dp, 0.0_dp)
          call read_blkmat(Sigma_ph_r(n+1,n),pnegf%scratch_path,'Sigma_ph_r_',n+1,n,pnegf%iE)
          ESH(n+1,n)%val = ESH(n+1,n)%val - Sigma_ph_r(n+1,n)%val
          call destroy(Sigma_ph_r(n+1,n))
       ENDDO
     ENDIF

     DEALLOCATE(Sigma_ph_r)

  END SUBROUTINE add_sigma_ph_r 


  !***********************************************************************
  !
  !  g_small right (gsmr) calculation - write on memory
  !
  !***********************************************************************

  SUBROUTINE calculate_gsmr_blocks(ESH,sbl,ebl)

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


  END SUBROUTINE calculate_gsmr_blocks

 


  !***********************************************************************
  !
  !  g_small left (gsml) calculation - write on memory
  !
  !***********************************************************************

  SUBROUTINE calculate_gsml_blocks(ESH,sbl,ebl)

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

  END SUBROUTINE calculate_gsml_blocks





  !***********************************************************************
  !
  !  Diagonal, Subdiagonal, Superdiagonal blocks of Green Retarded 
  !  Gr(nbl,nbl) - writing on memory
  !
  !*********************************************************************** 

  SUBROUTINE calculate_Gr_tridiag_blocks(ESH,sbl,ebl)

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

  END SUBROUTINE calculate_Gr_tridiag_blocks

  !**************************************************************************
  !
  !  Calculate Green Retarded column "n" - writing on memory 
  !
  !**************************************************************************

  SUBROUTINE calculate_Gr_column_blocks(ESH,n,indblk)

    !***********************************************************************
    !Input:
    !ESH: sparse matrices array ESH(nbl,nbl)
    !n: n umber of column to be calculated
    !
    !global variables needed:  
    !Gr diagonal, subadiagonal and superdiagonal, 
    !gsmr(:) for downgoing and gsml(:) for upgoing 
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
       ENDDO

    ENDIF

    if (debug) then
       WRITE(*,*) '********************'
       WRITE(*,*) 'Make_Grcol_mem done column',n
       WRITE(*,*) '********************'
    endif

  END SUBROUTINE calculate_Gr_column_blocks

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

  SUBROUTINE calculate_Gn_tridiag_blocks(ESH,SelfEneR,nX,ref,struct,Gn)

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
    REAL(dp), DIMENSION(:), intent(in) :: nX !bosons or fermions 
    INTEGER, intent(in) :: ref
 
    !Work
    Type(z_DNS) :: Gam
    TYPE(z_DNS) :: work1,Ga
    INTEGER :: i,j,cb
    INTEGER :: ncont, nbl
    INTEGER, DIMENSION(:), POINTER :: cblk
    COMPLEX(dp) :: nXdiff

    ncont = struct%num_conts    
    nbl = struct%num_PLs
    cblk => struct%cblk
    
    !*******************************************
    ! Contact Iteration
    !*******************************************
    DO j=1,ncont
  
       ! NOTE: this soubroutine uses prealloc_mult that performs 
       ! C = C + A*B 
       IF (j.NE.ref .AND. ABS(nX(j)-nX(ref)).GT.EPS) THEN
 
          cb=cblk(j) ! block corresponding to contact j

          CALL zspectral(SelfEneR(j),SelfEneR(j),0,Gam)

          nXdiff = cmplx(nX(j)-nX(ref),0.d0,dp)
          ! Computation of Gl(1,1) = Gr(1,cb) Gam(cb) Ga(cb,1)
          if (allocated(Gr(1,cb)%val)) then
            CALL prealloc_mult(Gr(1,cb),Gam,work1)    
            CALL zdagger(Gr(1,cb),Ga)
            CALL prealloc_mult(work1,Ga,nXdiff,Gn(1,1))
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
                CALL prealloc_mult(work1,Ga,nXdiff,Gn(i-1,i))
                CALL destroy(work1)
                
                CALL prealloc_mult(Gr(i,cb),Gam,work1)
                CALL prealloc_mult(work1,Ga,nXdiff,Gn(i,i))

                CALL destroy(work1, Ga)
               
                CALL prealloc_mult(Gr(i,cb),Gam,work1)
                CALL zdagger(Gr(i-1,cb),Ga)
                CALL prealloc_mult(work1,Ga,nXdiff,Gn(i,i-1))

                CALL destroy(work1, Ga)
             else
                Gn(i-1,i)%val=(0.d0,0.d0)
                Gn(i,i-1)%val=(0.d0,0.d0)
             endif      

          ENDDO
                  
          call destroy(Gam)

      ENDIF 

    ENDDO

  END SUBROUTINE calculate_Gn_tridiag_blocks


  !****************************************************************************
  !
  ! Calculate G_n contributions due to el-ph
  ! Writing on memory
  ! The subroutine assumes a block-tridiagonal structure for Sigma_ph_n
  ! Computes tridiagonal structure for Gn(i,j): 
  !           --       --
  ! Gn(i,j) = >        >       Gr(i,k) Sigma_ph(k,k+s) Ga(k+s,j)
  !           --       --
  !          k=1,nbl s=-1,0,1
  !****************************************************************************
  SUBROUTINE calculate_Gn_tridiag_elph_contributions(pnegf,ESH,Gn)

    TYPE(Tnegf), intent(in) :: pnegf
    TYPE(z_DNS), DIMENSION(:,:), intent(in) :: ESH
    TYPE(z_DNS), DIMENSION(:,:), intent(inout) :: Gn

    Type(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_ph_n
    Type(z_DNS) :: Ga, work1, work2
    INTEGER :: n, k, s, nbl, nrow, ncol, ierr

    nbl = pnegf%str%num_PLs
    ALLOCATE(Sigma_ph_n(nbl,nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_n'

    DO n = 1, nbl
      nrow = ESH(n,n)%nrow
      call create(Sigma_ph_n(n,n), nrow, nrow)
      Sigma_ph_n(n,n)%val = (0.0_dp, 0.0_dp)
      call read_blkmat(Sigma_ph_n(n,n),pnegf%scratch_path,'Sigma_ph_n_',n,n,pnegf%iE)
    END DO
    DO n = 1, nbl-1 
      nrow = ESH(n,n+1)%nrow
      ncol = ESH(n,n+1)%ncol
      call create(Sigma_ph_n(n,n+1), nrow, ncol)
      call create(Sigma_ph_n(n+1,n), ncol, nrow)
      Sigma_ph_n(n,n+1)%val = (0.0_dp, 0.0_dp)
      Sigma_ph_n(n+1,n)%val = (0.0_dp, 0.0_dp)
      call read_blkmat(Sigma_ph_n(n,n+1),pnegf%scratch_path,'Sigma_ph_n_',n,n+1,pnegf%iE)
      call read_blkmat(Sigma_ph_n(n+1,n),pnegf%scratch_path,'Sigma_ph_n_',n+1,n,pnegf%iE)
    END DO

    DO n = 1, nbl
      DO k = 1, nbl
        DO s = -1,+1 
          if (k+s.ge.1 .and. k+s.le.nbl) then
            if (Gr(n,k+s)%nrow.gt.0) then
              CALL prealloc_mult(Gr(n,k+s), Sigma_ph_n(k,k+s), work1)
              CALL zdagger(Gr(n,k+s),Ga)
              CALL prealloc_mult(work1, Ga, work2)
              Gn(n,n)%val = Gn(n,n)%val + work2%val
              call destroy(work2,Ga)
            endif
          endif
          if (k+s.ge.1 .and. k+s.le.nbl .and. n+1.le.nbl) then
            if (Gr(n+1,k+s)%nrow.gt.0) then
               CALL zdagger(Gr(n+1,k+s),Ga)
               CALL prealloc_mult(work1, Ga, work2)
               Gn(n,n+1)%val = Gn(n,n+1)%val + work2%val
               call destroy(work1,work2,Ga)
            endif
          endif
        END DO
      END DO    
      Gn(n+1,n)%val =  Gn(n+1,n)%val + conjg(transpose(Gn(n,n+1)%val))
    END DO
    
    DO n = 1, nbl-1
      CALL destroy(Sigma_ph_n(n,n))
      CALL destroy(Sigma_ph_n(n,n+1))
      CALL destroy(Sigma_ph_n(n+1,n))
    END DO
    CALL destroy(Sigma_ph_n(nbl,nbl))

    DEALLOCATE(Sigma_ph_n)

  END SUBROUTINE calculate_Gn_tridiag_elph_contributions
 

  !****************************************************************************
  !Subroutine used to convert blk-dense to csr format, masked by the matrix P
  !If only one block is present, concatenation is not needed and it's implemented in a
  !more trivial way
  !****************************************************************************
  SUBROUTINE blk2csr(G,struct,P,Gcsr)

    TYPE(z_DNS), DIMENSION(:,:) :: G
    TYPE(Tstruct_info), intent(in) :: struct
    TYPE(z_CSR), intent(in) :: P
    TYPE(z_CSR), intent(out) :: Gcsr

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

 
!  !****************************************************************************
!  !
!  ! Calculate G_p=iG> contributions due to el-ph
!  ! Writing on memory
!  !
!  !****************************************************************************
!  SUBROUTINE Make_Gp_ph(pnegf,ESH,iter,Gp)
!
!    TYPE(Tnegf) :: pnegf
!    TYPE(z_DNS), DIMENSION(:,:) :: ESH, Gp
!    INTEGER :: iter
!
!    Type(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_ph_p
!    Type(z_DNS) :: Ga, work1, work2 
!    INTEGER :: n, k, nbl, nrow, ierr
!
!    nbl = pnegf%str%num_PLs
!    ALLOCATE(Sigma_ph_p(nbl,nbl),stat=ierr)
!    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_n'
!
!
!    DO n = 1, nbl
!
!      nrow = ESH(n,n)%nrow
!
!      call create(Sigma_ph_p(n,n), nrow, nrow)
!
!      Sigma_ph_p(n,n)%val = (0.0_dp, 0.0_dp)
!      if (iter .gt. 0) then
!         call read_blkmat(Sigma_ph_p(n,n),pnegf%scratch_path,'Sigma_ph_p_',n,n,pnegf%iE)
!      else 
!         call write_blkmat(Sigma_ph_p(n,n),pnegf%scratch_path,'Sigma_ph_p_',n,n,pnegf%iE)
!      endif
!    
!    END DO
!
!    DO n = 1, nbl
!
!      DO k = 1, nbl
!
!         if (Gr(n,k)%nrow.gt.0) then
!            CALL zdagger(Gr(n,k),Ga)
!            CALL prealloc_mult(Gr(n,k), Sigma_ph_p(k,k), work1)
!            CALL prealloc_mult(work1, Ga, work2)
!            Gp(n,n)%val = Gp(n,n)%val + work2%val
!            call destroy(work1, work2, Ga)
!         endif
!
!      END DO    
!    
!    END DO
!    
!    DO n = 1, nbl
!      CALL destroy(Sigma_ph_p(n,n))
!    END DO
!
!    DEALLOCATE(Sigma_ph_p)
!
!  END SUBROUTINE Make_Gp_ph
!  ! ******************************************************************************
!  ! Computes Sigma_ph_n and save it file
!  ! ******************************************************************************
!  SUBROUTINE Sigma_ph_n(pnegf,Epnt)
!
!    TYPE(Tnegf) :: pnegf
!    REAL(dp), DIMENSION(:),allocatable :: Epnt
!
!    !Local variables
!    TYPE(z_DNS) :: work1,work2
!    TYPE(z_DNS) :: Mq_mat
!    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_n, G_n_interP, G_n_interN
!
!    REAL(dp), DIMENSION(:), POINTER :: Wq
!    REAL(dp), DIMENSION(:), POINTER :: Nq
!    REAL(dp), DIMENSION(:), POINTER :: Mq
!    LOGICAL, DIMENSION(:), POINTER :: selmodes
!
!    INTEGER :: i, m, iE, i1,i2
!    INTEGER :: nummodes, numselmodes, nbl
!    INTEGER, DIMENSION(:), pointer :: indblk
!    REAL(dp) :: E1, E2, En
!
!    nbl = pnegf%str%num_PLs
!    indblk => pnegf%str%mat_PL_start
!
!    selmodes => pnegf%elph%selmodes
!    Wq => pnegf%elph%Wq
!    Mq => pnegf%elph%Mq
!    Nq => pnegf%elph%Nq
!    numselmodes = pnegf%elph%numselmodes
!    nummodes = pnegf%elph%nummodes
!
!    iE = pnegf%iE
!
!    allocate(G_n_interP(nbl,nbl))
!    allocate(G_n_interN(nbl,nbl))
!    allocate(Sigma_n(nbl,nbl))
!
!    do i = 1, nbl
!      m = indblk(i+1)-indblk(i)
!      call create(Sigma_n(i,i), m, m)
!      Sigma_n(i,i)%val = (0.d0, 0.d0)
!      call create(G_n_interP(i,i), m, m)
!      call create(G_n_interN(i,i), m, m)
!    enddo
!    
!    do m = 1 , nummodes
!
!      if (.not.selmodes(m)) cycle
!      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
!      call search_points(pnegf,Wq(m),Epnt,i1,i2,E1,E2)
!      
!      En = real(pnegf%Epnt)+Wq(m)
!      
!      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
!                             'G_n_', G_n_interP)
!
!      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
!      call search_points(pnegf,-Wq(m),Epnt,i1,i2,E1,E2)
!      
!      En = real(pnegf%Epnt)-Wq(m)
!      
!      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
!                              'G_n_', G_n_interN)
!
!
!      do i = 1, nbl
!
!        i1 =  G_n_interN(i,i)%nrow
!
!        call create(work1, i1, i1)       
!
!        work1%val = (Nq(m)+1.0_dp)*G_n_interP(i,i)%val + Nq(m)* G_n_interN(i,i)%val
!
!        call create_id(Mq_mat,i1,Mq(m))
!
!        call prealloc_mult( Mq_mat, work1, work2)
!
!        call destroy(work1)
!
!        call prealloc_mult( work2, Mq_mat, work1)
!
!        call destroy(work2)
!
!        Sigma_n(i,i)%val = Sigma_n(i,i)%val + work1%val 
!
!        call destroy(work1, Mq_mat)
!
!     end do
!
!   end do
!
!   !Throw away all non diagonal parts
!   !do i = 1, nbl
!   !   i1 = Sigma_n(i,i)%nrow
!   !   do m = 1, i1 
!   !      Sigma_n(i,i)%val(1:m-1,m) = (0.0_dp,0.0_dp) 
!   !      Sigma_n(i,i)%val(m+1:i1,m) = (0.0_dp,0.0_dp) 
!   !   enddo
!   !enddo
!   !
!
!   do i = 1, nbl
!     call write_blkmat(Sigma_n(i,i),pnegf%scratch_path,'Sigma_ph_n_',i,i,iE)
!     call destroy(Sigma_n(i,i))
!     call destroy(G_n_interN(i,i))
!     call destroy(G_n_interP(i,i))
!   enddo
!
!   deallocate(Sigma_n,G_n_interN, G_n_interP)
!
! END SUBROUTINE Sigma_ph_n
!
!
!  ! ******************************************************************************
!  ! Computes Sigma_ph> and save it file
!  ! ******************************************************************************
!  SUBROUTINE Sigma_ph_p(pnegf,Epnt)
!
!    TYPE(Tnegf) :: pnegf
!    REAL(dp), DIMENSION(:),allocatable :: Epnt
!
!    !Local variables
!    TYPE(z_DNS) :: work1,work2
!    TYPE(z_DNS) :: Mq_mat
!    TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_p, G_p_interP, G_p_interN
!
!    REAL(dp), DIMENSION(:), POINTER :: Wq
!    REAL(dp), DIMENSION(:), POINTER :: Nq
!    REAL(dp), DIMENSION(:), POINTER :: Mq
!    LOGICAL, DIMENSION(:), POINTER :: selmodes
!
!    INTEGER :: i, m, iE, i1,i2
!    INTEGER :: nummodes, numselmodes, nbl
!    INTEGER, DIMENSION(:), pointer :: indblk
!    REAL(dp) :: E1, E2, En
!
!    nbl = pnegf%str%num_PLs
!    indblk => pnegf%str%mat_PL_start
!
!    selmodes => pnegf%elph%selmodes
!    Wq => pnegf%elph%Wq
!    Mq => pnegf%elph%Mq
!    Nq => pnegf%elph%Nq
!    numselmodes = pnegf%elph%numselmodes
!    nummodes = pnegf%elph%nummodes
!
!    iE = pnegf%iE
!
!    allocate(G_p_interP(nbl,nbl))
!    allocate(G_p_interN(nbl,nbl))
!    allocate(Sigma_p(nbl,nbl))
!
!    do i = 1, nbl
!      m = indblk(i+1)-indblk(i)
!      call create(Sigma_p(i,i), m, m)
!      Sigma_p(i,i)%val = (0.d0, 0.d0)
!      call create(G_p_interP(i,i), m, m)
!      call create(G_p_interN(i,i), m, m)
!    enddo
!    
!    do m = 1 , nummodes
!
!      if (.not.selmodes(m)) cycle
!      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
!      call search_points(pnegf,Wq(m),Epnt,i1,i2,E1,E2)
!    
!      En = real(pnegf%Epnt)+Wq(m)
!      
!      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
!                             'G_p_', G_p_interP)
!      
!      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
!      call search_points(pnegf,-Wq(m),Epnt,i1,i2,E1,E2)
!      En = real(pnegf%Epnt)-Wq(m)
!      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
!                              'G_p_', G_p_interN)
!
!      do i = 1, nbl
!
!        !print*
!        !print*,'(sigma_r) G_p+',maxval(abs(G_p_interP(i,i)%val))
!        !print*,'(sigma_r) G_p-',maxval(abs(G_p_interN(i,i)%val))
!
!        i1 =  G_p_interN(i,i)%nrow
!
!        call create(work1, i1, i1)       
!
!        work1%val = (Nq(m)+1.0_dp)*G_p_interN(i,i)%val + Nq(m)* G_p_interP(i,i)%val
!
!        call create_id(Mq_mat,i1,Mq(m))
!
!        call prealloc_mult( Mq_mat, work1, work2)
!
!        call destroy(work1)
!
!        call prealloc_mult ( work2, Mq_mat, work1)
!
!        call destroy(work2)
!
!        Sigma_p(i,i)%val = Sigma_p(i,i)%val + work1%val 
!
!        call destroy(work1, Mq_mat)
!
!     end do
!
!   end do
!   
!   ! throw away all non diagonal parts
!   !do i = 1, nbl
!   !   i1 = Sigma_p(i,i)%nrow
!   !   do m = 1, i1 
!   !      Sigma_p(i,i)%val(1:m-1,m) = (0.0_dp,0.0_dp) 
!   !      Sigma_p(i,i)%val(m+1:i1,m) = (0.0_dp,0.0_dp) 
!   !   enddo
!   !enddo
!   !
!
!   do i = 1, nbl
!     call write_blkmat(Sigma_p(i,i),pnegf%scratch_path,'Sigma_ph_p_',i,i,iE)
!     call destroy(Sigma_p(i,i))
!     call destroy(G_p_interN(i,i))
!     call destroy(G_p_interP(i,i))
!   enddo
!
!  deallocate(Sigma_p,G_p_interN, G_p_interP)
!
! END SUBROUTINE Sigma_ph_p
!
! !----------------------------------------------------------------------------------
! SUBROUTINE Sigma_ph_r(pnegf,Epnt)
!
!  TYPE(Tnegf) :: pnegf
!  REAL(dp), DIMENSION(:),allocatable :: Epnt
!
!  !Local variables
!  TYPE(z_DNS) :: work1,work2
!  TYPE(z_DNS) :: Mq_mat
!  TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_r, G_r_interP, G_r_interN
!  TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: G_n_interP, G_n_interN
!
!  REAL(dp), DIMENSION(:), POINTER :: Wq
!  REAL(dp), DIMENSION(:), POINTER :: Nq
!  REAL(dp), DIMENSION(:), POINTER :: Mq
!  LOGICAL, DIMENSION(:), POINTER :: selmodes
!  INTEGER :: i, m, iE, i1, i2
!  INTEGER :: nummodes, numselmodes, nbl
!  INTEGER, DIMENSION(:), pointer :: indblk
!  REAL(dp) :: E1, E2, En
!
!  nbl = pnegf%str%num_PLs
!  indblk => pnegf%str%mat_PL_start
!
!
!  selmodes => pnegf%elph%selmodes
!  Mq => pnegf%elph%Mq
!  Wq => pnegf%elph%Wq
!  Nq => pnegf%elph%Nq
!  nummodes = pnegf%elph%nummodes
!  numselmodes = pnegf%elph%numselmodes
!
!  iE=pnegf%iE
!
!  allocate(G_r_interP(nbl,nbl))
!  allocate(G_r_interN(nbl,nbl))
!  allocate(G_n_interP(nbl,nbl))
!  allocate(G_n_interN(nbl,nbl))
!  allocate(Sigma_r(nbl,nbl))
!  !allocate(G_r(nbl,nbl))
!
!    
!      
!  do i = 1, nbl
!      m = indblk(i+1)-indblk(i)
!      call create(Sigma_r(i,i), m, m)
!      Sigma_r(i,i)%val = (0.d0, 0.d0)
!      call create(G_r_interP(i,i), m, m)
!      call create(G_r_interN(i,i), m, m)
!      call create(G_n_interP(i,i), m, m)
!      call create(G_n_interN(i,i), m, m)
!      !call create(G_r(i,i), m, m)
!      !call read_blkmat(G_r(i,i),pnegf%scratch_path,'G_r_',i,i,iE)
!  enddo
!
!  do m = 1 , nummodes
!
!      if (.not.selmodes(m)) cycle
!
!      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
!      call search_points(pnegf,Wq(m),Epnt,i1,i2,E1,E2)
!      En = real(pnegf%Epnt)+Wq(m)
!      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
!                             'G_n_', G_n_interP)
!      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
!                             'G_r_', G_r_interP)
!
!      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
!      call search_points(pnegf,-Wq(m),Epnt,i1,i2,E1,E2)
!      En = real(pnegf%Epnt)-Wq(m)
!      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
!                             'G_n_', G_n_interN)
!      call interpolation(i1, i2, E1, E2, En, pnegf%scratch_path, &
!                             'G_r_', G_r_interN)
!
!      do i = 1, nbl
!         !print*
!         !print*,'(sigma_r) G_r+',maxval(abs(G_r_interP(i,i)%val))
!         !print*,'(sigma_r) G_r-',maxval(abs(G_r_interN(i,i)%val))
!         
!         !print*,'(sigma_r) G_n+',maxval(abs(G_n_interP(i,i)%val))
!         !print*,'(sigma_r) G_n-',maxval(abs(G_n_interN(i,i)%val))
!                         
!         i1 = Sigma_r(i,i)%nrow
!
!         call create(work1,i1,i1)
! 
!         work1%val = (0.d0,0.d0)
!
!         if (pnegf%elph%selfene_gr) then
!             !print*,'SelfEneR_Gr'    
!             ! Via I: should be exact for E >> Ef_max
!             work1%val = (Nq(m)+1.0_dp)*G_r_interN(i,i)%val + Nq(m)*G_r_interP(i,i)%val
!             ! Via II: should be exact for E << Ef_min
!             !work1%val = (Nq(m)+1.0_dp)*G_r_interP(i,i)%val + Nq(m)*G_r_interN(i,i)%val
!             ! Via III: should work as a compromise
!             !work1%val = Nq(m)*(G_r_interP(i,i)%val + G_r_interN(i,i)%val) + G_r(i,i)%val
!         endif   
!            
!         if (pnegf%elph%selfene_gless) then
!             !print*,'SelfEneR_G<'    
!             ! should be 1/2 [G<- - G<+] == i/2 [Gn- - Gn+]    
!             work1%val = work1%val + &
!                    (0.0_dp,0.5_dp) * (G_n_interN(i,i)%val - G_n_interP(i,i)%val)
!         endif
!
!         if (pnegf%elph%selfene_gless .or. pnegf%elph%selfene_gr) then
!
!             call create_id(Mq_mat,i1,Mq(m))
!                      
!             call prealloc_mult(Mq_mat, work1, work2)
!
!             call destroy(work1)
!
!             call prealloc_mult(work2, Mq_mat, work1)
!
!             call destroy(work2)
!
!             Sigma_r(i,i)%val = Sigma_r(i,i)%val + work1%val 
!
!             call destroy(work1, Mq_mat)
!
!          else
!
!             Sigma_r(i,i)%val = (0.d0,0.d0) 
!             call destroy(work1)
!
!          endif   
!
!      enddo
!
!  enddo
!   
!  ! throw away all non diagonal parts
!  !do i = 1, nbl
!  !    i1 = Sigma_r(i,i)%nrow
!  !    do m = 1, i1 
!  !       Sigma_r(i,i)%val(1:m-1,m) = (0.0_dp,0.0_dp) 
!  !       Sigma_r(i,i)%val(m+1:i1,m) = (0.0_dp,0.0_dp) 
!  !    enddo
!  ! enddo
!  !
!
!  do i = 1, nbl
!   !print*          
!   !print*,'(sigma_r) Sigma_ph_r',maxval(abs(Sigma_r(i,i)%val)), iE
!    call write_blkmat(Sigma_r(i,i),pnegf%scratch_path,'Sigma_ph_r_',i,i,iE)
!    call destroy(Sigma_r(i,i))
!    call destroy(G_r_interP(i,i))
!    call destroy(G_r_interN(i,i))
!    call destroy(G_n_interP(i,i))
!    call destroy(G_n_interN(i,i))
!    !call destroy(G_r(i,i))
!  enddo
!
!  deallocate(Sigma_r,G_r_interP,G_r_interN,G_n_interP,G_n_interN)
!  !deallocate(G_r)
!
! END SUBROUTINE Sigma_ph_r
!
! SUBROUTINE Sigma_ph_r_z(pnegf,z)
!
!  TYPE(Tnegf) :: pnegf
!  COMPLEX(dp), intent(in) :: z
!  
!  !Local variables
!  TYPE(z_DNS) :: work1,work2
!  TYPE(z_DNS) :: Mq_mat
!  TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_r, G_r
!  REAL(dp), DIMENSION(:), POINTER :: Wq
!  REAL(dp), DIMENSION(:), POINTER :: Nq
!  REAL(dp), DIMENSION(:), POINTER :: Mq
!  LOGICAL, DIMENSION(:), POINTER :: selmodes
!  INTEGER :: i, m, iE,i1
!  INTEGER :: nummodes, numselmodes, nbl
!  INTEGER, DIMENSION(:), pointer :: indblk
!
!  nbl = pnegf%str%num_PLs
!  indblk => pnegf%str%mat_PL_start
!
!
!  selmodes => pnegf%elph%selmodes
!  Mq => pnegf%elph%Mq
!  Wq => pnegf%elph%Wq
!  Nq => pnegf%elph%Nq
!  nummodes = pnegf%elph%nummodes
!  numselmodes = pnegf%elph%numselmodes
!
!  iE=pnegf%iE
!
!  allocate(Sigma_r(nbl,nbl))
!  allocate(G_r(nbl,nbl))
!
!  do i = 1, nbl
!      m = indblk(i+1)-indblk(i)
!      call create(Sigma_r(i,i), m, m)
!      Sigma_r(i,i)%val = (0.d0, 0.d0)
!      call create(G_r(i,i), m, m)
!  enddo 
!
!  do m = 1 , nummodes
!
!    if (.not.selmodes(m)) cycle
!    
!    do i = 1, nbl
!    
!      call read_blkmat(G_r(i,i),pnegf%scratch_path,'G_r_',i,i,iE)
!
!      i1 = Sigma_r(i,i)%nrow
!
!      if (pnegf%elph%selfene_gr) then
!         
!          call create(work1,i1,i1)
! 
!          work1%val = (0.d0,0.d0)
!          !print*,'SelfEneR_Gr'    
!          work1%val = (2.0_dp*Nq(m)+1.0_dp)*G_r(i,i)%val
!
!          call create_id(Mq_mat,i1,Mq(m))
!                   
!          call prealloc_mult(Mq_mat, work1, work2)
!
!          call destroy(work1)
!
!          call prealloc_mult(work2, Mq_mat, work1)
!
!          call destroy(work2)
!
!          Sigma_r(i,i)%val = Sigma_r(i,i)%val +  work1%val 
!
!          call destroy(work1, Mq_mat)
!       
!       else
!          Sigma_r(i,i)%val = (0.d0,0.d0) 
!       endif
!    
!     enddo  
!     
!   enddo  
!
!   do i = 1, nbl
!    !print*,'(sigma_r) Sigma_ph_r', pnegf%iE, maxval(abs(Sigma_r(i,i)%val))
!    call write_blkmat(Sigma_r(i,i),pnegf%scratch_path,'Sigma_ph_r_',i,i,iE)
!    call destroy(Sigma_r(i,i))
!    call destroy(G_r(i,i))
!  enddo
!
! END SUBROUTINE Sigma_ph_r_z
!
! ! ----------------------------------------------------------
! SUBROUTINE check_Gl_Gr(pnegf)
!  TYPE(Tnegf) :: pnegf
!  REAL(dp), DIMENSION(:),allocatable :: Epnt
!  integer :: ioffset
!  
!  TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: G_r, G_p, G_n
!  TYPE(z_DNS) :: A, T 
!  integer :: nbl, n, sizebl, i_start, i_stop, psize, maxpos(2)
!  real(dp) :: Wmax, maxdev, tmp, maxG
!  INTEGER, DIMENSION(:), pointer :: indblk
!
!  nbl = pnegf%str%num_PLs
!  indblk => pnegf%str%mat_PL_start
!
!  allocate(G_n(nbl,nbl))
!  allocate(G_p(nbl,nbl))
!  allocate(G_r(nbl,nbl))
!
!  maxdev = 0.0_dp
!  maxG = 0.0_dp
!  psize = 0
!
!  do n = 1, nbl
!    sizebl = indblk(n+1)-indblk(n)
!    call create(G_r(n,n),sizebl,sizebl)
!    call create(G_n(n,n),sizebl,sizebl)
!    call create(G_p(n,n),sizebl,sizebl)
!    call read_blkmat(G_r(n,n), pnegf%scratch_path, 'G_r_',n,n, pnegf%iE)
!    call read_blkmat(G_n(n,n),pnegf%scratch_path,'G_n_',n,n,pnegf%iE)
!    call read_blkmat(G_p(n,n),pnegf%scratch_path,'G_p_',n,n,pnegf%iE)
!
!    call zspectral(G_r(n,n),G_r(n,n),0,A)
!
!    call create(T,sizebl,sizebl)
!
!    !print*,'CHECK Sigma_ph_n+Sigma_ph_p',maxval(abs(Sigma_n(n,n)%val+Sigma_p(n,n)%val))
!  
!    T%val = G_n(n,n)%val + G_p(n,n)%val - A%val 
!    
!    tmp = maxval(abs(A%val))
!    if (tmp .gt. maxG) maxG=tmp
!
!    tmp = maxval(abs(T%val))/maxval(abs(A%val)) 
!    if (tmp .gt. maxdev) then
!        maxdev = tmp
!        maxpos = maxloc(abs(T%val)) + psize
!    endif     
!
!    !print*
!  
!    psize = psize + sizebl
!
!    call destroy(G_r(n,n))
!    call destroy(G_n(n,n))
!    call destroy(G_p(n,n))
!    call destroy(A)
!    call destroy(T)
!
!  enddo
!
!  print*,'CHECK Gn+Gp=Gr-Ga',pnegf%iE, maxG, maxdev
!
!  deallocate(G_n, G_p, G_r)
!
! END SUBROUTINE check_Gl_Gr
!
! ! ----------------------------------------------------------
! SUBROUTINE check_sigma_ph_r(pnegf)
!  TYPE(Tnegf) :: pnegf
!  
!  TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: Sigma_r, Sigma_p, Sigma_n
!  TYPE(z_DNS) :: Gam, T 
!  integer :: nbl, n, sizebl, psize, maxpos(2)
!  real(dp) :: maxdev, tmp, maxG
!  INTEGER, DIMENSION(:), pointer :: indblk
!
!  nbl = pnegf%str%num_PLs
!  indblk => pnegf%str%mat_PL_start
!
!  allocate(Sigma_n(nbl,nbl))
!  allocate(Sigma_p(nbl,nbl))
!  allocate(Sigma_r(nbl,nbl))
!
!  maxdev = 0.0_dp
!  maxG = 0.0_dp
!  psize = 0
!
!  do n = 1, nbl
!    sizebl = indblk(n+1)-indblk(n)
!    call create(Sigma_r(n,n),sizebl,sizebl)
!    call create(Sigma_n(n,n),sizebl,sizebl)
!    call create(Sigma_p(n,n),sizebl,sizebl)
!    call read_blkmat(Sigma_r(n,n), pnegf%scratch_path, 'Sigma_ph_r_',n,n, pnegf%iE)
!    call read_blkmat(Sigma_n(n,n),pnegf%scratch_path,'Sigma_ph_n_',n,n,pnegf%iE)
!    call read_blkmat(Sigma_p(n,n),pnegf%scratch_path,'Sigma_ph_p_',n,n,pnegf%iE)
!
!    call zspectral(Sigma_r(n,n),Sigma_r(n,n),0,Gam)
!
!    call create(T,sizebl,sizebl)
!
!    !print*,'CHECK Sigma_ph_n+Sigma_ph_p',maxval(abs(Sigma_n(n,n)%val+Sigma_p(n,n)%val))
!  
!    T%val = Sigma_n(n,n)%val + Sigma_p(n,n)%val - Gam%val 
!    
!    tmp = maxval(abs(Gam%val))
!    if (tmp .gt. maxG) maxG=tmp
!
!    tmp = maxval(abs(T%val))/maxval(abs(Gam%val)) 
!    if (tmp .gt. maxdev) then
!        maxdev = tmp
!        maxpos = maxloc(abs(T%val)) + psize
!    endif     
!
!    !print*
!  
!    psize = psize + sizebl
!
!    call destroy(Sigma_r(n,n))
!    call destroy(Sigma_n(n,n))
!    call destroy(Sigma_p(n,n))
!    call destroy(Gam)
!    call destroy(T)
!
!  enddo
!
!  print*,'CHECK Sigma_ph_r',pnegf%iE, maxG, maxdev
!
!  deallocate(Sigma_n, Sigma_p, Sigma_r)
!
!
! END SUBROUTINE check_sigma_ph_r
!
!
! ! ----------------------------------------------------------
! ! Search points for interpolations. 
! ! Wq has a sign (+/- Wq)
! !
! SUBROUTINE search_points(pnegf, Wq, Epnt, i1, i2, E1, E2)
!   use negf_energy_mesh, only : elem 
!   type(TNegf) :: pnegf
!   real(dp) :: Wq                
!   real(dp), dimension(:), allocatable :: Epnt
!   integer, intent(out) :: i1,i2 
!   real(dp), intent(out) :: E1, E2
!  
!   integer :: iE, iel, istart, iend, ip
!   Type(elem), pointer :: pel
!   real(dp) :: En
!  
!   if (Wq.eq.0) then
!      i1 = pnegf%iE
!      i2 = pnegf%iE
!      E1 = pnegf%Epnt
!      E2 = pnegf%Epnt
!      return
!   endif
!  
!   if (allocated(Epnt)) then
!     ! Remove offset such that search can work on Epnt(1..N)      
!     iE = pnegf%iE - pnegf%Np_n(1) - pnegf%Np_n(2) - pnegf%n_poles
!     En = real(pnegf%Epnt) + Wq !Wq carry the right sign      
!     !print*
!     !print*,'iE', iE, real(pnegf%Epnt) + Wq 
!  
!     if (sign(1.0_dp,Wq) .gt. 0) then
!       i2 = iE + 1
!       if (i2.gt.size(Epnt)) then
!          i2 = size(Epnt)
!       else    
!          do while (Epnt(i2) .lt. En)
!             i2 = i2 + 1
!          end do
!       endif
!       i1 = i2 - 1
!     else
!       i1 = iE - 1
!       if (i1.lt.1) then
!          i1 = 1 
!       else
!          do while (Epnt(i1) .gt. En)
!             i1 = i1 - 1
!          end do
!       endif   
!       i2 = i1 + 1
!     endif
!     E1 = Epnt(i1)
!     E2 = Epnt(i2)
!     ! add back the offset to the point
!     i1 = i1 + pnegf%Np_n(1) + pnegf%Np_n(2) + pnegf%n_poles
!     i2 = i2 + pnegf%Np_n(1) + pnegf%Np_n(2) + pnegf%n_poles
!   
!   else  
!   
!     !if (.not.allocated(pnegf%emesh%pactive)) STOP 'emesh not initialized'
!  
!     En = real(pnegf%Epnt) + Wq       
!           
!     if (sign(1.0_dp,Wq) .gt. 0) then
!       istart = pnegf%emesh%iactive
!       iend = pnegf%emesh%maxind
!  
!       elloop1: do iel = istart, iend
!          pel => pnegf%emesh%pactive(iel)%pelem
!          do ip = 1, 3
!             if (pel%pnt(ip) .gt. En) then
!                exit elloop1 
!             endif     
!          end do
!       end do elloop1   
!       i1 = pel%map(ip-1)
!       i2 = pel%map(ip)
!       E1 = pel%pnt(ip-1)
!       E2 = pel%pnt(ip)
!     else
!       istart = pnegf%emesh%iactive
!       iend = 1
!  
!       elloop2: do iel = istart, iend, -1
!          pel => pnegf%emesh%pactive(iel)%pelem
!          do ip = 3, 1, -1
!             if (pel%pnt(ip) .lt. En) then
!                exit elloop2 
!             endif     
!          end do
!       end do elloop2  
!       i1 = pel%map(ip)
!       i2 = pel%map(ip+1)
!       E1 = pel%pnt(ip)
!       E2 = pel%pnt(ip+1)
!     end if
!   end if
!   !print*
!   !print*,E1,En,E2
!   !print*,'interpolate between:',i1,i2
!
! END SUBROUTINE search_points 
!
!
! SUBROUTINE interpolation(i1,i2, E1, E2, E, path, name, G_interp)
!     INTEGER, intent(in) :: i1, i2
!     REAL(dp) :: E1, E2, E
!     TYPE(z_DNS), DIMENSION(:,:) :: G_interp
!     CHARACTER(*) :: path
!     CHARACTER(*) :: name
!
!     !local variables
!     TYPE(z_DNS) :: work1,work2
!     INTEGER :: i
!
!     do i = 1, size(G_interp,1)
!        
!        call create(work1,G_interp(i,i)%nrow,G_interp(i,i)%ncol)
!        work1%val = (0.d0,0.d0)
!        call read_blkmat(work1, path, name, i, i, i1)
!
!        if (E1.ne.E2 .and. i1.ne.i2) then
!         
!          call create(work2,G_interp(i,i)%nrow,G_interp(i,i)%ncol)
!          work2%val = (0.d0,0.d0)
!          call read_blkmat(work2, path, name, i, i, i2)
!
!          G_interp(i,i)%val = ((E-E1)*work2%val + (E2-E)*work1%val)/(E2-E1)
!
!          call destroy(work2)
!      
!        else
!
!          G_interp(i,i)%val = work1%val
!        
!        endif
!        
!        call destroy(work1)
!
!     end do
!
!
! END SUBROUTINE interpolation
!
!
! !*********************************************************************
! !ADD Hilbert-transform part to Sigma_ph_r  
! !Need to set back FFT transforms
! !********************************************************************
! SUBROUTINE complete_sigma_ph_r(pnegf, Epnt, ioffset)
!
!  TYPE(Tnegf) :: pnegf
!  REAL(dp), DIMENSION(:) :: Epnt
!  INTEGER :: ioffset
!
!  ! Locals
!  INTEGER :: i,j,n,m, iE, bl, sizebl,i_start,i_stop, ierr, nummodes
!  REAL(dp), DIMENSION(:), POINTER :: Wq
!  REAL(dp), DIMENSION(:), POINTER :: Mq
!  REAL(dp) :: Wmax, dE, tmp
!  
!  character(4) :: ofblki
!  complex(dp), DIMENSION(:), allocatable :: temp1 
!  type(z_DNS), DIMENSION(:), allocatable :: Gn_E, Sigma_r_E
!  !type(z_DNS) :: work1, work2, work3, Mq_mat 
!  LOGICAL, DIMENSION(:), POINTER :: selmodes
!
!  Mq => pnegf%elph%Mq
!  Wq => pnegf%elph%Wq
!  selmodes => pnegf%elph%selmodes
!
!  n = size(Epnt)
!  nummodes = size(Wq)
!
!  Wmax = maxval(Wq)
!  dE = Epnt(2)-Epnt(1)
!  ! ENERGY-INTERVAL FOR SELF-ENERGIES:
!  i_start = aint(Wmax/dE) + 1
!  i_stop =  (n - i_start) 
!  i_start = i_start + 1
!  
!  ! PROBLEM: 
!  ! The Hilb-Transf should be resticted on the sub interval
!  ! Emin+m*Wmax ... Emax-m*Wmax
!  ! However the number of points should be always 2^p
!  ! This condition is difficult to fulfill.
!  ! Possible solution: interpolation between grids.
!
!  ! CREATE THE ARRAY OF MATRICES (one for each E-point)
!  allocate(Gn_E(n),stat=ierr)
!  allocate(Sigma_r_E(n),stat=ierr)
!  if(ierr.ne.0) STOP 'ERROR in allocation of Gn_E'
!  call log_allocate(temp1,n)
!
!  print*
!  print*,'HILBERT TRANSFORM memory:',sizebl*sizebl*n*16
!  print*,'HILBERT TRANSFORM interval:',i_start,i_stop
! 
!  ! LOOP ON blocks
!  do bl = 1, pnegf%str%num_PLs
!
!    sizebl =  pnegf%str%mat_PL_start(bl+1) - pnegf%str%mat_PL_start(bl)
!    if (bl.le.9999) write(ofblki,'(i4.4)') bl 
!    if (bl.gt.9999) stop 'ERROR: too many blks (> 9999)'
! 
!    ! LOAD ALL G< FROM FILES and store in Gn_E(iE)
!    do iE = 1, n
!      call create(Gn_E(iE),sizebl,sizebl)
!      call read_blkmat(Gn_E(iE), pnegf%scratch_path, 'G_n_', bl, bl,iE+ioffset)
!    enddo
!    
!    ! LOAD ALL Sigma_r in the right interval
!    do iE = i_start, i_stop
!      call create(Sigma_r_E(iE),sizebl,sizebl)
!      call read_blkmat(Sigma_r_E(iE), pnegf%scratch_path, 'Sigma_ph_r_', bl, bl, iE+ioffset)
!    enddo
!
!    do m = 1, nummodes 
!      if (.not.selmodes(m)) cycle
!
!      ! COMPUTE   Mq G< Mq   (assume diagonal now) 
!      do iE = 1, n
!        !call create_id(Mq_mat,sizebl,Mq(m))
!        !call prealloc_mult(Mq_mat, Gn_E(iE), work1)
!        !call prealloc_mult(work1, Mq_mat, work2)
!        !call destroy(work1, Mq_mat)
!        !Gn_E(iE)%val = work2%val
!        !call destroy(work2)
!        tmp = Mq(m)*Mq(m)
!        Gn_E(iE)%val = tmp*Gn_E(iE)%val
!      enddo
!
!      ! PERFORM Hilbert in Energy for each (i,j)     
!      do j = 1, sizebl  
!        do i = 1, sizebl
!
!!          k = (j-1)*sizebl+i
!!          do iE = 1, n
!!             if (iE.le.99999) write(ofpnt,'(i5.5)') iE+ioffset
!!             filename = 'G_n_'//trim(ofblki)//'_'//trim(ofblki)//'_'//trim(ofpnt)//'.dat'
!!             open(100,file=trim(pnegf%scratch_path)//trim(filename), access='DIRECT', recl=4)
!!             READ(100, rec = k)  temp1(iE)
!!             close(100)
!!          enddo
!
!           ! SETUP a vector out of all G<_ij(E)
!           ! Here we could perform an efficient interpolation on a regular grid 2^p
!           do iE = 1, n
!             temp1(iE) = Gn_E(iE)%val(i,j)
!           enddo
! 
!           Wq(m) = Wq(m)*2.0_dp*pi/((n-1)*dE)
! 
!           !call Hilbert_shift(temp1, Wq(m))
!         
!           Wq(m) = Wq(m)*((n-1)*dE)/(2.0_dp*pi)
! 
!           ! UPDATE the self-energies with the Hilbert part. 
!           do iE = i_start, i_stop
!             ! Should be  -i/2 H[ G<+ - G<-  ] = 1/2 H[ Gn+ - Gn- ]
!             Sigma_r_E(iE)%val(i,j) =  Sigma_r_E(iE)%val(i,j) + (0.5_dp, 0.0)* temp1(iE) 
!           enddo
!
!        enddo ! Loop on block size   
!      enddo ! Loop on block size 
!
!    enddo !Loop on modes 
!
!    do iE = i_start, i_stop
!      call write_blkmat(Sigma_r_E(iE), pnegf%scratch_path, 'Sigma_ph_r_', bl, bl, iE+ioffset)
!    enddo
!
!    do iE = 1, n
!      call destroy(Gn_E(iE))
!    enddo
!    do iE = i_start, i_stop
!      call destroy(Sigma_r_E(iE))
!    enddo
!
!   enddo !Loop on blocks 
!  
!   call log_deallocate(temp1)
!   deallocate(Gn_E)
!   deallocate(Sigma_r_E)
!
! END SUBROUTINE complete_sigma_ph_r 

 !---------------------------------------------------------
 !---------------------------------------------------------
 SUBROUTINE create_scratch(nbl, npoints)
   integer :: nbl, npoints

   integer :: i,j,k,err
   
   call destroy_scratch(nbl, npoints)
   
   ALLOCATE(DDn(nbl,nbl),stat=err)
   ALLOCATE(DDr(nbl,nbl),stat=err)
   ALLOCATE(DDp(nbl,nbl),stat=err)
   ALLOCATE(Pin(nbl,nbl),stat=err)
   ALLOCATE(Pir(nbl,nbl),stat=err)
   ALLOCATE(Pip(nbl,nbl),stat=err)
  
   ! Initialize everything to 0 
   do j=1,nbl
     do i=1,nbl  
         DDn(i,j)%nrow=0
         DDn(i,j)%ncol=0
         DDn(i,j)%npoints=0
         DDr(i,j)%nrow=0
         DDr(i,j)%ncol=0
         DDr(i,j)%npoints=0
         DDp(i,j)%nrow=0
         DDp(i,j)%ncol=0
         DDp(i,j)%npoints=0
         Pin(i,j)%nrow=0
         Pin(i,j)%ncol=0
         Pin(i,j)%npoints=0
         Pir(i,j)%nrow=0
         Pir(i,j)%ncol=0
         Pir(i,j)%npoints=0
         Pip(i,j)%nrow=0
         Pip(i,j)%ncol=0
         Pip(i,j)%npoints=0
     enddo
   enddo
   
   if(err.ne.0) then
      STOP 'ERROR: Cannot allocate GG'
   endif     
   print*,'Created memory scratch',nbl,'x',nbl,'x',npoints

 END SUBROUTINE create_scratch
 !---------------------------------------------------------
 !---------------------------------------------------------

 SUBROUTINE destroy_scratch(nbl, npoints)
   integer :: nbl, npoints
   
   integer :: i, j, iE,  err
   
   err = 0
   
   
   do i = 1, nbl
     do j = 1, nbl
         
        if (allocated(DDn)) then
            if (allocated(DDn(i,j)%val)) call destroy(DDn(i,j))
        endif       
        if (allocated(DDr)) then
            if (allocated(DDr(i,j)%val)) call destroy(DDr(i,j))
        endif       
        if (allocated(DDp)) then
            if (allocated(DDp(i,j)%val)) call destroy(DDp(i,j))
        endif       
        if (allocated(Pin)) then
            if (allocated(Pin(i,j)%val)) call destroy(Pin(i,j))     
        endif       
        if (allocated(Pir)) then
            if (allocated(Pir(i,j)%val)) call destroy(Pir(i,j))
        endif       
        if (allocated(Pip)) then
            if (allocated(Pip(i,j)%val)) call destroy(Pip(i,j))     
        endif       
       
     end do
   end do
   
   if (allocated(DDn)) DEALLOCATE(DDn,stat=err)
   if (allocated(DDr)) DEALLOCATE(DDr,stat=err)
   if (allocated(DDp)) DEALLOCATE(DDp,stat=err)
   if (allocated(Pin)) DEALLOCATE(Pin,stat=err)
   if (allocated(Pir)) DEALLOCATE(Pir,stat=err)
   if (allocated(Pip)) DEALLOCATE(Pip,stat=err)
   
   if(err.ne.0) then
      STOP 'ERROR: Cannot deallocate GG'
   endif

 END SUBROUTINE destroy_scratch
 !---------------------------------------------------------
 !---------------------------------------------------------

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
          Matrix%val = DDn(i,j)%val(:,:,iE) 
       case('G_r_')                   
          Matrix%val = DDr(i,j)%val(:,:,iE)
       case('G_p_')                   
          Matrix%val = DDp(i,j)%val(:,:,iE)
       case('Sigma_ph_r_')          
          Matrix%val = Pir(i,j)%val(:,:,iE) 
       case('Sigma_ph_n_')          
          Matrix%val = Pin(i,j)%val(:,:,iE)
       case('Sigma_ph_p_')          
          Matrix%val = Pip(i,j)%val(:,:,iE)
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
  SUBROUTINE write_blkmat(Matrix, path, name, i, j, iE, npoints)
    TYPE(z_DNS), intent(in) :: Matrix
    CHARACTER(*), intent(in) :: path
    CHARACTER(*), intent(in) :: name
    INTEGER, intent(in) :: i, j, iE, npoints

    CHARACTER(64) :: filename
    character(4) :: ofblki, ofblkj
    character(10) :: ofpnt
    integer :: m,n,nrow,ncol

    if (memory) then
              
       nrow=Matrix%nrow
       ncol=Matrix%ncol

       select case(trim(name))
       case('G_n_')     
          if (.not.allocated(DDn(i,j)%val)) then
              call create(DDn(i,j),nrow,ncol,npoints)
          endif    
          DDn(i,j)%val(:,:,iE) = Matrix%val
       case('G_r_')     
          if (.not.allocated(DDr(i,j)%val)) then
              call create(DDr(i,j),nrow,ncol,npoints)
          endif    
          DDr(i,j)%val(:,:,iE) = Matrix%val
       case('G_p_')     
          if (.not.allocated(DDp(i,j)%val)) then
              call create(DDp(i,j),nrow,ncol,npoints)
          endif    
          DDp(i,j)%val(:,:,iE) = Matrix%val 
       case('Sigma_ph_r_')     
          if (.not.allocated(Pir(i,j)%val)) then
              call create(Pir(i,j),nrow,ncol,npoints)
          endif    
          Pir(i,j)%val(:,:,iE) = Matrix%val
       case('Sigma_ph_n_')     
          if (.not.allocated(Pin(i,j)%val)) then
              call create(Pin(i,j),nrow,ncol,npoints)
          endif    
          Pin(i,j)%val(:,:,iE) = Matrix%val
       case('Sigma_ph_p_')     
          if (.not.allocated(Pip(i,j)%val)) then
              call create(Pip(i,j),nrow,ncol,npoints)
          endif    
          Pip(i,j)%val(:,:,iE) = Matrix%val
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
  !  Calculate tunneling  
  !
  !****************************************************************************

  SUBROUTINE calculate_transmissions(H,S,Ec,SelfEneR,ni,nf,size_ni,str,TUN_MAT)

    implicit none

    Type(z_CSR) :: H
    Type(z_CSR) :: S           
    Complex(dp) :: Ec
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR
    !Type(z_DNS) :: SelfEner_d 
    Real(dp), Dimension(:) :: TUN_MAT
    Type(z_CSR) :: ESH_csr
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
    call prealloc_sum(H,S,(-1.d0, 0.d0),Ec,ESH_csr)    

    call allocate_blk_dns(ESH,nbl)
    
    call csr2blkdns(ESH_csr,ESH,str%mat_PL_start)
    call destroy(ESH_csr)

    !Inclusion of the contact Self-Energies to the relevant blocks
    do i=1,ncont
       ESH(str%cblk(i),str%cblk(i))%val = ESH(str%cblk(i),str%cblk(i))%val-SelfEneR(i)%val
    enddo
    call allocate_gsm(gsmr,nbl)

    !Iterative calculation up with gsmr 
    call calculate_gsmr_blocks(ESH,nbl,2)


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
    call calculate_Gr_tridiag_blocks(ESH,1)
    if (nt.gt.1) call calculate_Gr_tridiag_blocks(ESH,2,nt)
    
    if (ncont.eq.2) then
       call calculate_single_transmission_2_contacts(nit,nft,ESH,SelfEneR,str%cblk,tun) 
    else
       call calculate_single_transmission_N_contacts(nit,nft,ESH,SelfEneR,str%cblk,tun) 
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
          call calculate_Gr_tridiag_blocks(ESH,nt+1,nt1)
          nt = nt1
       endif

       call calculate_single_transmission_N_contacts(nit,nft,ESH,SelfEneR,str%cblk,tun) 
  
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
    call deallocate_gsm(gsmr)

    call destroy_ESH(ESH)
    call deallocate_blk_dns(ESH)

  end SUBROUTINE calculate_transmissions

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

  subroutine calculate_single_transmission_2_contacts(ni,nf,ESH,SelfEneR,cblk,TUN)

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
       write(*,*) "ERROR: calculate_single_transmission_2_contacts is valid only for 2 contacts"
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
   
  end subroutine calculate_single_transmission_2_contacts

  !************************************************************************
  !
  ! Subroutine for transmission calculation (GENERIC FOR N CONTACTS)
  !
  !************************************************************************ 
  subroutine calculate_single_transmission_N_contacts(ni,nf,ESH,SelfEneR,cblk,TUN)
    
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
    
  end subroutine calculate_single_transmission_N_contacts
 

  !---------------------------------------------------!
  !Subroutine for transmission and dos calculation    !
  !---------------------------------------------------!  

  subroutine calculate_transmissions_and_dos(H,S,Ec,SelfEneR,Gs,ni,nf,nLDOS,LDOS,size_ni,str,TUN_MAT,LEDOS)

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
    Type(z_CSR) :: ESH_csr, GrCSR
    Type(z_DNS), Dimension(:,:), allocatable :: ESH
    Type(r_CSR) :: Grm                          ! Green Retarded nella molecola           
    real(dp), dimension(:), allocatable :: diag
    Real(dp) :: tun
    Complex(dp) :: zc
    Integer :: nbl,ncont, ierr
    Integer :: nit, nft, icpl
    Integer :: iLDOS, i2, i, cb
    Character(1) :: Im

    nbl = str%num_PLs
    ncont = str%num_conts
    Im = 'I'

    !Calculation of ES-H and brak into blocks
    call prealloc_sum(H,S,(-1.d0, 0.d0),Ec,ESH_csr)    

    call allocate_blk_dns(ESH,nbl)
    call csr2blkdns(ESH_csr,ESH,str%mat_PL_start)
    call destroy(ESH_csr)

    !Inclusion of the contact Self-Energies to the relevant blocks
    do i=1,ncont
       cb = str%cblk(i)
       ESH(cb,cb)%val = ESH(cb,cb)%val-SelfEneR(i)%val
    enddo

    call allocate_gsm(gsmr,nbl)
    
    call calculate_gsmr_blocks(ESH,nbl,2)

    call allocate_blk_dns(Gr,nbl)

    call calculate_Gr_tridiag_blocks(ESH,1)
    call calculate_Gr_tridiag_blocks(ESH,2,nbl)

    !Computation of transmission(s) between contacts ni(:) -> nf(:)
    do icpl=1,size_ni

       nit=ni(icpl)
       nft=nf(icpl)
       
       if (ncont.eq.2) then
         call calculate_single_transmission_2_contacts(nit,nft,ESH,SelfEneR,str%cblk,tun) 
       else
         call calculate_single_transmission_N_contacts(nit,nft,ESH,SelfEneR,str%cblk,tun)
       endif

       TUN_MAT(icpl) = tun 
    
    enddo
   
    !Deallocate energy-dependent matrices
    call destroy_gsm(gsmr)
    call deallocate_gsm(gsmr)
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

  end subroutine calculate_transmissions_and_dos


  !---------------------------------------------------


  subroutine allocate_gsm(gsm,nbl)
    type(z_DNS), dimension(:), allocatable :: gsm 
    integer :: nbl, ierr

    allocate(gsm(nbl),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate gsm'
    
  end subroutine allocate_gsm
  
  !---------------------------------------------------

  subroutine allocate_blk_dns(blkM,nbl)
    type(z_DNS), dimension(:,:), allocatable :: blkM 
    integer :: nbl, ierr

    allocate(blkM(nbl,nbl),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate block-Matrix'

  end subroutine allocate_blk_dns

  !---------------------------------------------------

  subroutine deallocate_gsm(gsm)
    type(z_DNS), dimension(:), allocatable :: gsm 
    integer :: ierr

    deallocate(gsm,stat=ierr)
    if (ierr.ne.0) stop 'DEALLOCATION ERROR: could not deallocate gsmr'

  end subroutine deallocate_gsm

  !---------------------------------------------------

  subroutine deallocate_blk_dns(blkM)
    type(z_DNS), dimension(:,:), allocatable :: blkM 
    integer :: ierr

    deallocate(blkM,stat=ierr)
    if (ierr.ne.0) stop 'DEALLOCATION ERROR: could not deallocate block-Matrix'

  end subroutine deallocate_blk_dns

  !---------------------------------------------------


END  module negf_iterative_ph


! SUBROUTINE calls_neq_ph(pnegf,E,SelfEneR,Tlc,Tcl,gsurfR,frm,Glout,out)
!
! !****************************************************************************
! !
! !Input
! !H: sparse matrix contaning Device Hamiltonian
! !S: sparse matrix containing Device Overlap
! !SelfEneR: sparse matrices array containing contacts Self Energy
! !Tlc: sparse matrices array containing contacts-device interaction blocks (ES-H)
! !Tcl: sparse matrices array containing device-contacts interaction blocks (ES-H)
! !gsurfR: sparse matrices array containing contacts surface green
! !frm: array containing Fermi distribution values for all contacts
! !ref: reference contact excluded from summation
! !
! !Output:
! !Aout: G_n contributions (Device + Contacts overlap regions -> effective conductor)
! !
! !*****************************************************************************
! IMPLICIT NONE
!
! !In/Out
! TYPE(Tnegf) :: pnegf
! TYPE(z_CSR) :: Glout
! TYPE(z_DNS), DIMENSION(:) :: SelfEneR, gsurfR, Tlc, Tcl
! REAL(dp) :: E
! REAL(dp), DIMENSION(:) :: frm
! INTEGER :: out
! 
! !Work
! COMPLEX(dp) :: Ec
! INTEGER :: i,ierr,ncont,nbl, lbl
! INTEGER :: ref, iter
! INTEGER, DIMENSION(:), POINTER :: cblk, indblk
! TYPE(z_DNS), DIMENSION(:,:), ALLOCATABLE :: ESH, Gn, Gp
! TYPE(z_CSR) :: ESH_csr, Gl
! LOGICAL :: mask(MAXNCONT)
! TYPE(Tstruct_info) :: struct
! REAL(dp), DIMENSION(:), allocatable :: cfrm
! 
! struct = pnegf%str
! nbl = struct%num_PLs
! ncont = struct%num_conts
! indblk => struct%mat_PL_start
! cblk => struct%cblk
! ref = pnegf%refcont
! iter = pnegf%elph%scba_iter
! 
! Ec=cmplx(E,0.d0,dp)
! 
! !if (debug) write(*,*) '----------------------------------------------------'
! !if (debug) call writeMemInfo(6)
! !if (debug) call writePeakInfo(6)
! !if (debug) write(*,*) '----------------------------------------------------'
! 
!  !Costruiamo la matrice sparsa ESH
!  CALL prealloc_sum(pnegf%H,pnegf%S,(-1.d0, 0.d0),Ec,ESH_csr)
! 
!  call allocate_blk_dns(ESH,nbl)
! 
!  CALL csr2blkdns(ESH_csr,ESH,indblk)
! 
!  call destroy(ESH_csr)
! 
!  DO i=1,ncont
!        ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
!  ENDDO
! 
!  ! Reload and add Sigma_ph_r to ESH
!  call add_sigma_ph_r(pnegf, ESH)
! 
!  !Allocazione delle gsmr
!  call allocate_gsm(gsmr,nbl)
!  call allocate_gsm(gsml,nbl)
! 
!  !Chiamata di Make_gsmr_mem
!  CALL calculate_gsmr_blocks(ESH,nbl,2)
! 
!  !Chiamata di Make_gsml_mem solo per i blocchi 1..lbl dove  
!  ! lbl = maxval(cblk,mask) - 2 
!  lbl = nbl - 2    ! ALL BLOCKS are needed!!
! 
!  if( ncont.gt.1 ) then
!   CALL calculate_gsml_blocks(ESH,1,lbl)    
!  endif
! 
!  call allocate_blk_dns(Gr,nbl)
! 
!  CALL calculate_Gr_tridiag_blocks(ESH,1)
!  CALL calculate_Gr_tridiag_blocks(ESH,2,nbl)
! 
!
!  !Chiamata di Make_Grcol_mem per i contatti necessari 
!  !DO i=1,ncont 
!  !  IF (i.NE.ref) THEN
!  !    CALL calculate_Gr_column_blocks(ESH,cblk(i),indblk)
!  !  ENDIF
!  !ENDDO
! 
!  !With el-ph we need all columns
!  DO i=1,nbl
!    CALL calculate_Gr_column_blocks(ESH,i,indblk)
!  ENDDO  
! 
!  !Distruzione delle gsmall
!  call destroy_gsm(gsmr)
!  call deallocate_gsm(gsmr)
!  call destroy_gsm(gsml)
!  call deallocate_gsm(gsml)
! 
!  !Calcolo degli outer blocks 
!  !if (debug) write(*,*) 'Compute Outer_G_n' 
!  SELECT CASE (out)
!  CASE(0)
!  CASE(1)
!    CALL calculate_Gn_outer(Tlc,gsurfR,SelfEneR,pnegf%str,frm,ref,.false.,Glout)
!  CASE(2)
!    CALL calculate_Gn_outer(Tlc,gsurfR,SelfEneR,pnegf%str,frm,ref,.true.,Glout)
!  END SELECT
! 
!  !Calcolo della G_n nel device
!  !if (debug) write(*,*) 'Compute G_n' 
!  !Allocazione degli array di sparse
!  ALLOCATE(Gn(nbl,nbl),stat=ierr)
!  IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Gn(nbl,nbl)'
! 
!  call init_blkmat(Gn,ESH)
! 
!  CALL calculate_Gn_tridiag_blocks(ESH,SelfEneR,frm,ref,struct,Gn)
! 
!  call calculate_Gn_tridiag_elph_contributions(pnegf,ESH,iter,Gn)
! 
!  !print*
!  ! save diagonal blocks of G_n
!  DO i = 1, nbl
!    !print*,'(G_n) G_n',minval(abs(Gn(i,i)%val)), maxval(abs(Gn(i,i)%val))
!    call write_blkmat(Gn(i,i),pnegf%scratch_path,'G_n_',i,i,pnegf%iE)
!    !print*,'(G_r) G_r',minval(abs(Gn(i,i)%val)),maxval(abs(Gr(i,i)%val))
!    call write_blkmat(Gr(i,i),pnegf%scratch_path,'G_r_',i,i,pnegf%iE)
!  ENDDO
! 
! 
!  !open(5001,file=trim(pnegf%out_path)//'Lcurr.dat',position='APPEND')
!  !write(5001,'(ES15.5)',advance='NO') real(Ec)
!  !do i = 1, nbl
!  !  do i1 = 1, ESH(i,i)%nrow-1
!  !     Iloc = 2*aimag(Glsub(i,i)%val(i1,i1+1))*ESH(i,i)%val(i1+1,i1)
!  !     write(5001,'(ES15.5)',advance='NO') Iloc
!  !  enddo
!  !enddo    
!  !write(5001,*)
!  !close(5001)
! 
!  ! the following destroys Gn
!  call blk2csr(Gn,pnegf%str,pnegf%S,Gl)
! 
!  DEALLOCATE(Gn)
!  
!  !! Computation of G_p for Sigma_p
!  if (pnegf%elph%check) then
!     allocate(cfrm(10))
!     cfrm = 1.0_dp - frm
!     cfrm(ref) = 0.0_dp
!  
!     ALLOCATE(Gp(nbl,nbl),stat=ierr)
!  
!     IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Gp(nbl,nbl)'
!     call init_blkmat(Gp,ESH)
!     
!     CALL calculate_Gn_tridiag_blocks(ESH,SelfEneR,cfrm,ref,struct,Gp)
!     deallocate(cfrm)
!     
!     call Make_Gp_ph(pnegf,ESH,iter,Gp)
!     
!    !  save diagonal blocks of G_p
!    DO i = 1, nbl
!      !print*
!      !print*,'(G_p) G_p',minval(abs(Gn(i,i)%val)),maxval(abs(Gn(i,i)%val))
!      call write_blkmat(Gp(i,i),pnegf%scratch_path,'G_p_',i,i,pnegf%iE)
!    ENDDO
!     
!    call destroy_blk(Gp)
!    DEALLOCATE(Gp)
!  endif
!  !! ------------------------------------------------------------------------
! 
!  !Distruzione dell'array Gr
!  CALL destroy_blk(Gr)
!  DEALLOCATE(Gr)
! 
!  CALL destroy_ESH(ESH)
!  DEALLOCATE(ESH)
! 
!  !Concatenazione di Gl in Glout
!  SELECT CASE(out)
!  CASE(0)
!     call clone(Gl,Glout)
!  CASE(1:2)
!     call concat(Glout,Gl,1,1)
!  END SELECT
! 
!  call destroy(Gl)
!
!  END SUBROUTINE calls_neq_ph


