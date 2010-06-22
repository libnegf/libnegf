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




MODULE iterative

  USE allocation
  USE mat_def
  USE sparsekit_drv
  USE inversions
  USE structure, only : TStruct_Info

  !USE parameters, only : ncont, ncdim
  !use structure, only : nbl, indblk, cblk, cindblk
  !use clock

  private

  public :: calls_eq_mem
  public :: calls_neq_mem

  public :: sub_ESH
  public :: rebuild
  public :: Make_gsmr_mem
  public :: Make_gsml_mem
  public :: Make_Grdiag_mem
  public :: Make_Grcol_mem
  public :: Make_Spectral_mem
  public :: Make_GreenR_mem2  ! New. Best.
  public :: Make_Gl_mem
  public :: Outer_A_mem
  public :: Outer_GreenR_mem
  public :: Outer_Gl_mem

  LOGICAL, PARAMETER :: debug=.false. 
  !Dropout value
  REAL(8), PARAMETER :: drop=1e-20

  INTEGER :: t11,t21,cr1,cm1
  REAL(8) :: total_time,partial_time

  TYPE(z_CSR), DIMENSION(:), ALLOCATABLE :: gsmr
  INTEGER, DIMENSION(:), ALLOCATABLE :: gsmr_dim
  TYPE(z_CSR), DIMENSION(:), ALLOCATABLE :: gsml
  INTEGER, DIMENSION(:), ALLOCATABLE :: gsml_dim
  TYPE(z_CSR), DIMENSION(:,:), ALLOCATABLE :: Gr
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: Gr_dim

CONTAINS

  !****************************************************************************
  ! 
  ! Driver for computing Equilibrium Green Retarded (Outer blocks included)
  ! writing on memory
  !
  !****************************************************************************

  SUBROUTINE calls_eq_mem(H,S,E,SelfEneR,Tlc,Tcl,gsurfR,Aout,struct,lower_outer)

    !****************************************************************************
    !
    !Input
    !H: sparse matrix contaning Device Hamiltonian
    !S: sparse matrix containing Device Overlap
    !SelfEneR: sparse matrices array containing contacts Self Energy
    !Tlc: sparse matrices array containing contacts-device interaction blocks (ES-H)
    !Tcl: sparse matrices array containing device-contacts interaction blocks (ES-H)
    !gsurfR: sparse matrices array containing contacts surface green
    !min: collector index 
    !lower_outer: optional parameter. If defined (and true), Aout contains also 
    !the lower outer parts (needed for K-points calculations)
    !
    !Output:
    !Aout: Green Retarded (Device + Contacts overlap regions -> effective conductor)
    !
    !*****************************************************************************

    IMPLICIT NONE

    !In/Out
    TYPE(z_CSR) :: H, S, Aout
    TYPE(z_CSR), DIMENSION(:) :: SelfEneR, Tlc, Tcl, gsurfR 
    COMPLEX(8) :: E
    TYPE(Tstruct_info) :: struct
    LOGICAL, OPTIONAL :: lower_outer

    !Work
    TYPE(z_CSR), DIMENSION(:,:), ALLOCATABLE :: ESH
    TYPE(z_CSR) :: ESH_tot,A
    INTEGER :: i,ierr, nbl, ncont
    INTEGER, DIMENSION(:), POINTER :: cblk, indblk

    !if (debug) then 
    !   write(*,*) '----------------------------------------------------'
    !   call writeMemInfo(6)
    !   call writePeakInfo(6)
    !   write(*,*) '----------------------------------------------------'
    !end if
    nbl = struct%num_PLs
    ncont = struct%num_conts
    cblk => struct%cblk
    indblk => struct%mat_PL_start

    !Costruiamo la matrice sparsa ESH
    CALL prealloc_sum(H,S,(-1.d0, 0.d0),E,ESH_tot)

    !Costruiamo l'array di sparse ESH
    !Allocazione dell'array di sparse ESH
    ALLOCATE(ESH(nbl,nbl),stat=ierr)
    IF (ierr.EQ.1) STOP 'Error in ESH allocation'

    CALL sub_ESH(ESH_tot,ESH,indblk)
    !CALL destroy(ESH_tot)

    DO i=1,ncont
       CALL concat(ESH(cblk(i),cblk(i)),(-1.d0,0.d0),SelfEneR(i),1,1)
    ENDDO

    !if (debug) then
    !   write(*,*) '----------------------------------------------------'
    !   call writeMemInfo(6)
    !   call writePeakInfo(6)
    !   write(*,*) '----------------------------------------------------'
    !endif

    !Allocazione delle gsmr
    ALLOCATE(gsmr(nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate gsmr'

    !Chiamata di Make_gsmr_mem
    CALL Make_gsmr_mem(ESH,nbl,2)

    !if (debug) then
    !   write(*,*) '----------------------------------------------------'
    !   call writeMemInfo(6)
    !   call writePeakInfo(6)
    !   write(*,*) '----------------------------------------------------'
    !end if

    !Allocazione delle Gr
    ALLOCATE(Gr(nbl,nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Gr'

    CALL Make_Grdiag_mem(ESH,indblk)

    !if (debug) then
    !   write(*,*) '----------------------------------------------------'
    !   call writeMemInfo(6)
    !   call writePeakInfo(6)
    !   write(*,*) '----------------------------------------------------'
    !endif

    !Distruzione delle gsmall
    do i=2,nbl 
       call destroy(gsmr(i))
    enddo
    DEALLOCATE(gsmr)

    !Chiamata di Outer_A_mem
    !CALL Outer_A_mem(Tlc,Tcl,gsurfR,Aout) 

    IF ( (PRESENT (lower_outer)).AND.(lower_outer) ) THEN
       CALL Outer_GreenR_mem(Tlc,Tcl,gsurfR,struct,.TRUE.,Aout) 
    ELSE
       CALL Outer_GreenR_mem(Tlc,Tcl,gsurfR,struct,.FALSE.,Aout) 
    ENDIF

    !if (debug) then
    !   write(*,*) '----------------------------------------------------'
    !   call writePeakInfo(6)
    !   write(*,*) '----------------------------------------------------'
    !endif
    !Chiamata di Make_Spectral_mem
    !CALL Make_Spectral_mem(struct,0,A)
    !if (debug) write(*,*) 'Compute GreenR'   
    !Nota: Il flag 0 indica che non teniamo le Gr calcolate 
    !E' necessario metterlo a 1 se dopo il calcolo della spectral 
    !effettuiamo il calcolo dei contributi aggiunti della Gless
    !(omesso perchè è implementato un driver a parte)
    !CALL Make_GreenR_mem(H%nrow,0,A)
    call Make_GreenR_mem2(ESH_tot,nbl,indblk,0,A)

    !if (debug) then
    !   write(*,*) '----------------------------------------------------'
    !   call writeMemInfo(6)
    !   call writePeakInfo(6)
    !   write(*,*) '----------------------------------------------------'
    !endif

    CALL destroy(ESH_tot)

    !Distruzione dell'array Gr
    DEALLOCATE(Gr)
    DO i=1,nbl
       CALL destroy(ESH(i,i))
    ENDDO
    DO i=2,nbl
       CALL destroy(ESH(i-1,i))
       CALL destroy(ESH(i,i-1))
    ENDDO
    DEALLOCATE(ESH)

    !Concatenazioone di A in Aout

    call concat(Aout,A,1,1)
    call destroy(A)

    !if (debug) then
    !   write(*,*) '----------------------------------------------------'
    !   call writeMemInfo(6)
    !   call writePeakInfo(6)
    !   write(*,*) '----------------------------------------------------'
    !endif

  END SUBROUTINE calls_eq_mem




  !****************************************************************************
  !
  ! Driver for computing Gless contributions due to all contacts but collector
  ! writing on memory 
  !
  !****************************************************************************

  SUBROUTINE calls_neq_mem(H,S,E,SelfEneR,Tlc,Tcl,gsurfR,struct,frm,min,Glout)

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
    !min: collector index 
    !
    !Output:
    !Aout: Gless contributions (Device + Contacts overlap regions -> effective conductor)
    !
    !*****************************************************************************


    IMPLICIT NONE

    !In/Out
    TYPE(z_CSR) :: H,S,Glout
    TYPE(z_CSR), DIMENSION(:) :: SelfEneR, gsurfR, Tlc, Tcl
    COMPLEX(8) :: E
    TYPE(Tstruct_info) :: struct
    REAL(8), DIMENSION(:) :: frm
    INTEGER :: min

    !Work
    INTEGER :: i,ierr,i1,ncont,nbl
    INTEGER, DIMENSION(:), POINTER :: cblk, indblk
    TYPE(z_CSR), DIMENSION(:,:), ALLOCATABLE :: ESH
    TYPE(z_CSR) :: ESH_tot, Gl
    LOGICAL :: destr

    nbl = struct%num_PLs
    ncont = struct%num_conts
    indblk => struct%mat_PL_start
    cblk => struct%cblk

    !if (debug) write(*,*) '----------------------------------------------------'
    !if (debug) call writeMemInfo(6)
    !if (debug) call writePeakInfo(6)
    !if (debug) write(*,*) '----------------------------------------------------'

    !Costruiamo la matrice sparsa ESH
    CALL prealloc_sum(H,S,(-1.d0, 0.d0),E,ESH_tot)

    !Costruiamo l'array di sparse ESH
    !Allocazione dell'array di sparse ESH
    ALLOCATE(ESH(nbl,nbl),stat=ierr)
    IF (ierr.EQ.1) STOP 'Error in ESH allocation'

    CALL sub_ESH(ESH_tot,ESH,indblk)


    DO i=1,ncont
       CALL concat(ESH(cblk(i),cblk(i)),(-1.d0,0.d0),SelfEneR(i),1,1)
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
    CALL Make_gsmr_mem(ESH,nbl,2)

    !Chiamata di Make_gsml_mem solo se i contatti sono più di due o se il 
    !contatto a potenziale minimo è nel primo blocco
    IF (((ncont.GT.2).OR.(cblk(min).EQ.1)).AND.(nbl.gt.2)) THEN

       CALL Make_gsml_mem(ESH,1,nbl-1)    

    ENDIF

    !if (debug) write(*,*) '----------------------------------------------------'
    !if (debug) call writeMemInfo(6)
    !if (debug) call writePeakInfo(6)
    !if (debug) write(*,*) '----------------------------------------------------'

    !Allocazione delle Gr
    ALLOCATE(Gr(nbl,nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Gr'

    CALL Make_Grdiag_mem(ESH,indblk)

    !if (debug) write(*,*) '----------------------------------------------------'
    !if (debug) call writeMemInfo(6)
    !if (debug) call writePeakInfo(6)
    !if (debug) write(*,*) '----------------------------------------------------'


    do i=1,nbl
       destr=.true.

       do i1=1,ncont
          if(i1.eq.min) cycle
          if(i.eq.cblk(i1)) destr=.false.
       enddo

       if (destr) then
          call destroy(Gr(i,i)) 
          if (i.ne.1) call destroy(Gr(i-1,i))
          if (i.ne.nbl) call destroy(Gr(i+1,i))  
       endif
    enddo


    !Chiamata di Make_Grcol_mem per i contatti necessari 
    DO i=1,ncont
       IF (i.NE.min) THEN
          CALL Make_Grcol_mem(ESH,cblk(i),indblk)
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
    CALL Outer_Gl_mem(Tlc,gsurfR,SelfEneR,struct,frm,min,Glout)

    !if (debug) write(*,*) '----------------------------------------------------'
    !if (debug) call writePeakInfo(6)
    !if (debug) write(*,*) '----------------------------------------------------'

    !Calcolo della Gless nel device
    !if (debug) write(*,*) 'Compute Gless'   
    CALL Make_Gl_mem(ESH,SelfEneR,frm,min,0,ESH_tot,Gl)

    CALL destroy(ESH_tot)

    DO i=1,nbl
       DO i1=1,nbl
          IF (ALLOCATED(Gr(i,i1)%rowpnt)) THEN
             CALL destroy(Gr(i,i1))
             !if (debug) WRITE(*,*) 'Deallocating Gr out of Make_Gl_mem'
          ENDIF
       ENDDO
    ENDDO


    !Distruzione dell'array Gr
    DEALLOCATE(Gr)
    DO i=1,nbl
       CALL destroy(ESH(i,i))
    ENDDO
    DO i=2,nbl
       CALL destroy(ESH(i-1,i))
       CALL destroy(ESH(i,i-1))
    ENDDO
    DEALLOCATE(ESH)

    !Concatenazione di Gl in Glout
    CALL concat(Glout,Gl,1,1)
    CALL destroy(Gl)

    !if (debug) write(*,*) '----------------------------------------------------'
    !if (debug) call writeMemInfo(6)
    !if (debug) call writePeakInfo(6)
    !if (debug)  write(*,*) '----------------------------------------------------'

  END SUBROUTINE calls_neq_mem


  !**********************************************************************
  !
  !  Divides sparse matrix ES-H in the device region in a sparse matrices
  !  array ESH(nbl,nbl)  (needs global variable indblk)
  !
  !**********************************************************************

  SUBROUTINE sub_ESH(ESH_tot,ESH,indblk)

    !**********************************************************************
    !Input:
    !ESH_tot: sparse matrix ES-H related to device
    !
    !Output:
    !ESH(nbl,nbl): sparse matrix array -> single matrices allocated 
    !              internally, array ESH(nbl,nbl) allocated externally
    !**********************************************************************

    IMPLICIT NONE 

    INTEGER :: ierr, i
    TYPE(z_CSR) :: ESH_tot
    INTEGER :: nbl
    TYPE(z_CSR), DIMENSION(:,:) :: ESH
    INTEGER, DIMENSION(:) :: indblk

    nbl = size(ESH,1)

    DO i=1,nbl

       CALL zextract(ESH_tot,indblk(i),indblk(i+1)-1,indblk(i),indblk(i+1)-1,ESH(i,i))

    ENDDO

    DO i=2,nbl

       CALL zextract(ESH_tot,indblk(i-1),indblk(i)-1,indblk(i),indblk(i+1)-1,ESH(i-1,i))
       CALL zextract(ESH_tot,indblk(i),indblk(i+1)-1,indblk(i-1),indblk(i)-1,ESH(i,i-1))

    ENDDO

  END SUBROUTINE sub_ESH


  !***********************************************************************
  !
  !  Construct a sparse matrix starting from the sparse matrices array
  !  related to the blocks  
  !
  !***********************************************************************

  SUBROUTINE rebuild(Atot,A,n,indb)

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
    INTEGER :: i,j,Atot_nnz,Atot_nrow,i1,j1

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

  END SUBROUTINE rebuild

  !***********************************************************************
  !
  !  g_small right (gsmr) calculation - write on memory
  !
  !***********************************************************************

  SUBROUTINE Make_gsmr_mem(ESH,sbl,ebl)

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
    TYPE(z_CSR), DIMENSION(:,:) :: ESH
    integer :: sbl,ebl                            ! start block, end block

    !Work
    TYPE(z_CSR) :: work1, work2

    TYPE(z_DNS) :: gsmr_d
    INTEGER :: nrow,nrow_next,nnz
    INTEGER :: i, nbl

    nbl = size(ESH,1)

    if (sbl.lt.ebl) return

    !***
    !gsmr(sbl)
    !***

    if (nbl.gt.1) then

       nrow=ESH(sbl,sbl)%nrow    !indblk(nbl+1)-indblk(nbl)

       CALL create(gsmr_d,nrow,nrow)

#ifdef __PARDISO
       CALL zINV_PARDISO(ESH(sbl,sbl), nrow, gsmr_d%val)
#endif
#ifdef __SUPERLU
       CALL zINV_LU(ESH(sbl,sbl), gsmr_d%val)  
#endif 
#ifdef __LAPACK
       CALL zINV_LAPACK(ESH(sbl,sbl), gsmr_d%val)
#endif

       CALL create( gsmr(sbl),nrow,nrow,nzdrop(gsmr_d,drop) )

       CALL dns2csr(gsmr_d,gsmr(sbl))

       CALL destroy(gsmr_d)

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

       CALL create(gsmr_d,nrow,nrow)

#ifdef __PARDISO
       CALL zINV_PARDISO(work1, nrow, gsmr_d%val)
#endif 
#ifdef __SUPERLU
       CALL zINV_LU(work1, gsmr_d%val)  
#endif 
#ifdef __LAPACK
       CALL zINV_LAPACK(work1, gsmr_d%val)
#endif

       CALL destroy(work1)

       CALL create(gsmr(i),nrow,nrow,nzdrop(gsmr_d,drop) )

       CALL dns2csr(gsmr_d,gsmr(i))

       CALL destroy(gsmr_d)

    ENDDO

    if (debug) then
       WRITE(*,*)
       WRITE(*,*) '********************'
       WRITE(*,*) 'Make_gsmr_mem done'
       WRITE(*,*) '********************'
    endif

  END SUBROUTINE Make_gsmr_mem

 


  !***********************************************************************
  !
  !  g_small left (gsml) calculation - write on memory
  !
  !***********************************************************************

  SUBROUTINE Make_gsml_mem(ESH,sbl,ebl)

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
    TYPE(z_CSR), DIMENSION(:,:) :: ESH
    integer :: sbl,ebl                       ! start block, end block

    !Work
    TYPE(z_CSR) :: work1, work2
    TYPE(z_DNS) :: gsml_d
    INTEGER :: nrow,nrow_prev,nnz
    INTEGER :: i, nbl

    if (sbl.gt.ebl) return

    !***
    !gsml(sbl)
    !***
    nbl = size(ESH,1)
    nrow=ESH(sbl,sbl)%nrow  !indblk(2)-indblk(1)

    CALL create(gsml_d,nrow,nrow)

#ifdef __PARDISO
       CALL zINV_PARDISO(ESH(sbl,sbl), nrow, gsml_d%val)
#endif 
#ifdef __SUPERLU
       CALL zINV_LU(ESH(sbl,sbl), gsml_d%val)  
#endif 
#ifdef __LAPACK
       CALL zINV_LAPACK(ESH(sbl,sbl), gsml_d%val)
#endif

    CALL create(gsml(sbl),nrow,nrow,nzdrop(gsml_d,drop))

    CALL dns2csr(gsml_d,gsml(sbl))

    CALL destroy(gsml_d)

    !***
    !gsml(sbl+1):gsml(ebl)
    !***

    DO i=sbl+1,ebl

       nrow=ESH(i,i)%nrow   !indblk(i+1)-indblk(i)

       CALL prealloc_mult(ESH(i,i-1),gsml(i-1),(-1.d0, 0.d0),work1)

       CALL prealloc_mult(work1,ESH(i-1,i),work2)

       CALL destroy(work1)

       CALL prealloc_sum(ESH(i,i),work2,work1)

       CALL destroy(work2)

       CALL create(gsml_d,nrow,nrow)

#ifdef __PARDISO
       CALL zINV_PARDISO(work1, nrow, gsml_d%val)
#endif 
#ifdef __SUPERLU
       CALL zINV_LU(work1, gsml_d%val)  
#endif 
#ifdef __LAPACK
       CALL zINV_LAPACK(work1, gsml_d%val)
#endif

       CALL destroy(work1)

       CALL create(gsml(i),nrow,nrow,nzdrop(gsml_d,drop))

       CALL dns2csr(gsml_d,gsml(i))

       CALL destroy(gsml_d)

    ENDDO

    if (debug) then
       WRITE(*,*) '********************'
       WRITE(*,*) 'Make_gsml_mem done'
       WRITE(*,*) '********************'
    endif

  END SUBROUTINE Make_gsml_mem





  !***********************************************************************
  !
  !  Diagonal, Subdiagonal, Superdiagonal blocks of Green Retarded 
  !  Gr(nbl,nbl) - writing on memory
  !
  !*********************************************************************** 

  SUBROUTINE Make_Grdiag_mem(ESH,indblk)

    !***********************************************************************
    !Input:
    !ESH: sparse matrices array ESH(nbl,nbl)
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
    TYPE(z_CSR), DIMENSION(:,:) :: ESH
    INTEGER, DIMENSION(:), POINTER :: indblk

    !Work
    INTEGER :: i,nrow,nrow_prev,nrow_next,nnz,nbl
    TYPE(z_CSR) :: work1, work2, work3
    TYPE(z_DNS) :: Gr_d

    !***
    !Gr(1,1)
    !***

    nbl = size(ESH,1)
    nrow=indblk(2)-indblk(1)  

    if(nbl.gt.1) then

       CALL prealloc_mult(ESH(1,2),gsmr(2),work1)

       CALL prealloc_mult(work1,ESH(2,1),work2)

       CALL destroy(work1)

       CALL prealloc_sum(ESH(1,1),work2,(-1.d0, 0.d0),work1)

       CALL destroy(work2)

    else

       CALL clone(ESH(1,1), work1)

    endif

    CALL create(Gr_d,nrow,nrow)

#ifdef __PARDISO
       CALL zINV_PARDISO(work1, nrow, Gr_d%val)
#endif 
#ifdef __SUPERLU
       CALL zINV_LU(work1, Gr_d%val)  
#endif 
#ifdef __LAPACK
       CALL zINV_LAPACK(work1, Gr_d%val)
#endif

    CALL destroy(work1)

    CALL create(Gr(1,1),nrow,nrow,nzdrop(Gr_d,drop))

    CALL dns2csr(Gr_d,Gr(1,1))

    CALL destroy(Gr_d)

    !***
    !Diagonal, Subdiagonal and Superdiagonal blocks
    !***

    DO i=2,nbl

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

  END SUBROUTINE Make_Grdiag_mem





  !**************************************************************************
  !
  !  Calculate Green Retarded column "n" - writing on memory 
  !
  !**************************************************************************

  SUBROUTINE Make_Grcol_mem(ESH,n,indblk)

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
    TYPE(z_CSR), DIMENSION(:,:) :: ESH
    INTEGER :: n
    INTEGER, DIMENSION(:), POINTER :: indblk 

    !Work
    INTEGER :: i,j,nrow,nrow_prev,nrow_next,ncol,nbl
    TYPE(z_CSR) :: work1
    REAL(8) :: max

    nbl = size(ESH,1)

    IF (n.GT.nbl) THEN
       STOP 'Error in Make_Grcol_mem : n is greater than nbl'
    ENDIF

    !***
    !Downgoing (n<nbl-1)
    !**
    ncol=indblk(n+1)-indblk(n)

    IF (n.LT.(nbl-1)) THEN

       DO i=n+2,nbl

          nrow=indblk(i+1)-indblk(i)

          IF (Gr(i-1,n)%nnz.NE.0) THEN 
             max=MAXVAL(ABS(Gr(i-1,n)%nzval(:)))
          ELSE 
             max=0
          ENDIF

          IF (max.GT.drop) THEN

             CALL prealloc_mult(gsmr(i),ESH(i,i-1),(-1.d0, 0.d0),work1)
             CALL prealloc_mult(work1,Gr(i-1,n),Gr(i,n))
             CALL destroy(work1)

          ELSE

             CALL create(Gr(i,n),nrow,ncol,0)
             Gr(i,n)%rowpnt(:)=1

          ENDIF

       ENDDO

    ENDIF

    !***
    !Upgoing (n>2)
    !***

    IF (n.GT.2) THEN

       DO i=n-2,1,(-1)

          nrow=indblk(i+1)-indblk(i)

          IF (Gr(i+1,n)%nnz.NE.0) THEN 
             max=MAXVAL(ABS(Gr(i+1,n)%nzval(:)))
          ELSE 
             max=0
          ENDIF

          IF (max.GT.drop) THEN

             CALL prealloc_mult(gsml(i),ESH(i,i+1),(-1.d0, 0.d0),work1)
             CALL prealloc_mult(work1,Gr(i+1,n),Gr(i,n))
             CALL destroy(work1)

          ELSE

             CALL create(Gr(i,n),nrow,ncol,0)

             Gr(i,n)%rowpnt(:)=1

          ENDIF

       ENDDO

    ENDIF

    if (debug) then
       WRITE(*,*) '********************'
       WRITE(*,*) 'Make_Grcol_mem done column',n
       WRITE(*,*) '********************'
    endif

  END SUBROUTINE Make_Grcol_mem




  !****************************************************************************
  !
  !  Calculate Spectral Density - writing on memory
  !
  !****************************************************************************

  SUBROUTINE Make_Spectral_mem(struct,keep_Gr,A)

    !****************************************************************************
    !Input:
    !nrow_tot: total device matrix rows
    !keep_Gr (flag): if 0, destroy elements Gr(:,:) when no more necessary. If 1,
    !                keep them in memory
    !
    !global variable needed: nbl, indblk(nbl+1), Gr(:,:) 
    !
    !Output:
    !A: sparse matrix containing spectral density of device (allocated internally)
    !****************************************************************************

    IMPLICIT NONE

    !In/Out
    TYPE(z_CSR) :: A
    TYPE(Tstruct_info) :: struct
    INTEGER :: keep_Gr

    !Work
    INTEGER :: i,nrow,nrow_prev,i1,j1,ierr,nrow_tot,nbl
    TYPE(z_CSR), DIMENSION(:,:), ALLOCATABLE :: Asub
    INTEGER, DIMENSION(:), POINTER :: indblk

    nbl = struct%num_PLs
    indblk => struct%mat_PL_start
    nrow_tot = struct%total_dim


    IF ((keep_Gr.NE.1).AND.(keep_Gr.NE.0)) STOP 'keep_Gr must be 0 or 1'

    !Allocazione dell'array di sparse Asub e di A
    ALLOCATE(Asub(nbl,nbl),stat=ierr)
    IF (ierr.NE.0) THEN
       STOP 'ALLOCATION ERROR: could not allocate Asub(nbl,nbl)'
    ENDIF

    CALL create(A,nrow_tot,nrow_tot,0)
    A%rowpnt(:)=1

    !***
    !A(1,1) 
    !***
    nrow=indblk(2)-indblk(1)

    CALL zspectral(Gr(1,1),Gr(1,1),0,Asub(1,1))

    IF (keep_Gr.EQ.0) CALL destroy(Gr(1,1))

    CALL concat(A,Asub(1,1),1,1)
    CALL destroy(Asub(1,1))

    !***
    !Diagonal, Subdiagonal and Super diagonal blocks
    !***

    DO i=2,nbl

       nrow=indblk(i+1)-indblk(i)

       i1=indblk(i)
       j1=indblk(i)

       CALL zspectral(Gr(i,i),Gr(i,i),0,Asub(i,i))

       IF (keep_Gr.EQ.0) CALL destroy(Gr(i,i))

       CALL concat(A,Asub(i,i),i1,j1)
       CALL destroy(Asub(i,i))

       i1=indblk(i-1)
       j1=indblk(i)

       CALL zspectral(Gr(i-1,i),Gr(i,i-1),0,Asub(i-1,i))
       CALL concat(A,Asub(i-1,i),i1,j1)
       CALL destroy(Asub(i-1,i))  

       i1=indblk(i)
       j1=indblk(i-1)

       CALL zspectral(Gr(i,i-1),Gr(i-1,i),0,Asub(i,i-1))

       IF (keep_Gr.EQ.0) THEN 
          CALL destroy(Gr(i,i-1))
          CALL destroy(Gr(i-1,i))
       ENDIF

       !Blocks concatenation
       CALL concat(A,Asub(i,i-1),i1,j1)
       CALL destroy(Asub(i,i-1))

    ENDDO

    DEALLOCATE(Asub)

    if (debug) then    
       WRITE(*,*) '**********************'
       WRITE(*,*) 'Make_Spectral_mem done'
       WRITE(*,*) '**********************'
    endif

  END SUBROUTINE Make_Spectral_mem





  !****************************************************************************
  !
  !  Calculate Green Retarded - writing on memory
  !
  !****************************************************************************

!!$  SUBROUTINE Make_GreenR_mem(nrow_tot,keep_Gr,A)
!!$
!!$    !****************************************************************************
!!$    !Input:
!!$    !nrow_tot: total device matrix rows
!!$    !keep_Gr (flag): if 0, destroy elements Gr(:,:) when no more necessary. If 1,
!!$    !                keep them in memory
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
!!$    INTEGER :: nrow_tot, keep_Gr
!!$    TYPE(z_CSR) :: A
!!$
!!$    !TYPE(z_CSR), DIMENSION(nbl,nbl) :: ESH
!!$
!!$    !Work
!!$    INTEGER :: i,nrow,i1,j1,ierr
!!$
!!$    IF ((keep_Gr.NE.1).AND.(keep_Gr.NE.0)) STOP 'keep_Gr must be 0 or 1'
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
!!$    IF (keep_Gr.EQ.0) CALL destroy(Gr(1,1))
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
!!$       IF (keep_Gr.EQ.0) CALL destroy(Gr(i,i))
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
!!$       IF (keep_Gr.EQ.0) THEN 
!!$          CALL destroy(Gr(i,i-1))
!!$          CALL destroy(Gr(i-1,i))
!!$       ENDIF
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

  SUBROUTINE Make_GreenR_mem2(P,nbl,indblk,keep_Gr,A)

    !****************************************************************************
    !Input:
    !P: CSR matrix containing masking pattern
    !keep_Gr (flag): if 0, destroy elements Gr(:,:) when no more necessary. If 1,
    !                keep them in memory
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
    INTEGER :: keep_Gr
    TYPE(z_CSR) :: A, P

    !Work
    INTEGER :: i, j, j1, i1, ix, iy, x, y, col, oldx

    IF ((keep_Gr.NE.1).AND.(keep_Gr.NE.0)) STOP 'keep_Gr must be 0 or 1'

    !create A with same pattern of P
    CALL create(A,P%nrow,P%ncol,P%nnz)
    A%rowpnt(:)=P%rowpnt(:)
    A%colind(:)=P%colind(:)


    !If only one block is present, concatenation is not needed 
    !and it's implemented in a more trivial way
    IF (nbl.EQ.1) THEN

       CALL mask(Gr(1,1),P,A)

       CALL destroy(Gr(1,1))

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
                if (y.eq.0) then
                   write(*,*) 'ERROR in Make_GreenR_mem2: probably wrong PL size',x
                   stop
                endif
             ENDIF
             col = A%colind(j) - indblk(y) + 1
             A%nzval(j) = (0.0, 0.0)
             DO j1 = Gr(x,y)%rowpnt(i1), Gr(x,y)%rowpnt(i1 + 1) -1
                IF (Gr(x,y)%colind(j1).EQ.col)  then
                   A%nzval(j) = Gr(x,y)%nzval(j1)
                   exit
                ENDIF
             ENDDO
          ENDDO


          IF(oldx.NE.x) THEN

             IF ( (keep_Gr.EQ.0).AND.(x.eq.2) ) THEN 
                CALL destroy(Gr(x-1,x-1))
                CALL destroy(Gr(x-1,x))
             ELSEIF ( (keep_Gr.EQ.0).AND.(x.GT.2) ) THEN 
                CALL destroy(Gr(x-1,x-1))
                CALL destroy(Gr(x-1,x-2))
                CALL destroy(Gr(x-1,x))
             ENDIF


          ENDIF


       ENDDO

       CALL destroy(Gr(nbl, nbl))
       CALL destroy(Gr(nbl, nbl - 1))


    ENDIF

    !if (debug) call writePeakInfo(6)    
    if (debug) then
       WRITE(*,*) '**********************'
       WRITE(*,*) 'Make_GreenR_mem2 done'
       WRITE(*,*) '**********************'
    endif

  END SUBROUTINE Make_GreenR_mem2





  !****************************************************************************
  !
  ! Calculate Gless contributions for all contacts (except collector, defined
  ! as contact with minimum potential) - writing on memory
  !
  !****************************************************************************

  SUBROUTINE Make_Gl_mem(ESH,SelfEneR,frm,min,keep_Gr,P,Gl)

    !******************************************************************************
    !Input:  
    !ESH(nbl,nbl): sparse matrices array ES-H
    !SelfEneR(ncont): sparse matrices array with contacts Self Energy 
    !frm(ncont): Fermi diistribution value for each contact 
    !min:  collector index
    !keep_Gr (flag): if 0, destroy elements Gr(:,:) when no more necessary. If 1,
    !                keep them in memory 
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
    TYPE(z_CSR), DIMENSION(:,:) :: ESH
    TYPE(z_CSR), DIMENSION(:) :: SelfEneR
    TYPE(Tstruct_info) :: struct
    REAL(8), DIMENSION(:) :: frm
    INTEGER :: min,keep_Gr

    !Work
    TYPE(z_CSR), DIMENSION(:,:), ALLOCATABLE :: Glsub
    TYPE(z_CSR) :: Gam, P
    TYPE(z_CSR) :: work1,Ga
    INTEGER :: ierr,i,j,nrow,i1,j1,cb,nrow_tot,nrow_cb, nbl, ncont
    INTEGER :: oldx, col, iy, ix, x, y, ii, jj
    INTEGER, DIMENSION(:), POINTER :: indblk, cblk
    COMPLEX(dp) :: frmdiff

    ncont = struct%num_conts    
    nbl = struct%num_PLs
    cblk => struct%cblk
    indblk => struct%mat_PL_start

    IF ((keep_Gr.NE.1).AND.(keep_Gr.NE.0)) STOP 'keep_Gr must be 0 or 1'

    !Allocazione degli array di sparse
    ALLOCATE(Glsub(nbl,nbl),stat=ierr)
    IF (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Glsub(nbl,nbl)'


    !create A with same pattern of P
    CALL create(Gl,P%nrow,P%ncol,P%nnz)
    Gl%rowpnt(:)=P%rowpnt(:)
    Gl%colind(:)=P%colind(:)
    Gl%nzval(:) = (0.0 , 0.0) 


    !***
    !Iterazione sui contatti
    !***

    DO j=1,ncont

       IF (j.NE.min.AND.ABS(frm(j)-frm(min)).GT.drop) THEN

          !Calcolo della Gamma corrispondente
          cb=cblk(j)
          !nota:cb indica l'indice di blocco corrispondente al contatto j-esimo

          CALL zspectral(SelfEneR(j),SelfEneR(j),0,Gam)

          !***
          !Iterazione sui blocchi 
          !***
          frmdiff = cmplx(frm(j)-frm(min))
          !Calcolo del sottoblocco Gl(1,1) fuori iterazione
          nrow=indblk(2)-indblk(1)
          CALL prealloc_mult(Gr(1,cb),Gam,work1)
          CALL zdagacsr(Gr(1,cb),Ga)

          !G(1,j) va tenuta per la prima iterazione blocco sopradiagonale

          CALL prealloc_mult(work1,Ga,frmdiff,Glsub(1,1))
          CALL destroy(work1)
          CALL destroy(Ga)

          i1=indblk(1)
          j1=indblk(1)
          CALL zmask_realloc(Glsub(1,1), ESH(1,1))
          !CALL concat(Gl,Glsub(1,1),i1,j1)
          !CALL destroy(Glsub(1,1))

          !Calcolo sottoblocchi diagonali, sopradiagonali e sottodiagonali
          DO i=2,nbl

             nrow=indblk(i+1)-indblk(i)

             !Calcolo blocchi sopradiagonali
             CALL prealloc_mult(Gr(i-1,cb),Gam,work1)

             CALL zdagacsr(Gr(i,cb),Ga)

             CALL prealloc_mult(work1,Ga,frmdiff,Glsub(i-1,i))
             CALL zmask_realloc(Glsub(i-1,i), ESH(i-1,i))
             CALL destroy(work1)
             !Teniamo Ga per il calcolo del blocco diagonale

             i1=indblk(i-1)
             j1=indblk(i)
             !CALL concat(Gl,Glsub(i-1,i),i1,j1)

             !CALL destroy(Glsub(i-1,i))

             !Calcolo blocchi diagonali
             CALL prealloc_mult(Gr(i,cb),Gam,work1)

             CALL prealloc_mult(work1,Ga,frmdiff,Glsub(i,i))
             CALL zmask_realloc(Glsub(i,i), ESH(i,i))
             CALL destroy(work1)
             call destroy(Ga)

             i1=indblk(i)
             j1=indblk(i)
             !CALL concat(Gl,Glsub(i,i),i1,j1)
             !CALL destroy(Glsub(i,i))
             !Calcolo diretto blocchi sottodiagonali
             CALL prealloc_mult(Gr(i,cb),Gam,work1)
             CALL zdagacsr(Gr(i-1,cb),Ga)
             CALL prealloc_mult(work1,Ga,frmdiff,Glsub(i,i-1))
             CALL zmask_realloc(Glsub(i,i-1), ESH(i,i-1))

             CALL destroy(work1)
             call destroy(Ga)

             !Concateniamo Glsub come sottodiagonale
             i1=indblk(i)
             j1=indblk(i-1)
             !CALL concat(Gl,Glsub(i,i-1),i1,j1)
             !CALL destroy(Glsub(i,i-1))

             IF (keep_Gr.EQ.0) CALL destroy(Gr(i-1,cb)) 

          ENDDO

          IF (keep_Gr.EQ.0) CALL destroy(Gr(nbl,cb))                     

          call destroy(Gam)




          !Concatenation for every contact in Gless. Performs a sum on elements, not a replacement
          !Similar to Make_GreenR_mem2, except for sum of elements
          !Note: to backup old version zconcat calls (and Glsub deallocations) must be uncommented and all this part removed 

          !If only one block is present, concatenation is not needed and it's implemented in a more trivial way
          IF (nbl.EQ.1) THEN

             CALL zmask_realloc(Glsub(1,1), ESH(1,1))

             CALL concat(Gl,Glsub(1,1),1,1)

             CALL destroy(Glsub(1,1))

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
                   WRITE(*,*) 'Error in Make_GreenR_mem2: could not find row block of index ',ii,ix
                   STOP
                ENDIF

                !Offset: i1 is the index for separate blocks
                i1 = ii - indblk(x) + 1

                !Cycle upon columns 
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
                      IF (y.EQ.0) THEN
                         WRITE(*,*) 'ERROR in Make_Gl_mem: probably wrong PL size', x
                         STOP
                      ENDIF
                   ENDIF
                   col = Gl%colind(jj) - indblk(y) + 1
                   DO j1 = Glsub(x,y)%rowpnt(i1), Glsub(x,y)%rowpnt(i1 + 1) -1
                      IF (Glsub(x,y)%colind(j1).EQ.col)  then
                         Gl%nzval(jj) = Gl%nzval(jj) + Glsub(x,y)%nzval(j1)
                         exit
                      ENDIF
                   ENDDO
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

       ENDIF

    ENDDO


    if (debug) then
       WRITE(*,*) '********************'
       WRITE(*,*) 'Make_Gl_mem done'
       WRITE(*,*) '********************'
    endif

  END SUBROUTINE Make_Gl_mem




  !****************************************************************************
  !
  !  Calculate Gless contributions for all contacts but collector, in the 
  !  contacts regions, where overlap with device orbitals is non-zero
  !  writing on memory
  !
  !****************************************************************************

  SUBROUTINE Outer_A_mem(Tlc,Tcl,gsurfR,struct,Aout)

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

    IMPLICIT NONE 

    !In/Out
    TYPE(z_CSR), DIMENSION(:) :: Tlc,Tcl,gsurfR
    TYPE(Tstruct_Info), intent(in) :: struct
    TYPE(z_CSR), intent(out) :: Aout

    !Work
    TYPE(z_CSR) :: work1,Grlc,Grcl,Asub
    INTEGER :: i,cb,ierr,nrow_tot,i1,j1
    INTEGER :: ncont,nbl
    INTEGER, DIMENSION(:), POINTER :: indblk, cblk

    ncont = struct%num_conts
    nbl = struct%num_PLs
    nrow_tot = struct%total_dim
    indblk => struct%mat_PL_start
    cblk => struct%cblk

    !Righe totali del conduttore effettivo nrow_tot 
    !nrow_tot=indblk(nbl+1)-1
    !DO i=1,ncont
    !   nrow_tot=nrow_tot+ncdim(i)   !gsurfR(i)%nrow
    !ENDDO
    CALL create(Aout,nrow_tot,nrow_tot,0)
    Aout%rowpnt(:)=1

    DO i=1,ncont

       !Numero di blocco del contatto
       cb=cblk(i)

       !Calcolo di Grlc
       CALL prealloc_mult(Gr(cb,cb),Tlc(i),(-1.d0, 0.d0),work1)
       !Nota: numero colonne di Tlc = numero di righe di gsurf(i)

       CALL prealloc_mult(work1,gsurfR(i),Grlc)
       CALL destroy(work1)

       !Calcolo di Grcl
       CALL prealloc_mult(gsurfR(i),Tcl(i),(-1.d0, 0.d0),work1)

       CALL prealloc_mult(work1,Gr(cb,cb),Grcl)
       CALL destroy(work1)

       !Calcolo della spectral density
       CALL zspectral(Grlc,Grcl,0,Asub)

       !Dellocazione delle Green Retarded corrispondenti
       CALL destroy(Grlc)
       CALL destroy(Grcl)

       !Concatenazione di Asub nella posizione corrispondente
       i1=indblk(cb)
       j1=struct%mat_B_start(i)-struct%central_dim+indblk(nbl+1)-1
       CALL concat(Aout,Asub,i1,j1)
       CALL destroy(Asub)  

    ENDDO

    if (debug) then
       WRITE(*,*) '********************'
       WRITE(*,*) 'Outer_A_mem done'
       WRITE(*,*) '********************'
    endif

  END SUBROUTINE Outer_A_mem



  !****************************************************************************
  !
  !  Calculate Green Retarded in the 
  !  contacts regions, where overlap with device orbitals is non-zero
  !  (only the upper part, needed in both K=0 and K points calculations,
  !   or both upper and lower parts)
  !  writing on memory
  !
  !****************************************************************************

  SUBROUTINE Outer_GreenR_mem(Tlc,Tcl,gsurfR,struct,lower,Aout)

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
    TYPE(z_CSR), DIMENSION(:) :: Tlc,Tcl,gsurfR
    LOGICAL :: lower
    TYPE(Tstruct_info), intent(in) :: struct
    TYPE(z_CSR), intent(out) :: Aout


    !Work
    TYPE(z_CSR) :: work1,Grlc,Grcl,Asub
    INTEGER :: i,cb,ierr,nrow_tot,i1,j1
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

       !Calcolo di Grlc
       CALL prealloc_mult(Gr(cb,cb),Tlc(i),(-1.d0, 0.d0),work1)

       CALL prealloc_mult(work1,gsurfR(i),Grlc)

       CALL zmask_realloc(Grlc, Tlc(i))

       CALL destroy(work1)

       !Concatenazione di Asub nella posizione corrispondente
       i1=indblk(cb)
       j1=struct%mat_B_start(i)-struct%central_dim+indblk(nbl+1)-1

       CALL concat(Aout,Grlc,i1,j1)
       !CALL destroy(Asub) 
       CALL destroy(Grlc)
       !CALL destroy(Grcl)

       IF (lower) THEN

          CALL prealloc_mult(gsurfR(i),Tcl(i),(-1.d0, 0.d0), work1)

          CALL prealloc_mult(work1, Gr(cb,cb), Grcl)

          CALL zmask_realloc(Grcl, Tcl(i))

          CALL destroy(work1)

          i1 = struct%mat_B_start(i)-struct%central_dim+indblk(nbl+1)-1
          j1 = indblk(cb)

          CALL concat(Aout,Grcl,i1,j1)

          CALL destroy(Grcl)

       ENDIF


    ENDDO

    if (debug) then
       WRITE(*,*) '********************'
       WRITE(*,*) 'Outer_GreenR_mem done'
       WRITE(*,*) '********************'
    endif

  END SUBROUTINE Outer_GreenR_mem



  !****************************************************************************
  !
  !  Calculate Gless contributions for all contacts except collector, in the
  !  outer region where contacts-device overlap is non-zero - writing on 
  !  memory
  !
  !****************************************************************************

  SUBROUTINE Outer_Gl_mem(Tlc,gsurfR,SelfEneR,struct,frm,min,Glout)

    !****************************************************************************
    !Input:
    !Tlc: sparse matrices array containing contacts-device interaction blocks (ES-H)
    !gsurfR: sparse matrices array containing contacts surface green
    !SelfEneR: sparse matrices array containing contacts Self Energy
    !frm: array containing Fermi distribution values for all contacts
    !min: collector index
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
    TYPE(z_CSR), DIMENSION(:) :: Tlc, gsurfR, SelfEneR
    REAL(8), DIMENSION(:) :: frm
    TYPE(Tstruct_info), intent(in) :: struct
    INTEGER :: min  
    TYPE(z_CSR) :: Glout

    !Work
    TYPE(z_CSR) :: work1, work2, work3, Gam, gsurfA, Ga, Glsub
    INTEGER :: j,k,cb,cbj,i1,j1,nrow_tot,i,nrow_cb,nrow_cbj
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
    CALL create(Glout,nrow_tot,nrow_tot,0)
    Glout%rowpnt(:)=1

    !***
    !Iterazione su tutti i contatti "k" 
    !***
    DO k=1,ncont

       !Esegue le operazioni relative al contatto solo se è valida la condizione
       !sulle distribuzioni di Fermi e se non si tratta del contatto a potenziale minimo
       IF ((ABS(frm(k)-frm(min)).GT.drop).AND.(k.NE.min)) THEN

          !Calcolo della Gamma corrispondente
          cb=cblk(k)
          !nota:cb indica l'indice di blocco corrispondente al contatto j-esimo
          !nrow_cont=ncdim(k) !gsurfR(k)%nrow          
          !nrow_cb=indblk(cb+1)-indblk(cb)

          CALL zspectral(SelfEneR(k),SelfEneR(k),0,Gam)

          !***
          !Calcolo del contributo sulla proria regione
          !***          
          !Primo addendo
          frmdiff=cmplx(frm(min)-frm(k))

          !work1=j(gsurfR-gsurfA)
          CALL zspectral(gsurfR(k),gsurfR(k),0,work1)

          !work2=Tlc*work1=Tlc*j(gsurfR-gsurfA)
          CALL prealloc_mult(Tlc(k),work1,work2)
          CALL destroy(work1)

          !work1=-Gr(cb,cb)*work2=-Gr(cb,cb)*Tlc*j(gsurfR-gsurfA)
          CALL prealloc_mult(Gr(cb,cb),work2,frmdiff,work1)
          CALL destroy(work2)

          !if (debug) write(*,*) 'cont',k,'primo addendo ok'
          !Secondo addendo

          !work2=Tlc*gsurfA
          CALL zdagacsr(gsurfR(k),gsurfA)
          CALL prealloc_mult(Tlc(k),gsurfA,work2)
          CALL destroy(gsurfA)

          !work3=Ga*work2=Ga*Tlc*gsurfA           
          CALL zdagacsr(Gr(cb,cb),Ga)
          CALL prealloc_mult(Ga,work2,work3)

          CALL destroy(Ga)
          CALL destroy(work2)

          !work2=Gam*work3=Gam*Ga*Tlc*gsurfA          
          CALL prealloc_mult(Gam,work3,work2)
          CALL destroy(work3)

          !work3=-Gr*work2=-Gr*Gam*Ga*Tlc*gsurfA        
          CALL prealloc_mult(Gr(cb,cb),work2,frmdiff,work3)
          CALL destroy(work2)

          !if (debug) write(*,*) 'cont',k,'secondo addendo ok'

          !Contributo totale sulla propria regione
          CALL prealloc_sum(work3,work1,Glsub)
          CALL destroy(work1)
          CALL destroy(work3)

          CALL zmask_realloc(Glsub, Tlc(k)) 

          !Concatenazione di Glsub nella matrice globale Glout
          i1=indblk(cb)
          j1=struct%mat_B_start(k)-struct%central_dim+indblk(nbl+1)-1
          CALL concat(Glout,Glsub,i1,j1)
          CALL destroy(Glsub)

          !if (debug) write(*,*) 'cont',k,'concat ok'
          !***
          !Ciclo per il calcolo dei contributi sulle regioni degli altri contatti
          !***

          DO j=1,ncont

             cbj=cblk(j)
             !nrow_cbj=indblk(cbj+1)-indblk(cbj)
             !nrow_contj=ncdim(j) !gsurfR(j)%nrow
             !Esegue le operazioni del ciclo solo se il j.ne.k o se
             !il blocco colonna di Gr è non nullo (altrimenti il contributo è nullo) 

             IF ((j.NE.k).AND.(Gr(cbj,cb)%nnz.NE.0)) THEN

                !work1=Tlc*gsurfA  
                CALL zdagacsr(gsurfR(j),gsurfA)
                CALL prealloc_mult(Tlc(j),gsurfA,work1)
                CALL destroy(gsurfA)

                !if (debug) write(*,*) 'cont',j,'T*g'

                !work2=Ga*work1=Ga*Tlc*gsurfA  
                CALL zdagacsr(Gr(cbj,cb),Ga)
                CALL prealloc_mult(Ga,work1,work2)

                !if (debug) write(*,*) 'cont',j,'Ga*T*g'

                CALL destroy(Ga)
                CALL destroy(work1)

                !if (debug) write(*,*) 'Gam',Gam%nrow,Gam%ncol,Gam%nnz
                !if (debug) write(*,*) 'work2',work2%nrow,work2%ncol,work2%nnz
                !work1=Gam*work2=Gam*Ga*Tlc*gsurfA  
                CALL prealloc_mult(Gam,work2,work1)
                CALL destroy(work2)

                !if (debug) write(*,*) 'cont',j,'Gam*Ga*T*g'

                !if (debug) write(*,*) 'work1',work1%nrow,work1%ncol,work1%nnz
                !if (debug) write(*,*) 'Gr',Gr(cbj,cb)%nrow,Gr(cbj,cb)%ncol,Gr(cbj,cb)%nnz

                !Glsub=-Gr*work1=-Gr*Gam*Ga*Tlc*gsurfA  
                CALL prealloc_mult(Gr(cbj,cb),work1,frmdiff,Glsub)
                CALL destroy(work1)

                CALL zmask_realloc(Glsub, Tlc(j))

                !if (debug) write(*,*) 'cont',j,'Gr*Gam*Ga*T*g'

                !Concatenazione di Glsub nella posizione corrispondente al contatto "j"
                i1=indblk(cbj)
                j1=struct%mat_B_start(j)-struct%central_dim+indblk(nbl+1)-1
                CALL concat(Glout,Glsub,i1,j1)
                CALL destroy(Glsub)

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

  END SUBROUTINE Outer_Gl_mem







  !---------------------------------------------------

END MODULE iterative
