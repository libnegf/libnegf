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


!*************************************************************************
!  Copyright (c) 2004 by Univ. Rome 'Tor Vergata'. All rights reserved.  *  
!  Authors: A. Pecchia, L. Latessa, A. Di Carlo                          *
!                                                                        *
!  Permission is hereby granted to use or copy this program for any      *
!  purpose, provided the above notices are retained on all copies, and   *
!  a notice that the code was modifided is included with the above       *
!  copyright.                                                            *
!  Permission to redistribute the code to third parties is restricted    *
!  by the licence agreement.                                             *
!*************************************************************************
 module negf_greendftb

  use precision
  use constants 
  use negf_mpi_globals
  use allocation
  use negf_lib_param
  use negf_mat_def
  use negf_sparsekit_drv
  use structure, only : TStruct_Info
  use negf_contselfenergy
  use negf_iterative
  use fermi_dist
  use negf_clock

  
  implicit none
  private
  
  integer, PARAMETER :: VBT=80
  LOGICAL, PARAMETER :: mem=.true. 
  LOGICAL, PARAMETER :: profile=.false.
  

  public :: contour_int


contains
  
  subroutine contour_int(H,S,pnegf,struct,DensMat,EnMat)     

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !% Subroutine GreenDftb ver. 5.0
    !% The DFTB hamiltonian is used to calculate the density matrix of
    !% the device. The charge is used in the self-constent cycle (SCF).
    !% Completely rewritten for sparse matrices: by A. Pecchia 21/4/2006
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !% Input arrays and variables
    type(z_CSR), intent(in) :: H           ! hamiltonian
    type(z_CSR), intent(in) :: S	   ! overlap
    type(Tnegf), pointer :: pnegf
    type(TStruct_Info), intent(in) :: struct 

    !% Output arrays
    type(z_CSR), intent(inout) :: DensMat              ! Density Matrix
    type(z_CSR), intent(inout) :: EnMat                ! Energy wighted Density Matrix
    
    
    !% Local Array and variables    
    type(z_CSR) :: HM,SM                               ! molecular HAM and OVR
    type(z_CSR) :: TM(MAXNCONT),ST(MAXNCONT)           ! contact-mol HAM and OVR
    type(z_CSR) :: GS(MAXNCONT)                        ! surface Green functions
    type(z_CSR) :: SelfEneR(MAXNCONT)                  ! 
    type(z_CSR) :: Tlc(MAXNCONT),Tcl(MAXNCONT)         ! Temporary matrices E ST - TM
    type(z_CSR) :: GreenR                              ! Green's Function
    type(z_CSR) :: TpMt                                ! Temporary matrix    
#ifdef MPI
    real(dp), dimension(:), allocatable :: Mat_st      ! MPI mat storage
    real(dp), dimension(:), allocatable :: Mat_rcv     ! MPI mat storage
#endif
    
    type(z_DNS) :: HC_d(MAXNCONT),SC_d(MAXNCONT),GS_d(MAXNCONT),M_d           ! Contact H and S

    real(dp), DIMENSION(:), ALLOCATABLE :: wght,pnts   ! Gauss-quadrature points
    real(dp), DIMENSION(:), ALLOCATABLE :: t_wght,t_pnts ! Gauss-quadrature points

    integer :: cstart(MAXNCONT),cend(MAXNCONT),nmdim        ! Structure specifications
    integer :: ncdim(MAXNCONT)
    integer :: err, tmp,ncol, nc, npid, nn, ncont           ! Assistence variables
    integer :: i, k, l, i1, i2, j1, j2                      ! Assistence variables
    integer :: imin, imax, NumPoles, istart, iend, npc      ! Integration counters
    integer :: nc_vec(1), ibsize, npT, Np(4), N_omega, verbose
    real(dp) :: ncyc, avncyc                           ! Average num. decimations
    real(dp) :: Efermi(MAXNCONT), frm_f(MAXNCONT), mu(MAXNCONT)
    real(dp) :: E, dd, mumin, mumax, muref                
    real(dp) :: Rad, Centre, Lambda, Omega, Temp, Elow
    real(dp) :: c1,c2,T,teta,dt,alpha
    complex(kind=dp) :: Ec,Pc,z1,z2,z_diff,zt                  ! Integration variables

    INTEGER :: t1g,t2g,crg,cmg
    integer, parameter :: low = 2

    ! Trasferimento variabili locali dal contenitore
    ! ------------------------------------------------------------------------
    verbose = pnegf%verbose
    Np = pnegf%Np
    NumPoles = pnegf%n_poles
    Temp = pnegf%Temp
    Efermi = pnegf%Efermi
    mu = pnegf%mu
    N_omega = pnegf%N_omega
    Elow = pnegf%Emin

    ncont = struct%num_conts
    do i=1,ncont
      cstart(i) = struct%mat_B_start(i)
      cend(i)   = struct%mat_C_end(i)
      ncdim(i)  = cend(i)-cstart(i)+1

      pnegf%contdim(i) = ncdim(i)
      pnegf%surfdim(i) = struct%mat_C_start(i) - struct%mat_B_start(i)

      write(*,*) '(int)',i,cstart(i),cend(i),pnegf%surfdim(i),ncdim(i)
    enddo
    nmdim = struct%central_dim
  
    write(*,*) '(int) nmdim=',nmdim
    
    ! ------------------------------------------------------------------------
    if (Temp.ne.0.d0) then

      if (NumPoles.ne.0) then
        Lambda = 2.d0*NumPoles*Kb*Temp*pi
      else
        Lambda = 0.5d0*Kb*Temp*pi
      end if

    else   
      NumPoles = 0
      Lambda = 0.d0
    end if

    !**************************
    ! COMPUTE DENSITY MATRIX
    !**************************
    if (id0.and.verbose.gt.60) then
      write(*,'(73("="))')
      write(*,*) '          COMPUTING DENSITY MATRIX - Complex Integration'
      write(*,'(73("="))') 
    endif

    !-------------------------------------------------------
    ! Separates blocks of H and S
    !------------------------------------------------------- 
    call zextract(H,1,nmdim,1,nmdim,HM)
    call zextract(S,1,nmdim,1,nmdim,SM)

    ! -------------------------------------------------------------
    !  Extract Device-Contact blocks
    ! -------------------------------------------------------------

    do i=1,ncont
       print*, '(int) extract contact',i
       call zextract_dns(H,cstart(i),cend(i),cstart(i),cend(i),HC_d(i))       
       call zextract_dns(S,cstart(i),cend(i),cstart(i),cend(i),SC_d(i))
       
    enddo
    
    do i=1,ncont
       print*, '(int) extract central-contact',i
       i1=struct%mat_PL_start(struct%cblk(i))
       i2=struct%mat_PL_end(struct%cblk(i)) 
       j1=cstart(i); j2=j1+(ncdim(i)+pnegf%surfdim(i))/2-1
       print*, 'Interaction block:',i1,i2,j1,j2
       call zextract(H,i1,i2,j1,j2,TM(i))         
       call zextract(S,i1,i2,j1,j2,ST(i))       
    enddo

    !********************************************************************
    !****** MAIN LOOP OF GDFTB -> CHARGE or TUNNELING COMPUTATIONS ******
    !********************************************************************

    if (ncont.eq.0) then
      mumin=Efermi(1)
      mumax=mumin
    else  
      mumin=minval(Efermi(1:ncont)-mu(1:ncont))
      mumax=maxval(Efermi(1:ncont)-mu(1:ncont))
      nc_vec=maxloc(Efermi(1:ncont)-mu(1:ncont))
      nc=nc_vec(1)
    endif

    muref = mumax

#ifdef MPI 
   call MPI_BARRIER(mpi_comm,ierr)   
#endif
    !if(id0) open(47,file='partial.dat')

    !If at first SCC cycle save surf. green's else load
    !If ReadOldSGF is true then flag is forced to 0, so it reads

    if (pnegf%iteration.eq.1) then 
      pnegf%ReadOldSGF = 2
    else
      pnegf%ReadOldSGF = 0
    endif

    avncyc = 0.0   !average number of decimation iteration for SGFs computation     

    ! ***************************************************************************
    ! 1.  INTEGRATION OVER THE CIRCLE FROM PI...ALPHA  Np(1)
    ! ***************************************************************************
    ! NEW INTEGRATION FOR COMPLEX DENSITY:
    !---------------------------------------------------------
    !   2  [ /           ]     1  [ /         it   ] 
    !  --- [ | Gr(z) dz  ] =   -- [ | iGr(t)Re  dt ]  
    !  2pi [ /           ]     pi [ /              ]
    !
    ! The factor 2 is for spin, later we compute rho=i[Gr-Ga]
    !---------------------------------------------------------

    Omega = N_omega*Kb*Temp

    Centre = (Lambda**2-Elow**2+(muref-Omega)**2)/(2.d0*(muref-Omega-Elow))
    Rad = Centre - Elow

    if (Temp.ne.0.d0) then        
      alpha = datan(Lambda/(muref-Centre-Omega)) 
    else
      alpha = 0.1d0*pi             
    end if

    !Setting weights for gaussian integration 
    
    allocate(wght(max(Np(1),Np(2))),stat=err)
    IF (err.ne.0) STOP 'no space for allocation (wght)'
    allocate(pnts(max(Np(1),Np(2))),stat=err)
    IF (err.ne.0) STOP 'no space for allocation (pnts)'

    call gauleg(pi,alpha,pnts,wght,Np(1))

    !Computing complex integral (Common for T>=0)
    
    npid = int(Np(1)/numprocs)
    istart = id*npid+1
    if(id.ne.(numprocs-1)) then 
      iend = (id+1)*npid
    else
      iend = Np(1)
    end if

    IF (profile) then
      CALL SYSTEM_CLOCK(t1g,crg,cmg)
    endif

    ! -----------------------------------------------------------------------
    !  Integration loop starts here
    ! -----------------------------------------------------------------------

    do i = istart,iend

      if (verbose.gt.VBT) then 
        write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 1: point #',i,'/',iend,'  CPU=&
            &', id
      endif

      teta = pnts(i)

      Pc = Rad*exp(j*teta)
      Ec = Centre+Pc
      dt = 1.d0*wght(i)/pi 
      zt = dt*Pc*j
      if (low.lt.2) zt = zt*2.d0*j

      ! -----------------------------------------------------------------------
      !  Calculation of contact self-energies
      ! -----------------------------------------------------------------------
      ! For the time HC and SC are dense, GS is sparse (already allocated)
      ! TM and ST are sparse, SelfEneR is allocated inside SelfEnergy
      ! -----------------------------------------------------------------------
      do l=1,ncont
         pnegf%activecont=l
         call surface_green(Ec,HC_d(l),SC_d(l),pnegf,i,ncyc,GS_d(l),GS(l))
         avncyc = avncyc + ncyc
      enddo

      ! -- GF Calculation ----------------------------------------------------
      ! 
      !Tlc= Ec*ST - TM 
      !Tcl= (conjg(Ec)*ST - TM )
      !Array di GS sparse.

      !Tlc: matrici di interazione (ES-H) device-contatti (l=layer,c=contact)
      !Tcl: matrici di interazione (ES-H) contatti-device (l=layer,c=contact)

      do i1=1,ncont
        
        call prealloc_sum(TM(i1),ST(i1),(-1.d0, 0.d0),Ec,Tlc(i1))

        call prealloc_sum(TM(i1),ST(i1),(-1.d0, 0.d0),conjg(Ec),TpMt)

        call zdagacsr(TpMt,Tcl(i1))

        call destroy(TpMt)

      enddo
    
      if (id0.and.verbose.gt.VBT) call message_clock('Compute Green`s funct ')

      call SelfEnergies(Ec,ncont,GS,Tlc,Tcl,SelfEneR)

      if (mem) then 
         call calls_eq_mem(HM,SM,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,struct,low) 
      else
         !call calls_eq_dsk(HM,SM,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,struct)
      endif

      if (id0.and.verbose.gt.VBT) call write_clock
      !      
      do i1=1,ncont
        call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
      enddo

      ! ------------------------------------------------------------------------------
      if (id0.and.verbose.gt.VBT) call message_clock('Density matrix update ') 
 
      ! ------------------------------------------------------------------------------
      ! Note: the extra parts of G<, needed to correct for the charges 
      ! due to the overlap with the contacts are now included in GreenR
      ! ------------------------------------------------------------------------------
      ! GMk< =   G< Tk gk* + i Gr Tk gk - i Gr Tk gk*
      !      = ... = -2 Im{ Gr Tk gk }   
      ! ------------------------------------------------------------------------------
      if(pnegf%DorE.eq.'D'.or.pnegf%DorE.eq.'B') then
         CALL concat(DensMat,zt,GreenR,1,1)
      endif
      if(pnegf%DorE.eq.'E'.or.pnegf%DorE.eq.'B') then
         CALL concat(EnMat,zt*Ec,GreenR,1,1)
      endif

      if (id0.and.verbose.gt.VBT) call write_clock

      call destroy(GreenR)

      !if(id0) then
      !   write(*,*) '----------------------------------------------------'
      !   call writeMemInfo(6)
      !   call writePeakInfo(6)
      !   write(*,*) '----------------------------------------------------'
      !endif  
    end do
    ! *******************************************************************************
    !   END OF INTEGRATION OVER THE CIRCLE ...Np(1)
    ! *******************************************************************************

    IF (profile) then
      CALL SYSTEM_CLOCK(t2g,crg,cmg)
      WRITE(*,*) '-------- GREENDFTB: INTEGRATION OVER THE CIRCLE ----'
      CALL writeMemInfo(6)
      CALL writePeakInfo(6)
      WRITE(*,*) 'Operation time: ',(t2g-t1g)*1.d0/crg,'sec'
      WRITE(*,*) '----------------------------------------------------'
    ENDIF

    ! *******************************************************************************
    ! 2. INTEGRATION OVER THE SEGMENT [muref+Omega+j*Lambda,muref-Omega+j*Lambda]
    ! (Temp /= 0) OR OVER THE CIRCLE  ALPHA..0   (Temp == 0)          Np(2)
    ! *******************************************************************************
    ! NEW INTEGRATION FOR COMPLEX DENSITY (T>0):
    !----------------------------------------------------
    !   2  [ /           ]     1  [ /                  ] 
    !  --- [ |  Gr(z) dz ] =   -- [ | Gr(z)*(z2-z1)*dt ]  
    !  2pi [ /           ]     pi [ /                  ]
    !----------------------------------------------------

    if (Temp.eq.0.d0) then                        ! Circle integration T=0

      call  gauleg(alpha,0.d0,pnts,wght,Np(2))

    else                                          ! Segment integration T>0

      z1 = muref + Omega + j*Lambda
      z2 = muref - Omega + j*Lambda

      z_diff = z2 - z1

      call  gauleg(1.d0,0.d0,pnts,wght,Np(2))    !Setting weights for integration

    endif

    !Computing complex integral

    npid = int(Np(2)/numprocs)
    istart = id*npid+1
    if(id.ne.(numprocs-1)) then 
      iend = (id+1)*npid
    else
      iend = Np(2)
    end if

    IF (profile) then
      CALL SYSTEM_CLOCK(t1g,crg,cmg)
    ENDIF

    do i = istart,iend

      if (verbose.gt.VBT) then
        write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 2: point #',i,'/',iend,'  CPU=&
            &', id
      endif

      teta = pnts(i)

      dt = 1.d0*wght(i)/pi

      if (Temp.eq.0.d0) then                      ! Circle integration T=0            
         Pc = Rad*exp(j*teta)
         Ec = Centre+Pc
         zt = dt*Pc*j
      else                                         ! Segment integration T>0
         Ec = z1 + teta*z_diff
         zt = z_diff*fermi_fc(Ec,muref,Kb*Temp)*dt
      endif

      if(low.lt.2) zt = zt*2.d0*j
      

      do l=1,ncont
         pnegf%activecont=l
         npc = Np(1) + i
         call surface_green(Ec,HC_d(l),SC_d(l),pnegf,npc,ncyc,GS_d(l),GS(l))
         avncyc = avncyc + ncyc
      enddo

      ! -- GF Calculation ----------------------------------------------------
      ! 
      !Tlc= Ec*ST - TM 
      !Tcl= (conjg(Ec)*ST - TM )
      !Array di GS sparse.

      !Tlc: matrici di interazione (ES-H) device-contatti (l=layer,c=contact) 
      !Tcl: matrici di interazione (ES-H) contatti-device (l=layer,c=contact) 

      do i1=1,ncont

        call prealloc_sum(TM(i1),ST(i1),(-1.d0, 0.d0),Ec,Tlc(i1))

        call prealloc_sum(TM(i1),ST(i1),(-1.d0, 0.d0),conjg(Ec),TpMT)

        call zdagacsr(TpMt,Tcl(i1))

        call destroy(TpMt)

      enddo

      call SelfEnergies(Ec,ncont,GS,Tlc,Tcl,SelfEneR)

      if (id0.and.verbose.gt.VBT) call message_clock('Compute Green`s funct ') 

      if (mem) then
         call calls_eq_mem(HM,SM,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,struct,low)
      else
         !call calls_eq_dsk(HM,SM,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,struct)
      endif

      if (id0.and.verbose.gt.VBT) call write_clock

      do i1=1,ncont
        call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
      enddo

      if (id0.and.verbose.gt.VBT) call message_clock('Density matrix update ') 

      ! ------------------------------------------------------------------------------
      ! Integration of G< on the complex plane. 
      ! The real part comes from Im{Int{Gr(E)dE}}=Im{Int{Gr(z) R exp(iO) i dO}}=
      !                                          =Im{ Sum{ iGr(zi) Pc Wi } }   =
      !                                          =Sum{ Re{Gr(zi) Pc Wi} }  
      ! ------------------------------------------------------------------------------
      ! Note: the extra parts of G<, needed to correct for the charges 
      ! due to the overlap with the contacts are now included in GreenR
      ! ------------------------------------------------------------------------------
      ! GMk< =   G< Tk gk* + i Gr Tk gk - i Gr Tk gk*
      !      = ... = -2 Im{ Gr Tk gk }   
      ! ------------------------------------------------------------------------------  
      if(pnegf%DorE.eq.'D'.or.pnegf%DorE.eq.'B') then
         CALL concat(DensMat,zt,GreenR,1,1)
      endif
      if(pnegf%DorE.eq.'E'.or.pnegf%DorE.eq.'B') then
         CALL concat(EnMat,zt*Ec,GreenR,1,1)
      endif
 
      call destroy(GreenR)  

      if (id0.and.verbose.gt.VBT) call write_clock      


    end do

    IF (profile) then
      CALL SYSTEM_CLOCK(t2g,crg,cmg)
      WRITE(*,*) '----- GREENDFTB: INTEGRATION OVER THE SEGMENT (2) --, CPU #',id
      CALL writeMemInfo(6)
      CALL writePeakInfo(6)
      WRITE(*,*) 'Operation time: ',(t2g-t1g)*1.d0/crg,'sec'
      WRITE(*,*) '----------------------------------------------------'
    END IF


    ! *******************************************************************************
    ! 3. SUMMATION OVER THE POLES ENCLOSED IN THE CONTOUR  (NumPoles)
    ! *******************************************************************************          
    ! NEW INTEGRATION FOR COMPLEX DENSITY (T>=0):
    !---------------------------------------------------------------------
    !             [ 2                  ]    1      
    !  2 pi j* Res[ -- *Gr(z_k)f(z_k)  ] =  -- *2*pi*j*(-kb T)* Gr(z_k)     
    !             [ 2pi                ]    pi     
    !                                                  (-kb*T) <- Residue
    !---------------------------------------------------------------------
    npid = int(NumPoles/numprocs)
    istart = id*npid+1
    if(id.ne.(numprocs-1)) then 
      iend = (id+1)*npid
    else
      iend = NumPoles
    end if

    do i = istart,iend

      if (verbose.gt.VBT) then
        write(6,'(a17,i3,a1,i3,a6,i3)') 'POLES: point #',i,'/',iend,'  CPU=', id
      endif

      Ec = muref + j*Kb*Temp*pi* (2.d0*i - 1.d0)       
      
      do l=1,ncont
         pnegf%activecont=l
         npc =  Np(1)+Np(2)+i
         call surface_green(Ec,HC_d(l),SC_d(l),pnegf,npc,ncyc,GS_d(l),GS(l))
         avncyc = avncyc + ncyc
      enddo
      
      ! -- GF Calculation ----------------------------------------------------
      do i1=1,ncont

        call prealloc_sum(TM(i1),ST(i1),(-1.d0, 0.d0),Ec,Tlc(i1))

        call prealloc_sum(TM(i1),ST(i1),(-1.d0, 0.d0),conjg(Ec),TpMt)

        call zdagacsr(TpMt,Tcl(i1))

        call destroy(TpMt)

      enddo

      if (id0.and.verbose.gt.VBT) call message_clock('Compute Green`s funct ')

      call SelfEnergies(Ec,ncont,GS,Tlc,Tcl,SelfEneR)

      if (mem) then
         call calls_eq_mem(HM,SM,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,struct,low)
      else
         !call calls_eq_dsk(HM,SM,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,struct)
      endif

      if (id0.and.verbose.gt.VBT) call write_clock

      do i1=1,ncont
        call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
      enddo

      if (id0.and.verbose.gt.VBT) call message_clock('Density matrix update ') 

      zt= -2.d0*j*Kb*Temp*(1.d0,0.d0)
      if(low.lt.2) zt = zt*2.d0*j

      if(pnegf%DorE.eq.'D'.or.pnegf%DorE.eq.'B') then
         CALL concat(DensMat,zt,GreenR,1,1)
      endif
      if(pnegf%DorE.eq.'E'.or.pnegf%DorE.eq.'B') then
         CALL concat(EnMat,zt*Ec,GreenR,1,1)
      endif

      CALL destroy(GreenR) 

      if (id0.and.verbose.gt.VBT) call write_clock

    end do

    ! *******************************************************************************
    ! END OF SUMMATION OVER THE POLES ENCLOSED IN THE CONTOUR
    ! *******************************************************************************                  

    !if(id0) then
    !   write(*,*) '----------------------------------------------------'
    !   call writeMemInfo(6)
    !   call writePeakInfo(6)
    !   write(*,*) '----------------------------------------------------'
    !endif 
    ! -------------------------------------------------------------------------   
    ! Build Specral density, i(Gr-Ga), for the equilibrium part
    ! -------------------------------------------------------------------------
    
    if(pnegf%DorE.eq.'D'.or.pnegf%DorE.eq.'B') then
      call zspectral(DensMat,DensMat,0,TpMt)
      ! To do: check outer blocks of DensMat !
      call destroy(DensMat)
      call clone(TpMt,DensMat)
      call destroy(TpMt)
    end if
    if(pnegf%DorE.eq.'E'.or.pnegf%DorE.eq.'B') then
      call zspectral(EnMat,EnMat,0,TpMt)
      ! To do: check outer blocks of DensMat !
      call destroy(EnMat)
      call clone(TpMt,EnMat)
      call destroy(TpMt)
    end if



    ! *******************************************************************************
    ! 4. INTEGRATION OVER THE REAL SEGMENT mumin-Omega...mumax+Omega   Np(3)+2*npT
    ! *******************************************************************************
    ! NEW INTEGRATION FOR COMPLEX DENSITY (T>=0):
    !----------------------------------------------------------
    !    2    [ /           ]     1  /                  
    !  ------ [ | G<(E) dE  ] =   -- | Gr(E)*GAM_c(E)*Ga(E)*dE   
    !  2 pi i [ /           ]     pi /                  
    !----------------------------------------------------------
    
    if (mumax.gt.mumin) then 
      ! Real segment integration
      deallocate(pnts,wght)
       
      ! compute extended number of points due to kT.
      npT=nint(Np(3)/(mumax-mumin))*Omega
      allocate(pnts(Np(3)+2*npT))
      allocate(wght(Np(3)+2*npT))

      !Setting weights for gaussian integration
      call gauleg(mumin-omega,mumax+omega,pnts,wght,Np(3)+2*npT)

      !Computing real axis integral

      npid = int((Np(3)+2*npT)/numprocs)
      istart = id*npid+1
      if(id.ne.(numprocs-1)) then 
        iend = (id+1)*npid
      else
        iend = Np(3)+2*npT
      end if

      IF (profile) THEN
        CALL SYSTEM_CLOCK(t1g,crg,cmg)
      ENDIF

      do i = istart,iend

        if (verbose.gt.VBT) then
          write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 3: point #',i,'/',iend,'  CP&
              &U=', id
        endif

        Ec = cmplx(pnts(i),0.0,dp) 
        dt = wght(i)/pi
        zt = dt*(1.d0,0.d0)

        do j1=1,ncont
          frm_f(j1)=fermi_f(dreal(Ec),Efermi(j1)-mu(j1),Kb*Temp)
        enddo

        ! Real segment integration
        do l=1,ncont
           pnegf%activecont=l
           npc = Np(1)+Np(2)+NumPoles+i
           call surface_green(Ec+j*pnegf%delta,HC_d(l),SC_d(l),pnegf,npc,ncyc,&
                                                                  GS_d(l),GS(l))            
           avncyc = avncyc + ncyc
        enddo
        
        ! -- GF Calculation ----------------------------------------------------
        ! 
        !Tlc= Ec*ST - TM 
        !Tcl= (Ec*ST - TM )^+    Ec is real !

        do i1=1,ncont
          call prealloc_sum(TM(i1),ST(i1),(-1.d0, 0.d0),Ec,Tlc(i1))

          call zdagacsr(Tlc(i1),Tcl(i1))
       enddo

        if (id0.and.verbose.gt.VBT) call message_clock('Compute Green`s funct ')

        call SelfEnergies(Ec,ncont,GS,Tlc,Tcl,SelfEneR)

        if (mem) then
           call calls_neq_mem(HM,SM,Ec,SelfEneR,Tlc,Tcl,GS,struct,frm_f,nc,GreenR,low)
        else 
           !call calls_neq_dsk(HM,SM,Ec,SelfEneR,Tlc,Tcl,GS,struct,frm_f,nc,GreenR)
        endif

        if (id0.and.verbose.gt.VBT) call write_clock      

        do i1=1,ncont
           call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
        enddo

        if (id0.and.verbose.gt.VBT) call message_clock('Density matrix update ') 

        if(pnegf%DorE.eq.'D'.or.pnegf%DorE.eq.'B') then
           CALL concat(DensMat,zt,GreenR,1,1)
        endif
        if(pnegf%DorE.eq.'E'.or.pnegf%DorE.eq.'B') then
           CALL concat(DensMat,zt*Ec,GreenR,1,1)
        endif

        call destroy(GreenR)         

        if (id0.and.verbose.gt.VBT) call write_clock

        !Charge computation

      enddo

      IF (profile) then
        CALL SYSTEM_CLOCK(t2g,crg,cmg)
        WRITE(*,*) '------ GREENDFTB: INTEGRATION OVER REAL SEGMENT ----'
        CALL writeMemInfo(6)
        CALL writePeakInfo(6)
        WRITE(*,*) 'Operation time: ',(t2g-t1g)*1.d0/crg,'sec'
        WRITE(*,*) '----------------------------------------------------'
      ENDIF

    else   ! if mumax <= mumin

      Np(3) = 0

    end if ! if mumax <= mumin  


    !close(47)
    deallocate(pnts,wght) 

    call destroy(HM,SM)

    do i=1,ncont
       call destroy(HC_d(i),SC_d(i))
       call destroy(TM(i),ST(i))
    enddo


#ifdef MPI
    ! *******************************************************************************
    !  MPI REDUCE ALL  (REDUCE D.M. ON ALL NODES) 
    ! *******************************************************************************

    if ((id0).and.verbose.gt.VBT) call message_clock('MPI gather masking ') 
    !sorting DensMat for MPI 
    if(param%DorE.eq.'D'.or.param%DorE.eq.'B') then
       CALL  msort(DensMat)
    endif
    if(param%DorE.eq.'E'.or.param%DorE.eq.'B') then
       CALL  msort(EnMat)
    endif

    call MPI_BARRIER( mpi_comm, ierr)

    if(param%DorE.eq.'D'.or.param%DorE.eq.'B') then
       ibsize=DensMat%nnz
       call log_allocate(Mat_st,ibsize)
       call log_allocate(Mat_rcv,ibsize)

       Mat_st(1:ibsize)=DensMat%nzval(1:ibsize)
       Mat_rcv(1:ibsize)=0.d0
       
       call MPI_ALLREDUCE(Mat_st,Mat_rcv,ibsize,MPI_DOUBLE_PRECISION,&
                                               MPI_SUM,mpi_comm,ierr)

       DensMat%nzval(1:ibsize)=Mat_rcv(1:ibsize)

       call log_deallocate(Mat_st)
       call log_deallocate(Mat_rcv)
    endif

    if(param%DorE.eq.'E'.or.param%DorE.eq.'B') then
       ibsize=EnMat%nnz
       call log_allocate(Mat_st,ibsize)
       call log_allocate(Mat_rcv,ibsize)

       Mat_st(1:ibsize)=EnMat%nzval(1:ibsize)
       Mat_rcv(1:ibsize)=0.d0
       
       call MPI_ALLREDUCE(Mat_st,Mat_rcv,ibsize,MPI_DOUBLE_PRECISION,&
                                               MPI_SUM,mpi_comm,ierr)

       EnMat%nzval(1:ibsize)=Mat_rcv(1:ibsize)

       call log_deallocate(Mat_st)
       call log_deallocate(Mat_rcv)
    endif

    if (id0.and.verbose.gt.VBT) call write_clock
#endif

    if (id0.and.verbose.gt.VBT) then
      write(*,'(73("="))')
      write(*,'(A,f8.3)') 'Average number of decimation iter.:', &
           avncyc*1.0*numprocs/(Np(1)+Np(2)+Np(3)+2*npT+NumPoles)
    endif

    if(id0.and.verbose.gt.70) then
      call writeMemInfo(6)
      call writePeakInfo(6)
      write(*,*)
    endif

    return 

  end subroutine Contour_int

!***********************************************************************************
!
! Subroutine di mascheramento (zeri inclusi) necessario all'MPI
!
!***********************************************************************************

SUBROUTINE mask_dens(D,S)

  implicit none

  TYPE(z_CSR) :: S
  TYPE(r_CSR) :: D,DS
  INTEGER :: k,i,m

  call msort(D)
  call msort(S)

  CALL create(DS,S%nrow,S%ncol,S%nnz)
  DS%rowpnt(:)=S%rowpnt(:)
  DS%colind=S%colind(:)

  k=0

  DO k=1,S%nrow

     m=D%rowpnt(k)

     DO i=S%rowpnt(k),S%rowpnt(k+1)
     
        DO WHILE(m.LT.D%rowpnt(k+1))
        !DO j=D%rowpnt(k),D%rowpnt(k+1)-1
           IF (S%colind(i).EQ.D%colind(m)) THEN
              DS%nzval(i)=S%nzval(i)
              m=m+1
              exit
           ELSE IF (S%colind(i).GT.D%colind(m)) THEN
              DS%nzval(i)=0.d0
              exit 
           ELSE 
              m=m+1
           ENDIF
        ENDDO

     ENDDO

  ENDDO

  call destroy(D)
  CALL create(D,DS%nrow,DS%ncol,DS%nnz)
  D%rowpnt=DS%rowpnt
  D%colind=DS%colind
  D%nzval=DS%nzval
  call destroy(DS)

END SUBROUTINE mask_dens

!***********************************************************************************
!
! Subroutine di mascheramento (zeri inclusi) necessario all'MPI (di confronto,
! versione meno performante)
!
!***********************************************************************************

subroutine mask_dens2(D,S)

implicit none 

  TYPE(z_CSR) :: S
  TYPE(r_CSR) :: D,DS,SR
  COMPLEX(kind=dp), PARAMETER :: drop=(1.d-20,0.d0)
  INTEGER :: i

  !CALL create(SR,S%nrow,S%ncol,S%nnz)
  !SR%rowpnt(:)=S%rowpnt(:)
  !SR%colind(:)=S%colind(:)
  !SR%nzval(:)=dreal(S%nzval(:))
  
  !call resumcsrs(D,SR,drop,DS)
  !call destroy(SR)

  call concat(D,drop,S,'R',1,1)
  call clone(D,DS)

  call destroy(D)
  call create(D,S%nrow,S%ncol,S%nnz)

  call mask(DS,S,D)

  call destroy(DS)

end subroutine mask_dens2


!******************************************************************
!
! Subroutine for Gaussian integration. 
! Find zeros and wights of Gauss-Legendre polinomials
!
!******************************************************************
subroutine gauleg(x1,x2,x,w,n)

  real(kind=dp), PARAMETER :: ACC = 1d-15
  
  INTEGER n
  real(kind=dp) :: x1,x2,x(n),w(n)
  
  INTEGER i,k,m
  real(kind=dp) :: p1,p2,p3,pp,xl,xm,z,z1
  
  m=(n+1)/2
  
  xm=0.5d0*(x2+x1)
  xl=0.5d0*(x2-x1)
  
  do i=1,m
    
    z=cos(Pi*(i-0.25d0)/(n+0.5d0))
    
    do
      p1=1.d0
      p2=0.d0
      
      ! Legendre polynomial p1 evaluated by rec. relations:
      do k=1,n
        p3=p2
        p2=p1
        p1=((2.d0*k-1.d0)*z*p2-(k-1.d0)*p3)/k
      enddo
      ! Derivative pp using the relation of p1 and p2:
      pp=n*(z*p1-p2)/(z*z-1.d0)
      
      ! Newton method to refine the zeros:
      z1=z
      z=z1-p1/pp
      
      if(abs(z-z1).le.ACC) exit
    enddo
    
    ! Scale the interval to x1..x2:
    x(i)=xm-xl*z
    x(n+1-i)=xm+xl*z
    w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
    w(n+1-i)=w(i)
  enddo

  return
  
end subroutine Gauleg
  

end  module negf_greendftb
