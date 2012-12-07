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


module integrations 

 use ln_precision
 use ln_constants
 use ln_allocation
 use lib_param
 use mpi_globals, only : id, numprocs, id0
 use input_output
 use ln_structure
 use fermi_dist
 use sparsekit_drv
 use inversions
 use iterative
 use iterative_dns
 use mat_def
 use ln_extract
 use contselfenergy
 use clock
 use elph
 use energy_mesh 
  
 implicit none
 private

 public :: contour_int     ! standard contour integrations for DFT(B) 
 public :: real_axis_int   ! real-axis integration for DFT(B)
 public :: contour_int_n   ! contour integration for CB
 public :: real_axis_int_n ! real axis integration for CB
 public :: contour_int_p   ! contour integration for VB 
 public :: tunneling_and_current
 public :: integrate       ! integration of tunneling
 public :: compute_dos                ! compute local dos only

 ! ////////////////////////////////////////////////////////////
 ! Under development:
 public :: init_emesh, destroy_emesh
 private :: adaptive_int, trapez23 
 !public :: contour_int_ph, real_axis_int_ph, real_axis_int_ph2
 !
 type TG_pointer
   type(z_CSR),  pointer :: pG => null()   
   integer :: ind   
 end type TG_pointer
 ! ////////////////////////////////////////////////////////////

 integer, PARAMETER :: VBT=70


contains
  

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine compute_dos(negf) 
    type(Tnegf) :: negf

    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    Type(z_CSR) ::  Gr

    integer :: N, i, i1, it
    integer :: outer, nbl, ncont

    real(dp) :: ncyc
    complex(dp) :: Ec

    outer = 1

    it = negf%iteration
    nbl = negf%str%num_PLs
    ncont = negf%str%num_conts

    N = nint((negf%Emax-negf%Emin)/negf%Estep)

    !print*,'N=',N
    !print*,'delta=',negf%delta

    open(101,file='dos.dat')

    do i=1,N
  
       Ec=(negf%Emin+i*negf%Estep)+negf%delta*(0.d0,1.d0)

       call compute_contacts(Ec,negf,i,ncyc,Tlc,Tcl,SelfEneR,GS)
   
       call calls_eq_mem_dns(negf,Ec,SelfEneR,Tlc,Tcl,GS,Gr,negf%str,outer)

       do i1=1,ncont
          call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
       enddo

       negf%dos = -aimag( trace(Gr) )/pi

       write(101,*) real(Ec), negf%dos
    !   print*, real(Ec), negf%dos, Gr%nnz

       call destroy(Gr)

    enddo

    call writememinfo(6)

    close(101)   

  end subroutine compute_dos

!-----------------------------------------------------------------------
! Contour integration for density matrix 
! DOES INCLUDE FACTOR 2 FOR SPIN !! 
!-----------------------------------------------------------------------
  subroutine contour_int_n(negf)

    type(Tnegf) :: negf
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    type(z_CSR) :: GreenR, TmpMt 

    integer :: npid, istart, iend, NumPoles
    integer :: i, i1, outer, it, ncont, nbl

    real(dp), DIMENSION(:), ALLOCATABLE :: wght,pnts   ! Gauss-quadrature points
    real(dp) :: Omega, Lambda
    real(dp) :: muref, kbT
    real(dp) :: ncyc

    complex(dp) :: z1,z2,z_diff, zt
    complex(dp) :: Ec, ff

    it = negf%iteration
    ncont = negf%str%num_conts
    nbl = negf%str%num_PLs
    kbT = negf%kbT
    muref = negf%muref 

    Omega = negf%n_kt * kbT
    Lambda = 2.d0* negf%n_poles * KbT * pi
    NumPoles = negf%n_poles

    outer = negf%outer

    call create(TmpMt,negf%H%nrow,negf%H%ncol,negf%H%nrow)
    call initialize(TmpMt)


    ! *******************************************************************************
    ! 1. INTEGRATION OVER THE SEGMENT [Ec - dEc , Ec - dEc + j*Lambda]
    ! *******************************************************************************
    ! NEW INTEGRATION FOR COMPLEX DENSITY (T>0):
    !----------------------------------------------------
    !   g  [ /                ]      g   [ /                       ] 
    !  --- [ |  Gr(z)*f(z) dz ] =   ---  [ | Gr(z)*(z2-z1)*f(z)*dt ]  
    !  2pi [ /                ]     2*pi [ /                       ]
    !----------------------------------------------------

    z1 = negf%Ec + negf%DeltaEc
    z2 = negf%Ec + negf%DeltaEc + j*Lambda

    z_diff = z2 - z1

    allocate(wght(negf%Np_n(1)))
    allocate(pnts(negf%Np_n(1)))

    call gauleg(0.d0,1.d0,pnts,wght,negf%Np_n(1))

    npid = int(negf%Np_n(1)/numprocs)
    istart = id*npid+1
    if(id.ne.(numprocs-1)) then 
       iend = (id+1)*npid
    else
       iend = negf%Np_n(1)
    end if

    do i = istart,iend

       if (negf%verbose.gt.VBT) then
          write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 1: point #',i,'/',iend,'  CPU=&
               &', id
       endif

       if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Green`s funct ')

       Ec = z1 + pnts(i) * z_diff

       ff = fermi_fc(Ec,muref,KbT)

       zt = negf%spin * z_diff * ff * wght(i) / (2.d0 *pi)

       call compute_contacts(Ec,negf,i,ncyc,Tlc,Tcl,SelfEneR,GS)

       call calls_eq_mem_dns(negf,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer)

       if(negf%DorE.eq.'D') then
          call concat(TmpMt,zt,GreenR,1,1)
       endif
       if(negf%DorE.eq.'E') then
          call concat(TmpMt,zt*Ec,GreenR,1,1)
       endif

       call destroy(GreenR) 

       do i1=1,ncont
          call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
       enddo

       if (id0.and.negf%verbose.gt.VBT) call write_clock

    enddo

    deallocate(wght)
    deallocate(pnts)


    ! *******************************************************************************
    ! 2. INTEGRATION OVER THE SEGMENT [Ec + j*Lambda, mu(r) + Omega+j*Lambda]
    ! *******************************************************************************
    ! NEW INTEGRATION FOR COMPLEX DENSITY (T>0):
    !----------------------------------------------------
    !   g  [ /                ]      g   [ /                       ] 
    !  --- [ |  Gr(z)*f(z) dz ] =   ---  [ | Gr(z)*(z2-z1)*f(z)*dt ]  
    !  2pi [ /                ]     2*pi [ /                       ]
    !----------------------------------------------------


    allocate(wght(negf%Np_n(2)))
    allocate(pnts(negf%Np_n(2)))

    call gauleg(0.d0,1.d0,pnts,wght,negf%Np_n(2))    !Setting weights for integration

    z1 = negf%Ec + negf%DeltaEc  + j*Lambda
    z2 = muref + Omega + j*Lambda

    z_diff = z2 - z1

    npid = int(negf%Np_n(2)/numprocs)

    istart = id*npid+1
    if(id.ne.(numprocs-1)) then 
       iend = (id+1)*npid
    else
       iend = negf%Np_n(2)
    end if
  
    do i = istart,iend

       if (negf%verbose.gt.VBT) then
          write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 2: point #',i,'/',iend,'  CPU=&
               &', id
       endif

       if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Green`s funct ')

       Ec = z1 + pnts(i) * z_diff

       ff = fermi_fc(Ec,muref,KbT)

       zt = negf%spin *  z_diff * ff * wght(i) / (2.d0 *pi)

       call compute_contacts(Ec,negf,negf%Np_n(1)+i,ncyc,Tlc,Tcl,SelfEneR,GS)

       call calls_eq_mem_dns(negf,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer) 

       if(negf%DorE.eq.'D') then
          call concat(TmpMt,zt,GreenR,1,1)
       endif
       if(negf%DorE.eq.'E') then
          call concat(TmpMt,zt*Ec,GreenR,1,1)
       endif

       call destroy(GreenR)

       do i1=1,ncont
          call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
       enddo

       if (id0.and.negf%verbose.gt.VBT) call write_clock

    enddo

    deallocate(wght)
    deallocate(pnts)



    ! *******************************************************************************
    ! 3. SUMMATION OVER THE POLES ENCLOSED IN THE CONTOUR  (NumPoles)
    ! *******************************************************************************          
    ! NEW INTEGRATION FOR COMPLEX DENSITY (T>=0):
    !---------------------------------------------------------------------
    !                  [ 1                  ]          
    ! 2 pi * g * j* Res[ -- *Gr(z_k)f(z_k)  ] =  j*(-kb T)* Gr(z_k)     
    !                  [ 2pi                ]         
    !                                              (-kb*T) <- Residue
    !---------------------------------------------------------------------

    npid = int(NumPoles/numprocs)
    istart = id*npid+1
    if(id.ne.(numprocs-1)) then 
       iend = (id+1)*npid
    else
       iend = NumPoles
    end if

    do i = istart,iend

       if (negf%verbose.gt.VBT) then
          write(6,'(a17,i3,a1,i3,a6,i3)') 'POLES: point #',i,'/',iend,'  CPU=', id
       endif

       if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Green`s funct ')

       Ec = muref + j * KbT *pi* (2.d0*real(i,dp) - 1.d0)   

       zt= -j * KbT * negf%spin *(1.d0,0.d0) 

       call compute_contacts(Ec,negf,negf%Np_n(1)+negf%Np_n(2)+i,ncyc,Tlc,Tcl,SelfEneR,GS)

       call calls_eq_mem_dns(negf,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer)

       if(negf%DorE.eq.'D') then
          call concat(TmpMt,zt,GreenR,1,1)
       endif
       if(negf%DorE.eq.'E') then
          call concat(TmpMt,zt*Ec,GreenR,1,1)
       endif

       call destroy(GreenR)   

       do i1=1,ncont
          call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
       enddo

       if (id0.and.negf%verbose.gt.VBT) call write_clock

    enddo

    if(negf%DorE.eq.'D') then
       call zspectral(TmpMt,TmpMt,0,negf%rho)
    endif
    if(negf%DorE.eq.'E') then
       call zspectral(TmpMt,TmpMt,0,negf%rho_eps)
    endif


    call destroy(TmpMt)

  end subroutine contour_int_n
!--------------------------------------------------------------------------------
!-----------------------------------------------------------------------
! Contour integration for density matrix 
! DOES INCLUDE FACTOR 2 FOR SPIN !! 
!-----------------------------------------------------------------------
  subroutine contour_int_p(negf)

    type(Tnegf) :: negf
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    type(z_CSR) :: GreenR, TmpMt 

    integer :: npid, istart, iend, NumPoles
    integer :: i, i1, outer, it, ncont, nbl

    real(dp), DIMENSION(:), ALLOCATABLE :: wght,pnts   ! Gauss-quadrature points
    real(dp) :: Omega, Lambda
    real(dp) :: muref, kbT
    real(dp) :: ncyc

    complex(dp) :: z1,z2,z_diff, zt
    complex(dp) :: Ev, ff

    it = negf%iteration
    ncont = negf%str%num_conts
    nbl = negf%str%num_PLs
    kbT = negf%kbT
    muref = negf%muref

    Omega = negf%n_kt * kbT
    Lambda = 2.d0* negf%n_poles * KbT * pi
    NumPoles = negf%n_poles
    
    outer = negf%outer

    call create(TmpMt,negf%H%nrow,negf%H%ncol,negf%H%nrow)
    call initialize(TmpMt)

    !1. INTEGRATION OVER THE SEGMENT

    z1 = negf%Ev + negf%DeltaEv + j*Lambda
    z2 = negf%Ev + negf%DeltaEv

    z_diff = z2 - z1

    allocate(wght(negf%Np_p(1)))
    allocate(pnts(negf%Np_p(1)))

    call gauleg(0.d0,1.d0,pnts,wght,negf%Np_p(1))

    npid = int(negf%Np_p(1)/numprocs)
    istart = id*npid+1
    if(id.ne.(numprocs-1)) then 
       iend = (id+1)*npid
    else
       iend = negf%Np_p(1)
    end if


    do i = istart,iend

       !if (negf%verbose.gt.VBT) then
       !   write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 1: point #',i,'/',iend,'  CPU=&
       !        &', id
       !endif

       Ev = z1 + pnts(i) * z_diff

       ff = (1.d0,0.d0) - fermi_fc(Ev,muref,KbT)

       zt = negf%spin * z_diff * ff * wght(i) / (2.d0 *pi)

       call compute_contacts(Ev,negf,i,ncyc,Tlc,Tcl,SelfEneR,GS)

       call calls_eq_mem_dns(negf,Ev,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer)

       if(negf%DorE.eq.'D') then
          call concat(TmpMt,zt,GreenR,1,1)
       endif
       if(negf%DorE.eq.'E') then
          call concat(TmpMt,zt*Ev,GreenR,1,1)
       endif

       call destroy(GreenR) 

       do i1=1,ncont
          call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
       enddo

    enddo

    deallocate(wght)
    deallocate(pnts)

    ! 2. INTEGRATION OVER THE SEGMENT 

    allocate(wght(negf%Np_p(2)))
    allocate(pnts(negf%Np_p(2)))

    call gauleg(0.d0,1.d0,pnts,wght,negf%Np_p(2))    !Setting weights for integration

    z1 = muref - Omega + j*Lambda
    z2 = negf%Ev + negf%DeltaEv  + j*Lambda

    z_diff = z2 - z1

    npid = int(negf%Np_p(2)/numprocs)

    istart = id*npid+1
    if(id.ne.(numprocs-1)) then 
       iend = (id+1)*npid
    else
       iend = negf%Np_p(2)
    end if
  
    do i = istart,iend

       !if (negf%verbose.gt.VBT) then
       !   write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 2: point #',i,'/',iend,'  CPU=&
       !        &', id
       !endif

       Ev = z1 + pnts(i) * z_diff

       ff = (1.d0,0.d0) - fermi_fc(Ev,muref,KbT)

       zt = z_diff * negf%spin * ff * wght(i) / (2.d0 *pi)

       call compute_contacts(Ev,negf,negf%Np_p(1)+i,ncyc,Tlc,Tcl,SelfEneR,GS)

       call calls_eq_mem_dns(negf,Ev,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer) 

       if(negf%DorE.eq.'D') then
          call concat(TmpMt,zt,GreenR,1,1)
       endif
       if(negf%DorE.eq.'E') then
          call concat(TmpMt,zt*Ev,GreenR,1,1)
       endif

       call destroy(GreenR)

       do i1=1,ncont
          call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
       enddo

    enddo

    deallocate(wght)
    deallocate(pnts)
    
    ! SUMMATION OVER THE POLES ENCLOSED IN THE CONTOUR
    
    npid = int(NumPoles/numprocs)
    istart = id*npid+1
    if(id.ne.(numprocs-1)) then 
       iend = (id+1)*npid
    else
       iend = NumPoles
    end if

    do i = istart,iend

       !if (negf%verbose.gt.VBT) then
       !   write(6,'(a17,i3,a1,i3,a6,i3)') 'POLES: point #',i,'/',iend,'  CPU=', id
       !endif

       Ev =  muref + j * KbT *pi* (2.d0*real(i,dp) - 1.d0)   

       zt= j*negf%spin*KbT*(1.d0,0.d0) 

       call compute_contacts(Ev,negf,negf%Np_p(1)+negf%Np_p(2)+i,ncyc,Tlc,Tcl,SelfEneR,GS)

       call calls_eq_mem_dns(negf,Ev,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer)

       if(negf%DorE.eq.'D') then
          call concat(TmpMt,zt,GreenR,1,1)
       endif
       if(negf%DorE.eq.'E') then
          call concat(TmpMt,zt*Ev,GreenR,1,1)
       endif

       call destroy(GreenR)  

        do i1=1,ncont
          call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
       enddo

    enddo

    if(negf%DorE.eq.'D') then
       call zspectral(TmpMt,TmpMt,0,negf%rho)
    endif
    if(negf%DorE.eq.'E') then
       call zspectral(TmpMt,TmpMt,0,negf%rho_eps)
    endif

    call destroy(TmpMt)

  end subroutine contour_int_p

  !-----------------------------------------------------------------------
  ! Contour integration for density matrix 
  ! DOES INCLUDE FACTOR 2 FOR SPIN !! 
  !-----------------------------------------------------------------------
  ! Performs the complex contur integration
  !
  ! T=0                             T>0 
  !         * * * *                          * * *      +
  !      *          *                     *         *   + 
  !    *              *                 *             **+******
  !   *                *               *                + 
  !  -*----------------*------------- -*-----------------------------
  !  Elow              muref          Elow              muref       
  !
  !  muref is the energy of the reference contact
  !
  !  
  subroutine contour_int(negf)
  
     type(Tnegf) :: negf
     Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
     type(z_CSR) :: GreenR, TmpMt 
  
     integer :: npid, istart, iend, NumPoles
     integer :: i, i1, outer, ncont, nbl
  
     real(dp), DIMENSION(:), ALLOCATABLE :: wght,pnts   ! Gauss-quadrature points
     real(dp) :: Omega, Lambda, Rad, Centre
     real(dp) :: muref, mumin, kbT, nkT, alpha
     real(dp) :: ncyc, dt, Elow
    
     complex(dp) :: z1,z2,z_diff, zt
     complex(dp) :: Ec, ff, Pc
    
     ncont = negf%str%num_conts
     nbl = negf%str%num_PLs
     kbT = negf%kbT
    
     muref = negf%mu_n
  
     nkT = negf%n_kt * kbT
     Lambda = 2.d0* negf%n_poles * KbT * pi
     NumPoles = negf%n_poles
     mumin = muref - nkT
    
     outer = negf%outer 
    
     call create(TmpMt,negf%H%nrow,negf%H%ncol,negf%H%nrow)
     call initialize(TmpMt)
     ! -----------------------------------------------------------------------
     !  Integration loop starts here
     ! -----------------------------------------------------------------------
     ! ***********************************************************************
     ! 1. INTEGRATION OVER THE CIRCLE Pi..alpha    Np(1)
     ! ***********************************************************************
     ! NEW INTEGRATION FOR COMPLEX DENSITY:
     !----------------------------------------------------
     !   g  [ /           ]     g  [ /         it   ] 
     !  --- [ | Gr(z) dz  ] =  --- [ | iGr(t)Re  dt ]  
     !  2pi [ /           ]    2pi [ /              ]
     !----------------------------------------------------
    
     Elow = negf%Ec
     Centre = (Lambda**2-Elow**2+(mumin)**2)/(2.d0*(mumin-Elow))
     Rad = Centre - Elow
    
     if (negf%kbT.ne.0.d0) then        
        alpha = atan(Lambda/(mumin-Centre)) 
     else
        alpha = 0.1d0*pi             
     end if
    
     !Setting weights for gaussian integration 
     allocate(wght(negf%Np_n(1)))
     allocate(pnts(negf%Np_n(1)))
    
     call gauleg(pi,alpha,pnts,wght,negf%Np_n(1))
    
     !Computing complex integral (Common for T>=0)
    
     npid = int(negf%Np_n(1)/numprocs)
     istart = id*npid+1
     if(id.ne.(numprocs-1)) then 
        iend = (id+1)*npid
     else
        iend = negf%Np_n(1)
     end if
    
     do i = istart,iend
    
        if (negf%verbose.gt.VBT) then
           write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 1: point #',i,'/',iend,'  CPU=&
                &', id
        endif
    
        if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Green`s funct ')
    
        Pc = Rad*exp(j*pnts(i))
        Ec = Centre+Pc
        zt = j * Pc * negf%spin * wght(i)/(2.d0*pi)
        negf%iE = i

        call compute_contacts(Ec,negf,i,ncyc,Tlc,Tcl,SelfEneR,GS)
  
        call calls_eq_mem_dns(negf,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer)
  
        if(negf%DorE.eq.'D') then
           call concat(TmpMt,zt,GreenR,1,1)
        endif
        if(negf%DorE.eq.'E') then
           call concat(TmpMt,zt*Ec,GreenR,1,1)
        endif
    
        call destroy(GreenR)
    
        do i1=1,ncont
           call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
        enddo
    
        if (id0.and.negf%verbose.gt.VBT) call write_clock
    
     enddo
    
     deallocate(wght)
     deallocate(pnts)
  
     ! *******************************************************************************
     ! 2. INTEGRATION OVER THE SEGMENT [Ec + j*Lambda, mu(r) + Omega+j*Lambda]
     ! (Temp /= 0) OR OVER THE CIRCLE WITH TETA FROM ZERO TO ALPHA (Temp == 0)
     ! *******************************************************************************
     ! NEW INTEGRATION FOR COMPLEX DENSITY (T=0):
     !----------------------------------------------------
     !   2  [ /           ]     1  [ /         it   ] 
     !  --- [ |  Gr(z) dz ] =   -- [ | Gr(t)Rie  dt ]  
     !  2pi [ /           ]     pi [ /              ]
     !----------------------------------------------------
     ! NEW INTEGRATION FOR COMPLEX DENSITY (T>0):
     !----------------------------------------------------
     !   g  [ /                ]      g   [ /                       ] 
     !  --- [ |  Gr(z)*f(z) dz ] =   ---  [ | Gr(z)*(z2-z1)*f(z)*dt ]  
     !  2pi [ /                ]     2*pi [ /                       ]
     !----------------------------------------------------
    
     allocate(wght(negf%Np_n(2)))
     allocate(pnts(negf%Np_n(2)))
    
     if (negf%kbT.eq.0.d0) then                        ! Circle integration T=0
    
        call  gauleg(alpha,0.d0,pnts,wght,negf%Np_n(2))
        
     else                                          ! Segment integration T>0
        
        z1 = muref + nkT + j*Lambda
        z2 = muref - nkT + j*Lambda
        
        z_diff = z2 - z1
        
        call  gauleg(1.d0,0.d0,pnts,wght,negf%Np_n(2))    !Setting weights for integration
        
     endif
  
     npid = int(negf%Np_n(2)/numprocs)
 
     istart = id*npid+1
     if(id.ne.(numprocs-1)) then 
        iend = (id+1)*npid
     else
        iend = negf%Np_n(2)
     end if
  
     do i = istart,iend
  
        if (negf%verbose.gt.VBT) then
           write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 2: point #',i,'/',iend,'  CPU=&
                &', id
        endif
  
        if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Green`s funct ')
  
        if (negf%kbT.eq.0.d0) then                      ! Circle integration T=0            
           
           Pc = Rad*exp(j*pnts(i))
           Ec = Centre+Pc
           dt = negf%spin*wght(i)/(2.d0*pi)
           zt = dt*Pc*j
          
        else                                        ! Segment integration T>0
           
           Ec = z1 + pnts(i)*z_diff
           ff = fermi_fc(Ec,muref,KbT)
           zt = negf%spin * z_diff * ff * wght(i) / (2.d0 *pi)
  
        endif
  
        negf%iE = negf%Np_n(1)+i

        call compute_contacts(Ec,negf,negf%Np_n(1)+i,ncyc,Tlc,Tcl,SelfEneR,GS)
  
        call calls_eq_mem_dns(negf,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer) 
        
        if(negf%DorE.eq.'D') then
           call concat(TmpMt,zt,GreenR,1,1)
        endif
        if(negf%DorE.eq.'E') then
           call concat(TmpMt,zt*Ec,GreenR,1,1)
        endif
  
        call destroy(GreenR)
  
        do i1=1,ncont
           call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
        enddo
  
        if (id0.and.negf%verbose.gt.VBT) call write_clock
  
     enddo
  
     deallocate(wght)
     deallocate(pnts)
  
     ! *******************************************************************************
     ! 3. SUMMATION OVER THE POLES ENCLOSED IN THE CONTOUR  (NumPoles)
     ! *******************************************************************************          
     ! NEW INTEGRATION FOR COMPLEX DENSITY (T>=0):
     !---------------------------------------------------------------------
     !             [  g                  ]          
     !  2 pi j* Res[ --- *Gr(z_k)f(z_k)  ] =  j*(-kb T)* Gr(z_k)     
     !             [ 2pi                 ]         
     !                                         (-kb*T) <- Residue
     !---------------------------------------------------------------------
  
     npid = int(NumPoles/numprocs)
     istart = id*npid+1
     if(id.ne.(numprocs-1)) then 
        iend = (id+1)*npid
     else
        iend = NumPoles
     end if
  
     do i = istart,iend
  
        if (negf%verbose.gt.VBT) then
           write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 3: point #',i,'/',iend,'  CPU=&
                &', id
        endif
  
        if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Green`s funct ')
  
        Ec = muref + j * KbT *pi* (2.d0*real(i,dp) - 1.d0)   
  
        zt= -j*negf%spin*KbT

        negf%iE = negf%Np_n(1)+negf%Np_n(2)+i
  
        call compute_contacts(Ec,negf,negf%Np_n(1)+negf%Np_n(2)+i,ncyc,Tlc,Tcl,SelfEneR,GS)
  
        call calls_eq_mem_dns(negf,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer)
  
        if(negf%DorE.eq.'D') then
           call concat(TmpMt,zt,GreenR,1,1)
        endif
        if(negf%DorE.eq.'E') then
           call concat(TmpMt,zt*Ec,GreenR,1,1)
        endif
  
        call destroy(GreenR)   
  
        do i1=1,ncont
           call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
        enddo
  
        if (id0.and.negf%verbose.gt.VBT) call write_clock
  
     enddo
  
     if(negf%DorE.eq.'D') then
        call zspectral(TmpMt,TmpMt,0,negf%rho)
     endif
     if(negf%DorE.eq.'E') then
        call zspectral(TmpMt,TmpMt,0,negf%rho_eps)
     endif
  
  
     call destroy(TmpMt)


  end subroutine contour_int

  !--------------------------------------------!
  !--------------------------------------------!
  ! Non equilibrium integration over real axis !
  !--------------------------------------------!
  !--------------------------------------------!
  !-----------------------------------------------------------------------
  ! Contour integration for density matrix 
  ! DOES INCLUDE FACTOR 2 FOR SPIN !! 
  !-----------------------------------------------------------------------
  subroutine real_axis_int(negf)

    type(Tnegf) :: negf
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    type(z_CSR) :: GreenR, TmpMt 

    integer :: npid, istart, iend, ref
    integer :: i, i1, ioffset, outer, ncont, j1, npT

    real(dp), DIMENSION(:), allocatable :: wght,pnts   ! Gauss-quadrature points
    real(dp), DIMENSION(:), allocatable :: frm_f
    real(dp) :: Omega, mumin, mumax
    real(dp) :: ncyc, kbT, dt

    complex(dp) :: zt
    complex(dp) :: Ec

    ncont = negf%str%num_conts
    kbT = negf%kbT
    ref = negf%refcont
    ioffset = negf%Np_n(1) + negf%Np_n(2) + negf%n_poles
 
    mumin=minval(negf%Efermi(1:ncont)-negf%mu(1:ncont))
    mumax=maxval(negf%Efermi(1:ncont)-negf%mu(1:ncont))

    if (mumax.gt.mumin) then
       
       Omega = negf%n_kt * kbT

       outer = negf%outer
       
       call log_allocate(frm_f,ncont)
       
       call create(TmpMt,negf%H%nrow,negf%H%ncol,negf%H%nrow)
       call initialize(TmpMt)
       

       !Compute extended number of points due to kT
       !npT=nint(negf%Np_real(1)*Omega/(mumax-mumin))
       ! **** MODIFIED 4/04/2012 because for mumin-mumax -> 0 does not work !!!
       npT = 0
       
       allocate(pnts(negf%Np_real(1)+2*npT))
       allocate(wght(negf%Np_real(1)+2*npT))

       !Setting weights for gaussian integration
        
       call gauleg(mumin-Omega,mumax+Omega,pnts,wght,negf%Np_real(1)+2*npT)

       !Computing real axis integral       
       npid = int((negf%Np_real(1)+2*npT)/numprocs)
       istart = id*npid+1
       if(id.ne.(numprocs-1)) then 
          iend = (id+1)*npid
       else
          iend = negf%Np_real(1)+2*npT
       end if

       !---------------------------------------------------------
       !    g    --  [ /                                      ]
       ! ------  >   [ | j [f(i)-f(ref)] Gr(E) Gam_i Ga(E) dE  ]  
       ! 2*pi*j  --i [ /                                      ]
       !---------------------------------------------------------

       do i = istart,iend

          if (negf%verbose.gt.VBT) then
             write(6,'(a17,i3,a1,i3,a6,i3,f8.4)') 'INTEGRAL neq: pnt #',i,'/',iend,'  CPU=&
                  &', id, pnts(i)
          endif

          do j1 = 1,ncont
             frm_f(j1)=fermi_f(pnts(i),negf%Efermi(j1)-negf%mu(j1),KbT)
          enddo

          if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Green`s funct ')

          Ec = cmplx(pnts(i),negf%delta,dp)

          dt = negf%wght * negf%spin * wght(i)/(2*pi)

          zt = dt*(1.d0,0.d0)

          call compute_contacts(Ec,negf,ioffset+i,ncyc,Tlc,Tcl,SelfEneR,GS)

          call calls_neq_mem_dns(negf,real(Ec),SelfEneR,Tlc,Tcl,GS,negf%str,frm_f,GreenR,outer)

          do i1=1,ncont
             call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
          enddo

          if(negf%DorE.eq.'D') then
             call concat(TmpMt,zt,GreenR,1,1)
          endif
          if(negf%DorE.eq.'E') then
             call concat(TmpMt,zt*Ec,GreenR,1,1)
          endif

          call destroy(GreenR) 

          if (id0.and.negf%verbose.gt.VBT) call write_clock

       enddo

       deallocate(wght)
       deallocate(pnts)

       if(negf%DorE.eq.'D') then
          if(allocated(negf%rho%nzval)) then
             call concat(negf%rho,TmpMt,1,1)
          else
             call clone(TmpMt,negf%rho) 
          endif
       endif
       if(negf%DorE.eq.'E') then
          if(allocated(negf%rho_eps%nzval)) then
             call concat(negf%rho_eps,TmpMt,1,1)
          else
             call clone(TmpMt,negf%rho_eps) 
          endif
       endif
       
       call destroy(TmpMt)
       
       call log_deallocate(frm_f)

    end if


  end subroutine real_axis_int


!-----------------------------------------------------------------------
! Contour integration for density matrix 
! DOES INCLUDE FACTOR 2 FOR SPIN !! 
!-----------------------------------------------------------------------
  subroutine real_axis_int_n(negf)

    type(Tnegf) :: negf
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    type(z_CSR) :: GreenR, TmpMt 

    integer :: npid, istart, iend, ref
    integer :: i, i1, ioffset, outer, ncont, j1, npT

    real(dp), DIMENSION(:), allocatable :: wght,pnts   ! Gauss-quadrature points
    real(dp), DIMENSION(:), allocatable :: frm_f
    complex(dp), dimension(:), allocatable :: q_tmp
    real(dp) :: Omega, mumin, mumax
    real(dp) :: ncyc, kbT, dt

    complex(dp) :: zt
    complex(dp) :: Ec

    ncont = negf%str%num_conts
    kbT = negf%kbT
    ref = negf%refcont
    ioffset = negf%Np_n(1) + negf%Np_n(2) + negf%n_poles

    mumin=minval(negf%Efermi(1:ncont)-negf%mu(1:ncont))
    mumax=maxval(negf%Efermi(1:ncont)-negf%mu(1:ncont))

    if (negf%writeLDOS) then
       open(2001,file=trim(negf%out_path)//'LDOS.dat')
       open(2002,file=trim(negf%out_path)//'energy.dat')
       call log_allocate(q_tmp, negf%H%nrow)
    endif
       
       Omega = negf%n_kt * kbT
       outer = negf%outer
       
       call log_allocate(frm_f,ncont+1)
       frm_f = 0.d0
       
       call create(TmpMt,negf%H%nrow,negf%H%ncol,negf%H%nrow)
       call initialize(TmpMt)
       

       !Compute extended number of points due to kT
       !npT=nint(negf%Np_real(1)*Omega/(mumax-mumin))
       ! **** MODIFIED 4/04/2012 because for mumin-mumax -> 0 does not work !!!
       npT = 0
       
       allocate(pnts(negf%Np_real(1)+2*npT))
       allocate(wght(negf%Np_real(1)+2*npT))

       !Setting weights for gaussian integration
       call gauleg(mumin-Omega,mumax+Omega,pnts,wght,negf%Np_real(1)+2*npT)

       !Computing real axis integral       
       npid = int((negf%Np_real(1)+2*npT)/numprocs)
       istart = id*npid+1
       if(id.ne.(numprocs-1)) then 
          iend = (id+1)*npid
       else
          iend = negf%Np_real(1)+2*npT
       end if

       !---------------------------------------------------------------
       !    g    --  [ /                                              ]
       ! ------  >   [ | j [f(i) - f(ref)] Gr(E) Gam_i Ga(E) dE       ]  
       ! 2*pi*j  --i [ /                                              ]
       !---------------------------------------------------------------

       do i = istart,iend

          if (negf%verbose.gt.VBT) then
             write(6,'(a19,i3,a1,i3,a6,i3,f8.4)') 'INTEGRAL neq: pnt #',i,'/',iend,'  CPU=&
                  &', id, pnts(i)
          endif

          do j1 = 1,ncont
             frm_f(j1)=fermi_f(pnts(i),negf%Efermi(j1)-negf%mu(j1),KbT)
          enddo

          Ec = cmplx(pnts(i),negf%delta,dp)

          dt = negf%wght * negf%spin * wght(i)/(2*pi)

          zt = dt*(1.d0,0.d0)


          if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Green`s funct ')

          call compute_contacts(Ec,negf,ioffset+i,ncyc,Tlc,Tcl,SelfEneR,GS)

          call calls_neq_mem_dns(negf,real(Ec),SelfEneR,Tlc,Tcl,GS,negf%str,frm_f,GreenR,outer)
       
          do i1=1,ncont
             call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
          enddo

          if(negf%DorE.eq.'D') then
             call concat(TmpMt,zt,GreenR,1,1)
          endif
          if(negf%DorE.eq.'E') then
             call concat(TmpMt,zt*Ec,GreenR,1,1)
          endif

          if (negf%writeLDOS) then
             call zgetdiag(GreenR, q_tmp) 
             do i1 = 1,negf%str%central_dim
                write(2001,'((ES14.5))', advance='NO') -real(q_tmp(i1))
             enddo
             write(2001,*)
             write(2002,*) pnts(i)
          endif


          call destroy(GreenR) 

          if (id0.and.negf%verbose.gt.VBT) call write_clock

       enddo

       deallocate(wght)
       deallocate(pnts)

       if(negf%DorE.eq.'D') then
          if(allocated(negf%rho%nzval)) then
             call concat(negf%rho,TmpMt,1,1)
          else
             call clone(TmpMt,negf%rho) 
          endif
       endif
       if(negf%DorE.eq.'E') then
          if(allocated(negf%rho_eps%nzval)) then
             call concat(negf%rho_eps,TmpMt,1,1)
          else
             call clone(TmpMt,negf%rho_eps) 
          endif
       endif
       
       call destroy(TmpMt)
       
       call log_deallocate(frm_f)

    if (negf%writeLDOS) then
       close (2001)
       close (2002)
       call log_deallocate(q_tmp)
    endif

  end subroutine real_axis_int_n
  !------------------------------------------------------------------------
 
  subroutine init_emesh(negf, reflevel)
    type(Tnegf) :: negf
    integer :: reflevel

    real(dp) :: Emin, Emax, nkT, Wmax, mumin, mumax
    integer :: np,ncont,ioffset


    if (negf%elph%numselmodes .eq. 0) return

    ncont = 2 !TEMPORARY HACKING!!!!!
    !ncont = negf%str%num_conts
    
    !mumin = negf%Ec
    Wmax = maxval(negf%elph%Wq, negf%elph%selmodes)
    nkT = negf%n_kt * negf%kbT
    
    Emin=minval(negf%Efermi(1:ncont)-negf%mu(1:ncont)) - nkT - Wmax
    Emax=maxval(negf%Efermi(1:ncont)-negf%mu(1:ncont)) + nkT + Wmax


    mumin = Emin - negf%elph%scba_iterations * Wmax
    mumax = Emax + negf%elph%scba_iterations * Wmax

    np = negf%Np_real(1)
    
    if (np.eq.0) np=1 

    !CREATE energy mesh
    ioffset = negf%Np_n(1) + negf%Np_n(2) + negf%n_poles

    print*,'create emesh:',mumin,mumax,np

    call create_mesh(negf%emesh, reflevel, mumin, mumax, np, ioffset)

    negf%elph%Erange(1) = Emin 
    negf%elph%Erange(2) = Emax 

    negf%elph%self_range(1) = mumin + Wmax
    negf%elph%self_range(2) = mumax - Wmax

  end subroutine init_emesh

  !-----------------------------------------------------------------------------------------------------
  subroutine destroy_emesh(negf)
    type(Tnegf) :: negf

    call destroy_mesh(negf%emesh)

  end subroutine destroy_emesh

  !---------------------------------------------------------------------------
  subroutine test(negf)
    type(Tnegf) :: negf
    type(mesh) :: emesh
    type(elem) :: el
    type(z_CSR) :: Dens, Dens1
    type(z_CSR), target :: P
    type(z_CSR), pointer :: P2
    type(TG_pointer) :: G(3)

    integer :: i,j,lev,nelem, newelem

    print*,'Create mesh'
    call create_mesh(emesh,6, -1.0_dp, 0.0_dp, 5, 0)
    print*,'Done'

    i=1
    do while(i.le.emesh%maxind) 
      nelem = emesh%maxind
      print*,'elem',i, emesh%pactive(i)%pelem%map
      lev = 1
      
      call adaptive_int(lev,negf,emesh,emesh%pactive(i)%pelem,G,Dens1)

      newelem = emesh%maxind - nelem 
      print*,'new elements in mesh',newelem
      i = i + newelem + 1

      !print*,allocated(Gless(1)%nzval),allocated(Gless(2)%nzval),allocated(Gless(3)%nzval)
      call destroyp(G(1)%pG)
      call destroyp(G(2)%pG)
      G(1)%pG => G(3)%pG
      G(1)%ind = G(3)%ind
      G(3)%pG => null()
      G(3)%ind = 0

      call destroy(Dens1)
    enddo

    if (allocated(G(1)%pG%nzval)) then
       print*," Destroy",G(1)%ind
       call destroyp(G(1)%pG)
    endif   

  end subroutine test
  !---------------------------------------------------------------------------
  

  recursive subroutine adaptive_int(lev,negf,emesh,el,G,Dens)
    integer, intent(inout) :: lev
    type(TNegf) :: negf
    type(mesh) :: emesh
    type(elem), pointer :: el 
    type(TG_pointer) :: G(3)
    type(z_CSR) :: Dens
    
    ! local variables
    type(elem), pointer :: pel
    !type(z_CSR), target :: Gless2(3)
    type(TG_pointer) :: G2(3)
    type(z_CSR) :: Dens1, Dens2
    integer :: i
    real(dp) :: dE
    character(10) :: al
    logical :: dorefine

    al=''
    do i=1,lev
      al = trim(al)//'o'
    enddo

    print*,trim(al)//' Sub called at lev:',lev

    do i = 1,3
      if (associated(G(i)%pG) .and. allocated(G(i)%pG%nzval)) then
         G(i)%ind = el%map(i) 
         print *,trim(al)//" Gp: ",abs(trace(G(i)%pG)) !G(i)%ind,LOC(G(i)%pG)
      else   
         allocate(G(i)%pG)
         !!print*,'E=',el%pnt(i)*negf%eneconv
         call compute_gless(negf,el%pnt(i),el%map(i),G(i)%pG)
         G(i)%ind = el%map(i) 
         print *,trim(al)//" Create G: ",abs(trace(G(i)%pG))  !G(i)%ind,LOC(G(i)%pG)
      endif   
    enddo
    
    call trapez23(negf, G, el%pnt, Dens, el%error, dorefine)

    if (dorefine .and. lev.lt.emesh%maxreflev) then

       print*,'refine element', el%ind, 'error',el%error

       pel=>el !Workaround since after refine el points to el%child1 ?!?
       call refine(emesh,el)
       !print*, trim(al)//" (el,ch1,ch2):",LOC(pel),LOC(pel%child1),LOC(pel%child2)
 
       !print*,trim(al)//" New mesh"
       !do i = 1, emesh%maxind
       !    print*,trim(al),LOC(emesh%pactive(i)%pelem), emesh%pactive(i)%pelem%map
       !enddo
 
       !print*,trim(al)//" Go to elem",LOC(pel%child1)
       G2(1)%pG => G(1)%pG 
       G2(3)%pG => G(2)%pG 
       G2(2)%pG => null() 
 
       lev = lev + 1
       call adaptive_int(lev,negf,emesh,pel%child1,G2,Dens1)  
       lev = lev - 1

       call destroyp(G2(2)%pG)
       
       !print*,trim(al)//" Go to elem",LOC(pel%child2)
       G2(3)%pG => G(3)%pG
       G2(1)%pG => G(2)%pG
       G2(2)%pG => null() 
 
       lev = lev + 1
       call adaptive_int(lev,negf,emesh,pel%child2,G2,Dens2)  
       lev = lev - 1
  
       call destroyp(G2(2)%pG)
 
       call destroy(Dens)
       call create(Dens,negf%H%nrow,negf%H%ncol,negf%H%nrow)
       call initialize(Dens)

       call concat(Dens,Dens1,1,1)
       call concat(Dens,Dens2,1,1)

       call destroy(Dens1,Dens2)
    
       print*,'--------------------------------------------------------------'

    else

    !   if (lev.gt.1) then
    !      call destroyp(G(2)%pG)
    !      print*,trim(al)//" Destroy",G(2)%ind
    !   endif

    !if (negf%writeLDOS) then
    !   do j1 = 1, 2      
    !      call zgetdiag(Gless(j1), q_tmp)
    !      do i1 = 1,negf%str%central_dim
    !         write(2001,'((ES14.5))', advance='NO') real(q_tmp(i1))
    !      enddo
    !      write(2001,*)
    !      write(2002,*) pel%pnt(j1)
    !   enddo   
    endif
 
      

  end subroutine adaptive_int

  !---------------------------------------------------------------------------
  subroutine trapez23(negf,G, Epnt, Dens, error, refine)
    type(TNegf) :: negf
    type(TG_pointer) :: G(3)
    real(dp), intent(in) :: Epnt(3)
    type(z_CSR) :: Dens
    real(dp), intent(out) :: error
    logical, intent(out) :: refine
    
    type(z_CSR) :: I2
    real(dp) :: dt    
    complex(dp) :: zt, tr, tr2

    call create(I2,negf%H%nrow,negf%H%ncol,negf%H%nrow)
    call initialize(I2)

    dt = negf%wght * negf%spin * (Epnt(2)-Epnt(1)) / (2*pi)
    zt = dt*(1.d0,0.d0)
    ! Computing I2
    if(negf%DorE.eq.'D') then
       call concat(I2,zt,G(1)%pG,1,1)
       call concat(I2,zt,G(3)%pG,1,1)
    endif
    if(negf%DorE.eq.'E') then
       call concat(I2,zt*Epnt(1),G(1)%pG,1,1)
       call concat(I2,zt*Epnt(3),G(3)%pG,1,1)
    endif
    
    call clone(I2,Dens)

    ! Computing I3
    Dens%nzval = Dens%nzval/2.0_dp 

    if(negf%DorE.eq.'D') then
        call concat(Dens,zt,G(2)%pG,1,1)
    endif
    if(negf%DorE.eq.'E') then
        call concat(Dens,zt*Epnt(2),G(2)%pG,1,1)
    endif
    
    tr2 = trace(I2)
    tr = trace(Dens) 

    if (abs(tr2).lt.negf%int_acc) then
            error = abs(tr - tr2)
    else
            error = abs(tr - tr2)/abs(tr2)
    endif


    call destroy(I2)
    refine = .false.
    if (error > negf%int_acc) then
       refine = .true.
    endif

  end subroutine trapez23
  !---------------------------------------------------------------------------
  
  subroutine compute_gless(negf,E,pnt,Gless)
    type(TNegf) ::  negf
    real(dp) :: E
    integer :: pnt
    type(z_CSR) :: Gless

    real(dp), DIMENSION(:), allocatable :: frm_f
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    integer :: i1, j1, ref, ncont, outer, iter
    complex(dp) :: Ec
    real(dp) :: ncyc


    call log_allocate(frm_f,ncont+1)
    frm_f = 0.d0
    ncont = negf%str%num_conts
    outer = negf%outer 
    iter = negf%elph%scba_iter

    if (negf%iteration.eq.1 .and. iter.eq.0) then
        negf%ReadOldSGF = 2 ! Compute & save Surface GF
    else
        negf%ReadOldSGF = 0 ! Reload Surface GF
    endif

    Ec = cmplx(E,negf%delta,dp)

    negf%iE = pnt    

    call  compute_contacts(Ec,negf,pnt,ncyc,Tlc,Tcl,SelfEneR,GS)

    do j1 = 1,ncont
       frm_f(j1)=fermi_f(E,negf%Efermi(j1)-negf%mu(j1),negf%KbT)
    enddo
    
    if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Green`s funct ')

    ! Compute Gr and G<
    call calls_neq_mem_dns(negf,E,SelfEneR,Tlc,Tcl,GS,negf%str,frm_f,Gless,outer)

    if (id0.and.negf%verbose.gt.VBT) call write_clock

    do i1=1,ncont
       call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
    enddo

    call log_deallocate(frm_f)

  end subroutine compute_gless
  
  !----------------------------------------------------------------------------
  ! Gauss - Legendre quadrature weights and points
  !----------------------------------------------------------------------------
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
    
  end subroutine gauleg

  !--------------------------------------------------------------------  
  subroutine trapez(x1,x2,x,w,n)

    real(kind=dp), PARAMETER :: ACC = 1d-15

    INTEGER n
    real(kind=dp) :: x1,x2,x(n),w(n)

    real(dp) :: d
    integer :: i

    d = (x2-x1)/(n-1)

    w = d * 1.0_dp
    w(1) = d * 0.5_dp
    w(n) = d * 0.5_dp

    do i = 1, n
       x(i) = ( x1*(n-i) + x2*(i-1) ) / (n-1) 
    enddo
 
  end subroutine trapez
  !--------------------------------------------------------------------  

  !------------------------------------------------------------------------------- 
  subroutine tunneling_and_current(negf)
    
    type(Tnegf) :: negf

    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    
    Real(dp), Dimension(:), allocatable :: TUN_MAT
    !Real(kind=dp), Dimension(:,:,:), allocatable :: TUN_PMAT, TUN_TOT_PMAT      
    
    Real(dp), Dimension(:), allocatable :: mumin_array, mumax_array
    Real(dp), Dimension(:), allocatable :: LEDOS
    Real(dp) :: ncyc, Ec_min, mumin, mumax
    
    Complex(dp) :: Ec
    
    Integer, Dimension(:), pointer :: cblk, indblk
    Integer :: i, Nstep, npid, istart, iend, i1
    Integer :: size_ni, size_nf, icpl, ncont, icont
    Integer :: nbl
    Integer :: iLDOS
    
    Character(6) :: ofKP
    Logical :: do_LEDOS, lex

    do_LEDOS = .false.
    if(negf%nLDOS.gt.0) do_LEDOS=.true.
    
    nbl = negf%str%num_PLs
    ncont = negf%str%num_conts
    cblk => negf%str%cblk
    indblk => negf%str%mat_PL_start
    
    Nstep=NINT((negf%Emax-negf%Emin)/negf%Estep)
    npid = int((Nstep+1)/numprocs)
    
    !Get out if Emax<Emin and Nstep<0
    if (Nstep.lt.0) then
       if(id0) write(*,*) '0 tunneling points;  current = 0.0'
       call log_allocatep(negf%tunn_mat,0,0)
       call log_allocatep(negf%currents,1)
       negf%currents(1) = 0.0_dp 
       if (do_ledos) call log_allocatep(negf%ldos_mat,0,0)
       return
    endif
    
    !Extract Contacts in main
    !Tunneling set-up
    do i=1,size(negf%ni)
       if (negf%ni(i).eq.0) then
          size_ni=i-1
          exit
       endif
    enddo
    
    do i=1,size(negf%nf)
       if (negf%nf(i).eq.0) then
          size_nf=i-1
          exit
       endif
    enddo

    !check size_ni .ne. size_nf
    if (size_ni.ne.size_nf) then 
       size_ni=min(size_ni,size_nf)
       size_nf=min(size_ni,size_nf)
    endif
    
    call log_allocate(mumin_array,size_ni)
    call log_allocate(mumax_array,size_ni)
    
    ! find bias window for each contact pair
    do icpl=1,size_ni
       mumin_array(icpl)=min(negf%Efermi(negf%ni(icpl))-negf%mu(negf%ni(icpl)),&
            negf%Efermi(negf%nf(icpl))-negf%mu(negf%nf(icpl)))
       mumax_array(icpl)=max(negf%Efermi(negf%ni(icpl))-negf%mu(negf%ni(icpl)),&
            negf%Efermi(negf%nf(icpl))-negf%mu(negf%nf(icpl)))
    enddo
    
    ncyc=0
    istart = 1
    iend = npid
    
    call log_allocate(TUN_MAT,size_ni)
    call log_allocatep(negf%tunn_mat,Nstep+1,size_ni)   
    !call log_allocate(TUN_PMAT,npid,size_ni,num_channels) 
    !call log_allocate(TUN_TOT_PMAT,Nstep+1,size_ni,num_channels) 
    negf%tunn_mat = 0.0_dp 

    if (do_LEDOS) then
       call log_allocatep(negf%ldos_mat,Nstep+1,negf%nLDOS)
       call log_allocate(LEDOS,negf%nLDOS)          
       negf%ldos_mat(:,:)=0.d0
    endif
    
    
    !Loop on energy points: tunneling 
    do i1 = istart,iend
       
       Ec_min = negf%Emin + id*npid*negf%Estep
       Ec = (Ec_min + negf%Estep*(i1-1))*(1.d0,0.d0) !+negf%delta*(0.d0,1.d0) 

       if (negf%verbose.gt.VBT) then
          write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL:   point #',i1,'/',iend,'  CPU=&
             &', id
       endif

       if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Contact SE ')       

       call compute_contacts(Ec+negf%delta*(0.d0,1.d0),negf,i1,ncyc,Tlc,Tcl,SelfEneR,GS)

       if (id0.and.negf%verbose.gt.VBT) call write_clock
       
       do icont=1,ncont
          call destroy(Tlc(icont))
          call destroy(Tcl(icont))
       enddo

       if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Tunneling ') 
       
       if (.not.do_LEDOS) then
          
          call tunneling_dns(negf%HM,negf%SM,Ec,SelfEneR,negf%ni,negf%nf,size_ni, &
               negf%str,TUN_MAT)
          negf%tunn_mat(i1,:) = TUN_MAT(:) * negf%wght
          
       else
          
          LEDOS(:) = 0.d0
          
          call tun_and_dos(negf%HM,negf%SM,Ec,SelfEneR,GS,negf%ni,negf%nf,negf%nLDOS, &
               negf%LDOS,size_ni,negf%str,TUN_MAT,LEDOS)
          
          negf%tunn_mat(i1,:) = TUN_MAT(:) * negf%wght
          negf%ldos_mat(i1,:) = LEDOS(:) * negf%wght
          
       endif

       if (id0.and.negf%verbose.gt.VBT) call write_clock
       
       do icont=1,ncont
          call destroy(SelfEneR(icont))
          call destroy(GS(icont))
       enddo
       
    enddo !Loop on energy i1 = istart,iend

    !---------------------------------------------------------------------------
    !   COMPUTATION OF CURRENTS 
    !---------------------------------------------------------------------------
    call log_allocatep(negf%currents,size_ni)
    negf%currents(:)=0.d0
    
    do icpl=1,size_ni
       
       mumin=mumin_array(icpl)
       mumax=mumax_array(icpl)

       !checks if the energy interval is appropriate
       !if (id0.and.mumin.lt.mumax) then
               
        !if (negf%Emin.gt.mumin-10*negf%kbT .or. negf%Emax.lt.mumin-10*negf%kbT) then
        !  write(*,*) 'WARNING: the interval Emin..Emax is smaller than the bias window'
        !  write(*,*) 'Emin=',negf%emin*negf%eneconv, 'Emax=',negf%emax*negf%eneconv
        !  write(*,*) 'kT=',negf%kbT*negf%eneconv    
        !  write(*,*) 'Suggested interval:', &
        !        (mumin-10*negf%kbT)*negf%eneconv,(mumax+10*negf%kbT)*negf%eneconv
        ! endif
       !endif
       
       negf%currents(icpl)= integrate(negf%tunn_mat(:,icpl),mumin,mumax,negf%kbT, &
            negf%Emin,negf%Emax,negf%Estep)
       
    enddo
   
    if (negf%verbose.gt.VBT) then 
     do i1=1,size_ni
       write(*,'(1x,a,i3,i3,a,i3,a,ES14.5,a,ES14.5,a)') 'contacts:',negf%ni(i1),negf%nf(i1), &
            '; k-point:',negf%kpoint,'; current:', negf%currents(i1),' A'
     enddo
    endif

    call log_deallocate(TUN_MAT)
    !call log_deallocate(TUN_PMAT)
    !call log_deallocate(TUN_TOT_PMAT)  
    
    call log_deallocate(mumin_array)
    call log_deallocate(mumax_array)
    
    if(do_LEDOS) then
       call log_deallocate(LEDOS)
    endif
    
    
  end subroutine tunneling_and_current
 

!////////////////////////////////////////////////////////////////////////
!************************************************************************
!
! Function to integrate the tunneling and get the current
!
!************************************************************************

function integrate(TUN_TOT,mumin,mumax,kT,emin,emax,estep)

  implicit none

  real(dp) :: integrate
  real(dp), intent(in) :: mumin,mumax,emin,emax,estep
  real(dp), dimension(:), intent(in) :: TUN_TOT
  real(dp), intent(in) :: kT 

  REAL(dp) :: destep,kbT,TT1,TT2,E3,E4,TT3,TT4
  REAL(dp) :: E1,E2,c1,c2,curr
  INTEGER :: i,i1,N,Nstep,imin,imax
  
  curr=0.d0
  N=0
  destep=1.0d10 
  Nstep=NINT((emax-emin)/estep);

  if (kT.lt.3.d-6) then
    kbT = 1.d-5
  else
    kbT = kT     
  endif
  ! Find initial step for integration

  imin=0
  do i=0,Nstep
     E1=emin+estep*i     
     imin=i-1
     if(E1.ge.mumin-10*kbT) then 
        exit
     endif
  enddo

  ! Find final step for integration 
  imax=0
  do i=Nstep,imin,-1    
     E1=emin+estep*i 
     imax=i+1
     if(E1.le.mumax+10*kbT) then 
        exit
     endif
  enddo


  !rest the min and max to the actual interval
  if (imin.lt.0) imin=0
  if (imax.gt.Nstep) imax=Nstep 

  ! performs the integration with simple trapezium rule. 
  do i=imin,imax-1
     
     E1=emin+estep*i  
     TT1=TUN_TOT(i+1)     
     E2=emin+estep*(i+1)
     TT2=TUN_TOT(i+2)
     
     ! Each step is devided into substeps in order to
     ! smooth out the Fermi function
     do while (destep.ge.2*kbT) 
        N=N+1
        destep=(E2-E1)/N
     enddo
     
     ! within each substep the tunneling is linearly 
     ! interpolated
     do i1=0,N-1
        
        E3=E1+(E2-E1)*i1/N
        E4=E3+(E2-E1)/N
        TT3=( TT2-TT1 )*i1/N + TT1
        TT4=TT3 + (TT2-TT1)/N
        
        c1=2.d0*eovh*(fermi_f(E3,mumax,KbT)-fermi_f(E3,mumin,KbT))*TT3
        c2=2.d0*eovh*(fermi_f(E4,mumax,KbT)-fermi_f(E4,mumin,KbT))*TT4
        
        curr=curr+(c1+c2)*(E4-E3)/2.d0
        
     enddo
     
  enddo
  
  integrate = curr
  
end function integrate
         
end module integrations 
