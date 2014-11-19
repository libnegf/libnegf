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
 use distributions 
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

 public :: contour_int       ! generalized contour integrator 
 public :: real_axis_int     ! real-axis integrator
 public :: ldos_int          ! ldos only integrator

 public :: contour_int_def   ! contour integration for DFT(B)
 public :: contour_int_n_def ! contour integration for CB
 public :: contour_int_p_def ! contour integration for VB
 public :: real_axis_int_def ! real axis integration 
 public :: real_axis_int_n_def ! integration of CB on real axis
 public :: real_axis_int_p_def ! integration of VB on real axis

 public :: tunneling_int_def  !
 public :: tunneling_and_dos  ! computes of T(E) & LDOS(E)
 public :: electron_current   ! computes terminal currents

 public :: phonon_tunneling   ! computes T(E) for phonons
 public :: phonon_current     ! computes heat currents
 public :: thermal_conductance ! computes thermal conductance

 public :: integrate_el       ! integration of tunneling (el)
 public :: integrate_ph       ! integration of tunneling (ph)
 !!public :: compute_dos      ! compute local dos only

 !public :: menage_scratch
 ! ////////////////////////////////////////////////////////////
 ! Under development:
 !public :: init_emesh, destroy_emesh
 !private :: adaptive_int, trapez23 
 !public :: contour_int_ph, real_axis_int_ph, real_axis_int_ph2
 !
 !type TG_pointer
 !  type(z_CSR),  pointer :: pG => null()   
 !  integer :: ind   
 !end type TG_pointer
 ! ////////////////////////////////////////////////////////////

 integer, PARAMETER :: VBT=70

 !! Structure used to define energy points for the integration
 !! For every point we define
 !!     path (1,2 or 3): the energy point belongs to a real axis 
 !!     integration (1), a complex plane integration (2) or a 
 !!     pole summation (3)
 !!     pt_path: relative point number within a single path
 !!     pt: absolute point number along the whole integration path
 !!     cpu: cpu assigned to the calculation of the given energy point
 !!     Ec: energy value
 !!     wght: a weight used in final summation to evaluate integrals
 type TEnGrid   
     integer :: path
     integer :: pt_path
     integer :: pt
     integer :: cpu
     complex(dp) :: Ec
     complex(dp) :: wght
 end type TEnGrid
     
 type(TEnGrid), dimension(:), allocatable :: en_grid

contains
  
  subroutine destroy_en_grid()
    if (allocated(en_grid)) deallocate(en_grid)
  end subroutine    

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine write_info(verbose,message,Npoints)
    integer, intent(in) :: verbose
    character(*), intent(in) :: message
    integer, intent(in) :: Npoints

     if (id0 .and. verbose.gt.30) then
       write(6,'(a,a,i0,a,i0)') message,', CPU ',id,' points ',Npoints
     end if

  end subroutine write_info
  !-----------------------------------------------------------------------
  subroutine write_point(verbose,gridpn,Npoints)
    integer, intent(in) :: verbose
    type(TEnGrid), intent(in) :: gridpn
    integer, intent(in) :: Npoints

    if (verbose.gt.VBT) then
      write(6,'(3(a,i0),a,ES15.8)') 'INTEGRAL: point # ',gridpn%pt, &
          &'/',Npoints,'  CPU= ', id, '  E=',real(gridpn%Ec)
    endif

  end subroutine write_point
               
  !-----------------------------------------------------------------------
  ! Projected DOS on atoms or obitals
  !-----------------------------------------------------------------------
  subroutine ldos_int(negf) 
    type(Tnegf) :: negf

    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    Type(z_CSR) ::  Gr
    complex(dp), Dimension(:), ALLOCATABLE :: diag

    integer :: Nstep, i, i1, l, kb, ke
    integer :: outer, ncont

    real(dp) :: ncyc
    complex(dp) :: Ec
    character(6) :: ofKP
    character(1) :: ofSp

    outer = 1
    ncont = negf%str%num_conts

    Nstep = size(en_grid)
    
    call log_allocatep(negf%ldos_mat,Nstep,negf%nLDOS)
    negf%ldos_mat(:,:)=0.d0
    
    do i = 1, Nstep
  
       if (en_grid(i)%cpu /= id) cycle
      
       Ec = en_grid(i)%Ec+(0.d0,1.d0)*negf%dos_delta
       negf%iE = en_grid(i)%pt

       call write_point(negf%verbose,en_grid(i), size(en_grid))

       call compute_contacts(Ec,negf,ncyc,Tlc,Tcl,SelfEneR,GS)
    
       call calls_eq_mem_dns(negf,Ec,SelfEneR,Tlc,Tcl,GS,Gr,negf%str,outer)
       
       do i1=1,ncont
          call destroy(Tlc(i1),Tcl(i1))
       enddo

       do i1=1,ncont
          call destroy(SelfEneR(i1),GS(i1))
       enddo

       call log_allocate(diag, Gr%nrow)
       call getdiag(Gr,diag)


       do i1 = 1, size(negf%LDOS)
           negf%ldos_mat(i, i1) = - aimag( sum(diag(negf%LDOS(i1)%indexes)) )/pi
       enddo
        
       call destroy(Gr)
       call log_deallocate(diag)

    enddo
      
    call destroy_en_grid()

  end subroutine ldos_int
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Contour integration for density matrix 
  ! DOES INCLUDE FACTOR 2 FOR SPIN !! 
  !-----------------------------------------------------------------------
  ! Performs the complex contur integration
  !
  !      T>=0 
  !                +
  !                + 
  !       * * * * * * * * * * * 
  !       *        +                 
  !  --- -*---========================
  !      Elow Ec  muref       
  ! 
  !-----------------------------------------------------------------------
  ! Contour integration for density matrix 
  ! DOES INCLUDE FACTOR 2 FOR SPIN !! 
  !-----------------------------------------------------------------------

  subroutine contour_int_n_def(negf)
    type(Tnegf) :: negf  
    
    integer :: i, Ntot, Npoles, ioffs
    real(dp), DIMENSION(:), ALLOCATABLE :: wght,pnts   ! Gauss-quadrature points
    real(dp) :: Omega, Lambda
    real(dp) :: muref, kbT, Emin

    complex(dp) :: z1,z2,z_diff, zt
    complex(dp) :: Ec, ff

    kbT = negf%kbT(negf%refcont)
    muref = negf%muref
    
    Omega = negf%n_kt * kbT

    if (negf%n_poles.eq.0) then
      Lambda = 0.5d0* kbT * pi
    else
      Lambda = 2.d0* negf%n_poles * KbT * pi
    endif

    Emin = negf%Ec + negf%DeltaEc

    if ((Emin < (muref + 1.d-3)) .and. &
        (Emin > (muref - 1.d-3))) then
       Emin = muref - kbT
    endif

    Npoles = negf%n_poles
    if (Emin > muref) then
      Npoles = 0
    endif

    Ntot = negf%Np_n(1) + negf%Np_n(2) + Npoles
    allocate(en_grid(Ntot))

    ! *******************************************************************************
    ! 1. INTEGRATION OVER THE SEGMENT [Ec - dEc , Ec - dEc + j*Lambda]
    ! *******************************************************************************
    ! NEW INTEGRATION FOR COMPLEX DENSITY (T>0):
    !----------------------------------------------------
    !   g  [ /                ]      g   [ /                       ] 
    !  --- [ |  Gr(z)*f(z) dz ] =   ---  [ | Gr(z)*(z2-z1)*f(z)*dt ]  
    !  2pi [ /                ]     2*pi [ /                       ]
    !----------------------------------------------------

    z1 = Emin
    z2 = Emin + j*Lambda

    z_diff = z2 - z1

    allocate(wght(negf%Np_n(1)))
    allocate(pnts(negf%Np_n(1)))

    call gauleg(0.d0,1.d0,pnts,wght,negf%Np_n(1))
    !call trapez(0.d0,1.d0,pnts,wght,negf%Np_n(1))

    do i = 1, negf%Np_n(1)
      Ec = z1 + pnts(i) * z_diff
      ff = fermi(Ec,muref,KbT)
      zt = negf%g_spin * z_diff * ff * wght(i) / (2.d0 *pi)

      en_grid(i)%path = 1
      en_grid(i)%pt = i
      en_grid(i)%pt_path = i
      en_grid(i)%Ec = Ec
      en_grid(i)%wght = zt
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
    !call trapez(0.d0,1.d0,pnts,wght,negf%Np_n(2))    !Setting weights for integration

    z1 = z2
    z2 = muref + Omega + j*Lambda

    z_diff = z2 - z1
    
    ioffs = negf%Np_n(1)
    
    do i = 1, negf%Np_n(2)
      Ec = z1 + pnts(i) * z_diff
      ff = fermi(Ec,muref,KbT)
      zt = negf%g_spin * z_diff * ff * wght(i) / (2.d0 *pi)

      en_grid(ioffs+i)%path = 2
      en_grid(ioffs+i)%pt = ioffs + i
      en_grid(ioffs+i)%pt_path = i
      en_grid(ioffs+i)%Ec = Ec
      en_grid(ioffs+i)%wght = zt
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
    ioffs = negf%Np_n(1)+negf%Np_n(2)
    
    do i = 1, Npoles
      Ec = muref + j * KbT *pi* (2.d0*i - 1.d0)
      zt= -j * KbT * negf%g_spin *(1.d0,0.d0)

      en_grid(ioffs+i)%path = 3
      en_grid(ioffs+i)%pt = ioffs + i
      en_grid(ioffs+i)%pt_path = i
      en_grid(ioffs+i)%Ec = Ec
      en_grid(ioffs+i)%wght = zt
    enddo

    ! *******************************************************************************          
    ! Distribution of Energy grid 
    ! pts 1 2 3 4 5 6 7 8 9 ...
    ! cpu 0 1 2 3 0 1 2 3 0 ...  
    ! *******************************************************************************          
    do i = 0, Ntot-1
       en_grid(i+1)%cpu = mod(i,numprocs) 
    enddo

  end subroutine contour_int_n_def



  !-----------------------------------------------------------------------
  ! Contour integration for density matrix, holes
  ! DOES INCLUDE FACTOR 2 FOR SPIN !!
  !-----------------------------------------------------------------------
  ! Performs the complex contur integration
  !
  !      T>=0
  !                +
  !                +
  !       * * * * * * * * * * *
  !       *        +
  !  --- -*---========================
  !      Elow Ec  muref
  !
  !-----------------------------------------------------------------------
  ! Contour integration for density matrix
  ! DOES INCLUDE FACTOR 2 FOR SPIN !!
  !-----------------------------------------------------------------------

  subroutine contour_int_p_def(negf)
    type(Tnegf) :: negf

    integer :: i, Ntot, Npoles, ioffs
    real(dp), DIMENSION(:), ALLOCATABLE :: wght,pnts   ! Gauss-quadrature points
    real(dp) :: Omega, Lambda
    real(dp) :: muref, kbT, Emax

    complex(dp) :: z1,z2,z_diff, zt
    complex(dp) :: Ec, ff

    kbT = negf%kbT(negf%refcont)
    muref = negf%muref
    Omega = negf%n_kt * kbT

    if (negf%n_poles.eq.0) then
      Lambda = 0.5d0* kbT * pi
    else
      Lambda = 2.d0* negf%n_poles * KbT * pi
    endif

    Emax = negf%Ev - negf%DeltaEv

    if ((Emax < (muref + 1.d-3)) .and. &
        (Emax > (muref - 1.d-3))) then
       Emax = muref + kbT
    endif

    Npoles = negf%n_poles
    if (Emax < muref) then
      Npoles = 0
    endif

    Ntot=negf%Np_p(1)+negf%Np_p(2)+Npoles
    allocate(en_grid(Ntot))

    ! *******************************************************************************
    ! 1. INTEGRATION OVER THE SEGMENT [Ec - dEc , Ec - dEc + j*Lambda]
    ! *******************************************************************************
    ! NEW INTEGRATION FOR COMPLEX DENSITY (T>0):
    !----------------------------------------------------
    !   g  [ /                ]      g   [ /                       ]
    !  --- [ |  Gr(z)*f(z) dz ] =   ---  [ | Gr(z)*(z2-z1)*f(z)*dt ]
    !  2pi [ /                ]     2*pi [ /                       ]
    !----------------------------------------------------

    z1 = Emax
    z2 = Emax + j*Lambda

    z_diff = z2 - z1

    allocate(wght(negf%Np_p(1)))
    allocate(pnts(negf%Np_p(1)))

    call gauleg(0.d0,1.d0,pnts,wght,negf%Np_p(1))

    do i = 1, negf%Np_p(1)
      Ec = z1 + pnts(i) * z_diff
      ff = fermi(-Ec,-muref,KbT)
      zt = negf%g_spin * z_diff * ff * wght(i) / (2.d0 *pi)

      en_grid(i)%path = 1
      en_grid(i)%pt = i
      en_grid(i)%pt_path = i
      en_grid(i)%Ec = Ec
      en_grid(i)%wght = zt
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

    allocate(wght(negf%Np_p(2)))
    allocate(pnts(negf%Np_p(2)))

    call gauleg(0.d0,1.d0,pnts,wght,negf%Np_p(2))    !Setting weights for integration

    z1 = z2 
    z2 = muref - Omega + j*Lambda

    z_diff = z2 - z1

    ioffs = negf%Np_p(1)

    do i = 1, negf%Np_p(2)
      Ec = z1 + pnts(i) * z_diff
      ff = fermi(-Ec,-muref,KbT)
      zt = negf%g_spin * z_diff * ff * wght(i) / (2.d0 *pi)

      en_grid(ioffs+i)%path = 2
      en_grid(ioffs+i)%pt = ioffs + i
      en_grid(ioffs+i)%pt_path = i
      en_grid(ioffs+i)%Ec = Ec
      en_grid(ioffs+i)%wght = zt
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
    ioffs = negf%Np_p(1)+negf%Np_p(2)

    do i = 1, Npoles
      Ec = muref + j * KbT *pi* (2.d0*i - 1.d0)
      zt= -j * KbT * negf%g_spin *(1.d0,0.d0)

      en_grid(ioffs+i)%path = 3
      en_grid(ioffs+i)%pt = ioffs + i
      en_grid(ioffs+i)%pt_path = i
      en_grid(ioffs+i)%Ec = Ec
      en_grid(ioffs+i)%wght = zt
    enddo

    ! *******************************************************************************
    ! Distribution of Energy grid
    ! pts 1 2 3 4 5 6 7 8 9 ...
    ! cpu 0 1 2 3 0 1 2 3 0 ...
    ! *******************************************************************************
    do i = 0, Ntot-1
       en_grid(i+1)%cpu = mod(i,numprocs)
    enddo

  end subroutine contour_int_p_def

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
  subroutine contour_int_def(negf)
     type(Tnegf) :: negf

     integer :: i, Ntot, ioffs
     real(dp) :: Lambda, Rad, Centre
     real(dp), DIMENSION(:), ALLOCATABLE :: wght,pnts   ! Gauss-quadrature points
     real(dp) :: dt, Elow
     real(dp) :: muref, mumin, kbT, nkT, alpha
     complex(dp) :: z1,z2,z_diff, zt
     complex(dp) :: Ec, ff, Pc
     
     kbT = negf%kbT(negf%refcont)
     muref = negf%muref
     nkT = negf%n_kt * kbT
     Lambda = 2.d0* negf%n_poles * KbT * pi
     mumin = muref - nkT
     
     Ntot=negf%Np_n(1)+negf%Np_n(2)+negf%n_poles
     allocate(en_grid(Ntot))
 
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
     if (kbT.ne.0.d0) then        
        alpha = atan(Lambda/(mumin-Centre)) 
     else
        alpha = 0.1d0*pi             
     end if
    
     !Setting weights for gaussian integration 
     allocate(wght(negf%Np_n(1)))
     allocate(pnts(negf%Np_n(1)))
    
     call gauleg(pi,alpha,pnts,wght,negf%Np_n(1))
     do i = 1, negf%Np_n(1)
        Pc = Rad*exp(j*pnts(i))
        Ec = Centre+Pc
        zt = j * Pc * negf%g_spin * wght(i)/(2.d0*pi)
        en_grid(i)%path=1
        en_grid(i)%pt_path=i
        en_grid(i)%pt=i
        en_grid(i)%Ec=Ec
        en_grid(i)%wght=zt
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
    
     if (kbT.eq.0.d0) then                        ! Circle integration T=0
       call  gauleg(alpha,0.d0,pnts,wght,negf%Np_n(2))
     else                                          ! Segment integration T>0
       z1 = muref + nkT + j*Lambda
       z2 = muref - nkT + j*Lambda
       z_diff = z2 - z1
       call  gauleg(1.d0,0.d0,pnts,wght,negf%Np_n(2))    !Setting weights for integration
     endif

     ioffs = negf%Np_n(1)
     
     do i = 1, negf%Np_n(2)
        if (kbT.eq.0.d0) then                  ! Circle integration T=0            
           Pc = Rad*exp(j*pnts(i))
           Ec = Centre+Pc
           dt = negf%g_spin*wght(i)/(2.d0*pi)
           zt = dt*Pc*j
        else                                        ! Segment integration T>0
           Ec = z1 + pnts(i)*z_diff
           ff = fermi(Ec,muref,KbT)
           zt = negf%g_spin * z_diff * ff * wght(i) / (2.d0 *pi)
        endif
        en_grid(ioffs+i)%path=2
        en_grid(ioffs+i)%pt_path=i
        en_grid(ioffs+i)%pt=ioffs+i
        en_grid(ioffs+i)%Ec=Ec
        en_grid(ioffs+i)%wght=zt
     end do
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
     ioffs = negf%Np_n(1)+negf%Np_n(2)
     do i = 1, negf%n_poles
        Ec = muref + j * KbT *pi* (2.d0*real(i,dp) - 1.d0)   
        zt= -j*negf%g_spin*KbT
        en_grid(ioffs+i)%path=3
        en_grid(ioffs+i)%pt_path=i
        en_grid(ioffs+i)%pt=ioffs+i
        en_grid(ioffs+i)%Ec=Ec
        en_grid(ioffs+i)%wght=zt
     enddo

     ! *******************************************************************************          
     ! Distribution of Energy grid 
     ! pts 1 2 3 4 5 6 7 8 9 ...
     ! cpu 0 1 2 3 0 1 2 3 0 ...  
     ! *******************************************************************************          
     do i = 0, Ntot-1
        en_grid(i+1)%cpu = mod(i,numprocs) 
     enddo
     
  end subroutine contour_int_def
  !-----------------------------------------------------------------------

  subroutine contour_int(negf)
     type(Tnegf) :: negf

     Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
     type(z_CSR) :: GreenR, TmpMt 
     integer :: i, i1, ncont, Ntot, outer
     real(dp) :: ncyc
     complex(dp) :: Ec, zt
    
     ncont = negf%str%num_conts
     outer = negf%outer 
     Ntot = size(en_grid)  
     call create(TmpMt,negf%H%nrow,negf%H%ncol,negf%H%nrow)
     call initialize(TmpMt)
    
     call write_info(negf%verbose,'CONTOUR INTEGRAL',Ntot)
     
     do i = 1, Ntot
   
        if (en_grid(i)%cpu .ne. id) cycle

        call write_point(negf%verbose,en_grid(i), Ntot) 
    
        if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Green`s funct ')
    
        Ec = en_grid(i)%Ec 
        zt = en_grid(i)%wght
        negf%iE = en_grid(i)%pt 
        
        call compute_contacts(Ec,negf,ncyc,Tlc,Tcl,SelfEneR,GS)
  
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
 
     ! The MPI gather here

     call destroy(TmpMt)

     call destroy_en_grid()

  end subroutine contour_int

  !--------------------------------------------!
  !--------------------------------------------!
  ! Non equilibrium integration over real axis !
  !--------------------------------------------!
  !-----------------------------------------------------------------------
  !    g    --  [ /                                         ]
  ! ------  >   [ | j [f(i) - f(ref)] Gr(E) Gam_i Ga(E) dE  ]  
  ! 2*pi*j  --i [ /                                         ]
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  ! Contour integration for density matrix 
  ! DOES INCLUDE FACTOR 2 FOR SPIN !! 
  !-----------------------------------------------------------------------
  subroutine real_axis_int_def(negf)
    type(Tnegf) :: negf

    integer :: i, ioffset, ncont, Ntot
    real(dp), DIMENSION(:), allocatable :: wght,pnts   ! Gauss-quadrature points
    real(dp) :: Omega, mumin, mumax
    
    ncont = negf%str%num_conts
    ioffset = negf%Np_n(1) + negf%Np_n(2) + negf%n_poles
    ! Omega considers maximum kT so interval is always large enough
    Omega = negf%n_kt * maxval(negf%kbT) 
    
    if (ncont.gt.0) then 
       mumin=minval(negf%mu(1:ncont))
       mumax=maxval(negf%mu(1:ncont))
    else
       mumin=negf%mu(1)
       mumax=negf%mu(1)
    endif  

    Ntot = negf%Np_real(1)
    allocate(en_grid(Ntot))

    allocate(pnts(Ntot))
    allocate(wght(Ntot))
     
    call gauleg(mumin-Omega,mumax+Omega,pnts,wght,Ntot)
    
    do i = 1, Ntot
       en_grid(i)%path = 1
       en_grid(i)%pt = ioffset + i
       en_grid(i)%pt_path = i
       en_grid(i)%Ec = cmplx(pnts(i),negf%delta,dp)
       en_grid(i)%wght = negf%wght * negf%g_spin * wght(i)/(2.d0 *pi)
    enddo

    deallocate(wght)
    deallocate(pnts)

    ! distribute energy grid
    do i = 0, Ntot-1
       en_grid(i+1)%cpu = mod(i,numprocs)
    enddo

  end subroutine real_axis_int_def
  !-----------------------------------------------------------------------
  !    g    --  [ /                                         ]
  ! ------  >   [ | j [f(i) - f(ref)] Gr(E) Gam_i Ga(E) dE  ]  
  ! 2*pi*j  --i [ /                                         ]
  !-----------------------------------------------------------------------
  ! note: refcont is passed to calls_neq via TNegf
  !-----------------------------------------------------------------------

  subroutine real_axis_int(negf)
    type(Tnegf) :: negf

    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    type(z_CSR) :: GreenR, TmpMt 

    integer :: ref, Npoints
    integer :: i, i1, j1, outer, ncont

    real(dp), DIMENSION(:), allocatable :: frm_f
    real(dp) :: ncyc, Er

    complex(dp) :: zt
    complex(dp) :: Ec

    ncont = negf%str%num_conts
    outer = negf%outer
    Npoints = size(en_grid)

    call log_allocate(frm_f,ncont)
    
    call create(TmpMt,negf%H%nrow,negf%H%ncol,negf%H%nrow)
    call initialize(TmpMt)
    
    call write_info(negf%verbose,'REAL AXIS INTEGRAL',Npoints)

    do i = 1, Npoints

       if (en_grid(i)%cpu .ne. id) cycle

       Ec = en_grid(i)%Ec
       Er = real(Ec)
       zt = en_grid(i)%wght
       negf%iE = en_grid(i)%pt

       call write_point(negf%verbose,en_grid(i),Npoints)

       do j1 = 1,ncont
          frm_f(j1)=fermi(Er,negf%mu(j1),negf%kbT(j1))
       enddo

       if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Green`s funct ')

       call compute_contacts(Ec,negf,ncyc,Tlc,Tcl,SelfEneR,GS)

       call calls_neq_mem_dns(negf,Er,SelfEneR,Tlc,Tcl,GS,negf%str,frm_f,GreenR,outer)

       do i1=1,ncont
          call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
       enddo

       if(negf%DorE.eq.'D') then
          call concat(TmpMt,zt,GreenR,1,1)
       endif
       if(negf%DorE.eq.'E') then
          call concat(TmpMt,zt*Er,GreenR,1,1)
       endif

       call destroy(GreenR) 

       if (id0.and.negf%verbose.gt.VBT) call write_clock

    enddo

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
     
    call destroy_en_grid()

  end subroutine real_axis_int
  
  !--------------------------------------------!
  !--------------------------------------------!
  ! Equilibrium integration over real axis     !
  !--------------------------------------------!
  !-----------------------------------------------------------------------
  !    g    --  [ /                                         ]
  ! ------  >   [ | j [f(E) Gr(E) dE ]                      ]  
  ! 2*pi*j  --i [ /                                         ]
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  ! Contour integration for density matrix 
  ! DOES INCLUDE FACTOR 2 FOR SPIN !! 
  !-----------------------------------------------------------------------
  subroutine real_axis_int_n_def(negf)
    type(Tnegf) :: negf

    integer :: i, ioffset, ncont, Ntot
    real(dp), DIMENSION(:), allocatable :: wght,pnts   ! Gauss-quadrature points
    real(dp) :: Omega, mumin, mumax, muref, ff, kbT
    
    ncont = negf%str%num_conts
    ioffset = negf%Np_n(1) + negf%Np_n(2) + negf%n_poles
    kbT = maxval(negf%kbT)
    Omega = negf%n_kt * kbT 
    muref = negf%muref
    
    if (ncont.gt.0) then 
       mumin=minval(negf%mu_n(1:ncont))
       mumax=maxval(negf%mu_n(1:ncont))
    else
       mumin=negf%mu_n(1)
       mumax=negf%mu_n(1)
    endif  

    Ntot = negf%Np_real(1)
    allocate(en_grid(Ntot))

    allocate(pnts(Ntot))
    allocate(wght(Ntot))
    
    call gauleg(mumin+negf%Ec+negf%DeltaEc, mumax+Omega,pnts,wght,Ntot)
    
    do i = 1, Ntot
       en_grid(i)%path = 1
       en_grid(i)%pt = ioffset + i
       en_grid(i)%pt_path = i
       en_grid(i)%Ec = cmplx(pnts(i),negf%delta,dp)
       ff = fermi(pnts(i),muref,KbT)
       en_grid(i)%wght = negf%g_spin * negf%wght * ff * wght(i) / (2.d0 * pi)
    enddo

    deallocate(wght)
    deallocate(pnts)

    ! distribute energy grid
    do i = 0, Ntot-1
       en_grid(i+1)%cpu = mod(i,numprocs)
    enddo

  end subroutine real_axis_int_n_def

  !--------------------------------------------!
  !--------------------------------------------!
  ! Equilibrium integration over real axis     !
  !--------------------------------------------!
  !-----------------------------------------------------------------------
  !    g    --  [ /                                         ]
  ! ------  >   [ | j [1-f(E) Gr(E) dE ]                    ]  
  ! 2*pi*j  --i [ /                                         ]
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  ! Contour integration for density matrix 
  ! DOES INCLUDE FACTOR 2 FOR SPIN !! 
  !-----------------------------------------------------------------------
  subroutine real_axis_int_p_def(negf)
    type(Tnegf) :: negf

    integer :: i, ioffset, ncont, Ntot
    real(dp), DIMENSION(:), allocatable :: wght,pnts   ! Gauss-quadrature points
    real(dp) :: Omega, mumin, mumax, muref, ff, kbT
    
    ncont = negf%str%num_conts
    ioffset = negf%Np_n(1) + negf%Np_n(2) + negf%n_poles
    kbT = maxval(negf%kbT) 
    Omega = negf%n_kt * kbT
    muref = negf%muref
    
    if (ncont.gt.0) then 
       mumin=minval(negf%mu_p(1:ncont))
       mumax=maxval(negf%mu_p(1:ncont))
    else
       mumin=negf%mu_p(1)
       mumax=negf%mu_p(1)
    endif  

    Ntot = negf%Np_real(1)
    allocate(en_grid(Ntot))

    allocate(pnts(Ntot))
    allocate(wght(Ntot))
    
    mumax = negf%Ev+negf%DeltaEv
    
    call gauleg(mumin-Omega, mumax, pnts, wght, Ntot)
    
    do i = 1, Ntot
       en_grid(i)%path = 1
       en_grid(i)%pt = ioffset + i
       en_grid(i)%pt_path = i
       en_grid(i)%Ec = cmplx(pnts(i),negf%delta,dp)
       ff = fermi(-pnts(i),-muref,KbT)
       en_grid(i)%wght = negf%g_spin * negf%wght * ff * wght(i) / (2.d0 * pi)
    enddo

    deallocate(wght)
    deallocate(pnts)

    ! distribute energy grid
    do i = 0, Ntot-1
       en_grid(i+1)%cpu = mod(i,numprocs)
    enddo

  end subroutine real_axis_int_p_def

  !-----------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! Gauss - Legendre quadrature weights and points
  !----------------------------------------------------------------------------
  subroutine gauleg(x1,x2,x,w,n)

    real(kind=dp), PARAMETER :: ACC = 1d-15

    INTEGER n
    real(kind=dp) :: x1,x2,x(n),w(n)

    INTEGER i,k,m
    real(kind=dp) :: p0,p1,p2,pp,xl,xm,z,z1

    m=(n+1)/2

    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)

    do i=1,m

       ! Approssimazione degli zeri dei polinomi di Legendre:
       z=cos(Pi*(i-0.25d0)/(n+0.5d0))

       ! Legendre polynomial, p1, evaluated by rec. relations:
       ! P(0)=1; P(-1)=0
       ! P(n) = (2n-1)*x*P(n-1) - (n-1)*P(n-2)
       !
       ! Derivative pp using the relation of p1 and p2:
       ! P'(n) = (2n-1)*P(n-1) + (2n-1)*x*P'(n-1) - (n-1)*P'(n-2)
       !
       ! Newton method is used to refine the zeros
       !
       do
          p0=1.d0  !p(0)
          p1=0.d0  !p(-1)

          ! Legendre polynomial p1 evaluated by rec. relations:
          do k=1,n
             p2=p1 !p(-2)=p(-1)
             p1=p0 !p(-1)=p(0)
             p0=((2.d0*k-1.d0)*z*p1-(k-1.d0)*p2)/k 
          enddo
          
          pp=n*(z*p0-p1)/(z*z-1.d0)

          ! Newton method to refine the zeros:
          z1=z
          z=z1-p0/pp

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

    d = (x2-x1)/(1.0_dp*(n-1))

    w = d * 1.0_dp
    w(1) = d * 0.5_dp
    w(n) = d * 0.5_dp

    do i = 1, n
       x(i) = ( x1*(n-i) + x2*(i-1) ) / (n-1) 
    enddo
 
  end subroutine trapez
  !--------------------------------------------------------------------  
  subroutine simpsons(x1,x2,x,w,n)
    real(kind=dp), PARAMETER :: ACC = 1d-15

    INTEGER n
    real(kind=dp) :: x1,x2,x(n),w(n)

    real(dp) :: d
    integer :: i

    if (mod(n-1,2).ne.0) STOP 'ERROR: N is not multiple of 2'

    d = (x2-x1)/((n-1)*1.0_dp)

    w = d * 4.0_dp/3.0_dp
    w(1) = d * 1.0_dp/3.0_dp
    w(n) = d * 1.0_dp/3.0_dp
    do i = 3,n-1,2
      w(i) = d* 2.0_dp/3.0_dp
    enddo

    do i = 1, n
       x(i) = ( x1*(n-i) + x2*(i-1) ) / (n-1) 
    enddo


  end subroutine simpsons
  !--------------------------------------------------------------------  
  subroutine three_eigth(x1,x2,x,w,n)

    real(kind=dp), PARAMETER :: ACC = 1d-15

    INTEGER n
    real(kind=dp) :: x1,x2,x(n),w(n)

    real(dp) :: d
    integer :: i

    if (mod(n-1,3).ne.0) STOP 'ERROR: N-1 is not multiple of 3'

    d = (x2-x1)/((n-1)*1.0_dp)

    w = d * 9.0_dp/8.0_dp
    w(1) = d * 3.0_dp/8.0_dp
    w(n) = d * 3.0_dp/8.0_dp
    do i = 4,n-1,3
      w(i) = d* 6.0_dp/8.0_dp
    enddo

    do i = 1, n
       x(i) = ( x1*(n-i) + x2*(i-1) ) / (n-1) 
    enddo
 
  end subroutine three_eigth
  !--------------------------------------------------------------------  

  !-----------------------------------------------------------------------
  !  --  [ /                           ]
  !  >   [ |  [f(i) - f(ref)] T(E) dE  ]  
  !  --i [ /                           ]
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  ! Contour integration for density matrix 
  ! DOES INCLUDE FACTOR 2 FOR SPIN !! 
  !-----------------------------------------------------------------------
  subroutine tunneling_int_def(negf)
    type(Tnegf) :: negf

    integer :: i, ncont, Nsteps
    
    Nsteps=NINT((negf%Emax-negf%Emin)/negf%Estep) + 1
    allocate(en_grid(Nsteps))

    do i = 1, Nsteps
       en_grid(i)%path = 1
       en_grid(i)%pt = i
       en_grid(i)%pt_path = i
       en_grid(i)%Ec = cmplx(negf%Emin + negf%Estep*(i-1), 0.0, dp) 
       en_grid(i)%wght = negf%wght * negf%g_spin 
    enddo

    ! distribute energy grid
    do i = 0, Nsteps-1
       en_grid(i+1)%cpu = mod(i,numprocs)
    enddo

  end subroutine tunneling_int_def

  !-----------------------------------------------------------------------
  !  Routine to compute T(E) and (optionally) LDOS(E)
  !  PDOS is computed if negf%nLDOS > 0
  !  When only T(E) is needed, a fast algorithm is used (reduction to one block)
  !------------------------------------------------------------------------------- 
  subroutine tunneling_and_dos(negf)
    type(Tnegf) :: negf

    ! Local Variables
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    Real(dp), Dimension(:), allocatable :: TUN_MAT
    Real(dp), Dimension(:), allocatable :: LEDOS
    Real(dp) :: mu1, mu2   ! contact potentials
    Real(dp) :: ncyc       ! stores average number of iters in decimation 

    Integer :: i, icont, icpl      ! dummy counters
    Integer :: ncont               ! number of contacts
    Integer :: size_ni, size_nf    ! emitter-collector contacts
    
    Integer :: Nstep               ! number of integration points
    Complex(dp) :: Ec              ! Energy point
    
    Logical :: do_LEDOS            ! performs or not LDOS
    
    ! Get out immediately if Emax<Emin 
    if (negf%Emax.le.negf%Emin) then
       if(id0) write(*,*) '0 tunneling points;  current = 0.0'
       !call log_allocatep(negf%tunn_mat,0,0)
       !if (do_ledos) call log_allocatep(negf%ldos_mat,0,0)
       call destroy_en_grid()
       call log_allocatep(negf%currents,1) 
       negf%currents = 0.0_dp 
       return
    endif
    !-------------------------------------------------------
    
    do_LEDOS = .false. 
    if(negf%nLDOS.gt.0) do_LEDOS=.true.
    ncont = negf%str%num_conts
    Nstep = size(en_grid)
    ncyc=0
    
    !Extract emitter-collector contacts -------------------
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
    !-------------------------------------------------------
    
    call log_allocate(TUN_MAT,size_ni)
    call log_allocatep(negf%tunn_mat,Nstep,size_ni)   
    negf%tunn_mat = 0.0_dp 

    if (do_LEDOS) then
       call log_allocatep(negf%ldos_mat,Nstep,negf%nLDOS)
       call log_allocate(LEDOS,negf%nLDOS)          
       negf%ldos_mat(:,:)=0.d0
    endif
    !-------------------------------------------------------
    
    
    !Loop on energy points: tunneling 
    do i = 1, Nstep
      
       if (en_grid(i)%cpu /= id) cycle
      
       Ec = en_grid(i)%Ec
       negf%iE = en_grid(i)%pt

       call write_point(negf%verbose,en_grid(i), size(en_grid))

       if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Contact SE ')       
       call compute_contacts(Ec+(0.d0,1.d0)*negf%delta,negf,ncyc,Tlc,Tcl,SelfEneR,GS)
       if (id0.and.negf%verbose.gt.VBT) call write_clock
      

       if (.not.do_LEDOS) then
          if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Tunneling ') 

          call tunneling_dns(negf%H,negf%S,Ec,SelfEneR,negf%ni,negf%nf,size_ni, &
                             & negf%str,TUN_MAT)

          negf%tunn_mat(i,:) = TUN_MAT(:) * negf%wght
       else
          if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Tunneling and DOS') 
          LEDOS(:) = 0.d0
          
          call tun_and_dos(negf%H,negf%S,Ec,SelfEneR,GS,negf%ni,negf%nf,negf%nLDOS, &
                           & negf%LDOS,size_ni,negf%str,TUN_MAT,LEDOS)
          
          negf%tunn_mat(i,:) = TUN_MAT(:) * negf%wght
          negf%ldos_mat(i,:) = LEDOS(:) * negf%wght
       endif

       if (id0.and.negf%verbose.gt.VBT) call write_clock
       
       do icont=1,ncont
          call destroy(Tlc(icont))
          call destroy(Tcl(icont))
          call destroy(SelfEneR(icont))
          call destroy(GS(icont))
       enddo
       
    enddo !Loop on energy 

    call destroy_en_grid()
    call log_deallocate(TUN_MAT)
    if(do_LEDOS) call log_deallocate(LEDOS)
  
  end subroutine tunneling_and_dos

  !---------------------------------------------------------------------------
  !   COMPUTATION OF CURRENTS 
  !---------------------------------------------------------------------------
  subroutine electron_current(negf)
    type(Tnegf) :: negf
    
    integer :: size_ni, ii
    real(dp) :: mu1, mu2

    size_ni = size(negf%tunn_mat,2)

    if (.not.associated(negf%currents)) then
      call log_allocatep(negf%currents,size_ni)
    end if
    negf%currents=0.d0

    do ii=1,size_ni
       mu1=negf%mu(negf%ni(ii))
       mu2=negf%mu(negf%nf(ii))
       
       negf%currents(ii)= integrate_el(negf%tunn_mat(:,ii), mu1, mu2, &
                            & negf%kbT(negf%ni(ii)), negf%kbT(negf%nf(ii)), &
                            & negf%Emin, negf%Emax, negf%Estep, negf%g_spin)
    enddo

  end subroutine electron_current

  !-----------------------------------------------------------------------
  !  Routine to compute T(E) and (optionally) LDOS(E)
  !  PDOS is computed if negf%nLDOS > 0
  !  When only T(E) is needed, a fast algorithm is used (reduction to one block)
  !------------------------------------------------------------------------------- 
  subroutine phonon_tunneling(negf)
    type(Tnegf) :: negf

    ! Local Variables
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    Real(dp), Dimension(:), allocatable :: TUN_MAT
    Real(dp), Dimension(:), allocatable :: LEDOS
    Real(dp) :: mu1, mu2   ! contact potentials
    Real(dp) :: ncyc       ! stores average number of iters in decimation 

    Integer :: i, icont, icpl      ! dummy counters
    Integer :: ncont               ! number of contacts
    Integer :: size_ni, size_nf    ! emitter-collector contacts
    
    Integer :: Nstep               ! number of integration points
    Complex(dp) :: Ec              ! Energy point
    Complex(dp) :: delta
    Logical :: do_LEDOS            ! performs or not LDOS
    
    ! Get out immediately if Emax<Emin 
    if (negf%Emax.le.negf%Emin) then
       if(id0) write(*,*) '0 tunneling points;  current = 0.0'
       !call log_allocatep(negf%tunn_mat,0,0)
       !if (do_ledos) call log_allocatep(negf%ldos_mat,0,0)
       call destroy_en_grid()
       call log_allocatep(negf%currents,1) 
       negf%currents = 0.0_dp 
       return
    endif
    !-------------------------------------------------------
    
    do_LEDOS = .false. 
    if(negf%nLDOS.gt.0) do_LEDOS=.true.
    ncont = negf%str%num_conts
    Nstep = size(en_grid)
    ncyc=0
    
    !Extract emitter-collector contacts -------------------
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
    !-------------------------------------------------------
    
    call log_allocate(TUN_MAT,size_ni)
    call log_allocatep(negf%tunn_mat,Nstep,size_ni)   
    negf%tunn_mat = 0.0_dp 

    if (do_LEDOS) then
       call log_allocatep(negf%ldos_mat,Nstep,negf%nLDOS)
       call log_allocate(LEDOS,negf%nLDOS)          
       negf%ldos_mat(:,:)=0.d0
    endif
    !-------------------------------------------------------
    
    
    !Loop on energy points: tunneling 
    do i = 1, Nstep
      
       if (en_grid(i)%cpu /= id) cycle
      
       Ec = en_grid(i)%Ec * en_grid(i)%Ec
       negf%iE = en_grid(i)%pt
       !delta = negf%delta * negf%delta 
       delta = negf%delta * (1.0_dp - real(en_grid(i)%Ec)/(negf%Emax+1d-12)) * Ec 

       call write_point(negf%verbose,en_grid(i), size(en_grid))

       if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Contact SE ')       
       call compute_contacts(Ec+(0.d0,1.d0)*delta,negf,ncyc,Tlc,Tcl,SelfEneR,GS)
       if (id0.and.negf%verbose.gt.VBT) call write_clock
      

       if (.not.do_LEDOS) then
          if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Tunneling ') 

          call tunneling_dns(negf%H,negf%S,Ec,SelfEneR,negf%ni,negf%nf,size_ni, &
                             & negf%str,TUN_MAT)

          negf%tunn_mat(i,:) = TUN_MAT(:) * negf%wght
       else
          if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Tunneling and DOS') 
          LEDOS(:) = 0.d0
          
          call tun_and_dos(negf%H,negf%S,Ec,SelfEneR,GS,negf%ni,negf%nf,negf%nLDOS, &
                           & negf%LDOS,size_ni,negf%str,TUN_MAT,LEDOS)
          
          negf%tunn_mat(i,:) = TUN_MAT(:) * negf%wght
          negf%ldos_mat(i,:) = LEDOS(:) * negf%wght
       endif

       if (id0.and.negf%verbose.gt.VBT) call write_clock
       
       do icont=1,ncont
          call destroy(Tlc(icont))
          call destroy(Tcl(icont))
          call destroy(SelfEneR(icont))
          call destroy(GS(icont))
       enddo
       
    enddo !Loop on energy 

    call destroy_en_grid()
    call log_deallocate(TUN_MAT)
    if(do_LEDOS) call log_deallocate(LEDOS)
  
  end subroutine phonon_tunneling
  
  !---------------------------------------------------------------------------
  subroutine phonon_current(negf)
    type(Tnegf) :: negf
    
    integer :: size_ni, ii

    size_ni = size(negf%tunn_mat,2)

    call log_allocatep(negf%currents,size_ni)
    negf%currents=0.d0

    do ii=1,size_ni
       
       negf%currents(ii)= integrate_ph(negf%tunn_mat(:,ii),  &
                            & negf%kbT(negf%ni(ii)), negf%kbT(negf%nf(ii)), &
                            & negf%Emin, negf%Emax, negf%Estep)
    enddo
        
    
  end subroutine phonon_current
 
  !////////////////////////////////////////////////////////////////////////
  !************************************************************************
  ! Function to integrate the tunneling and get the current
  ! The function resolves fermi(E) on a fine grid interpolating linearly T(E)  
  ! In this way a more precise integration is obtained when T ~ constant
  !************************************************************************
  function integrate_el(TUN_TOT,mu1,mu2,kT1,kT2,emin,emax,estep,spin_g)

    implicit none
 
    real(dp) :: integrate_el
    real(dp), intent(in) :: mu1,mu2,emin,emax,estep
    real(dp), dimension(:), intent(in) :: TUN_TOT
    real(dp), intent(in) :: kT1, kT2 
    real(dp), intent(in) :: spin_g 
 
    REAL(dp) :: destep,kbT1,kbT2,TT1,TT2,E3,E4,TT3,TT4
    REAL(dp) :: E1,E2,c1,c2,curr
    INTEGER :: i,i1,N,Nstep,imin,imax
 
    curr=0.d0
    N=0
    destep=1.0d10 
    Nstep=NINT((emax-emin)/estep);
 
    if (kT1.lt.0.01_dp*Kb) then
      kbT1 = Kb*0.01_dp
    else
      kbT1 = kT1     
    endif
    if (kT2.lt.0.01_dp*Kb) then
      kbT2 = Kb*0.01_dp
    else
      kbT2 = kT2     
    endif
 
    if (mu2<mu1) then
      call swap(mu1,mu2)
      call swap(kbT1,kbT2)
    end if
 
    imin=0
    imax=Nstep 
 
    ! performs the integration with simple trapezium rule. 
    do i=imin,imax-1
       
       E1=emin+estep*i  
       TT1=TUN_TOT(i+1)     
       E2=emin+estep*(i+1)
       TT2=TUN_TOT(i+2)
       
       ! Each step is devided into substeps in order to
       ! smooth out the Fermi function
       do while (destep.ge.(kbT1+kbT2)/20.0_dp) 
          N=N+1
          destep=(E2-E1)/N
       enddo
       
       ! Within each substep the tunneling is linearly interpolated
       ! Possibly perform a cubic-spline interpolation in future 
       do i1=0,N-1
          
          E3=E1+(E2-E1)*i1/N
          E4=E3+(E2-E1)/N
          TT3=( TT2-TT1 )*i1/N + TT1
          TT4=TT3 + (TT2-TT1)/N
          
          c1=(fermi(E3,mu2,KbT2)-fermi(E3,mu1,KbT1))*TT3
          c2=(fermi(E4,mu2,KbT2)-fermi(E4,mu1,KbT1))*TT4
          
          curr=curr+spin_g*(c1+c2)*(E4-E3)/2.d0
          
       enddo
       
    enddo
 
    integrate_el = curr
  
  end function integrate_el
 
  !////////////////////////////////////////////////////////////////////////
  !************************************************************************
  ! Function to integrate the tunneling and get the current
  ! The function resolves fermi(E) on a fine grid interpolating linearly T(E)  
  ! In this way a more precise integration is obtained when T ~ constant
  !************************************************************************
  function integrate_ph(TUN_TOT,kT1,kT2,emin,emax,estep)

    implicit none
 
    real(dp) :: integrate_ph
    real(dp), intent(in) :: emin,emax,estep
    real(dp), dimension(:), intent(in) :: TUN_TOT
    real(dp), intent(in) :: kT1, kT2 
 
    REAL(dp) :: destep,kbT1,kbT2,TT1,TT2,E3,E4,TT3,TT4
    REAL(dp) :: E1,E2,c1,c2,curr
    INTEGER :: i,i1,N,Nstep,imin,imax
 
    curr=0.d0
    N=0
    destep=1.0d10 
    Nstep=NINT((emax-emin)/estep);
 
    if (kT1.lt.0.01_dp*Kb) then
      kbT1 = Kb*0.01_dp
    else
      kbT1 = kT1     
    endif
    if (kT2.lt.0.01_dp*Kb) then
      kbT2 = Kb*0.01_dp
    else
      kbT2 = kT2     
    endif
 
    imin=0
    imax=Nstep 
 
    ! performs the integration with simple trapezium rule. 
    do i=imin,imax-1
       
       E1=emin+estep*i  
       TT1=TUN_TOT(i+1)     
       E2=emin+estep*(i+1)
       TT2=TUN_TOT(i+2)
       
       ! Each step is devided into substeps in order to
       ! smooth out the Fermi function
       do while (destep.ge.(kbT1+kbT2)/20.0_dp) 
          N=N+1
          destep=(E2-E1)/N
       enddo
       
       ! Within each substep the tunneling is linearly interpolated
       ! Possibly perform a cubic-spline interpolation in future 
       do i1=0,N-1
          
          E3=E1+(E2-E1)*i1/N
          E4=E3+(E2-E1)/N
          TT3=( TT2-TT1 )*i1/N + TT1
          TT4=TT3 + (TT2-TT1)/N
          
          c1=(bose(E3,KbT2)-bose(E3,KbT1))*TT3
          c2=(bose(E4,KbT2)-bose(E4,KbT1))*TT4
          
          curr=curr+(c1+c2)*(E4-E3)*(E4-E3)/2.d0
          
       enddo
       
    enddo
 
    integrate_ph = curr  
  
  end function integrate_ph

  !/////////////////////////////////////////////////////////////////////////
  function thermal_conductance(TUN_TOT,kbT,emin,emax,estep)
    implicit none

    real(dp) :: thermal_conductance
    real(dp), intent(in) :: emin,emax,estep
    real(dp), dimension(:), intent(in) :: TUN_TOT
    real(dp), intent(in) :: kbT  ! temperature

    REAL(dp) :: destep,TT1,TT2
    REAL(dp) :: E1,E2,c1,c2,curr
    INTEGER :: i,i1,N,Nstep,imin,imax

    curr=0.d0
    Nstep=NINT((emax-emin)/estep);

    TT1=TUN_TOT(1)
    do i = 0, 9
      E1=emin*i/10
      E2=emin*(i+1)/10           
      c1=diff_bose(E1,kbT)*TT1
      c2=diff_bose(E2,kbT)*TT1
      curr=curr+(c1+c2)*emin/20.d0 
    end do

    ! performs the integration with simple trapezium rule. 
!    do i=1,100

!       TT1=10.d0+3.0*(i-1)

       ! Within each substep the tunneling is linearly interpolated
       ! Possibly perform a cubic-spline interpolation in future 
       do i=0,Nstep-1

         E1=emin+estep*i
         TT1=TUN_TOT(i+1)
         E2=emin+estep*(i+1)           
         TT2=TUN_TOT(i+2)

         c1=diff_bose(E1,kbT)*TT1
         c2=diff_bose(E2,kbT)*TT2

         curr=curr+(c1+c2)*estep/2.d0 
       enddo
!    enddo

    thermal_conductance = curr

  end function thermal_conductance

  !////////////////////////////////////////////////////////////////////////
  subroutine partial_charge(negf,DensMat,qmulli,qtot)
    type(TNegf) :: negf
    type(z_CSR) :: DensMat
    real(dp), dimension(:) :: qmulli
    real(dp) :: qtot
    
    integer nrow,ii,jj,jcol,ka,kb
    real(dp) :: dd

    nrow = negf%S%nrow 
    qmulli = 0.d0
    
    ! Partial Sum_j[G_ij S_ji]
    do ii=1, nrow 
      do ka=DensMat%rowpnt(ii), DensMat%rowpnt(ii+1)-1 
        dd = DensMat%nzval(ka)
        jj = DensMat%colind(ka)
        
        do kb=negf%S%rowpnt(jj),negf%S%rowpnt(jj+1)-1
          jcol = negf%S%colind(kb)
          if (jcol .eq. ii) then
            qmulli(jcol) = qmulli(jcol) + real(dd*negf%S%nzval(kb))
          endif
        enddo
      enddo
    enddo
    
    !.............................................................
    ! Calculation of total charge 
    
    qtot = 0.d0
    do ii = 1,nrow
      qtot = qtot+qmulli(ii)
    enddo
    
  end subroutine partial_charge 
  
  subroutine swap(x1,x2)
    real(dp) :: x1,x2

    real(dp) :: tmp

    x1=tmp
    x1=x2
    x2=tmp

  end subroutine swap

end module integrations 
