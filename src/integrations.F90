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
 !use iterative_ph
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
 public :: quasiEq_int_p     ! real-axis integration - quasi-equilibrium - holes     
 public :: quasiEq_int_n     ! real-axis integration - quasi-equilibrium - electrons    
 public :: ldos_int          ! ldos only integrator

 public :: contour_int_def   ! contour integration for DFT(B)
 public :: contour_int_n_def ! contour integration for CB
 public :: contour_int_p_def ! contour integration for VB
 public :: real_axis_int_def ! real axis integration
 public :: real_axis_int_n_def ! integration of CB on real axis
 public :: real_axis_int_p_def ! integration of VB on real axis

 public :: tunneling_int_def  !
 public :: tunneling_and_dos  ! computes of T(E) & dos_proj(E)
 public :: meir_wingreen      ! computes effective T(E) with el-ph
 public :: layer_current      ! computes the current layer-to-layer
 public :: electron_current   ! computes terminal currents
 public :: electron_current_meir_wingreen

 public :: phonon_tunneling   ! computes T(E) for phonons
 public :: phonon_current     ! computes heat currents
 public :: thermal_conductance ! computes thermal conductance

 public :: integrate_el       ! integration of tunneling (el)
 public :: integrate_el_meir_wingreen
 public :: integrate_ph       ! integration of tunneling (ph)


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



contains

  subroutine destroy_en_grid(en_grid)
    type(TEnGrid), dimension(:), allocatable :: en_grid
    if (allocated(en_grid)) deallocate(en_grid)
  end subroutine

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine write_info(verbose,message,Npoints)
    integer, intent(in) :: verbose
    character(*), intent(in) :: message
    integer, intent(in) :: Npoints

     if (id0 .and. verbose.gt.30) then
       write(6,'(a,a,i0,a,a,i0,a)') message,': ',Npoints,' points ', &
            & ' parallelized on ',numprocs,' processes'
     end if

  end subroutine write_info
  !-----------------------------------------------------------------------
  subroutine write_point(verbose,gridpn,Npoints)
    integer, intent(in) :: verbose
    type(TEnGrid), intent(in) :: gridpn
    integer, intent(in) :: Npoints

    if (id0 .and. verbose.gt.VBT) then
      write(6,'(3(a,i0),a,ES15.8)') 'INTEGRAL: point # ',gridpn%pt_path, &
          &'/',Npoints,'  CPU= ', gridpn%cpu, '  E=',real(gridpn%Ec)
    endif

  end subroutine write_point
  !-----------------------------------------------------------------------

  subroutine write_message_clock(verbose,message)
    integer, intent(in) :: verbose
    character(*), intent(in) :: message

    if (id0 .and. verbose.gt.VBT) then
        call message_clock(message)
    end if

  end subroutine write_message_clock
  !-----------------------------------------------------------------------

  subroutine write_end_clock(verbose)
    integer, intent(in) :: verbose

    if (id0 .and. verbose.gt.VBT) call write_clock()

  end subroutine write_end_clock

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

    Nstep = size(negf%en_grid)

    call log_allocate(negf%ldos_mat, Nstep, negf%ndos_proj)
    negf%ldos_mat(:,:)=0.d0

    do i = 1, Nstep

       call write_point(negf%verbose,negf%en_grid(i), size(negf%en_grid))

       if (negf%en_grid(i)%cpu /= id) cycle

       Ec = negf%en_grid(i)%Ec + j*negf%dos_delta
       negf%iE = negf%en_grid(i)%pt

       call compute_Gr(negf, outer, ncont, Ec, Gr)

       call log_allocate(diag, Gr%nrow)
       call getdiag(Gr,diag)
       diag = - aimag(diag)/pi

       do i1 = 1, size(negf%dos_proj)
           negf%ldos_mat(i, i1) = real(sum(diag(negf%dos_proj(i1)%indexes)))
       enddo

       call destroy(Gr)
       call log_deallocate(diag)

    enddo

    !call destroy_en_grid()

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

    if (negf%str%num_conts > 0) then
      kbT = negf%cont(negf%refcont)%kbT_dm
    else
      kbT = negf%kbT
    end if
    muref = negf%muref
    Omega = negf%n_kt * kbT

    if (negf%n_poles.eq.0) then
      Lambda = 0.5d0* kbT * pi
    else
      Lambda = 2.d0* negf%n_poles * KbT * pi
    endif

    Emin = negf%Ec - negf%DeltaEc

    if ((Emin < (muref + 1.d-3)) .and. &
        (Emin > (muref - 1.d-3))) then
       Emin = muref - kbT
    endif

    Npoles = negf%n_poles
    if (Emin > muref) then
      Npoles = 0
    endif

    Ntot = negf%Np_n(1) + negf%Np_n(2) + Npoles
    call destroy_en_grid(negf%en_grid)
    allocate(negf%en_grid(Ntot))

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

    do i = 1, negf%Np_n(1)
      Ec = z1 + pnts(i) * z_diff
      ff = fermi(Ec,muref,KbT)
      zt = negf%g_spin * z_diff * ff * wght(i) / (2.d0 *pi)

      negf%en_grid(i)%path = 1
      negf%en_grid(i)%pt = i
      negf%en_grid(i)%pt_path = i
      negf%en_grid(i)%Ec = Ec
      negf%en_grid(i)%wght = zt
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

    z1 = z2
    z2 = muref + Omega + j*Lambda

    z_diff = z2 - z1

    ioffs = negf%Np_n(1)

    do i = 1, negf%Np_n(2)
      Ec = z1 + pnts(i) * z_diff
      ff = fermi(Ec,muref,KbT)
      zt = negf%g_spin * z_diff * ff * wght(i) / (2.d0 *pi)

      negf%en_grid(ioffs+i)%path = 2
      negf%en_grid(ioffs+i)%pt = ioffs + i
      negf%en_grid(ioffs+i)%pt_path = ioffs + i
      negf%en_grid(ioffs+i)%Ec = Ec
      negf%en_grid(ioffs+i)%wght = zt
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

      negf%en_grid(ioffs+i)%path = 3
      negf%en_grid(ioffs+i)%pt = ioffs + i
      negf%en_grid(ioffs+i)%pt_path = ioffs + i
      negf%en_grid(ioffs+i)%Ec = Ec
      negf%en_grid(ioffs+i)%wght = zt
    enddo

    ! *******************************************************************************
    ! Distribution of Energy grid
    ! pts 1 2 3 4 5 6 7 8 9 ...
    ! cpu 0 1 2 3 0 1 2 3 0 ...
    ! *******************************************************************************
    do i = 0, Ntot-1
       negf%en_grid(i+1)%cpu = mod(i,numprocs)
    enddo

  end subroutine contour_int_n_def



  !-----------------------------------------------------------------------
  ! Contour integration for density matrix, holes
  ! DOES INCLUDE FACTOR 2 FOR SPIN !!
  !-----------------------------------------------------------------------
  ! Performs the complex contur integration
  !
  !      T>=0
  !                          +
  !                          +
  !        * * * * * * * * * * * *
  !                          +   *
  !  ========================+---*---
  !                        Ev mu Emax
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

    if (negf%str%num_conts > 0) then
      kbT = negf%cont(negf%refcont)%kbT_dm
    else
      kbT = negf%kbT
    end if
    muref = negf%muref
    Omega = negf%n_kt * kbT

    if (negf%n_poles.eq.0) then
      Lambda = 0.5d0* kbT * pi
    else
      Lambda = 2.d0* negf%n_poles * KbT * pi
    endif

    Emax = negf%Ev + negf%DeltaEv

    if ((Emax < (muref + 1.d-3)) .and. &
        (Emax > (muref - 1.d-3))) then
       Emax = muref + kbT
    endif

    Npoles = negf%n_poles
    if (Emax < muref) then
      Npoles = 0
    endif

    Ntot=negf%Np_p(1)+negf%Np_p(2)+Npoles
    call destroy_en_grid(negf%en_grid)
    allocate(negf%en_grid(Ntot))

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
      ff = fermi(-Ec,-muref,KbT)   ! 1-f(E-muref)
      zt = - negf%g_spin * z_diff * ff * wght(i) / (2.d0 *pi) !zt is with minus sign because the integration of holes is in
                                                              !the opposite direction compared to the one of electrons
      negf%en_grid(i)%path = 1
      negf%en_grid(i)%pt = i
      negf%en_grid(i)%pt_path = i
      negf%en_grid(i)%Ec = Ec
      negf%en_grid(i)%wght = zt
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
      zt = - negf%g_spin * z_diff * ff * wght(i) / (2.d0 *pi)  !zt is with minus sign because the integration of holes is in  
                                                               !the opposite direction compared to the one of electrons 
      negf%en_grid(ioffs+i)%path = 2
      negf%en_grid(ioffs+i)%pt = ioffs + i
      negf%en_grid(ioffs+i)%pt_path = ioffs + i
      negf%en_grid(ioffs+i)%Ec = Ec
      negf%en_grid(ioffs+i)%wght = zt
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
      zt = j * KbT * negf%g_spin *(1.d0,0.d0)  !zt is with plus sign because the integration of holes is in the opposite
                                               !direction compared to the one of electrons  
      negf%en_grid(ioffs+i)%path = 3
      negf%en_grid(ioffs+i)%pt = ioffs + i
      negf%en_grid(ioffs+i)%pt_path = ioffs + i
      negf%en_grid(ioffs+i)%Ec = Ec
      negf%en_grid(ioffs+i)%wght = zt
    enddo

    ! *******************************************************************************
    ! Distribution of Energy grid
    ! pts 1 2 3 4 5 6 7 8 9 ...
    ! cpu 0 1 2 3 0 1 2 3 0 ...
    ! *******************************************************************************
    do i = 0, Ntot-1
       negf%en_grid(i+1)%cpu = mod(i,numprocs)
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

     if (negf%str%num_conts > 0) then
       kbT = negf%cont(negf%refcont)%kbT_dm
     else
       kbT = negf%kbT
     end if
     muref = negf%muref
     nkT = negf%n_kt * kbT
     Lambda = 2.d0* negf%n_poles * KbT * pi
     mumin = muref - nkT
     Elow = negf%Ec

     Ntot=negf%Np_n(1)+negf%Np_n(2)+negf%n_poles
     !! destroy previously defined grids, if any
     call destroy_en_grid(negf%en_grid)
     allocate(negf%en_grid(Ntot))

     ! ***********************************************************************
     ! 1. INTEGRATION OVER THE CIRCLE Pi..alpha    Np(1)
     ! ***********************************************************************
     ! NEW INTEGRATION FOR COMPLEX DENSITY:
     !----------------------------------------------------
     !   g  [ /           ]     g  [ /         it   ]
     !  --- [ | Gr(z) dz  ] =  --- [ | iGr(t)Re  dt ]
     !  2pi [ /           ]    2pi [ /              ]
     !----------------------------------------------------
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
        negf%en_grid(i)%path=1
        negf%en_grid(i)%pt_path=i
        negf%en_grid(i)%pt=i
        negf%en_grid(i)%Ec=Ec
        negf%en_grid(i)%wght=zt
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
        negf%en_grid(ioffs+i)%path=2
        negf%en_grid(ioffs+i)%pt_path=ioffs+i
        negf%en_grid(ioffs+i)%pt=ioffs+i
        negf%en_grid(ioffs+i)%Ec=Ec
        negf%en_grid(ioffs+i)%wght=zt
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
        negf%en_grid(ioffs+i)%path=3
        negf%en_grid(ioffs+i)%pt_path=ioffs+i
        negf%en_grid(ioffs+i)%pt=ioffs+i
        negf%en_grid(ioffs+i)%Ec=Ec
        negf%en_grid(ioffs+i)%wght=zt
     enddo

     ! *******************************************************************************
     ! Distribution of Energy grid
     ! pts 1 2 3 4 5 6 7 8 9 ...
     ! cpu 0 1 2 3 0 1 2 3 0 ...
     ! *******************************************************************************
     do i = 0, Ntot-1
        negf%en_grid(i+1)%cpu = mod(i,numprocs)
     enddo

  end subroutine contour_int_def

  !-----------------------------------------------------------------------
  subroutine contour_int(negf)
     type(Tnegf) :: negf

     type(z_CSR) :: GreenR, TmpMt
     integer :: i, i1, ncont, Ntot, outer
     real(dp) :: ncyc
     complex(dp) :: Ec, zt

     ncont = negf%str%num_conts
     outer = negf%outer
     Ntot = size(negf%en_grid)
     call create(TmpMt,negf%H%nrow,negf%H%ncol,negf%H%nrow)
     call initialize(TmpMt)

     call write_info(negf%verbose,'CONTOUR INTEGRAL',Ntot)

     do i = 1, Ntot

        call write_point(negf%verbose,negf%en_grid(i), Ntot)
        if (negf%en_grid(i)%cpu .ne. id) cycle

        Ec = negf%en_grid(i)%Ec
        zt = negf%en_grid(i)%wght
        negf%iE = negf%en_grid(i)%pt

        if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Green`s funct ')

        call compute_Gr(negf, outer, ncont, Ec, GreenR)

        if (id0.and.negf%verbose.gt.VBT) call write_clock

        if (allocated(negf%inter)) then
          if (negf%inter%scba_iter > 0) then
            call write_info(negf%verbose,'SCBA iterations',negf%inter%scba_iter)
          end if
        end if

        if(negf%DorE.eq.'E') zt = zt * Ec

        call concat(TmpMt,zt,GreenR,1,1)

        call destroy(GreenR)

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

    integer :: i, i1, np, ioffset, ncont, Ntot
    real(dp), DIMENSION(:), allocatable :: wght,pnts   ! Gauss-quadrature points
    real(dp) :: Omega, mumin, mumax

    ncont = negf%str%num_conts
    ioffset = negf%Np_n(1) + negf%Np_n(2) + negf%n_poles
    ! Omega considers maximum kT so interval is always large enough
    Omega = negf%n_kt * maxval(negf%cont(:)%kbT_dm)

    if (ncont.gt.0) then
       mumin=minval(negf%cont(:)%mu)
       mumax=maxval(negf%cont(:)%mu)
    else
       mumin=negf%muref
       mumax=negf%muref
    endif

    Ntot = negf%Np_real
    !! destroy en_grid from previous calculation, if any
    call destroy_en_grid(negf%en_grid)
    allocate(negf%en_grid(Ntot))

    allocate(pnts(Ntot))
    allocate(wght(Ntot))

    call gauleg(mumin-Omega,mumax+Omega,pnts,wght,Ntot)

    do i = 1, Ntot
       negf%en_grid(i)%path = 1
       negf%en_grid(i)%pt = ioffset + i
       negf%en_grid(i)%pt_path = i
       negf%en_grid(i)%Ec = cmplx(pnts(i),negf%delta,dp)
       negf%en_grid(i)%wght = negf%kwght * negf%g_spin * wght(i)/(2.d0 *pi)
    enddo

    deallocate(wght)
    deallocate(pnts)

    ! distribute energy grid:  0 1 2 3 ... 0 1 2 .... 0 1 2 ....
    do i = 0, Ntot-1
       negf%en_grid(i+1)%cpu = mod(i,numprocs)
    enddo

    ! distribute energy grid:  0 0 0 0 ... 1 1 1 .... 2 2 2 ....
    !np=Ntot/numprocs
    !i1 = 0
    !do i = 1, Ntot, np
    !  negf%en_grid(i:i+np-1)%cpu = i1
    !  i1 = i1 + 1
    !end do
    !
    !if (id .ne. numprocs-1) then
    !  negf%local_en_points = np
    !else
    !  negf%local_en_points = np + mod(Ntot,numprocs)
    !endif

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

    type(z_CSR) :: Gn, TmpMt

    integer :: ref, Npoints
    integer :: i, i1, j1, outer, ncont

    real(dp), DIMENSION(:), allocatable :: frm_f
    real(dp) :: ncyc, Er

    complex(dp) :: zt
    complex(dp) :: Ec

    ncont = negf%str%num_conts
    outer = negf%outer
    Npoints = size(negf%en_grid)

    call log_allocate(frm_f,ncont)

    call create(TmpMt,negf%H%nrow,negf%H%ncol,negf%H%nrow)
    call initialize(TmpMt)

    call write_info(negf%verbose,'REAL AXIS INTEGRAL',Npoints)

    do i = 1, Npoints

       call write_point(negf%verbose,negf%en_grid(i),Npoints)
       if (negf%en_grid(i)%cpu .ne. id) cycle

       Ec = negf%en_grid(i)%Ec
       Er = real(Ec)
       zt = negf%en_grid(i)%wght
       negf%iE = negf%en_grid(i)%pt

       do j1 = 1,ncont
          frm_f(j1)=fermi(Er,negf%cont(j1)%mu,negf%cont(j1)%kbT_dm)
       enddo

       if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Green`s funct ')

       call compute_Gn(negf, outer, ncont, Ec, frm_f, Gn)

       if (id0.and.negf%verbose.gt.VBT) call write_clock

       if (allocated(negf%inter)) then
         if (negf%inter%scba_iter > 0) then
           call write_info(negf%verbose,'SCBA iterations',negf%inter%scba_iter)
         end if
       end if

       if(negf%DorE.eq.'E') zt = zt * Er

       call concat(TmpMt,zt,Gn,1,1)

       call destroy(Gn)


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
    kbT = maxval(negf%cont(:)%kbT_dm)
    Omega = negf%n_kt * kbT
    muref = negf%muref

    if (ncont.gt.0) then
       mumin=minval(negf%cont(:)%mu_n)
       mumax=maxval(negf%cont(:)%mu_n)
    else
       mumin=negf%muref
       mumax=negf%muref
    endif

    Ntot = negf%Np_real
    call destroy_en_grid(negf%en_grid)
    allocate(negf%en_grid(Ntot))

    allocate(pnts(Ntot))
    allocate(wght(Ntot))

    call gauleg(negf%Ec-negf%DeltaEc, max(mumax, negf%Ec) + Omega,pnts,wght,Ntot)

    do i = 1, Ntot
       negf%en_grid(i)%path = 1
       negf%en_grid(i)%pt = ioffset + i
       negf%en_grid(i)%pt_path = i
       negf%en_grid(i)%Ec = cmplx(pnts(i),negf%delta,dp)
       ff = fermi(pnts(i),muref,KbT)
       negf%en_grid(i)%wght = negf%g_spin * negf%kwght * ff * wght(i) / (2.d0 * pi)
    enddo

    deallocate(wght)
    deallocate(pnts)

    ! distribute energy grid
    do i = 0, Ntot-1
       negf%en_grid(i+1)%cpu = mod(i,numprocs)
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
    kbT = maxval(negf%cont(:)%kbT_dm)
    Omega = negf%n_kt * kbT
    muref = negf%muref

    if (ncont.gt.0) then
       mumin=minval(negf%cont(:)%mu_p)
       mumax=maxval(negf%cont(:)%mu_p)
    else
       mumin=negf%muref
       mumax=negf%muref
    endif

    Ntot = negf%Np_real
    call destroy_en_grid(negf%en_grid)
    allocate(negf%en_grid(Ntot))

    allocate(pnts(Ntot))
    allocate(wght(Ntot))

    mumax = negf%Ev+negf%DeltaEv

    call gauleg(mumin-Omega, mumax, pnts, wght, Ntot)

    do i = 1, Ntot
       negf%en_grid(i)%path = 1
       negf%en_grid(i)%pt = ioffset + i
       negf%en_grid(i)%pt_path = i
       negf%en_grid(i)%Ec = cmplx(pnts(i),negf%delta,dp)
       ff = fermi(-pnts(i),-muref,KbT)
       negf%en_grid(i)%wght = negf%g_spin * negf%kwght * ff * wght(i) / (2.d0 * pi)
    enddo

    deallocate(wght)
    deallocate(pnts)

    ! distribute energy grid
    do i = 0, Ntot-1
       negf%en_grid(i+1)%cpu = mod(i,numprocs)
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
  ! Real axis integration for Landauer or Meir-Wingreen
  !-----------------------------------------------------------------------
  subroutine tunneling_int_def(negf)
    type(Tnegf) :: negf

    integer :: i, ncont, Nsteps

    Nsteps=NINT((negf%Emax-negf%Emin)/negf%Estep) + 1
    !! Destroy en_grid from previous calculation, if any
    call destroy_en_grid(negf%en_grid)
    allocate(negf%en_grid(Nsteps))

    do i = 1, Nsteps
       negf%en_grid(i)%path = 1
       negf%en_grid(i)%pt = i
       negf%en_grid(i)%pt_path = i
       negf%en_grid(i)%Ec = cmplx(negf%Emin + negf%Estep*(i-1), 0.0, dp)
       negf%en_grid(i)%wght = negf%kwght
    enddo

    ! distribute energy grid
    do i = 0, Nsteps-1
       negf%en_grid(i+1)%cpu = mod(i,numprocs)
    enddo

  end subroutine tunneling_int_def

  !-----------------------------------------------------------------------------
  !  Routine to compute T(E) and (optionally) dos_proj(E)
  !  PDOS is computed if negf%ndos_proj > 0
  !  When only T(E) is needed, a fast algorithm is used (reduction to one block)
  !-----------------------------------------------------------------------------
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
    Integer :: size_ni             ! emitter-collector contacts

    Integer :: Nstep               ! number of integration points
    Complex(dp) :: Ec              ! Energy point

    Logical :: do_LEDOS            ! performs or not dos_proj


    ! Get out immediately if Emax<Emin
    if (negf%Emax.le.negf%Emin) then
       if(id0) write(*,*) '0 tunneling points;  current = 0.0'
       if (allocated(negf%currents)) call log_deallocate(negf%currents)
       call log_allocate(negf%currents,1)
       negf%currents = 0.0_dp
       return
    endif

    !-------------------------------------------------------

    do_LEDOS = .false.
    if(negf%ndos_proj.gt.0) do_LEDOS=.true.
    ncont = negf%str%num_conts
    Nstep = size(negf%en_grid)
    ncyc=0
    size_ni = size(negf%ni)
    negf%readOldSGF = negf%readOldT_SGFs

    !-------------------------------------------------------

    call log_allocate(TUN_MAT,size_ni)
    !If previous calculation is there, destroy output
    if (allocated(negf%tunn_mat)) then
       call log_deallocate(negf%tunn_mat)
    end if
    call log_allocate(negf%tunn_mat,Nstep,size_ni)

    negf%tunn_mat = 0.0_dp

    if (do_LEDOS) then
       !If previous calculation is there, destroy output
       if (allocated(negf%ldos_mat)) then
         call log_deallocate(negf%ldos_mat)
       end if
       call log_allocate(negf%ldos_mat,Nstep,negf%ndos_proj)
       call log_allocate(LEDOS,negf%ndos_proj)
       negf%ldos_mat(:,:)=0.d0
    endif

    !-------------------------------------------------------
    call write_info(negf%verbose,'CALCULATION OF COHERENT TRANSMISSION',Nstep)

    !Loop on energy points: tunneling
    do i = 1, Nstep

       call write_point(negf%verbose,negf%en_grid(i), size(negf%en_grid))
       if (negf%en_grid(i)%cpu /= id) cycle

       Ec = negf%en_grid(i)%Ec
       negf%iE = negf%en_grid(i)%pt

       if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Contact SE ')
       call compute_contacts(Ec+j*negf%delta,negf,ncyc,Tlc,Tcl,SelfEneR,GS)
       if (id0.and.negf%verbose.gt.VBT) call write_clock

       if (.not.do_LEDOS) then
          if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Tunneling ')

          call calculate_transmissions(negf%H,negf%S,Ec,SelfEneR,negf%ni,negf%nf, &
                             & negf%str, negf%tun_proj, TUN_MAT)

          negf%tunn_mat(i,:) = TUN_MAT(:) * negf%kwght
       else
          if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Tunneling and DOS')
          LEDOS(:) = 0.d0

          call calculate_transmissions_and_dos(negf%H,negf%S,Ec,SelfEneR,GS,negf%ni,negf%nf, &
                             & negf%str, negf%tun_proj, TUN_MAT, negf%dos_proj, LEDOS)

          negf%tunn_mat(i,:) = TUN_MAT(:) * negf%kwght
          negf%ldos_mat(i,:) = LEDOS(:) * negf%kwght
       endif

       if (id0.and.negf%verbose.gt.VBT) call write_clock

       do icont=1,ncont
         call destroy(Tlc(icont),Tcl(icont),SelfEneR(icont),GS(icont))
       enddo

    enddo !Loop on energy

    !call destroy_en_grid()
    call log_deallocate(TUN_MAT)
    if (do_LEDOS) then
       call log_deallocate(LEDOS)
    end if

  end subroutine tunneling_and_dos

  !---------------------------------------------------------------------------
  !>
  !  Calculate the contact current per unit energy according to the
  !  Meir-Wingreen formula on the energy points specified by tunneling_int_def.
  !
  !    I_i(E) = Tr[Sigma^n_i(E)*G^p(E)-Sigma^p_i(E)*G^n(E)]
  !             G^p = A - G^n 
  !    I_i(E) = Tr[Sigma^n_i(E)*A(E)-Gamma_i(E)*G^n(E)]
  !
  !  The solution is calculated on an arbitrary number of
  !  leads and stored on negf%curr_mat. The leads are specified in negf%ni
  !  We don't use the collector negf%nf because we need to specify only the
  !  lead for integration
  !
  ! The Fermi distribution for the contacts can be specified as optional input.
  ! In this case the chemical potential and temperature are ignored and the
  ! values must be constant in energy. This can be used for example to
  ! calculate the effective transmission by forcing the value to be 0 or 1.
  !---------------------------------------------------------------------------
  subroutine meir_wingreen(negf, fixed_occupations)
    type(Tnegf) :: negf
    real(dp), dimension(:), optional :: fixed_occupations

    integer :: scba_iter, i1
    real(dp) :: ncyc
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    Real(dp), Dimension(:), allocatable :: curr_mat
    real(dp), dimension(:), allocatable :: frm
    integer :: size_ni, ii, Nstep, outer, ncont, npl, j1, icont, jj, ref_bk
    complex(dp) :: Ec
    real(dp) :: scba_error
    Type(z_CSR) :: Gn, Gn_previous

    ncont = negf%str%num_conts
    Nstep = size(negf%en_grid)
    size_ni = size(negf%ni)
    negf%readOldSGF = negf%readOldT_SGFs
    outer = 0

    if (.not. allocated(negf%curr_mat)) then
      call log_allocate(negf%curr_mat,Nstep,ncont)
    end if
    negf%curr_mat = 0.0_dp
    call log_allocate(curr_mat, ncont)
    
    ! Create Fermi array. Set reference such that f(ref)=0.
    ref_bk = negf%refcont
    negf%refcont = ncont + 1
    call log_allocate(frm, ncont+1)
    frm = 0.0_dp   
    
    ! Fixed occupations (e.g. 1.0, 0.0) can be used to get a transmission
    if (present(fixed_occupations)) then
      frm(1:ncont) = fixed_occupations
    end if

    call write_info(negf%verbose,'CALCULATION OF MEIR-WINGREEN FORMULA',Nstep)

    !! Loop on energy points
    do ii = 1, Nstep

      call write_point(negf%verbose, negf%en_grid(ii), size(negf%en_grid))

      if (negf%en_grid(ii)%cpu /= id) cycle
      Ec = negf%en_grid(ii)%Ec
      negf%iE = negf%en_grid(ii)%pt
      if (.not.present(fixed_occupations)) then
        do j1 = 1,ncont
           frm(j1)=fermi(real(Ec), negf%cont(j1)%mu, negf%cont(j1)%kbT_t)
        enddo
      end if

      if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Contact SE ')
      call compute_contacts(Ec+j*negf%delta, negf, ncyc, Tlc, Tcl, SelfEneR, GS)
      if (id0.and.negf%verbose.gt.VBT) call write_clock

      ! Calculate the SCBA before meir-wingreen current so el-ph self-energies are stored
      if (allocated(negf%inter)) then

         do scba_iter = 0, negf%inter%scba_niter
            negf%inter%scba_iter = scba_iter
            negf%tDestroyGr = .true.; negf%tDestroyGn = .true.
            call destroy_all_blk(negf)
            negf%tDestroyGr = .false.; negf%tDestroyGn = .false.

            call calculate_Gn_neq_components(negf,real(Ec),SelfEneR,Tlc,Tcl,GS,frm,Gn,outer)
            
            if (negf%inter%scba_iter.ne.0) then
               scba_error = maxval(abs(Gn%nzval - Gn_previous%nzval))
               if (scba_error < negf%inter%scba_tol) then
                  call destroy(Gn)
                  exit
               end if
               call destroy(Gn_previous)
            end if
            call clone(Gn,Gn_previous)
            call destroy(Gn)
         enddo

         call destroy(Gn_previous)

         if (id0 .and. negf%verbose .gt. VBT) then
           if (scba_error < negf%inter%scba_tol) then
              !write(*,*) "SCBA loop converged in",negf%inter%scba_iter,&
              !      & " iterations with error",scba_error
           else
              write(*,*) "WARNING: SCBA exit with error ",scba_error, &
                    & "  > ",negf%inter%scba_tol
           end if
         end if
      endif
            
      negf%tDestroyGr = .true.; negf%tDestroyGn = .true.
      call iterative_meir_wingreen(negf,real(Ec),SelfEneR,frm,curr_mat)
      negf%curr_mat(ii,:) = curr_mat(:) * negf%kwght

      if (id0.and.negf%verbose.gt.VBT) call write_clock
      do icont=1,ncont
        call destroy(Tlc(icont),Tcl(icont),SelfEneR(icont),GS(icont))
      enddo

    enddo
    call log_deallocate(curr_mat)
    call log_deallocate(frm)
    negf%refcont = ref_bk

  end subroutine meir_wingreen

  !---------------------------------------------------------------------------
  !>
  !  Calculate the layer current per unit energy
  !
  !    I_LL'(E) = Tr[(ES-H)_LL' * Gn_L'L(E)-(ES-H)_L'L * Gn_LL'(E)]
  !
  !  The solution is calculated on adjecent layers and stored
  !  negf%tunn_mat. The leads are specified in negf%ni
  !  We don't use the collector negf%nf because we need to specify only the
  !  lead for integration
  !
  !---------------------------------------------------------------------------
  subroutine layer_current(negf)
    type(Tnegf) :: negf

    integer :: nbl, scba_iter, scba_niter, Nstep
    integer :: ii, i1, j1, iK, icont, ncont, ref_bk
    real(dp) :: ncyc, scba_error
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    real(dp), dimension(:), allocatable :: curr_mat, frm
    complex(dp) :: Ec

    ncont = negf%str%num_conts
    Nstep = size(negf%en_grid)
    nbl = negf%str%num_PLs
    negf%readOldSGF = negf%readOldT_SGFs

    ! Allocating curr_mat to nbl-1 so we compute L->L+1 currents
    if (.not. allocated(negf%curr_mat)) then
      call log_allocate(negf%curr_mat, Nstep, nbl-1)
    end if
    negf%curr_mat = 0.0_dp
    call log_allocate(curr_mat, nbl-1)

    ! Create Fermi array. Set reference such that f(ref)=0.
    ref_bk = negf%refcont
    negf%refcont = ncont + 1
    call log_allocate(frm, ncont+1)
    frm = 0.0_dp
    do j1 = 1,ncont
       frm(j1)=fermi(real(Ec), negf%cont(j1)%mu, negf%cont(j1)%kbT_t)
    enddo

    !! Loop over energy points
    enloop: do ii = 1, Nstep

      call write_point(negf%verbose, negf%en_grid(ii), size(negf%en_grid))
      if (negf%en_grid(ii)%cpu /= id) cycle
      Ec = negf%en_grid(ii)%Ec
      negf%iE = negf%en_grid(ii)%pt  ! global energy index
    
      ! ---------------------------------------------------------------------
      ! Compute contact GF
      ! ---------------------------------------------------------------------
      if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Contact SE ')
      call compute_contacts(Ec+j*negf%delta, negf, ncyc, Tlc, Tcl, SelfEneR, GS)
      if (id0.and.negf%verbose.gt.VBT) call write_clock
      !call write_int_info(negf%verbose, VBT, 'Average number of iterations', int(ncyc))

      ! ---------------------------------------------------------------------
      ! Compute block tri-diagonal Gr and Gn
      ! ---------------------------------------------------------------------
      ! Avoids cleanup of Gn and ESH components for later use
      negf%tDestroyGn = .false.; negf%tDestroyESH = .false.
      if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Gn ')
      call calculate_Gn_neq_components(negf,real(Ec),SelfEneR,Tlc,Tcl,GS,frm)
      if (id0.and.negf%verbose.gt.VBT) call write_clock

      ! ---------------------------------------------------------------------
      ! Compute layer-to-layer currents and release memory
      ! ---------------------------------------------------------------------
      negf%tDestroyGn = .true.; negf%tDestroyESH = .true.
      if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Jn,n+1 ')
      call iterative_layer_current(negf,real(Ec),curr_mat)
      if (id0.and.negf%verbose.gt.VBT) call write_clock

      negf%curr_mat(ii,:) = curr_mat(:) * negf%kwght 
    
      do icont=1,ncont
        call destroy(Tlc(icont),Tcl(icont),SelfEneR(icont),GS(icont))
      end do

    end do enloop

      
    print*,'call compute_electron_current  CPU#', id
    call electron_current_meir_wingreen(negf)


    call log_deallocate(curr_mat)
    call log_deallocate(frm)
    negf%refcont = ref_bk

  end subroutine layer_current

  !---------------------------------------------------------------------------
  !>
  !  Calculate the equilibrium Retarded Green's function (extended diagonal)
  !  on a single energy point
  !  It groups calculation of leads, scba loop if any and deallocations of
  !  working arrays. This routine is used in contour integration and DOS and
  !
  !---------------------------------------------------------------------------
  subroutine compute_Gr(negf, outer, ncont, Ec, Gr)
    type(Tnegf), intent(inout) :: negf
    Type(z_CSR), intent(out) :: Gr
    complex(dp), intent(in) :: Ec
    integer, intent(in) :: outer, ncont

    integer :: scba_iter, i1
    real(dp) :: ncyc
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS

    real(dp) :: scba_error
    Type(z_CSR) :: Gr_previous

    negf%readOldSGF = negf%readOldDM_SGFs

    call compute_contacts(Ec,negf,ncyc,Tlc,Tcl,SelfEneR,GS)

    call calculate_Gr(negf,Ec,SelfEneR,Tlc,Tcl,GS,Gr,outer)

    if (allocated(negf%inter)) then
      if (negf%inter%scba_niter /= 0) then
        call clone(Gr,Gr_previous)

        do scba_iter = 1, negf%inter%scba_niter
          negf%inter%scba_iter = scba_iter
          call destroy(Gr)
          call calculate_Gr(negf,Ec,SelfEneR,Tlc,Tcl,GS,Gr,outer)

          scba_error = maxval(abs(Gr%nzval - Gr_previous%nzval))

          if (scba_error .lt. negf%inter%scba_tol) then
            exit
          end if

          call destroy(Gr_previous)
          call clone(Gr,Gr_previous)
        end do
      end if
    end if

    do i1=1,ncont
      call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
    enddo

  end subroutine compute_Gr

  !------------------------------------------------------------------------------
  !>
  !  Calculates the non equilibrium LESSER Green function (extended diagonal)
  !  on a single energy point on real axis.
  !  It groups calculation of leads, scba loop if any and deallocations of
  !  working arrays.
  !
  !-----------------------------------------------------------------------------
  subroutine compute_Gn(negf, outer, ncont, Ec, frm, Gn)
    type(Tnegf), intent(inout) :: negf
    Type(z_CSR), intent(out) :: Gn
    complex(dp), intent(in) :: Ec
    real(dp), dimension(:), intent(in) :: frm

    integer, intent(in) :: outer, ncont
    integer :: scba_iter, i1, max_scba_iter
    real(dp) :: ncyc
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    real(dp) :: Er

    !DAR begin - compute_Gr
    real(dp) :: scba_error
    Type(z_CSR) :: Gn_previous
    !DAR end

    negf%readOldSGF = negf%readOldDM_SGFs
    Er = real(Ec,dp)
    call compute_contacts(Ec,negf,ncyc,Tlc,Tcl,SelfEneR,GS)

    call calculate_Gn_neq_components(negf, Er, SelfEneR, Tlc, Tcl, GS, frm, Gn, outer)

    ! In case of interactions (only elastic supported now) we go into
    ! the Self Consistent Born Approximation loop.
    if (.not.allocated(negf%inter)) then
      max_scba_iter = 0
    else
      max_scba_iter = negf%inter%scba_niter
      negf%inter%scba_iter = 0
    end if

    do scba_iter = 1, max_scba_iter
      negf%inter%scba_iter = scba_iter
      call clone(Gn,Gn_previous)
      call destroy(Gn)
      call calculate_Gn_neq_components(negf, Er, SelfEneR, Tlc, Tcl, GS, frm, Gn, outer)
      scba_error = maxval(abs(Gn%nzval - Gn_previous%nzval))
      call destroy(Gn_previous)
      if (scba_error .lt. negf%inter%scba_tol) then
        exit
      end if
    enddo

    do i1=1,ncont
      call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
    enddo

  end subroutine compute_Gn

  !---------------------------------------------------------------------------
  !   COMPUTATION OF CURRENTS - INTEGRATION OF T(E)
  !---------------------------------------------------------------------------
  subroutine electron_current(negf)
    type(Tnegf) :: negf

    integer :: size_ni, ii, ni, nf
    real(dp) :: mu1, mu2

    if (.not.allocated(negf%tunn_mat)) then
      write(*,*) 'Internal error: electron_current must be invoked'
      write(*,*) 'after tunneling calculation'
      stop
    end if

    size_ni = size(negf%tunn_mat,2)

    !print *, 'negf%ni',negf%ni
    !print *, 'negf%nf',negf%nf
    !print *, 'negf%ref',negf%refcont
    !print *, 'negf%cont',size(negf%cont)

    ! If previous calculation is there, destroy it
    if (allocated(negf%currents)) call log_deallocate(negf%currents)
    call log_allocate(negf%currents,size_ni)

    negf%currents=0.d0

    if (size(negf%cont) < 2) then
      return
    end if

    do ii=1,size_ni
       ni = negf%ni(ii); nf = negf%nf(ii)
       mu1=negf%cont(ni)%mu; mu2=negf%cont(nf)%mu
       negf%currents(ii)= integrate_el(negf%tunn_mat(:,ii), mu1, mu2, &
                          & negf%cont(ni)%kbT_t, negf%cont(nf)%kbT_t, &
                          & negf%Emin, negf%Emax, negf%Estep) * negf%g_spin
    enddo

  end subroutine electron_current

  !---------------------------------------------------------------------------
  !   COMPUTATION OF CURRENTS - INTEGRATION OF I_i(E)
  !---------------------------------------------------------------------------
  subroutine electron_current_meir_wingreen(negf)
    type(Tnegf) :: negf

    integer :: size_ni, ii
    real(dp) :: mu1, mu2

    if (.not.allocated(negf%curr_mat)) then
      write(*,*) 'Internal error: electron_current_meir_wingreen must be invoked'
      write(*,*) 'after tunneling calculation'
      stop
    end if

    size_ni = size(negf%curr_mat,2)

    ! If previous calculation is there, destroy it
    if (allocated(negf%currents)) call log_deallocate(negf%currents)
    call log_allocate(negf%currents,size_ni)

    negf%currents=0.d0
    do ii=1,size_ni
       negf%currents(ii)= integrate_el_meir_wingreen(negf%curr_mat(:,ii), &
                          & negf%Emin, negf%Emax, negf%Estep) * negf%g_spin
    enddo

  end subroutine electron_current_meir_wingreen


  !-----------------------------------------------------------------------
  !  Routine to compute T(E) and (optionally) dos_proj(E)
  !  PDOS is computed if negf%ndos_proj > 0
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
    Integer :: size_ni             ! emitter-collector contacts

    Integer :: Nstep               ! number of integration points
    Complex(dp) :: Ec              ! Energy point
    Complex(dp) :: delta
    Logical :: do_LEDOS            ! performs or not dos_proj

    ! Get out immediately if Emax<Emin
    if (negf%Emax.le.negf%Emin) then
       if(id0) write(*,*) '0 tunneling points;  current = 0.0'
       call log_allocate(negf%currents,1)
       negf%currents = 0.0_dp
       return
    endif
    !-------------------------------------------------------

    do_LEDOS = .false.
    if(negf%ndos_proj.gt.0) do_LEDOS=.true.
    ncont = negf%str%num_conts
    Nstep = size(negf%en_grid)
    ncyc=0
    size_ni = size(negf%ni)
    negf%readOldSGF = negf%readOldT_SGFs

    !-------------------------------------------------------

    call log_allocate(TUN_MAT,size_ni)
    call log_allocate(negf%tunn_mat,Nstep,size_ni)
    negf%tunn_mat = 0.0_dp

    if (do_LEDOS) then
       call log_allocate(negf%ldos_mat,Nstep,negf%ndos_proj)
       call log_allocate(LEDOS,negf%ndos_proj)
       negf%ldos_mat(:,:)=0.d0
    endif
    !-------------------------------------------------------


    !Loop on energy points: tunneling
    do i = 1, Nstep

       call write_point(negf%verbose,negf%en_grid(i), size(negf%en_grid))
       if (negf%en_grid(i)%cpu /= id) cycle

       Ec = negf%en_grid(i)%Ec * negf%en_grid(i)%Ec
       negf%iE = negf%en_grid(i)%pt

       select case(negf%deltaModel)
       case(DELTA_SQ)
         delta = negf%delta * negf%delta
       case(DELTA_W)
         delta = negf%delta * real(negf%en_grid(i)%Ec)
       case(DELTA_MINGO)
         delta = negf%delta * (1.0_dp - real(negf%en_grid(i)%Ec)/(negf%wmax+1d-12)) * Ec
       end select

       if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Contact SE ')
       call compute_contacts(Ec+j*delta,negf,ncyc,Tlc,Tcl,SelfEneR,GS)
       if (id0.and.negf%verbose.gt.VBT) call write_clock


       if (.not.do_LEDOS) then
          if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Tunneling ')

          call calculate_transmissions(negf%H,negf%S,Ec,SelfEneR,negf%ni,negf%nf, &
                             & negf%str, negf%tun_proj, TUN_MAT)

          negf%tunn_mat(i,:) = TUN_MAT(:) * negf%kwght
       else
          if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Tunneling and DOS')
          LEDOS(:) = 0.d0

          call calculate_transmissions_and_dos(negf%H,negf%S,Ec,SelfEneR,GS,negf%ni,negf%nf, &
                             & negf%str, negf%tun_proj, TUN_MAT, negf%dos_proj, LEDOS)

          negf%tunn_mat(i,:) = TUN_MAT(:) * negf%kwght
          negf%ldos_mat(i,:) = LEDOS(:) * negf%kwght
       endif

       if (id0.and.negf%verbose.gt.VBT) call write_clock

       do icont=1,ncont
         call destroy(Tlc(icont),Tcl(icont),SelfEneR(icont),GS(icont))
       enddo

    enddo !Loop on energy

    !call destroy_en_grid()
    call log_deallocate(TUN_MAT)
    if(do_LEDOS) call log_deallocate(LEDOS)

  end subroutine phonon_tunneling

  !---------------------------------------------------------------------------
  subroutine phonon_current(negf)
    type(Tnegf) :: negf

    integer :: size_ni, ii, ni, nf

    size_ni = size(negf%tunn_mat,2)

    call log_allocate(negf%currents,size_ni)
    negf%currents=0.d0

    do ii=1,size_ni
       ni = negf%ni(ii); nf = negf%nf(ii)
       negf%currents(ii)= integrate_ph(negf%tunn_mat(:,ii),  &
                            & negf%cont(ni)%kbT_t, negf%cont(nf)%kbT_t, &
                            & negf%Emin, negf%Emax, negf%Estep)
    enddo


  end subroutine phonon_current



  !////////////////////////////////////////////////////////////////////////
  !************************************************************************
  ! Function to integrate the tunneling and get the current
  ! The function resolves fermi(E) on a fine grid interpolating linearly T(E)
  ! In this way a more precise integration is obtained when T ~ constant
  !************************************************************************
  function integrate_el(TUN_TOT,mu1,mu2,kT1,kT2,emin,emax,estep)

    implicit none

    real(dp) :: integrate_el
    real(dp), intent(inout) :: mu1, mu2
    real(dp), intent(in) :: emin,emax,estep
    real(dp), dimension(:), intent(in) :: TUN_TOT
    real(dp), intent(in) :: kT1, kT2

    REAL(dp) :: destep,kbT1,kbT2,TT1,TT2,E3,E4,TT3,TT4
    REAL(dp) :: E1,E2,c1,c2,curr
    INTEGER :: i,i1,N,Nstep,imin,imax
    logical :: swapped

    curr=0.d0
    N=0
    destep=1.0d10
    Nstep=NINT((emax-emin)/estep);

    !We set a minimum possible value T = 1.0 K to avoid
    !numerical issues
    if (kT1 < 1.0_dp*Kb) then
      kbT1 = Kb*1.0_dp
    else
      kbT1 = kT1
    endif
    if (kT2 < 1.0_dp*Kb) then
      kbT2 = Kb*1.0_dp
    else
      kbT2 = kT2
    endif

    swapped = .false.
    if (mu2<mu1) then
      !We have to remember we swapped to get the right sign at the end
      swapped = .true.
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

          curr=curr+(c1+c2)*(E4-E3)/2.d0

       enddo

    enddo

    if (swapped) curr = -1.d0*curr
    integrate_el = curr

  end function integrate_el

  !************************************************************************
  ! Function to integrate the current density I(E) !!! and get the current
  ! for meir_wingreen
  !************************************************************************
  function integrate_el_meir_wingreen(TUN_TOT,emin,emax,estep)

    implicit none

    real(dp) :: integrate_el_meir_wingreen
    real(dp), intent(in) :: emin,emax,estep
    real(dp), dimension(:), intent(in) :: TUN_TOT

    REAL(dp) :: TT1,TT2,E3,E4,TT3,TT4
    REAL(dp) :: E1,E2,c1,c2,curr
    INTEGER :: i,i1,N,Nstep,imin,imax

    curr=0.d0
    N=0
    Nstep=NINT((emax-emin)/estep);

    imin=0
    imax=Nstep

    ! performs the integration with simple trapezium rule.
    do i=imin,imax-1

       E1=emin+estep*i
       TT1=TUN_TOT(i+1)
       E2=emin+estep*(i+1)
       TT2=TUN_TOT(i+2)

       curr=curr+(TT1+TT2)*(E2-E1)/2.d0

    enddo

    integrate_el_meir_wingreen = curr

  end function integrate_el_meir_wingreen

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
        dd = real(DensMat%nzval(ka))
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

    tmp=x1
    x1=x2
    x2=tmp

  end subroutine swap

  subroutine quasiEq_int_n(negf, mu_n, E_half, rho)
    !In/Out
    type(Tnegf), intent(inout) :: negf
    real(dp), dimension(:), intent(inout) :: rho
    real(dp), dimension(:), intent(in) :: E_half, mu_n

    !Work
    integer :: i, nr
    integer ::  ncont, outer
    integer ::  Ntot, N_coarse, N_refined, Nz
    complex(dp) :: Ez
    type(z_CSR) :: Gr 
    real(dp), dimension(:), allocatable :: wght, E
    real(dp), dimension(:), allocatable :: minE, maxE
    complex(dp), dimension(:), allocatable :: diag   
    real(dp) :: Omega, ff, ww, kbT 

    kbT = maxval(negf%cont(:)%kbT_dm)
    ncont = negf%str%num_conts
    outer = negf%outer
    
    ! Omega considers maximum kT so interval is always large enough
    Omega = negf%n_kt * kbT
    Nz = size(mu_n)

    allocate(minE(Nz))
    allocate(maxE(Nz))

    !First range for electrons: from E_half to mu_n (if E_half < mu_n); coarse integration
    if (all(E_half < mu_n)) then
       minE(:) = E_half(:)
       maxE(:) = mu_n(:)
    else
       minE(:) = mu_n(:)
       maxE(:) = E_half(:)
    endif
    N_coarse = nint(abs(maxval(maxE)-minval(minE)) / negf%Estep_coarse) 

    allocate(E(N_coarse))
    allocate(wght(N_coarse))

    call gauleg(minval(minE), maxval(maxE), E, wght, N_coarse)

    do i = 1, N_coarse 
 
       if (mod(i-1,numprocs) .ne. id) cycle
       
       Ez = cmplx(E(i),negf%delta,dp)
       call compute_Gr(negf, outer, ncont, Ez, Gr)
       call log_allocate(diag, Gr%nrow)
       call getdiag(Gr,diag)

       do nr = 1,Nz
          if (E(i) > minE(nr) .and. E(i) < maxE(nr)) then
             ff = fermi(E(i), mu_n(nr), kbT)
             ww = negf%g_spin * negf%kwght * wght(i) * ff/pi   
             rho(nr) = rho(nr) - aimag(diag(nr)) * ww
          endif
       enddo

       call log_deallocate(diag)
       call destroy(Gr)

    enddo

    deallocate(wght)
    deallocate(E)

    !Second range, from maxE of previous integration to mu_n + nKT; refined integration
    !Redefinition of minE,maxE
    minE(:) = maxE(:)
    maxE(:) = mu_n(:) + Omega

    N_refined = nint(abs(maxval(maxE) - minval(minE)) / negf%Estep)  
    Ntot = N_coarse + N_refined  

    allocate(E(N_refined))
    allocate(wght(N_refined))

    call gauleg(minval(minE), maxval(maxE), E, wght, N_refined)

    do i = 1, N_refined 

       if (mod(i-1,numprocs) .ne. id) cycle

       Ez = cmplx(E(i),negf%delta,dp)
       call compute_Gr(negf, outer, ncont, Ez, Gr)
       call log_allocate(diag, Gr%nrow)
       call getdiag(Gr,diag)

       do nr = 1,Nz
          if (E(i) > minE(nr) .and. E(i) < maxE(nr)) then
             ff = fermi(E(i), mu_n(nr), kbT)
             ww = negf%g_spin * negf%kwght * wght(i) * ff/pi   
             rho(nr) = rho(nr) - aimag(diag(nr)) * ww
          endif
       enddo

       call log_deallocate(diag)
       call destroy(Gr)

    enddo

    deallocate(E)
    deallocate(wght)

    deallocate(minE)
    deallocate(maxE)
  end subroutine quasiEq_int_n

  subroutine quasiEq_int_p(negf, mu_p, E_half, rho)
    !In/Out
    type(Tnegf), intent(inout) :: negf
    real(dp), dimension(:), intent(inout) :: rho
    real(dp), dimension(:), intent(in) :: E_half, mu_p

    !Work
    integer :: i, nr
    integer ::  ncont, outer
    integer ::  Ntot, N_coarse, N_refined, Nz
    complex(dp) :: Ez
    type(z_CSR) :: Gr 
    real(dp), dimension(:), allocatable :: wght, E
    real(dp), dimension(:), allocatable :: minE, maxE
    complex(dp), dimension(:), allocatable :: diag   
    real(dp) :: Omega, ff, ww, kbT 

    kbT = maxval(negf%cont(:)%kbT_dm)
    ncont = negf%str%num_conts
    outer = negf%outer
    
    ! Omega considers maximum kT so interval is always large enough
    Omega = negf%n_kt * kbT
    Nz = size(mu_p)

    allocate(minE(Nz))
    allocate(maxE(Nz))

    !First range for holes: from mu_p to E_half (if E_half > mu_p); coarse integration
    if (all(E_half > mu_p)) then
       maxE(:) = E_half(:)
       minE(:) = mu_p(:)
    else
       maxE(:) = mu_p(:)
       minE(:) = E_half(:)
    endif
    N_coarse = nint(abs(maxval(maxE)-minval(minE)) / negf%Estep_coarse) 

    allocate(E(N_coarse))
    allocate(wght(N_coarse))

    call gauleg(minval(minE), maxval(maxE), E, wght, N_coarse)

    do i = 1, N_coarse 
 
       if (mod(i-1,numprocs) .ne. id) cycle
       
       Ez = cmplx(E(i),negf%delta,dp)
       call compute_Gr(negf, outer, ncont, Ez, Gr)
       call log_allocate(diag, Gr%nrow)
       call getdiag(Gr,diag)

       do nr = 1,Nz
          if (E(i) > minE(nr) .and. E(i) < maxE(nr)) then
             ff = fermi(-E(i), -mu_p(nr), kbT)
             ww = negf%g_spin * negf%kwght * wght(i) * ff/pi   
             rho(nr) = rho(nr) - aimag(diag(nr)) * ww
          endif
       enddo

       call log_deallocate(diag)
       call destroy(Gr)

    enddo

    deallocate(wght)
    deallocate(E)

    !Second range, from maxE of previous integration to mu_n + nKT; refined integration
    !Redefinition of minE,maxE
    maxE(:) = minE(:)
    minE(:) = mu_p(:) - Omega

    N_refined = nint(abs(maxval(maxE) - minval(minE)) / negf%Estep)  
    Ntot = N_coarse + N_refined  

    allocate(E(N_refined))
    allocate(wght(N_refined))

    call gauleg(minval(minE), maxval(maxE), E, wght, N_refined)

    do i = 1, N_refined 

       if (mod(i-1,numprocs) .ne. id) cycle

       Ez = cmplx(E(i),negf%delta,dp)
       call compute_Gr(negf, outer, ncont, Ez, Gr)
       call log_allocate(diag, Gr%nrow)
       call getdiag(Gr,diag)

       do nr = 1,Nz
          if (E(i) > minE(nr) .and. E(i) < maxE(nr)) then
             ff = fermi(-E(i), -mu_p(nr), kbT)
             ww = negf%g_spin * negf%kwght * wght(i) * ff/pi   
             rho(nr) = rho(nr) - aimag(diag(nr)) * ww
          endif
       enddo

       call log_deallocate(diag)
       call destroy(Gr)

    enddo

    deallocate(E)
    deallocate(wght)

    deallocate(minE)
    deallocate(maxE)
  end subroutine quasiEq_int_p

end module integrations
