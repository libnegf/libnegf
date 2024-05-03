!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Greens functions.         !
!! Copyright (C) 2012 - 2026                                                !
!!                                                                          !
!! This file is part of libNEGF: a library for                              !
!! Non Equilibrium Green's Functions calculations                           !
!!                                                                          !
!! Developers: Alessandro Pecchia, Daniele Soccodato                        !
!! Former Contributors: Gabriele Penazzi, Luca Latessa, Aldo Di Carlo       !
!!                                                                          !
!! libNEGF is free software: you can redistribute and/or modify it          !
!! under the terms of the GNU Lesser General Public License as published    !
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
 use mat_def
 use ln_extract
 use contselfenergy
 use clock
 use energy_mesh
 use interactions
 use ln_elastic
 use ln_inelastic
 use elphinel
 use ln_enums, only : integration_type
#:if defined("MPI")
  use libmpifx_module
#:endif

 implicit none

 private

 public :: contour_int       ! generalized contour integrator
 public :: contour_int_inel  ! generalized contour integrator with inelasti scatt.
 public :: real_axis_int     ! real-axis integrator
 public :: real_axis_int_inel! real-axis integrator with inelastic scattering
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
  subroutine write_info(verbose,vbmin,message,int_num)
    integer, intent(in) :: verbose
    integer, intent(in) :: vbmin
    character(*), intent(in) :: message
    integer, optional, intent(in) :: Int_num
    if (id0 .and. verbose.gt.vbmin) then
      if (present(int_num)) then
        write(6,'(a,a,i0)') message,': ',int_num
      else
        write(6,'(a)') message
      end if
    end if
  end subroutine write_info
  !-----------------------------------------------------------------------

  subroutine write_info_parallel(verbose,vbmin,message,Npoints)
    integer, intent(in) :: verbose
    integer, intent(in) :: vbmin
    character(*), intent(in) :: message
    integer, intent(in) :: Npoints

     if (id0 .and. verbose.gt.vbmin) then
       write(6,'(a,a,i0,a,a,i0,a)') message,': ',Npoints,' points ', &
            & ' parallelized on ',numprocs,' processes'
     end if

  end subroutine write_info_parallel
  !-----------------------------------------------------------------------

  subroutine write_real_info(verbose,vbmin,message,rr)
    integer, intent(in) :: verbose
    integer, intent(in) :: vbmin
    character(*), intent(in) :: message
    real(dp), intent(in) :: rr

     if (id0 .and. verbose.gt.vbmin) then
       write(6,'(a,a,es12.3)') message,': ',rr
     end if

  end subroutine write_real_info

  subroutine write_int_info(verbose,vbmin,message,ii)
    integer, intent(in) :: verbose
    integer, intent(in) :: vbmin
    character(*), intent(in) :: message
    integer, intent(in) :: ii

     if (id0 .and. verbose.gt.vbmin) then
       write(6,'(a,a,i3)') message,': ',ii
     end if

  end subroutine write_int_info

  !-----------------------------------------------------------------------
  subroutine write_Epoint(verbose,gridpn,Npoints)
    integer, intent(in) :: verbose
    type(TEnGrid), intent(in) :: gridpn
    integer, intent(in) :: Npoints

    if (id0 .and. verbose.gt.VBT) then
      write(6,'(3(a,i0),a,ES15.8)') 'INTEGRAL: point # ',gridpn%pt, &
          &'/',Npoints,'  CPU= ', gridpn%cpu, '  E=',real(gridpn%Ec)
    endif

  end subroutine write_Epoint
  !-----------------------------------------------------------------------

  subroutine write_kpoint(verbose,local_point,global_point)
    integer, intent(in) :: verbose
    integer, intent(in) :: local_point
    integer, intent(in) :: global_point

    if (id0 .and. verbose.gt.VBT) then
      write(6,'(a,i0,a,i0)') 'K-Point: local point # ', local_point, &
          &' -> global ',global_point
    endif

  end subroutine write_kpoint
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
    Type(z_CSR) ::  Gr, TmpMt
    complex(dp), Dimension(:), ALLOCATABLE :: diag

    integer :: Npoints, i, i1, l, kb, ke
    integer :: outer, ncont

    real(dp) :: ncyc, scba_error
    complex(dp) :: Ec
    character(6) :: ofKP
    character(1) :: ofSp

    outer = 1
    ncont = negf%str%num_conts

    Npoints = size(negf%en_grid)

    call log_allocate(negf%ldos_mat, Npoints, negf%ndos_proj)
    negf%ldos_mat(:,:)=0.0_dp

    do i = 1, Npoints

       call write_Epoint(negf%verbose,negf%en_grid(i), size(negf%en_grid))
       if (negf%en_grid(i)%cpu /= id) cycle

       Ec = negf%en_grid(i)%Ec + j*negf%dos_delta
       negf%iE = negf%en_grid(i)%pt
       negf%iE_path = negf%en_grid(i)%pt_path

       call compute_Gr(negf, outer, ncont, Ec, scba_error, Gr)

       call prealloc_mult(Gr, negf%S, TmpMt)
       call log_allocate(diag, Gr%nrow)
       call getdiag(TmpMt, diag)
       diag = - negf%g_spin * aimag(diag)/pi

       do i1 = 1, size(negf%dos_proj)
           negf%ldos_mat(i, i1) = real(sum(diag(negf%dos_proj(i1)%indexes)))
       enddo

       call destroy(Gr)
       call destroy(TmpMt)
       call log_deallocate(diag)

    enddo

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

    integer :: i, Npoints, Npoles, ioffs
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
      Lambda = 0.5_dp* kbT * pi
    else
      Lambda = 2.0_dp* negf%n_poles * KbT * pi
    endif

    Emin = negf%Ec - negf%DeltaEc

    if ((Emin < (muref + 0.001_dp)) .and. &
        (Emin > (muref - 0.001_dp))) then
       Emin = muref - kbT
    endif

    Npoles = negf%n_poles
    if (Emin > muref) then
      Npoles = 0
    endif

    Npoints = negf%Np_n(1) + negf%Np_n(2) + Npoles
    call destroy_en_grid(negf%en_grid)
    allocate(negf%en_grid(Npoints))

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

    call gauleg(0.0_dp,1.0_dp,pnts,wght,negf%Np_n(1))

    do i = 1, negf%Np_n(1)
      Ec = z1 + pnts(i) * z_diff
      ff = fermi(Ec,muref,KbT)
      zt = negf%g_spin * z_diff * ff * wght(i) / (2.0_dp *pi)

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

    call gauleg(0.0_dp,1.0_dp,pnts,wght,negf%Np_n(2))    !Setting weights for integration

    z1 = z2
    z2 = muref + Omega + j*Lambda

    z_diff = z2 - z1

    ioffs = negf%Np_n(1)

    do i = 1, negf%Np_n(2)
      Ec = z1 + pnts(i) * z_diff
      ff = fermi(Ec,muref,KbT)
      zt = negf%g_spin * z_diff * ff * wght(i) / (2.0_dp *pi)

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
      Ec = muref + j * KbT *pi* (2.0_dp*i - 1.0_dp)
      zt= -j * KbT * negf%g_spin *(1.0_dp,0.0_dp)

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
    do i = 0, Npoints-1
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

    integer :: i, Npoints, Npoles, ioffs
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
      Lambda = 2.0_dp* negf%n_poles * KbT * pi
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

    Npoints=negf%Np_p(1)+negf%Np_p(2)+Npoles
    call destroy_en_grid(negf%en_grid)
    allocate(negf%en_grid(Npoints))

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

    call gauleg(0.0_dp,1.0_dp,pnts,wght,negf%Np_p(1))

    do i = 1, negf%Np_p(1)
      Ec = z1 + pnts(i) * z_diff
      ff = fermi(-Ec,-muref,KbT)   ! 1-f(E-muref)
      zt = - negf%g_spin * z_diff * ff * wght(i) / (2.0_dp *pi) !zt is with minus sign because the integration of holes is in
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

    call gauleg(0.0_dp,1.0_dp,pnts,wght,negf%Np_p(2))    !Setting weights for integration

    z1 = z2
    z2 = muref - Omega + j*Lambda

    z_diff = z2 - z1

    ioffs = negf%Np_p(1)

    do i = 1, negf%Np_p(2)
      Ec = z1 + pnts(i) * z_diff
      ff = fermi(-Ec,-muref,KbT)
      !zt is with minus sign because the integration of holes is in
      !the opposite direction compared to the one of electrons
      zt = - negf%g_spin * z_diff * ff * wght(i) / (2.0_dp *pi)
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
      Ec = muref + j * KbT *pi* (2.0_dp*i - 1.0_dp)
      zt = j * KbT * negf%g_spin *(1.0_dp,0.0_dp)  !zt is with plus sign because the integration of holes is in the opposite
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
    do i = 0, Npoints-1
       negf%en_grid(i+1)%cpu = mod(i,numprocs)
    enddo

  end subroutine contour_int_p_def

  !-----------------------------------------------------------------------
  ! Contour integration for density matrix
  ! DOES INCLUDE FACTOR 2 FOR SPIN !!
  !-----------------------------------------------------------------------
  ! Performs the complex contur integration.
  ! Coherent transport:
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
  !  With inelastic scattering (wq > 0) the contour goes as:
  !
  !         * * * *
  !      *          *
  !    *              *
  !   *                *
  !  -*----------------*-------
  !  Elow         mumin - nkT - wq
  !
  !-----------------------------------------------------------------------
  subroutine contour_int_def(negf)
     type(Tnegf) :: negf

     integer :: i, Npoints, ioffs, nPoles
     real(dp) :: Lambda, Rad, Centre
     real(dp), DIMENSION(:), ALLOCATABLE :: wght,pnts   ! Gauss-quadrature points
     real(dp) :: dt, Elow
     real(dp) :: muref, mumin, kbT, nkT, alpha, wqmax = 0.0_dp
     complex(dp) :: z1,z2,z_diff, zt
     complex(dp) :: Ec, ff, Pc

     if (negf%str%num_conts > 0) then
       kbT = negf%cont(negf%refcont)%kbT_dm
       mumin = minval(negf%cont(:)%mu)
     else
       kbT = negf%kbT
       mumin = negf%mu
     end if
     nkT = negf%n_kt * kbT
     Elow = negf%Ec
     muref = negf%muref
     ! Get maximum wq for inelastic interactions
     wqmax = get_max_wq(negf%interactList)
     if (wqmax > 0.0_dp) then
       mumin = mumin - nkT - wqmax
       Lambda = 0.0_dp
       nPoles = 0
     else
       mumin = muref - nkT
       nPoles = negf%n_poles
       Lambda = 2.0_dp* negf%n_poles * KbT * pi
     end if

     Npoints=negf%Np_n(1) + negf%Np_n(2) + nPoles
     !! destroy previously defined grids, if any
     call destroy_en_grid(negf%en_grid)
     allocate(negf%en_grid(Npoints))

     ! ***********************************************************************
     ! 1. INTEGRATION OVER THE CIRCLE Pi..alpha    Np(1)
     ! ***********************************************************************
     ! NEW INTEGRATION FOR COMPLEX DENSITY:
     !----------------------------------------------------
     !   g  [ /           ]     g  [ /         it   ]
     !  --- [ | Gr(z) dz  ] =  --- [ | iGr(t)Re  dt ]
     !  2pi [ /           ]    2pi [ /              ]
     !----------------------------------------------------
     Centre = (Lambda**2-Elow**2+(mumin)**2)/(2.0_dp*(mumin-Elow))
     Rad = Centre - Elow
     if (kbT.eq.0.0_dp .or. wqmax>0.0_dp) then
        alpha = 0.1_dp*pi
     else
        alpha = atan(Lambda/(mumin-Centre))
     end if

     !Setting weights for gaussian integration
     allocate(wght(negf%Np_n(1)))
     allocate(pnts(negf%Np_n(1)))

     call gauleg(pi,alpha,pnts,wght,negf%Np_n(1))
     do i = 1, negf%Np_n(1)
        Pc = Rad*exp(j*pnts(i))
        Ec = Centre+Pc
        zt = j * Pc * negf%g_spin * wght(i)/(2.0_dp*pi)
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

     if (kbT.eq.0.0_dp .or. wqmax>0.0_dp) then     ! Circle integration T=0
       call  gauleg(alpha,0.0_dp,pnts,wght,negf%Np_n(2))
     else                                          ! Segment integration T>0
       z1 = muref + nkT + j*Lambda
       z2 = muref - nkT + j*Lambda
       z_diff = z2 - z1
       call  gauleg(1.0_dp,0.0_dp,pnts,wght,negf%Np_n(2))    !Setting weights for integration
     endif

     ioffs = negf%Np_n(1)

     do i = 1, negf%Np_n(2)
        if (kbT.eq.0.0_dp .or. wqmax>0.0_dp) then   ! Circle integration T=0
           Pc = Rad*exp(j*pnts(i))
           Ec = Centre+Pc
           dt = negf%g_spin*wght(i)/(2.0_dp*pi)
           zt = dt*Pc*j
        else                                        ! Segment integration T>0
           Ec = z1 + pnts(i)*z_diff
           ff = fermi(Ec,muref,KbT)
           zt = negf%g_spin * z_diff * ff * wght(i) / (2.0_dp *pi)
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
     do i = 1, nPoles
        Ec = muref + j * KbT *pi* (2.0_dp*real(i,dp) - 1.0_dp)
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
     do i = 0, Npoints-1
        negf%en_grid(i+1)%cpu = mod(i,numprocs)
     enddo

  end subroutine contour_int_def

  !----------------------------------------------------------------------------------
  ! Performs a generic contour integration + eventual poles.
  ! Must be called after contour_int_def (see)
  !
  ! NOTE: For the moment it is implicitly assumed that contour_int is called BEFORE
  !       real-axis integration (any existing matrix is overwritten)
  !----------------------------------------------------------------------------------
  subroutine contour_int(negf)
     type(Tnegf) :: negf

     type(z_CSR) :: GreenR, TmpMt
     integer :: i, i1, ncont, Npoints, outer
     real(dp) :: ncyc, scba_error
     complex(dp) :: Ec, zt

     ncont = negf%str%num_conts
     outer = negf%outer_blocks
     Npoints = size(negf%en_grid)
     call create(TmpMt,negf%H%nrow,negf%H%ncol,negf%H%nrow)
     call initialize(TmpMt)

     call write_info_parallel(negf%verbose,30,'CONTOUR INTEGRAL',Npoints)

     do i = 1, Npoints

        call write_Epoint(negf%verbose,negf%en_grid(i), Npoints)
        if (negf%en_grid(i)%cpu .ne. id) cycle

        Ec = negf%en_grid(i)%Ec
        zt = negf%en_grid(i)%wght
        negf%iE = negf%en_grid(i)%pt
        negf%iE_path = negf%en_grid(i)%pt_path

        call compute_Gr(negf, outer, ncont, Ec, scba_error, GreenR)

        if(negf%DorE.eq.'E') zt = zt * Ec

        call concat(TmpMt,zt,GreenR,1,1)

        call destroy(GreenR)

     enddo

     ! NOTE: The following assumes contour_int is called BEFORE any
     !       real-axis integration.
     if(negf%DorE.eq.'D') then
        call zspectral(TmpMt,TmpMt,0,negf%rho)
     endif
     if(negf%DorE.eq.'E') then
        call zspectral(TmpMt,TmpMt,0,negf%rho_eps)
     endif

     call destroy(TmpMt)

  end subroutine contour_int

  !----------------------------------------------------------------------------------
  ! Performs a generic contour integration. Must be called after contour_int_def
  ! Currently it is similar to contour_int, with the k-sum performed internally
  ! in order to be consistent with real_axis_int_inel()
  ! Therefore at the moment the inelastic Sigma^r are not computed
  !
  ! NOTE: For the moment it is implicitly assumed that contour_int is called BEFORE
  !       real-axis integration (any existing matrix is overwritten)
  !----------------------------------------------------------------------------------
  subroutine contour_int_inel(negf)
     type(Tnegf) :: negf

     type(z_CSR) :: GreenR, TmpMt
     integer :: i, i1, ncont, Npoints, outer
     real(dp) :: ncyc, scba_error
     complex(dp) :: Ec, zt

     ncont = negf%str%num_conts
     outer = negf%outer_blocks
     Npoints = size(negf%en_grid)
     call create(TmpMt,negf%H%nrow,negf%H%ncol,negf%H%nrow)
     call initialize(TmpMt)

     call write_info_parallel(negf%verbose,30,'CONTOUR INTEGRAL',Npoints)

     ! Loop over local k-points
     kloop: do iK = 1, size(negf%local_k_index)

       negf%iKloc = iK
       negf%iKpoint = negf%local_k_index(iK) ! global k index
       negf%kwght = negf%kweights(negf%iKpoint)
       call write_kpoint(negf%verbose, iK, negf%iKpoint)

       ! Set the working pointers
       negf%H => negf%HS(iK)%H
       negf%S => negf%HS(iK)%S
       call destroy_contact_matrices(negf)
       call extract_cont(negf)


       enloop:do i = 1, Npoints

         call write_Epoint(negf%verbose,negf%en_grid(i), Npoints)
         if (negf%en_grid(i)%cpu .ne. id) cycle

         Ec = negf%en_grid(i)%Ec
         zt = negf%en_grid(i)%wght
         negf%iE = negf%en_grid(i)%pt
         negf%iE_path = negf%en_grid(i)%pt_path

         call compute_Gr(negf, outer, ncont, Ec, scba_error, GreenR)

         if(negf%DorE.eq.'E') zt = zt * Ec

         call concat(TmpMt,zt,GreenR,1,1)

         call destroy(GreenR)

       enddo enloop

     enddo kloop

     ! NOTE: The following assumes contour_int is called BEFORE any
     !       real-axis integration.
     if(negf%DorE.eq.'D') then
        call zspectral(TmpMt,TmpMt,0,negf%rho)
#:if defined("MPI")
      if (id0.and.negf%verbose.gt.VBT) call message_clock('Gather MPI results ')
      call mpifx_allreduceip(negf%energyComm, negf%rho%nzval, MPI_SUM)
      call mpifx_allreduceip(negf%kComm, negf%rho%nzval, MPI_SUM)
      if (id0.and.negf%verbose.gt.VBT) call write_clock
#:endif
     endif

     if(negf%DorE.eq.'E') then
        call zspectral(TmpMt,TmpMt,0,negf%rho_eps)
#:if defined("MPI")
      if (id0.and.negf%verbose.gt.VBT) call message_clock('Gather MPI results ')
      call mpifx_allreduceip(negf%energyComm, negf%rho_eps%nzval, MPI_SUM)
      call mpifx_allreduceip(negf%kComm, negf%rho_eps%nzval, MPI_SUM)
      if (id0.and.negf%verbose.gt.VBT) call write_clock
#:endif
     endif

     call destroy(TmpMt)

  end subroutine contour_int_inel

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

    integer :: i, i1, np, ioffset, ncont, Npoints
    real(dp), DIMENSION(:), allocatable :: wght, pnts
    real(dp) :: Omega, mumin, mumax, wqmax
    integer, allocatable :: seq(:)

    ncont = negf%str%num_conts
    Npoints = negf%Np_real

    ! With inelastic scattering equalize the energy points/process
    if (negf%interactList%counter /= 0) then
      if (mod(Npoints,numprocs) .ne. 0) then
        do while (mod(Npoints,numprocs) .ne. 0)
          Npoints = Npoints + 1
        end do
        negf%Np_real = Npoints
      end if
    end if

    ! Get maximum wq for inelastic interactions
    wqmax = get_max_wq(negf%interactList)

    ! Omega considers maximum kT so interval is always large enough
    Omega = negf%n_kt * maxval(negf%cont(:)%kbT_dm) + wqmax

    if (ncont.gt.0) then
       mumin=minval(negf%cont(:)%mu)
       mumax=maxval(negf%cont(:)%mu)
    else
       mumin=negf%muref
       mumax=negf%muref
    endif

    ! Balance N points over the energy grid
    if (mod(Npoints,numprocs) .ne. 0) then
      do while (mod(Npoints,numprocs) .ne. 0)
        Npoints = Npoints + 1
      end do
    end if

    !! destroy en_grid from previous calculation, if any
    call destroy_en_grid(negf%en_grid)
    allocate(negf%en_grid(Npoints))

    allocate(pnts(Npoints))
    allocate(wght(Npoints))

    if (wqmax == 0.0_dp) then
      call gauleg(mumin-Omega,mumax+Omega,pnts,wght,Npoints)
    else
      ! Trapezoidal is used as it gives a regular grid for inelastic S.E.
      call trapez(mumin-Omega,mumax+Omega,pnts,wght,Npoints)
    end if

    ! Offset due to the contour integration (for contact storage)
    ioffset = negf%Np_n(1) + negf%Np_n(2) + negf%n_poles
    do i = 1, Npoints
       negf%en_grid(i)%path = 1
       negf%en_grid(i)%pt = i
       negf%en_grid(i)%pt_path = ioffset + i
       negf%en_grid(i)%Ec = cmplx(pnts(i),negf%delta,dp)
       negf%en_grid(i)%wght = negf%g_spin * wght(i)/(2.0_dp *pi)
    enddo

    deallocate(wght)
    deallocate(pnts)

    np=Npoints/numprocs
    if (id .ne. numprocs-1) then
      negf%local_en_points = np
    else
      negf%local_en_points = np + mod(Npoints,numprocs)
    end if

    if (wqmax == 0.0_dp) then
      ! distribute energy grid round robin scheme:
      !     0 1 2 3 ... 0 1 2 .... 0 1 2 ....
      do i = 0, Npoints-1
        negf%en_grid(i+1)%cpu = mod(i,numprocs)
      enddo
    else
      ! With inelastic points must be contigous
      !    0 0 0 0 ... 1 1 1 .... 2 2 2 ....
      allocate(seq(np))
      do i1 = 1, np
        seq(i1) = i1
      end do
      do i = 0, numprocs-1
        negf%en_grid(i*np+1:(i+1)*np)%cpu = i
        negf%en_grid(i*np+1:(i+1)*np)%pt_cpu = seq
      end do
    end if

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

    type(z_CSR) :: G, TmpMt

    integer :: ref, Npoints
    integer :: iE, i1, j1, iK, outer, ncont
    integer :: particle

    real(dp), DIMENSION(:), allocatable :: frm_f
    complex(dp), DIMENSION(:), allocatable :: diag
    real(dp) :: ncyc, Er, scba_error, scba_elastic_tol

    complex(dp) :: zt
    complex(dp) :: Ec

    ncont = negf%str%num_conts
    outer = negf%outer_blocks
    Npoints = size(negf%en_grid)
    particle = negf%particle

    call log_allocate(frm_f,ncont)

    call create(TmpMt,negf%H%nrow,negf%H%ncol,negf%H%nrow)
    call initialize(TmpMt)

    if (negf%cartComm%rank == 0) then
      call write_info_parallel(negf%verbose,30,'REAL AXIS INTEGRAL',Npoints)
    end if

    !! Loop over energy points
    enloop: do iE = 1, Npoints

      call write_Epoint(negf%verbose,negf%en_grid(iE),Npoints)
      if (negf%en_grid(iE)%cpu .ne. id) cycle
      Ec = negf%en_grid(iE)%Ec
      Er = real(Ec)
      negf%iE = negf%en_grid(iE)%pt
      negf%iE_path = negf%en_grid(iE)%pt_path

      frm_f = 0.0_dp
      do j1 = 1,ncont
        frm_f(j1)=fermi(Er,negf%cont(j1)%mu,negf%cont(j1)%kbT_dm)
      enddo

      if (particle == 1) then
        call compute_Gn(negf, outer, ncont, Ec, frm_f, G, scba_error)
      else
        call compute_Gp(negf, outer, ncont, Ec, frm_f, G, scba_error)
      endif

      zt = negf%en_grid(iE)%wght !!* negf%kwght
      if (negf%DorE.eq.'E') then
         zt = zt * Er
      end if

      call concat(TmpMt,zt,G,1,1)

      call destroy(G)

    enddo enloop

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

  !-----------------------------------------------------------------------
  ! Real axis integral with inelastic scattering
  !-----------------------------------------------------------------------
  !    g     /
  ! ------   | Gn(E) dE
  ! 2*pi*j   /
  !-----------------------------------------------------------------------
  subroutine real_axis_int_inel(negf)
    type(Tnegf) :: negf

    integer :: ref, Npoints
    integer :: scba_iter, scba_elastic_iter, scba_niter_ela, scba_niter_inela
    integer :: iE, i1, j1, iK, outer, ncont, ref_bk
    integer :: particle
    type(z_CSR) :: G, TmpMt

    real(dp), DIMENSION(:), allocatable :: frm_f
    complex(dp), DIMENSION(:), allocatable :: diag
    real(dp) :: ncyc, Er, scba_elastic_error, scba_inelastic_error
    complex(dp) :: zt, Ec

    ncont = negf%str%num_conts
    outer = negf%outer_blocks
    Npoints = size(negf%en_grid)
    particle = negf%particle

    ref_bk = negf%refcont
    negf%refcont = ncont + 1
    call log_allocate(frm_f,ncont+1)
    frm_f = 0.0_dp

    call create(TmpMt,negf%H%nrow,negf%H%ncol,negf%H%nrow)
    call initialize(TmpMt)

    if (negf%cartComm%rank == 0) then
      call write_info_parallel(negf%verbose,30,'REAL AXIS INTEGRAL (inel)',Npoints)
    end if
    ! ---------------------------------------------------------------------
    ! SCBA Iteration
    ! ---------------------------------------------------------------------
    call interaction_prepare(negf)
    call negf%scbaDriverInelastic%init(negf%scba_inelastic_tol, dowrite=.false.)
    scba_niter_inela = get_max_niter_inelastic(negf%interactList)
    scba_niter_ela = get_max_niter_elastic(negf%interactList)
    if (negf%cartComm%rank == 0) then
      call write_info(negf%verbose, 30, 'NUMBER OF SCBA INELASTIC ITERATIONS', scba_niter_inela)
      call write_info(negf%verbose, 30, 'NUMBER OF SCBA ELASTIC ITERATIONS', scba_niter_ela)
    end if
    scba_iter = 0

    scba: do while (.not.negf%scbaDriverInelastic%is_converged() .and. scba_iter <= scba_niter_inela)

      if (negf%cartComm%rank == 0) then
        call write_info(negf%verbose, 30, ' INELASTIC SCBA ITERATION', scba_iter)
      end if

      call negf%scbaDriverInelastic%set_scba_iter(scba_iter, negf%interactList)

      ! Density matrix initialization
      TmpMt%nzval = (0.0_dp, 0.0_dp)

      ! Loop over local k-points
      kloop: do iK = 1, size(negf%local_k_index)

        negf%iKloc = iK
        negf%iKpoint = negf%local_k_index(iK) ! global k index
        negf%kwght = negf%kweights(negf%iKpoint)
        call write_kpoint(negf%verbose, iK, negf%iKpoint)

        ! Set the working pointers
        negf%H => negf%HS(iK)%H
        negf%S => negf%HS(iK)%S
        call destroy_contact_matrices(negf)
        call extract_cont(negf)

        !! Loop over energy points
        enloop: do iE = 1, Npoints

          call write_Epoint(negf%verbose,negf%en_grid(iE),Npoints)
          if (negf%en_grid(iE)%cpu .ne. id) cycle
          Ec = negf%en_grid(iE)%Ec
          Er = real(Ec)
          negf%iEloc = negf%en_grid(iE)%pt_cpu
          negf%iE = negf%en_grid(iE)%pt
          negf%iE_path = negf%en_grid(iE)%pt_path

          frm_f = 0.0_dp
          do j1 = 1,ncont
            frm_f(j1)=fermi(Er,negf%cont(j1)%mu,negf%cont(j1)%kbT_dm)
          enddo

          if (particle == 1) then
            call compute_Gn(negf, outer, ncont, Ec, frm_f, G, scba_elastic_error)
          else
            call compute_Gp(negf, outer, ncont, Ec, frm_f, G, scba_elastic_error)
          endif

          zt = negf%en_grid(iE)%wght * negf%kwght
          if (negf%DorE.eq.'E') then
             zt = zt * Er
          end if

          call concat(TmpMt,zt,G,1,1)

          call destroy(G)

        enddo enloop

      enddo kloop

      ! ---------------------------------------------------------------------
      ! COMPUTE SELF-ENERGIES
      ! ---------------------------------------------------------------------
      if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute self-energies ')
      call compute_Sigmas_inelastic(negf)
      if (id0.and.negf%verbose.gt.VBT) call write_clock

      !Check SCBA convergence on layer currents.
      call negf%scbaDriverInelastic%check_Mat_convergence(TmpMt)

      if (negf%cartComm%rank == 0) then
        scba_inelastic_error = negf%scbaDriverInelastic%scba_error()
        call write_real_info(negf%verbose, 30, 'SCBA inelastic error',scba_inelastic_error)
      end if

#:if defined("MPI")
      !In MPI runs only root has a meaningful result => Bcast the result
      call mpifx_bcast(negf%cartComm, negf%scbaDriverInelastic%converged)
#:endif

      scba_iter = scba_iter + 1

    end do scba

    if(negf%DorE.eq.'D') then
       if(allocated(negf%rho%nzval)) then
          call concat(negf%rho,TmpMt,1,1)
       else
          call clone(TmpMt,negf%rho)
       endif
       call destroy(TmpMt)
#:if defined("MPI")
      if (id0.and.negf%verbose.gt.VBT) call message_clock('Gather MPI results ')
      call mpifx_allreduceip(negf%energyComm, negf%rho%nzval, MPI_SUM)
      call mpifx_allreduceip(negf%kComm, negf%rho%nzval, MPI_SUM)
      if (id0.and.negf%verbose.gt.VBT) call write_clock
#:endif
    endif

    if(negf%DorE.eq.'E') then
       if(allocated(negf%rho_eps%nzval)) then
          call concat(negf%rho_eps,TmpMt,1,1)
       else
          call clone(TmpMt,negf%rho_eps)
       endif
       call destroy(TmpMt)
#:if defined("MPI")
      if (id0.and.negf%verbose.gt.VBT) call message_clock('Gather MPI results ')
      call mpifx_allreduceip(negf%energyComm, negf%rho_eps%nzval, MPI_SUM)
      call mpifx_allreduceip(negf%kComm, negf%rho_eps%nzval, MPI_SUM)
      if (id0.and.negf%verbose.gt.VBT) call write_clock
#:endif
    endif

    call negf%scbaDriverInelastic%destroy()
    call log_deallocate(frm_f)

  end subroutine real_axis_int_inel


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

    integer :: i, ioffset, ncont, Npoints
    real(dp), DIMENSION(:), allocatable :: wght,pnts   ! Gauss-quadrature points
    real(dp) :: Omega, mumax, muref, ff, kbT

    ncont = negf%str%num_conts
    ioffset = negf%Np_n(1) + negf%Np_n(2) + negf%n_poles
    kbT = maxval(negf%cont(:)%kbT_dm)
    Omega = negf%n_kt * kbT
    muref = negf%muref

    if (ncont.gt.0) then
       mumax=maxval(negf%cont(:)%mu_n)
    else
       mumax=negf%muref
    endif

    Npoints = negf%Np_real
    call destroy_en_grid(negf%en_grid)
    allocate(negf%en_grid(Npoints))

    allocate(pnts(Npoints))
    allocate(wght(Npoints))

    call trapez(negf%Ec-negf%DeltaEc, mumax + Omega,pnts,wght,Npoints)

    do i = 1, Npoints
       negf%en_grid(i)%path = 1
       negf%en_grid(i)%pt = i
       negf%en_grid(i)%pt_path = ioffset + i
       negf%en_grid(i)%Ec = cmplx(pnts(i),negf%delta,dp)

       !IMPORTANT: there is no fermi function multiplied in negf%en_grid(i)%wght anymore,
       !so we cannot use this routine in combination with contour_int
       negf%en_grid(i)%wght = negf%g_spin * wght(i) / (2.d0 * pi)
    enddo

    deallocate(wght)
    deallocate(pnts)

    ! distribute energy grid
    do i = 0, Npoints-1
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

    integer :: i, ioffset, ncont, Npoints
    real(dp), DIMENSION(:), allocatable :: wght,pnts   ! Gauss-quadrature points
    real(dp) :: Omega, mumin, Emax, muref, ff, kbT

    ncont = negf%str%num_conts
    ioffset = negf%Np_n(1) + negf%Np_n(2) + negf%n_poles
    kbT = maxval(negf%cont(:)%kbT_dm)
    Omega = negf%n_kt * kbT
    muref = negf%muref

    if (ncont.gt.0) then
       mumin=minval(negf%cont(:)%mu_p)
    else
       mumin=negf%muref
    endif

    Npoints = negf%Np_real
    call destroy_en_grid(negf%en_grid)
    allocate(negf%en_grid(Npoints))

    allocate(pnts(Npoints))
    allocate(wght(Npoints))

    Emax = negf%Ev+negf%DeltaEv

    call trapez(mumin-Omega, Emax, pnts, wght, Npoints)

    do i = 1, Npoints
       negf%en_grid(i)%path = 1
       negf%en_grid(i)%pt = i
       negf%en_grid(i)%pt_path = ioffset + i
       negf%en_grid(i)%Ec = cmplx(pnts(i),negf%delta,dp)

       !IMPORTANT: there is no fermi function multiplied in negf%en_grid(i)%wght anymore,
       !so we cannot use this routine in combination with contour_int
       negf%en_grid(i)%wght = negf%g_spin * wght(i) / (2.d0 * pi)
    enddo

    deallocate(wght)
    deallocate(pnts)

    ! distribute energy grid
    do i = 0, Npoints-1
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

    xm=0.5_dp*(x2+x1)
    xl=0.5_dp*(x2-x1)

    do i=1,m

       ! Approssimazione degli zeri dei polinomi di Legendre:
       z=cos(Pi*(i-0.25_dp)/(n+0.5_dp))

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
          p0=1.0_dp  !p(0)
          p1=0.0_dp  !p(-1)

          ! Legendre polynomial p1 evaluated by rec. relations:
          do k=1,n
             p2=p1 !p(-2)=p(-1)
             p1=p0 !p(-1)=p(0)
             p0=((2.0_dp*k-1.0_dp)*z*p1-(k-1.0_dp)*p2)/k
          enddo

          pp=n*(z*p0-p1)/(z*z-1.0_dp)

          ! Newton method to refine the zeros:
          z1=z
          z=z1-p0/pp

          if(abs(z-z1).le.ACC) exit
       enddo

       ! Scale the interval to x1..x2:
       x(i)=xm-xl*z
       x(n+1-i)=xm+xl*z
       w(i)=2.0_dp*xl/((1.0_dp-z*z)*pp*pp)
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

    integer :: i, k, ncont, Npoints, np
    integer, allocatable :: seq(:)

    Npoints=nint((negf%Emax-negf%Emin)/negf%Estep) + 1

    ! With inelastic scattering equalize the energy points/process
    if (negf%interactList%counter /= 0) then
      if (mod(Npoints,numprocs) .ne. 0) then
        do while (mod(Npoints,numprocs) .ne. 0)
          Npoints = Npoints + 1
        end do
        negf%Emax = negf%Emin + (Npoints-1) * negf%Estep
      end if
    end if

    call destroy_en_grid(negf%en_grid)
    allocate(negf%en_grid(Npoints))

    do i = 1, Npoints
       negf%en_grid(i)%path = 1
       negf%en_grid(i)%pt = i
       negf%en_grid(i)%pt_path = i
       negf%en_grid(i)%Ec = cmplx(negf%Emin + negf%Estep*(i-1), 0.0, dp)
       negf%en_grid(i)%wght = negf%kwght
    enddo

    np=Npoints/numprocs
    if (id .ne. numprocs-1) then
      negf%local_en_points = np
    else
      negf%local_en_points = np + mod(Npoints,numprocs)
    end if

    if ( get_max_wq(negf%interactList) == 0.0_dp) then
      ! distribute energy grid round robin scheme
      !     0 1 2 3 ... 0 1 2 .... 0 1 2 ....
      do i = 0, Npoints-1
        negf%en_grid(i+1)%cpu = mod(i,numprocs)
      enddo
    else
      ! With inelastic points must be contigous
      !    0 0 0 0 ... 1 1 1 .... 2 2 2 ....
      allocate(seq(np))
      do k = 1, np
        seq(k) = k
      end do
      do i = 0, numprocs-1
        negf%en_grid(i*np+1:(i+1)*np)%cpu = i
        negf%en_grid(i*np+1:(i+1)*np)%pt_cpu = seq
      end do
    end if

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
    Real(dp), Dimension(:), allocatable :: tun_mat
    Real(dp), Dimension(:), allocatable :: ledos
    Real(dp) :: mu1, mu2   ! contact potentials
    Real(dp) :: ncyc       ! stores average number of iters in decimation

    Integer :: i, icont, icpl      ! dummy counters
    Integer :: ncont               ! number of contacts
    Integer :: size_ni             ! emitter-collector contacts

    Integer :: Npoints               ! number of integration points
    Complex(dp) :: Ec              ! Energy point

    Logical :: do_ledos            ! performs or not dos_proj


    ! Get out immediately if Emax<Emin
    if (negf%Emax.le.negf%Emin) then
       if(id0) write(*,*) '0 tunneling points;  current = 0.0'
       if (allocated(negf%currents)) call log_deallocate(negf%currents)
       call log_allocate(negf%currents,1)
       negf%currents = 0.0_dp
       return
    endif

    !-------------------------------------------------------

    do_ledos = .false.
    if(negf%ndos_proj.gt.0) do_ledos=.true.
    ncont = negf%str%num_conts
    Npoints = size(negf%en_grid)
    ncyc=0
    size_ni = size(negf%ni)
    negf%readOldSGF = negf%readOldT_SGFs

    !-------------------------------------------------------

    call log_allocate(tun_mat,size_ni)
    !If previous calculation is there, destroy output
    if (allocated(negf%tunn_mat)) then
       call log_deallocate(negf%tunn_mat)
    end if
    call log_allocate(negf%tunn_mat,Npoints,size_ni)

    negf%tunn_mat = 0.0_dp

    if (do_ledos) then
       !If previous calculation is there, destroy output
       if (allocated(negf%ldos_mat)) then
         call log_deallocate(negf%ldos_mat)
       end if
       call log_allocate(negf%ldos_mat,Npoints,negf%ndos_proj)
       call log_allocate(ledos,negf%ndos_proj)
       negf%ldos_mat(:,:)=0.0_dp
    endif

    !-------------------------------------------------------
    call write_info_parallel(negf%verbose,30,'COHERENT TRANSMISSION; Npoints',Npoints)

    !Loop on energy points: tunneling
    do i = 1, Npoints

       call write_Epoint(negf%verbose,negf%en_grid(i), size(negf%en_grid))
       if (negf%en_grid(i)%cpu /= id) cycle

       Ec = negf%en_grid(i)%Ec
       negf%iE = negf%en_grid(i)%pt
       negf%iE_path = negf%en_grid(i)%pt_path

       if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Contact SE ')
       call compute_contacts(Ec+j*negf%delta,negf,ncyc,Tlc,Tcl,SelfEneR,GS)
       if (id0.and.negf%verbose.gt.VBT) call write_clock
       call write_int_info(negf%verbose, VBT, 'Average number of iterations', int(ncyc))

       if (.not.do_ledos) then
          if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Tunneling ')

          call calculate_transmissions(negf, Ec, SelfEneR, negf%tun_proj, tun_mat)

          negf%tunn_mat(i,:) = tun_mat(:) * negf%kwght
       else
          if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Tunneling and DOS')
          ledos(:) = 0.0_dp

          call calculate_transmissions_and_dos(negf, Ec, SelfEneR, GS, negf%tun_proj, tun_mat, &
                                             &  negf%dos_proj, ledos)

          negf%tunn_mat(i,:) = tun_mat(:) * negf%kwght
          negf%ldos_mat(i,:) = ledos(:) * negf%kwght
       endif

       if (id0.and.negf%verbose.gt.VBT) call write_clock

       do icont=1,ncont
         call destroy(Tlc(icont),Tcl(icont),SelfEneR(icont),GS(icont))
       enddo

    enddo !Loop on energy

    !call destroy_en_grid()
    call log_deallocate(tun_mat)
    if (do_ledos) then
       call log_deallocate(ledos)
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

    integer :: scba_iter, scba_niter_inela, scba_niter_ela, size_ni, ref_bk
    integer :: iE, iK, jj, i1, j1, Npoints, outer, ncont, icont, scba_elastic_iter
    real(dp) :: ncyc, scba_elastic_tol, scba_elastic_error, scba_inelastic_error
    Type(z_DNS), dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    real(dp), dimension(:), allocatable :: curr_mat, frm
    complex(dp) :: Ec

    ncont = negf%str%num_conts
    Npoints = size(negf%en_grid)
    size_ni = size(negf%ni)
    negf%readOldSGF = negf%readOldT_SGFs
    outer = 0

    if (.not. allocated(negf%curr_mat)) then
      call log_allocate(negf%curr_mat,Npoints,ncont)
    end if
    call log_allocate(curr_mat, ncont)

    ! Clean up caches of G_r and G_n
    call negf%G_r%destroy()
    call negf%G_n%destroy()

    ! Create Fermi array. Set reference such that f(ref)=0.
    ref_bk = negf%refcont
    negf%refcont = ncont + 1
    call log_allocate(frm, ncont+1)
    frm = 0.0_dp

    ! Fixed occupations (e.g. 1.0, 0.0) can be used to get a transmission
    if (present(fixed_occupations)) then
      frm(1:ncont) = fixed_occupations
    end if

    call write_info_parallel(negf%verbose,30,'MEIR-WINGREEN FORMULA; Npoints',Npoints)
    ! ---------------------------------------------------------------------
    ! SCBA Iteration
    ! ---------------------------------------------------------------------
    call interaction_prepare(negf)
    call negf%scbaDriverInelastic%init(tol = negf%scba_inelastic_tol, dowrite = .true.)
    scba_niter_inela = get_max_niter_inelastic(negf%interactList)
    scba_niter_ela = get_max_niter_elastic(negf%interactList)
    scba_elastic_tol = negf%scba_elastic_tol
    call write_info(negf%verbose, 30, 'NUMBER OF SCBA INELASTIC ITERATIONS', scba_niter_inela)
    call write_info(negf%verbose, 30, 'NUMBER OF SCBA ELASTIC ITERATIONS', scba_niter_ela)
    scba_iter = 0

    scba: do while (.not.negf%scbaDriverInelastic%is_converged() .and. scba_iter <= scba_niter_inela)

      if (negf%cartComm%rank == 0) then
        call write_info(negf%verbose, 30, ' INELASTIC SCBA ITERATION', scba_iter)
      end if

      call negf%scbaDriverInelastic%set_scba_iter(scba_iter, negf%interactList)

      negf%curr_mat = 0.0_dp
      ! Loop over local k-points
      kloop: do iK = 1, size(negf%local_k_index)

        negf%iKloc = iK
        negf%iKpoint = negf%local_k_index(iK) ! global k index
        negf%kwght = negf%kweights(negf%iKpoint)
        call write_kpoint(negf%verbose, iK, negf%iKpoint)

        ! Set the working pointers
        negf%H => negf%HS(iK)%H
        negf%S => negf%HS(iK)%S
        call extract_cont(negf)

        !! Loop over energy points
        enloop: do iE = 1, Npoints

          call write_Epoint(negf%verbose, negf%en_grid(iE), size(negf%en_grid))
          if (negf%en_grid(iE)%cpu /= id) cycle
          Ec = negf%en_grid(iE)%Ec
          negf%iEloc = negf%en_grid(iE)%pt_cpu
          negf%iE = negf%en_grid(iE)%pt
          negf%iE_path = negf%en_grid(iE)%pt_path

          if (.not.present(fixed_occupations)) then
            do j1 = 1,ncont
              frm(j1)=fermi(real(Ec), negf%cont(j1)%mu, negf%cont(j1)%kbT_t)
            enddo
          end if

          ! ---------------------------------------------------------------------
          ! Compute contact GF
          ! ---------------------------------------------------------------------
          if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Contact SE ')
          call compute_contacts(Ec+j*negf%delta, negf, ncyc, Tlc, Tcl, SelfEneR, GS)
          if (id0.and.negf%verbose.gt.VBT) call write_clock
          call write_int_info(negf%verbose, VBT, 'Average number of iterations', int(ncyc))

          ! ---------------------------------------------------------------------
          ! Compute block tri-diagonal Gr and Gn
          ! ---------------------------------------------------------------------
          negf%tDestroyGr = .false.; negf%tDestroyGn = .false.
          if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute elastic SCBA ')
          call calculate_elastic_scba(negf,real(Ec),SelfEneR,Tlc,Tcl,GS,frm,scba_niter_ela, &
                scba_elastic_tol, scba_elastic_iter, scba_elastic_error)
          if (id0.and.negf%verbose.gt.VBT) call write_clock

          !call write_real_info(negf%verbose, VBT, 'scba elastic error',scba_elastic_error)

          negf%tDestroyGr = .true.; negf%tDestroyGn = .true.
          if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Meir-Wingreen ')
          call iterative_meir_wingreen(negf,real(Ec),SelfEneR,frm,curr_mat)
          if (id0.and.negf%verbose.gt.VBT) call write_clock

          ! Recursive sum adds up k-dependent partial results
          negf%curr_mat(iE,:) = negf%curr_mat(iE,:) + curr_mat(:) * negf%kwght

          do icont=1,ncont
            call destroy(Tlc(icont),Tcl(icont),SelfEneR(icont),GS(icont))
          enddo

        end do enloop

        call destroy_contact_matrices(negf)

      end do kloop

      ! ---------------------------------------------------------------------
      ! COMPUTE SELF-ENERGIES
      ! ---------------------------------------------------------------------
      if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute self-energies ')
      call compute_Sigmas_inelastic(negf)
      if (id0.and.negf%verbose.gt.VBT) call write_clock

      call electron_current_meir_wingreen(negf)

#:if defined("MPI")
      if (id0.and.negf%verbose.gt.VBT) call message_clock('Gather MPI results ')
      call mpifx_reduceip(negf%energyComm, negf%currents, MPI_SUM)
      call mpifx_reduceip(negf%kComm, negf%currents, MPI_SUM)
      if (id0.and.negf%verbose.gt.VBT) call write_clock
#:endif

      !Check SCBA convergence on layer currents.
      !In MPI runs only root has a meaningful result => Bcast the result
      call negf%scbaDriverInelastic%check_J_convergence(negf%currents)
      if (negf%cartComm%rank == 0) then
        scba_inelastic_error = negf%scbaDriverInelastic%scba_error()
        call write_real_info(negf%verbose, VBT, 'scba inelastic error',scba_inelastic_error)
      end if
#:if defined("MPI")
      call mpifx_bcast(negf%cartComm, negf%scbaDriverInelastic%converged)
#:endif

      scba_iter = scba_iter + 1

    end do scba

    !call interaction_cleanup(negf)
    call negf%scbaDriverInelastic%destroy()
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

    integer :: nbl, Npoints, outer
    integer :: iE, i1, j1, iK, icont, ncont, ref_bk
    integer :: scba_iter, scba_elastic_iter, scba_niter_ela
    integer :: scba_inelastic_iter, scba_niter_inela
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    real(dp) :: Er, ncyc, scba_elastic_tol, scba_elastic_error
    !real(dp) :: scba_elastic_error, scba_inelastic_error, scba_elastic_tol
    real(dp) :: scba_inelastic_error
    real(dp), dimension(:), allocatable :: curr_mat, ldos_mat, frm
    complex(dp), DIMENSION(:), allocatable :: diag
    complex(dp) :: Ec, zt
    type(z_CSR) :: G, TmpMt

    ncont = negf%str%num_conts
    Npoints = size(negf%en_grid)
    nbl = negf%str%num_PLs
    negf%readOldSGF = negf%readOldT_SGFs
    outer = negf%outer_blocks

    ! Allocating curr_mat to nbl-1 so we compute L->L+1 currents
    if (.not. allocated(negf%curr_mat)) then
      call log_allocate(negf%curr_mat, Npoints, nbl-1)
    end if
    call log_allocate(curr_mat, nbl-1)

    ! Allocating ldos_mat to nbl
    if (.not. allocated(negf%ldos_mat)) then
      call log_allocate(negf%ldos_mat, Npoints, nbl)
    end if
    call log_allocate(ldos_mat, nbl)

    call create(TmpMt,negf%H%nrow,negf%H%ncol,negf%H%nrow)
    call initialize(TmpMt)

    ! Clean up caches of G_r and G_n
    call negf%G_r%destroy()
    call negf%G_n%destroy()

    ! Create Fermi array. Set reference such that f(ref)=0.
    ref_bk = negf%refcont
    negf%refcont = ncont + 1
    call log_allocate(frm, ncont+1)
    frm = 0.0_dp

    if (negf%cartComm%rank == 0) then
      call write_info_parallel(negf%verbose,30,'CALCULATION OF LAYER CURRENTS; Npoints=',Npoints)
    end if

    !if (id0) then
    !  do iK = 1, size(negf%local_k_index)
    !    negf%iKpoint = negf%local_k_index(iK) ! global k index
    !    print*,'cpu',negf%kComm%rank,':  k=',negf%kpoints(:,negf%iKpoint)
    !  end do
    !end if
    ! ---------------------------------------------------------------------
    ! SCBA Iteration
    ! ---------------------------------------------------------------------
    call interaction_prepare(negf)
    call negf%scbaDriverInelastic%init(negf%scba_inelastic_tol, dowrite=.false.)
    scba_niter_inela = get_max_niter_inelastic(negf%interactList)
    scba_niter_ela = get_max_niter_elastic(negf%interactList)
    scba_elastic_tol = negf%scba_elastic_tol
    if (negf%cartComm%rank == 0) then
      call write_info(negf%verbose, 30, 'NUMBER OF SCBA INELASTIC ITERATIONS', scba_niter_inela)
      call write_info(negf%verbose, 30, 'NUMBER OF SCBA ELASTIC ITERATIONS', scba_niter_ela)
    end if
    scba_iter = 0

    scba: do while (.not.negf%scbaDriverInelastic%is_converged() .and. scba_iter <= scba_niter_inela)

      if (negf%cartComm%rank == 0) then
        call write_info(negf%verbose, 30, ' INELASTIC SCBA ITERATION', scba_iter)
      end if

      call negf%scbaDriverInelastic%set_scba_iter(scba_iter, negf%interactList)

      ! Density matrix initialization
      TmpMt%nzval = (0.0_dp, 0.0_dp)
      negf%curr_mat = 0.0_dp
      negf%ldos_mat = 0.0_dp
      ! Loop over local k-points
      kloop: do iK = 1, size(negf%local_k_index)

        negf%iKloc = iK
        negf%iKpoint = negf%local_k_index(iK) ! global k index
        negf%kwght = negf%kweights(negf%iKpoint)
        call write_kpoint(negf%verbose, iK, negf%iKpoint)

        negf%H => negf%HS(iK)%H
        negf%S => negf%HS(iK)%S
        call extract_cont(negf)

        !! Loop over energy points
        enloop: do iE = 1, Npoints

          if (negf%en_grid(iE)%cpu /= id) cycle
          call write_Epoint(negf%verbose, negf%en_grid(iE), size(negf%en_grid))
          Ec = negf%en_grid(iE)%Ec + j*negf%delta
          Er = real(Ec)
          negf%iEloc = negf%en_grid(iE)%pt_cpu
          negf%iE = negf%en_grid(iE)%pt  ! global energy index
          negf%iE_path = negf%en_grid(iE)%pt_path

          frm = 0.0_dp
          do j1 = 1,ncont
             frm(j1)=fermi(real(Ec), negf%cont(j1)%mu, negf%cont(j1)%kbT_t)
          enddo

          ! Avoids cleanup of block-dense Gn for layer currents (ESH is recomputed)
          negf%tDestroyGn = .false.

          call compute_Gn(negf, outer, ncont, Ec, frm, G, scba_elastic_error)

          ! ---------------------------------------------------------------------
          ! Compute layer-to-layer currents and release memory
          ! ---------------------------------------------------------------------
          negf%tDestroyGn = .true.
          if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Jn,n+1 ')
          call iterative_layer_current(negf, Er, curr_mat, ldos_mat)
          if (id0.and.negf%verbose.gt.VBT) call write_clock

          ! Recursive sum adds up k-dependent partial results
          negf%curr_mat(iE,:) = negf%curr_mat(iE,:) + curr_mat(:) * negf%kwght
          negf%ldos_mat(iE,:) = negf%ldos_mat(iE,:) + ldos_mat(:) * negf%kwght

          call destroy(G)

        end do enloop

        call destroy_contact_matrices(negf)

      end do kloop

      ! ---------------------------------------------------------------------
      ! COMPUTE SELF-ENERGIES
      ! ---------------------------------------------------------------------
      if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute self-energies ')
      call compute_Sigmas_inelastic(negf)
      if (id0.and.negf%verbose.gt.VBT) call write_clock

      call electron_current_meir_wingreen(negf)

#:if defined("MPI")
      if (id0.and.negf%verbose.gt.VBT) call message_clock('Gather MPI results ')
      call mpifx_reduceip(negf%energyComm, negf%currents, MPI_SUM)
      call mpifx_reduceip(negf%kComm, negf%currents, MPI_SUM)
      call mpifx_reduceip(negf%energyComm, negf%ldos_mat, MPI_SUM)
      call mpifx_reduceip(negf%kComm, negf%ldos_mat, MPI_SUM)
      if (id0.and.negf%verbose.gt.VBT) call write_clock
#:endif

      !Check SCBA convergence on layer currents.
      call negf%scbaDriverInelastic%check_J_convergence(negf%currents)

      if (negf%cartComm%rank == 0) then
        scba_inelastic_error = negf%scbaDriverInelastic%scba_error()
        call write_real_info(negf%verbose, 30, 'SCBA inelastic error',scba_inelastic_error)
      end if

#:if defined("MPI")
      !In MPI runs only root has a meaningful result => Bcast the result
      call mpifx_bcast(negf%cartComm, negf%scbaDriverInelastic%converged)
#:endif

      scba_iter = scba_iter + 1

    end do scba

    !call interaction_cleanup(negf)
    call negf%scbaDriverInelastic%destroy()
    call log_deallocate(curr_mat)
    call log_deallocate(ldos_mat)
    call log_deallocate(frm)
    negf%refcont = ref_bk

  end subroutine layer_current


  !---------------------------------------------------------------------------
  subroutine interaction_prepare(negf)
    type(TNegf) :: negf

    real(dp) :: deltaE
    type(TInteractionNode), pointer :: it
    it => negf%interactList%first

    do while (associated(it))
      select type(pInter => it%inter)
      class is(ElPhonInel)
        deltaE = real(negf%en_grid(2)%Ec - negf%en_grid(1)%Ec)
        call pInter%set_EnGrid(deltaE, size(negf%en_grid), negf%local_en_points)
        if (allocated(negf%equivalent_kpoints)) then
          call pInter%set_kpoints(negf%kpoints, negf%kweights, negf%local_k_index, &
                               & negf%equivalent_kpoints)
        else
          call pInter%set_kpoints(negf%kpoints, negf%kweights, negf%local_k_index)
        end if
        call pInter%prepare()
      end select
      it => it%next
    end do

  end subroutine interaction_prepare

  !---------------------------------------------------------------------------
  !subroutine interaction_cleanup(negf)
  !  type(TNegf) :: negf
  !
  !  type(TInteractionNode), pointer :: it
  !  it => negf%interactList%first
  !
  !  do while (associated(it))
  !    select type(pInter => it%inter)
  !    class is(ElPhonInel)
  !       call pInter%destroy()
  !    end select
   !   it => it%next
  !  end do

  !end subroutine interaction_cleanup

  !---------------------------------------------------------------------------
  subroutine compute_sigmas_inelastic(negf)
    type(TNegf) :: negf

    type(TInteractionNode), pointer :: it
    it => negf%interactList%first

    do while (associated(it))
      select type(pInter => it%inter)
      class is (TInelastic)
        call it%inter%destroy_Sigma_n()
        call it%inter%compute_Sigma_n(spin=negf%spin)
        call it%inter%destroy_Sigma_r()
        call it%inter%compute_Sigma_r(spin=negf%spin)
      end select
      it => it%next
    end do
  end subroutine compute_sigmas_inelastic

  !---------------------------------------------------------------------------
  subroutine compute_sigmas_elastic(negf)
    type(TNegf) :: negf

    type(TInteractionNode), pointer :: it
    it => negf%interactList%first

    do while (associated(it))
      select type(pInter => it%inter)
      class is (TElastic)
        call it%inter%compute_Sigma_r(spin=negf%spin)
        call it%inter%compute_Sigma_n(spin=negf%spin)
      end select
      it => it%next
    end do
  end subroutine compute_sigmas_elastic

  !---------------------------------------------------------------------------
  function get_max_niter_elastic(interactList) result (maxiter)
    type(TInteractionList), intent(in) :: interactList
    integer :: maxiter

    type(TInteractionNode), pointer :: it
    it => interactList%first

    maxiter = 0
    do while (associated(it))
      select type(pInter => it%inter)
      class is(TElastic)
        if (pInter%scba_niter > maxiter) then
           maxiter = pInter%scba_niter
        end if
      end select
      it => it%next
    end do
  end function get_max_niter_elastic

  !---------------------------------------------------------------------------
  function get_max_niter_inelastic(interactList) result (maxiter)
    type(TInteractionList), intent(in) :: interactList
    integer :: maxiter

    type(TInteractionNode), pointer :: it
    it => interactList%first

    maxiter = 0
    do while (associated(it))
      select type(pInter => it%inter)
      class is(TInelastic)
        if (pInter%scba_niter > maxiter) then
           maxiter = pInter%scba_niter
        end if
      end select
      it => it%next
    end do
  end function get_max_niter_inelastic

  !---------------------------------------------------------------------------
  !>
  !  Calculate the equilibrium Retarded Green's function (extended diagonal)
  !  on a single energy point
  !  It groups calculation of leads, scba loop if any and deallocations of
  !  working arrays. This routine is used in contour integration and DOS and
  !
  !---------------------------------------------------------------------------
  subroutine compute_Gr(negf, outer, ncont, Ec, scba_error, Gr)
    type(Tnegf), intent(inout) :: negf
    complex(dp), intent(in) :: Ec
    integer, intent(in) :: outer, ncont
    Type(z_CSR), intent(out) :: Gr
    real(dp), intent(out) :: scba_error

    integer :: scba_iter, scba_niter, i1
    real(dp) :: ncyc, scba_tol
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS

    negf%readOldSGF = negf%readOldDM_SGFs

    call compute_contacts(Ec,negf,ncyc,Tlc,Tcl,SelfEneR,GS)

    scba_error = 0.0_dp
    scba_niter = get_max_niter_elastic(negf%interactList)

    !   Note: scba loop is only for elastic scattering because the
    !         on the contour integrations there are issues with z+/-wq
    call negf%scbaDriverElastic%init(tol = negf%scba_elastic_tol, dowrite = .false.)
    scba_iter = 0

    do while (.not.negf%scbaDriverElastic%is_converged() .and. scba_iter <= scba_niter)

      call negf%scbaDriverElastic%set_scba_iter(scba_iter, negf%interactList)

      if (allocated(Gr%nzval)) call destroy(Gr)

      if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Green`s funct ')
      call calculate_Gr(negf,Ec,SelfEneR,Tlc,Tcl,GS,Gr,outer)
      if (id0.and.negf%verbose.gt.VBT) call write_clock

      call negf%scbaDriverElastic%check_Mat_convergence(Gr)

      scba_error = negf%scbaDriverElastic%scba_error()

      if (negf%cartComm%rank == 0 .and. scba_niter>0) then
        call write_int_info(negf%verbose, VBT, 'scba elastic iteration',scba_iter)
        call write_real_info(negf%verbose, VBT, 'scba elastic error',scba_error)
      end if

      scba_iter = scba_iter + 1

    enddo

    call negf%scbaDriverElastic%destroy()

    do i1=1,ncont
      call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
    enddo


  end subroutine compute_Gr

  !------------------------------------------------------------------------------
  !
  !  Calculates the non equilibrium Green function Gn = -i G<
  !  on a single energy point on real axis.
  !  It groups calculation of leads, scba loop if any and deallocations of
  !  working arrays.
  !
  !-----------------------------------------------------------------------------
  subroutine compute_Gn(negf, outer, ncont, Ec, frm, Gn, scba_error, Gr)
    type(Tnegf), intent(inout) :: negf
    integer, intent(in) :: outer, ncont
    complex(dp), intent(in) :: Ec
    real(dp), dimension(:), intent(in) :: frm
    Type(z_CSR), intent(out) :: Gn
    real(dp), intent(out) :: scba_error
    Type(z_CSR), intent(out), optional :: Gr

    integer :: i1
    integer :: scba_iter, scba_elastic_iter, scba_niter_ela
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    real(dp) :: Er, ncyc, scba_elastic_tol, scba_elastic_error

    Er = real(Ec,dp)
    scba_elastic_tol = negf%scba_elastic_tol
    scba_niter_ela = get_max_niter_elastic(negf%interactList)

    ! ---------------------------------------------------------------------
    ! Compute contact GF
    ! ---------------------------------------------------------------------
    if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Contact SE ')
    call compute_contacts(Ec, negf, ncyc, Tlc, Tcl, SelfEneR, GS)
    if (id0.and.negf%verbose.gt.VBT) call write_clock
    call write_int_info(negf%verbose, VBT, 'Average number of iterations', int(ncyc))

    ! ---------------------------------------------------------------------
    ! Compute block tri-diagonal Gn and optionally Gr
    ! ---------------------------------------------------------------------
    if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute elastic SCBA ')
    call calculate_elastic_scba(negf, Er, SelfEneR, Tlc, Tcl, GS, frm, scba_niter_ela, &
                 & scba_elastic_tol, scba_elastic_iter, scba_elastic_error, outer, Gn, Gr)
    if (id0.and.negf%verbose.gt.VBT) call write_clock

    if (negf%cartComm%rank == 0 .and. negf%interactList%counter /= 0) then
      call write_int_info(negf%verbose, VBT, 'scba elastic iterations',scba_elastic_iter)
      call write_real_info(negf%verbose, VBT, 'scba elastic error',scba_elastic_error)
    end if

    scba_error = scba_elastic_error

    do i1=1,ncont
      call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
    enddo

  end subroutine compute_Gn

  !------------------------------------------------------------------------------
  !
  !  Calculates the non equilibrium Green function Gp = i G> = A - Gn
  !  on a single energy point on real axis.
  !  It groups calculation of leads, scba loop if any and deallocations of
  !  working arrays.
  !
  !-----------------------------------------------------------------------------
  subroutine compute_Gp(negf, outer, ncont, Ec, frm, Gp, scba_error)
    type(Tnegf), intent(inout) :: negf
    integer, intent(in) :: outer, ncont
    complex(dp), intent(in) :: Ec
    real(dp), dimension(:), intent(in) :: frm
    Type(z_CSR), intent(out) :: Gp
    real(dp), intent(out) :: scba_error

    Type(z_CSR) :: Gr, Gn, A
    complex(dp) :: minusOne = (-1.0_dp, 0.0_dp)

    call compute_Gn(negf, outer, ncont, Ec, frm, Gn, scba_error, Gr)
    call zspectral(Gr, Gr, 0, A)
    call prealloc_sum(A, Gn, minusOne, Gp)

    call destroy(Gr, Gn, A)

  end subroutine compute_Gp

  !---------------------------------------------------------------------------
  !   COMPUTATION OF CURRENTS - INTEGRATION OF T(E)
  !---------------------------------------------------------------------------
  subroutine electron_current(negf)
    type(Tnegf) :: negf

    integer :: size_ni, ii, ni, nf
    real(dp) :: mu1, mu2

    if (.not.allocated(negf%tunn_mat)) then
       return
    end if

    size_ni = size(negf%tunn_mat,2)

    ! If previous calculation is there, destroy it
    if (allocated(negf%currents)) call log_deallocate(negf%currents)
    call log_allocate(negf%currents,size_ni)

    negf%currents=0.0_dp

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

    negf%currents=0.0_dp
    do ii=1,size_ni
       negf%currents(ii)= integrate_el_meir_wingreen(negf%integration, &
             negf%curr_mat(:,ii), negf%Emin, negf%Emax, negf%Estep) * negf%g_spin
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
    Real(dp), Dimension(:), allocatable :: tun_mat
    Real(dp), Dimension(:), allocatable :: ledos
    Real(dp) :: mu1, mu2   ! contact potentials
    Real(dp) :: ncyc       ! stores average number of iters in decimation

    Integer :: i, icont, icpl      ! dummy counters
    Integer :: ncont               ! number of contacts
    Integer :: size_ni             ! emitter-collector contacts

    Integer :: Npoints             ! number of integration points
    Complex(dp) :: Ec              ! Energy point
    Complex(dp) :: delta
    Logical :: do_ledos            ! performs or not dos_proj

    ! Get out immediately if Emax<Emin
    if (negf%Emax.le.negf%Emin) then
       if(id0) write(*,*) '0 tunneling points;  current = 0.0'
       call log_allocate(negf%currents,1)
       negf%currents = 0.0_dp
       return
    endif
    !-------------------------------------------------------

    do_ledos = .false.
    if(negf%ndos_proj.gt.0) do_ledos=.true.
    ncont = negf%str%num_conts
    Npoints = size(negf%en_grid)
    ncyc=0
    size_ni = size(negf%ni)
    negf%readOldSGF = negf%readOldT_SGFs

    !-------------------------------------------------------

    call log_allocate(tun_mat,size_ni)
    call log_allocate(negf%tunn_mat,Npoints,size_ni)
    negf%tunn_mat = 0.0_dp

    if (do_ledos) then
       call log_allocate(negf%ldos_mat,Npoints,negf%ndos_proj)
       call log_allocate(ledos,negf%ndos_proj)
       negf%ldos_mat(:,:)=0.0_dp
    endif
    !-------------------------------------------------------


    !Loop on energy points: tunneling
    do i = 1, Npoints

       call write_Epoint(negf%verbose,negf%en_grid(i), size(negf%en_grid))
       if (negf%en_grid(i)%cpu /= id) cycle

       Ec = negf%en_grid(i)%Ec * negf%en_grid(i)%Ec
       negf%iE = negf%en_grid(i)%pt
       negf%iE_path = negf%en_grid(i)%pt_path

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
       call write_int_info(negf%verbose, VBT, 'Average number of iterations', int(ncyc))


       if (.not.do_ledos) then
          if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Tunneling ')

          call calculate_transmissions(negf, Ec, SelfEneR, negf%tun_proj, tun_mat)

          negf%tunn_mat(i,:) = tun_mat(:) * negf%kwght
       else
          if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Tunneling and DOS')
          ledos(:) = 0.0_dp

          call calculate_transmissions_and_dos(negf, Ec, SelfEneR, GS, negf%tun_proj, tun_mat, &
                                             &  negf%dos_proj, ledos)

          negf%tunn_mat(i,:) = tun_mat(:) * negf%kwght
          negf%ldos_mat(i,:) = ledos(:) * negf%kwght
       endif

       if (id0.and.negf%verbose.gt.VBT) call write_clock

       do icont=1,ncont
         call destroy(Tlc(icont),Tcl(icont),SelfEneR(icont),GS(icont))
       enddo

    enddo !Loop on energy

    !call destroy_en_grid()
    call log_deallocate(tun_mat)
    if(do_ledos) call log_deallocate(ledos)

  end subroutine phonon_tunneling

  !---------------------------------------------------------------------------
  subroutine phonon_current(negf)
    type(Tnegf) :: negf

    integer :: size_ni, ii, ni, nf

    size_ni = size(negf%tunn_mat,2)

    call log_allocate(negf%currents,size_ni)
    negf%currents=0.0_dp

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
    INTEGER :: i,i1,N,Npoints,imin,imax
    logical :: swapped

    curr=0.0_dp
    N=0
    destep=1.0d10
    Npoints=NINT((emax-emin)/estep);

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
    imax=Npoints

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

          curr=curr+(c1+c2)*(E4-E3)/2.0_dp

       enddo

    enddo

    if (swapped) curr = -1.0_dp*curr
    integrate_el = curr

  end function integrate_el

  !************************************************************************
  ! Function to integrate the current density I(E) !!! and get the current
  ! for meir_wingreen
  !************************************************************************
  function integrate_el_meir_wingreen(integration,TUN_TOT,emin,emax,estep)

    implicit none

    real(dp) :: integrate_el_meir_wingreen
    integer, intent(in) :: integration
    real(dp), intent(in) :: emin,emax,estep
    real(dp), dimension(:), intent(in) :: TUN_TOT

    real(dp), dimension(:), allocatable :: w
    REAL(dp) :: curr
    !REAL(dp) :: TT1,TT2
    !REAL(dp) :: E1,E2,curr
    INTEGER :: i,Npoints

    Npoints = size(TUN_TOT)

    ! performs the integration with Simpson's rule
    ! w = (1,4,2,4,2,4,...,4,1)/3*h
    ! h = (Emax - Emin)/(N-1)
    ! 1  2  3  4  5  6  7  8  9  10 11 12 13
    ! |  |  |  |  |  |  |  |  |  |  |  |  |
    ! 1  2  2  2  2  2  2  2  2  2  2  2  1  h/2
    ! 1  4  2  4  2  4  2  4  2  4  2  4  1  h/3   (N-1)/2  15h/45
    ! 1  3  3  2  3  3  2  3  3  2  3  3  1  3h/8  (N-1)/3
    ! 7  32 12 32 14 32 12 32 14 32 12 32 7  2h/45 (N-1)/4

    allocate(w(Npoints))
    select case(integration)
    case(integration_type%trapezoidal)
       w = estep*1.0_dp
       w(1) = estep*0.5_dp
       w(Npoints) = estep*0.5_dp
    case(integration_type%simpson13)
      do i = 1, Npoints
        w(i) = estep/3.0_dp*(mod(i-1,2)+1)**2
      end do
      do i = 3, Npoints-1, 2
        w(i) = w(i) + estep/3.0_dp
      end do
    case(integration_type%simpson38)
      w=9.0_dp/8.0_dp*estep
      w(1)=3.0_dp/8.0_dp*estep
      do i = 4, Npoints, 3
        w(i)=6.0_dp/8.0_dp*estep
      end do
    end select

    curr=0.0_dp

    do i = 1, Npoints
       curr=curr+TUN_TOT(i)*w(i)
    enddo

    deallocate(w)
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
    INTEGER :: i,i1,N,Npoints,imin,imax

    curr=0.0_dp
    N=0
    destep=1.0e10_dp
    Npoints=NINT((emax-emin)/estep);

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
    imax=Npoints

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

          curr=curr+(c1+c2)*(E4-E3)*(E4-E3)/2.0_dp

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
    INTEGER :: i,i1,N,Npoints,imin,imax

    curr=0.0_dp
    Npoints=NINT((emax-emin)/estep);

    TT1=TUN_TOT(1)
    do i = 0, 9
      E1=emin*i/10
      E2=emin*(i+1)/10
      c1=diff_bose(E1,kbT)*TT1
      c2=diff_bose(E2,kbT)*TT1
      curr=curr+(c1+c2)*emin/20.0_dp
    end do

    ! performs the integration with simple trapezium rule.
    ! Within each substep the tunneling is linearly interpolated
    ! Possibly perform a cubic-spline interpolation in future
    do i=0,Npoints-1

      E1=emin+estep*i
      TT1=TUN_TOT(i+1)
      E2=emin+estep*(i+1)
      TT2=TUN_TOT(i+2)

      c1=diff_bose(E1,kbT)*TT1
      c2=diff_bose(E2,kbT)*TT2

      curr=curr+(c1+c2)*estep/2.0_dp
    enddo

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
    qmulli = 0.0_dp

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

    qtot = 0.0_dp
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

  subroutine quasiEq_int_n(negf, mu_n, Ec, rho)
    !In/Out
    type(Tnegf), intent(inout) :: negf
    real(dp), dimension(:), intent(inout) :: rho
    real(dp), dimension(:), intent(in) :: Ec, mu_n

    !Work
    integer :: i, nr, ioffs
    integer ::  ncont, outer
    integer ::  Nz, Npoles
    complex(dp) :: Ez, ww, z1, z2, z_diff, ff
    type(z_CSR) :: Gr
    real(dp), dimension(:), allocatable :: wght, pnts, minE, maxE
    complex(dp), dimension(:), allocatable :: diag, temp
    real(dp) :: Omega, kbT, Lambda, scba_error

    kbT = maxval(negf%cont(:)%kbT_dm)
    ncont = negf%str%num_conts
    outer = negf%outer_blocks

    if (negf%n_poles.eq.0) then
      Lambda = 0.5d0* kbT * pi
    else
      Lambda = 2.d0* negf%n_poles * KbT * pi
    endif

    ! Omega considers maximum kT so interval is always large enough
    Omega = negf%n_kt * kbT
    Nz = size(mu_n)
    Npoles = negf%n_poles

    allocate(temp(negf%H%nrow))
    allocate(minE(Nz))
    allocate(maxE(Nz))

    minE(:) = Ec(:) - negf%deltaEc
    maxE(:) = Ec(:) + Omega
    temp = 0.0_dp

    z1 = minval(minE) + j*Lambda
    z2 = maxval(maxE) + j*Lambda
    z_diff = z2 - z1

    allocate(pnts(negf%Np_n(2)))
    allocate(wght(negf%Np_n(2)))

    call gauleg(0.0_dp, 1.0_dp, pnts, wght, negf%Np_n(2))
    do i = 1, negf%Np_n(2)
       if (mod(i-1,numprocs) .ne. id) cycle

       Ez = z1 + pnts(i) * z_diff
       call compute_Gr(negf, outer, ncont, Ez, scba_error, Gr)
       call log_allocate(diag, Gr%nrow)
       call getdiag(Gr,diag)

       do nr = 1,Nz
          !if (real(Ez) > minE(nr) .and. real(Ez) < maxE(nr)) then
             ff = fermi(Ez, mu_n(nr), kbT)
             ww = negf%g_spin * negf%kwght * wght(i) * z_diff * ff/(2.0_dp*pi)
             temp(nr) = temp(nr) + diag(nr)*ww
          !endif
       enddo

       call log_deallocate(diag)
       call destroy(Gr)
    enddo

    deallocate(pnts)
    deallocate(wght)


    allocate(wght(negf%Np_n(1)))
    allocate(pnts(negf%Np_n(1)))

    call gauleg(0.0_dp, 1.0_dp, pnts, wght, negf%Np_n(1))
    z1 = minval(minE(:))
    z2 = minval(minE(:)) + j*Lambda


    z_diff = z2 - z1

    ioffs = negf%Np_n(2)
    do i = 1, negf%Np_n(1)
         if (mod(i-1+ioffs, numprocs) .ne. id) cycle
         Ez = z1 + pnts(i) * z_diff

         call compute_Gr(negf, outer, ncont, Ez, scba_error, Gr)
         call log_allocate(diag, Gr%nrow)
         call getdiag(Gr,diag)
         do nr = 1,Nz
            ff = fermi(Ez, mu_n(nr), kbT)
            ww = negf%g_spin * ff * negf%kwght * wght(i) * z_diff / (2.0_dp*pi)
            temp(nr) = temp(nr) + diag(nr)*ww
         end do

         call log_deallocate(diag)
         call destroy(Gr)
    end do
    deallocate(wght)
    deallocate(pnts)

    ioffs = negf%Np_n(2) + negf%Np_n(1)
    if (Npoles.ne.0) then
      do nr = 1,Nz
          do i = 1,Npoles
             if (mod(i-1+ioffs, numprocs) .ne. id) cycle
             Ez = mu_n(nr) + j * kbT * pi * (2.0_dp*i - 1.0_dp)
             ww = -j * kbT * negf%g_spin *(1.0_dp, 0.0_dp)

             call compute_Gr(negf, outer, ncont, Ez, scba_error, Gr)
             call log_allocate(diag, Gr%nrow)
             call getdiag(Gr,diag)

             temp(nr) = temp(nr) + diag(nr) * ww
             call log_deallocate(diag)
             call destroy(Gr)
          end do
      end do
    endif

    rho(:) = real(j*(temp(:)-conjg(temp(:))),dp)

    deallocate(minE)
    deallocate(maxE)
    deallocate(temp)

  end subroutine quasiEq_int_n

  subroutine quasiEq_int_p(negf, mu_p, Ev, rho)
    !In/Out
    type(Tnegf), intent(inout) :: negf
    real(dp), dimension(:), intent(inout) :: rho
    real(dp), dimension(:), intent(in) :: Ev, mu_p

    !Work
    integer :: i, nr, ioffs
    integer ::  ncont, outer
    integer ::  Nz, Npoles
    complex(dp) :: Ez, ww, z1, z2, z_diff, ff
    type(z_CSR) :: Gr
    real(dp), dimension(:), allocatable :: wght, pnts, minE, maxE
    complex(dp), dimension(:), allocatable :: diag, temp
    real(dp) :: Omega, kbT, Lambda, scba_error

    kbT = maxval(negf%cont(:)%kbT_dm)
    ncont = negf%str%num_conts
    outer = negf%outer_blocks

    if (negf%n_poles.eq.0) then
      Lambda = 0.5d0* kbT * pi
    else
      Lambda = 2.d0* negf%n_poles * KbT * pi
    endif

    ! Omega considers maximum kT so interval is always large enough
    Omega = negf%n_kt * kbT
    Nz = size(mu_p)
    Npoles = negf%n_poles

    allocate(temp(negf%H%nrow))
    allocate(minE(Nz))
    allocate(maxE(Nz))

    minE(:) = Ev(:) - Omega
    maxE(:) = Ev(:) + negf%deltaEv
    temp = 0.0_dp

    !TODO: extremes and weights are correct but a bit convoluted
    !      (following contour_int_p_def): should evetually fix it
    z1 = maxval(maxE) + j*Lambda
    z2 = minval(minE) + j*Lambda
    z_diff = z2 - z1

    allocate(pnts(negf%Np_p(2)))
    allocate(wght(negf%Np_p(2)))

    call gauleg(0.0_dp, 1.0_dp, pnts, wght, negf%Np_p(2))
    do i = 1, negf%Np_p(2)
       if (mod(i-1,numprocs) .ne. id) cycle

       Ez = z1 + pnts(i) * z_diff
       call compute_Gr(negf, outer, ncont, Ez, scba_error, Gr)
       call log_allocate(diag, Gr%nrow)
       call getdiag(Gr,diag)

       do nr = 1,Nz
          !if (real(Ez) > minE(nr) .and. real(Ez) < maxE(nr)) then
             ff = fermi(-Ez, -mu_p(nr), kbT)
             ww = - negf%g_spin * negf%kwght * wght(i) * z_diff * ff/(2.0_dp*pi)
             temp(nr) = temp(nr) + diag(nr)*ww
          !endif
       enddo

       call log_deallocate(diag)
       call destroy(Gr)
    enddo

    deallocate(pnts)
    deallocate(wght)


    allocate(wght(negf%Np_p(1)))
    allocate(pnts(negf%Np_p(1)))

    call gauleg(0.0_dp, 1.0_dp, pnts, wght, negf%Np_p(1))
    z1 = maxval(maxE(:))
    z2 = maxval(maxE(:)) + j*Lambda


    z_diff = z2 - z1

    ioffs = negf%Np_p(2)
    do i = 1, negf%Np_p(1)
         if (mod(i-1+ioffs,numprocs) .ne. id) cycle
         Ez = z1 + pnts(i) * z_diff

         call compute_Gr(negf, outer, ncont, Ez, scba_error, Gr)
         call log_allocate(diag, Gr%nrow)
         call getdiag(Gr,diag)
         do nr = 1,Nz
            ff = fermi(-Ez, -mu_p(nr), kbT)
            ww = - negf%g_spin * ff * negf%kwght * wght(i) * z_diff / (2.0_dp*pi)
            temp(nr) = temp(nr) + diag(nr)*ww
         end do

         call log_deallocate(diag)
         call destroy(Gr)
    end do
    deallocate(wght)
    deallocate(pnts)

    ioffs = negf%Np_p(2) + negf%Np_p(1)
    if (Npoles.ne.0) then
      do nr = 1,Nz
          do i = 1,Npoles
             if (mod(i-1+ioffs,numprocs) .ne. id) cycle
             Ez = mu_p(nr) + j * kbT * pi * (2.0_dp*i - 1.0_dp)
             ww = j * kbT * negf%g_spin *(1.0_dp, 0.0_dp)

             call compute_Gr(negf, outer, ncont, Ez, scba_error, Gr)
             call log_allocate(diag, Gr%nrow)
             call getdiag(Gr,diag)

             temp(nr) = temp(nr) + diag(nr) * ww
             call log_deallocate(diag)
             call destroy(Gr)
          end do
      end do
    endif

    rho(:) = real(j*(temp(:)-conjg(temp(:))),dp)

    deallocate(minE)
    deallocate(maxE)
    deallocate(temp)

  end subroutine quasiEq_int_p

end module integrations

