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

!!--------------------------------------------------------------------------!
!! Unit conversion utility for length and energy units. 
!! Example:
!! unitin%name = "eV"
!! unitout%name = "H"
!! unitin%value = val
!! call convertUnits(unitin, unitout)
!! unitout%value  contains converted value

!! Current conversion utility:
!! Define your energy units and desired current units ("A", "mA", ...)
!! I = value * convertCurrent(UnitsofEnergy, currUnits)
!! For Heat current, define desired output units ("W", "mW", ...)
!! IF unitsofEnergy is "unknown" (default), then convert returns 1.0
!!--------------------------------------------------------------------------!

module ln_constants
  
  use ln_precision
  
  COMPLEX(dp), PARAMETER ::    j = (0.d0,1.d0)  ! CMPX unity
  
  real(dp), parameter    :: pi = 3.1415926535897932384_dp ! Greek p real
  real(dp), parameter    :: hh = 4.135667525931100d-15  ! eV s
  real(dp), parameter    :: ee = 1.60217656d-19         ! C  
  real(dp), parameter    :: a0 = 0.529177249_dp         ! Ang

  real(dp), parameter    :: eovh= 1.054181532603995d-3  ! A/H
  real(dp), parameter    :: Kb = 3.166811389002312d-6   ! H/K
  real(dp), parameter    :: oneovh = 2.8685739606828d-2 ! W/H^2
 !real(dp), parameter    :: eovh= 3.874045846176399d-5  ! A/eV
 !real(dp), parameter    :: Kb = 8.617332411853536d-5   ! eV/K
 !real(dp), parameter    :: oneovh = 3.8740458461764d-5 ! W/eV^2

  ! conversion factors
  real(dp), parameter    :: HAR = 27.21138506_dp        ! H/eV
  real(dp), parameter    :: ATU = a0                    ! a.u./Ang
  real(dp), parameter    :: AA_Bohr = 1.0_dp/ATU        ! Ang/Bohr
  real(dp), parameter    :: e2 = HAR*ATU                ! eV * Ang
  real(dp), parameter    :: J_eV = 1.0_dp/ee            ! J/eV 

  ! Heat conductance quantum:  
  ! Note, g0= 9.464309602837370d-013 W/K^2*T with these units
  ! This contrast literature: g0 = 9.456d-13 W/K^2 T
  real(dp), parameter    :: g0= (pi**2)*Kb*Kb/(3.d0*hh) *ee ! W/K^2*T

  ! electric conductance quantum:  
  real(dp), parameter    :: Go= 2*ee/hh                   ! A/V 
 

  !!* Contains name of a unit and its conversion factor
  type unit
    character(8) :: name = "unknown "
    real(dp) :: value = 1.0_dp
  end type unit
   

  !!* Length units ---------------------------------------------
  integer, parameter :: nLengthUnits = 10

  type(unit), parameter :: lengthUnits(nLengthUnits) = (/ &
      &unit("unknown ", 1.0_dp), &
      &unit("Bohr    ", 1.0_dp), &
      &unit("au      ", 1.0_dp), &
      &unit("AA      ", 1.0_dp/ATU ), &
      &unit("Ang     ", 1.0_dp/ATU ), &
      &unit("pm      ", 0.01_dp/ATU), &
      &unit("nm      ", 10.0_dp/ATU), &
      &unit("um      ", 1.0e4_dp/ATU), &
      &unit("cm      ", 1.0e8_dp/ATU ), &
      &unit("m       ", 1.0e10_dp/ATU) &
      /)

  
  !!* Energy units ---------------------------------------------
  integer, parameter :: nEnergyUnits = 9 

  type(unit), parameter :: energyUnits(nEnergyUnits) = (/ &
      &unit("unknown ", 1.0_dp), &
      &unit("au      ", 1.0_dp), &
      &unit("H       ", 1.0_dp), &
      &unit("Ry      ", 0.5_dp), &
      &unit("eV      ", 1.0_dp/HAR), &
      &unit("kcal/mol", 0.0433634_dp/HAR), &
      &unit("K       ", Kb), &
      &unit("cm^-1   ", 1.239841930e-4_dp/HAR), &
      &unit("J       ", J_eV/HAR) &
      &/)
  
   ! Internally current is computed in e/h * Energy 
   integer, parameter :: nCurrentUnits = 4
   type(unit), parameter :: currentUnits(4) = (/ &
      &unit("unknown ", 1.0_dp), &
      &unit("A       ", 1.0_dp), &
      &unit("mA      ", 1.0e-3_dp), &
      &unit("nA      ", 1.0e-9_dp) &
      &/)
     
   ! Heat current is computed in 1/h * Energy^2
   integer, parameter :: nHeatCurrentUnits = 4
   type(unit), parameter :: heatCurrentUnits(4) = (/ &
      &unit("unknown ", 1.0_dp), &
      &unit("W       ", 1.0_dp), &
      &unit("mW      ", 1.0e-3_dp), &
      &unit("nW      ", 1.0e-9_dp) &
      &/)

   ! Thermal conductance is computed in kb/h * Energy
   integer, parameter :: nHeatCondUnits = 4
   type(unit), parameter :: heatCondUnits(4) = (/ &
      &unit("unknown ", 1.0_dp), &
      &unit("W/K     ", 1.0_dp), &
      &unit("mW/K    ", 1.0e-3_dp), &
      &unit("nW/K    ", 1.0e-9_dp) &
      &/)

contains 
 
  subroutine convertUnits(unitin,unitout)
    type(unit), intent(in) :: unitin
    type(unit),intent(inout) :: unitout

    integer :: ii, jj
    
    do ii = 1, nEnergyUnits
      if (unitin%name .eq. energyUnits(ii)%name) then
        do jj = 1, nEnergyUnits
          if (unitout%name .eq. energyUnits(jj)%name) then
             unitout%value = unitin%value * &
               &  energyUnits(ii)%value / energyUnits(jj)%value
             return 
          endif
        end do
      end if
    end do

    do ii = 1, nLengthUnits
      if (unitin%name .eq. lengthUnits(ii)%name) then
        do jj = 1, nLengthUnits
          if (unitout%name .eq. lengthUnits(jj)%name) then
             unitout%value = unitin%value * &
               &  lengthUnits(ii)%value / lengthUnits(jj)%value
             return 
          endif
        end do
      end if
    end do

    unitout%name="unknown "

  end subroutine
 
  
  ! Compute conversion factor of current units
  ! [j] = e/h * [Energy]
  ! Current output is converted to Amp
  function convertCurrent(unitsOfH,currUnits) result(curr)
    type(unit), intent(in) :: unitsOfH
    type(unit), intent(in) :: currUnits
    real(dp) :: curr

    integer :: ii
    
    if (unitsOfH%name .eq. "unknown ") then
      curr = 1.0_dp
      return
    end if

    do ii = 1, nEnergyUnits
      if (unitsOfH%name .eq. energyUnits(ii)%name) then
         ! Current is transformed to Ampere
         curr = eovh * energyUnits(ii)%value   
      end if
    end do
    do ii = 1, nCurrentUnits
      if (currUnits%name .eq. currentUnits(ii)%name) then
        curr = curr * currentUnits(ii)%value
      end if
    end do

  end function convertCurrent

  function convertHeatCurrent(unitsOfH,currUnits) result(curr)
    type(unit), intent(in) :: unitsOfH
    type(unit), intent(in) :: currUnits
    real(dp) :: curr

    integer :: ii
    
    if (unitsOfH%name .eq. "unknown ") then
      curr = 1.0_dp
      return
    end if

    do ii = 1, nEnergyUnits
      if (unitsOfH%name .eq. energyUnits(ii)%name) then
         ! Heat Current is transformed to W
         curr = oneovh * (energyUnits(ii)%value)**2   
      end if
    end do
    do ii = 1, nHeatCurrentUnits
      if (currUnits%name .eq. heatCurrentUnits(ii)%name) then
        curr = curr * heatCurrentUnits(ii)%value
      end if
    end do

  end function convertHeatCurrent

  function convertHeatConductance(unitsOfH,condUnits) result(cond)
    type(unit), intent(in) :: unitsOfH
    type(unit), intent(in) :: condUnits
    real(dp) :: cond

    integer :: ii
  
    if (unitsOfH%name .eq. "unknown ") then
      cond = 1.0_dp
      return
    end if

    do ii = 1, nEnergyUnits
      if (unitsOfH%name .eq. energyUnits(ii)%name) then
         ! Thermal Conductance is transformed into W/K
         cond = kb*oneovh * energyUnits(ii)%value   
      end if
    end do

    do ii = 1, nHeatCondUnits
      if (condUnits%name .eq. heatCondUnits(ii)%name) then
        cond = cond * heatCondUnits(ii)%value
      end if
    end do

  end function convertHeatConductance

end module ln_constants

