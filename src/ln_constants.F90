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
!! Conversion utility for length and energy units.
!! Example:
!! unitsin%name = "eV"
!! unitsout%name = "H"
!! unitsin%value = val
!! call convertUnits(unitin, unitout)
!! unitsout%value  contains converted value

!! Current conversion utility:
!! Define your energy units and desired current units ("A", "mA", ...)
!! I = value * convertCurrent(unitsOfEnergy, currUnits)
!! For Heat current, define desired output units ("W", "mW", ...)
!! IF unitsofEnergy is "unknown" (default), then convert returns 1.0
!!--------------------------------------------------------------------------!

module ln_constants
  
  use ln_precision
  
  complex(dp), parameter ::  j = (0.0_dp,1.0_dp)     ! CMPX unity
  complex(dp), parameter :: one = (1.0_dp,0.0_dp)        ! some 
  complex(dp), parameter :: minusone = (-1.0_dp,0.0_dp)  ! shortcuts
  complex(dp), parameter :: zero = (0.0_dp,0.0_dp)       !
  
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
  real(dp), parameter    :: g0= (pi**2)*Kb*Kb/(3.0_dp*hh) *ee ! W/K^2*T

  ! electric conductance quantum:  
  real(dp), parameter    :: Go= 2*ee/hh                   ! A/V 
 
  integer, parameter :: DELTA_SQ = 0
  integer, parameter :: DELTA_W = 1
  integer, parameter :: DELTA_MINGO = 2

  !!* Contains name of a units and its conversion factor
  type units
    character(8) :: name = "unknown "
    real(dp) :: value = 1.0_dp
  end type units
   

  !!* Length units ---------------------------------------------
  integer, parameter :: nLengthunits = 10

  type(units), parameter :: lengthUnits(nLengthUnits) = (/ &
      &units("unknown ", 1.0_dp), &
      &units("Bohr    ", 1.0_dp), &
      &units("au      ", 1.0_dp), &
      &units("AA      ", 1.0_dp/ATU ), &
      &units("Ang     ", 1.0_dp/ATU ), &
      &units("pm      ", 0.01_dp/ATU), &
      &units("nm      ", 10.0_dp/ATU), &
      &units("um      ", 1.0e4_dp/ATU), &
      &units("cm      ", 1.0e8_dp/ATU ), &
      &units("m       ", 1.0e10_dp/ATU) &
      /)

  
  !!* Energy units ---------------------------------------------
  integer, parameter :: nEnergyUnits = 9 

  type(units), parameter :: energyUnits(nEnergyUnits) = (/ &
      &units("unknown ", 1.0_dp), &
      &units("au      ", 1.0_dp), &
      &units("H       ", 1.0_dp), &
      &units("Ry      ", 0.5_dp), &
      &units("eV      ", 1.0_dp/HAR), &
      &units("kcal/mol", 0.0433634_dp/HAR), &
      &units("K       ", Kb), &
      &units("cm^-1   ", 1.239841930e-4_dp/HAR), &
      &units("J       ", J_eV/HAR) &
      &/)
  
   ! Internally current is computed in e/h * Energy 
   integer, parameter :: nCurrentUnits = 4
   type(units), parameter :: currentUnits(4) = (/ &
      &units("unknown ", 1.0_dp), &
      &units("A       ", 1.0_dp), &
      &units("mA      ", 1.0e-3_dp), &
      &units("nA      ", 1.0e-9_dp) &
      &/)
     
   ! Heat current is computed in 1/h * Energy^2
   integer, parameter :: nHeatCurrentUnits = 4
   type(units), parameter :: heatCurrentUnits(4) = (/ &
      &units("unknown ", 1.0_dp), &
      &units("W       ", 1.0_dp), &
      &units("mW      ", 1.0e-3_dp), &
      &units("nW      ", 1.0e-9_dp) &
      &/)

   ! Thermal conductance is computed in kb/h * Energy
   integer, parameter :: nHeatCondUnits = 4
   type(units), parameter :: heatCondUnits(4) = (/ &
      &units("unknown ", 1.0_dp), &
      &units("W/K     ", 1.0_dp), &
      &units("mW/K    ", 1.0e-3_dp), &
      &units("nW/K    ", 1.0e-9_dp) &
      &/)

contains 
 
  subroutine convertUnits(unitsin,unitsout)
    type(units), intent(in) :: unitsin
    type(units),intent(inout) :: unitsout

    integer :: ii, jj
    
    do ii = 1, nEnergyUnits
      if (unitsin%name .eq. energyUnits(ii)%name) then
        do jj = 1, nEnergyUnits
          if (unitsout%name .eq. energyUnits(jj)%name) then
             unitsout%value = unitsin%value * &
               &  energyUnits(ii)%value / energyUnits(jj)%value
             return 
          endif
        end do
      end if
    end do

    do ii = 1, nLengthUnits
      if (unitsin%name .eq. lengthUnits(ii)%name) then
        do jj = 1, nLengthUnits
          if (unitsout%name .eq. lengthUnits(jj)%name) then
             unitsout%value = unitsin%value * &
               &  lengthUnits(ii)%value / lengthUnits(jj)%value
             return 
          endif
        end do
      end if
    end do

    unitsout%name="unknown "

  end subroutine
 
  
  ! Compute conversion factor of current unitss
  ! [j] = e/h * [Energy]
  ! Current output is converted to Amp
  function convertCurrent(unitsOfH,currUnits) result(curr)
    type(units), intent(in) :: unitsOfH
    type(units), intent(in) :: currUnits
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
      if (currunits%name .eq. currentUnits(ii)%name) then
        curr = curr * currentUnits(ii)%value
      end if
    end do

  end function convertCurrent

  function convertHeatCurrent(unitsOfH,currUnits) result(curr)
    type(units), intent(in) :: unitsOfH
    type(units), intent(in) :: currUnits
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
      if (currunits%name .eq. heatCurrentUnits(ii)%name) then
        curr = curr * heatCurrentUnits(ii)%value
      end if
    end do

  end function convertHeatCurrent

  function convertHeatConductance(unitsOfH,condUnits) result(cond)
    type(units), intent(in) :: unitsOfH
    type(units), intent(in) :: condUnits
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

