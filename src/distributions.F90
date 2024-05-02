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


module distributions

  use ln_precision
  use ln_constants, only : zero, one 
  implicit none
  private

  public :: bose
  public :: fermi
  public :: diff_bose

  interface bose
    module procedure bose_r
    module procedure bose_c
  end interface
  
  interface fermi
    module procedure fermi_f
    module procedure fermi_fc
  end interface


contains
  
  real(kind=dp) function fermi_f(E,Ef,kT) 

    implicit none

    real(kind=dp), intent(in) :: E, Ef, kT

    ! the check over 0 is important otherwise the next fails
    if (kT.eq.0.0_dp) then
      if(E.gt.Ef) then
        fermi_f = 0.0_dp
      else
        fermi_f = 1.0_dp
      end if
      return
    endif

    if (abs((E-Ef)/kT).gt.30) then
      if(E.gt.Ef) then
        fermi_f = exp(-(E-Ef)/kT) 
      else
        fermi_f = 1.0_dp
      end if
      return
    else        
      fermi_f = 1.0_dp/(exp((E-Ef)/kT)+1.0_dp);
      return
    endif

  end function fermi_f


  complex(kind=dp) function fermi_fc(Ec,Ef,kT)

    implicit none
    
    complex(kind=dp), intent(in) :: Ec
    real(kind=dp), intent(in) :: Ef, kT

    complex(kind=dp) :: Efc,kTc

    Efc=Ef*ONE
    kTc=kT*ONE

    if (kT.eq.0.0_dp) then
      if(real(Ec).gt.Ef) then
        fermi_fc = zero
      else
        fermi_fc = one
      end if
      return
    endif

    if (abs( (real(Ec)-Ef)/kT ).gt.30) then
      if(real(Ec).gt.Ef) then
        fermi_fc = exp( -(Ec-Efc)/kTc )
      else
        fermi_fc = one
      end if
      return
    else        

      fermi_fc = one/(exp( (Ec-Efc)/kTc ) + one);
      return
    endif

  end function fermi_fc
  
  !////////////////////////////////////////////////////////////////////////
  real(dp) function bose_r(E,kT) 
    real(dp), intent(in) :: E, kT

    ! the check over 0 is important otherwise the next fails
    if (kT.eq.0.0_dp) then
      bose_r = 0.0_dp
      return
    endif

    if (abs(E/kT).gt.30.0_dp) then
      bose_r = exp(-E/kT);
    else        
      bose_r = 1.0_dp/(exp(E/kT) - 1.0_dp);
    endif
     
  end function bose_r   
  
  complex(kind=dp) function bose_c(Ec,Ef,kT)
    complex(kind=dp), intent(in) :: Ec
    real(kind=dp), intent(in) :: Ef, kT

    complex(kind=dp) :: Efc,kTc

    Efc=Ef*ONE
    kTc=kT*ONE

    if (kT.eq.0.0_dp) then
      bose_c = zero
      return
    endif

    if (abs( real(Ec)/kT ).gt.30.0_dp) then
      bose_c = exp( -Ec/kTc )
    else        
      bose_c = ONE/(exp( Ec/kTc ) - ONE);
    endif

  end function bose_c
  
  !////////////////////////////////////////////////////////////////////////
  real(dp) function diff_bose(E,kT) 
    real(dp), intent(in) :: E, kT

    ! the check over 0 is important otherwise the next fails
    if (kT.eq.0.0_dp) then
      diff_bose = 0.0_dp
      return
    endif

    if (E.eq.0.0_dp) then
      diff_bose = 1.0_dp
      return
    endif

    if (abs(E/kT).gt.30.0_dp) then
      diff_bose = exp(-E/kT)*(E/kT)**2
    else        
      diff_bose = exp(E/kT)/((exp(E/kT)-1.0_dp)**2)*(E/kT)**2
    endif
     
  end function diff_bose
  
end module distributions
