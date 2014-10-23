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


module fermi_dist

  use ln_precision
  
  implicit none
  private

  public :: fermi
  
  interface fermi
    module procedure fermi_f
    module procedure fermi_fc
  end interface

contains
  
  real(kind=dp) function fermi_f(E,Ef,kT) 

    implicit none

    real(kind=dp), intent(in) :: E, Ef, kT

    ! the check over 0 is important otherwise the next fails
    if (kT.eq.0.d0) then
      if(E.gt.Ef) then
        fermi_f = 0.D0
      else
        fermi_f = 1.D0
      end if
      return
    endif

    if (abs((E-Ef)/kT).gt.30) then
      if(E.gt.Ef) then
        fermi_f = 0.D0
      else
        fermi_f = 1.D0
      end if
      return
    else        
      fermi_f = 1.d0/(dexp((E-Ef)/kT)+1.d0);
      return
    endif

  end function fermi_f


  complex(kind=dp) function fermi_fc(Ec,Ef,kT)

    implicit none
    
    complex(kind=dp), intent(in) :: Ec
    real(kind=dp), intent(in) :: Ef, kT

    complex(kind=dp) :: Efc,kTc,ONE=(1.d0,0.d0)

    Efc=Ef*ONE
    kTc=kT*ONE

    if (kT.eq.0.d0) then
      if(real(Ec).gt.Ef) then
        fermi_fc = (0.D0,0.D0)
      else
        fermi_fc = (1.D0,0.D0)
      end if
      return
    endif

    if (abs( (real(Ec)-Ef)/kT ).gt.30) then
      if(real(Ec).gt.Ef) then
        fermi_fc = (0.D0,0.D0)
      else
        fermi_fc = (1.D0,0.D0)
      end if
      return
    else        

      fermi_fc = ONE/(exp( (Ec-Efc)/kTc ) + ONE);
      return
    endif

  end function fermi_fc
  
end module fermi_dist
