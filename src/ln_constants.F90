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


module ln_constants
  
  use ln_precision
  
  real(dp), parameter    :: pi =  3.14159265358979323844_dp ! Greek p real
  real(dp), parameter    :: eovh = (1.05420882d-3)      ! A/H
  real(dp), parameter    :: HAR = 27.2113845_dp         ! H/eV
  real(dp), parameter    :: ATU = 0.529177249_dp        ! a.u./Ang
  real(dp), parameter    :: Kb = (3.166830814d-6)       ! H/K
  real(dp), parameter    :: oneovh = 0.180237804832302  ! W/H^2
  real(dp), parameter    :: hh =  1.054571726d-34       ! J s
  real(dp), parameter    :: HovJ =  4.3597441775d-18    ! H/J 
  

   
  COMPLEX(dp), PARAMETER ::    j = (0.d0,1.d0)  ! CMPX unity
  
  
end module ln_constants

