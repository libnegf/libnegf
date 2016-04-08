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

!> The module implements an abstract class to interface different
!! many body interactions. 

module interactions

  use globals, only : LST
  use ln_precision, only : dp

  implicit none
  private

  public :: interaction

  type, abstract :: Interaction

    character(len=LST) :: descriptor
    !> Maximum number of SCBA iterations. 
    !! corresponds to no iterations (self energy is not calculated)
    integer :: scba_niter = 0
    !> Keep track of SCBA iteration 
    integer :: scba_iter = 0
    !> SCBA Tolerance
    real(dp) :: scba_tol = 1.0d-7

  contains

    procedure(abst_destroy), deferred :: destroy

  end type Interaction

  abstract interface
    subroutine abst_destroy(self)
      import interaction
      class(interaction) :: self
    end subroutine abst_destroy
  end interface

end module interactions
