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


module ln_precision
  use :: iso_fortran_env, only : int64
  implicit none
  private

  public :: sp, dp
  public :: EPS, EPS10, EPS12, EPS15, SGFACC
  public :: get_machine_prec
  public :: set_drop
  public :: long

  integer, parameter :: sp = selected_real_kind(6,30)
  integer, parameter :: dp = selected_real_kind(14,100)

  real(dp) :: EPS
  real(dp), parameter :: EPS10=1e-10_dp
  real(dp), parameter :: EPS12=1e-12_dp
  real(dp), parameter :: EPS15=1e-15_dp
  real(dp), parameter :: SGFACC=EPS10*10.0_dp
  integer, parameter :: long = int64

contains

  function get_machine_prec() result(racc)
    interface
      real(dp) function DLAMCH(C)
        import :: dp
        character C
      end function DLAMCH
    end interface
    real(dp) :: racc
    racc=DLAMCH('Precision')
  end function get_machine_prec

  subroutine set_drop(drop)
     real(dp) :: drop
     EPS = drop
  end subroutine set_drop

end module ln_precision
