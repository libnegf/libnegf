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

!--------------------------------------------------------------------
!>  Regression test for issue #6. Verify that we can initalize
!!  an instance of libnegf and immediately destroy. It was not
!!  previously possible due to problems in the inizalization of
!!  contacts, and this test would crash.
!--------------------------------------------------------------------
program main

  use libnegf
  use lib_param

  implicit none

  Type(Tnegf) :: pnegf
  Type(lnParams) :: params

  call init_negf(pnegf)
  call get_params(pnegf, params)
  call destroy_negf(pnegf)

end program main
