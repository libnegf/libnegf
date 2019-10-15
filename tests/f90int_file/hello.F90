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

program hello

  use libnegf
  use lib_param

  implicit none

  Type(Tnegf), target :: negf
  Type(Tnegf), pointer :: pnegf
  real(kind(1.d0)), dimension(:,:), pointer :: transmission

  pnegf => negf

  write(*,*) 'Libnegf hello world'
  write(*,*) 'Init...'
  call init_negf(pnegf)
  call init_contacts(pnegf, 2)
  call read_negf_in(negf)

  write(*,*) 'Compute landauer tunneling and current'
  call compute_current(pnegf)
  call associate_transmission(pnegf, transmission)
  ! The above passes the transmission, but we write to file for debug
  call write_tunneling_and_dos(pnegf)
  write(*,*) 'Release libnegf'
  call destroy_negf(pnegf)

end program hello
