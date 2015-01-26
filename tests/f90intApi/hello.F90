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

include 'libnegf.h'
program hello

  implicit none
  integer :: handler_size
  integer, allocatable, dimension(:) :: handler1, handler2
  double precision :: current
  character(500) :: realmat
  character(500) :: imagmat
  character(500) :: unitsH, unitsJ

  realmat="H_real.dat"
  imagmat="H_imm.dat"
  unitsH = "H"
  unitsH = "J"
  write(*,*) 'Libnegf API hello world'

  call negf_gethandlersize(handler_size)

  write(*,*) 'Handler size', handler_size
  allocate(handler1(handler_size))
  allocate(handler2(handler_size))
  write(*,*) 'Create two libnegf instances'
  call negf_init_session(handler1)
  call negf_init_session(handler2)
  call negf_init(handler1)
  call negf_init(handler2)
  write(*,*) 'Instance 1 at address',handler1
  write(*,*) 'Instance 2 at address',handler2
  write(*,*) 'Instance 1 will calculate tunneling'
  call negf_read_hs(handler1, realmat, imagmat, 0) 
  call negf_set_S_id(handler1, 100)
  call negf_read_input(handler1)
  call negf_current(handler1, current, unitsH, unitsJ)
  write(*,*) 'Instance 2 will calculate tunneling'
  call negf_read_hs(handler2, realmat, imagmat, 0) 
  call negf_set_S_id(handler2, 100)
  call negf_read_input(handler2)
  call negf_current(handler2, current, unitsH, unitsJ)
  write(*,*) 'Destroy negf'
  call negf_destruct_libnegf(handler1)
  call negf_destruct_session(handler1)
  call negf_destruct_libnegf(handler2)
  call negf_destruct_session(handler2)
  write(*,*) 'Done'

end program hello
