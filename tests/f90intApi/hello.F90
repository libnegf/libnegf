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

  write(*,*) 'Libnegf API hello world'

  call negf_gethandlersize(handler_size)

  write(*,*) 'Handler size', handler_size
  allocate(handler1(handler_size))
  allocate(handler2(handler_size))
  write(*,*) 'I will create two libnegf instances'
  call negf_init_session(handler1)
  call negf_init_session(handler2)

  write(*,*) 'Done'

end program hello
