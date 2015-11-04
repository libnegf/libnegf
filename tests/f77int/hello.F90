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

  implicit none
  integer :: handler_size, en_size, ii
  integer, allocatable, dimension(:) :: handler1, handler2
  double precision, allocatable, dimension(:) :: trans, energies
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
  call negf_solve_landauer(handler1)
  call negf_get_current(handler1, 1, unitsH, unitsJ, current)
  write(*,*) 'Writing down tunneling for instance 1'
  call negf_write_tunneling_and_dos(handler1)
  write(*,*) 'Instance 2 will calculate tunneling'
  call negf_read_hs(handler2, realmat, imagmat, 0) 
  call negf_set_S_id(handler2, 100)
  call negf_read_input(handler2)
  call negf_solve_landauer(handler2)
  call negf_get_current(handler2, 1, unitsH, unitsJ, current)
  call negf_get_energygrid_size(handler2, en_size)
  write(*,*) 'en_size is',en_size
  allocate(energies(en_size))
  allocate(trans(en_size))
  call negf_get_transmission(handler2, 1, en_size, energies, trans)
  do ii=1,en_size
    write(*,*) 'Instance 2 energy and transmission: ',energies(ii), '  ', trans(ii)
  enddo
  write(*,*) 'Destroy negf'
  call negf_destruct_libnegf(handler1)
  call negf_destruct_session(handler1)
  call negf_destruct_libnegf(handler2)
  call negf_destruct_session(handler2)
  deallocate(energies)
  deallocate(trans)
  write(*,*) 'Done'

end program hello
