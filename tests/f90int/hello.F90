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

  implicit none

  Type(Tnegf), target :: negf
  Type(Tnegf), pointer :: pnegf
  Type(lnParams) :: params
  integer, allocatable :: surfstart(:), surfend(:), contend(:), plend(:), cblk(:)
  real(kind(1.d0)), dimension(2) :: mu, kt
  real(kind(1.d0)), dimension(:,:), pointer :: transmission

  surfstart = [61,81]
  surfend = [60,80]
  contend = [80,100]
  plend = [60]
  mu = [0.d0, 0.d0]
  kt = [1.0d-3, 1.0d-3]
  cblk = [1,1]

  pnegf => negf

  write(*,*) 'Libnegf hello world'
  write(*,*) 'Init...'
  call init_negf(pnegf)
  call init_contacts(pnegf, 2)
  write(*,*) 'Import Hamiltonian'
  call read_HS(pnegf, "H_real.dat", "H_imm.dat", 0)
  call set_S_id(pnegf, 100)
  call init_structure(pnegf, 2, surfstart, surfend, contend, 1, plend, cblk)

  ! Here we set the parameters, only the ones different from default
  call get_params(pnegf, params)
  params%Emin = -3.d0
  params%Emax = 3.d0
  params%Estep = 1.d-2
  call set_params(pnegf, params)

  write(*,*) 'Compute landauer tunneling and current'
  call compute_current(pnegf)
  call associate_transmission(pnegf, transmission)
  ! The above passes the transmission, but we write to file for debug
  call write_tunneling_and_dos(pnegf)
  write(*,*) 'Destroy negf'
  call destroy_negf(pnegf)

end program hello
