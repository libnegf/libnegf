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
  use integrations

  implicit none

  Type(Tnegf), target :: negf
  Type(Tnegf), pointer :: pnegf

  pnegf => negf

  write(*,*) 'Libnegf hello world'
  write(*,*) 'Init...'
  call init_negf(pnegf)
  write(*,*) 'Import Hamiltonian'
  call read_HS(pnegf, "H_real.dat", "H_imm.dat", "H")
  call set_S_id(pnegf, 100)
  write(*,*) 'Import input file'
  call read_negf_in(pnegf)
  call negf_partition_info(pnegf)
  !!write(*,*) 'Compute density'
  !!call compute_density_dft(pnegf)
  write(*,*) 'Compute landauer tunneling and current'
  call compute_current(pnegf)
  call write_tunneling_and_dos(pnegf)
  write(*,*) 'Destroy negf'
  call destroy_negf(pnegf)
  write(*,*) 'Done'

end program hello
