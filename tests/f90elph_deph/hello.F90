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
  Type(lnParams) :: params
  integer, allocatable :: surfend(:), contend(:), plend(:), cblk(:)
  real(kind(1.d0)), allocatable :: mu(:), kt(:), coupling(:)
  real(kind(1.d0)), dimension(:,:), pointer :: transmission

  surfend = [60, 80]
  contend = [80, 100]
  plend = [10, 20, 30, 40, 50, 60]
  mu = [0.5d0, -0.5d0]
  kt = [1.0d-3, 1.0d-3]
  cblk = [6, 1]

  pnegf => negf

  write(*,*) 'Libnegf hello world'
  write(*,*) 'Init...'
  call init_negf(pnegf)
  call init_contacts(pnegf, 2)
  write(*,*) 'Import Hamiltonian'
  call read_HS(pnegf, "H_real.dat", "H_imm.dat", 0)
  call set_S_id(pnegf, 100)

  call init_structure(pnegf, 2, contend, surfend, 6, plend, cblk)

  ! Here we set the parameters, only the ones different from default
  call get_params(pnegf, params)
  params%Emin = -2.d0
  params%Emax = 2.d0
  params%Estep = 1.d-1
  params%mu = mu
  params%kbT_t = kt
  call set_params(pnegf, params)

  ! Check values for 0 and finite coupling.
  allocate(coupling(60))
  coupling = 0.0
  call set_elph_dephasing(pnegf, coupling, 5)
  write(*,*) 'Compute current'
  call compute_current(pnegf)
  ! The current should be 2: energy window is 1 and
  ! spin degeneracy is 2.
  write(*,*) 'Current ', pnegf%currents(1)
  if (abs(pnegf%currents(1) - 2.0) > 1e-4) error stop "Wrong current for zero coupling"

  ! Compute for finite coupling, we should have a smaller current.
  coupling = 0.05
  call destroy_elph_model(pnegf)
  call set_elph_dephasing(pnegf, coupling, 5)
  write(*,*) 'Compute current'
  call compute_current(pnegf)
  write(*,*) 'Current ', pnegf%currents(1)
  ! The current with dephasing is smaller than the ballistic current.
  if (pnegf%currents(1) > 1.95 .or. pnegf%currents(1) < 1.90) error stop "Wrong current for finite coupling"
  deallocate(coupling)
  write(*,*) 'Destroy negf'
  call destroy_negf(pnegf)
  write(*,*) 'Done'

end program hello
