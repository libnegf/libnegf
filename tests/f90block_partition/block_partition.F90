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

program test_block_partition

  use negf_libnegf
  use negf_lib_param

  implicit none

  Type(Tnegf), target :: negf
  Type(Tnegf), pointer :: pnegf
  Type(lnParams) :: params
  integer, allocatable :: surfend(:), contend(:), plend(:), cblk(:)
  real(kind(1.d0)), allocatable :: mu(:), kt(:)
  real(kind(1.d0)), dimension(:,:), pointer :: transmission

  ! Use the the partition algorithm on a linear chain.
  ! We should get blocks of the minimum allowed size (2).
  surfend = [60,80]
  contend = [80,100]
  mu = [0.d0, 0.d0]
  kt = [1.0d-3, 1.0d-3]

  pnegf => negf

  call init_negf(pnegf)
  call init_contacts(pnegf, 2)
  write(*,*) 'Import Hamiltonian'
  call read_HS(pnegf, "H_real.dat", "H_imm.dat", 0)
  call set_S_id(pnegf, 100)
  call init_structure(pnegf, 2, contend, surfend, 0, plend, cblk)

  ! Verify that the structure is correct. We should have
  ! 30 blocks and interaction with the first and last.
  write(*,*) pnegf%str%cblk(:)
  if (pnegf%str%num_PLs .ne. 30) error stop "Wrong number of PLs"
  if (pnegf%str%cblk(1) .ne. 30) error stop "Wrong contact interaction"
  if (pnegf%str%cblk(2) .ne. 1) error stop "Wrong contact interaction"

  ! Calculate also a transmission as sanity check.
  call get_params(pnegf, params)
  params%Emin = -1.d0
  params%Emax = 1.d0
  params%Estep = 1.d-1
  call set_params(pnegf, params)

  call compute_current(pnegf)
  call destroy_negf(pnegf)

end program test_block_partition
