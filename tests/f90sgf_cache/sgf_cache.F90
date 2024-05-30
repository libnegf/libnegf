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

program main

  use libnegf

  implicit none

  Type(Tnegf), target :: negf
  Type(Tnegf), pointer :: pnegf
  Type(lnParams) :: params
  integer, allocatable :: surfstart(:), surfend(:), contend(:), plend(:), cblk(:)
  real(kind(1.d0)), allocatable :: mu(:), kt(:)

  type(z_csr) :: dm1, dm2, dm3, dm4
  !integer :: nnz, nrow
  !integer, allocatable :: rowpnt(:), colind(:)
  real(kind(1.d0)), allocatable :: nzval(:) !nzval1(:), nzval2(:), nzval3(:), nzval4(:)

  surfstart = [61,81]
  surfend = [60,80]
  contend = [80,100]
  plend = [60]
  mu = [0.d0, 0.d0]
  kt = [1.0d-3, 1.0d-3]
  cblk = [1,1]

  pnegf => negf

  call create(dm1, 60, 60, 60)
  call create(dm2, 60, 60, 60)
  call create(dm3, 60, 60, 60)
  call create(dm4, 60, 60, 60)
  allocate(nzval(60))

  call init_negf(pnegf)
  call init_contacts(pnegf, 2)
  call read_HS(pnegf, "H_real.dat", "H_imm.dat", 0)
  call set_S_id(pnegf, 100)
  call init_structure(pnegf, 2, surfstart, surfend, contend, 1, plend, cblk)

  ! Calculate a density with disk cache.
  call get_params(pnegf, params)
  params%readOldDM_SGFs = 2
  params%SGFcache = 0
  params%ec = -5.0
  params%kbT_dm = kt
  params%verbose = 101
  call set_params(pnegf, params)

  call create_scratch(pnegf)

  write(*,*) 'Compute density and save SGF'
  call compute_density_dft(pnegf)
  !call get_DM(pnegf, nnz, nrow, rowpnt, colind, nzval1)
  call get_DM(pnegf, dm1)
  call destroy_DM(pnegf)

  call set_readOldDMsgf(pnegf, 0)
  call compute_density_dft(pnegf)
  call get_DM(pnegf, dm2)
  call destroy_DM(pnegf)

  nzval = real(dm1%nzval - dm2%nzval)
  if (norm2(nzval) .gt. 1e-10) then
    error stop "Mismatch between density matrix after reloading surface green's functions"
  end if

  write(*,*) 'Compute and save SGF in memory'
  params%readOldDM_SGFs = 2
  params%SGFcache = 1
  call set_params(pnegf, params)
  call compute_density_dft(pnegf)
  call get_DM(pnegf, dm3)
  call destroy_DM(pnegf)
  nzval = real(dm2%nzval - dm3%nzval)

  if (norm2(nzval) .gt. 1e-10) then
    error stop "Mismatch between density matrix with SGF from disk and re-computed"
  end if

  call set_readOldDMsgf(pnegf, 0)
  call compute_density_dft(pnegf)
  call get_DM(pnegf, dm4)
  call destroy_DM(pnegf)
  nzval = real(dm3%nzval - dm4%nzval)

  if (norm2(nzval) .gt. 1e-10) then
    error stop "Mismatch between density matrix with SGF reading from memory"
  end if

  call destroy_negf(pnegf)

end program main
