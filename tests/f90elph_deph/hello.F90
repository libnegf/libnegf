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
  use mpi_f08

  implicit none

  Type(Tnegf), target :: negf
  Type(Tnegf), pointer :: pnegf
  Type(lnParams) :: params
  integer, allocatable :: surfstart(:), surfend(:), contend(:), plend(:), cblk(:)
  real(kind(1.d0)), allocatable :: mu(:), kt(:), coupling(:)
  real(kind(1.d0)) :: current, sendbuff, recvbuff
  integer :: ierr, rank
  type(MPI_Comm) :: cartComm
  type(MPI_Comm) :: kComm, enComm

  call MPI_Init(ierr);

  surfstart = [61, 81]
  surfend = [60, 80]
  contend = [80, 100]
  plend = [10, 20, 30, 40, 50, 60]
  mu = [-0.5d0, 0.5d0]
  kt = [1.0d-3, 1.0d-3]
  cblk = [6, 1]

  pnegf => negf

  write(*,*) 'Initializing libNEGF'
  call init_negf(pnegf)

  write(*,*) 'Setup MPI communicator'
  !call set_mpi_bare_comm(pnegf, MPI_COMM_WORLD) 
  call set_cartesian_bare_comms(pnegf, MPI_COMM_WORLD, 1, cartComm, kComm, enComm)
  call MPI_Comm_rank(enComm, rank, ierr);

  write(*,*) 'Import Hamiltonian'
  call read_HS(pnegf, "H_real.dat", "H_imm.dat", 0)
  call set_S_id(pnegf, 100)

  call init_contacts(pnegf, 2)
  call init_structure(pnegf, 2, surfstart, surfend, contend, 6, plend, cblk)

  ! Here we set the parameters, only the ones different from default
  call get_params(pnegf, params)
  params%verbose = 100
  params%Emin = -2.d0
  params%Emax = 2.d0
  params%Estep = 1.d-1
  params%mu(1:2) = mu
  params%kbT_t(1:2) = kt
  call set_params(pnegf, params)

  ! Check values for 0 and finite coupling.
  if (rank==0) then
    write(*,*) '------------------------------------------------------------------ '
    write(*,*) 'Test 1 - current with coupling = 0'
    write(*,*) '------------------------------------------------------------------ '
  end if
  allocate(coupling(60))
  coupling = 0.0
  call set_elph_dephasing(pnegf, coupling, 5)
  if (rank==0) then
    write(*,*) 'Compute current'
  end if
  call compute_current(pnegf)
  ! The current should be 2: energy window is 1 and
  ! spin degeneracy is 2.
  sendbuff = pnegf%currents(1)
  call MPI_reduce(sendbuff, current, 1, MPI_DOUBLE, MPI_SUM, 0, encomm, ierr)

  if (rank==0) then
    write(*,*) 'Current ', current
    if (abs(current - 2.0) > 1e-4) then
       error stop "Wrong current for zero coupling"
    end if
  end if

  ! Compute for finite coupling, we should have a smaller current.
  if (rank==0) then
    write(*,*) '------------------------------------------------------------------ '
    write(*,*) 'Test 2 - current with coupling = 0.05'
    write(*,*) '------------------------------------------------------------------ '
  end if
  coupling = 0.05
  call destroy_interactions(pnegf)
  call set_elph_dephasing(pnegf, coupling, 5)
  call compute_current(pnegf)
  current = abs(pnegf%currents(1))
  sendbuff = current 
  call MPI_reduce(sendbuff, current, 1, MPI_DOUBLE, MPI_SUM, 0, encomm, ierr)
  if (rank==0) then
    write(*,*) 'Current with Meir-Wingreen', current 
    if (current > 1.95 .or. current < 1.90) then
       error stop "Wrong current for finite coupling"
    end if
  end if
  ! The current with dephasing is smaller than the ballistic current.

  !write(*,*) '------------------------------------------------------------------ '
  !write(*,*) 'Test 3 - cross-check using effective transmission'
  !write(*,*) '------------------------------------------------------------------ '
  ! Check that the effective transmission produces the same current.
  !call compute_dephasing_transmission(pnegf)
  !write(*,*) 'Current from effective transmission', pnegf%currents(1)
  !if (abs(pnegf%currents(1) - current) .gt. 1e-6) then
  !  error stop "Current evaluated with effective transmission does not match Meir Wingreen"
  !end if

  deallocate(coupling)
  call destroy_negf(pnegf)
  
  if (rank==0) then
    call writePeakInfo(6)
    call writeMemInfo(6)
  end if
  
  call MPI_finalize(ierr);

end program hello
