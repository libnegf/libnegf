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

!> Atom block diagonal elastic dephasing model

 module negf_elphdb

  use negf_ln_precision, only : dp
  use negf_interactions, only : interaction
  use negf_ln_allocation, only : log_allocate, log_deallocate
  use negf_ln_structure, only : TStruct_info
  use negf_mat_def, only : z_dns, create
  
  implicit none
  private

  public :: ElPhonDephB, ElPhonDephB_create

  type, extends(interaction) :: ElPhonDephB

    private
    !> Coupling squared per each atomic block, dimension energy^2
    type(z_DNS), allocatable, dimension(:) :: coupling
    !> Local block diagonal representation of retarded self energy
    !! Each z_dns is an atomic block
    type(z_DNS), allocatable, dimension(:) :: sigma_r
    !> Local block diagonal representation of lesser particle self energy
    !! Each z_dns is an atomic block
    type(z_DNS), allocatable, dimension(:) :: sigma_n
    !> Number of vibrational modes
    integer :: nummodes
    !> An array specifying how many orbital per atom (assumed contiguous)
    integer, allocatable, dimension(:) :: orbsperatm
    !> From orbsperatom, a work array containing starting orbital index
    !  for each atom, to speed up some patchworking
    integer, allocatable, dimension(:) :: atmorbstart
    !> For each atom, determine in which PL it sits, to accelerate block-sparse
    !  assignment. Used in model 2,3
    integer, allocatable, dimension(:) :: atmpl

  contains

    procedure :: add_sigma_r
    procedure :: get_sigma_n
    procedure :: set_Gr
    procedure :: set_Gn

  end type ElPhonDephB

contains

  !>
  ! Factory for el-ph dephasing diagonal model
  ! @param struct : contact/device partitioning infos
  ! @param coupling: coupling (energy units) 
  ! @param orbsperatom: number of orbitals per each atom
  ! @param niter: fixed number of scba iterations
  ! @param tol: scba tolerance
  subroutine ElPhonDephB_create(this, struct, coupling, orbsperatm, niter, tol)
    
    type(ElPhonDephB), intent(inout) :: this
    type(TStruct_info), intent(in) :: struct
    real(dp), dimension(:), intent(in) :: coupling
    integer, dimension(:), intent(in) :: orbsperatm
    integer, intent(in) :: niter
    real(dp), intent(in) :: tol

    integer :: ii, jj, ierr, natm

    this%descriptor = &
        & "Electron-Phonon dephasing model in atom block diagonal model"

    !Check input size
    if (size(coupling).ne.sum(orbsperatm)) then
      stop 'Error: coupling and orbsperatom not compatible'
    end if
    
    this%scba_niter = niter
    this%scba_tol = tol
    this%struct = struct
    this%orbsperatm = orbsperatm

    natm = size(this%orbsperatm)
    call log_allocate(this%atmorbstart, natm)
    call log_allocate(this%atmpl, natm)
    this%atmorbstart(1) = 1
    do ii = 2,natm
      this%atmorbstart(ii) = sum(this%orbsperatm(1:ii-1)) + 1
    enddo
    allocate(this%sigma_r(natm),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate sigma_r'
    allocate(this%sigma_n(natm),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate sigma_n'
    allocate(this%coupling(natm),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate coupling'
    do ii = 1,natm
      call create(this%sigma_r(ii),orbsperatm(ii),orbsperatm(ii))
      this%sigma_r(ii)%val = (0.0d0, 0.0d0)
      call create(this%sigma_n(ii),orbsperatm(ii),orbsperatm(ii))
      this%sigma_n(ii)%val = (0.0d0, 0.0d0)
      call create(this%coupling(ii),orbsperatm(ii),orbsperatm(ii))
      this%coupling(ii)%val = 0.d0
    end do
    this%nummodes = 1  !Single 0eV mode (Actually n localized modes, but we 
                       !treat them all contemporary)

    !Assign couplings on atom block diagonal
    do ii = 1,natm
      do jj = 1,this%orbsperatm(ii)
      this%coupling(ii)%val(jj,jj) = coupling(jj + this%atmorbstart(ii) - 1)
      end do
    end do
    ! Determine atmpl
    this%atmpl = 0
    do ii = 1,natm
      associate(pl_start=>this%struct%mat_PL_start)
      do jj = 1, size(pl_start) - 1
        if (this%atmorbstart(ii).ge.pl_start(jj).and. &
            this%atmorbstart(ii).lt.pl_start(jj + 1)) then
          this%atmpl(ii) = jj
        end if
      end do
      end associate
    end do
    ! Check that they are all assigned
    do ii = 1,natm
      if (this%atmpl(ii).eq.0) then
        write(*,*) this%atmpl
        stop 'atmpl not correctly set'
      end if
    end do

  end subroutine ElPhonDephB_create


  !> This interface should append
  !  the retarded self energy to ESH
  subroutine add_sigma_r(this, esh)
    class(ElPhonDephB) :: this
    type(z_dns), dimension(:,:), intent(inout) :: esh

    integer :: n, nbl, ii, indend, indstart, natm, norbs

    nbl = this%struct%num_PLs

    if (this%scba_iter .eq. 0) return

    natm = size(this%orbsperatm)

    do ii = 1,natm
      n = this%atmpl(ii)
      norbs = this%orbsperatm(ii)
      indstart = this%atmorbstart(ii) - this%struct%mat_PL_start(n) + 1
      indend = indstart + norbs - 1
      ESH(n,n)%val(indstart:indend, indstart:indend) = &
          ESH(n,n)%val(indstart:indend, indstart:indend) - &
          this%sigma_r(ii)%val(:,:)
    end do

  end subroutine add_sigma_r
  

  !> Returns the lesser (n) Self Energy in block format
  !  
  subroutine get_sigma_n(this, blk_sigma_n, en_index)
    class(ElPhonDephB) :: this
    type(z_dns), dimension(:,:), intent(inout) :: blk_sigma_n
    integer, intent(in) :: en_index

    integer :: n, nbl, ii, nrow, indend, indstart, natm, norbs

    nbl = this%struct%num_PLs
    natm = size(this%orbsperatm)

    if (this%scba_iter .eq. 0) return
    
    do n = 1, nbl
      nrow = this%struct%mat_PL_end(n) - this%struct%mat_PL_start(n) + 1
      if (.not.allocated(blk_sigma_n(n,n)%val)) then
        call create(blk_sigma_n(n,n), nrow, nrow)
      end if  
      blk_sigma_n(n,n)%val = (0.0_dp, 0.0_dp)
    end do
    do ii = 1,natm
      n = this%atmpl(ii)
      norbs = this%orbsperatm(ii)
      indstart = this%atmorbstart(ii) - this%struct%mat_PL_start(n) + 1
      indend = indstart + norbs - 1
      blk_sigma_n(n,n)%val(indstart:indend, indstart:indend) = &
          this%sigma_n(ii)%val(:,:)
    enddo

  end subroutine get_sigma_n

  !> Give the Gr at given energy point to the interaction
  subroutine set_Gr(this, Gr, en_index)
    class(ElPhonDephB) :: this
    type(z_dns), dimension(:,:), intent(in) :: Gr
    integer :: en_index

    integer :: n, npl, ii, indstart, indend, natm, norbs

    ! This model does not keep track of Gr at different energies because
    ! the model is local in energy. We directly calculate sigma_r
    npl = this%struct%num_PLs
    natm = size(this%orbsperatm)
    
    do ii = 1,natm
      n = this%atmpl(ii)
      norbs = this%orbsperatm(ii)
      indstart = this%atmorbstart(ii) - this%struct%mat_PL_start(n) + 1
      indend = indstart + norbs - 1
      this%sigma_r(ii)%val = matmul(matmul( &
          & this%coupling(ii)%val, &
          & Gr(n,n)%val(indstart:indend, indstart:indend)), &
          & this%coupling(ii)%val)
    end do

  end subroutine set_Gr

  !> Give the Gn at given energy point to the interaction
  subroutine set_Gn(this, Gn, en_index)
    class(ElPhonDephB) :: this
    type(z_dns), dimension(:,:), intent(in) :: Gn
    integer :: en_index

    integer :: n, npl, ii, indstart, indend, natm, norbs
    ! This model does not keep track of Gr at different energies because
    ! the model is local in energy. We directly calculate sigma_n
    npl = this%struct%num_PLs
    natm = size(this%orbsperatm)
    
    do ii = 1,natm
      n = this%atmpl(ii)
      norbs = this%orbsperatm(ii)
      indstart = this%atmorbstart(ii) - this%struct%mat_PL_start(n) + 1
      indend = indstart + norbs - 1
      this%sigma_n(ii)%val = matmul(matmul( &
          & this%coupling(ii)%val, &
          & Gn(n,n)%val(indstart:indend, indstart:indend)), &
          & this%coupling(ii)%val)
    end do

  end subroutine set_Gn


end  module negf_elphdb
