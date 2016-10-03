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

!> Overlap mask elastic dephasing model

module elphds

  use ln_precision, only : dp
  use interactions, only : interaction
  use ln_allocation, only : log_allocate, log_deallocate
  use ln_structure, only : TStruct_info
  use mat_def, only : z_csr, z_dns, create, destroy
  use sparsekit_drv, only : extract, zcsr2blk_sod, nzdrop, &
      & prealloc_mult, dns2csr, csr2dns, prealloc_sum
  
  implicit none
  private

  public :: ElPhonDephS, ElPhonDephS_create

  type, extends(interaction) :: ElPhonDephS

    private
    !> Coupling for each atomic mode in CSR form, dimension energy^2
    type(z_CSR), allocatable, dimension(:) :: couplings
    !> CSR retarded self energy for each mode
    type(z_CSR), allocatable, dimension(:) :: sigma_r
    !> CSR lesser self energy for each mode
    type(z_CSR), allocatable, dimension(:) :: sigma_n
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

  end type ElPhonDephS

contains

  !>
  ! Factory for el-ph dephasing diagonal model
  ! @param struct : contact/device partitioning infos
  ! @param coupling: coupling (energy units) 
  ! @param orbsperatom: number of orbitals per each atom
  ! @param over: overlap matrix in csr
  ! @param niter: fixed number of scba iterations
  ! @param tol: scba tolerance
  subroutine ElPhonDephS_create(this, struct, coupling, orbsperatm, over, niter, tol)
    
    type(ElPhonDephS), intent(inout) :: this
    type(TStruct_info), intent(in) :: struct
    real(dp), dimension(:), allocatable :: coupling
    integer, dimension(:), allocatable, intent(in) :: orbsperatm
    integer, intent(in) :: niter
    real(dp), intent(in) :: tol
    type(z_CSR), pointer, intent(in) :: over

    integer :: ii, jj, natm, ierr, norbs, offset, iimode
    type(z_CSR) :: work1, work2, work3, over_device
    real(dp), parameter :: coupling_tol = 1.0d-8

    this%descriptor = &
        & "Electron-Phonon dephasing model in overlap masked local oscillator model"
    !Check input size
    if (size(coupling).ne.sum(orbsperatm)) then
      stop 'Error: coupling and orbsperatom not compatible'
    end if
    
    this%scba_niter = niter
    this%scba_tol = tol
    this%struct = struct
    this%orbsperatm = orbsperatm
    natm = size(this%orbsperatm)
    norbs = size(coupling)
    this%nummodes = natm !One localized oscillator per atom
    call extract(over, 1, norbs, 1, norbs, over_device)
    !! Redefine number of modes taking away atoms with zero (or almost) coupling
    do ii = 1,natm
      offset = sum(orbsperatm(1:ii-1))
      !! if all couplings are zero, just skip it.
      if (all(coupling(offset+1:offset+orbsperatm(ii)) .lt. coupling_tol)) then
        this%nummodes = this%nummodes - 1
      end if
    end do
    
    allocate(this%sigma_r(this%nummodes), stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate csr_sigma_r'
    allocate(this%sigma_n(this%nummodes), stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate csr_sigma_n'
    allocate(this%couplings(this%nummodes), stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate csr_couplings'

    
    iimode = 0
    do ii = 1, natm
     offset = sum(orbsperatm(1:ii-1))
      !! if all couplings are zero, just skip it.
      if (all(coupling(offset+1:offset+orbsperatm(ii)) .lt. coupling_tol)) then
        cycle
      end if
      iimode = iimode + 1
      !! I assemble directly the (very) sparse CSR
      call create(work1, norbs, norbs, orbsperatm(ii))
      work1%rowpnt(1:offset) = 1
      work1%rowpnt(offset+1:norbs+1) = orbsperatm(ii) + 1
      do jj = 1,orbsperatm(ii)
        work1%nzval(jj) = coupling(jj + offset)
        work1%rowpnt(jj + offset) = jj
        work1%colind = jj + offset
      end do
      !! M -> M*S/2 + S*M/2
      call prealloc_mult(work1, over_device, (0.5d0, 0.0d0), work2)
      call prealloc_mult(over_device, work1, (0.5d0, 0.0d0), work3)
      call destroy(work1)
      call prealloc_sum(work2, work3, this%couplings(iimode))
      call destroy(work2)
      call destroy(work3)
    end do
    call destroy(over_device)

  end subroutine ElPhonDephS_create


  !> This interface should append
  !  the retarded self energy to ESH
  subroutine add_sigma_r(this, esh)
    class(ElPhonDephS) :: this
    type(z_dns), dimension(:,:), allocatable, intent(inout) :: esh

    type(z_dns), dimension(:,:), allocatable :: tmp_blk
    integer :: n, npl, ii, ierr, jj

    if (this%scba_iter .eq. 0) return
    npl = this%struct%num_PLs
    if (npl .ne. 1) then
      write(*,*) 'ElphPhonDephB works only with single PL'
      stop 0
    end if

    ! This could be done more performant, but this model won't see much use
    ! so I am leaving this way
    do ii = 1,this%nummodes
      allocate(tmp_blk(npl,npl),stat=ierr)
      if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate block-Matrix'
      call zcsr2blk_sod(this%sigma_r(ii), tmp_blk, this%struct%mat_PL_start)
      do jj = 1,npl
        ESH(jj, jj)%val = ESH(jj, jj)%val - tmp_blk(jj, jj)%val
        if (jj .lt. npl) then
          ESH(jj, jj + 1)%val = ESH(jj, jj + 1)%val - tmp_blk(jj, jj + 1)%val
          ESH(jj + 1, jj)%val = ESH(jj + 1, jj)%val - tmp_blk(jj + 1, jj)%val
        end if
      end do
      deallocate(tmp_blk)
    end do

  end subroutine add_sigma_r
  

  !> Returns the lesser (n) Self Energy in block format
  !  
  subroutine get_sigma_n(this, blk_sigma_n, en_index)
    class(ElPhonDephS) :: this
    type(z_dns), dimension(:,:), allocatable, intent(inout) :: blk_sigma_n
    integer, intent(in) :: en_index

    type(z_dns), dimension(:,:), allocatable :: tmp_blk
    integer :: n, npl, ii, nrow, ierr, jj

    if (this%scba_iter .eq. 0) return
    npl = this%struct%num_PLs
    if (npl .ne. 1) then
      write(*,*) 'ElphPhonDephB works only with single PL now'
      stop 0
    end if
    
    do n = 1, npl
      nrow = this%struct%mat_PL_end(n) - this%struct%mat_PL_start(n) + 1
      call create(blk_sigma_n(n,n), nrow, nrow)
      blk_sigma_n(n,n)%val = (0.0_dp, 0.0_dp)
    end do
      
    do ii = 1, this%nummodes
      allocate(tmp_blk(npl,npl),stat=ierr)
      if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate block-Matrix'
      call zcsr2blk_sod(this%sigma_n(ii), tmp_blk, this%struct%mat_PL_start)
      do jj = 1,npl
        blk_sigma_n(jj, jj)%val = blk_sigma_n(jj, jj)%val + tmp_blk(jj, jj)%val
        if (jj .lt. npl) then
          blk_sigma_n(jj, jj + 1)%val =  blk_sigma_n(jj, jj + 1)%val + &
              & tmp_blk(jj, jj + 1)%val
          blk_sigma_n(jj + 1, jj)%val = blk_sigma_n(jj + 1, jj)%val + &
              & tmp_blk(jj + 1, jj)%val
        end if
      end do
      deallocate(tmp_blk)
    end do

  end subroutine get_sigma_n

  !> Give the Gr at given energy point to the interaction
  subroutine set_Gr(this, Gr, en_index)
    class(ElPhonDephS) :: this
    type(z_dns), dimension(:,:), allocatable, intent(in) :: Gr
    integer :: en_index

    type(z_dns) :: work1, work2, work3
    integer :: n, npl, ii, natm, nnz

    !! Implement sigma = M*G*M, assuming that PL structure is not only
    !! preserved, but that only the corresponding Gr blocks are used
    !! Now dirty, only working without PL
    npl = this%struct%num_PLs
    if (npl .ne. 1) then
      write(*,*) 'ElphPhonDephB works only with single PL'
      stop 0
    end if  
    do ii=1,this%nummodes
      if (allocated(this%sigma_r(ii)%rowpnt)) then
        call destroy(this%sigma_r(ii))
      end if
      call create(work1, this%couplings(ii)%nrow, this%couplings(ii)%ncol)
      call csr2dns(this%couplings(ii), work1)
      call prealloc_mult(work1, Gr(1,1), work2)
      work1%val = conjg(transpose(work1%val))
      call prealloc_mult(work2, work1, work3)
      call destroy(work1)
      call destroy(work2)
      nnz = nzdrop(work3, 1.0d-12)
      call create(this%sigma_r(ii), Gr(1,1)%nrow, Gr(1,1)%ncol, nnz)
      call dns2csr(work3, this%sigma_r(ii))
      call destroy(work3)
    end do

  end subroutine set_Gr

  !> Give the Gn at given energy point to the interaction
  subroutine set_Gn(this, Gn, en_index)
    class(ElPhonDephS) :: this
    type(z_dns), dimension(:,:), allocatable, intent(in) :: Gn
    integer :: en_index

    type(z_dns) :: work1, work2, work3
    integer :: n, npl, ii, nnz
    !! Implement sigma = M*G*M^dagger, assuming that PL structure is not only
    !! preserved, but that only the corresponding Gr blocks are used
    !! Now dirty, only working without PL
    npl = this%struct%num_PLs
    if (npl .ne. 1) then
      write(*,*) 'ElphPhonDephB works only with single PL'
      stop 0
    end if   
    do ii=1,this%nummodes    
      if (allocated(this%sigma_n(ii)%rowpnt)) then
        call destroy(this%sigma_n(ii))
      end if 
      call create(work1, this%couplings(ii)%nrow, this%couplings(ii)%ncol)  
      call csr2dns(this%couplings(ii), work1)
      call prealloc_mult(work1, Gn(1,1), work2)
      work1%val = conjg(transpose(work1%val))
      call prealloc_mult(work2, work1, work3)
      call destroy(work1)
      call destroy(work2)
      nnz = nzdrop(work3, 1.0d-12)   
      call create(this%sigma_n(ii), Gn(1,1)%nrow, Gn(1,1)%ncol, nnz)  
      call dns2csr(work3, this%sigma_n(ii))
      call destroy(work3)   
    end do  
  end subroutine set_Gn


end module elphds
