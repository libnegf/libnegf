!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Greens functions.         !
!! Copyright (C) 2012 - 2026                                                !
!!                                                                          !
!! This file is part of libNEGF: a library for                              !
!! Non Equilibrium Green's Functions calculations                           !
!!                                                                          !
!! Developers: Alessandro Pecchia, Daniele Soccodato                        !
!! Former Contributors: Gabriele Penazzi, Luca Latessa, Aldo Di Carlo       !
!!                                                                          !
!! libNEGF is free software: you can redistribute and/or modify it          !
!! under the terms of the GNU Lesser General Public License as published    !
!! by the Free Software Foundation, either version 3 of the License, or     !
!! (at your option) any later version.                                      !
!!                                                                          !
!!  You should have received a copy of the GNU Lesser General Public        !
!!  License along with libNEGF.  If not, see                                !
!!  <http://www.gnu.org/licenses/>.                                         !
!!--------------------------------------------------------------------------!

!> Overlap mask elastic dephasing model
#:include "types.fypp"

module elphds

  use ln_precision, only : sp, dp
  use interactions, only : TInteraction
  use ln_elastic, only : TElastic
  use ln_allocation, only : log_allocate, log_deallocate
  use ln_structure, only : TStruct_info
  use mat_def, only : c_csr, c_dns, z_csr, z_dns, create, destroy
  use sparsekit_drv, only : extract, zcsr2blk_sod, nzdrop, &
      & prealloc_mult, dns2csr, csr2dns, prealloc_sum

  implicit none
  private

  public :: ElPhonDephS, ElPhonDephS_create
  public :: ElPhonDephS_init

  type, extends(TElastic) :: ElPhonDephS

    private
    !> Electron-phonon Coupling for each atomic mode in CSR form, dimension energy
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

    procedure :: add_sigma_r_sp
    procedure :: add_sigma_r_dp
    procedure :: add_sigma_n_sp
    procedure :: add_sigma_n_dp
    procedure :: get_sigma_n_blk
    procedure :: get_sigma_n_mat
    procedure :: set_Gr
    procedure :: set_Gn
    procedure :: compute_Sigma_r
    procedure :: compute_Sigma_n
    procedure :: destroy_Sigma_r
    procedure :: destroy_Sigma_n
    procedure :: destroy => destroy_all

  end type ElPhonDephS

contains

  subroutine ElPhonDephS_create(this)
    class(TInteraction), allocatable :: this
    allocate(ElPhonDephS::this)
  end subroutine ElPhonDephS_create

  !>
  ! Factory for el-ph dephasing diagonal model
  ! @param struct : contact/device partitioning infos
  ! @param coupling: coupling (energy units)
  ! @param orbsperatom: number of orbitals per each atom
  ! @param over: overlap matrix in csr
  ! @param niter: fixed number of scba iterations
  ! @param tol: scba tolerance
  subroutine ElPhonDephS_init(this, struct, coupling, orbsperatm, over, niter)
    type(ElPhonDephS) :: this
    type(TStruct_info), intent(in) :: struct
    real(dp), dimension(:), intent(in) :: coupling
    integer, dimension(:), intent(in) :: orbsperatm
    integer, intent(in) :: niter
    type(z_CSR), pointer, intent(in) :: over

    integer :: ii, jj, natm, ierr, norbs, offset, iimode
    type(z_CSR) :: work1, work2, work3, over_device
    real(dp), parameter :: coupling_tol = 1.0d-8

    this%descriptor = &
        & "Electron-Phonon dephasing model in overlap masked local oscillator model"
    !Check input size
    if (size(coupling).ne.sum(orbsperatm)) then
      error stop 'Error: coupling and orbsperatom not compatible'
    end if

    this%scba_niter = niter
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

    this%wq = 0.0_dp   ! Zero energy mode

    allocate(this%sigma_r(this%nummodes), stat=ierr)
    if (ierr.ne.0) error stop 'ALLOCATION ERROR: could not allocate csr_sigma_r'
    allocate(this%sigma_n(this%nummodes), stat=ierr)
    if (ierr.ne.0) error stop 'ALLOCATION ERROR: could not allocate csr_sigma_n'
    allocate(this%couplings(this%nummodes), stat=ierr)
    if (ierr.ne.0) error stop 'ALLOCATION ERROR: could not allocate csr_couplings'

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

  end subroutine ElPhonDephS_init


  !> This interface should append
  !  the retarded self energy to ESH
#:def add_sigma_r_template(KIND,MTYPE)  
  subroutine add_sigma_r_${KIND}$(this, esh, en_index, k_index, spin)
    class(ElPhonDephS) :: this
    type(${MTYPE}$), dimension(:,:), intent(inout) :: esh
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    type(z_dns), dimension(:,:), allocatable :: tmp_blk
    integer :: npl, ii, ierr, jj

    if (this%scba_iter .eq. 0) return
    npl = this%struct%num_PLs
    if (npl .ne. 1) then
      error stop 'ElphPhonDephB works only with single PL'
    end if

    ! This could be done more performant, but this model won't see much use
    ! so I am leaving this way
    allocate(tmp_blk(npl,npl),stat=ierr)
    do ii = 1,this%nummodes
      if (ierr.ne.0) error stop 'ALLOCATION ERROR: could not allocate block-Matrix'
      call zcsr2blk_sod(this%sigma_r(ii), tmp_blk, this%struct%mat_PL_start)
      do jj = 1,npl
        ESH(jj, jj)%val = ESH(jj, jj)%val - tmp_blk(jj, jj)%val
        call destroy(tmp_blk(jj,jj))
        if (jj .lt. npl) then
          ESH(jj, jj + 1)%val = ESH(jj, jj + 1)%val - tmp_blk(jj, jj + 1)%val
          call destroy(tmp_blk(jj,jj+1))
          ESH(jj + 1, jj)%val = ESH(jj + 1, jj)%val - tmp_blk(jj + 1, jj)%val
          call destroy(tmp_blk(jj+1,jj))
        end if
      end do
    end do
    deallocate(tmp_blk)

  end subroutine add_sigma_r_${KIND}$
#:enddef add_sigma_r_template  

  !--------------------------------------------------------------------------
  !> This interface should append
  !  sigma_n to a passed self energy, sigma
#:def add_sigma_n_template(KIND,MTYPE)  
  subroutine add_sigma_n_${KIND}$(this, sigma, en_index, k_index, spin)
    class(ElPhonDephS) :: this
    type(${MTYPE}$), dimension(:,:), intent(inout) :: sigma
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    type(z_dns), dimension(:,:), allocatable :: tmp_blk
    integer :: npl, ii, ierr, jj

    if (this%scba_iter .eq. 0) return
    npl = this%struct%num_PLs
    if (npl .ne. 1) then
      error stop 'ElphPhonDephB works only with single PL'
    end if

    ! This could be done more performant, but this model won't see much use
    ! so I am leaving this way
    allocate(tmp_blk(npl,npl),stat=ierr)
    do ii = 1,this%nummodes
      if (ierr.ne.0) error stop 'ALLOCATION ERROR: could not allocate block-Matrix'
      call zcsr2blk_sod(this%sigma_r(ii), tmp_blk, this%struct%mat_PL_start)
      do jj = 1,npl
        sigma(jj, jj)%val = sigma(jj, jj)%val + tmp_blk(jj, jj)%val
        call destroy(tmp_blk(jj,jj))
        if (jj .lt. npl) then
          sigma(jj, jj + 1)%val = sigma(jj, jj + 1)%val + tmp_blk(jj, jj + 1)%val
          call destroy(tmp_blk(jj,jj+1))
          sigma(jj + 1, jj)%val = sigma(jj + 1, jj)%val + tmp_blk(jj + 1, jj)%val
          call destroy(tmp_blk(jj+1,jj))
        end if
      end do
    end do
    deallocate(tmp_blk)

  end subroutine add_sigma_n_${KIND}$
#:enddef add_sigma_n_template  

#:for PREC in PRECISIONS
     #:set KIND = PREC_ABBREVS[PREC]
     #:set MTYPE = MAT_TYPES['complex'][PREC]

     $:add_sigma_r_template(KIND,MTYPE)
     
     $:add_sigma_n_template(KIND,MTYPE)
#:endfor


  !> Returns the lesser (n) Self Energy in block format
  !
  subroutine get_sigma_n_blk(this, blk_sigma_n, en_index, k_index, spin)
    class(ElPhonDephS) :: this
    type(z_dns), dimension(:,:), intent(inout) :: blk_sigma_n
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    type(z_dns), dimension(:,:), allocatable :: tmp_blk
    integer :: n, npl, ii, nrow, ierr, jj

    if (this%scba_iter .eq. 0) return
    npl = this%struct%num_PLs
    if (npl .ne. 1) then
      error stop 'ElphPhonDephB works only with single PL now'
    end if

    do n = 1, npl
      nrow = this%struct%mat_PL_end(n) - this%struct%mat_PL_start(n) + 1
      if (.not.allocated(blk_sigma_n(n,n)%val)) then
        call create(blk_sigma_n(n,n), nrow, nrow)
      end if
      blk_sigma_n(n,n)%val = (0.0_dp, 0.0_dp)
    end do

    allocate(tmp_blk(npl,npl),stat=ierr)
    do ii = 1, this%nummodes
      if (ierr.ne.0) error stop 'ALLOCATION ERROR: could not allocate block-Matrix'
      call zcsr2blk_sod(this%sigma_n(ii), tmp_blk, this%struct%mat_PL_start)
      do jj = 1,npl
        blk_sigma_n(jj, jj)%val = blk_sigma_n(jj, jj)%val + tmp_blk(jj, jj)%val
        call destroy(tmp_blk(jj,jj))
        if (jj .lt. npl) then
          blk_sigma_n(jj, jj + 1)%val =  blk_sigma_n(jj, jj + 1)%val + &
              & tmp_blk(jj, jj + 1)%val
          call destroy(tmp_blk(jj,jj+1))
          blk_sigma_n(jj + 1, jj)%val = blk_sigma_n(jj + 1, jj)%val + &
              & tmp_blk(jj + 1, jj)%val
          call destroy(tmp_blk(jj+1,jj))
        end if
      end do
    end do
    deallocate(tmp_blk)

  end subroutine get_sigma_n_blk

  subroutine get_sigma_n_mat(this, sigma_n, ii, jj, en_index, k_index, spin)
    class(ElPhonDephS) :: this
    type(z_dns), intent(inout) :: sigma_n
    integer, intent(in) :: ii
    integer, intent(in) :: jj
    integer, intent(in), optional  :: en_index
    integer, intent(in), optional  :: k_index
    integer, intent(in), optional  :: spin
  end subroutine get_sigma_n_mat

  !> Give the Gr at given energy point to the interaction
  subroutine set_Gr(this, Gr, en_index, k_index, spin)
    class(ElPhonDephS) :: this
    type(z_dns), dimension(:,:), intent(in) :: Gr
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    type(z_dns) :: work1, work2, work3
    integer :: npl, ii, nnz

    !! Implement sigma = M*G*M, assuming that PL structure is not only
    !! preserved, but that only the corresponding Gr blocks are used
    !! Now dirty, only working with 1 PL
    npl = this%struct%num_PLs
    if (npl .ne. 1) then
      error stop 'ElphPhonDephB works only with single PL'
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
  subroutine set_Gn(this, Gn, en_index, k_index, spin)
    class(ElPhonDephS) :: this
    type(z_dns), dimension(:,:), intent(in) :: Gn
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    type(z_dns) :: work1, work2, work3
    integer :: npl, ii, nnz
    !! Implement sigma = M*G*M^dagger, assuming that PL structure is not only
    !! preserved, but that only the corresponding Gr blocks are used
    !! Now dirty, only working without PL
    npl = this%struct%num_PLs
    if (npl .ne. 1) then
      error stop 'ElphPhonDephB works only with single PL'
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

  !>  Compute Sigma_r : necessary for inelastic
  subroutine compute_Sigma_r(this, en_index, k_index, spin)
    class(ElPhonDephS) :: this
    integer, intent(in), optional  :: en_index
    integer, intent(in), optional  :: k_index
    integer, intent(in), optional  :: spin
  end subroutine compute_Sigma_r

  !>  Compute Sigma_n : necessary for inelastic
  subroutine compute_Sigma_n(this, en_index, k_index, spin)
    class(ElPhonDephS) :: this
    integer, intent(in), optional  :: en_index
    integer, intent(in), optional  :: k_index
    integer, intent(in), optional  :: spin
  end subroutine compute_Sigma_n

  !> Destroy Sigma_r
  subroutine destroy_Sigma_r(this)
    class(ElPhonDephS) :: this
    integer :: ii
    do ii=1,this%nummodes
      if (allocated(this%sigma_r(ii)%rowpnt)) then
        call destroy(this%sigma_r(ii))
      end if
    end do
    deallocate(this%sigma_r)
  end subroutine destroy_Sigma_r

  !> Destroy Sigma_n
  subroutine destroy_Sigma_n(this)
    class(ElPhonDephS) :: this
    integer :: ii
    do ii=1,this%nummodes
      if (allocated(this%sigma_n(ii)%rowpnt)) then
        call destroy(this%sigma_n(ii))
      end if
    end do
    deallocate(this%sigma_n)
  end subroutine destroy_Sigma_n

  !> Destroy All
  subroutine destroy_all(this)
    class(ElPhonDephS) :: this
    integer :: ii
    do ii=1,this%nummodes
      if (allocated(this%couplings(ii)%rowpnt)) then
        call destroy(this%couplings(ii))
      end if
    end do
    deallocate(this%couplings)
    call destroy_Sigma_r(this)
    call destroy_Sigma_n(this)
  end subroutine destroy_all

end module elphds
