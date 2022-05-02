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

! Inelastic electron-phonon interactions
! This is the base class for inelastic interactions

module elphinel

  use ln_precision, only : dp
  use interactions, only : TInteraction
  use ln_inelastic, only : TInelastic
  use ln_allocation, only : log_allocate, log_deallocate, writeMemInfo
  use ln_structure, only : TStruct_info
  use ln_structure, only : TBasisCenters, create_TBasis, destroy_TBasis
  use mat_def, only : z_csr, z_dns, create, destroy
  use distributions, only : bose
  use ln_cache
  use iso_c_binding
  use mpi_globals
#:if defined("MPI")
  use libmpifx_module
#:endif

  implicit none
  private

  public :: ElPhonInel, ElPhonInel_create
  public :: ElPhonInel_init, ElPhonInel_setkpoints, ElPhonInel_setEnGrid
  public :: ElPhonInel_destroy

  type, extends(TInelastic) :: ElPhonInel
    private
    !> communicator of the cartesian grid
    integer(c_int) :: cart_comm

    !> holds atomic structure: Only the scattering region
    type(TBasisCenters) :: basis

    !> binning resolution for z-coordinates
    real(dp) :: dz

    !> Bose-Einstein phonon occupation
    real(dp) :: Nq

    !> Kpoints
    real(dp), allocatable :: kpoint(:,:)
    real(dp), allocatable :: kweight(:)
    integer, allocatable :: local_kindex(:)
    real(dp) :: cell_area
    real(dp) :: Ce

    !> Energy grid. global num of energy points
    integer :: nE_global
    !> Energy grid. Local num of energy points
    integer :: nE_local

    !> Paramters for the Kappa function
    real(dp) :: eps0
    real(dp) :: eps_inf
    real(dp) :: q0

    !> Energy grid spacing
    real(dp) :: dE

    !> Matrix KK
    real(dp), allocatable :: Kmat(:,:,:)


  contains

    procedure :: add_sigma_r
    procedure :: add_sigma_n
    procedure :: get_sigma_n_blk
    procedure :: get_sigma_n_mat
    procedure :: set_Gr
    procedure :: set_Gn
    procedure :: compute_sigma_r
    procedure :: compute_sigma_n

  end type ElPhonInel

  !> Interface of C function that perform all the MPI communications in order to compute
  !  sigma_r and sigma_n
  interface
   integer (c_int) function self_energy(cart_comm, Nrow, Ncol, NK, NE, NKloc, NEloc, iEshift, &
        & GG, Sigma, sbuff1, sbuff2, rbuff1, rbuff2, rbuffH, sbuffH, fac_min, fac_plus, iz, KK, Ndz) &
        & bind(C, name='self_energy')
     use iso_c_binding
     !> Cartesian MPI_communicator (negf%cartgrid%id)
     integer(c_int) :: cart_comm
     !> matrix size Nrow x Ncol
     integer(c_int), value :: Nrow
     !> matrix size Nrow x Ncol
     integer(c_int), value :: Ncol
     !> Total Number of K-points
     integer(c_int), value :: NK
     !> Total Number of E-points
     integer(c_int), value :: NE
     !> Total Number of local K-points
     integer(c_int), value :: NKloc
     !> Total Number of local E-points
     integer(c_int), value :: NEloc
     !> grid shift of hwq  iEshif
     integer(c_int), value :: iEshift
     !> Green Function
     type(c_ptr) :: GG
     !> self energy
     type(c_ptr) :: Sigma
     !> 6 Buffers of size Np x Np
     complex(c_double_complex) :: sbuff1(*)
     complex(c_double_complex) :: sbuff2(*)
     complex(c_double_complex) :: rbuff1(*)
     complex(c_double_complex) :: rbuff2(*)
     complex(c_double_complex) :: sbuffH(*)
     complex(c_double_complex) :: rbuffH(*)
     !> prefactor of GG(E-wq) (e.g. nq+1, nq, i/2)
     complex(c_double_complex), value :: fac_min
     !> prefactor of GG(E+wq)
     complex(c_double_complex), value :: fac_plus
     !> index of z coordinates on the coarse grid
     integer(c_int) :: iz(*)
     !> look-up array to store KK(iQ, iK, |zi - zj|)
     real(c_double) :: KK(*)
     !> leading dimension of KK
     integer(c_int), value :: Ndz
   end function self_energy
  end interface

contains

  subroutine ElPhonInel_create(this)
    class(TInteraction), allocatable :: this
    allocate(ElPhonInel::this)
  end subroutine ElPhonInel_create

  !>
  ! Factory for el-ph inelastic model based on polar-optical modes
  ! @param struct : contact/device partitioning infos
  ! @param coupling: coupling (energy units)
  ! @param wq: phonon frequence
  ! @param niter: fixed number of scba iterations
  subroutine ElPhonInel_init(this, comm, struct, basis, coupling, wq, Temp, &
             & dz, eps0, eps_inf, q0, cell_area, niter)
    type(ElPhonInel) :: this
    integer, intent(in) :: comm
    type(TStruct_Info), intent(in) :: struct
    type(TBasisCenters), intent(in) :: basis
    real(dp), dimension(:), intent(in) :: coupling
    real(dp), intent(in) :: wq
    real(dp), intent(in) :: Temp
    real(dp), intent(in) :: dz
    real(dp), intent(in) :: eps0
    real(dp), intent(in) :: eps_inf
    real(dp), intent(in) :: q0
    real(dp), intent(in) :: cell_area
    integer, intent(in) :: niter

    this%descriptor = &
        & "Electron-Phonon inelastic model for optical phonons"

    this%cart_comm = comm
    this%scba_niter = niter
    this%struct = struct
    this%basis = basis
    this%wq = wq
    this%Nq = bose(wq, Temp)
    ! Parameters for function Kappa
    this%dz = dz
    this%eps0 = eps0
    this%eps_inf = eps_inf
    this%q0 = q0
    this%cell_area = cell_area
    this%Ce = coupling(1)*this%wq/2.0_dp/this%cell_area* &
          & (1.0_dp/this%eps_inf - 1.0_dp/this%eps0)

    ! Initialize the cache space
    if (allocated(this%sigma_r)) call this%sigma_r%destroy()
    this%sigma_r = TMatrixCacheMem(tagname='Sigma_r')
    if (allocated(this%sigma_n)) call this%sigma_n%destroy()
    this%sigma_n = TMatrixCacheMem(tagname='Sigma_n')

  end subroutine ElPhonInel_init

  !--------------------------------------------------------------------------
  ! This function should be called before the SCBA loop
  subroutine ElPhonInel_setkpoints(this, kpoints, kweights, kindex)
    type(ElPhonInel), intent(inout) :: this
    real(dp), intent(in) :: kpoints(:,:)
    real(dp), intent(in) :: kweights(:)
    integer, intent(in) :: kindex(:)

    integer :: nDeltaZ, nCentralAtoms, iZ, iQ, iK, fu
    real(dp) :: kq(3), kk(3), QQ(3), Q2, bb, z_mn
    real(dp) :: zmin, zmax, Kf

    this%kpoint = kpoints
    this%kweight = kweights
    this%local_kindex = kindex

    ! Compute the matrix Kmat as lookup table
    nCentralAtoms = this%basis%nCentralAtoms

    zmin = minval(this%basis%x(3,1:nCentralAtoms))
    zmax = maxval(this%basis%x(3,1:nCentralAtoms))

    nDeltaZ = nint(zmax - zmin)/this%dz

    call log_allocate(this%Kmat, nDeltaZ+1, size(kweights), size(kweights))

    do iQ = 1, size(kweights)
      kq = this%kpoint(:, iQ)
      do iK = 1, size(kweights)
        kk = this%kpoint(:, iK)
        QQ = kk - kq
        Q2 = dot_product(QQ, QQ)
        bb = sqrt(this%q0*this%q0 + Q2)
        do iZ = 0, nDeltaZ
          z_mn = iZ * this%dz
          Kf = (2.0_dp*Q2 + this%q0*this%q0*(1.0_dp-bb*z_mn))*exp(-bb*z_mn)/ (4.0_dp*bb**3)
          this%Kmat(iZ+1,iK,iQ) = this%Ce * this%kweight(iQ) * Kf
        end do
      end do
    end do
    !open(newunit=fu, file='Kmat.dat')
    !do iZ = 0, nDeltaZ
    !  write(fu,*) iZ*this%dz, this%Kmat(iZ+1,1,1)
    !end do
    !close(fu)
  end subroutine ElPhonInel_setkpoints

  !--------------------------------------------------------------------------
  subroutine ElPhonInel_setEnGrid(this, deltaE, nE_global, nE_local)
    type(ElPhonInel), intent(inout) :: this
    real(dp), intent(in) :: deltaE
    integer, intent(in) :: nE_global
    integer, intent(in) :: nE_local

    this%dE = deltaE
    this%nE_global = nE_global
    this%nE_local = nE_local

  end subroutine ElPhonInel_SetEnGrid

  !--------------------------------------------------------------------------
  subroutine ElPhonInel_destroy(this)
    type(ElPhonInel), intent(inout) :: this

    if (allocated(this%kpoint)) deallocate(this%kpoint)
    if (allocated(this%kweight)) deallocate(this%kweight)
    if (allocated(this%local_kindex)) deallocate(this%local_kindex)
    if (allocated(this%Kmat)) call log_deallocate(this%Kmat)
    if (allocated(this%sigma_r)) call this%sigma_r%destroy()
    if (allocated(this%sigma_n)) call this%sigma_n%destroy()

  end subroutine ElPhonInel_destroy

  !--------------------------------------------------------------------------
  !> This interface should append
  !  the retarded self energy to ESH
  subroutine add_sigma_r(this, esh, en_index, k_index, spin)
    class(ElPhonInel) :: this
    type(z_dns), dimension(:,:), intent(inout) :: esh
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    integer :: n, npl, ii, ierr, jj
    type(z_dns), pointer :: tmp_blk
    type(TMatLabel) :: label
    !print*,'inel%add_sigma_r'
    npl = this%struct%num_PLs

    if (this%scba_iter .eq. 0) return

    label%kpoint = k_index
    label%energy_point = en_index
    label%spin = 1

    do jj = 1, npl
      label%row_block = jj
      label%col_block = jj
      !call print_label(label)
      tmp_blk => this%sigma_r%retrieve_pointer(label)
      ESH(jj, jj)%val = ESH(jj, jj)%val - tmp_blk%val
      ! 3-diagonal blocks
      if (jj .lt. npl) then
        label%row_block = jj
        label%col_block = jj + 1
        tmp_blk => this%sigma_r%retrieve_pointer(label)
        ESH(jj, jj + 1)%val = ESH(jj, jj + 1)%val - tmp_blk%val
        ESH(jj + 1, jj)%val = ESH(jj + 1, jj)%val - tmp_blk%val
      end if
    end do

  end subroutine add_sigma_r

  !--------------------------------------------------------------------------
  !> This interface should append
  !  sigma_n to a passed self energy, sigma
  subroutine add_sigma_n(this, sigma, en_index, k_index, spin)
    class(ElPhonInel) :: this
    type(z_dns), dimension(:,:), intent(inout) :: sigma
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    integer :: n, npl, ii, ierr, jj
    type(z_dns), pointer :: tmp_blk
    type(TMatLabel) :: label
    !print*,'inel%add_sigma_n'
    npl = this%struct%num_PLs

    if (this%scba_iter .eq. 0) return

    label%kpoint = k_index
    label%energy_point = en_index
    label%spin = 1

    do jj = 1, npl
      label%row_block = jj
      label%col_block = jj
      !call print_label(label)
      tmp_blk => this%sigma_n%retrieve_pointer(label)
      sigma(jj, jj)%val = sigma(jj, jj)%val + tmp_blk%val
      ! 3-diagonal blocks
      if (jj .lt. npl) then
        label%row_block = jj
        label%col_block = jj + 1
        tmp_blk => this%sigma_r%retrieve_pointer(label)
        sigma(jj, jj + 1)%val = sigma(jj, jj + 1)%val - tmp_blk%val
        sigma(jj + 1, jj)%val = sigma(jj + 1, jj)%val - tmp_blk%val
      end if
    end do

  end subroutine add_sigma_n

  !--------------------------------------------------------------------------
  !> Returns the lesser (n) Self Energy in block format
  !
  subroutine get_sigma_n_blk(this, blk_sigma_n, en_index, k_index, spin)
    class(ElPhonInel) :: this
    type(z_dns), dimension(:,:), intent(inout) :: blk_sigma_n
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    type(z_dns), pointer :: tmp_blk
    type(TMatLabel) :: label
    integer :: ii, jj

    if (this%scba_iter .eq. 0) return

    label%kpoint = 0
    label%energy_point = 0
    label%spin = 0
    if (present(k_index)) then
      label%kpoint = k_index
    end if
    if (present(en_index)) then
      label%energy_point = en_index
    end if
    if (present(spin)) then
      label%spin = spin
    end if

    ! Retrieve the diagonal blocks
    do ii = 1, size(blk_sigma_n, 1)
      label%row_block = ii
      label%col_block = ii
      call this%sigma_n%retrieve(blk_sigma_n(ii,ii), label)
    end do

  end subroutine get_sigma_n_blk

  !--------------------------------------------------------------------------
  !> Returns the lesser (n) Self Energy for a given block
  !
  subroutine get_sigma_n_mat(this, sigma_n, ii, jj, en_index, k_index, spin)
    class(ElPhonInel) :: this
    type(z_dns), intent(inout) :: sigma_n
    integer, intent(in) :: ii
    integer, intent(in) :: jj
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    type(TMatLabel) :: label

    if (this%scba_iter .eq. 0) return

    label%row_block = ii
    label%col_block = jj
    label%kpoint = 0
    label%energy_point = 0
    label%spin = 0

    if (present(k_index)) then
      label%kpoint = k_index
    end if
    if (present(en_index)) then
      label%energy_point = en_index
    end if
    if (present(spin)) then
      label%spin = spin
    end if

    call this%sigma_n%retrieve(sigma_n, label)

  end subroutine get_sigma_n_mat


  !> Dummy subroutine to set Gr
  subroutine set_Gr(this, Gr, en_index, k_index, spin)
    class(ElPhonInel) :: this
    type(z_dns), dimension(:,:), intent(in) :: Gr
    integer, intent(in), optional  :: en_index
    integer, intent(in), optional  :: k_index
    integer, intent(in), optional  :: spin
  end subroutine set_Gr

  !> Dummy subroutine to set Gn
  subroutine set_Gn(this, Gn, en_index, k_index, spin)
    class(ElPhonInel) :: this
    type(z_dns), dimension(:,:), intent(in) :: Gn
    integer, intent(in), optional  :: en_index
    integer, intent(in), optional  :: k_index
    integer, intent(in), optional  :: spin
  end subroutine set_Gn


  !--------------------------------------------------------------------------
  !> Give the Gn at given energy point to the interaction
  !>
  subroutine compute_Sigma_r(this, en_index, k_index, spin)
    class(ElPhonInel) :: this
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    integer :: ii, ibl, nbl, Np, NK, NE, NKloc, NEloc, iEshift, err
    integer :: iK, iE, Ndz, PL_start
    complex(c_double_complex) :: fac_min, fac_plus
    complex(c_double_complex), allocatable :: sbuff1(:,:), rbuff1(:,:)
    complex(c_double_complex), allocatable :: sbuff2(:,:), rbuff2(:,:)
    complex(c_double_complex), allocatable :: sbuffH(:,:), rbuffH(:,:)

    type(C_PTR), allocatable :: pGG(:,:), pSigma(:,:)
    type(z_DNS), target :: Sigma_r
    type(z_DNS), pointer :: pMat
    integer, allocatable :: iz(:)
    type(TMatLabel) :: label
    logical :: buff

    print*,'....Compute Sigma_r diagonal blocks'

    nbl = this%struct%num_PLs
    NK = size(this%kpoint,2)
    NKloc = size(this%local_kindex)
    NE = this%nE_global
    NEloc = this%nE_local
    iEshift = nint(this%wq/this%dE)
    Ndz = size(this%Kmat,1)
    label%spin = 0
    if (present(spin)) then
       label%spin = spin
    end if
    allocate(pGG(NEloc,NKloc))
    allocate(pSigma(NEloc,NKloc))

    ! Compute the diagonal blocks
    do ibl = 1, nbl
      ! block dimension
      PL_start = this%struct%mat_PL_start(ibl)
      Np = this%struct%mat_PL_end(ibl) - PL_start + 1

      ! create buffers
      call log_allocate(sbuff1, Np, Np)
      call log_allocate(rbuff1, Np, Np)
      call log_allocate(sbuff2, Np, Np)
      call log_allocate(rbuff2, Np, Np)
      call log_allocate(sbuffH, Np, Np)
      call log_allocate(rbuffH, Np, Np)

      ! Project atom position on the coarser grid
      call log_allocate(iz,Np)
      do ii = 1, Np
        iz(ii) = nint(this%basis%x(3, this%basis%matrixToBasis(PL_start+ii-1))/this%dz)
      end do
      label%row_block = ibl
      label%col_block = ibl

      ! serialize printing
      !if (id > 0) then
      !  call mpifx_recv(debug_comm, buff)
      !end if

      ! setup the array of pointers to G_r and sigma_r
      do iK = 1, NKloc
        do iE = 1, NEloc
          ! map local to global index for label
          label%kpoint = this%local_kindex(iK)
          label%energy_point = iE + id*NEloc
          !call print_label(label)
          pGG(iE,iK) = this%G_r%retrieve_loc(label)
          if (.not.this%Sigma_r%is_cached(label)) then
            !print*,'create Sigma_r', iK, iE
            call create(Sigma_r,Np,Np)
            Sigma_r%val = (0.0_dp, 0.0_dp)
            call this%Sigma_r%add(Sigma_r, label)
            call destroy(Sigma_r)
          end if
          pSigma(iE,iK) = this%Sigma_r%retrieve_loc(label)
        end do
      end do

      !if ( id /= numprocs-1 ) then
      !   call mpifx_send(debug_comm, buff, id+1)
      !end if
      !call mpifx_barrier(debug_comm)

      ! Compute the retarded part
      fac_min = cmplx(this%Nq + 1, 0.0, dp)
      fac_plus = cmplx(this%Nq, 0.0, dp)

      err = self_energy(this%cart_comm, Np, Np, NK, NE, NKloc, NEloc, iEshift, pGG, pSigma, &
            & sbuff1, sbuff2, rbuff1, rbuff2, rbuffH, sbuffH, fac_min, fac_plus, iz, this%Kmat, Ndz)

      ! setup the array of pointers to G_n
      do iK = 1, NKloc
        do iE = 1, NEloc
          label%kpoint = this%local_kindex(iK)
          label%energy_point = iE + id*NEloc
          !print*,'retrieve address of Gn', label%kpoint, label%energy_point, ibl, ibl
          pGG(iE,iK) = this%G_n%retrieve_loc(label)
        end do
      end do

      ! Compute the Gn part to Sigma_r
      fac_min = (0.0_dp, 0.5_dp)
      fac_plus = (0.0_dp, -0.5_dp)

      err = self_energy(this%cart_comm, Np, Np, NK, NE, NKloc, NEloc, iEshift, pGG, pSigma, &
            & sbuff1, sbuff2, rbuff1, rbuff2, rbuffH, sbuffH, fac_min, fac_plus, iz, this%Kmat, Ndz)

      call log_deallocate(sbuff1)
      call log_deallocate(rbuff1)
      call log_deallocate(sbuff2)
      call log_deallocate(rbuff2)
      call log_deallocate(sbuffH)
      call log_deallocate(rbuffH)
      call log_deallocate(iz)

    end do

    deallocate(pGG)
    deallocate(pSigma)

  end subroutine compute_Sigma_r


  !> Give the Gn at given energy point to the interaction
  subroutine compute_Sigma_n(this, en_index, k_index, spin)
    class(ElPhonInel) :: this
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    integer :: ii, ibl, nbl, Np, NK, NE, NKloc, NEloc, iEshift, err
    integer :: iK, iE, Ndz, PL_start
    complex(c_double_complex) :: fac_min, fac_plus
    complex(c_double_complex), allocatable :: sbuff1(:,:), rbuff1(:,:)
    complex(c_double_complex), allocatable :: sbuff2(:,:), rbuff2(:,:)
    complex(c_double_complex), allocatable :: sbuffH(:,:), rbuffH(:,:)

    type(C_PTR), allocatable :: pGG(:,:), pSigma(:,:)
    type(z_DNS), pointer :: pMat
    type(z_DNS) :: Sigma_n
    integer, allocatable :: iz(:)
    type(TMatLabel) :: label

    print*,'....Compute Sigma_n diagonal blocks'

    nbl = this%struct%num_PLs
    NK = size(this%kpoint,2)
    NKloc = size(this%local_kindex)
    NE = this%nE_global
    NEloc = this%nE_local
    iEshift = nint(this%wq/this%dE)
    Ndz = size(this%Kmat, 1)
    label%spin = 1
    if (present(spin)) then
      label%spin = spin
    end if
    allocate(pGG(NEloc,NKloc))
    allocate(pSigma(NEloc,NKloc))

    ! Compute the diagonal blocks
    do ibl = 1, nbl
      ! block dimension
      PL_start = this%struct%mat_PL_start(ibl)
      Np = this%struct%mat_PL_end(ibl) - PL_start + 1
      ! create buffers
      call log_allocate(sbuff1, Np, Np)
      call log_allocate(rbuff1, Np, Np)
      call log_allocate(sbuff2, Np, Np)
      call log_allocate(rbuff2, Np, Np)
      call log_allocate(sbuffH, Np, Np)
      call log_allocate(rbuffH, Np, Np)

      ! Project atom position on the coarser grid
      call log_allocate(iz,Np)
      do ii = 1, Np
        iz(ii) = nint(this%basis%x(3, this%basis%matrixToBasis(PL_start+ii-1))/this%dz)
      end do
      label%row_block = ibl
      label%col_block = ibl

      ! setup the array of pointers to G_n and sigma_n
      do iK = 1, NKloc
        do iE = 1, NEloc

          label%kpoint = this%local_kindex(iK)
          label%energy_point = iE + id*NEloc

          !call print_label(label)
          pGG(iE,iK) = this%G_n%retrieve_loc(label)
          if (.not.this%Sigma_n%is_cached(label)) then
            ! print*,'create Sigma_n', iK, iE
            call create(Sigma_n,Np,Np)
            Sigma_n%val = (0.0_dp, 0.0_dp)
            call this%Sigma_n%add(Sigma_n, label)
            call destroy(Sigma_n)
          end if
          pSigma(iE,iK) = this%Sigma_n%retrieve_loc(label)
        end do
      end do

      ! Compute the retarded part
      fac_min = cmplx(this%Nq, 0.0, dp)
      fac_plus = cmplx(this%Nq + 1, 0.0, dp)

      err = self_energy(this%cart_comm, Np, Np, NK, NE, NKloc, NEloc, iEshift, pGG, pSigma, &
            & sbuff1, sbuff2, rbuff1, rbuff2, rbuffH, sbuffH, fac_min, fac_plus, iz, this%Kmat, Ndz)

      call log_deallocate(sbuff1)
      call log_deallocate(rbuff1)
      call log_deallocate(sbuff2)
      call log_deallocate(rbuff2)
      call log_deallocate(sbuffH)
      call log_deallocate(rbuffH)
      call log_deallocate(iz)

    end do

    deallocate(pGG)
    deallocate(pSigma)

  end subroutine compute_Sigma_n


  !> Integral over qz of the electron-phonon polar optical couping.
  !
  !  Kappa = Sum_G [ w(q) Ce  I(k, q+G, |z_mu - z_nu|) ]
  !
  !  Ce = e^2/eps0 hbar wq/2 (eps0/eps_inf - 1)
  !  Q = k - q
  !  b = sqrt(q0^2 + Q^2)
  !  I = (2Q^2 + q0^2 (1-b |z_mu - z_nu|)) exp(-b |z_mu - z_nu|)/ 4 b^3
  !
  ! |U|^2 = En*L En * L^2   <= OK
  ! Ce ~ En^2 L;  b ~ 1/L;  I ~ L;  =>  kappa ~ En^2 L^2
  !
  ! =>> Must divide by the Area of the cell
  !
  ! LOOKUP TABLE: KK(|zi - zj|, iQ, iK)
  !
  ! NK ~ 256 => NK*NK = 65563
  ! Natom = 1000 * 10 (npl) = 10,000
  ! Option to approximately map zi's onto a regular grid:
  ! Essentially setup a binning parameter. e.g. Si unit cell is 5.4 AA.
  ! Can be pre-computed for every block just beforehand
  ! Could place a resolution of 0.01 AA => 540 points
  ! => 65563 * 1000 * 8 ~ 500 Mb
  ! We need a mapping: matrix index(mu) -> bin(iz)

  real(c_double) function kappa(handler, iK, iQ, mu, nu) bind(c, name='kappa')
    implicit none
    !> handle as c_ptr
    type(c_ptr), intent(in), value :: handler
    !> Index of k (final k-vector)
    integer(c_int), intent(in), value :: iK
    !> Index of q (initial k-vector)
    integer(c_int), intent(in), value :: iQ
    !> Index of mu in the array
    integer(c_int), intent(in), value :: mu
    !> Index of nu in the array
    integer(c_int), intent(in) :: nu

    real(dp) :: Ce, kk(3), kq(3), QQ(3), bb, z_mu, z_nu, z_mn, Kf, Q2
    integer :: atm_index, loc_iK, loc_iQ
    type(ElPhonInel), pointer :: this

    !> cast the c_ptr to the instance Telph
    call c_f_pointer(handler, this)

    ! Forget for the moment Umklapp. G = 0

    atm_index = this%basis%matrixToBasis(mu)
    z_mu = this%basis%x(3, atm_index)
    atm_index = this%basis%matrixToBasis(nu)
    z_nu = this%basis%x(3, atm_index)
    z_mn = abs(z_mu - z_nu)

    kk = this%kpoint(:, iK)
    kq = this%kpoint(:, iQ)

    QQ = kk - kq
    Q2 = dot_product(QQ, QQ)
    bb = sqrt(this%q0*this%q0 + Q2)

    Kf = (2.0_dp*Q2 + this%q0*this%q0*(1.0_dp-bb*z_mn))*exp(-bb*z_mn)/ (4.0_dp*bb**3)

    kappa = this%Ce * this%kweight(iQ) * Kf

  end function kappa

end module elphinel
