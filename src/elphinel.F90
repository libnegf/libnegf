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
  use ln_constants, only : pi
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
  use inversions, only : inverse
#:if defined("MPI")
  use libmpifx_module
#:endif
  use clock
  implicit none
  private

  public :: ElPhonInel, ElPhonPolarOptical, ElPhonNonPolarOptical
  public :: ElPhonPO_create, ElPhonNonPO_create
  public :: ElPhonPO_init, ElPhonNonPO_init

  type, abstract, extends(TInelastic) :: ElPhonInel
    private
    !> communicator of the cartesian grid
    integer(c_int) :: cart_comm

    !> holds atomic structure: Only the scattering region
    type(TBasisCenters) :: basis

    !> Bose-Einstein phonon occupation
    real(dp) :: Nq

    !> A general scalar coupling
    real(dp) :: coupling

    !> Kpoints
    real(dp), allocatable :: kpoint(:,:)
    real(dp), allocatable :: kweight(:)
    integer, allocatable :: local_kindex(:)
    !> Energy grid. global num of energy points
    integer :: nE_global
    !> Energy grid. Local num of energy points
    integer :: nE_local
    !> Energy grid spacing
    real(dp) :: dE
    !> Matrix KK
    real(dp), allocatable :: Kmat(:,:,:)
    !> binning resolution for z-coordinates
    real(dp) :: dz
    !> Supercell area for integrations
    real(dp) :: cell_area

  contains

    procedure :: add_sigma_r
    procedure :: add_sigma_n
    procedure :: get_sigma_n_blk
    procedure :: get_sigma_n_mat
    procedure :: set_Gr
    procedure :: set_Gn
    procedure :: compute_sigma_r
    procedure :: compute_sigma_n
    procedure :: destroy_sigma_r
    procedure :: destroy_sigma_n
    procedure :: set_kpoints
    procedure :: set_EnGrid
    procedure(abstract_prepare), deferred :: prepare
    procedure :: destroy => destroy_elph

  end type ElPhonInel

  abstract interface
    subroutine abstract_prepare(this)
      import ElPhonInel
      class(ElPhonInel) :: this
    end subroutine abstract_prepare
  end interface


  type, extends(ElPhonInel) :: ElPhonPolarOptical
    private
    !> Parameters for Polar-optical K-function
    real(dp) :: Ce
    real(dp) :: eps0
    real(dp) :: eps_inf
    real(dp) :: q0
    contains
    procedure :: prepare => prepare_POKmat
  end type ElPhonPolarOptical

  type, extends(ElPhonInel) :: ElPhonNonPolarOptical
    private
    !> Parameters for Non-polar-optical K-function
    real(dp) :: D0
    contains
    procedure :: prepare => prepare_NonPOKmat
  end type ElPhonNonPolarOptical

  !type, extends(TInelastic) :: ElPhoton
  !  private
  !  !> Parameters for Polar-optical K-function
  !  real(dp) :: I0
  !  real(dp) :: eps_inf
  !end type ElPhoton

  !> Interface of C function that perform all the MPI communications in order to compute
  !  sigma_r and sigma_n
  interface
   integer (c_int) function self_energy(cart_comm, Nrow, Ncol, NK, NE, NKloc, NEloc, iEshift, &
        & GG, Sigma, sbuff1, sbuff2, rbuff1, rbuff2, rbuffH, sbuffH, fac_min, fac_plus, izr, izc, &
        & KK, Ndz) bind(C, name='self_energy')
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
     !> index of z coordinates on the coarse grid corresponding to the block rows
     integer(c_int) :: izr(*)
     !> index of z coordinates on the coarse grid corresponding to the block cols
     integer(c_int) :: izc(*)
     !> look-up array to store KK(iQ, iK, |zi - zj|)
     real(c_double) :: KK(*)
     !> leading dimension of KK
     integer(c_int), value :: Ndz
   end function self_energy
  end interface

contains

  subroutine ElPhonPO_create(this)
    class(TInteraction), allocatable :: this
    allocate(ElPhonPolarOptical::this)
  end subroutine ElPhonPO_create

  subroutine ElPhonNonPO_create(this)
    class(TInteraction), allocatable :: this
    allocate(ElPhonNonPolarOptical::this)
  end subroutine ElPhonNonPO_create

  !>
  ! Factory for el-ph inelastic model based on polar-optical modes
  ! @param struct : contact/device partitioning infos
  ! @param coupling: coupling (energy units)
  ! @param wq: phonon frequence
  ! @param niter: fixed number of scba iterations
  subroutine ElPhonPO_init(this, comm, struct, basis, coupling, wq, Temp, &
             & dz, eps0, eps_inf, q0, cell_area, niter, tridiag)
    type(ElPhonPolarOptical) :: this
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
    logical, intent(in) :: tridiag

    this%descriptor = &
        & "Electron-Phonon inelastic model for polar-optical phonons"

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
    this%coupling = coupling(1)
    this%tTridiagonal = tridiag

    ! Initialize the cache space
    if (allocated(this%sigma_r)) call this%sigma_r%destroy()
    this%sigma_r = TMatrixCacheMem(tagname='Sigma_r')
    if (allocated(this%sigma_n)) call this%sigma_n%destroy()
    this%sigma_n = TMatrixCacheMem(tagname='Sigma_n')

  end subroutine ElPhonPO_init

  !>
  ! Factory for el-ph inelastic model based on polar-optical modes
  ! @param struct : contact/device partitioning infos
  ! @param coupling: coupling (energy units)
  ! @param wq: phonon frequence
  ! @param niter: fixed number of scba iterations
  subroutine ElPhonNonPO_init(this, comm, struct, basis, coupling, wq, Temp, &
             & dz, D0, cell_area, niter, tridiag)
    type(ElPhonNonPolarOptical) :: this
    integer, intent(in) :: comm
    type(TStruct_Info), intent(in) :: struct
    type(TBasisCenters), intent(in) :: basis
    real(dp), dimension(:), intent(in) :: coupling
    real(dp), intent(in) :: wq
    real(dp), intent(in) :: Temp
    real(dp), intent(in) :: dz
    real(dp), intent(in) :: D0
    real(dp), intent(in) :: cell_area
    integer, intent(in) :: niter
    logical, intent(in) :: tridiag

    this%descriptor = &
        & "Electron-Phonon inelastic model for polar-optical phonons"

    this%cart_comm = comm
    this%scba_niter = niter
    this%struct = struct
    this%basis = basis
    this%wq = wq
    this%Nq = bose(wq, Temp)
    ! Parameters for function Kappa
    this%coupling = coupling(1)
    this%dz = dz
    this%D0 = D0
    this%cell_area = cell_area
    this%tTridiagonal = tridiag

    ! Initialize the cache space
    if (allocated(this%sigma_r)) call this%sigma_r%destroy()
    this%sigma_r = TMatrixCacheMem(tagname='Sigma_r')
    if (allocated(this%sigma_n)) call this%sigma_n%destroy()
    this%sigma_n = TMatrixCacheMem(tagname='Sigma_n')

  end subroutine ElPhonNonPO_init
  !--------------------------------------------------------------------------
  ! This function should be called before the SCBA loop
  subroutine set_kpoints(this, kpoints, kweights, kindex)
    class(ElPhonInel) :: this
    real(dp), intent(in) :: kpoints(:,:)
    real(dp), intent(in) :: kweights(:)
    integer, intent(in) :: kindex(:)

    integer :: nCentralAtoms

    this%kpoint = kpoints
    this%kweight = kweights
    this%local_kindex = kindex

  end subroutine set_kpoints

  !--------------------------------------------------------------------------
  subroutine set_EnGrid(this, deltaE, nE_global, nE_local)
    class(ElPhonInel) :: this
    real(dp), intent(in) :: deltaE
    integer, intent(in) :: nE_global
    integer, intent(in) :: nE_local

    this%dE = deltaE
    this%nE_global = nE_global
    this%nE_local = nE_local

  end subroutine set_EnGrid

  !--------------------------------------------------------------------------
  ! PO phonons: see derivation slides:
  ! Sigma_mn(k,E) = Sum_q,G K(|zm - zn|, k, q+G) [B(-)*Gmn(q,E-wq) + B(+)*Gmn(q,E+wq)]
  subroutine prepare_POKmat(this)
    class(ElPhonPolarOptical) :: this

    integer :: iZ, iQ, iK, fu, nDeltaZ, nCentralAtoms
    real(dp) :: kq(3), kk(3), QQ(3), Q2, bb, z_mn, Kf
    real(dp) :: zmin, zmax, recVecs2p(3,3)
    real(dp), allocatable :: kpoint(:,:)

    ! Compute the matrix Kmat as lookup table
    nCentralAtoms = this%basis%nCentralAtoms

    zmin = minval(this%basis%x(3,1:nCentralAtoms))
    zmax = maxval(this%basis%x(3,1:nCentralAtoms))

    nDeltaZ = nint(zmax - zmin)/this%dz
    call log_allocate(this%Kmat, nDeltaZ+1, size(this%kweight), size(this%kweight))

    this%Ce = this%coupling*this%wq/2.0_dp/this%cell_area* &
          & (1.0_dp/this%eps_inf - 1.0_dp/this%eps0)

    !compute the absolute k-points
    allocate(kpoint(3,size(this%kweight)))
    if (all(this%basis%lattVecs == 0.0_dp)) then
       kpoint = 0.0_dp
    else
       recVecs2p = 0.0_dp
       recVecs2p(1,1) = 2.0_dp*pi/this%basis%lattVecs(1,1)
       recVecs2p(2,2) = 2.0_dp*pi/this%basis%lattVecs(2,2)
       recVecs2p(3,3) = 2.0_dp*pi/this%basis%lattVecs(3,3)
       kpoint = matmul(recVecs2p, this%kpoint)
    end if

    do iQ = 1, size(this%kweight)
      kq = kpoint(:, iQ)
      do iK = 1, size(this%kweight)
        kk = kpoint(:, iK)
        QQ = kk - kq
        Q2 = dot_product(QQ, QQ)
        bb = sqrt(this%q0*this%q0 + Q2)
        !print*,"QQ=",kk-kq
        !print*,"bb=",bb
        do iZ = 0, nDeltaZ
          z_mn = iZ * this%dz
          Kf = (2.0_dp*Q2 + this%q0*this%q0*(1.0_dp-bb*z_mn))*exp(-bb*z_mn)/ (4.0_dp*bb**3)
          this%Kmat(iZ+1,iK,iQ) = this%Ce * this%kweight(iQ) * Kf
        end do
      end do
    end do

    !print*,'debug print Kmat'
    !open(newunit=fu, file='Kmat.dat')
    !do iZ = 0, nDeltaZ
    !  write(fu,*) iZ*this%dz, this%Kmat(iZ+1,1,1)
    !end do
    !close(fu)

  end subroutine prepare_POKmat

  !--------------------------------------------------------------------------
  ! Non PO phonons: see derivation slides:
  ! Sigma_mn(k,E) = Sum_q,G K(|zm - zn|, k, q+G) [B(-)*Gmn(q,E-wq) + B(+)*Gmn(q,E+wq)]
  subroutine prepare_NonPOKmat(this)
    class(ElPhonNonPolarOptical) :: this

    integer :: iZ, iQ, iK, fu, nDeltaZ, nCentralAtoms
    real(dp) :: kq(3), kk(3), QQ, Q2, kq2, kk2, z_mn, Kf
    real(dp) :: zmin, zmax, recVecs2p(3,3)
    real(dp), allocatable :: kpoint(:,:)

    ! Compute the matrix Kmat as lookup table
    nCentralAtoms = this%basis%nCentralAtoms

    zmin = minval(this%basis%x(3,1:nCentralAtoms))
    zmax = maxval(this%basis%x(3,1:nCentralAtoms))

    nDeltaZ = nint(zmax - zmin)/this%dz
    call log_allocate(this%Kmat, nDeltaZ+1, size(this%kweight), size(this%kweight))

    !compute the absolute k-points
    allocate(kpoint(3,size(this%kweight)))
    if (all(this%basis%lattVecs == 0.0_dp)) then
       kpoint = 0.0_dp
    else
       recVecs2p = 0.0_dp
       recVecs2p(1,1) = 2.0_dp*pi/this%basis%lattVecs(1,1)
       recVecs2p(2,2) = 2.0_dp*pi/this%basis%lattVecs(2,2)
       recVecs2p(3,3) = 2.0_dp*pi/this%basis%lattVecs(3,3)
       kpoint = matmul(recVecs2p, this%kpoint)
    end if

    do iQ = 1, size(this%kweight)
      kq = kpoint(:, iQ)
      kq2 = dot_product(kq,kq)
      if (kq2==0.0_dp) then
         this%Kmat(:,:,iQ) = 0.0_dp
         cycle
      end if
      do iK = 1, size(this%kweight)
        kk = kpoint(:, iK)
        kk2 = dot_product(kk,kk)
        if (kk2==0.0_dp) then
           this%Kmat(:,iK,iQ) = 0.0_dp
           cycle
        end if
        QQ = dot_product(kk, kq)
        Q2 = QQ*QQ
        do iZ = 0, nDeltaZ
          z_mn = iZ * this%dz
          Kf = this%D0*Q2*exp(-sqrt(kq2)*z_mn)/(2.0_dp*kq2*kk2)
          !if (isNan(Kf)) then
          !   print*,kq2,kk2,QQ
          !   stop
          !end if
          this%Kmat(iZ+1,iK,iQ) = this%coupling * this%kweight(iQ) * Kf
        end do
      end do
    end do
    !if (any(isNan(this%Kmat))) then
    !  print*,'Kmat= NaN'
    !  stop
    !end if
    !open(newunit=fu, file='Kmat.dat')
    !do iZ = 0, nDeltaZ
    !  write(fu,*) iZ*this%dz, this%Kmat(iZ+1,1,1)
    !end do
    !close(fu)

  end subroutine prepare_NonPOKmat


  !--------------------------------------------------------------------------
  subroutine destroy_elph(this)
    class(ElPhonInel) :: this

    if (allocated(this%kpoint)) deallocate(this%kpoint)
    if (allocated(this%kweight)) deallocate(this%kweight)
    if (allocated(this%local_kindex)) deallocate(this%local_kindex)
    if (allocated(this%Kmat)) call log_deallocate(this%Kmat)
    if (allocated(this%sigma_r)) call this%sigma_r%destroy()
    if (allocated(this%sigma_n)) call this%sigma_n%destroy()

  end subroutine destroy_elph

  !--------------------------------------------------------------------------
  !> Retrieve matrix pointer from cache and add it to ESH
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
      if (this%tTridiagonal) then
        if (jj .lt. npl) then
          label%row_block = jj
          label%col_block = jj + 1
          tmp_blk => this%sigma_r%retrieve_pointer(label)
          ESH(jj, jj + 1)%val = ESH(jj, jj + 1)%val - tmp_blk%val
          label%row_block = jj + 1
          label%col_block = jj
          tmp_blk => this%sigma_r%retrieve_pointer(label)
          ESH(jj + 1, jj)%val = ESH(jj + 1, jj)%val - tmp_blk%val
        end if
      end if
    end do

  end subroutine add_sigma_r

  !--------------------------------------------------------------------------
  !> Retrieve matrix pointer from cache and add it to sigma
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
      ! Gn^dag = Gn => Sigma_n(i,j) = Sigma_n(j,i)^dag
      if (this%tTridiagonal) then
        if (jj .lt. npl) then
          label%row_block = jj
          label%col_block = jj + 1
          tmp_blk => this%sigma_n%retrieve_pointer(label)
          sigma(jj, jj + 1)%val = sigma(jj, jj + 1)%val + tmp_blk%val
          sigma(jj + 1, jj)%val = sigma(jj + 1, jj)%val + conjg(transpose(tmp_blk%val))
        end if
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

#:if defined("MPI")
    integer :: ii, ibl, nbl, Np, Mp, NK, NE, NKloc, NEloc, iEshift, err
    integer :: iK, iE, Ndz, PL_start
    complex(c_double_complex) :: fac_min, fac_plus
    complex(c_double_complex), allocatable :: sbuff1(:,:), rbuff1(:,:)
    complex(c_double_complex), allocatable :: sbuff2(:,:), rbuff2(:,:)
    complex(c_double_complex), allocatable :: sbuffH(:,:), rbuffH(:,:)

    type(C_PTR), allocatable :: pGG(:,:), pSigma(:,:)
    type(z_DNS), target :: Sigma_r
    type(z_DNS), pointer :: pMat
    integer, allocatable :: izr(:), izc(:)
    type(TMatLabel) :: label
    logical :: buff

    real(dp) :: maxvalue

    nbl = this%struct%num_PLs
    NK = size(this%kpoint,2)
    NKloc = size(this%local_kindex)
    NE = this%nE_global
    NEloc = this%nE_local
    iEshift = nint(this%wq/this%dE)
    !if (iEshift > NEloc) then
    !  stop "ERROR: iEshift>NEloc. More points are needed in the energy grid"
    !end if
    Ndz = size(this%Kmat,1)
    label%spin = 0
    if (present(spin)) then
       label%spin = spin
    end if
    allocate(pGG(NEloc,NKloc))
    allocate(pSigma(NEloc,NKloc))

    ! Compute the diagonal blocks
    do ibl = 1, nbl
      ! ==================================================================================================
      !  DIAGONAL BLOCKS
      ! ==================================================================================================
      ! block dimension
      PL_start = this%struct%mat_PL_start(ibl)
      Np = this%struct%mat_PL_end(ibl) - PL_start + 1
      ! create buffers
      call allocate_buff(Np,Np)

      ! Project atom position on the coarser grid
      call log_allocate(izr,Np)
      do ii = 1, Np
        izr(ii) = nint(this%basis%x(3, this%basis%matrixToBasis(PL_start+ii-1))/this%dz)
      end do
      label%row_block = ibl
      label%col_block = ibl

      call setup_pointers_Gr()
      call setup_pointers_Sigma_r(Np,Np)

      ! Compute the retarded part
      fac_min = cmplx(this%Nq + 1, 0.0, dp)
      fac_plus = cmplx(this%Nq, 0.0, dp)

      err = self_energy(this%cart_comm, Np, Np, NK, NE, NKloc, NEloc, iEshift, pGG, pSigma, &
            & sbuff1, sbuff2, rbuff1, rbuff2, rbuffH, sbuffH, fac_min, fac_plus, izr, izr, this%Kmat, Ndz)

      call setup_pointers_Gn()
      ! Compute the Gn part to Sigma_r
      fac_min = (0.0_dp, 0.5_dp)
      fac_plus = (0.0_dp, -0.5_dp)

      err = self_energy(this%cart_comm, Np, Np, NK, NE, NKloc, NEloc, iEshift, pGG, pSigma, &
            & sbuff1, sbuff2, rbuff1, rbuff2, rbuffH, sbuffH, fac_min, fac_plus, izr, izr, this%Kmat, Ndz)

      call deallocate_buff()

      ! ==================================================================================================
      !  UPPER/LOWER TRI-DIAGONAL BLOCKS
      ! ==================================================================================================
      if (this%tTridiagonal .and. ibl < nbl) then
        PL_start = this%struct%mat_PL_start(ibl+1)
        Mp = this%struct%mat_PL_end(ibl+1) - PL_start + 1
        ! Project atom position on the coarser grid (two indep. arrays for rows and cols)
        call log_allocate(izc,Mp)
        do ii = 1, Mp
          izc(ii) = nint(this%basis%x(3, this%basis%matrixToBasis(PL_start+ii-1))/this%dz)
        end do

        ! -------------------- supradiagonal blocks -----------------------------------------------------
        label%row_block = ibl
        label%col_block = ibl + 1

        call allocate_buff(Np,Mp)

        call setup_pointers_Gr()
        call setup_pointers_Sigma_r(Np,Mp)

        ! Compute the retarded part
        fac_min = cmplx(this%Nq + 1, 0.0, dp)
        fac_plus = cmplx(this%Nq, 0.0, dp)

        err = self_energy(this%cart_comm, Np, Mp, NK, NE, NKloc, NEloc, iEshift, pGG, pSigma, &
              & sbuff1, sbuff2, rbuff1, rbuff2, rbuffH, sbuffH, fac_min, fac_plus, izr, izc, this%Kmat, Ndz)

        call setup_pointers_Gn()
        ! Compute the Gn part to Sigma_r
        fac_min = (0.0_dp, 0.5_dp)
        fac_plus = (0.0_dp, -0.5_dp)

        err = self_energy(this%cart_comm, Np, Mp, NK, NE, NKloc, NEloc, iEshift, pGG, pSigma, &
              & sbuff1, sbuff2, rbuff1, rbuff2, rbuffH, sbuffH, fac_min, fac_plus, izr, izc, this%Kmat, Ndz)

        call deallocate_buff()

        ! -------------------- subdiagonal blocks -----------------------------------------------------
        label%row_block = ibl + 1
        label%col_block = ibl

        call allocate_buff(Mp,Np)

        call setup_pointers_Gr()
        call setup_pointers_Sigma_r(Mp,Np)

        ! Compute the retarded part
        fac_min = cmplx(this%Nq + 1, 0.0, dp)
        fac_plus = cmplx(this%Nq, 0.0, dp)

        err = self_energy(this%cart_comm, Mp, Np, NK, NE, NKloc, NEloc, iEshift, pGG, pSigma, &
              & sbuff1, sbuff2, rbuff1, rbuff2, rbuffH, sbuffH, fac_min, fac_plus, izc, izr, this%Kmat, Ndz)

        call setup_pointers_Gn()
        ! Compute the Gn part to Sigma_r
        fac_min = (0.0_dp, 0.5_dp)
        fac_plus = (0.0_dp, -0.5_dp)

        err = self_energy(this%cart_comm, Mp, Np, NK, NE, NKloc, NEloc, iEshift, pGG, pSigma, &
              & sbuff1, sbuff2, rbuff1, rbuff2, rbuffH, sbuffH, fac_min, fac_plus, izc, izr, this%Kmat, Ndz)

        call deallocate_buff()
        call log_deallocate(izc)
      end if
      ! ==================================================================================================

      call log_deallocate(izr)
    end do

    deallocate(pGG)
    deallocate(pSigma)
#:else
    stop "inelastic scattering requires MPI compilation"
#:endif
    contains

    subroutine allocate_buff(Nr,Nc)
      integer, intent(in) :: Nr, Nc

      call log_allocate(sbuff1, Nr, Nc)
      call log_allocate(rbuff1, Nr, Nc)
      call log_allocate(sbuff2, Nr, Nc)
      call log_allocate(rbuff2, Nr, Nc)
      call log_allocate(sbuffH, Nr, Nc)
      call log_allocate(rbuffH, Nr, Nc)
    end subroutine allocate_buff

    subroutine deallocate_buff()
      call log_deallocate(sbuff1)
      call log_deallocate(rbuff1)
      call log_deallocate(sbuff2)
      call log_deallocate(rbuff2)
      call log_deallocate(sbuffH)
      call log_deallocate(rbuffH)
    end subroutine deallocate_buff

    ! setup the array of pointers to G_r and sigma_r
    subroutine setup_pointers_Gr()

      do iK = 1, NKloc
        do iE = 1, NEloc
          ! map local to global index for label
          label%kpoint = this%local_kindex(iK)
          label%energy_point = iE + id*NEloc
          pGG(iE,iK) = this%G_r%retrieve_loc(label)
        end do
      end do
    end subroutine setup_pointers_Gr

    ! setup the array of pointers to G_n
    subroutine setup_pointers_Gn()
      do iK = 1, NKloc
        do iE = 1, NEloc
          label%kpoint = this%local_kindex(iK)
          label%energy_point = iE + id*NEloc
          pGG(iE,iK) = this%G_n%retrieve_loc(label)
        end do
      end do
    end subroutine setup_pointers_Gn


    subroutine setup_pointers_Sigma_r(Nr,Nc)
      integer, intent(in) :: Nr, Nc

      do iK = 1, NKloc
        do iE = 1, NEloc
          label%kpoint = this%local_kindex(iK)
          label%energy_point = iE + id*NEloc
          if (.not.this%Sigma_r%is_cached(label)) then
            !print*,'create Sigma_r', label%kpoint, label%energy_point
            call create(Sigma_r,Nr,Nc)
            Sigma_r%val = (0.0_dp, 0.0_dp)
            call this%Sigma_r%add(Sigma_r, label)
            call destroy(Sigma_r)
          end if
          pSigma(iE,iK) = this%Sigma_r%retrieve_loc(label)
        end do
      end do

    end subroutine setup_pointers_Sigma_r

    subroutine check_elements(Mat,arg)
       class(TMatrixCache) :: Mat
       character(*) :: arg
       do iK = 1, NKloc
        do iE = 1, NEloc
          label%kpoint = this%local_kindex(iK)
          label%energy_point = iE + id*NEloc
          pMat => Mat%retrieve_pointer(label)
          if (any(isNaN(abs(pMat%val)))) then
             print*,arg//'=NaN',iE,iK
          end if
        end do
      end do
    end subroutine check_elements

  end subroutine compute_Sigma_r


  !> Give the Gn at given energy point to the interaction
  subroutine compute_Sigma_n(this, en_index, k_index, spin)
    class(ElPhonInel) :: this
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

#:if defined("MPI")
    integer :: ii, ibl, nbl, Np, Mp, NK, NE, NKloc, NEloc, iEshift, err
    integer :: iK, iE, Ndz, PL_start
    complex(c_double_complex) :: fac_min, fac_plus
    complex(c_double_complex), allocatable :: sbuff1(:,:), rbuff1(:,:)
    complex(c_double_complex), allocatable :: sbuff2(:,:), rbuff2(:,:)
    complex(c_double_complex), allocatable :: sbuffH(:,:), rbuffH(:,:)

    type(C_PTR), allocatable :: pGG(:,:), pSigma(:,:)
    type(z_DNS), pointer :: pMat
    type(z_DNS) :: Sigma_n
    integer, allocatable :: izr(:), izc(:)
    type(TMatLabel) :: label

    real(dp) :: maxvalue
    !print*,'....Compute Sigma_n diagonal blocks'

    nbl = this%struct%num_PLs
    NK = size(this%kpoint,2)
    NKloc = size(this%local_kindex)
    NE = this%nE_global
    NEloc = this%nE_local
    iEshift = nint(this%wq/this%dE)
    !if (iEshift > NEloc) then
    !  stop "ERROR: iEshift>NEloc. More points are needed in the energy grid"
    !end if
    Ndz = size(this%Kmat, 1)
    label%spin = 1
    if (present(spin)) then
      label%spin = spin
    end if
    allocate(pGG(NEloc,NKloc))
    allocate(pSigma(NEloc,NKloc))

    do ibl = 1, nbl
      ! ==================================================================================================
      !  Compute the diagonal blocks
      ! ==================================================================================================
      ! block dimension
      PL_start = this%struct%mat_PL_start(ibl)
      Np = this%struct%mat_PL_end(ibl) - PL_start + 1
      ! create buffers
      call allocate_buff(Np,Np)

      ! Project atom position on the coarser grid
      call log_allocate(izr,Np)
      do ii = 1, Np
        izr(ii) = nint(this%basis%x(3, this%basis%matrixToBasis(PL_start+ii-1))/this%dz)
      end do
      label%row_block = ibl
      label%col_block = ibl

      call setup_pointers_Gn()
      call setup_pointers_Sigma_n(Np,Np)

      ! Compute the retarded part
      fac_min = cmplx(this%Nq, 0.0, dp)
      fac_plus = cmplx(this%Nq + 1, 0.0, dp)

      err = self_energy(this%cart_comm, Np, Np, NK, NE, NKloc, NEloc, iEshift, pGG, pSigma, &
            & sbuff1, sbuff2, rbuff1, rbuff2, rbuffH, sbuffH, fac_min, fac_plus, izr, izr, this%Kmat, Ndz)

      call deallocate_buff()

      ! ==================================================================================================
      !   UPPER TRI- DIAGONAL BLOCKS
      ! ==================================================================================================
      if (this%tTridiagonal .and. ibl < nbl) then
        PL_start = this%struct%mat_PL_start(ibl+1)
        Mp = this%struct%mat_PL_end(ibl+1) - PL_start + 1
        ! Project atom position on the coarser grid (two indep. arrays for rows and cols)
        call log_allocate(izc,Mp)
        do ii = 1, Mp
          izc(ii) = nint(this%basis%x(3, this%basis%matrixToBasis(PL_start+ii-1))/this%dz)
        end do

        ! -------------------- supradiagonal blocks -----------------------------------------------------
        label%row_block = ibl
        label%col_block = ibl + 1
        ! create buffers
        call allocate_buff(Np,Mp)

        call setup_pointers_Gn()
        call setup_pointers_Sigma_n(Np,Mp)
        !call check_elements(this%G_n,"G_n")

        ! Compute the retarded part
        fac_min = cmplx(this%Nq, 0.0, dp)
        fac_plus = cmplx(this%Nq + 1, 0.0, dp)

        err = self_energy(this%cart_comm, Np, Mp, NK, NE, NKloc, NEloc, iEshift, pGG, pSigma, &
              & sbuff1, sbuff2, rbuff1, rbuff2, rbuffH, sbuffH, fac_min, fac_plus, izr, izc, this%Kmat, Ndz)

       ! call check_elements(this%Sigma_n,"Sigma_n")

        call deallocate_buff()
        call log_deallocate(izc)
      end if
      ! ==================================================================================================

      call log_deallocate(izr)
    end do

    deallocate(pGG)
    deallocate(pSigma)
#:else
    stop "inelastic scattering requires MPI compilation"
#:endif

    contains

    subroutine allocate_buff(Nr,Nc)
      integer, intent(in) :: Nr, Nc

      call log_allocate(sbuff1, Nr, Nc)
      call log_allocate(rbuff1, Nr, Nc)
      call log_allocate(sbuff2, Nr, Nc)
      call log_allocate(rbuff2, Nr, Nc)
      call log_allocate(sbuffH, Nr, Nc)
      call log_allocate(rbuffH, Nr, Nc)
    end subroutine allocate_buff

    subroutine deallocate_buff()
      call log_deallocate(sbuff1)
      call log_deallocate(rbuff1)
      call log_deallocate(sbuff2)
      call log_deallocate(rbuff2)
      call log_deallocate(sbuffH)
      call log_deallocate(rbuffH)
    end subroutine deallocate_buff

    ! setup the array of pointers to G_n and sigma_n
    subroutine setup_pointers_Gn()

      do iK = 1, NKloc
        do iE = 1, NEloc
          label%kpoint = this%local_kindex(iK)
          label%energy_point = iE + id*NEloc
          pGG(iE,iK) = this%G_n%retrieve_loc(label)
        end do
      end do
    end subroutine setup_pointers_Gn

    subroutine setup_pointers_Sigma_n(Nr,Nc)
      integer, intent(in) :: Nr, Nc

      do iK = 1, NKloc
        do iE = 1, NEloc
          label%kpoint = this%local_kindex(iK)
          label%energy_point = iE + id*NEloc
          if (.not.this%Sigma_n%is_cached(label)) then
            !print*,'create Sigma_r', label%kpoint, label%energy_point
            call create(Sigma_n,Nr,Nc)
            Sigma_n%val = (0.0_dp, 0.0_dp)
            call this%Sigma_n%add(Sigma_n, label)
            call destroy(Sigma_n)
          end if
          pSigma(iE,iK) = this%Sigma_n%retrieve_loc(label)
        end do
      end do

    end subroutine setup_pointers_Sigma_n

    subroutine check_elements(Mat,arg)
      class(TMatrixCache) :: Mat
      character(*) :: arg
      do iK = 1, NKloc
        do iE = 1, NEloc
          label%kpoint = this%local_kindex(iK)
          label%energy_point = iE + id*NEloc
          pMat => Mat%retrieve_pointer(label)
          if (any(isNaN(abs(pMat%val)))) then
             print*,arg//'=NaN',iE,iK
          end if
        end do
      end do
    end subroutine check_elements

  end subroutine compute_Sigma_n

  !> Destroy Sigma_r
  subroutine destroy_Sigma_r(this)
    class(ElPhonInel) :: this
    if (allocated(this%sigma_r)) call this%sigma_r%destroy()
  end subroutine destroy_Sigma_r

  !> Destroy Sigma_n
  subroutine destroy_Sigma_n(this)
    class(ElPhonInel) :: this
    if (allocated(this%sigma_n)) call this%sigma_n%destroy()
  end subroutine destroy_Sigma_n

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

 ! real(c_double) function kappa(handler, iK, iQ, mu, nu) bind(c, name='kappa')
 !   implicit none
 !   !> handle as c_ptr
 !   type(c_ptr), intent(in), value :: handler
 !   !> Index of k (final k-vector)
 !   integer(c_int), intent(in), value :: iK
 !   !> Index of q (initial k-vector)
 !   integer(c_int), intent(in), value :: iQ
 !   !> Index of mu in the array
 !   integer(c_int), intent(in), value :: mu
 !   !> Index of nu in the array
 !   integer(c_int), intent(in) :: nu

 !   real(dp) :: Ce, kk(3), kq(3), QQ(3), bb, z_mu, z_nu, z_mn, Kf, Q2
 !   integer :: atm_index, loc_iK, loc_iQ
 !   type(ElPhonInel), pointer :: this

 !   !> cast the c_ptr to the instance Telph
 !   call c_f_pointer(handler, this)

 !   ! Forget for the moment Umklapp. G = 0

 !   atm_index = this%basis%matrixToBasis(mu)
 !   z_mu = this%basis%x(3, atm_index)
 !   atm_index = this%basis%matrixToBasis(nu)
 !   z_nu = this%basis%x(3, atm_index)
 !   z_mn = abs(z_mu - z_nu)

 !   kk = this%kpoint(:, iK)
 !   kq = this%kpoint(:, iQ)

 !   QQ = kk - kq
 !   Q2 = dot_product(QQ, QQ)
 !   bb = sqrt(this%q0*this%q0 + Q2)

 !   Kf = (2.0_dp*Q2 + this%q0*this%q0*(1.0_dp-bb*z_mn))*exp(-bb*z_mn)/ (4.0_dp*bb**3)

 !   kappa = this%Ce * this%kweight(iQ) * Kf

 ! end function kappa

end module elphinel
