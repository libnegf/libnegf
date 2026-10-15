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

! Inelastic electron-photon interaction
#:include "types.fypp"

module elptinel
  use ieee_arithmetic, only : isnan => ieee_is_nan
  use ln_precision, only : sp, dp, lp => dp
  use mat_def, only : c_csr, z_csr, c_dns, z_dns, x_dns => z_dns, x_csr => z_csr
  use mat_def, only : create, destroy, print_mat, assignment(=)
  use block_3_diag, only : init_tridiag_blk, destroy_tridiag_blk, blk2dns
  use block_3_diag, only : csr2blk3diag
  use ln_math, only : cross_product, volume
  use sparsekit_drv, only : prealloc_mult, getdiag
  use ln_constants, only : c, esq_ov_eps0, alpha, pi
  use interactions, only : TInteraction
  use ln_inelastic, only : TInelastic
  use ln_allocation, only : log_allocate, log_deallocate, writeMemInfo
  use ln_structure, only : TStruct_info
  use ln_structure, only : TBasisCenters, create_TBasis, destroy_TBasis
  use ln_cache
  use ln_types, only : TMatrixArrayHS
  use ln_enums, only : interaction_models
  use iso_c_binding
  use mpi_f08, MPI_XX_COMPLEX => MPI_DOUBLE_COMPLEX
  use mpi_globals
  use clock
  use ln_messages, only : error_msg
  implicit none
  private

  public :: ElPhoton
  public :: ElPhoton_create
  public :: ElPhoton_init
  public :: set_Pmat

  interface set_Pmat
    module procedure :: set_Pmat_array
    module procedure :: set_Pmat_struct
  end interface set_Pmat

  type, public :: TMatPointer
    type(x_DNS), pointer :: pMat
  end type TMatPointer

  type, extends(TInelastic) :: ElPhoton
    private
    !> communicator of the cartesian grid
    type(MPI_Comm) :: cart_comm
    !> holds atomic structure: Only the scattering region
    type(TBasisCenters) :: basis

    !> Number of photons
    real(lp) :: Nq

    !> A general scalar coupling
    real(lp) :: coupling

    !> The refraction index
    real(lp) :: nr

    !> Supercell volume for integrations in units of nm
    real(lp) :: cell_volume

    !> The first dimension of Pbkl and total number of Kpoints
    integer :: nK

    !> Polarization direction
    integer::poldir

    !> prefactor (squared) for the el-photon coupling
    real(lp) :: prefactor

    !> The polarization matrices (R_i - R_j) * H_ij in units of nm * ev
    !Dimension: nK, 3 polarazations
    type(x_CSR), allocatable, dimension(:, :) :: Pcsr

    !> Whether the Pcsr matrix has been set
    logical :: tPcsr_set = .false.

    !> Energy grid. global num of energy points
    integer :: nE_global
    !> Energy grid. Local num of energy points
    integer :: nE_local
    !> Energy grid spacing
    real(lp) :: dE
    !> Number of frequency points for the blackbody spectrum
    integer :: nfreq
    !> Blackbody temperature in energy units (eV)
    real(lp) :: kbT
    !> Whether spontaneous emission is included in the calculation
    logical :: tSpontEmission

    integer :: count_iter
    !> Model type (laser, blackbody, debug)
    integer :: model_type

  contains

    procedure :: add_sigma_r_sp
    procedure :: add_sigma_r_dp
    procedure :: add_sigma_n_sp
    procedure :: add_sigma_n_dp
    procedure :: get_sigma_n_blk
    procedure :: get_Gn_blk
    procedure :: get_sigma_n_mat
    procedure :: set_Gr
    procedure :: set_Gn
    procedure :: init_sigma_r
    procedure :: init_sigma_n
    procedure :: compute_sigma_r
    procedure :: compute_sigma_n
    procedure :: destroy_sigma_r
    procedure :: destroy_sigma_n
    procedure :: clear_sigma_r
    procedure :: clear_sigma_n
    procedure :: set_EnGrid
    !procedure :: set_Pmat => set_Pmat_i
    procedure :: is_Pmat_set
    procedure :: prepare => prepare_Kmat
    procedure :: destroy => destroy_elphot

  end type ElPhoton

contains

  subroutine ElPhoton_create(this)
    class(TInteraction), allocatable :: this
    allocate(ElPhoton::this)
  end subroutine ElPhoton_create

  !>
  ! Factory for el-phot inelastic model
  ! @param coupling: coupling (energy units)
  subroutine ElPhoton_init(this, comm, basis, coupling, cell_vol, nr, Ephoton, &
              & intensity, poldir, niter, kbT, nfreq, tSpontEmission, tridiag, modeltype)
    type(ElPhoton) :: this
    type(MPI_Comm), intent(in) :: comm
    type(TBasisCenters), intent(in) :: basis
    real(lp), dimension(:), intent(in) :: coupling
    real(lp), intent(in) :: cell_vol
    real(lp), intent(in) :: nr
    real(lp), intent(in) :: Ephoton
    real(lp), intent(in) :: intensity
    integer, intent(in) :: poldir
    integer, intent(in) :: niter
    real(lp), intent(in) :: kbT
    integer, intent(in) :: nfreq
    logical, intent(in) :: tSpontEmission
    logical, intent(in) :: tridiag
    integer, intent(in) :: modeltype

    integer :: i
    real(lp) :: c__au

    this%model_type = modeltype
    this%count_iter = 0

    this%descriptor = &
        & "Electron-Photon inelastic model"

    this%cart_comm = comm
    this%scba_niter = niter
    this%basis = basis

    ! Nq = I_au*Vol/(Ephoton*c_au)
    ! This is the number of photons in the volume Vol, given the intensity I_au and photon energy Ephoton
    c__au = 1.0_lp / alpha
    this%Nq = intensity * cell_vol / (Ephoton * c__au)
    this%nr = nr
    this%wq = Ephoton
    this%cell_volume = cell_vol
    ! Temperature of blackbody radiation in energy units (eV)
    this%kbT = kbT
    this%nfreq = nfreq

    this%poldir = poldir
    this%coupling = coupling(1)
    this%tTridiagonal = tridiag
    this%tSpontEmission = tSpontEmission

    ! Initialize the cache space
    if (associated(this%sigma_r)) then
      call error_msg('Initialization of Photon called twice')
    end if
    allocate(TMatrixCacheMem::this%sigma_r)
    select type(p => this%sigma_r)
    type is(TMatrixCacheMem)
      p%tagname = 'Sigma_r'
    end select

    if (associated(this%sigma_n)) then
      call error_msg('Initialization of Photon called twice')
    end if
    allocate(TMatrixCacheMem::this%sigma_n)
    select type(p => this%sigma_n)
    type is(TMatrixCacheMem)
      p%tagname = 'Sigma_n'
    end select

    ! electron-photon coupling factor
    ! M_uv = i sqrt[ esq / (2 eps0 hbar w nr Vol) ] xi * (R_u - R_v) H_uv exp[i k (R_u - R_v)]
    !
    ! We later define P_uv = i (R_u - R_v) H_uv, and the prefactor is then
    !
    ! this%prefactor = esq_ov_eps0/(2.0 * cell_volume * nr * hbar wq)
    ! In atomic units, esq/4 pi eps0 = 1.0, wq is already in a.u.
    !
    this%prefactor = this%coupling * 2.0_dp*pi/(this%cell_volume * this%nr * this%wq)

    if (id0) then
      print*, "DEBUG: Initializing ElPhoton inelastic model"
      print*, "DEBUG: coupling = ", this%coupling
      print*, "DEBUG: I = ", intensity
      print*, "DEBUG: prefactor = ", this%prefactor
      print*, "DEBUG: Ephoton = ", this%wq
      print*, "DEBUG: Volume = ", this%cell_volume
      print*, "DEBUG: Nq = ", this%Nq
      !print*, "DEBUG: Nk = ", this%nk
      !print*, "DEBUG: Nbl = ", nbl
      print*, "DEBUG: tridiag = ", this%tTridiagonal
    end if

  end subroutine ElPhoton_init

  !--------------------------------------------------------------------------
  subroutine set_EnGrid(this, deltaE, nE_global, nE_local)
    class(ElPhoton) :: this
    real(lp), intent(in) :: deltaE
    integer, intent(in) :: nE_global
    integer, intent(in) :: nE_local

    this%dE = deltaE
    this%nE_global = nE_global
    this%nE_local = nE_local

  end subroutine set_EnGrid
  !--------------------------------------------------------------------------

  subroutine set_Pmat_array(this, Pcsr)
    type(ElPhoton) :: this
    type(z_CSR), intent(in) :: Pcsr(:,:)


    this%nk = size(Pcsr, 1)
    if (size(Pcsr, 2) .ne. 3) then
      call error_msg('set_Pmat_array: Pcsr must have 3 columns for polarizations')
    end if
    if (.not.allocated(this%Pcsr)) then
      allocate(this%Pcsr(this%nk, 3))
    else
      if (size(this%Pcsr) .ne. this%nk) then
         deallocate(this%Pcsr)
         allocate(this%Pcsr(this%nk, 3))
      end if
    end if

    this%Pcsr = Pcsr
    this%tPcsr_set = .true.

  end subroutine set_Pmat_array

  !--------------------------------------------------------------------------
  subroutine set_Pmat_struct(this, matHS)
    type(ElPhoton) :: this
    type(TMatrixArrayHS) :: matHS(:)

    integer :: iK, ii, ic, ir, nbl, j_col

    if (this%model_type == interaction_models%debugelpt) then
      call create_debug_Pmat(this, matHS(1)%H)
      return
    end if

    this%nk = size(matHS)
    if (.not.allocated(this%Pcsr)) then
      allocate(this%Pcsr(this%nk, 3))
    else
      if (size(this%Pcsr, 1) .ne. this%nk) then
         deallocate(this%Pcsr)
         allocate(this%Pcsr(this%nk, 3))
      end if
    end if
    if (.not.allocated(this%basis%matrixToBasis)) then
      call error_msg('set_Pmat_struct: matrixToBasis not allocated')
    end if

    do ii = 1, 3
      do iK = 1, this%nk
        this%Pcsr(iK, ii) = matHS(iK)%H
      end do
    end do

    do iK = 1, this%nk
      ! Now we compute the polarization matrix P = i (R_i - R_j) * H_ij
      associate(Hcsr => matHS(iK)%H, coord => this%basis%x, &
                  & matrixToBasis => this%basis%matrixToBasis)
      do ir = 1, Hcsr%nrow
        do ic = Hcsr%rowpnt(ir), Hcsr%rowpnt(ir+1)-1
          if (size(this%basis%matrixToBasis) .ne. Hcsr%nrow) then
            call error_msg('set_Pmat_struct: matrixToBasis size does not match matHS%nrow')
          end if
          j_col = Hcsr%colind(ic)
          do ii = 1, 3
            this%Pcsr(iK, ii)%nzval(ic) = (coord(ii, matrixToBasis(ir))-coord(ii, matrixToBasis(j_col))) * &
                             &  Hcsr%nzval(ic) * (0.0_lp, 1.0_lp)
            if (j_col == ir) then
              this%Pcsr(iK, ii)%nzval(ic) = (0.0_lp, 0.0_lp)
            end if
          end do
        end do
      end do
      end associate
    end do
    this%tPcsr_set = .true.

  end subroutine set_Pmat_struct

  !--------------------------------------------------------------------------
  subroutine create_debug_Pmat(this, Hcsr)
    type(ElPhoton) :: this
    type(z_CSR), intent(in) :: Hcsr

    integer :: ir, ic

    this%nk = 1
    if (.not.allocated(this%Pcsr)) then
      allocate(this%Pcsr(this%nk, 3))
    else
      if (size(this%Pcsr, 1) .ne. this%nk) then
         deallocate(this%Pcsr)
         allocate(this%Pcsr(this%nk, 3))
      end if
    end if

    this%Pcsr(1, 1) = Hcsr
    this%Pcsr(1, 1)%nzval = (0.0_lp, 0.0_lp)

    associate(Pcsr => this%Pcsr(1, 1))
    do ir = 1, Pcsr%nrow
       do ic = Pcsr%rowpnt(ir), Pcsr%rowpnt(ir+1)-1
          if (mod(ir, 2) /= 0) then
            if (Pcsr%colind(ic) == ir+1) then
              Pcsr%nzval(ic) = (0.0_lp, 1.0_lp)
            end if
          end if
          if (mod(ir, 2) == 0) then
            if (Pcsr%colind(ic) == ir-1) then
              Pcsr%nzval(ic) = -(0.0_lp, 1.0_lp)
            end if
          end if
       end do
    end do
    end associate
    this%tPcsr_set = .true.

  end subroutine create_debug_Pmat

  !--------------------------------------------------------------------------
  function is_Pmat_set(this) result(is_set)
    class(ElPhoton) :: this
    logical :: is_set
    is_set = this%tPcsr_set
  end function is_Pmat_set

  !--------------------------------------------------------------------------
  subroutine prepare_Kmat(this)
    class(ElPhoton) :: this

    integer :: nbl, k, ii, ir


  end subroutine prepare_Kmat

  !--------------------------------------------------------------------------
  !> Retrieve matrix pointer from cache and add it to ESH
#:def add_sigma_r_template(KIND,MTYPE)
  subroutine add_sigma_r_${KIND}$(this, esh, en_index, k_index, spin)
    class(ElPhoton) :: this
    type(${MTYPE}$), dimension(:,:), intent(inout) :: esh
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    integer :: npl, jj
    type(x_dns), pointer :: tmp_blk
    type(TMatLabel) :: label
    npl = this%struct%num_PLs

    if (this%scba_iter .eq. 0) return

    label%kpoint = k_index
    label%energy_point = en_index
    label%spin = 1
    if (present(spin)) then
      label%spin = spin
    end if

    do jj = 1, npl
      label%row_block = jj
      label%col_block = jj
      call this%sigma_r%retrieve_pointer(tmp_blk, label)
      ESH(jj, jj)%val = ESH(jj, jj)%val - tmp_blk%val
      ! 3-diagonal blocks
      if (this%tTridiagonal) then
        if (jj .lt. npl) then
          label%row_block = jj
          label%col_block = jj + 1
          call this%sigma_r%retrieve_pointer(tmp_blk, label)
          ESH(jj, jj + 1)%val = ESH(jj, jj + 1)%val - tmp_blk%val
          label%row_block = jj + 1
          label%col_block = jj
          call this%sigma_r%retrieve_pointer(tmp_blk, label)
          ESH(jj + 1, jj)%val = ESH(jj + 1, jj)%val - tmp_blk%val
        end if
      end if
    end do

  end subroutine add_sigma_r_${KIND}$
#:enddef add_sigma_r_template

  !--------------------------------------------------------------------------
  !> Retrieve matrix pointer from cache and add it to sigma
#:def add_sigma_n_template(KIND,MTYPE)
  subroutine add_sigma_n_${KIND}$(this, sigma, en_index, k_index, spin)
    class(ElPhoton) :: this
    type(${MTYPE}$), dimension(:,:), intent(inout) :: sigma
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    integer :: npl, jj
    type(x_dns), pointer :: tmp_blk
    type(TMatLabel) :: label
    npl = this%struct%num_PLs

    if (this%scba_iter .eq. 0) return

    label%kpoint = k_index
    label%energy_point = en_index
    label%spin = 1
    if (present(spin)) then
      label%spin = spin
    end if

    do jj = 1, npl
      label%row_block = jj
      label%col_block = jj
      call this%sigma_n%retrieve_pointer(tmp_blk, label)
      sigma(jj, jj)%val = sigma(jj, jj)%val + tmp_blk%val
      ! 3-diagonal blocks
      ! Gn^dag = Gn => Sigma_n(i,j) = Sigma_n(j,i)^dag
      if (this%tTridiagonal) then
        if (jj .lt. npl) then
          label%row_block = jj
          label%col_block = jj + 1
          call this%sigma_n%retrieve_pointer(tmp_blk, label)
          sigma(jj, jj + 1)%val = sigma(jj, jj + 1)%val + tmp_blk%val
          sigma(jj + 1, jj)%val = sigma(jj + 1, jj)%val + conjg(transpose(tmp_blk%val))
        end if
      end if
    end do

  end subroutine add_sigma_n_${KIND}$
#:enddef add_sigma_n_template

#:for PREC in PRECISIONS
     #:set KIND = PREC_ABBREVS[PREC]
     #:set MTYPE = MAT_TYPES['complex'][PREC]

     $:add_sigma_r_template(KIND,MTYPE)

     $:add_sigma_n_template(KIND,MTYPE)
#:endfor

  !--------------------------------------------------------------------------
  !> Returns the lesser (n) Self Energy in block format
  !
  subroutine get_sigma_n_blk(this, blk_sigma_n, en_index, k_index, spin)
    class(ElPhoton) :: this
    type(z_dns), dimension(:,:), intent(inout) :: blk_sigma_n
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    type(TMatLabel) :: label
    integer :: ii

    if (this%scba_iter .eq. 0) return

    label%kpoint = 0
    label%energy_point = 0
    label%spin = 1
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

  subroutine get_Gn_blk(this, blk_Gn, en_index, k_index)
    class(ElPhoton) :: this
    type(z_dns), dimension(:,:), intent(inout) :: blk_Gn
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index

    type(TMatLabel) :: label
    integer :: ii

    if (this%scba_iter .eq. 0) return

    label%kpoint = 0
    label%energy_point = 0
    label%spin = 1
    if (present(k_index)) then
      label%kpoint = k_index
    end if
    if (present(en_index)) then
      label%energy_point = en_index
    end if

    ! Retrieve the diagonal blocks
    do ii = 1, size(blk_Gn, 1)
      label%row_block = ii
      label%col_block = ii
      call this%G_n%retrieve(blk_Gn(ii,ii), label)
    end do

  end subroutine get_Gn_blk

  !--------------------------------------------------------------------------
  !> Returns the lesser (n) Self Energy for a given block
  !
  subroutine get_sigma_n_mat(this, sigma_n, ii, jj, en_index, k_index, spin)
    class(ElPhoton) :: this
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
    label%spin = 1

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
    class(ElPhoton) :: this
    type(z_dns), dimension(:,:), intent(in) :: Gr
    integer, intent(in), optional  :: en_index
    integer, intent(in), optional  :: k_index
    integer, intent(in), optional  :: spin
  end subroutine set_Gr

  !> Dummy subroutine to set Gn
  subroutine set_Gn(this, Gn, en_index, k_index, spin)
    class(ElPhoton) :: this
    type(z_dns), dimension(:,:), intent(in) :: Gn
    integer, intent(in), optional  :: en_index
    integer, intent(in), optional  :: k_index
    integer, intent(in), optional  :: spin
  end subroutine set_Gn

  subroutine init_sigma_r(this, spin)
    class(ElPhoton) :: this
    integer, intent(in), optional :: spin

    integer :: nbl, nEloc, nKloc, iK, iE, ibl
    integer :: PL_start, Np, Mp
    type(TMatLabel) :: label
    type(x_DNS) :: tmp

    nbl = this%struct%num_PLs
    NEloc = this%nE_local
    NKloc = this%nK
    label%spin = 1
    if (present(spin)) then
       label%spin = spin
    end if

    select type(p => this%sigma_r)
    class is(TMatrixCacheMem)
       if (.not.p%isInitialized) then
          if (this%tTridiagonal) then
            call p%init(nEloc, nKloc, nbl, 3, 1)
          else
            call p%init(nEloc, nKloc, nbl, 1, 1)
          end if
       end if
    end select

    ! Init blocks of Sigma_r
    do iK = 1, NKloc
      do iE = 1, NEloc
        label%kpoint = iK
        label%energy_point = iE

        do ibl = 1, nbl
          PL_start = this%struct%mat_PL_start(ibl)
          Np = this%struct%mat_PL_end(ibl) - PL_start + 1
          associate(Sigma_r => this%Sigma_r)

          ! Diag blocks
          label%row_block = ibl
          label%col_block = ibl
          if (.not.Sigma_r%is_cached(label)) then
            call create(tmp,Np,Np)
            tmp%val = (0.0_lp, 0.0_lp)
            call Sigma_r%add(tmp, label)
            call destroy(tmp)
          end if

          if (this%tTridiagonal .and. ibl < nbl) then
            PL_start = this%struct%mat_PL_start(ibl+1)
            Mp = this%struct%mat_PL_end(ibl+1) - PL_start + 1

            ! Subdiag blocks
            label%row_block = ibl
            label%col_block = ibl + 1
            if (.not.Sigma_r%is_cached(label)) then
              call create(tmp,Np,Mp)
              tmp%val = (0.0_lp, 0.0_lp)
              call Sigma_r%add(tmp, label)
              call destroy(tmp)
            end if

            ! Overdiag blocks
            label%row_block = ibl + 1
            label%col_block = ibl
            if (.not.Sigma_r%is_cached(label)) then
              call create(tmp,Mp,Np)
              tmp%val = (0.0_lp, 0.0_lp)
              call Sigma_r%add(tmp, label)
              call destroy(tmp)
            end if

          end if
          end associate
        end do
      end do
    end do

  end subroutine init_sigma_r

  subroutine init_sigma_n(this, spin)
    class(ElPhoton) :: this
    integer, intent(in), optional :: spin

    integer :: nbl, nEloc, nKloc, iK, iE, ibl
    integer :: PL_start, Np, Mp
    type(TMatLabel) :: label
    type(x_DNS) :: tmp

    nbl = this%struct%num_PLs
    NEloc = this%nE_local
    NKloc = this%nK

    label%spin = 1
    if (present(spin)) then
       label%spin = spin
    end if

    select type(p => this%sigma_n)
    class is(TMatrixCacheMem)
       if (.not.p%isInitialized) then
          if (this%tTridiagonal) then
            call p%init(nEloc, nKloc, nbl, 3, 1)
          else
            call p%init(nEloc, nKloc, nbl, 1, 1)
          end if
       end if
    end select

    ! Init blocks of Sigma_n
    do iK = 1, NKloc
      do iE = 1, NEloc
        label%kpoint = iK
        label%energy_point = iE

        do ibl = 1, nbl
          PL_start = this%struct%mat_PL_start(ibl)
          Np = this%struct%mat_PL_end(ibl) - PL_start + 1
          associate(Sigma_n => this%Sigma_n)

          ! Diag blocks
          label%row_block = ibl
          label%col_block = ibl
          if (.not.Sigma_n%is_cached(label)) then
            call create(tmp,Np,Np)
            tmp%val = (0.0_lp, 0.0_lp)
            call Sigma_n%add(tmp, label)
            call destroy(tmp)
          end if

          if (this%tTridiagonal .and. ibl < nbl) then
            PL_start = this%struct%mat_PL_start(ibl+1)
            Mp = this%struct%mat_PL_end(ibl+1) - PL_start + 1

            ! Subdiag blocks
            label%row_block = ibl
            label%col_block = ibl + 1
            if (.not.Sigma_n%is_cached(label)) then
              call create(tmp,Np,Mp)
              tmp%val = (0.0_lp, 0.0_lp)
              call Sigma_n%add(tmp, label)
              call destroy(tmp)
            end if

            ! Overdiag blocks
            label%row_block = ibl + 1
            label%col_block = ibl
            if (.not.Sigma_n%is_cached(label)) then
              call create(tmp,Mp,Np)
              tmp%val = (0.0_lp, 0.0_lp)
              call Sigma_n%add(tmp, label)
              call destroy(tmp)
            end if

          end if
          end associate
        end do
      end do
    end do

  end subroutine init_sigma_n

  !--------------------------------------------------------------------------
  !> Give the Gn at given energy point to the interaction
  !>
  subroutine compute_Sigma_r(this, en_index, k_index, spin)
    class(ElPhoton) :: this
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

#:if defined("MPI")
    integer :: ibl, nbl, NE, NKloc, NEloc, iEshift
    integer :: iKK, PL_start, ir, ic
    complex(lp) :: fac_min, fac_plus
    type(TMatPointer), allocatable :: pGG(:,:), pSigma(:,:), pAA(:,:)
    type(x_DNS), dimension(:,:), allocatable :: Pblk
    type(TMatLabel) :: label

    nbl = this%struct%num_PLs
    NKloc = this%nK
    NE = this%nE_global
    NEloc = this%nE_local
    iEshift = nint(this%wq/this%dE)
    label%spin = 1
    if (present(spin)) then
       label%spin = spin
    end if

    allocate(pGG(NEloc,NKloc))
    allocate(pAA(NEloc,NKloc))
    allocate(pSigma(NEloc,NKloc))
    allocate(Pblk(nbl,nbl))

    ! k-loop for the self energy
    do iKK = 1, this%nk

      call csr2blk3diag(this%Pcsr(iKK,1), Pblk, this%struct%mat_PL_start)

      do ibl = 1, nbl
        !print*, "Debug Pblk", ibl
        !print*, Pblk(ibl,ibl)%val
        call compute_diagonal_contributions()

        if (this%tTridiagonal .and. ibl < nbl) then
          call compute_supradiagonal_contributions()
          call compute_subdiagonal_contributions()
        end if

      end do
      call destroy_tridiag_blk(Pblk)
    end do

    deallocate(pAA)
    deallocate(pGG)
    deallocate(pSigma)
    deallocate(Pblk)
#:else
    call error_msg("inelastic scattering requires MPI compilation")
#:endif
    contains

    !> Generic subroutine performs the product between blocks:
    !         P(iM,jM) * G(i,j) * P(iMdag, jMdag)
    ! The correct pointers for Gr(iE,iK,i,j) are set so the final calculation is
    ! the full Sigma^r expression:
    !
    ! Sigma^r = -i/2 (Sigma^n + Sigma^p)
    ! Sigma^r = -i/2 P ( Nq A(+) + (Nq+1) A(-) - G^n(-) + G^n(+) ) P*

    subroutine perform_block_product(iM, jM, iG, jG, iMdag, jMdag)
      integer, intent(in) :: iM, jM, iG, jG, iMdag, jMdag

      if (maxval(abs(Pblk(iM,jM)%val)) > 0.0_lp .and. &
        & maxval(abs(Pblk(IMdag,jMdag)%val)) > 0.0_lp) then

        call setup_pointers_A(iG, jG)

        fac_min = cmplx(0.0, -0.5*(this%Nq+1.0), lp)
        fac_plus = cmplx(0.0, -0.5*this%Nq, lp)

        call selfenergy(this%cart_comm, pAA, fac_min, fac_plus, iEshift, Pblk(iM,jM), &
            & Pblk(iMdag,jMdag), iKK, NE, NEloc, pSigma, this%prefactor)

        call release_pointers_A()

        call setup_pointers_Gn(iG, jG)

        fac_min = (0.0_lp, 0.5_dp)
        fac_plus = (0.0_lp, -0.5_dp)

        call selfenergy(this%cart_comm, pGG, fac_min, fac_plus, iEshift, Pblk(iM,jM), &
            & Pblk(iMdag,jMdag), iKK, NE, NEloc, pSigma, this%prefactor)

      end if

    end subroutine perform_block_product

    subroutine compute_diagonal_contributions()

      call setup_pointers_Sigma_r(ibl,ibl)

      ! M(i,i)G(i,i)M(i,i)
      call perform_block_product(ibl,ibl, ibl,ibl, ibl,ibl)

      if (ibl .ne. nbl) then
        ! M(i,i+1) G(i+1,i+1) M(i+1,i)
        call perform_block_product(ibl,ibl+1, ibl+1,ibl+1, ibl+1,ibl)
        ! M(i,i) G(i,i+1) M(i+1,i)
        call perform_block_product(ibl,ibl, ibl,ibl+1, ibl+1,ibl)
        ! M(i,i+1) G(i+1,i) M(i,i)
        call perform_block_product(ibl,ibl+1, ibl+1,ibl, ibl,ibl)
      end if

      if (ibl .ne. 1) then
        ! M(i,i-1) G(i-1,i-1) M(i-1,i)
        call perform_block_product(ibl,ibl-1, ibl-1,ibl-1, ibl-1,ibl)
        ! M(i,i) G(i,i-1) M(i-1,i)
        call perform_block_product(ibl,ibl, ibl,ibl-1, ibl-1,ibl)
        ! M(i,i-1) G(i-1,i) M(i,i)
        call perform_block_product(ibl,ibl-1, ibl-1,ibl, ibl,ibl)
      end if

    end subroutine compute_diagonal_contributions

    !> blocks for Sigma(ibl, ibl+1)
    subroutine compute_supradiagonal_contributions()

      call setup_pointers_Sigma_r(ibl, ibl+1)

      ! M(i,i) G(i,i) M(i,i+1)
      call perform_block_product(ibl,ibl, ibl,ibl, ibl,ibl+1)
      ! M(i,i+1) G(i+1,i+1) M(i+1,i+1)
      call perform_block_product(ibl,ibl+1, ibl+1,ibl+1, ibl+1,ibl+1)
      ! M(i,i) G(i,i+1) M(i+1,i+1)
      call perform_block_product(ibl,ibl, ibl,ibl+1, ibl+1,ibl+1)
      ! M(i,i+1) G(i+1,i) M(i,i+1)
      call perform_block_product(ibl,ibl+1, ibl+1,ibl, ibl,ibl+1)

      if (ibl > 1) then
        ! M(i,i-1) G(i-1,i) M(i,i+1)
        call perform_block_product(ibl,ibl-1, ibl-1,ibl, ibl,ibl+1)
      end if

      if (ibl < nbl-1) then
        ! M(i,i+1) G(i+1,i+2) M(i+2,i+1)
        call perform_block_product(ibl,ibl+1, ibl+1,ibl+2, ibl+2,ibl+1)
      end if


    end subroutine compute_supradiagonal_contributions

    !> blocks for Sigma(ibl+1, ibl)
    subroutine compute_subdiagonal_contributions()

      call setup_pointers_Sigma_r(ibl+1,ibl)

      ! M(i+1,i)G(i,i)M(i,i)
      call perform_block_product(ibl+1,ibl, ibl,ibl, ibl,ibl)
      ! M(i+1,i+1)G(i+1,i+1)M(i+1,i)
      call perform_block_product(ibl+1,ibl+1, ibl+1,ibl+1, ibl+1,ibl)
      ! M(i+1,i+1)G(i+1,i)M(i,i)
      call perform_block_product(ibl+1,ibl+1, ibl+1,ibl, ibl,ibl)
      ! M(i+1,i)G(i,i+1)M(i+1,i)
      call perform_block_product(ibl+1,ibl, ibl,ibl+1, ibl+1,ibl)

      if (ibl < nbl-1) then
        ! M(i+1,i+2)G(i+2,i+1)M(i+1,i)
        call perform_block_product(ibl+1,ibl+2, ibl+2,ibl+1, ibl+1,ibl)
      end if

      if (ibl > 1) then
        ! M(i+1,i)G(i,i-1)M(i-1,i)
        call perform_block_product(ibl+1,ibl, ibl,ibl-1, ibl-1,ibl)
      end if

    end subroutine compute_subdiagonal_contributions


    ! setup the array of pointers to G_r and sigma_r
    subroutine setup_pointers_Gr(Nr, Nc)
      integer, intent(in) :: Nr, Nc

      integer :: iK, iE

      label%row_block = Nr
      label%col_block = Nc

      do iK = 1, NKloc
        do iE = 1, NEloc
          label%kpoint = iK
          label%energy_point = iE
          call this%G_r%retrieve_pointer(pGG(iE,iK)%pMat, label)
        end do
      end do
    end subroutine setup_pointers_Gr

    ! setup the array of pointers to G_n
    ! Only the current k-point iKK is set, as selfenergy uses only that column
    subroutine setup_pointers_Gn(Nr, Nc)
      integer, intent(in) :: Nr, Nc

      integer :: iE

      label%row_block = Nr
      label%col_block = Nc
      label%kpoint = iKK

      do iE = 1, NEloc
        label%energy_point = iE
        call this%G_n%retrieve_pointer(pGG(iE,iKK)%pMat, label)
      end do
    end subroutine setup_pointers_Gn

    ! Computes the spectral function blocks A(Nr,Nc) = i [G_r(Nr,Nc) - G_r(Nc,Nr)^dag]
    ! Note: zspectral cannot be used here since it computes i [G - G^dag] from a
    !       single block, which is correct only for diagonal blocks. For off-diagonal
    !       blocks the advanced part must be taken from the transposed block (Nc,Nr).
    ! Only the current k-point iKK is computed, as selfenergy uses only that column
    subroutine setup_pointers_A(Nr, Nc)
      integer, intent(in) :: Nr, Nc

      integer :: iE
      type(x_DNS), pointer :: pGji

      label%kpoint = iKK

      do iE = 1, NEloc
        label%energy_point = iE
        label%row_block = Nc
        label%col_block = Nr
        call this%G_r%retrieve_pointer(pGji, label)
        label%row_block = Nr
        label%col_block = Nc
        call this%G_r%retrieve_pointer(pGG(iE,iKK)%pMat, label)
        allocate(pAA(iE,iKK)%pMat)
        call create(pAA(iE,iKK)%pMat, pGG(iE,iKK)%pMat%nrow, pGG(iE,iKK)%pMat%ncol)
        pAA(iE,iKK)%pMat%val = (0.0_lp, 1.0_lp) * &
            & (pGG(iE,iKK)%pMat%val - conjg(transpose(pGji%val)))
      end do
    end subroutine setup_pointers_A

    subroutine release_pointers_A()
      integer :: iE

      do iE = 1, NEloc
        call destroy(pAA(iE,iKK)%pMat)
        deallocate(pAA(iE,iKK)%pMat)
      end do
    end subroutine release_pointers_A

    ! setup the array of pointers to Sigma_r
    ! Only the current k-point iKK is set, as selfenergy uses only that column
    subroutine setup_pointers_Sigma_r(Nr,Nc)
      integer, intent(in) :: Nr, Nc

      integer :: iE

      label%row_block = Nr
      label%col_block = Nc
      label%kpoint = iKK

      associate(Sigma_r => this%Sigma_r)
      do iE = 1, NEloc
        label%energy_point = iE
        call Sigma_r%retrieve_pointer(pSigma(iE,iKK)%pMat, label)
      end do
      end associate

    end subroutine setup_pointers_Sigma_r

    subroutine check_elements(Mat,arg)
      class(TMatrixCache), intent(in) :: Mat
      character(*), intent(in) :: arg

      type(x_DNS), pointer :: pMat
      integer :: iK, iE

      do iK = 1, NKloc
        do iE = 1, NEloc
          label%kpoint = iK
          label%energy_point = iE
          call Mat%retrieve_pointer(pMat, label)
          if (any(isNaN(abs(pMat%val)))) then
             print*,arg//'=NaN',iE,iK
          end if
        end do
      end do
    end subroutine check_elements

  end subroutine compute_Sigma_r


  !> Give the Gn at given energy point to the interaction
  subroutine compute_Sigma_n(this, en_index, k_index, spin)
    class(ElPhoton) :: this
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

#:if defined("MPI")
    integer :: ibl, nbl, NE, NKloc, NEloc, iEshift
    integer :: iKK, PL_start
    complex(lp) :: fac_min, fac_plus
    type(TMatPointer), allocatable :: pGG(:,:), pSigma(:,:)
    type(x_DNS), dimension(:,:), allocatable :: Pblk
    type(TMatLabel) :: label

    nbl = this%struct%num_PLs
    NKloc = this%nK
    NE = this%nE_global
    NEloc = this%nE_local
    iEshift = nint(this%wq/this%dE)
    label%spin = 1
    if (present(spin)) then
      label%spin = spin
    end if

    allocate(pGG(NEloc,NKloc))
    allocate(pSigma(NEloc,NKloc))
    allocate(Pblk(nbl,nbl))

    ! k-loop for the self energy
    do iKK = 1, this%nk

      call csr2blk3diag(this%Pcsr(iKK,1), Pblk, this%struct%mat_PL_start)

      do ibl = 1, nbl
        call compute_diagonal_contributions()

        if (this%tTridiagonal .and. ibl < nbl) then
          call compute_supradiagonal_contributions()
        end if

      end do
      call destroy_tridiag_blk(Pblk)
    end do

    deallocate(pGG)
    deallocate(pSigma)
    deallocate(Pblk)
#:else
    error stop "inelastic scattering requires MPI compilation"
#:endif

    contains

    !> Generic subroutine performs the product between blocks:
    !         P(iM,jM) * G(i,j) * P(iMdag, jMdag)
    ! The correct pointers for Gr(iE,iK,i,j) are set so the final calculation is
    ! the full Sigma_n expression
    subroutine perform_block_product(iM, jM, iG, jG, iMdag, jMdag)
      integer, intent(in) :: iM, jM, iG, jG, iMdag, jMdag

      if (maxval(abs(Pblk(iM,jM)%val)) > 0.0_lp .and. &
        & maxval(abs(Pblk(IMdag,jMdag)%val)) > 0.0_lp) then

        call setup_pointers_Gn(iG, jG)

        fac_min = cmplx(this%Nq, 0.0, lp)
        fac_plus = cmplx(this%Nq + 1, 0.0, lp)

        call selfenergy(this%cart_comm, pGG, fac_min, fac_plus, iEshift, Pblk(iM,jM), &
            & Pblk(iMdag,jMdag), iKK, NE, NEloc, pSigma, this%prefactor)
      end if

    end subroutine perform_block_product

    !> blocks for Sigma(ibl, ibl)
    subroutine compute_diagonal_contributions()

      call setup_pointers_Sigma_n(ibl, ibl)

      ! M(i,i)G(i,i)M(i,i)
      call perform_block_product(ibl,ibl, ibl,ibl, ibl,ibl)

      if (ibl .ne. nbl) then
        ! M(i,i+1)G(i+1,i+1)M(i+1,i)
        call perform_block_product(ibl,ibl+1, ibl+1,ibl+1, ibl+1,ibl)
        ! M(i,i)G(i,i+1)M(i+1,i)
        call perform_block_product(ibl,ibl, ibl,ibl+1, ibl+1,ibl)
        ! M(i,i+1)G(i+1,i)M(i,i)
        call perform_block_product(ibl,ibl+1, ibl+1,ibl, ibl,ibl)
      end if

      if (ibl .ne. 1) then
        ! M(i,i-1)G(i-1,i-1)M(i-1,i)
        call perform_block_product(ibl,ibl-1, ibl-1,ibl-1, ibl-1,ibl)
        ! M(i,i)G(i,i-1)M(i-1,i)
        call perform_block_product(ibl,ibl, ibl,ibl-1, ibl-1,ibl)
        ! M(i,i-1)G(i-1,i)M(i,i)
        call perform_block_product(ibl,ibl-1, ibl-1,ibl, ibl,ibl)
      end if

    end subroutine compute_diagonal_contributions

    !> blocks for Sigma(ibl, ibl+1)
    subroutine compute_supradiagonal_contributions()

      call setup_pointers_Sigma_n(ibl, ibl+1)

      ! M(i,i)G(i,i)M(i,i+1)
      call perform_block_product(ibl,ibl, ibl,ibl, ibl,ibl+1)
      ! M(i,i+1)G(i+1,i+1)M(i+1,i+1)
      call perform_block_product(ibl,ibl+1, ibl+1,ibl+1, ibl+1,ibl+1)
      ! M(i,i)G(i,i+1)M(i+1,i+1)
      call perform_block_product(ibl,ibl, ibl,ibl+1, ibl+1,ibl+1)
      ! M(i,i+1)G(i+1,i)M(i,i+1)
      call perform_block_product(ibl,ibl+1, ibl+1,ibl, ibl,ibl+1)

      if (ibl > 1) then
        ! M(i,i-1)G(i-1,i)M(i,i+1)
        call perform_block_product(ibl,ibl-1, ibl-1,ibl, ibl,ibl+1)
      end if

      if (ibl < nbl-1) then
        ! M(i,i+1)G(i+1,i+2)M(i+2,i+1)
        call perform_block_product(ibl,ibl+1, ibl+1,ibl+2, ibl+2,ibl+1)
      end if

    end subroutine compute_supradiagonal_contributions

    ! setup the array of pointers to G_n
    ! Only the current k-point iKK is set, as selfenergy uses only that column
    subroutine setup_pointers_Gn(Nr, Nc)
      integer, intent(in) :: Nr, Nc

      integer :: iE

      label%row_block = Nr
      label%col_block = Nc
      label%kpoint = iKK

      do iE = 1, NEloc
        label%energy_point = iE
        call this%G_n%retrieve_pointer(pGG(iE,iKK)%pMat, label)
      end do
    end subroutine setup_pointers_Gn

    ! setup the array of pointers to Sigma_n
    ! Only the current k-point iKK is set, as selfenergy uses only that column
    subroutine setup_pointers_Sigma_n(Nr,Nc)
      integer, intent(in) :: Nr, Nc

      integer :: iE

      label%row_block = Nr
      label%col_block = Nc
      label%kpoint = iKK

      associate(Sigma_n => this%Sigma_n)
      do iE = 1, NEloc
        label%energy_point = iE
        call Sigma_n%retrieve_pointer(pSigma(iE,iKK)%pMat, label)
      end do
      end associate
    end subroutine setup_pointers_Sigma_n

    subroutine check_elements(Mat,arg)
      class(TMatrixCache), intent(in) :: Mat
      character(*), intent(in) :: arg

      integer :: iK, iE
      type(x_DNS), pointer :: pMat

      do iK = 1, NKloc
        do iE = 1, NEloc
          label%kpoint = iK
          label%energy_point = iE
          call Mat%retrieve_pointer(pMat, label)
          if (any(isNaN(abs(pMat%val)))) then
             print*,arg//'=NaN',iE,iK
          end if
        end do
      end do
    end subroutine check_elements

  end subroutine compute_Sigma_n

  subroutine selfenergy(comm2d, GG, fac_minus, fac_plus, iEhbaromega, PP, PP_dag, &
              & iK, NE, NEloc, Sigma, prefactor)
    type(MPI_Comm), intent(in) :: comm2d
    type(TMatPointer), intent(in) :: GG(0:,1:)
    complex(lp), intent(in) :: fac_minus
    complex(lp), intent(in) :: fac_plus
    type(x_DNS), intent(in) :: PP, PP_dag
    integer, intent(in) :: iEhbaromega
    integer, intent(in) :: iK
    integer, intent(in) :: NE
    integer, intent(in) :: NEloc
    type(TMatPointer), intent(in) :: Sigma(0:,1:)
    real(lp), intent(in) :: prefactor

    ! locals
    integer :: Np, Mp, nu
    integer :: iEminus, iEplus, iE, iEglo
    integer :: myid, commsize, ndims = 2
    integer :: dims(2), coords(2)
    type(x_DNS), target :: rbuff1, rbuff2, sbuffH, rbuffH, support_Sigma
    type(x_DNS), pointer :: pbuff1, pbuff2, pGG1, pGG2, pSigma
    type(x_DNS) :: work1, work2, Kdag
    integer :: ierr
    type(MPI_Request) :: rqE1,rqE2,rqE3,rqE4
    type(MPI_Status) :: statE(4)

    integer :: msource, mdest, psource, pdest
    integer :: ndiff

    nullify(pbuff1)
    nullify(pbuff2)

    dims(1) = 0 ! We don't need the kGrid dimension, as there is no k-communication in this case
    dims(2) = NE/NEloc

    call MPI_comm_size(comm2d, commsize, ierr)
    call MPI_comm_rank(comm2d, myid, ierr)
    call MPI_cart_coords(comm2d, myid, ndims, coords, ierr)
    call MPI_barrier(comm2d, ierr)

    Np = size(GG(0,iK)%pMat%val,1)
    Mp = size(GG(0,iK)%pMat%val,2)

    call create(rbuff1, Np, Mp)
    call create(rbuff2, Np, Mp)
    call create(sbuffH, Np, Mp)

    ! Sigma_ij(iK, iE) = KK_il(iK) * (fac_minus * GG_lm(iK, E-wq) + fac_plus * GG_lm(iK, E+wq))*K_mj

      eloop: do iE = 0, NEloc-1
        iEglo = iE + coords(2)*NEloc
        !////////////////////////////////////////////////////////////////////////////////////
        !// Communications of G(k, E-hwq)
        !////////////////////////////////////////////////////////////////////////////////////
        !//  012345 012345 012345 012345    iE < NEloc = 6
        !//  -------------------------------------------------------------------------------
        !//  000000 000011 111111 112222
        !//  012345 678901 234567 890123    iEglo = 1 + 6 = 7
        !// |oooooo|oooooo|oooooo|oooooo|
        !//  ^-----|-E----|--^              ihbarOmega = 7; ndiff = 1 iMinus = 0
        !//     ^--|----E-|-----^                           ndiff = 1 iMinus = 3
        !//       ^|------|E-----|--^                       ndiff = 2 iMinus = 5
        !//^|------|E-----|-^                               ndiff = 1 iMinus = 0
        !////////////////////////////////////////////////////////////////////////////////////
        nullify(pbuff1)
        nullify(pbuff2)
        !// pbuff1 => G(k,E-wq)
        !// checks if iE-iEhbaromega is on the same processor => no communication
        !//          for E<0 take G(E<0) = 0
        mdest = MPI_PROC_NULL
        msource = MPI_PROC_NULL
        iEminus = iE - iEhbaromega
        if (iEminus >= 0) then
          pbuff1 => GG(iEminus, iK)%pMat
        else !(iEminus < 0)
          rbuff1%val=(0.0_lp, 0.0_lp)
          pbuff1 => rbuff1
        end if

        ! MPI communications
        if (dims(2) > 1) then
          if (iE < iEhbaromega) then
            if (iEglo - iEhbaromega < 0 ) then
              rbuff1%val=(0.0_lp, 0.0_lp)
              pbuff1 => rbuff1
              ndiff = (iEhbaromega-iEglo)/NEloc + iEglo/NEloc
              if (mod(iEhbaromega-iEglo, NEloc) /= 0) then
                ndiff = ndiff + 1
              end if
            else
              ! get processor distance between iEglo and iEglo-ihbarOmega
              ndiff = iEglo/NEloc - (iEglo-iEhbaromega)/NEloc
            end if

            iEminus = mod(iE + (ndiff+1)*NEloc - iEhbaromega, NEloc)

            call MPI_Cart_shift( comm2d, 1, ndiff, msource, mdest, ierr)

            if (mdest /= MPI_PROC_NULL) then
               pGG1 => GG(iEminus, iK)%pMat
               call MPI_Isend(pGG1%val, Mp*Np, MPI_XX_COMPLEX, mdest, 41, comm2d, rqE1, ierr)
            end if

            if (msource /= MPI_PROC_NULL) then
               call MPI_Irecv(rbuff1%val, Mp*Np, MPI_XX_COMPLEX, msource, 41, comm2d, rqE2, ierr)
               pbuff1 => rbuff1
            end if

            ! Wait for complete communications
            if (mdest/=MPI_PROC_NULL) call MPI_Wait(rqE1, statE(1), ierr);
            if (msource/=MPI_PROC_NULL) call MPI_Wait(rqE2, statE(2), ierr);

          end if
        end if

        !////////////////////////////////////////////////////////////////////////////////////
        !// Communications of G(k, E+hwq)
        !////////////////////////////////////////////////////////////////////////////////////
        !// checks if iE+iEhbaromega is on the same processor => no communication
        !// pbuff2 points to G(k,E+wq)
        pdest = MPI_PROC_NULL
        psource = MPI_PROC_NULL
        iEplus = iE + iEhbaromega
        if (iEplus < NEloc) then
          pbuff2 => GG(iEplus, iK)%pMat
        else
          rbuff2%val=(0.0_lp, 0.0_lp)
          pbuff2 => rbuff2
        end if

        if (dims(2) > 1) then
          if (iE + iEhbaromega >= NEloc) then
            if (iEglo + iEhbaromega >= NE ) then
              rbuff2%val=(0.0_lp, 0.0_lp)
              pbuff2 => rbuff2
            end if

            ! get processor distance between iEglo+ihbarOmega and iEglo
            ndiff = (iEglo+iEhbaromega)/NEloc - iEglo/NEloc
            iEplus = mod(iEglo + iEhbaromega, NEloc)

            call MPI_Cart_shift( comm2d, 1, -ndiff, psource, pdest, ierr)

            if (pdest /= MPI_PROC_NULL) then
               pGG2 => GG(iEplus, iK)%pMat
               call MPI_Isend(pGG2%val, Mp*Np, MPI_XX_COMPLEX, pdest, 42, comm2d, rqE3, ierr)
            end if

            if (psource /= MPI_PROC_NULL) then
               call MPI_Irecv(rbuff2%val, Mp*Np, MPI_XX_COMPLEX, psource, 42, comm2d, rqE4, ierr)
               pbuff2 => rbuff2
            end if

            if(pdest/=MPI_PROC_NULL) call MPI_Wait(rqE3, statE(3), ierr);
            if(psource/=MPI_PROC_NULL) call MPI_Wait(rqE4, statE(4), ierr);

          end if
        end if

        !//   Compute:
        !//   sbuffH_ij(iK) = (fac_minus * GG_ij(iK, E-wq) + fac_plus * GG_ij(iK, E+wq))
        !$OMP PARALLEL DO private (nu)
        do nu = 1, Mp
          sbuffH%val(:,nu) = (fac_minus*pbuff1%val(:,nu) + fac_plus*pbuff2%val(:,nu))
        end do
        !$OMP END PARALLEL DO

        pSigma => Sigma(iE, iK)%pMat
        call prealloc_mult(PP, sbuffH, work1)
        call prealloc_mult(work1, PP_dag, support_Sigma)
        pSigma%val = pSigma%val + prefactor * support_Sigma%val
        call destroy(support_Sigma)
        call destroy(work1)

      end do eloop

    call destroy(rbuff1)
    call destroy(rbuff2)
    call destroy(sbuffH)

  end subroutine selfenergy

  !--------------------------------------------------------------------------
  !> Destroy Sigma_r
  subroutine destroy_Sigma_r(this)
    class(ElPhoton) :: this
    if (associated(this%sigma_r)) then
      call this%sigma_r%destroy()
    end if
  end subroutine destroy_Sigma_r

  !> Destroy Sigma_n
  subroutine destroy_Sigma_n(this)
    class(ElPhoton) :: this
    if (associated(this%sigma_n)) then
      call this%sigma_n%destroy()
    end if
  end subroutine destroy_Sigma_n

    !> Zero-out Sigma_r
  subroutine clear_Sigma_r(this)
    class(ElPhoton) :: this
    if (associated(this%sigma_r)) then
      call this%sigma_r%clear()
    end if
  end subroutine clear_Sigma_r

  !> Zero-out Sigma_n
  subroutine clear_Sigma_n(this)
    class(ElPhoton) :: this
    if (associated(this%sigma_n)) then
      call this%sigma_n%clear()
    end if
  end subroutine clear_Sigma_n

  subroutine destroy_elphot(this)
    class(ElPhoton) :: this

    integer :: iK, ix

    if (allocated(this%Pcsr)) then
      do iK = 1, size(this%Pcsr,1)
        do ix = 1, size(this%Pcsr,2)
           if (allocated(this%Pcsr(iK,ix)%nzval)) then
              call destroy(this%Pcsr(iK,ix))
           end if
        end do
      end do
      deallocate(this%Pcsr)
    end if
    call this%destroy_sigma_r()
    call this%destroy_sigma_n()

  end subroutine destroy_elphot


  subroutine dump_mat(matrix, iter)
    type(z_DNS), intent(in) :: matrix
    integer :: iter


    character(len=200) :: filename
    type(z_DNS), dimension(:,:), allocatable :: blk
    real(lp) :: real_part, imag_part

    complex(lp), dimension(:), allocatable :: diag

    ! Local variables
    integer :: i, j, iE, ibl, nbl, PL_start, Np
    integer :: unit_number_imag, unit_number_real, rows, cols

    write(filename, '("Pmatrix")')
    unit_number_real = 10 ! Assign a file unit number (can be any number not in use)
    unit_number_imag = 20 ! Assign a file unit number (can be any number not in use)

    rows = matrix%nrow
    cols = matrix%ncol

    open(unit=unit_number_imag, file=trim(filename)//"_imag.dat", status='replace')
    open(unit=unit_number_real, file=trim(filename)//"_real.dat", status='replace')
    ! Write the matrix to the file
    do i = 1, rows
        do j = 1, cols

            real_part = real(matrix%val(i, j)) ! Extract real part
            imag_part = aimag(matrix%val(i, j)) ! Extract imaginary part

            write(unit_number_imag, '(E20.8)', advance='no') imag_part
            write(unit_number_real, '(E20.8)', advance='no') real_part
            if (j < cols) then
                write(unit_number_imag, '(A)', advance='no') ', ' ! Space between columns
                write(unit_number_real, '(A)', advance='no') ', '
            endif
        end do
        write(unit_number_imag, *) ! New line after each row
        write(unit_number_real, *)
    end do
    ! Close the file
    close(unit_number_real)
    close(unit_number_imag)

  end subroutine dump_mat

end module elptinel
