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

! Inelastic electron-phonon interactions
! This is the base class for inelastic interactions
#:include "types.fypp"

module elphinel
  use ieee_arithmetic, only : isnan => ieee_is_nan
  use ln_precision, only : sp, dp, lp => dp
  use mat_def, only : z_csr, c_dns, z_dns, x_dns => z_dns, create, destroy, cross_product, volume
  use ln_constants, only : pi
  use interactions, only : TInteraction
  use ln_inelastic, only : TInelastic
  use ln_allocation, only : log_allocate, log_deallocate, writeMemInfo
  use ln_structure, only : TStruct_info
  use ln_structure, only : TBasisCenters, create_TBasis, destroy_TBasis
  use distributions, only : bose
  use ln_cache
  use iso_c_binding
  use mpi_f08
  use mpi_globals
  use inversions, only : inverse
  use self_energy, only: TMatPointer, selfenergy
  use equiv_kpoints, only : TEqPoint, TEqPointsArray, create, destroy
  use clock
  use ln_messages, only : error_msg
  implicit none
  private

  public :: ElPhonInel, ElPhonPolarOptical, ElPhonNonPolarOptical
  public :: ElPhonPO_create, ElPhonNonPO_create
  public :: ElPhonPO_init, ElPhonNonPO_init

  type, abstract, extends(TInelastic) :: ElPhonInel
    private
    !> communicator of the cartesian grid
    type(MPI_Comm) :: cart_comm

    !> holds atomic structure: Only the scattering region
    type(TBasisCenters) :: basis

    !> Bose-Einstein phonon occupation
    real(dp) :: Nq

    !> A general scalar coupling
    real(dp) :: coupling

    !> Kpoints
    real(dp) :: recVecs2p(3,3)
    real(dp), allocatable :: kpoint(:,:)
    real(dp), allocatable :: kweight(:)
    integer, allocatable :: local_kindex(:)
    !> Equivalent Kpoints. For extending the irreducible wedge
    type(TEqPointsArray), allocatable :: equivalent_kpoints
    !> Map of (iK, cpuID) -> iKglobal
    integer, allocatable :: kindices(:,:)
    !> Energy grid. indeces of local points
    integer, allocatable :: local_Eindex(:)

    !> Energy grid. global num of energy points
    integer :: nE_global
    !> Energy grid. Local num of energy points
    integer :: nE_local
    !> Energy grid spacing
    real(dp) :: dE
    !> Matrix KK
    real(lp), allocatable :: Kmat(:,:,:)
    !> binning resolution for z-coordinates
    real(dp) :: dz
    !> Supercell area for integrations
    real(dp) :: cell_area

  contains

    procedure :: add_sigma_r_sp
    procedure :: add_sigma_r_dp
    procedure :: add_sigma_n_sp
    procedure :: add_sigma_n_dp
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
    type(MPI_Comm), intent(in) :: comm
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

    call set_reciprocal_vectors(basis, this%recVecs2p)

    ! Initialize the cache space
    if (associated(this%sigma_r)) then
      call error_msg('Initialization of PO Phonon called twice') 
    end if
    allocate(TMatrixCacheMem::this%sigma_r)
    select type(p => this%sigma_r)
    type is(TMatrixCacheMem)
      p%tagname = 'Sigma_r'
    end select   

    if (associated(this%sigma_n)) then
      call error_msg('Initialization of PO Phonon called twice') 
    end if
    allocate(TMatrixCacheMem::this%sigma_n)
    select type(p => this%sigma_n)
    type is(TMatrixCacheMem)
      p%tagname = 'Sigma_n'
    end select   

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
    type(MPI_Comm), intent(in) :: comm
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
    if (associated(this%sigma_r)) then
      call error_msg('Initialization of non-PO Phonon called twice') 
    end if
    allocate(this%sigma_r, source=TMatrixCacheMem(tagname='Sigma_r'))
    
    if (associated(this%sigma_n)) then
      call error_msg('Initialization of non-PO Phonon called twice') 
    end if
    allocate(this%sigma_n, source=TMatrixCacheMem(tagname='Sigma_n'))

  end subroutine ElPhonNonPO_init
  !--------------------------------------------------------------------------
  ! This function should be called before the SCBA loop
  subroutine set_kpoints(this, kpoints, kweights, kindex, equiv_kpoints)
    class(ElPhonInel) :: this
    real(dp), intent(in) :: kpoints(:,:)
    real(dp), intent(in) :: kweights(:)
    integer, intent(in) :: kindex(:)
    type(TEqPointsArray), intent(in), optional :: equiv_kpoints

    type(MPI_Comm) :: kComm

    integer :: nKloc, nKprocs, mpierr, myid
    logical :: remain_dims(2) = [.true., .false.]

    this%kpoint = kpoints
    this%kweight = kweights
    this%local_kindex = kindex
    nKloc = size(kindex)
    nKprocs = size(kweights)/nKloc

    if (present(equiv_kpoints)) then
      call create(this%equivalent_kpoints, equiv_kpoints)
    else
      if (allocated(this%equivalent_kpoints)) then
         call destroy(this%equivalent_kpoints)
      end if
    endif

    if (allocated(this%kindices)) call log_deallocate(this%kindices)
    call log_allocate(this%kindices, nKloc, Nkprocs)

    call MPI_cart_sub(this%cart_comm, remain_dims, kComm, mpierr)
    call MPI_allgather(this%local_kindex,nKloc,MPI_INTEGER, &
                     & this%kindices,nKloc,MPI_INT,kComm,mpierr)
    call MPI_comm_rank(kComm, myid, mpierr)

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

    integer :: iZ, iQ, eQ, iK, nDeltaZ, nCentralAtoms, n_eq
    real(dp) :: kq(3), kk(3), ekp(3), QQ(3), Q2, bb, z_mn, Kf
    real(dp) :: zmin, zmax
    real(dp), allocatable :: kpoint(:,:), ekpoints(:,:)

    integer :: transDir
    type(TEqPoint) :: eqv_points_iQ

    ! Compute the matrix Kmat as lookup table
    nCentralAtoms = this%basis%nCentralAtoms

    transDir = this%basis%transportDirection
    zmin = minval(this%basis%x(transDir,1:nCentralAtoms))
    zmax = maxval(this%basis%x(transDir,1:nCentralAtoms))

    nDeltaZ = nint((zmax - zmin)/this%dz)
    if (allocated(this%Kmat)) then
       call log_deallocate(this%Kmat)
    end if
    call log_allocate(this%Kmat, nDeltaZ+1, size(this%kweight), size(this%kweight))

    this%Ce = this%coupling*this%wq/2.0_dp/this%cell_area* &
          & (1.0_dp/this%eps_inf - 1.0_dp/this%eps0)

    !compute the absolute k-points
    allocate(kpoint(3,size(this%kweight)))
    if (all(this%basis%lattVecs == 0.0_dp)) then
       kpoint = 0.0_dp
    else
       kpoint = matmul(this%recVecs2p, this%kpoint)
    end if

    iQloop: do iQ = 1, size(this%kweight)
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

      ! Add contribution of kpoints equivalent to iQ
      if (allocated(this%equivalent_kpoints)) then
        !Select the set of points equivalent to iQ
        eqv_points_iQ = this%equivalent_kpoints%EqPoints(iQ)
        n_eq = size(eqv_points_iQ%points, 2)
        !Compute the absolute k-points
        call log_allocate(ekpoints, 3, n_eq)
        ekpoints = matmul(this%recVecs2p, eqv_points_iQ%points)

        do eQ = 1, n_eq
          ekp = ekpoints(:,eQ)
          do iK = 1, size(this%kweight)
            kk = kpoint(:, iK)
            QQ = kk - ekp
            Q2 = dot_product(QQ, QQ)
            bb = sqrt(this%q0*this%q0 + Q2)

            do iZ = 0, nDeltaZ
              z_mn = iZ * this%dz
              Kf = (2.0_dp*Q2 + this%q0*this%q0*(1.0_dp-bb*z_mn))*exp(-bb*z_mn)/ (4.0_dp*bb**3)
              this%Kmat(iZ+1,iK,iQ) = this%Kmat(iZ+1,iK,iQ) + this%Ce * this%kweight(iQ) * Kf
            end do
          end do
        enddo
        call log_deallocate(ekpoints)
      endif

    end do iQloop

    deallocate(kpoint)
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

    integer :: iZ, iQ, iK, nDeltaZ, nCentralAtoms, eQ, n_eq
    real(dp) :: kq(3), kk(3), QQ, Q2, kq2, kk2, z_mn, Kf, ekp(3), ekp2
    real(dp) :: zmin, zmax
    real(dp), allocatable :: kpoint(:,:), ekpoints(:,:)

    integer :: transDir
    type(TEqPoint) :: eqv_points_iQ

    ! Compute the matrix Kmat as lookup table
    nCentralAtoms = this%basis%nCentralAtoms

    transDir = this%basis%transportDirection
    zmin = minval(this%basis%x(transDir,1:nCentralAtoms))
    zmax = maxval(this%basis%x(transDir,1:nCentralAtoms))

    nDeltaZ = nint((zmax - zmin)/this%dz)
    if (allocated(this%Kmat)) then
       call log_deallocate(this%Kmat)
    end if
    call log_allocate(this%Kmat, nDeltaZ+1, size(this%kweight), size(this%kweight))

    !compute the absolute k-points
    allocate(kpoint(3,size(this%kweight)))
    if (all(this%basis%lattVecs == 0.0_dp)) then
       kpoint = 0.0_dp
    else
       kpoint = matmul(this%recVecs2p, this%kpoint)
    end if

    iQloop: do iQ = 1, size(this%kweight)
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

      ! Add contribution of kpoints equivalent to iQ
      if (allocated(this%equivalent_kpoints)) then
        !Select the set of points equivalent to iQ
        eqv_points_iQ = this%equivalent_kpoints%EqPoints(iQ)
        n_eq = size(eqv_points_iQ%points, 2)
        !Compute the absolute k-points
        call log_allocate(ekpoints, 3, n_eq)
        ekpoints = matmul(this%recVecs2p, eqv_points_iQ%points)

        do eQ = 1, n_eq
          ekp = ekpoints(:,eQ)
          ekp2 = dot_product(ekp,ekp)
          if (ekp2==0.0_dp) then
            this%Kmat(:,:,iQ) = this%Kmat(:,:,iQ) + 0.0_dp
            cycle
         end if

         do iK = 1, size(this%kweight)
          kk = kpoint(:, iK)
          kk2 = dot_product(kk,kk)
          if (kk2==0.0_dp) then
             this%Kmat(:,iK,iQ) = this%Kmat(:,iK,iQ) + 0.0_dp
             cycle
          end if
          QQ = dot_product(kk, ekp)
          Q2 = QQ*QQ
          do iZ = 0, nDeltaZ
            z_mn = iZ * this%dz
            Kf = this%D0*Q2*exp(-sqrt(ekp2)*z_mn)/(2.0_dp*ekp2*kk2)

            this%Kmat(iZ+1,iK,iQ) = this%Kmat(iZ+1,iK,iQ) + this%coupling * this%kweight(iQ) * Kf
          end do
         end do

        end do
        call log_deallocate(ekpoints)
      endif

    end do iQloop
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
    if (allocated(this%kindices)) call log_deallocate(this%kindices)
    if (allocated(this%Kmat)) call log_deallocate(this%Kmat)
    call this%destroy_sigma_r()
    call this%destroy_sigma_n()
    call destroy(this%equivalent_kpoints)

  end subroutine destroy_elph

  !--------------------------------------------------------------------------
  !> Retrieve matrix pointer from cache and add it to ESH
#:def add_sigma_r_template(KIND,MTYPE)  
  subroutine add_sigma_r_${KIND}$(this, esh, en_index, k_index, spin)
    class(ElPhonInel) :: this
    type(${MTYPE}$), dimension(:,:), intent(inout) :: esh
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    integer :: npl, jj
    type(x_dns), pointer :: tmp_blk
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
    class(ElPhonInel) :: this
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
    class(ElPhonInel) :: this
    type(z_dns), dimension(:,:), intent(inout) :: blk_sigma_n
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    type(TMatLabel) :: label
    integer :: ii

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
    integer :: ii, ibl, nbl, Np, Mp, NK, NE, NKloc, NEloc, iEshift
    integer :: iK, iE, Ndz, PL_start
    complex(lp) :: fac_min, fac_plus
    type(TMatPointer), allocatable :: pGG(:,:), pSigma(:,:)
    integer, allocatable :: izr(:), izc(:)
    type(TMatLabel) :: label
    integer :: transDir

    transDir = this%basis%transportDirection
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
    label%spin = 1
    if (present(spin)) then
       label%spin = spin
    end if

    allocate(pGG(NEloc,NKloc))
    allocate(pSigma(NEloc,NKloc))

    select type(p => this%sigma_r)
    class is(TMatrixCacheMem)
       if (.not.p%isInitialized) then
          call p%init(nEloc, nKloc, nbl, 3, 1)
       end if
    end select

    ! Compute the diagonal blocks
    do ibl = 1, nbl
      ! ==================================================================================================
      !  DIAGONAL BLOCKS
      ! ==================================================================================================
      ! block dimension
      PL_start = this%struct%mat_PL_start(ibl)
      Np = this%struct%mat_PL_end(ibl) - PL_start + 1

      ! Project atom position on the coarser grid
      call log_allocate(izr,Np)
      do ii = 1, Np
        izr(ii) = nint(this%basis%x(transDir, this%basis%matrixToBasis(PL_start+ii-1))/this%dz)
      end do
      label%row_block = ibl
      label%col_block = ibl

      call setup_pointers_Gr()
      call setup_pointers_Sigma_r(Np,Np)

      ! Compute the retarded part
      fac_min = cmplx(this%Nq + 1, 0.0, dp)
      fac_plus = cmplx(this%Nq, 0.0, dp)

      call selfenergy(this%cart_comm, pGG, fac_min, fac_plus, iEshift, &
                          & izr, izr, this%Kmat, NK, NKloc, NE, NEloc, this%kindices, pSigma)

      call setup_pointers_Gn()
      ! Compute the Gn part to Sigma_r
      fac_min = (0.0_dp, 0.5_dp)
      fac_plus = (0.0_dp, -0.5_dp)

      call selfenergy(this%cart_comm, pGG, fac_min, fac_plus, iEshift, &
                          & izr, izr, this%Kmat, NK, NKloc, NE, NEloc, this%kindices, pSigma)

      ! ==================================================================================================
      !  UPPER/LOWER TRI-DIAGONAL BLOCKS
      ! ==================================================================================================
      if (this%tTridiagonal .and. ibl < nbl) then
        PL_start = this%struct%mat_PL_start(ibl+1)
        Mp = this%struct%mat_PL_end(ibl+1) - PL_start + 1
        ! Project atom position on the coarser grid (two indep. arrays for rows and cols)
        call log_allocate(izc,Mp)
        do ii = 1, Mp
          izc(ii) = nint(this%basis%x(transDir, this%basis%matrixToBasis(PL_start+ii-1))/this%dz)
        end do

        ! -------------------- supradiagonal blocks -----------------------------------------------------
        label%row_block = ibl
        label%col_block = ibl + 1

        call setup_pointers_Gr()
        call setup_pointers_Sigma_r(Np,Mp)

        ! Compute the retarded part
        fac_min = cmplx(this%Nq + 1, 0.0, dp)
        fac_plus = cmplx(this%Nq, 0.0, dp)

        call selfenergy(this%cart_comm, pGG, fac_min, fac_plus, iEshift, &
                          & izr, izc, this%Kmat, NK, NKloc, NE, NEloc, this%kindices, pSigma)

        call setup_pointers_Gn()
        ! Compute the Gn part to Sigma_r
        fac_min = (0.0_dp, 0.5_dp)
        fac_plus = (0.0_dp, -0.5_dp)

        call selfenergy(this%cart_comm, pGG, fac_min, fac_plus, iEshift, &
                          & izr, izc, this%Kmat, NK, NKloc, NE, NEloc, this%kindices, pSigma)

        ! -------------------- subdiagonal blocks -----------------------------------------------------
        label%row_block = ibl + 1
        label%col_block = ibl

        call setup_pointers_Gr()
        call setup_pointers_Sigma_r(Mp,Np)

        ! Compute the retarded part
        fac_min = cmplx(this%Nq + 1, 0.0, dp)
        fac_plus = cmplx(this%Nq, 0.0, dp)

        call selfenergy(this%cart_comm, pGG, fac_min, fac_plus, iEshift, &
                          & izc, izr, this%Kmat, NK, NKloc, NE, NEloc, this%kindices, pSigma)

        call setup_pointers_Gn()
        ! Compute the Gn part to Sigma_r
        fac_min = (0.0_dp, 0.5_dp)
        fac_plus = (0.0_dp, -0.5_dp)

        call selfenergy(this%cart_comm, pGG, fac_min, fac_plus, iEshift, &
                          & izc, izr, this%Kmat, NK, NKloc, NE, NEloc, this%kindices, pSigma)

        call log_deallocate(izc)
      end if
      ! ==================================================================================================

      call log_deallocate(izr)
    end do

    deallocate(pGG)
    deallocate(pSigma)
#:else
    call error_msg("inelastic scattering requires MPI compilation")
#:endif
    contains

    ! setup the array of pointers to G_r and sigma_r
    subroutine setup_pointers_Gr()

      do iK = 1, NKloc
        do iE = 1, NEloc
          label%kpoint = iK
          label%energy_point = iE
          call this%G_r%retrieve_pointer(pGG(iE,iK)%pMat, label)
        end do
      end do
    end subroutine setup_pointers_Gr

    ! setup the array of pointers to G_n
    subroutine setup_pointers_Gn()
      do iK = 1, NKloc
        do iE = 1, NEloc
          label%kpoint = iK
          label%energy_point = iE
          call this%G_n%retrieve_pointer(pGG(iE,iK)%pMat, label)
        end do
      end do
    end subroutine setup_pointers_Gn


    subroutine setup_pointers_Sigma_r(Nr,Nc)
      integer, intent(in) :: Nr, Nc

      type(x_DNS) :: tmp

      associate(Sigma_r => this%Sigma_r)
      do iK = 1, NKloc
        do iE = 1, NEloc
          label%kpoint = iK
          label%energy_point = iE
          if (.not.Sigma_r%is_cached(label)) then
            call create(tmp,Nr,Nc)
            tmp%val = (0.0_dp, 0.0_dp)
            call Sigma_r%add(tmp, label)
            call destroy(tmp)
          end if
          call Sigma_r%retrieve_pointer(pSigma(iE,iK)%pMat, label)
        end do
      end do
      end associate

    end subroutine setup_pointers_Sigma_r

    subroutine check_elements(Mat,arg)
       class(TMatrixCache) :: Mat
       character(*) :: arg
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

  end subroutine compute_Sigma_r


  !> Give the Gn at given energy point to the interaction
  subroutine compute_Sigma_n(this, en_index, k_index, spin)
    class(ElPhonInel) :: this
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

#:if defined("MPI")
    integer :: ii, ibl, nbl, Np, Mp, NK, NE, NKloc, NEloc, iEshift
    integer :: iK, iE, Ndz, PL_start
    complex(lp) :: fac_min, fac_plus
    type(TMatPointer), allocatable :: pGG(:,:), pSigma(:,:)
    !type(z_DNS), pointer :: pMat
    !type(z_DNS) :: Sigma_n
    integer, allocatable :: izr(:), izc(:)
    type(TMatLabel) :: label
    integer :: transDir

    transDir = this%basis%transportDirection
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

    select type(p => this%sigma_n)
    class is(TMatrixCacheMem)
       if (.not.p%isInitialized) then
          call p%init(nEloc, nKloc, nbl, 3, 1)
       end if
    end select

    do ibl = 1, nbl
      ! ==================================================================================================
      !  Compute the diagonal blocks
      ! ==================================================================================================
      ! block dimension
      PL_start = this%struct%mat_PL_start(ibl)
      Np = this%struct%mat_PL_end(ibl) - PL_start + 1

      ! Project atom position on the coarser grid
      call log_allocate(izr,Np)
      do ii = 1, Np
        izr(ii) = nint(this%basis%x(transDir, this%basis%matrixToBasis(PL_start+ii-1))/this%dz)
      end do
      label%row_block = ibl
      label%col_block = ibl

      call setup_pointers_Gn()
      call setup_pointers_Sigma_n(Np,Np)

      ! Compute the retarded part
      fac_min = cmplx(this%Nq, 0.0, dp)
      fac_plus = cmplx(this%Nq + 1, 0.0, dp)

      call selfenergy(this%cart_comm, pGG, fac_min, fac_plus, iEshift, &
                          & izr, izr, this%Kmat, NK, NKloc, NE, NEloc, this%kindices, pSigma)

      !call check_elements(this%Sigma_n,"Sigma_n")
      ! ==================================================================================================
      !   UPPER TRI- DIAGONAL BLOCKS
      ! ==================================================================================================
      if (this%tTridiagonal .and. ibl < nbl) then
        PL_start = this%struct%mat_PL_start(ibl+1)
        Mp = this%struct%mat_PL_end(ibl+1) - PL_start + 1
        ! Project atom position on the coarser grid (two indep. arrays for rows and cols)
        call log_allocate(izc,Mp)
        do ii = 1, Mp
          izc(ii) = nint(this%basis%x(transDir, this%basis%matrixToBasis(PL_start+ii-1))/this%dz)
        end do

        ! -------------------- supradiagonal blocks -----------------------------------------------------
        label%row_block = ibl
        label%col_block = ibl + 1

        call setup_pointers_Gn()
        call setup_pointers_Sigma_n(Np,Mp)

        ! Compute the retarded part
        fac_min = cmplx(this%Nq, 0.0, dp)
        fac_plus = cmplx(this%Nq + 1, 0.0, dp)

        call selfenergy(this%cart_comm, pGG, fac_min, fac_plus, iEshift, &
                          & izr, izc, this%Kmat, NK, NKloc, NE, NEloc, this%kindices, pSigma)

        call log_deallocate(izc)
      end if
      ! ==================================================================================================

      call log_deallocate(izr)
    end do

    deallocate(pGG)
    deallocate(pSigma)
#:else
    call error_msg("inelastic scattering requires MPI compilation")
#:endif

    contains

    subroutine setup_pointers_Gn()
      do iK = 1, NKloc
        do iE = 1, NEloc
          label%kpoint = iK
          label%energy_point = iE
          call this%G_n%retrieve_pointer(pGG(iE,iK)%pMat, label)
        end do
      end do
    end subroutine setup_pointers_Gn

    subroutine setup_pointers_Sigma_n(Nr,Nc)
      integer, intent(in) :: Nr, Nc
      type(x_DNS) :: tmp

      associate(Sigma_n => this%Sigma_n)
      do iK = 1, NKloc
        do iE = 1, NEloc
          label%kpoint = iK
          label%energy_point = iE
          if (.not.Sigma_n%is_cached(label)) then
            call create(tmp,Nr,Nc)
            tmp%val = (0.0_dp, 0.0_dp)
            call Sigma_n%add(tmp, label)
            call destroy(tmp)
          end if
          call Sigma_n%retrieve_pointer(pSigma(iE,iK)%pMat, label)
        end do
      end do
      end associate
    end subroutine setup_pointers_Sigma_n

    subroutine check_elements(Mat,arg)
      class(TMatrixCache) :: Mat
      character(*) :: arg

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

  !> Destroy Sigma_r
  subroutine destroy_Sigma_r(this)
    class(ElPhonInel) :: this
    if (associated(this%sigma_r)) then
      call this%sigma_r%destroy()
    end if
  end subroutine destroy_Sigma_r

  !> Destroy Sigma_n
  subroutine destroy_Sigma_n(this)
    class(ElPhonInel) :: this
    if (associated(this%sigma_n)) then
      call this%sigma_n%destroy()
    end if
  end subroutine destroy_Sigma_n

  subroutine set_reciprocal_vectors(basis, reciprocal_mat)
    type(TBasisCenters), intent(in) :: basis
    real(dp), intent(inout) :: reciprocal_mat(3,3)

    real(dp), dimension(3) :: latt_vec1, latt_vec2, latt_vec3, rec_vec1, rec_vec2, rec_vec3
    real(dp), dimension(3) :: cross_prod
    real(dp) :: vol

    latt_vec1 = basis%lattVecs(:, 1)
    latt_vec2 = basis%lattVecs(:, 2)
    latt_vec3 = basis%lattVecs(:, 3)

    vol = volume(latt_vec1, latt_vec2, latt_vec3)

    call cross_product(latt_vec2, latt_vec3, cross_prod)
    rec_vec1 = (2.0_dp*pi*cross_prod) / vol
    call cross_product(latt_vec3, latt_vec1, cross_prod)
    rec_vec2 = (2.0_dp*pi*cross_prod) / vol
    call cross_product(latt_vec1, latt_vec2, cross_prod)
    rec_vec3 = (2.0_dp*pi*cross_prod) / vol

    reciprocal_mat(:, 1) = rec_vec1
    reciprocal_mat(:, 2) = rec_vec2
    reciprocal_mat(:, 3) = rec_vec3

  end subroutine set_reciprocal_vectors

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

end module elphinel
