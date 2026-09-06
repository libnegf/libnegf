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
  use mat_def, only : c_DNS, z_DNS
  use mat_def, only : x_DNS => z_DNS, create, destroy, assignment(=)
  use ln_math, only : cross_product, volume
  use ln_constants, only : pi
  use interactions, only : TInteraction
  use ln_inelastic, only : TInelastic
  use ln_allocation, only : log_allocate, log_deallocate, writeMemInfo
  use ln_structure, only : TStruct_info
  use ln_structure, only : TBasisCenters, create_TBasis, destroy_TBasis
  use distributions, only : bose
  use ln_cache
  use iso_c_binding
  use mpi_f08, MPI_XX_COMPLEX => MPI_DOUBLE_COMPLEX
  use mpi_globals
  use inversions, only : inverse
  use equiv_kpoints, only : TEqPoint, TEqPointsArray, create, destroy
  use clock
  use ln_messages, only : error_msg
  implicit none
  private

  public :: ElPhonInel, ElPhonPolarOptical, ElPhonNonPolarOptical
  public :: ElPhonPO_create, ElPhonNonPO_create
  public :: ElPhonPO_init, ElPhonNonPO_init
  
  type, public :: TMatPointer
    type(x_DNS), pointer :: pMat
  end type TMatPointer


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
    procedure :: init_sigma_r
    procedure :: init_sigma_n
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
        & "Electron-Phonon inelastic model for non-PO phonons"

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

    call set_reciprocal_vectors(basis, this%recVecs2p)
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
  ! PO phonons: see also derivation slides:
  !
  ! Sigma_mn(k,E) = Sum_qt',G Sum_qz |U(Q+G)|^2 e^(i*qz*|zm - zn|) * Gmn(kt-qt', E-E')*D(E')
  ! kt-qt' in the BZ  => kt-qt'+G' with G' that does the appropriate folding into the BZ
  ! |U(q)|^2 = C*q^2/(q^2 + q0^2)^2
  ! q = qt' + qz z^  (q is a 3D wave-vector, qt' is 2D)
  ! Set qt=kt-qt'+G'  =>  qt'=kt-qt+G' => q^2 = kt^2 + (qt-G)^2 - 2kt*(qt-G) + qz^2
  ! The sum over G absorbs G' (contains G'). Also note that Sum_{G} == Sum_{-G}
  !
  ! First perform F = Sum_qz  (Q^2 + qz^2)/[ Q^2 + q0^2 + qz^2 ] exp(i*qz*|zm-zn|)  Q=k-q+G
  ! F = F(|zm-zn|, Q^2) = K(|zm-zn|, k, q+G)
  !
  ! Sigma_mn(k,E) = Sum_q Sum_G K(|zm - zn|, k, q+G) [B(-)*Gmn(q,E-wq) + B(+)*Gmn(q,E+wq)]
  !
  ! k <-> -k   symmetry:
  ! H(k)=H(-k) => G(q)=G(-q)  Since F depends on Q^2, F(Q)=F(-Q)
  ! => K(k,q)=K(-k,-q)  and  K(k,-q) = K(-k, q)
  ! Sigma(+k) = Sum_{q\2} K(k,q) G(q) + Sum_{-q\2} K(k,-q) G(q)
  !                       =                          =
  ! Sigma(-k) = Sum_{-q\2} K(-k,-q) G(q) + Sum_{q\2} K(-k,q) G(q)
  !
  ! => Define:  K(|zm-zn|, k, q)  = Sum_G [F(|zm-zn|, (k-q+G)^2) + F(|zm-zn|, (k+q+G)^2)]
  ! ^
  ! |o o x x|   4x4 k-sampling on the [0..1] x [0..1]
  ! |o o x x|   The x are what DFTB+ removes for symmetry. Notice it is NOT k->-k symmetry
  ! |o o x x|   The BZ does not extend -Pi/2..Pi/2 but [0..Pi]
  ! |o o x x|
  ! ---------->
  ! 0       1
  !--------------------------------------------------------------------------
  subroutine prepare_POKmat(this)
    class(ElPhonPolarOptical) :: this

    integer :: iZ, iQ, eQ, iK, nDeltaZ, nCentralAtoms, n_eq
    real(dp) :: kq(3), dk(3), QQ(3), Q2, bb, z_mn, Kf
    real(dp) :: zmin, zmax
    real(dp), allocatable :: myKmat(:,:)
    integer :: transDir, nKpoints, iG
    type(TEqPoint) :: eqv_points_iQ

    ! Compute the matrix Kmat as lookup table
    nCentralAtoms = this%basis%nCentralAtoms

    transDir = this%basis%transportDirection
    zmin = minval(this%basis%x(transDir,1:nCentralAtoms))
    zmax = maxval(this%basis%x(transDir,1:nCentralAtoms))

    nDeltaZ = nint((zmax - zmin)/this%dz)
    nKpoints = size(this%kweight)
    print*,'DEBUG:  #kpoints', nKpoints
    do iG = 1, nKpoints
      print*,'K= ',this%kpoint(:,iG)
    end do

    if (allocated(this%Kmat)) then
       call log_deallocate(this%Kmat)
    end if
    call log_allocate(this%Kmat, nDeltaZ+1, nKpoints, nKpoints)
    call log_allocate(myKmat, nDeltaZ+1, nKpoints)
    myKmat = 0.0_dp

    this%Ce = this%coupling*this%wq/2.0_dp/this%cell_area* &
          & (1.0_dp/this%eps_inf - 1.0_dp/this%eps0)

    ! K(|zm-zn|, k, q) = Sum_G [I(|zm-zn|,k+G,q)+Sum_qe I(|zm-zn|,k+G,qe)]*w(q)
    ! {G} = rVec2Pi * {(1,0,0),(0,1,0),(-1,0,0),(0,-1,0),(1,1,0),(1,-1,0),...}
    ! ISSUE: The G summation does not converge, so we restrict to first BZ
    iQloop: do iQ = 1, nKpoints
       iKloop: do iK = 1, nKpoints
         dk(:) = this%kpoint(:,iK)-this%kpoint(:,iQ)
         call fold_into_BZ(dk, kQ)
         ! We now get QQ in absolute coordinates
         QQ = matmul(this%recVecs2p, kQ)
         Q2 = dot_product(QQ, QQ)
         bb = sqrt(this%q0*this%q0 + Q2)
         if (id0) print*,'iK, bb:',iK,bb
         do iZ = 0, nDeltaZ
           z_mn = iZ * this%dz
           myKmat(iZ+1,iK) = myKmat(iZ+1,iK) + & 
                 & (2.0_dp*Q2 + this%q0*this%q0*(1.0_dp-bb*z_mn))*exp(-bb*z_mn)/ (4.0_dp*bb**3)
           !Kf = (2.0_dp*Q2 + this%q0*this%q0*(1.0_dp-bb*z_mn))*exp(-bb*z_mn)/ (4.0_dp*bb**3)
           !this%Kmat(iZ+1,iK,iQ) = this%Ce * this%kweight(iQ) * Kf
         end do

         ! Add contribution of kpoints equivalent to iQ
         if (allocated(this%equivalent_kpoints)) then
           !Select the set of points equivalent to iQ
           eqv_points_iQ = this%equivalent_kpoints%EqPoints(iQ)
           n_eq = size(eqv_points_iQ%points, 2)

           do eQ = 1, n_eq
             dk(:) = this%kpoint(:,iK)-eqv_points_iQ%points(:,eQ)
             call fold_into_BZ(dk, kQ)
             QQ = matmul(this%recVecs2p, kQ)
             Q2 = dot_product(QQ, QQ)
             bb = sqrt(this%q0*this%q0 + Q2)
             do iZ = 0, nDeltaZ
               z_mn = iZ * this%dz
               myKmat(iZ+1,iK) = myKmat(iZ+1,iK) + & 
                  & (2.0_dp*Q2 + this%q0*this%q0*(1.0_dp-bb*z_mn))*exp(-bb*z_mn)/ (4.0_dp*bb**3)
             end do
           end do
         end if
       end do iKloop 
       this%Kmat(:,:,iQ) = myKmat(:,:) * this%Ce * this%kweight(iQ)
    end do iQloop

    call log_deallocate(myKmat)

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
  ! This routine generates 2d integer G vectors of size r2=m**2+n**2
  ! It does not take into account the actual ||G|| of recLatticeVecs
  ! Potential issue: The internal predefined size of myG is 128
  !                  which is found ok even for r2=1,000,000
  subroutine set_shell_G(r2, transDir, GG, dummyDir)
    integer, intent(in) :: r2
    integer, intent(in) :: transDir
    integer, intent(in), optional :: dummyDir
    integer, allocatable :: GG(:,:)

    integer :: mmax, iG, m, n, ortDir1, ortDir2
    integer :: counts
    integer :: myG(2,128)
    real :: nr

    ! Special handling of 2D system with dummy periodicity
    ! and G-vectors in 1 dimension only:
    if (present(dummyDir)) then
      if (dummyDir == transDir) then
        error stop "dummyDir == transport Dir ?"
      end if
      if (transDir==1 .and. dummyDir==2) ortDir1=3
      if (transDir==1 .and. dummyDir==3) ortDir1=2
      if (transDir==2 .and. dummyDir==1) ortDir1=3
      if (transDir==2 .and. dummyDir==3) ortDir1=1
      if (transDir==3 .and. dummyDir==1) ortDir1=2
      if (transDir==3 .and. dummyDir==2) ortDir1=1
      nr = int(sqrt(1.0*r2))
      if (nr == 1.0*int(nr)) then
        n = int(nr)
        allocate(GG(3,2))
        GG(transDir,:) = 0
        GG(dummyDir,:) = 0
        GG(ortDir1,1) = n
        GG(ortDir1,2) = -n
      end if
      return
    end if

    select case(transDir)
    case(1)
      ortDir1=2; ortDir2=3
    case(2)
      ortDir1=3; ortDir2=1
    case(3)
      ortDir1=1; ortDir2=2
    end select

    if (allocated(GG)) deallocate(GG)

    mmax=int(sqrt(1.0*r2))
    counts=0

    do m = 0, mmax
      nr = sqrt(1.0*(r2-m**2))
      if (nr == 1.0*int(nr)) then
        n = int(nr)
        counts = counts + 1
        myG(1,counts) = m
        myG(2,counts) = n
        if (n /= 0) then
          counts = counts + 1
          myG(1,counts) = m
          myG(2,counts) = -n
        end if
        if (m /= 0) then
          counts = counts + 1
          myG(1,counts) = -m
          myG(2,counts) = n
          if (n /= 0) then
            counts = counts + 1
            myG(1,counts) = -m
            myG(2,counts) = -n
          end if
        end if
      end if
    end do
    if (counts>0) then
      allocate(GG(3,counts))
      GG(transDir,:) = 0
      do iG = 1, counts
        GG(ortDir1,iG) = myG(1,iG)
        GG(ortDir2,iG) = myG(2,iG)
      end do
    end if

  end subroutine set_shell_G

  !--------------------------------------------------------------------------
  ! Folds a kpoint into the first BZ [-pi/2..pi/2] x [-pi/2..pi/2]
  !
  subroutine fold_into_BZ(kpoint, kpoint_fold)
    real(dp), intent(in) :: kpoint(3)
    real(dp), intent(out) :: kpoint_fold(3)

    if (kpoint(1) > pi/2) then
      kpoint_fold(1) = kpoint(1) - pi
    end if
    if (kpoint(2) > pi/2) then
      kpoint_fold(2) = kpoint(2) - pi
    end if
    if (kpoint(3) > pi/2) then
      kpoint_fold(3) = kpoint(3) - pi
    end if

  end subroutine fold_into_BZ

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
    type(x_DNS), pointer :: tmp_blk
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
    type(x_DNS), pointer :: tmp_blk
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
    type(x_DNS), dimension(:,:), intent(inout) :: blk_sigma_n
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
    type(x_DNS), intent(inout) :: sigma_n
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

  subroutine init_sigma_r(this, spin)
    class(ElPhonInel) :: this
    integer, intent(in), optional :: spin

    integer :: nbl, nEloc, nKloc, iK, iE, ibl
    integer :: PL_start, Np, Mp
    type(TMatLabel) :: label
    type(x_DNS) :: tmp

    nbl = this%struct%num_PLs
    NEloc = this%nE_local
    NKloc = size(this%local_kindex)

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
            tmp%val = (0.0_dp, 0.0_dp)
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
              tmp%val = (0.0_dp, 0.0_dp)
              call Sigma_r%add(tmp, label)
              call destroy(tmp)
            end if

            ! Overdiag blocks
            label%row_block = ibl + 1
            label%col_block = ibl
            if (.not.Sigma_r%is_cached(label)) then
              call create(tmp,Mp,Np)
              tmp%val = (0.0_dp, 0.0_dp)
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
    class(ElPhonInel) :: this
    integer, intent(in), optional :: spin

    integer :: nbl, nEloc, nKloc, iK, iE, ibl
    integer :: PL_start, Np, Mp
    type(TMatLabel) :: label
    type(x_DNS) :: tmp

    nbl = this%struct%num_PLs
    NEloc = this%nE_local
    NKloc = size(this%local_kindex)

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
            tmp%val = (0.0_dp, 0.0_dp)
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
              tmp%val = (0.0_dp, 0.0_dp)
              call Sigma_n%add(tmp, label)
              call destroy(tmp)
            end if

            ! Overdiag blocks
            label%row_block = ibl + 1
            label%col_block = ibl
            if (.not.Sigma_n%is_cached(label)) then
              call create(tmp,Mp,Np)
              tmp%val = (0.0_dp, 0.0_dp)
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

    if (vol > 0.0_dp) then
      call cross_product(latt_vec2, latt_vec3, cross_prod)
      rec_vec1 = (2.0_dp*pi*cross_prod) / vol
      call cross_product(latt_vec3, latt_vec1, cross_prod)
      rec_vec2 = (2.0_dp*pi*cross_prod) / vol
      call cross_product(latt_vec1, latt_vec2, cross_prod)
      rec_vec3 = (2.0_dp*pi*cross_prod) / vol
 
      reciprocal_mat(:, 1) = rec_vec1
      reciprocal_mat(:, 2) = rec_vec2
      reciprocal_mat(:, 3) = rec_vec3
   else
      reciprocal_mat = 0.0_dp
   end if   
  
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

  subroutine selfenergy(comm2d, GG, fac_minus, fac_plus, iEhbaromega, &
                               &  izr, izc, KK, NK, NKloc, NE, NEloc, kindices_map, Sigma)
    type(MPI_Comm), intent(in) :: comm2d
    type(TMatPointer), intent(in) :: GG(0:,0:)
    complex(lp), intent(in) :: fac_minus
    complex(lp), intent(in) :: fac_plus
    integer, intent(in) :: izr(:), izc(:)
    real(lp), intent(in) :: KK(0:,0:,0:)
    integer, intent(in) :: iEhbaromega
    integer, intent(in) :: NK
    integer, intent(in) :: NKloc
    integer, intent(in) :: NE
    integer, intent(in) :: NEloc
    integer, intent(in) :: kindices_map(0:,0:)
    type(TMatPointer), intent(in) :: Sigma(0:,0:)

    ! locals
    integer :: Np, Mp, iQ, iK, iQglo, iQglo2, iKglo, nu
    integer :: iEminus, iEplus, iE, iEglo, iin
    integer :: myid, commsize, ndims = 2
    integer :: dims(2), coords(2), coordsH(2)
    type(x_DNS), target :: rbuff1, rbuff2, sbuffK, rbuffK
    type(x_DNS), pointer :: pbuff1, pbuff2, pGG1, pGG2, pSigma
    real(lp), allocatable :: KKbuf(:)
    integer :: ierr
    type(MPI_Request) :: rqE1,rqE2,rqE3,rqE4, rqK1, rqK2
    type(MPI_Status) :: statE(4), statK(2)

    integer :: msource, mdest, psource, pdest
    integer :: hsource, hdest
    integer :: ndiff

    dims(1) = NK/NKloc
    dims(2) = NE/NEloc

    call MPI_comm_size(comm2d, commsize, ierr)
    call MPI_comm_rank(comm2d, myid, ierr)
    call MPI_cart_coords(comm2d, myid, ndims, coords, ierr)
    call MPI_barrier(comm2d, ierr)

    Np = size(GG(0,0)%pMat%val,1)
    Mp = size(GG(0,0)%pMat%val,2)

    call create(rbuff1, Np, Mp)
    call create(rbuff2, Np, Mp)
    call create(sbuffK, Np, Mp)
    call create(rbuffK, Np, Mp)

    ! Sigma_ij(iK, iE) = Sum_iQ   KK(|z_i-z_j|, iK, iQ) *
    !               * (fac_minus * GG_ij(iQ, E-wq) + fac_plus * GG_ij(iQ, E+wq))
    qloop:do iQ = 0, NKloc-1
      iQglo = kindices_map(iQ,coords(1))-1
      eloop:do iE = 0, NEloc-1
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
          pbuff1 => GG(iEminus, iQ)%pMat
        else !(iEminus < 0)
          rbuff1%val=(0.0_lp, 0.0_lp)
          pbuff1 => rbuff1
        end if

        ! -------------------------------------------------------------------------------
        ! MPI communications.
        ! Communications take place if (iE < iEhbaromega .and. iEglo >= iEhbaromega)
        ! The latter is necessary to avoid 1-to-many communications
        ! If we take G(0) for all G(E<0) it is possible that rank 0 has to send
        ! G's to more than 1 rank, whenever ihbaromega > NEloc
        ! Example behaviour for ihbaromega=7 and NEloc=6, for iE=0:
        ! iEglo=0; ndiff=0 => no communication
        ! iEglo=6; ihbaromega=7; NEloc=6 => ndiff=1; iEminus=0
        ! => rank=0 SEND TO rank=1 GG(0,iQ)
        ! iEglo=12; ihbaromega=7; ndiff=12/6 - 5/6 = 2 => iEminus=(18-7)%6=5
        ! => rank=0 SEND TO rank=2 GG(5,iQ)
        ! We have a 1 TO many communication
        ! -------------------------------------------------------------------------------
        if (dims(2) > 1) then
          if (iE < iEhbaromega) then
            if (iEglo - iEhbaromega < 0 ) then
              rbuff1%val=(0.0_lp, 0.0_lp)
              pbuff1 => rbuff1
              ! takes care of the negative side where -1 turns into 0
              ! so +1 is added at the end if (mod(iEhbaromega-iEglo, NEloc) /= 0)
              !  ______________________
              !  333333222222221111111100000000111111112222222233333333
              ! --ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo----------
              !   345670123456701234567012345670123456701234567012345670
              !        ^-------------------^
              ! -----------------------------------------------------------------------------
              !rank=           0 iE           4 iEglo           4 iEminus         -16     -16
              !ndiff = (ihwq-iEglo)/NEloc + iEglo/NEloc = 16/8 + 0 = 2
              !rank=           2 iE           4 iEglo          20 iEminus         -16      0
              !ndiff = iEglo/NEloc - (iEglo-ihwq)/NEloc = 20/8 - (0) = 2
              ! -----------------------------------------------------------------------------
              !rank=           0 iE           5 iEglo           5 iEminus         -15     -15
              !ndiff = (ihwq-iEglo)/NEloc + iEglo/NEloc + 1 = 15/8 + 0 + 1 = 2
              !rank=           2 iE           5 iEglo          21 iEminus         -15      1
              !ndiff = iEglo/NEloc - (iEglo-ihwq)/NEloc = 21/8 - 1/8  = 2
              ! -----------------------------------------------------------------------------
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
               pGG1 => GG(iEminus, iQ)%pMat
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
        !// pbuff2 => G(k,E+wq)
        pdest = MPI_PROC_NULL
        psource = MPI_PROC_NULL
        iEplus = iE + iEhbaromega
        if (iEplus < NEloc) then
          pbuff2 => GG(iEplus, iQ)%pMat
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
               pGG2 => GG(iEplus, iQ)%pMat
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
        !//   sbuffK_ij(iQ) = (fac_minus * GG_ij(iQ, E-wq) + fac_plus * GG_ij(iQ, E+wq))

        !$OMP PARALLEL DO private (nu)
        do nu = 1, Mp
          sbuffK%val(:,nu) = (fac_minus*pbuff1%val(:,nu) + fac_plus*pbuff2%val(:,nu))
        end do
        !$OMP END PARALLEL DO

        !// The Outer loop performs the summation over iQ
        !// This Inner loop updates Sigma(iK):
        !// Sigma_ij(iK, iE) = Sum_iQ   KK(|z_i-z_j|, iK, iQ) * sbuffK_ij(iQ)
        do iK = 0, NKloc-1
          iKglo = kindices_map(iK,coords(1))-1
          pSigma => Sigma(iE, iK)%pMat

          !$OMP PARALLEL private (nu, KKbuf)
          allocate(KKbuf(Np))
          !$OMP DO
          do nu = 1, Mp
            KKbuf(:) = KK(abs(izr(1:Np)-izc(nu)), iKglo, iQglo)
            pSigma%val(:,nu) = pSigma%val(:,nu) + KKbuf(:) * sbuffK%val(:,nu)
          end do
          !$OMP END DO
          deallocate(KKbuf)
          !$OMP END PARALLEL
        end do


        !////////////////////////////////////////////////////////////////////////////////////
        !// Communications over the k-grid (Note k-grid is periodic)
        !// Suppose iQglo = 0, 1, 2, 3, 4, 5
        !// dims(1) == 3 =>
        !// cpu#0: iQ = 0, 1; coords(1)=0 => iQglo = 0, 1
        !// cpu#1: iQ = 0, 1; coords(1)=1 => iQglo = 2, 3
        !// cpu#2: iQ = 0, 1; coords(1)=2 => iQglo = 4, 5
        !//
        !// ii = 1: #0->#1, #1->#2,  #2->#0
        !// cpu#0: (recv sbuffK(iQglo=4,5) from cpu#2  coordsH(1)=2) => iQglo2 = 4,5  ()
        !// cpu#1: (recv sbuffK(iQglo=0,1) from cpu#0  coordsH(1)=0) => iQglo2 = 0,1  ()
        !// cpu#2: (recv sbuffK(iQglo=2,3) from cpu#1  coordsH(1)=1) => iQglo2 = 2,3  ()
        !// ii = 2: #0->#2, #1->#0,  #2->#1
        !// cpu#0: (recv sbuffK(iQglo=2,3) from cpu#1  coordsH(1)=1) => iQglo2 = 2,3  ()
        !// cpu#1: (recv sbuffK(iQglo=4,5) from cpu#2  coordsH(1)=2) => iQglo2 = 4,5  ()
        !// cpu#2: (recv sbuffK(iQglo=0,1) from cpu#0  coordsH(1)=0) => iQglo2 = 0,1  ()
        !////////////////////////////////////////////////////////////////////////////////////
        if (dims(1) > 1) then
          do iin = 1, dims(1)-1
            ! Compute the communication points applying a shift in dimension 0 (k-grid)
            ! NOTE: this means a communication of 1->2,3,4,5... 2->3,4,5,6...
            !       maybe a nearest-neighbours communication 1->2->3->4 ...
            !       is more efficient ? => sbuff = rbuff and send to neighbour
            ! First shift used to compute iQglo2:
            call MPI_Cart_shift(comm2d, 0, iin, hsource, hdest, ierr)
            call MPI_Cart_coords(comm2d, hsource, ndims, coordsH, ierr)
            iQglo2 = kindices_map(iQ,coordsH(1))-1

            ! Actual shift used in send/recv:
            call MPI_Cart_shift(comm2d, 0, 1, hsource, hdest, ierr)

            call MPI_Isend(sbuffK%val, Mp*Np, MPI_XX_COMPLEX, hdest, 43, comm2d, rqK1, ierr)
            call MPI_Irecv(rbuffK%val, Mp*Np, MPI_XX_COMPLEX, hsource, 43, comm2d, rqK2, ierr)

            call MPI_Wait(rqK2, statK(2), ierr)

            ! Sigma_ij(iK, iE) = Sum_iQ   KK(|z_i-z_j|, iK, iQ2) * rbuffK_ij(iQ2)
            do iK = 0, NKloc-1
              iKglo = kindices_map(iK,coords(1))-1
              pSigma => Sigma(iE, iK)%pMat
              !$OMP PARALLEL private (nu, KKbuf)
              allocate(KKbuf(Np))
              !$OMP DO
              do nu = 1, Mp
                KKbuf(:) = KK(abs(izr(1:Np)-izc(nu)), iKglo, iQglo2)
                pSigma%val(:,nu) = pSigma%val(:,nu) + KKbuf(:) * rbuffK%val(:,nu)
              end do
              !$OMP END DO
              deallocate(KKbuf)
              !$OMP END PARALLEL
            end do

            call MPI_Wait(rqK1, statK(1), ierr)

            !$OMP PARALLEL DO
            do nu = 1, Mp
              sbuffK%val(:,nu) = rbuffK%val(:,nu)
            end do
            !$OMP END PARALLEL DO

          end do
        end if

      end do eloop
    end do qloop

    call destroy(rbuff1)
    call destroy(rbuff2)
    call destroy(sbuffK)
    call destroy(rbuffK)


  end subroutine selfenergy

end module elphinel
