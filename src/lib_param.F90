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


module lib_param

  use iso_c_binding
  use ln_precision, only : dp
  use globals
  use mat_def
  use ln_structure, only : TStruct_info, TBasisCenters, TNeighbourMap
  use ln_messages
  use input_output
  use energy_mesh, only : mesh
  use interactions, only : TInteraction, TInteractionList, TInteractionNode
  use elphdd, only : ElPhonDephD, ElPhonDephD_create, ElPhonDephD_init
  use elphdb, only : ElPhonDephB, ElPhonDephB_create, ElPhonDephB_init
  use elphds, only : ElPhonDephS, ElPhonDephS_create, ElPhonDephS_init
  use elphinel, only : ElPhonPolarOptical, ElPhonPO_create, ElPhonPO_init
  use elphinel, only : ElPhonNonPolarOptical, ElPhonNonPO_create, ElPhonNonPO_init
  use equiv_kpoints, only : TEqPointsArray
  use scba
  use ln_cache
  use ln_enums, only : integration_type
#:if defined("MPI")
  use libmpifx_module, only : mpifx_comm
#:endif
  implicit none
  private

  public :: Tnegf, lnParams, intArray, complArray, realArray, TEnGrid, TInteractionList
  public :: TMatrixArrayHS, TMatrixArrayDM
  public :: set_defaults, print_all_vars
  public :: destroy_interactions
  public :: init_cache_space, destroy_cache_space
  public :: set_bp_dephasing
  public :: set_elph_dephasing
  public :: set_elph_block_dephasing
  public :: set_elph_s_dephasing
  public :: set_elph_polaroptical
  public :: set_elph_nonpolaroptical
  integer, public, parameter :: MAXNCONT=10
  integer, public, parameter :: MAXNLAYERS=10000

  type intArray
    integer, dimension(:), allocatable :: indexes
  end type intArray

  type complArray
    complex(dp), dimension(:), allocatable :: array
  end type complArray

  type realArray
    real(dp), dimension(:), allocatable :: array
  end type realArray

  ! Array used to store the k-dependent matrices
  type TMatrixArrayHS
    type(z_CSR), pointer :: H => null()    ! Points to externally allocated H
    type(z_CSR), pointer :: S => null()
    logical    :: isSid = .false.          ! True if overlap S == Id
    logical    :: internalHS = .true.      ! tells HS are internally allocated
  end type TMatrixArrayHS

  ! Array used to store the k-dependent matrices
  type TMatrixArrayDM
    type(z_CSR), pointer :: rho => null()      ! Holding output Matrix
    logical    :: internalDM = .true.          ! tells DM is internally allocated
  end type TMatrixArrayDM

#:if defined("GPU")
  !Type definition of cublas and cusolver handles (no module defs?)
  type, bind(C) :: cublasHandle
    type(c_ptr) :: handle
  end type cublasHandle

  type, bind(C) :: cusolverDnHandle
    type(c_ptr) :: handle
  end type cusolverDnHandle

  public :: cusolverDnHandle
  public :: cublasHandle
#:endif

  !! Structure used to define energy points for the integration
  !! For every point we define
  !!     path (1,2 or 3): the energy point belongs to a real axis
  !!     integration (1), a complex plane integration (2) or a
  !!     pole summation (3)
  !!     pt_path: relative point number within a single path
  !!     pt: absolute point number along the whole integration path
  !!     cpu: cpu assigned to the calculation of the given energy point
  !!     Ec: energy value
  !!     wght: a weight used in final summation to evaluate integrals
  type TEnGrid
     integer :: path
     integer :: pt_path
     integer :: pt
     integer :: pt_cpu
     integer :: cpu
     complex(dp) :: Ec
     complex(dp) :: wght
  end type TEnGrid

  !-----------------------------------------------------------------------------
  ! Buttiker probe dephasing
  type :: Tdeph_bp
     real(dp), allocatable, dimension(:) :: coupling
  end type Tdeph_bp

  !-----------------------------------------------------------------------------
  ! Structure to hold contact and reservoir data
  type Tcontact
    character(132) :: name
    logical :: tWriteSelfEnergy = .false.
    logical :: tReadSelfEnergy = .false.
    logical :: tWriteSurfaceGF = .false.
    logical :: tReadSurfaceGF = .false.
    complex(dp), dimension(:,:,:), allocatable :: SelfEnergy ! Electrode Self-Energy
    complex(dp), dimension(:,:,:), allocatable :: SurfaceGF  ! Electrode Surface Green Function
    real(dp) :: mu_n        ! Electrochemical potential (el)
    real(dp) :: mu_p        ! Electrochemical potential (hl)
    real(dp) :: mu          ! Electrochemical Potential (dft calculation)
    real(dp) :: contact_DOS ! Ficticious contact DOS
    logical  :: FictCont    ! Ficticious contact
    real(dp) :: kbT_dm      ! Electronic temperature
    real(dp) :: kbT_t       ! Electronic temperature
    type(z_DNS) :: HC       ! Contact Hamiltonian
    type(z_DNS) :: SC       ! Contact Overlap
    type(z_DNS) :: HMC      ! Device-Contact Hamiltonian
    type(z_DNS) :: SMC      ! Device-Contact Overlap
  end type Tcontact

  !> General libnegf container
  !! Contains input data, runtime quantities and output data
  type Tnegf
    !! Input parameters: set by library user
    !! General
    integer :: verbose
#:if defined("MPI")
    type(mpifx_comm) :: globalComm
    type(mpifx_comm) :: cartComm
    type(mpifx_comm) :: energyComm
    type(mpifx_comm) :: kComm
#:endif
#:if defined("GPU")
   type(cublasHandle) :: hcublas
   type(cusolverDnHandle) :: hcusolver
   integer :: devnum            ! Device Number
   integer(C_SIZE_T) :: freemem   ! Amount of free mem on GPU
   integer(C_SIZE_T) :: totalmem  ! Total memory on GPU
#:endif
    integer  :: ReadoldDM_SGFs    ! 0: Read 1: compute 2: comp & save
    integer  :: ReadoldT_SGFs     ! 0: Read 1: compute 2: comp & save
    character(len=LST) :: scratch_path    ! Folder for scratch work
    character(len=LST) :: out_path        ! Folder for output data
    real(dp) :: g_spin            ! spin degeneracy
    real(dp) :: delta             ! delta for G.F.
    real(dp) :: dos_delta         ! additional delta to force more broadening in the DOS
    integer  :: deltaModel        ! Used for phonon G.F. ! delta**2, 2*delta*w, Mingo's
    real(dp) :: wmax              ! Maximum frequency in Mingo's model
                                  ! See 'Numerical Heat Transfer, Part B', 51:333, 2007
    real(dp) :: eneconv           ! Energy conversion factor
    integer  :: spin = 1          ! spin component
    character(1) :: DorE          ! Density or En.Density

    !! Contacts info
    type(Tcontact), dimension(:), allocatable :: cont

    !! Contour integral
    integer :: Np_n(2)            ! Number of points for n
    integer :: Np_p(2)            ! Number of points for p
    integer :: Np_real            ! Number of points for integration over real axis
    integer :: n_kt               ! Number of kT extending integrations
    integer :: n_poles            ! Number of poles
    real(dp) :: Ec                ! conduction band edge
    real(dp) :: Ev                ! valence band edge

    !! Real axis
    real(dp) :: Emin              ! Tunneling or dos interval
    real(dp) :: Emax              !
    real(dp) :: Estep             ! Tunneling or dos E step
    real(dp) :: Estep_coarse      ! dos E step for coarse integration (quasiEq integral)

    !Particle info for holes/electrons integration
    integer :: particle

    !! Emitter and collector for transmission or Meir-Wingreen
    !! (only emitter in this case)
    integer, allocatable :: ni(:) ! ni: emitter contact list
    integer, allocatable :: nf(:) ! nf: collector contact list

    integer :: ndos_proj              ! number of LDOS interval ranges
    type(intArray), dimension(:), allocatable :: DOS_proj ! PDOS descriptor
                                                          ! contain matrix index
                                                          ! for PDOS projection
    type(intArray) :: TUN_proj    ! Array of TUN projection indices

    real(dp) :: DeltaEc           ! safe guard energy below Ec
    real(dp) :: DeltaEv           ! safe guard energy above Ev

    !! Runtime variables: used internally by the library
    type(format) :: form          ! Form of file-Hamiltonian
    logical  :: dumpHS            ! Used for debug
    real(dp) :: muref             ! reference elec.chem potential
    real(dp) :: E                 ! Holding variable
    real(dp) :: dos               ! Holding variable
    integer :: activecont         ! contact selfenergy
    integer :: min_or_max         ! in input: 0 take minimum, 1 take maximum mu
    integer :: refcont            ! reference contact (for non equilib)
    integer :: outer_blocks       ! flag switching computation of
                                  ! the Device/Contact DM
                                  ! 0 none; 1 upper block; 2 all

    real(dp) :: mu                !    chem potential used without contacts
    real(dp) :: mu_n              ! el chem potential used without contacts
    real(dp) :: mu_p              ! hl chem potential used without contacts
    real(dp) :: kbT               ! temperature used without contacts

    !! Note: H,S are partitioned immediately after input, therefore they are
    !! built runtime from input variable
    class(TMatrixArrayHS), allocatable :: HS(:)
    class(TMatrixArrayDM), allocatable :: DM(:)
    ! Ponters of current working Hamiltonian
    type(z_CSR), pointer :: H => null()
    type(z_CSR), pointer :: S => null()
    ! Ponters of current working density matrix
    type(z_CSR), pointer :: rho => null()
    !logical :: internalDM

    !Bulk contact density calculation
    logical :: bulk_cont_density                    ! Flag to trigger the calculation
    type(z_DNS) :: cont_bulkG(MAXNCONT)             ! Collection of contact bulk Green's functions
    type(complArray) :: bulk_diags(MAXNCONT)        ! Diagonals of the bulk Green that contain
                                                    !   the result of the energy integration (one for each contact)
    type(realArray) :: contact_density(MAXNCONT)    ! Final result of the contact density calculation

    type(TStruct_Info) :: str     ! system structure
    type(TBasisCenters) :: basis  ! local basis centers
    type(TNeighbourMap), dimension(:), allocatable :: neighbour_map
    integer :: iE                 ! Currently processed En point (global index)
    integer :: iE_path            ! Currently processed En point (on path)
    integer :: iEloc              ! Currently processed En point (local index)
    complex(dp) :: Epnt           ! Currently processed En point (value)
    real(dp) :: kwght             ! currently processed k-point weight
    integer :: iKpoint            ! currently processed k-point (global index)
    integer :: iKloc              ! currently processed k-point (local index)

    ! Energy grid information
    type(TEnGrid), dimension(:), allocatable :: en_grid
    integer :: local_en_points    ! Local number of energy points

    ! Array to store all k-points and k-weights
    real(dp), allocatable, dimension(:,:) :: kpoints
    real(dp), allocatable, dimension(:) :: kweights
    ! Array of local k-point indices
    integer, allocatable, dimension(:) :: local_k_index

    ! Information about equivalent kpoints
    type(TEqPointsArray), allocatable :: equivalent_kpoints

    type(mesh) :: emesh           ! energy mesh for adaptive Simpson
    real(dp) :: int_acc           ! adaptive integration accuracy

    ! Many Body Interactions as array of pointers
    type(TInteractionList)  :: interactList
    type(TScbaDriverElastic) :: scbaDriverElastic
    type(TScbaDriverInelastic) :: scbaDriverInelastic

    real(dp) :: scba_elastic_tol = 1.0e-7_dp
    real(dp) :: scba_inelastic_tol = 1.0e-7_dp

    !! Output variables: these arrays are filled by internal subroutines to store
    !! library outputs
    real(dp), dimension(:,:), allocatable :: tunn_mat
    real(dp), dimension(:,:), allocatable :: curr_mat
    real(dp), dimension(:,:), allocatable :: ldos_mat
    real(dp), dimension(:), allocatable :: currents

    ! Integration type for layercurrents or meirwingreen
    integer :: integration = integration_type%trapezoidal

    ! These variables need to be done to clean up
    logical :: tOrthonormal = .false.
    logical :: tOrthonormalDevice = .false.
    integer :: numStates = 0
    logical :: tDestroyGr = .true.
    logical :: tDestroyGn = .true.
    logical :: tDestroyESH = .true.
    ! Buttiker Probes dephasing
    logical :: tDephasingBP = .false.
    logical :: tDephasingVE = .false.
    type(Tdeph_bp) :: bp_deph

    ! internal use only
    integer :: readOldSGF

    ! Work variable: surface green cache.
    class(TMatrixCache), allocatable :: surface_green_cache
    class(TMatrixCache), allocatable :: ESH
    ! These are pointers so they can be passed to inelastic
    class(TMatrixCache), pointer :: G_r => null()
    class(TMatrixCache), pointer :: G_n => null()

    contains

    !procedure :: set_defaults => set_defaults
    !procedure :: print_all_vars => print_all_vars
    !procedure :: destroy_interactions => destroy_interactions
    !procedure :: init_cache_space => init_cache_space
    !procedure :: destroy_cache_space => destroy_cache_space

  end type Tnegf

  !-----------------------------------------------------------------------------
  !> Contains all the general parameters to be passed as input to library
  !! which are compatible with iso_c_binding representations
  !! 1-1 correspondance to input parameter in type Tnegf
  type, bind(c) :: lnparams
    !! General
    !> verbosity, > 100 is maximum
    integer(c_int) :: verbose
    !> Managing SGF readwrite for DM: 0: Read 1: compute 2: comp & save
    integer(c_int)  :: readOldDM_SGFs
    !> Managing SGF readwrite for Tunn: 0: Read 1: compute 2: comp & save
    integer(c_int)  :: readOldT_SGFs
    !> SGF cache destination: 0 for disk, 1 for memory, 2 for a dummy cache (no save performed)
    integer(c_int) :: SGFcache
    !> Spin component (for io)
    integer(c_int)  :: spin
    !> k-point index (for io)
    integer(c_int) :: ikpoint
    !> Spin degeneracy
    real(c_double) :: g_spin
    !> Imaginary delta
    real(c_double) :: delta
    !> delta model for phonon GF
    integer(c_int) :: deltaModel
    !> Maximal energy in Mingo delta model
    real(c_double) :: wmax
    !> Additional delta to force more broadening in the DOS
    real(c_double) :: dos_delta
    !> Energy conversion factor
    real(c_double) :: eneconv
    !> Weight for k-point integration, output integrals are scaled accordingly
    real(c_double) :: kwght
    !> Conduction band edge
    real(c_double) :: ec
    !> Valence band edge
    real(c_double) :: ev
    !> Safe guard energy below Ec
    real(c_double) :: deltaec
    !> Safe guard energy above Ev
    real(c_double) :: deltaev
    !! Real axis integral
    !> Minimum energy for real axis (current integration, DOS, tunneling)
    real(c_double) :: emin
    !> Maximum energy for real axis
    real(c_double) :: emax
    !> Energy step for real axis
    real(c_double) :: estep
    !> Energy step for coarse integrations
    real(c_double) :: estep_coarse
    !! Contacts info
    !> Electron electrochemical potential
    real(c_double) :: mu_n(MAXNCONT)
    !> Hole electrochemical potential
    real(c_double) :: mu_p(MAXNCONT)
    !> Electron electrochemical Potential (for dft)
    real(c_double) :: mu(MAXNCONT)
    !> Contact DOS for WBA
    real(c_double) :: contact_dos(MAXNCONT)
    !> Logical value: is the contact WB?
    logical(c_bool)  :: fictcont(MAXNCONT)
    !> Electronic temperature for each contact (Density Matrix)
    real(c_double) :: kbT_dm(MAXNCONT)
    !> Electronic temperature for each contact (Transmission)
    real(c_double) :: kbT_t(MAXNCONT)
    !> SCBA tolerance for inelastic loop
    real(c_double) :: scba_inelastic_tol
    !> SCBA tolerance for elastic loop
    real(c_double) :: scba_elastic_tol
    !! Contour integral
    !> Number of points for n
    integer(c_int) :: np_n(2)
    !> Number of points for p
    integer(c_int) :: np_p(2)
    !> Number of real axis points
    integer(c_int) :: np_real
    !> ! Number of kT extending integrations
    integer(c_int) :: n_kt
    !> Number of poles
    integer(c_int) :: n_poles
    !! Emitter and collector for transmission or Meir-Wingreen
    !! (only emitter in this case)
    !> Emitter contact list (or lead for integration in MW)
    integer(c_int) :: ni(MAXNLAYERS)
    !> Collector contact list
    integer(c_int) :: nf(MAXNLAYERS)
    !> Should I calculate the density ("D") or the energy weighted density ("E")?
    character(kind=c_char, len=1) :: dore  ! Density or En.Density
    !> Reference contact is set to maximum or minimum Fermi level
    integer(c_int) :: min_or_max
    logical(c_bool) :: is_s_is;
  end type lnparams

contains

  !> Set buttiker probe dephasing
  subroutine set_bp_dephasing(negf, coupling)
    type(Tnegf) :: negf
    real(dp),  dimension(:), intent(in) :: coupling

    if (.not.allocated(negf%bp_deph%coupling)) then
       allocate(negf%bp_deph%coupling(size(coupling)))
    end if
    negf%bp_deph%coupling = coupling

  end subroutine set_bp_dephasing

  !> Set values for the local electron phonon dephasing model
  !! (elastic scattering only)
  subroutine set_elph_dephasing(negf, coupling, niter)
    type(Tnegf), intent(inout) :: negf
    real(dp),  dimension(:), allocatable, intent(in) :: coupling
    integer, intent(in) :: niter

    type(TInteractionNode), pointer :: node

    call negf%interactList%add(node)
    call elphondephd_create(node%inter)
    select type(pInter => node%inter)
    type is(ElPhonDephD)
      call elphondephd_init(pInter, negf%str, coupling, niter)
    class default
      call error_msg( 'ERROR: error of type downcast to ElPhonDephD')
    end select

  end subroutine set_elph_dephasing

  !> Set values for the semi-local electron phonon dephasing model
  !! (elastic scattering only)
  subroutine set_elph_block_dephasing(negf, coupling, orbsperatom, niter)
    type(Tnegf), intent(inout) :: negf
    real(dp),  dimension(:), allocatable, intent(in) :: coupling
    integer,  dimension(:), allocatable, intent(in) :: orbsperatom
    integer, intent(in) :: niter

    type(TInteractionNode), pointer :: node

    call negf%interactList%add(node)
    call elphondephb_create(node%inter)
    select type(pInter => node%inter)
    type is(ElPhonDephB)
      call elphondephb_init(pInter, negf%str, coupling, orbsperatom, niter)
    class default
      call error_msg( 'ERROR: error of type downcast to ElPhonDephB')
    end select

  end subroutine set_elph_block_dephasing

  !> Set values for the semi-local electron phonon dephasing model
  !! (elastic scattering only)
  subroutine set_elph_s_dephasing(negf, coupling, orbsperatom, niter)
    type(Tnegf), intent(inout) :: negf
    real(dp),  dimension(:), allocatable, intent(in) :: coupling
    integer,  dimension(:), allocatable, intent(in) :: orbsperatom
    integer, intent(in) :: niter

    type(TInteractionNode), pointer :: node
    call negf%interactList%add(node)
    call elphondephs_create(node%inter)
    select type(pInter => node%inter)
    type is(ElPhonDephS)
      call elphondephs_init(pInter, negf%str, coupling, orbsperatom, negf%S, niter)
    class default
      call error_msg( 'ERROR: error of type downcast to ElPhonDephS')
    end select

  end subroutine set_elph_s_dephasing

  subroutine set_elph_polaroptical(negf, coupling, wq, Temp, dz, eps0, eps_inf, q0, area, niter, &
              & tridiag)
    type(Tnegf), intent(inout) :: negf
    real(dp),  dimension(:), allocatable, intent(in) :: coupling
    real(dp), intent(in) :: wq
    real(dp), intent(in) :: Temp
    real(dp), intent(in) :: dz
    real(dp), intent(in) :: eps0
    real(dp), intent(in) :: eps_inf
    real(dp), intent(in) :: q0
    real(dp), intent(in) :: area
    integer, intent(in) :: niter
    logical, intent(in) :: tridiag

    type(TInteractionNode), pointer :: node
    call negf%interactList%add(node)
    call elphonpo_create(node%inter)
    select type(pInter => node%inter)
    type is(ElPhonPolarOptical)
#:if defined("MPI")
      call ElPhonPO_init(pInter, negf%cartComm%comm, negf%str, negf%basis, coupling, &
          &  wq, Temp, dz, eps0, eps_inf, q0, area, niter, tridiag)
#:else
      call error_msg( "Inelastic scattering requires MPI")
#:endif
    class default
      call error_msg( 'ERROR: error of type downcast to ElPhonInel')
    end select

  end subroutine set_elph_polaroptical

  subroutine set_elph_nonpolaroptical(negf, coupling, wq, Temp, dz, D0, area, niter, &
              & tridiag)
    type(Tnegf), intent(inout) :: negf
    real(dp),  dimension(:), allocatable, intent(in) :: coupling
    real(dp), intent(in) :: wq
    real(dp), intent(in) :: Temp
    real(dp), intent(in) :: dz
    real(dp), intent(in) :: D0
    real(dp), intent(in) :: area
    integer, intent(in) :: niter
    logical, intent(in) :: tridiag

    type(TInteractionNode), pointer :: node
    call negf%interactList%add(node)
    call elphonnonpo_create(node%inter)
    select type(pInter => node%inter)
    type is(ElPhonNonPolarOptical)
#:if defined("MPI")
    call ElPhonNonPO_init(pInter, negf%cartComm%comm, negf%str, negf%basis, coupling, &
            &  wq, Temp, dz, D0, area, niter, tridiag)
#:else
      call error_msg( "Inelastic scattering requires MPI")
#:endif
    class default
      call error_msg( 'ERROR: error of type downcast to ElPhonInel')
    end select

  end subroutine set_elph_nonpolaroptical


  !> clean interactions objects
  subroutine destroy_interactions(negf)
    type(Tnegf) :: negf
    call negf%interactList%destroy()
  end subroutine destroy_interactions

  subroutine init_cache_space(negf, matrix)
    type(Tnegf) :: negf
    character(3), intent(in) :: matrix

    if (matrix == 'G_r' .and. .not.associated(negf%G_r)) then
       allocate(TMatrixCacheMem::negf%G_r)
       select type( p => negf%G_r)
       type is(TMatrixCacheMem)
         p%tagname='G_r'
       end select
    end if

    if (matrix == 'G_n' .and. .not.associated(negf%G_n)) then
       allocate(TMatrixCacheMem::negf%G_n)
       select type( p => negf%G_n)
       type is(TMatrixCacheMem)
         p%tagname='G_n'
       end select
    end if
  end subroutine init_cache_space

  subroutine destroy_cache_space(negf)
    type(Tnegf) :: negf

    if (associated(negf%G_r)) then
       call negf%G_r%destroy()
       deallocate(negf%G_r)
    end if

    if (associated(negf%G_n)) then
       call negf%G_n%destroy()
       deallocate(negf%G_n)
    end if
  end subroutine destroy_cache_space


  subroutine set_defaults(negf)
     type(Tnegf) :: negf

     negf%verbose = 10

     negf%scratch_path = './GS/'
     negf%out_path = './'
     negf%DorE = 'D'           ! Density or En.Density

     negf%ReadOldDM_SGFs = 1   ! Compute Surface G.F. do not save
     negf%ReadOldT_SGFs = 1    ! Compute Surface G.F. do not save

     negf%kwght = 1.0_dp
     negf%ikpoint = 1

     negf%Ec = 0.0_dp
     negf%Ev = 0.0_dp
     negf%DeltaEc = 0.0_dp
     negf%DeltaEv = 0.0_dp

     negf%E = 0.0_dp            ! Holding variable
     negf%dos = 0.0_dp          ! Holding variable
     negf%eneconv = 1.0_dp      ! Energy conversion factor

     !negf%internalDM = .true.   ! Internal allocation of D.M.

     negf%delta = 1.d-4      ! delta for G.F.
     negf%dos_delta = 0.0_dp  ! delta for DOS
     negf%deltaModel = 1     ! deltaOmega model
     negf%wmax = 0.009_dp     ! about 2000 cm^-1 cutoff
     negf%Emin = 0.0_dp        ! Tunneling or dos interval
     negf%Emax = 0.0_dp        !
     negf%Estep = 0.0_dp       ! Tunneling or dos E step
     negf%Estep_coarse = 0.0_dp       ! Tunneling or dos E step
     negf%g_spin = 2.0_dp      ! spin degeneracy

     negf%Np_n = (/20, 20/)  ! Number of points for n
     negf%Np_p = (/20, 20/)  ! Number of points for p
     negf%n_kt = 10          ! Numero di kT per l'integrazione
     negf%n_poles = 3        ! Numero di poli
     negf%activecont = 0     ! contact selfenergy
     allocate(negf%ni(1))
     allocate(negf%nf(1))
     negf%ni(1) = 1
     negf%nf(1) = 2          !
     negf%min_or_max = 1     ! Set reference cont to max(mu)
     negf%refcont = 1        ! call set_ref_cont()
     negf%outer_blocks = 2   ! Compute full D.M. L,U extra
     negf%dumpHS = .false.
     negf%int_acc = 1.d-3    ! Integration accuracy
                             ! Only in adaptive refinement
     negf%scba_elastic_tol = 1.d-7
     negf%scba_inelastic_tol = 1.d-7
     negf%ndos_proj = 0
     negf%particle = 1       ! Used for setting correct fermi function in real_axis_int.
                             ! 1: electrons
                             !-1: holes
                             ! Can be -1 only in compute_density_efa

     negf%bulk_cont_density = .false. ! Flag to calculate the density of the bulk contacts

     negf%surface_green_cache = TMatrixCacheDisk(scratch_path=negf%scratch_path)

     ! Initialize the cache space for Gr and Gn
     call init_cache_space(negf, 'G_r')
     call init_cache_space(negf, 'G_n')

   end subroutine set_defaults



   subroutine print_all_vars(negf,io)
     type(Tnegf) :: negf
     integer, intent(in) :: io

     integer :: ii

     write(io,*) 'verbose=',negf%verbose

     write(io,*) 'scratch= "'//trim(negf%scratch_path)//'"'
     write(io,*) 'output= "'//trim(negf%out_path)//'"'

     write(io,*) 'mu= ', negf%mu
     write(io,*) 'mu_n= ', negf%mu_n
     write(io,*) 'mu_p= ', negf%mu_p
     write(io,*) 'T= ',negf%kbT

     do ii = 1, negf%str%num_conts
       write(io,*) 'Contact Parameters for ',trim(negf%cont(ii)%name)
       write(io,*) 'mu= ', negf%cont(ii)%mu
       write(io,*) 'WideBand= ', negf%cont(ii)%FictCont
       write(io,*) 'DOS= ', negf%cont(ii)%contact_DOS
       write(io,*) 'mu_n= ', negf%cont(ii)%mu_n
       write(io,*) 'mu_p= ', negf%cont(ii)%mu_p
       write(io,*) 'Density Matrix T= ',negf%cont(ii)%kbT_dm
       write(io,*) 'Transmission T= ',negf%cont(ii)%kbT_t
     end do

     write(io,*) 'Contour Parameters:'
     write(io,*) 'Ec= ', negf%Ec
     write(io,*) 'Ev= ', negf%Ev
     write(io,*) 'DEc= ', negf%DeltaEc
     write(io,*) 'DEv= ', negf%DeltaEv
     write(io,*) 'ReadoldDM_SGFs= ',negf%ReadoldDM_SGFs
     write(io,*) 'ReadoldT_SGF= ',negf%ReadoldT_SGFs
     write(io,*) 'Np_n= ',negf%Np_n
     write(io,*) 'Np_p= ', negf%Np_p
     write(io,*) 'nkT= ', negf%n_kt
     write(io,*) 'nPoles= ', negf%n_poles
     write(io,*) 'min_or_max= ', negf%min_or_max
     write(io,*) 'refcont= ', negf%refcont

     write(io,*) 'Transmission Parameters:'
     write(io,*) 'Emin= ',negf%Emin
     write(io,*) 'Emax= ',negf%Emax
     write(io,*) 'Estep= ',negf%Estep
     write(io,*) 'g_spin= ',negf%g_spin
     write(io,*) 'delta= ',negf%delta
     write(io,*) 'dos_delta= ',negf%dos_delta
     write(io,*) 'ni= ', negf%ni
     write(io,*) 'nf= ', negf%nf
     write(io,*) 'DorE= ',negf%DorE

     write(io,*) 'Internal variables:'
     write(io,*) 'kp= ', negf%ikpoint
     write(io,*) 'wght= ', negf%kwght
     write(io,*) 'E= ',negf%E
     select case(negf%outer_blocks)
     case(0)
        write(io,*) 'outer_blocks= No'
     case(1)
        write(io,*) 'outer_blocks= upper diagonal'
     case(2)
        write(io,*) 'outer_blocks= upper/lower diagonal'
     end select
     write(io,*) 'DOS= ',negf%dos
     write(io,*) 'Eneconv= ',negf%eneconv
     write(io,*) 'activecont= ', negf%activecont
     !write(io,*) 'internalDM= ', negf%internalDM
  end subroutine print_all_vars

end module lib_param
