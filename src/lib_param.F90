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


module lib_param

  use ln_precision, only : dp
  use globals
  use mat_def
  use ln_structure, only : TStruct_info, print_Tstruct
  use input_output
  use elph, only : init_elph_1, Telph, destroy_elph, init_elph_2, init_elph_3
  use phph
  use energy_mesh, only : mesh
  use interactions, only : Interaction
  use elphdd, only : ElPhonDephD, ElPhonDephD_create 
  use elphdb, only : ElPhonDephB, ElPhonDephB_create
  use elphds, only : ElPhonDephS, ElPhonDephS_create
#:if defined("MPI")  
  use libmpifx_module, only : mpifx_comm
#:endif
  implicit none
  private

  public :: Tnegf, intArray, TEnGrid
  public :: set_defaults, print_all_vars

  public :: set_bp_dephasing
  public :: set_elph_dephasing, destroy_elph_model
  public :: set_elph_block_dephasing, set_elph_s_dephasing
  public :: set_phph
  integer, public, parameter :: MAXNCONT=10

  type intArray
    integer, dimension(:), allocatable :: indexes
  end type intArray 


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
     integer :: cpu
     complex(dp) :: Ec
     complex(dp) :: wght
  end type TEnGrid

  !DAR begin - in Tnegf
  !-----------------------------------------------------------------------------

  type :: Tdeph_bp
     real(dp), allocatable, dimension(:) :: coupling
  end type Tdeph_bp

  type :: Tdephasing
     type(Tdeph_bp) :: bp
  end type Tdephasing

  type Tcontact
     character(132) :: name
     logical :: tWriteSelfEnergy = .false.
     logical :: tReadSelfEnergy = .false.
     complex(kind=dp), dimension(:,:,:), allocatable :: SelfEnergy ! Electrode Self-Energy
     logical :: tWriteSurfaceGF = .false.
     logical :: tReadSurfaceGF = .false.
     complex(kind=dp), dimension(:,:,:), allocatable :: SurfaceGF ! Electrode Surface Green Function
  end type Tcontact

  type Telectrons
     integer :: IndexEnergy
  end type Telectrons

  type Toutput
     logical :: tWriteDOS = .false.
     logical :: tDOSwithS = .false.
  end type Toutput

  type Ttranas
     type(Telectrons) :: e
     type(Tcontact), dimension(:), allocatable :: cont
     type(Toutput) :: out 
  end type Ttranas

  !-----------------------------------------------------------------------------
  !DAR - end

 !> General libnegf container
 !! Contains input data, runtime quantities and output data
 type Tnegf
   !! Input parameters: set by library user
   !! General
   integer :: verbose
#:if defined("MPI")  
   type(mpifx_comm) :: mpicomm
#:endif   
   integer  :: ReadoldSGF            ! 0: Read 1: compute 2: comp & save
   character(len=LST) :: scratch_path    ! Folder for scratch work
   character(len=LST) :: out_path        ! Folder for output data
   real(dp) :: g_spin            ! spin degeneracy
   real(dp) :: delta             ! delta for G.F. 
   real(dp) :: dos_delta         ! additional delta to force more broadening in the DOS 
   real(dp) :: eneconv           ! Energy conversion factor
   integer :: iteration          ! Number of current SCC itaration
   integer  :: spin              ! spin component
   real(dp) :: wght              ! k-point weight 
   integer :: kpoint             ! k-point index
   character(1) :: DorE              ! Density or En.Density

   !! Contacts info
   real(dp) :: mu_n(MAXNCONT)    ! electrochemical potential (el)
   real(dp) :: mu_p(MAXNCONT)    ! electrochemical potential (hl)
   real(dp) :: mu(MAXNCONT)          ! Electrochemical Potential (dft calculation)   
   real(dp) :: contact_DOS(MAXNCONT) ! Ficticious contact DOS
   logical  :: FictCont(MAXNCONT)    ! Ficticious contact 
   real(dp) :: kbT_dm(MAXNCONT) ! Electronic temperature
   real(dp) :: kbT_t(MAXNCONT)    ! Electronic temperature

   !! Contour integral
   integer :: Np_n(2)            ! Number of points for n 
   integer :: Np_p(2)            ! Number of points for p 
   integer :: Np_real(11)        ! Number of points for integration over real axis
   integer :: n_kt               ! Number of kT extending integrations
   integer :: n_poles            ! Number of poles 
   real(dp) :: Ec                ! conduction band edge 
   real(dp) :: Ev                ! valence band edge

   !! Real axis
   real(dp) :: Emin              ! Tunneling or dos interval
   real(dp) :: Emax              ! 
   real(dp) :: Estep             ! Tunneling or dos E step

   !! Emitter and collector for transmission or Meir-Wingreen 
   !! (only emitter in this case)
   integer :: ni(MAXNCONT)       ! ni: emitter contact list 
   integer :: nf(MAXNCONT)       ! nf: collector contact list


   integer  :: nldos                 ! Number of LDOS intervals
   type(intArray), dimension(:), allocatable :: LDOS !Array of LDOS descriptor 
                                                     !(contain only index of atoms 
                                                     !for LDOS projection)  
   real(dp) :: DeltaEc           ! safe guard energy below Ec
   real(dp) :: DeltaEv           ! safe guard energy above Ev

   !! Runtime variables: used internally by the library
   type(format) :: form              ! Form of file-Hamiltonian
   logical  :: dumpHS                ! Used for debug
   real(dp) :: muref             ! reference elec.chem potential
   real(dp) :: E                 ! Holding variable 
   real(dp) :: dos               ! Holding variable
   integer :: activecont         ! contact selfenergy
   integer :: min_or_max         ! in input: 0 take minimum, 1 take maximum mu  
   integer :: refcont            ! reference contact (for non equilib)
   integer :: outer              ! flag switching computation of     
                                 ! the Device/Contact DM
                                 ! 0 none; 1 upper block; 2 all


   !! Note: H,S are partitioned immediately after input, therefore they are 
   !! built runtime from input variable
   type(z_CSR), pointer :: H => null()    ! Points to externally allocated H
   type(z_CSR), pointer :: S => null()
   type(z_DNS) :: HC(MAXNCONT)
   type(z_DNS) :: SC(MAXNCONT)
   type(z_DNS) :: HMC(MAXNCONT)
   type(z_DNS) :: SMC(MAXNCONT)
   type(z_CSR), pointer :: rho => null()      ! Holding output Matrix
   type(z_CSR), pointer :: rho_eps => null()  ! Holding output Matrix
   logical    :: isSid           ! True if overlap S == Id
   logical    :: intHS           ! tells HS are internally allocated
   logical    :: intDM           ! tells DM is internally allocated

   type(TStruct_Info) :: str     ! system structure
   integer :: iE                 ! Energy point (integer point)
   complex(dp) :: Epnt           ! Energy point (complex)
   integer :: local_en_points    ! Local number of energy points
   type(TEnGrid), dimension(:), allocatable :: en_grid
   real(dp) :: int_acc           ! integration accuracy
   real(dp), dimension(:), pointer :: E_singular => null()
   real(dp) :: delta_singular
   type(Telph) :: elph           ! electron-phonon data
   type(Tphph) :: phph           ! phonon-phonon data

   type(mesh) :: emesh           ! energy mesh for adaptive Simpson

   ! Many Body Interactions
   class(Interaction), allocatable :: inter

   !! Output variables: these are filled by internal subroutines to stor
   !! library output
   real(dp), dimension(:,:), pointer :: tunn_mat => null()
   !$real(dp), dimension(:,:), pointer :: tunn_mat_bp => null()               !DAR
   real(dp), dimension(:,:), pointer :: ldos_mat => null()
   real(dp), dimension(:), pointer :: currents => null() ! value of contact currents 

   !DAR begin - in Tnegf
   !----------------------------------------------------------------------------
   logical :: tNoGeometry = .false.
   logical :: tOrthonormal = .false.
   logical :: tOrthonormalDevice = .false.
   logical :: tTrans = .false.
   logical :: tCalcSelfEnergies = .true.
   integer :: NumStates       
   character(len=LST) :: FileName
   logical :: tManyBody = .false.
   logical :: tElastic = .true.
   logical :: tDephasingVE = .false.
   logical :: tDephasingBP = .false.
   logical :: tZeroCurrent = .false.
   integer :: MaxIter = 1000
   logical :: tWrite_ldos = .false.
   logical :: tWrite_negf_params = .false.
   type(Ttranas) :: tranas
   type(Tdephasing) :: deph
   !----------------------------------------------------------------------------
   !DAR end

 end type Tnegf

contains
 
  !> Set buttiker probe dephasing
  subroutine set_bp_dephasing(negf, coupling)
    type(Tnegf) :: negf
    real(dp),  dimension(:), intent(in) :: coupling

    if (.not.allocated(negf%deph%bp%coupling)) then
       allocate(negf%deph%bp%coupling(size(coupling)))
    end if
    negf%deph%bp%coupling = coupling

  end subroutine set_bp_dephasing
  
  !> Set values for the local electron phonon dephasing model
  !! (elastic scattering only)
  subroutine set_elph_dephasing(negf, coupling, niter)

    type(Tnegf) :: negf
    type(ElPhonDephD) :: elphdd_tmp
    real(dp),  dimension(:), allocatable, intent(in) :: coupling
    integer :: niter
    
    call elphondephd_create(elphdd_tmp, negf%str, coupling, niter, 1.0d-7)
    if(.not.allocated(negf%inter)) allocate(negf%inter, source=elphdd_tmp)

  end subroutine set_elph_dephasing

  !> Set values for the semi-local electron phonon dephasing model
  !! (elastic scattering only)
  subroutine set_elph_block_dephasing(negf, coupling, orbsperatom, niter)

    type(Tnegf) :: negf
    type(ElPhonDephB) :: elphdb_tmp
    real(dp),  dimension(:), allocatable, intent(in) :: coupling
    integer,  dimension(:), allocatable, intent(in) :: orbsperatom
    integer :: niter
    
    call elphondephb_create(elphdb_tmp, negf%str, coupling, orbsperatom, niter, 1.0d-7)
    if(.not.allocated(negf%inter)) allocate(negf%inter, source=elphdb_tmp)

  end subroutine set_elph_block_dephasing

 !> Set values for the semi-local electron phonon dephasing model
  !! (elastic scattering only)
  subroutine set_elph_s_dephasing(negf, coupling, orbsperatom, niter)

    type(Tnegf) :: negf
    type(ElPhonDephS) :: elphds_tmp
    real(dp),  dimension(:), allocatable, intent(in) :: coupling
    integer,  dimension(:), allocatable, intent(in) :: orbsperatom
    integer :: niter

    call elphondephs_create(elphds_tmp, negf%str, coupling, orbsperatom, negf%S, niter, 1.0d-7)
    if(.not.allocated(negf%inter)) allocate(negf%inter, source=elphds_tmp)

  end subroutine set_elph_s_dephasing



  !> Destroy elph model. This routine is accessible from interface as 
  !! it can be meaningful to "switch off" elph when doing different
  !! task (density or current)
  subroutine destroy_elph_model(negf)
    type(Tnegf) :: negf
    
    if (allocated(negf%inter)) deallocate(negf%inter)
    call destroy_elph(negf%elph)

  end subroutine destroy_elph_model

  subroutine set_defaults(negf)
    type(Tnegf) :: negf    

     negf%verbose = 10

     negf%scratch_path = './GS/'
     negf%out_path = './'
     negf%DorE = 'D'           ! Density or En.Density

     negf%ReadoldSGF = 1       ! Compute Surface G.F. do not save
     negf%FictCont = .false.   ! Ficticious contact 

     negf%mu = 0.d0            ! Potenziale elettrochimico
     negf%contact_DOS = 0.d0   ! Ficticious contact DOS

     negf%wght = 1.d0
     negf%kpoint = 1

     negf%mu_n = 0.d0
     negf%mu_p = 0.d0
     negf%Ec = 0.d0
     negf%Ev = 0.d0
     negf%DeltaEc = 0.d0
     negf%DeltaEv = 0.d0

     negf%E = 0.d0            ! Holding variable 
     negf%dos = 0.d0          ! Holding variable
     negf%eneconv = 1.d0      ! Energy conversion factor

     negf%isSid = .false.         
     negf%intHS = .true.
     negf%intDM = .true.

     negf%delta = 1.d-4      ! delta for G.F. 
     negf%dos_delta = 1.d-4  ! delta for DOS 
     negf%Emin = 0.d0        ! Tunneling or dos interval
     negf%Emax = 0.d0        ! 
     negf%Estep = 0.d0       ! Tunneling or dos E step
     negf%kbT_dm = 0.d0         ! electronic temperature
     negf%kbT_t = 0.d0         ! electronic temperature
     negf%g_spin = 2.d0      ! spin degeneracy

     negf%Np_n = (/20, 20/)  ! Number of points for n 
     negf%Np_p = (/20, 20/)  ! Number of points for p 
     negf%n_kt = 10          ! Numero di kT per l'integrazione
     negf%n_poles = 3        ! Numero di poli 
     negf%iteration = 1      ! Iterazione (SCC)
     negf%activecont = 0     ! contact selfenergy
     negf%ni = 0             ! ni
     negf%ni(1) = 1
     negf%nf = 0             ! nf
     negf%nf(1) = 2          !
     negf%min_or_max = 1     ! Set reference cont to max(mu)  
     negf%refcont = 1        ! call set_ref_cont()
     negf%outer = 2          ! Compute full D.M. L,U extra
     negf%dumpHS = .false.
     negf%int_acc = 1.d-3    ! Integration accuracy 
                             ! Only in adaptive refinement 
     negf%nldos = 0

   end subroutine set_defaults


   subroutine set_phph(negf, order, filename)
      type(Tnegf) :: negf
      integer, intent(in) :: order
      character(*), intent(in) :: filename
 
      print*,'(set_phph) init_phph' 
      call init_phph(negf%phph, negf%str%central_dim, order, negf%str%mat_PL_start, negf%str%mat_PL_end) 
      
      print*,'(set_phph) load_phph_coupl' 
      call load_phph_couplings(negf%phph, filename) 

   end subroutine set_phph



   subroutine print_all_vars(negf,io)
     type(Tnegf) :: negf    
     integer :: io

     call print_Tstruct(negf%str)

     write(io,*) 'verbose=',negf%verbose

     write(io,*) 'scratch= "'//trim(negf%scratch_path)//'"'
     write(io,*) 'output= "'//trim(negf%out_path)//'"'

     write(io,*) 'isSid= ',negf%isSid

     write(io,*) 'Contact Parameters:'
     write(io,*) 'mu= ', negf%mu
     write(io,*) 'WideBand= ', negf%FictCont
     write(io,*) 'DOS= ', negf%contact_DOS

     write(io,*) 'Contour Parameters:'
     write(io,*) 'mu_n= ', negf%mu_n
     write(io,*) 'mu_p= ', negf%mu_p
     write(io,*) 'Ec= ', negf%Ec 
     write(io,*) 'Ev= ', negf%Ev
     write(io,*) 'DEc= ', negf%DeltaEc 
     write(io,*) 'DEv= ', negf%DeltaEv
     write(io,*) 'ReadoldSGF= ',negf%ReadoldSGF
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
     write(io,*) 'Density Matrix T= ',negf%kbT_dm
     write(io,*) 'Transmission T= ',negf%kbT_t
     write(io,*) 'g_spin= ',negf%g_spin 
     write(io,*) 'delta= ',negf%delta
     write(io,*) 'dos_delta= ',negf%dos_delta
     write(io,*) 'ni= ', negf%ni
     write(io,*) 'nf= ', negf%nf
     write(io,*) 'DorE= ',negf%DorE

     write(io,*) 'Internal variables:'
     write(io,*) 'intHS= ',negf%intHS
     write(io,*) 'intDM= ',negf%intDM
     write(io,*) 'kp= ', negf%kpoint
     write(io,*) 'wght= ', negf%wght
     write(io,*) 'E= ',negf%E
     write(io,*) 'outer= ', negf%outer
     write(io,*) 'DOS= ',negf%dos
     write(io,*) 'Eneconv= ',negf%eneconv
     write(io,*) 'iter= ', negf%iteration
     write(io,*) 'activecont= ', negf%activecont 


  end subroutine print_all_vars

end module lib_param
 




