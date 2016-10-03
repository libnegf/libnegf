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
  use energy_mesh, only : mesh
  use interactions, only : Interaction
  use elphdd, only : ElPhonDephD, ElPhonDephD_create 
  use elphdb, only : ElPhonDephB, ElPhonDephB_create
  use elphds, only : ElPhonDephS, ElPhonDephS_create

  implicit none
  private

  public :: Tnegf, intarray, TEnGrid
  public :: fill_parameters, pass_HS, pass_DM
  public :: set_convfactor, set_fermi, set_potentials, set_fictcont
  public :: set_readoldsgf, set_computation, set_iteration, set_defaults
  public :: print_all_vars, set_elph_dephasing, destroy_elph_model
  public :: set_elph_block_dephasing, set_elph_s_dephasing
  integer, public, parameter :: MAXNCONT=10

  type intarray
    integer, dimension(:), allocatable :: indexes
  end type intarray 




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


 !> General libnegf container
 !! Contains input data, runtime quantities and output data
 type Tnegf
   !! Input parameters: set by library user
   !! General
   integer :: verbose
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
   real(dp) :: kbT(MAXNCONT)         ! Electronic temperature

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

 

   type(intarray), dimension(:), allocatable :: LDOS !Array of LDOS descriptor 
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
   integer :: minmax             ! in input: 0 take minimum, 1 take maximum mu  
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
   type(TEnGrid), dimension(:), allocatable :: en_grid
   real(dp) :: int_acc           ! integration accuracy
   real(dp), dimension(:), pointer :: E_singular => null()
   real(dp) :: delta_singular
   type(Telph) :: elph           ! electron-phonon data
   type(mesh) :: emesh           ! energy mesh for adaptive Simpson

   ! Many Body Interactions
   class(Interaction), allocatable :: inter

   !! Output variables: these are filled by internal subroutines to stor
   !! library output
   real(dp), dimension(:,:), pointer :: tunn_mat => null()
   real(dp), dimension(:,:), pointer :: ldos_mat => null()
   real(dp), dimension(:), pointer :: currents => null() ! value of contact currents 
   

 end type Tnegf

contains
  
 ! -------------------------------------------------------------------
  ! -----------------------------------------------------
  !  Pass an externally allocated density matrix
  ! -----------------------------------------------------
  subroutine pass_DM(negf,rho, rhoE)
    type(Tnegf) :: negf    
    type(z_CSR), optional, target :: rho
    type(z_CSR), optional, target :: rhoE
 
    if (present(rho)) then
       negf%rho => rho
       if(allocated(negf%rho%nzval)) then
          call destroy(negf%rho)
       endif   
    endif 

    if (present(rhoE)) then
       negf%rho_eps => rhoE
       if (allocated(negf%rho_eps%nzval)) then
          call destroy(negf%rho_eps)
       endif 
    end if
    
    negf%intDM = .false.

  end subroutine pass_DM
 
  ! -----------------------------------------------------
  !  Allocate and copy H,S 
  ! -----------------------------------------------------
  subroutine copy_HS(negf,H,S)
    type(Tnegf) :: negf    
    type(z_CSR), target :: H
    type(z_CSR), optional, target :: S

    call create(negf%H,H%nrow,H%ncol,H%nnz)
    negf%H%nzval = H%nzval
    negf%H%colind = H%colind
    negf%H%rowpnt = H%rowpnt
    negf%H%sorted = H%sorted   
    
    if (present(S)) then
       negf%isSid=.false.
       call create(negf%S,S%nrow,S%ncol,S%nnz)
       negf%S%nzval = S%nzval
       negf%S%colind = S%colind
       negf%S%rowpnt = S%rowpnt
       negf%S%sorted = S%sorted
    else
       negf%isSid=.true.
       call create_id(negf%S,negf%H%nrow) 
    endif
    
    negf%intHS = .true.

  end subroutine copy_HS
  
  !> Set values for the local electron phonon dephasing model
  !! (elastic scattering only)
  subroutine set_elph_dephasing(negf, coupling, niter)
    type(Tnegf) :: negf
    type(ElPhonDephD) :: elphdd_tmp
    real(dp),  dimension(:), allocatable, intent(in) :: coupling
    integer :: niter
    
    !call init_elph_1(negf%elph, coupling, niter)
    call elphondephd_create(elphdd_tmp, negf%str, coupling, niter, 1.0d-7)
    allocate(negf%inter, source=elphdd_tmp)

  end subroutine set_elph_dephasing

  !> Set values for the semi-local electron phonon dephasing model
  !! (elastic scattering only)
  subroutine set_elph_block_dephasing(negf, coupling, orbsperatom, niter)
    type(Tnegf) :: negf
    type(ElPhonDephB) :: elphdb_tmp
    real(dp),  dimension(:), allocatable, intent(in) :: coupling
    integer,  dimension(:), allocatable, intent(in) :: orbsperatom
    integer :: niter
    
    !! Verify that the size of the coupling fits with the Hamiltonian
    !! of the device region
    !if (size(coupling).ne.negf%H%nrow) then
    !  write(*,*) 'Elph dephasing model coupling size does not match '
    !endif
    !call init_elph_2(negf%elph, coupling, orbsperatom, niter, &
    !    negf%str%mat_PL_start)
    call elphondephb_create(elphdb_tmp, negf%str, coupling, orbsperatom, niter, 1.0d-7)
    allocate(negf%inter, source=elphdb_tmp)

  end subroutine set_elph_block_dephasing

 !> Set values for the semi-local electron phonon dephasing model
  !! (elastic scattering only)
  subroutine set_elph_s_dephasing(negf, coupling, orbsperatom, niter)
    type(Tnegf) :: negf
    type(ElPhonDephS) :: elphds_tmp
    real(dp),  dimension(:), allocatable, intent(in) :: coupling
    integer,  dimension(:), allocatable, intent(in) :: orbsperatom
    integer :: niter
    
    !! Verify that the size of the coupling fits with the Hamiltonian
    !! of the device region
    !if (size(coupling).ne.negf%H%nrow) then
    !  write(*,*) 'Elph dephasing model coupling size does not match '
    !endif
    !call init_elph_3(negf%elph, coupling, orbsperatom, niter, &
    !    negf%str%mat_PL_start, negf%S)
    call elphondephs_create(elphds_tmp, negf%str, coupling, orbsperatom, negf%S, niter, 1.0d-7)
    allocate(negf%inter, source=elphds_tmp)

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
     negf%dos_delta = 1.d-4      ! delta for DOS 
     negf%Emin = 0.d0        ! Tunneling or dos interval
     negf%Emax = 0.d0        ! 
     negf%Estep = 0.d0       ! Tunneling or dos E step
     negf%kbT = 0.d0        ! electronic temperature
     negf%g_spin = 2.d0        ! spin degeneracy

     negf%Np_n = (/20, 20/)   ! Number of points for n 
     negf%Np_p = (/20, 20/)   ! Number of points for p 
     negf%n_kt = 10          ! Numero di kT per l'integrazione
     negf%n_poles = 3        ! Numero di poli 
     negf%iteration = 1      ! Iterazione (SCC)
     negf%activecont = 0     ! contact selfenergy
     negf%ni = 0             ! ni
     negf%ni(1) = 1
     negf%nf = 0             ! nf
     negf%nf(1) = 2          !
     negf%minmax = 1         ! Set reference cont to max(mu)  
     negf%refcont = 1        ! call set_ref_cont()
     negf%outer = 2          ! Compute full D.M. L,U extra
     negf%dumpHS = .false.
     negf%int_acc = 1.d-3    ! Integration accuracy 
                             ! Only in adaptive refinement 
     negf%nldos = 0

   end subroutine set_defaults

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
     write(io,*) 'minmax= ', negf%minmax
     write(io,*) 'refcont= ', negf%refcont

     write(io,*) 'Transmission Parameters:'
     write(io,*) 'Emin= ',negf%Emin
     write(io,*) 'Emax= ',negf%Emax
     write(io,*) 'Estep= ',negf%Estep
     write(io,*) 'T= ',negf%kbT
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
 




