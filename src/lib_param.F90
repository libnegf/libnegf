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
  use ln_structure, only : TStruct_info
  use input_output
  use elph

  implicit none
  private

  public :: Tnegf
  public :: fill_parameters, pass_HS
  public :: set_convfactor, set_fermi, set_potentials, set_fictcont
  public :: set_readoldsgf, set_computation, set_iteration, set_defaults
  integer, public, parameter :: MAXNCONT=10

  type Tnegf
     
     integer :: verbose

     character(LST) :: file_re_H
     character(LST) :: file_im_H
     character(LST) :: file_re_S
     character(LST) :: file_im_S
     character(LST) :: file_struct
     character(LST) :: scratch_path
     character(LST) :: out_path
     character(1) :: DorE              ! Density or En.Density
     type(format) :: form              ! Form of file-Hamiltonian

     integer   :: ReadoldSGF           ! 0: Read 1: compute 2: comp & save
     logical   :: FictCont(MAXNCONT)   ! Ficticious contact 
     logical   :: dumpHS               ! Used for debug
     logical   :: writeLDOS            ! Write the LDOS for every k-point 
     logical   :: writeTunn            ! Write T(E) for every k-point

     real(dp) :: mu(MAXNCONT)          ! Potenziale elettrico
     real(dp) :: Efermi(MAXNCONT)      ! Energia di Fermi dei contatti
     real(dp) :: contact_DOS(MAXNCONT) ! Ficticious contact DOS

     integer  :: nLdos                 ! Number of LDOS intervals
     integer,Dimension(:,:), pointer :: LDOS => null()    ! LDOS intervals

     real(dp) :: mu_n
     real(dp) :: mu_p
     real(dp) :: muref             ! reference elec.chem potential
     real(dp) :: Ec
     real(dp) :: Ev
     real(dp) :: DeltaEc
     real(dp) :: DeltaEv

     real(dp) :: E                 ! Holding variable 
     real(dp) :: dos               ! Holding variable
     real(dp) :: eneconv           ! Energy conversion factor

     type(z_CSR) :: H
     type(z_CSR) :: S
     type(z_CSR) :: HM
     type(z_CSR) :: SM     
     type(z_DNS) :: HC(MAXNCONT)
     type(z_DNS) :: SC(MAXNCONT)
     type(z_DNS) :: HMC(MAXNCONT)
     type(z_DNS) :: SMC(MAXNCONT)
     type(z_CSR) :: Gr             ! Holding output Matrix
     type(z_CSR) :: rho            ! Holding output Matrix
     type(z_CSR) :: rho_eps        ! Holding output Matrix
     logical    :: isSid          

     type(TStruct_Info) :: str     ! system structure

     real(dp) :: delta       ! delta for G.F. 
     real(dp) :: Emin        ! Tunneling or dos interval
     real(dp) :: Emax        ! 
     real(dp) :: Estep       ! Tunneling or dos E step
     real(dp) :: kbT         ! electronic temperature
     real(dp) :: spin        ! spin degeneracy

     real(dp) :: wght        ! kp weight 
     integer :: kpoint       ! kp index
     integer :: Epnt         ! Energy point
     real(dp), dimension(:,:), pointer :: tunn_mat => null()
     real(dp), dimension(:,:), pointer :: ldos_mat => null()
     real(dp), dimension(:), pointer :: currents => null() ! value of contact currents 

     integer :: Np_n(2)      ! Number of points for n 
     integer :: Np_p(2)      ! Number of points for p 
     integer :: Np_real      ! Number of points for integration over real axis
     integer :: n_kt         ! Numero di kT per l'integrazione
     integer :: n_poles      ! Numero di poli 
     integer :: iteration    ! Iterazione (SCC)
     integer :: activecont   ! contact selfenergy
     integer :: ni(MAXNCONT) ! ni
     integer :: nf(MAXNCONT) ! nf
     integer :: refcont      ! reference contact (for non equilib)
                             ! in input: 0 take minimum, 1 take maximum mu  
     integer :: outer        ! flag switching computation of  
                             ! the Device/Contact DM
                             ! 0 none; 1 upper block; 2 all

	 real(dp), dimension(:), pointer :: E_singular => null()
	 real(dp) :: delta_singular

	 type(Telph) :: elph     ! electron-phonon data

  end type Tnegf

contains
  
  subroutine set_convfactor(negf, eneconv)
    type(Tnegf) :: negf
    real(dp) :: eneconv
    
    negf%eneconv=eneconv

  end subroutine set_convfactor
 ! -------------------------------------------------------------------
  
  subroutine set_fictcont(negf,cont,dos)
    type(Tnegf) :: negf
    integer :: cont
    real(dp) :: DOS
 
    negf%FictCont(cont) = .true. 
    negf%contact_DOS(cont) = DOS

  end subroutine set_fictcont
 ! -------------------------------------------------------------------

  subroutine set_iteration(negf,iter)
    type(Tnegf) :: negf
    integer :: iter

    negf%iteration = iter
  end subroutine set_iteration      
 ! -------------------------------------------------------------------

  subroutine set_computation(negf,DorE) 
    type(Tnegf) :: negf
    character(1) :: DorE           !Density or En.Density

    negf%DorE=DorE
  end subroutine set_computation      
 ! -------------------------------------------------------------------

  subroutine set_readOldSGF(negf,flag) 
    type(Tnegf) :: negf
    integer :: flag

    negf%ReadoldSGF=flag
  end subroutine set_readoldsgf    
 ! -------------------------------------------------------------------

  subroutine set_fermi(negf,ncont,efermi)
    type(Tnegf) :: negf
    integer :: ncont
    real(dp) :: efermi(*)
    
    negf%Efermi(1:ncont) = Efermi(1:ncont)
  end subroutine set_fermi

  subroutine set_potentials(negf,ncont,mu)
    type(Tnegf) :: negf
    integer :: ncont
    real(dp) :: mu(*)

    negf%mu(1:ncont) = mu(1:ncont)
  end subroutine set_potentials
 ! -------------------------------------------------------------------


  subroutine fill_parameters(negf, verbose, mu_n, mu_p, Ec, Ev, &
        DeltaEc, DeltaEv, delta, Emin, Emax, Estep, &
        kbT, Np_n, Np_p, n_kt, n_poles, spin)

    type(Tnegf) :: negf
    integer :: verbose

    real(dp) :: mu_n
    real(dp) :: mu_p
    real(dp) :: Ec
    real(dp) :: Ev
    real(dp) :: DeltaEc
    real(dp) :: DeltaEv     

    real(dp) :: delta
    real(dp) :: Emin
    real(dp) :: Emax
    real(dp) :: Estep
    real(dp) :: kbT

    integer :: Np_n(2)
    integer :: Np_p(2)   
    integer :: n_kt
    integer :: n_poles
    real(dp) :: spin   

    negf%verbose = verbose

    negf%mu_n = mu_n
    negf%mu_p = mu_p
    negf%Ec = Ec 
    negf%Ev = Ev
    negf%DeltaEc = DeltaEc 
    negf%DeltaEv = DeltaEv 

    negf%delta = delta
    negf%Emin = Emin
    negf%Emax = Emax
    negf%Estep = Estep
    negf%kbT = kbT

    negf%Np_n = Np_n
    negf%Np_p = Np_p
    negf%n_kt = n_kt
    negf%n_poles = n_poles
    negf%spin = spin

  end subroutine fill_parameters

  subroutine pass_HS(negf,H,S)
    type(Tnegf) :: negf    
    type(z_CSR) :: H
    type(z_CSR), optional :: S

    call create(negf%H,H%nrow,H%ncol,H%nnz)
    negf%H%nzval = H%nzval
    negf%H%colind = H%colind
    negf%H%rowpnt = H%rowpnt
          

    if (present(S)) then
       negf%isSid=.false.
       call create(negf%S,S%nrow,S%ncol,S%nnz)
       negf%S%nzval = S%nzval
       negf%S%colind = S%colind
       negf%S%rowpnt = S%rowpnt
    else
       negf%isSid=.true.
       call create_id(negf%S,negf%H%nrow) 
    endif

  end subroutine pass_HS


  subroutine set_defaults(negf)
    type(Tnegf) :: negf    

     negf%verbose = 10

     negf%file_re_H = ''
     negf%file_im_H = ''
     negf%file_re_S = ''
     negf%file_im_S = ''
     negf%file_struct = ''
     negf%scratch_path = './GS/'
     negf%out_path = './'
     negf%DorE = 'D'           ! Density or En.Density

     negf%ReadoldSGF = 1       ! Compute Surface G.F. do not save
     negf%FictCont = .false.   ! Ficticious contact 

     negf%mu = 0.d0            ! Potenziale elettrico
     negf%efermi= 0.d0         ! Energia di Fermi dei contatti
     negf%contact_DOS = 0.d0   ! Ficticious contact DOS

    ! negf%nLdos = 0.d0                ! Number of LDOS intervals
    ! negf%LDOS(:,:) = 0.d0    ! LDOS intervals
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
     negf%writeLDOS = .false. 

     negf%delta = 1.d-4      ! delta for G.F. 
     negf%Emin = 0.d0        ! Tunneling or dos interval
     negf%Emax = 0.d0        ! 
     negf%Estep = 0.d0       ! Tunneling or dos E step
     negf%kbT = 0.d0        ! electronic temperature
     negf%spin = 2.d0        ! spin degeneracy

     negf%Np_n = (/20, 20/)   ! Number of points for n 
     negf%Np_p = (/20, 20/)   ! Number of points for p 
     negf%n_kt = 10          ! Numero di kT per l'integrazione
     negf%n_poles = 3        ! Numero di poli 
     negf%iteration = 1      ! Iterazione (SCC)
     negf%activecont = 0     ! contact selfenergy
     negf%ni = 0             ! ni
     negf%ni(1) = 1
     negf%nf = 0             ! nf
     negf%nf(1) = 2
     negf%refcont = 1        ! reference contact (?)
     negf%outer = 2          ! Compute full D.M.
     
   end subroutine set_defaults

end module lib_param
 




