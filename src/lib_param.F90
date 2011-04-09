module lib_param

  use ln_precision, only : dp
  use globals
  use mat_def
  use ln_structure, only : TStruct_info
  use input_output

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
     character(1) :: DorE              ! Density or En.Density
     type(format) :: form              ! Form of file-Hamiltonian

     logical   :: ReadoldSGF           ! Read computed Surface G.F.
     logical   :: FictCont(MAXNCONT)   ! Ficticious contact 

     real(dp) :: mu(MAXNCONT)          ! Potenziale elettrico
     real(dp) :: Efermi(MAXNCONT)      ! Energia di Fermi dei contatti
     real(dp) :: contact_DOS(MAXNCONT) ! Ficticious contact DOS

     integer  :: nLdos                 ! Number of LDOS intervals
     integer,Dimension(:,:), pointer :: LDOS(:,:) => null()    ! LDOS intervals

     real(dp) :: mu_n
     real(dp) :: mu_p
     real(dp) :: Ec
     real(dp) :: Ev
     real(dp) :: DeltaEc
     real(dp) :: DeltaEv

     real(dp) :: E                 ! Holding variable 
     real(dp) :: dos               ! Holding variable
     real(dp) :: eneconv           ! Energy conversion factor

     type(z_CSR), pointer :: H
     type(z_CSR), pointer :: S
     type(z_CSR) :: HM
     type(z_CSR) :: SM     
     type(z_DNS) :: HC(MAXNCONT)
     type(z_DNS) :: SC(MAXNCONT)
     type(z_CSR) :: HMC(MAXNCONT)
     type(z_CSR) :: SMC(MAXNCONT)
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

  subroutine set_readOldSGF(negf,logic) 
    type(Tnegf) :: negf
    logical :: logic

    negf%ReadoldSGF=logic
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
    type(z_CSR), target :: H
    type(z_CSR), optional, target :: S

    !call create(negf%H,H%nrow,H%ncol,H%nnz)
    !negf%H%nzval = H%nzval
    !negf%H%colind = H%colind
    !negf%H%rowpnt = H%rowpnt
    negf%H => H      

    if (present(S)) then
       negf%isSid=.false.
       !call create(negf%S,S%nrow,S%ncol,S%nnz)
       !negf%S%nzval = S%nzval
       !negf%S%colind = S%colind
       !negf%S%rowpnt = S%rowpnt
       negf%S => S
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
     negf%DorE = 'D'             ! Density or En.Density

     negf%ReadoldSGF = .false.   ! Read computed Surface G.F.
     negf%FictCont = .false.     ! Ficticious contact 

     negf%mu = 0.d0            ! Potenziale elettrico
     negf%efermi= 0.d0         ! Energia di Fermi dei contatti
     negf%contact_DOS = 0.d0   ! Ficticious contact DOS

    ! negf%nLdos = 0.d0                ! Number of LDOS intervals
    ! negf%LDOS(:,:) = 0.d0    ! LDOS intervals

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
     
   end subroutine set_defaults

end module lib_param
 




