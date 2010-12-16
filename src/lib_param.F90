module lib_param

  use precision, only : dp
  use globals
  use mat_def
  use structure, only : TStruct_info
  use input_output

  implicit none
  private

  public :: Tnegf
  public :: fill_parameters
  integer, public, parameter :: MAXNCONT=10

  type Tnegf
     
     integer :: verbose

     character(LST) :: file_re_H
     character(LST) :: file_im_H
     character(LST) :: file_re_S
     character(LST) :: file_im_S
     character(LST) :: file_struct
     character(1) :: DorE           !Density or En.Density
     type(format) :: form

     logical   :: ReadoldSGF
     logical   :: FictCont(MAXNCONT)

     real(dp) :: mu(MAXNCONT)    ! Potenziale elettrico
     real(dp) :: Efermi(MAXNCONT)! Energia di Fermi dei contatti
     real(dp) :: contact_DOS(MAXNCONT)

     real(dp) :: mu_n
     real(dp) :: mu_p
     real(dp) :: Ec
     real(dp) :: Ev
     real(dp) :: DeltaEc
     real(dp) :: DeltaEv     
     real(dp) :: E
     real(dp) :: dos
     real(dp) :: eneconv

     type(z_CSR) :: H
     type(z_CSR) :: S
     type(z_CSR) :: HM
     type(z_CSR) :: SM     
     type(z_DNS) :: HC(MAXNCONT)
     type(z_DNS) :: SC(MAXNCONT)
     type(z_CSR) :: HMC(MAXNCONT)
     type(z_CSR) :: SMC(MAXNCONT)
     type(z_CSR) :: Gr
     type(z_CSR) :: rho
     type(z_CSR) :: rho_eps
     logical    :: isSid

     type(TStruct_Info) :: str

     real(dp) :: delta
     real(dp) :: Emin
     real(dp) :: Emax
     real(dp) :: Estep
     real(dp) :: Temp


     integer :: Np_n(2)
     integer :: Np_p(2)
     !integer :: Np(4)   
     integer :: n_kt
     integer :: n_poles
     integer :: iteration
     integer :: activecont 
     integer :: N_omega         ! Numero di kT per l'integrazione
     real(dp) :: spin   

  

  end type Tnegf

contains
  
  subroutine fill_parameters(negf, verbose, DorE, form, mu_n, mu_p, Ec, Ev, DeltaEc, DeltaEv, E, dos, eneconv, &
       ReadoldSGF, FictCont, mu, Efermi, contact_DOS, H, S, HM, SM, HC, SC, HMC, SMC, Gr, rho, rho_eps, isSid, &
       str, delta, Emin, Emax, Estep, Temp, Np_n, Np_p, n_kt, n_poles, iteration, activecont, N_omega, spin)


    type(Tnegf) :: negf
    integer :: verbose
    character(1) :: DorE           !Density or En.Density
    type(format) :: form
    logical   :: ReadoldSGF
    logical   :: FictCont(MAXNCONT)
    real(dp) :: mu(MAXNCONT)    ! Potenziale elettrico
    real(dp) :: Efermi(MAXNCONT)! Energia di Fermi dei contatti
    real(dp) :: contact_DOS(MAXNCONT)

    real(dp) :: mu_n
    real(dp) :: mu_p
    real(dp) :: Ec
    real(dp) :: Ev
    real(dp) :: DeltaEc
    real(dp) :: DeltaEv     
    real(dp) :: E
    real(dp) :: dos
    real(dp) :: eneconv

    type(z_CSR) :: H
    type(z_CSR) :: S
    type(z_CSR) :: HM
    type(z_CSR) :: SM     
    type(z_DNS) :: HC(MAXNCONT)
    type(z_DNS) :: SC(MAXNCONT)
    type(z_CSR) :: HMC(MAXNCONT)
    type(z_CSR) :: SMC(MAXNCONT)
    type(z_CSR) :: Gr
    type(z_CSR) :: rho
    type(z_CSR) :: rho_eps
    logical    :: isSid

    type(TStruct_Info) :: str

    real(dp) :: delta
    real(dp) :: Emin
    real(dp) :: Emax
    real(dp) :: Estep
    real(dp) :: Temp


    integer :: Np_n(2)
    integer :: Np_p(2)   
    integer :: n_kt
    integer :: n_poles
    integer :: iteration
    integer :: activecont 
    integer :: N_omega         ! Numero di kT per l'integrazione
    real(dp) :: spin   


    negf%verbose = verbose
    negf%DorE = DorE
    negf%form = form
    negf%ReadoldSGF = ReadoldSGF
    negf%FictCont = FictCont
    negf%mu(MAXNCONT) = mu(MAXNCONT)
    negf%Efermi(MAXNCONT) = Efermi(MAXNCONT)
    negf%contact_DOS(MAXNCONT) = contact_DOS(MAXNCONT)

    negf%mu_n = mu_n
    negf%mu_p = mu_p
    negf%Ec = Ec 
    negf%Ev = Ev
    negf%DeltaEc = DeltaEc 
    negf%DeltaEv = DeltaEv 
    negf%E = E
    negf%dos = dos
    negf%eneconv = eneconv 

    negf%H = H
    negf%S = S
    negf%HM = HM
    negf%SM = SM    
    negf%HC = HC 
    negf%SC = SC 
    negf%HMC = HMC
    negf%SMC = SMC
    negf%Gr = Gr
    negf%rho = rho
    negf%rho_eps = rho_eps
    negf%isSid = isSid

    negf%str = str
    negf%delta = delta
    negf%Emin = Emin
    negf%Emax = Emax
    negf%Estep = Estep
    negf%Temp = Temp

    negf%Np_n(2) = Np_n(2)
    negf%Np_p(2) = Np_p(2)
    negf%n_kt = n_kt
    negf%n_poles = n_poles
    negf%iteration = iteration
    negf%activecont = activecont
    negf%N_omega = N_omega
    negf%spin = spin


  end subroutine fill_parameters

end module lib_param
