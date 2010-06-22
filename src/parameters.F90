module parameters

  use precision
  use constants, only : HAR,ATU

implicit none
private

  public :: init_defaults, Tparam

  integer, public, parameter :: MAXNCONT=10

  type Tparam
     real(dp)  :: hartree         ! 
     real(dp)  :: a_u             ! 
     integer   :: verbose         ! livello di verbosita`
     integer   :: Np(4)           ! Numero di punti di quadratura
     integer   :: nPoles          ! Numero di poli inclusi
     integer   :: N_omega         ! Numero di kT per l'integrazione (10)  
     integer   :: activecont      ! Active contact
     integer   :: iter

     real(dp)  :: Temp            ! Temperatura (Fermi)
     real(dp)  :: Efermi(MAXNCONT)! Energia di Fermi dei contatti
     real(dp)  :: mu(MAXNCONT)    ! Potenziale elettrico
     real(dp)  :: DOS(MAXNCONT)
     real(dp)  :: delta           ! delta per le G.F. (1e-5)
     real(dp)  :: Elow            ! Lowest energy (-50 eV)
     
     logical   :: cluster         
     logical   :: ReadoldSGF         
     logical   :: FictCont(MAXNCONT)
     character(1) :: DorE           !Density or En.Density
   
     integer :: contdim(MAXNCONT)  !contact dimension
     integer :: surfdim(MAXNCONT)  !surface dimension

  end type Tparam



contains

subroutine init_defaults(param)
  
  type(Tparam) :: param

 param%verbose=0
 param%Np(:)=20
 param%N_omega=10
 param%nPoles=0
 param%iter=1

 param%hartree=HAR      ! this sets the conversion factors as set 
 param%a_u=ATU           ! in constants.F90
 param%Temp=0.0_dp
 param%Efermi(:)=0.0_dp
 param%Elow=-50.0_dp/HAR
 param%DOS(:)=0.05_dp
 param%delta=1e-4_dp/HAR
 param%mu(:)=0.0_dp

 param%FictCont(:)=.false.
 param%ReadoldSGF=.false.
 param%cluster=.false.
 param%DorE='D'

end subroutine



end module parameters
