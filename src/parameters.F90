module parameters

  use precision
  use constants, only : HAR,ATU

implicit none
private

  public :: init_defaults, Tparam

  integer, public, parameter :: MAXNCONT=10

  type Tparam
     real(dp)  :: hartree
     real(dp)  :: a_u
     integer   :: verbose
     integer   :: ncont
     integer   :: iatm(2)
     integer   :: N_omega
     integer   :: Np(4)
     integer   :: nPoles
     integer   :: niter

     integer   :: iatc(3,MAXNCONT)
     integer   :: ncdim(MAXNCONT)
     integer   :: mbound_end(MAXNCONT)

     real(dp)  :: Temp   
     real(dp)  :: Efermi(MAXNCONT)
     real(dp)  :: mu(MAXNCONT)
     real(dp)  :: DOS(MAXNCONT)
     real(dp)  :: LmbMax
     real(dp)  :: delta
     real(dp)  :: Elow
     
     logical   :: cluster
     logical   :: Readold
     logical   :: FictCont(MAXNCONT)   
  end type Tparam



contains

subroutine init_defaults(param)
  
  type(Tparam) :: param

 param%verbose=0
 param%iatc(:,:)=0
 param%iatm(:)=0
 param%ncdim(:)=0
 param%mbound_end(:)=0
 param%Np(:)=20
 param%N_omega=35
 param%nPoles=0

 param%hartree=HAR      ! this sets the conversion factors as set 
 param%a_u=ATU           ! in constants.F90
 param%Temp=0.0_dp
 param%LmbMax=0.50_dp
 param%Efermi(:)=0.0_dp
 param%Elow=-50.0_dp
 param%DOS(:)=0.0_dp
 param%delta=1e-4_dp
 param%mu(:)=0.0_dp

 param%FictCont(:)=.false.
 param%Readold=.false.
 param%cluster=.false.

end subroutine



end module parameters
