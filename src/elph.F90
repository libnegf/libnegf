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


module elph

  use ln_precision, only : dp
  use globals
  use ln_allocation

  implicit none
  private

  public :: Telph
  public :: init_elph_1


  !> This type contains information describing different electron phonon 
  !! models: input parameters, temporary data and output quantities
  type Telph
    !> Describe the model implemented. Currently supported:
    !! 0 : dummy model, no electron-phonon interactions
    !! 1 : electron phonon dephasing limit
    !!     Assumes elastic scattering and local coupling
    integer :: model = 0
    !> Diagonal coupling. Used in local coupling models (1)
    !! Note: it is assumed to be squared (units energy^2)
    real(dp), allocatable, dimension(:) :: coupling_array
    !> Diagonal elelents of retarded self energy. Used only in model (1)
    complex(dp), allocatable, dimension(:) :: diag_sigma_r


    !> Specify if the model is diagonal only. Depending on the model
    !! this will imply different approximations
    logical :: diagonal
    !> Number of active modes
    integer :: nummodes
     
    !> SCBA option: number of iterations
    !! 0 corresponds to no iterations (self energy is not calculated)
    integer :: scba_niter = 0
    
    !> Keep track of SCBA iteration 
    integer :: scba_iter = 0

    integer :: numselmodes
    logical, dimension(:), pointer :: selmodes => null()
    real(dp), dimension(:), pointer :: Wq => null()
    real(dp), dimension(:), pointer :: Nq => null()
 
    real(dp), dimension(:), pointer :: Mq => null()

    real(dp) :: self_range(2) !Energy interval for Sigma< 
    real(dp) :: Erange(2)     !Integration interval

    integer :: scba_iterations
    logical :: Selfene_Gr
    logical :: Selfene_Gless
    logical :: Selfene_Hilb
    logical :: memory
    logical :: check
  end type Telph

contains

  subroutine init_elph_1(elph, coupling, niter)
    Type(Telph) :: elph
    real(dp), dimension(:), allocatable, intent(in) :: coupling
    integer :: niter

    elph%model = 1
    elph%coupling_array = coupling
    elph%diagonal = .true.
    elph%scba_niter = niter

  end subroutine init_elph_1

  subroutine init_elph(elph,nummodes)
    Type(Telph) :: elph
    integer, optional :: nummodes

    if (present(nummodes)) then
       elph%nummodes = nummodes
       elph%numselmodes = nummodes
    else
       elph%nummodes = 0
       elph%numselmodes = 0
    endif
    elph%scba_iterations = 0 ! starts from 0  
    elph%scba_iter = 0       ! initialize at 0

    if (elph%nummodes .eq. 0) then
      
      elph%diagonal = .false.
      elph%Selfene_Gr = .false.
      elph%Selfene_Gless =.false.
      elph%Selfene_Hilb = .false.
     
    else

      call log_allocatep(elph%selmodes, elph%nummodes)
      call log_allocatep(elph%Wq, elph%nummodes)
      call log_allocatep(elph%Nq, elph%nummodes)
      call log_allocatep(elph%Mq, elph%nummodes)
 
      elph%Wq = 0.0_dp/27.2114_dp  !0.050_dp  ! 50 meV phonon
      elph%Nq = 0.0_dp              ! should be bose-einstain
      elph%Mq = 0.0_dp/27.2114_dp    ! 100 meV phonon coupling
      elph%selmodes = .true.
 
      elph%diagonal = .true.
      elph%Selfene_Gr = .true.
      elph%Selfene_Gless = .true.
      elph%Selfene_Hilb = .true.
      
      elph%memory = .true.
      elph%check = .false. 
    endif

  end subroutine init_elph



end module elph


!*******
!*******

!subroutine pippo(n)
!   integer, intent(in) :: n
!
!   integer, dimension(:), allocatable :: A
!   real(dp), dimension(:), pointer :: P
!   double DEPRECATED
!
!   allocate(A(n))

!   call log_allocate(A,n)     ! safe allocation in libNEGF

!   operazioni su A
!   call sub(A)

!   call log_deallocate(A)

!end subroutine pippo


!subroutine sub(B)
!   integer, dimension(:) :: B

!   integer :: i

!   do i = 1, size(B)
!      B(i) = 5
!   end do

! end subroutine sub


!subroutine sub2(B,n)
!   integer :: B(*)
!   integer :: n

!   integer :: i

!   do i = 1, n
!      B(i) = 5
!   end do

! end subroutine sub2

