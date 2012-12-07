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
  public :: init_elph

  type Telph
    integer :: nummodes
    integer :: numselmodes
    logical, dimension(:), pointer :: selmodes => null()
    real(dp), dimension(:), pointer :: Wq => null()
    real(dp), dimension(:), pointer :: Nq => null()
 
    real(dp), dimension(:), pointer :: Mq => null()

    real(dp) :: self_range(2) !Energy interval for Sigma< 
    real(dp) :: Erange(2)     !Integration interval

    integer :: scba_iterations
    integer :: scba_iter
    logical :: diagonal
 
    logical :: Selfene_Gr
    logical :: Selfene_Gless
    logical :: Selfene_Hilb
  end type Telph

contains

  subroutine init_elph(elph)
    Type(Telph) :: elph
    
    elph%nummodes = 0
    elph%numselmodes = 0
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
 
      elph%Wq(1) = 0.0_dp/27.2114_dp  !0.050_dp  ! 50 meV phonon
      elph%Nq(1) = 0.0_dp              ! should be bose-einstain
      elph%Mq(1) = 0.0_dp/27.2114_dp    ! 100 meV phonon coupling
      elph%selmodes(1) = .true.
 
      elph%diagonal = .true.
      elph%Selfene_Gr = .true.
      elph%Selfene_Gless = .true.
      elph%Selfene_Hilb = .true.
    
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

