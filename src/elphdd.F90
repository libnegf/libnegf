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

module elphdd

  use ln_precision, only : dp
  use interactions, only : interaction
  use ln_allocation, only : log_allocate, log_deallocate
  
  implicit none
  private

  public :: ElPhonDephD

  type, extends(interaction) :: ElPhonDephD

    !> Coupling squared per each orbital , dimension energy^2
    real(dp), allocatable, dimension(:) :: coupling
    complex(dp), allocatable, dimension(:) :: sigma_r
    complex(dp), allocatable, dimension(:) :: sigma_n
    integer :: nummodes
  contains
    procedure :: destroy  

  end type ElPhonDephD

  interface ElPhonDephD
    module procedure init_ElPhonDephD
  end interface ElPhonDephD

contains

  !>
  ! Constructor for el-ph dephasing diagonal model
  ! @param coupling: coupling (energy units) 
  ! @param niter: fixed number of scba iterations
  ! @param tol: scba tolerance
  function init_ElPhonDephD(coupling, niter, tol) result(this)

    type(ElPhonDephD) :: this
    real(dp), dimension(:), allocatable :: coupling
    real(dp), dimension(:), allocatable :: sigma_r
    real(dp), dimension(:), allocatable :: sigma_n
    integer, intent(in) :: niter
    real(dp), intent(in) :: tol

    this%descriptor = &
        & "Electron-Phonon dephasing model in fully diagonal model"

    call log_allocate(this%coupling, size(coupling))
    call log_allocate(this%sigma_r, size(coupling))
    call log_allocate(this%sigma_n, size(coupling))
    this%coupling = coupling * coupling
    this%scba_niter = niter
    this%scba_tol = tol
    this%sigma_r = 0.d0
    this%sigma_n = 0.d0
    this%nummodes = 1  !Single 0eV mode (Actually n localized modes, but we 
                       !treat them all contemporary)

  end function init_ElPhonDephD

  subroutine destroy(self)
    class(ElPhonDephD) :: self

    call log_deallocate(self%coupling)
    call log_deallocate(self%sigma_r)
    call log_deallocate(self%sigma_n)

  end subroutine destroy

end module elphdd
