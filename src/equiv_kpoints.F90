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

! Inelastic electron-phonon interactions
! This is the base class for inelastic interactions

module equiv_kpoints
  use ln_precision, only : dp
  use ln_allocation
  implicit none
  private
  
  public :: TEqPoint
  public :: TEqPointsArray
  public :: create
  public :: destroy
  public :: set

  !> For each kpoint, one set of equivalent points consists of a 3 x n_eq matrix
  type TEqPoint
    real(dp), dimension(:,:), allocatable :: points
  end type  TEqPoint

  !> Structure containing information about symmetrically equivalent points, 
  !> for each kpoint in kpoints
  type TEqPointsArray
    type(TEqPoint), dimension(:), allocatable :: EqPoints
  end type TEqPointsArray
 
  interface create
    module procedure :: create_equivalent_points
  end interface

  interface destroy
    module procedure ::  destroy_equivalent_points
  end interface
  
  interface set
    module procedure :: set_equivalent_points
  end interface

  contains

  !> Create the equivalent kpoint container 
  subroutine create_equivalent_points(mycontainer, equiv_kpoints)
    type(TEqPointsArray), allocatable, intent(inout) :: mycontainer
    type(TEqPointsArray), intent(in) :: equiv_kpoints

    integer :: ii, nk, n_eq

    if (allocated(mycontainer)) then
       error stop "Equivalent k-points container already created"
    end if   
    
    nk = size(equiv_kpoints%EqPoints)

    allocate(mycontainer)
    allocate(mycontainer%EqPoints(nk))

    do ii = 1, nk
      n_eq = size(equiv_kpoints%EqPoints(ii)%points, 2)
      call log_allocate(mycontainer%EqPoints(ii)%points, 3, n_eq)
      mycontainer%EqPoints(ii)%points = equiv_kpoints%EqPoints(ii)%points
    end do
        
  end subroutine create_equivalent_points

  !> Destroy the container 
  subroutine destroy_equivalent_points(eq_points)
    type(TEqPointsArray), allocatable :: eq_points

    integer :: ii, npoints

    if (allocated(eq_points)) then
      if (allocated(eq_points%EqPoints)) then
        npoints = size(eq_points%EqPoints)
        do ii = 1, npoints
          if (allocated(eq_points%EqPoints(ii)%points)) then
            call log_deallocate(eq_points%EqPoints(ii)%points)
          end if  
        end do
        deallocate(eq_points%EqPoints)
      end if
      deallocate(eq_points)
    end if

  end subroutine destroy_equivalent_points

  !> Define the equivalent k-point structure from an external array
  subroutine set_equivalent_points(mycontainer, equiv_kpoints, nkpoints, equiv_mult)
    type(TEqPointsArray), allocatable, intent(inout) :: mycontainer 
    real(dp), intent(in) :: equiv_kpoints(:,:)
    integer, intent(in) :: nkpoints
    integer, intent(in) :: equiv_mult(:)

    integer :: i, j, n_eq, begin, last

    if (allocated(mycontainer)) then
       call destroy(mycontainer)
    end if   
    allocate(mycontainer)
    allocate(mycontainer%EqPoints(nkpoints))

    !Eq. points of the first kpoint must be initialized outside the loop 
    n_eq = equiv_mult(1)
    !mycontainer%EqPoints(1)%n_eq = n_eq
    !if (allocated(mycontainer%EqPoints(1)%points)) then 
    !  call log_deallocate(mycontainer%EqPoints(1)%points)
    !endif
    call log_allocate(mycontainer%EqPoints(1)%points, 3, n_eq)

    begin = 1
    last = n_eq
    mycontainer%EqPoints(1)%points(:,:) = equiv_kpoints(:, begin:last)

    ! Now initialize all the remaining kpoints
    do i = 2, nkpoints
      n_eq = equiv_mult(i)
      !mycontainer%EqPoints(i)%n_eq = n_eq
      !if (allocated(mycontainer%EqPoints(i)%points)) then
      !  call log_deallocate(mycontainer%EqPoints(i)%points)
      !endif
      call log_allocate(mycontainer%EqPoints(i)%points, 3, n_eq)

      begin = last + 1
      last = begin + n_eq - 1
      mycontainer%EqPoints(i)%points(:,:) = equiv_kpoints(:, begin:last)
    enddo

    !mycontainer%present = .true.

  end subroutine set_equivalent_points

end module equiv_kpoints

