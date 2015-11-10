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


module phph

  use ln_precision, only : dp
  use globals
  use ln_allocation
  use mat_def

  implicit none
  private

  public :: Tphph
  public :: init_phph

  type Tphph
    integer :: numatoms
    type(r_CSR) :: T3    ! T[i,i,j]=T[i,j]
    type(r_CSR), dimension(:), allocatable :: T4  !T[i,i,j,k]=T(i)[j,k]

    integer :: scba_iterations
    integer :: scba_iter

    logical :: include_phph
    logical :: cubic
    logical :: quartic
    logical :: Selfene_Gr
    logical :: Selfene_Gless
    logical :: Selfene_Hilb
    logical :: memory
    logical :: check
  end type Tphph

contains

  subroutine init_phph(phph,numatoms,order)
    Type(Tphph) :: phph
    integer, intent(in) :: numatoms
    integer, intent(in) :: order

    integer :: error, i 

    phph%numatoms = numatoms
    phph%scba_iterations = 0 ! starts from 0  
    phph%scba_iter = 0       ! initialize at 0
    phph%include_phph =.false.
    phph%cubic = .false.
    phph%quartic = .false.
    
    if (order == 3 .or. order == 34) then 
       phph%include_phph =.true.
       phph%cubic = .true.
       call create(phph%T3, 3*numatoms, 3*numatoms, 0)
    end if 
    
    if (order == 4 .or. order == 34) then 
       phph%include_phph =.true.
       phph%quartic = .true.
       allocate(phph%T4(3*numatoms), stat=error)
       if (error /= 0) then
         write(*,*) "ALLOCATION ERROR"; STOP 
       end if
       do i = 1, 3*numatoms
         call create(phph%T4(i), 3*numatoms, 3*numatoms, 0)
       enddo
    end if

 
    phph%Selfene_Gr = .true.
    phph%Selfene_Gless = .true.
    phph%Selfene_Hilb = .true.
    
    phph%memory = .true.
    phph%check = .false. 

  end subroutine init_phph



end module phph


