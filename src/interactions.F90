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

!> The module implements an abstract class to interface different
!! many body interactions. 

module interactions

  use globals, only : LST
  use ln_precision, only : dp
  use mat_def, only : z_dns
  use ln_structure, only : TStruct_info

  implicit none
  private

  public :: interaction

  type, abstract :: Interaction

    character(len=LST) :: descriptor
    !> Maximum number of SCBA iterations. 
    !! corresponds to no iterations (self energy is not calculated)
    integer :: scba_niter = 0
    !> Keep track of SCBA iteration 
    integer :: scba_iter = 0
    !> SCBA Tolerance
    real(dp) :: scba_tol = 1.0d-7
    !> Number of energy points from integration grid
    !integer :: en_npoints = 0
    !> Buffer for Gr

    !> System partitioning (as in TNEGF)
    type(TStruct_info) :: struct

  contains

    procedure(abst_add_sigma_r), deferred :: add_sigma_r
    procedure(abst_get_sigma_n), deferred :: get_sigma_n
    procedure(abst_set_Gr), deferred :: set_Gr
    procedure(abst_set_Gn), deferred :: set_Gn

  end type Interaction

  abstract interface

    !> This interface should append
    !! the retarded self energy to ESH
    subroutine abst_add_sigma_r(this, esh)
      import :: interaction
      import :: z_dns
      class(interaction) :: this
      type(z_dns), dimension(:,:), allocatable, intent(inout) :: esh
    end subroutine abst_add_sigma_r

    !> Returns the lesser (n) Self Energy in block format
    !! @param [in] this: calling instance
    !! @param [in] struct: system structure
    !! @param [inout] blk_sigma_n: block dense sigma_n
    !! @param [in] ie: index of energy point
    subroutine abst_get_sigma_n(this, blk_sigma_n, en_index)
      import :: interaction
      import :: z_dns
      class(interaction) :: this
      type(z_dns), dimension(:,:), allocatable, intent(inout) :: blk_sigma_n
      integer, intent(in) :: en_index
    end subroutine abst_get_sigma_n

    !> Give the Gr at given energy point to the interaction
    subroutine abst_set_Gr(this, Gr, en_index)
      import :: interaction
      import :: z_dns
      class(interaction) :: this
      type(z_dns), dimension(:,:), allocatable, intent(in) :: Gr
      integer :: en_index
    end subroutine abst_set_Gr

    !> Give the Gn at given energy point to the interaction
    subroutine abst_set_Gn(this, Gn, en_index)
      import :: interaction
      import :: z_dns
      class(interaction) :: this
      type(z_dns), dimension(:,:), allocatable, intent(in) :: Gn
      integer :: en_index
    end subroutine abst_set_Gn


  end interface

contains

    !> Initialize information needed for buffering G on memory or disk
    !  Now it only pass the number of energy grid points but it 
    !  could turn in something more complicated (e.g. an energy path object)
!!$    subroutine init_Gbuffer(this, en_npoints)
!!$      class(interaction) :: this
!!$      integer, intent(in) :: en_npoints
!!$
!!$      this%en_npoints = en_npoints
!!$      
!!$      
!!$    end subroutine init_Gbuffer

end module interactions
