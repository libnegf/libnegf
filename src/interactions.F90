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

  public :: TInteraction
  public :: TInteractionList, TInteractionNode
  public :: get_max_wq
  public :: get_max_niter

  type, abstract :: TInteraction
    !> Textual descriptor for output
    character(len=LST) :: descriptor
    !> Maximum number of SCBA iterations
    !! corresponds to no iterations (self energy is not calculated)
    integer :: scba_niter = 0
    !> SCBA iteration (set from outside)
    integer :: scba_iter = 0
    !> Energy of the mode (what about wq(k) ??)
    real(dp) :: wq = 0.0_dp

    !> System partitioning (as in TNEGF)
    type(TStruct_info) :: struct

  contains

    procedure, non_overridable :: set_scba_iter
    procedure(abst_add_sigma_r), deferred :: add_sigma_r
    procedure(abst_add_sigma_n), deferred :: add_sigma_n
    procedure(abst_get_sigma_n_blk), deferred, private :: get_sigma_n_blk
    procedure(abst_get_sigma_n_mat), deferred, private :: get_sigma_n_mat
    generic :: get_sigma_n => get_sigma_n_blk, get_sigma_n_mat
    procedure(abst_set_Gr), deferred :: set_Gr
    procedure(abst_set_Gn), deferred :: set_Gn
    procedure(abst_comp_Sigma_r), deferred :: compute_Sigma_r
    procedure(abst_comp_Sigma_n), deferred :: compute_Sigma_n
    procedure(abst_destroy_Sigma_r), deferred :: destroy_Sigma_r
    procedure(abst_destroy_Sigma_n), deferred :: destroy_Sigma_n
    procedure(abst_destroy), deferred :: destroy

  end type TInteraction

  !-----------------------------------------------------------------------------
  ! derived type to create list of interaction objects
  type TInteractionList
    integer :: counter = 0
    type(TInteractionNode), pointer :: first => null()
    type(TInteractionNode), pointer :: curr => null()
    contains
    procedure :: add => add
    procedure :: destroy => destroy
  end type TInteractionList

  type TInteractionNode
    class(TInteraction), allocatable :: inter
    type(TInteractionNode), pointer :: next => null()
  end type TInteractionNode

  abstract interface

    !> This interface should append
    !! the retarded self energy to ESH
    subroutine abst_add_sigma_r(this, esh, en_index, k_index, spin)
      import :: TInteraction
      import :: z_dns
      class(TInteraction) :: this
      type(z_dns), dimension(:,:), intent(inout) :: esh
      integer, intent(in), optional  :: en_index
      integer, intent(in), optional  :: k_index
      integer, intent(in), optional  :: spin
    end subroutine abst_add_sigma_r

    !> This interface should append
    !! the retarded self energy to sigma
    subroutine abst_add_sigma_n(this, sigma, en_index, k_index, spin)
      import :: TInteraction
      import :: z_dns
      class(TInteraction) :: this
      type(z_dns), dimension(:,:), intent(inout) :: sigma
      integer, intent(in), optional  :: en_index
      integer, intent(in), optional  :: k_index
      integer, intent(in), optional  :: spin
    end subroutine abst_add_sigma_n

    !> Returns the lesser (n) Self Energy in block format
    !! @param [in] this: calling instance
    !! @param [in] struct: system structure
    !! @param [inout] blk_sigma_n: block dense sigma_n
    !! @param [in] ie: index of energy point
    !! @param [in] ie: index of k point
    subroutine abst_get_sigma_n_blk(this, blk_sigma_n, en_index, k_index, spin)
      import :: TInteraction
      import :: z_dns
      class(TInteraction) :: this
      type(z_dns), dimension(:,:), intent(inout) :: blk_sigma_n
      integer, intent(in), optional :: en_index
      integer, intent(in), optional :: k_index
      integer, intent(in), optional :: spin
    end subroutine abst_get_sigma_n_blk

    !> Returns the lesser (n) Self Energy in block format
    !! @param [in] this: calling instance
    !! @param [in] struct: system structure
    !! @param [inout] blk_sigma_n: block dense sigma_n
    !! @param [in] ie: index of energy point
    !! @param [in] ie: index of k point
    subroutine abst_get_sigma_n_mat(this, sigma_n, ii, jj, en_index, k_index, spin)
      import :: TInteraction
      import :: z_dns
      class(TInteraction) :: this
      type(z_dns), intent(inout) :: sigma_n
      integer, intent(in) :: ii
      integer, intent(in) :: jj
      integer, intent(in), optional :: en_index
      integer, intent(in), optional :: k_index
      integer, intent(in), optional :: spin
    end subroutine abst_get_sigma_n_mat


    !> Give the Gr at given energy point to the TInteraction
    subroutine abst_set_Gr(this, Gr, en_index, k_index, spin)
      import :: TInteraction
      import :: z_dns
      class(TInteraction) :: this
      type(z_dns), dimension(:,:), intent(in) :: Gr
      integer, intent(in), optional  :: en_index
      integer, intent(in), optional  :: k_index
      integer, intent(in), optional  :: spin
    end subroutine abst_set_Gr

    !> Give the Gn at given energy point to the TInteraction
    subroutine abst_set_Gn(this, Gn, en_index, k_index, spin)
      import :: TInteraction
      import :: z_dns
      class(TInteraction) :: this
      type(z_dns), dimension(:,:), intent(in) :: Gn
      integer, intent(in), optional  :: en_index
      integer, intent(in), optional  :: k_index
      integer, intent(in), optional  :: spin
    end subroutine abst_set_Gn

    !>  Compute Sigma_n : necessary for inelastic
    subroutine abst_comp_Sigma_n(this, en_index, k_index, spin)
      import :: TInteraction
      class(TInteraction) :: this
      integer, intent(in), optional  :: en_index
      integer, intent(in), optional  :: k_index
      integer, intent(in), optional  :: spin
    end subroutine abst_comp_Sigma_n

    !>  Compute Sigma_r : necessary for inelastic
    subroutine abst_comp_Sigma_r(this, en_index, k_index, spin)
      import :: TInteraction
      class(TInteraction) :: this
      integer, intent(in), optional  :: en_index
      integer, intent(in), optional  :: k_index
      integer, intent(in), optional  :: spin
    end subroutine abst_comp_Sigma_r

    !>  Destroy Sigma_r : cleanup memory
    subroutine abst_destroy_Sigma_r(this)
      import :: TInteraction
      class(TInteraction) :: this
    end subroutine abst_destroy_Sigma_r

    !>  Destroy Sigma_n : cleanup memory
    subroutine abst_destroy_Sigma_n(this)
      import :: TInteraction
      class(TInteraction) :: this
    end subroutine abst_destroy_Sigma_n

    !>  Destroy object : cleanup memory
    subroutine abst_destroy(this)
      import :: TInteraction
      class(TInteraction) :: this
    end subroutine abst_destroy

  end interface

  contains

  subroutine set_scba_iter(this, scba_iter)
    class(TInteraction) :: this
    integer, intent(in) :: scba_iter
    this%scba_iter = scba_iter
  end subroutine set_scba_iter

  function get_max_wq(list) result(maxwq)
    type(TInteractionList), intent(in) :: list
    real(dp) :: maxwq

    type(TInteractionNode), pointer :: it
    it => list%first
    maxwq = 0.0_dp
    do while (associated(it))
      if (it%inter%wq > maxwq) then
         maxwq = it%inter%wq
      end if
      it => it%next
    end do
  end function get_max_wq

  function get_max_niter(list) result (max_niter)
    type(TInteractionList), intent(in) :: list
    integer :: max_niter
    type(TInteractionNode), pointer :: it
    max_niter = 0
    it => list%first
    do while (associated(it))
      if (it%inter%scba_niter > max_niter) then
         max_niter = it%inter%scba_niter
      end if
      it => it%next
    end do
  end function get_max_niter

  ! Interaction list methods
  subroutine add(this, node)
    class(TInteractionList) :: this
    type(TInteractionNode), pointer, intent(out) :: node

    allocate(node)
    if (.not.associated(this%first)) then
       this%first => node
       this%curr => node
       this%counter = 1
    else
       this%curr%next => node
       this%curr => node
       this%counter = this%counter + 1
    end if
  end subroutine add

  subroutine destroy(this)
    class(TInteractionList) :: this

    type(TInteractionNode), pointer :: it

    if (.not.associated(this%first)) then
       return
    else
       this%curr => this%first
       do while (associated(this%curr))
          it => this%curr
          this%curr => this%curr%next
          if (allocated(it%inter)) then
             call it%inter%destroy()
             deallocate(it%inter)
          end if
          deallocate(it)
       end do
    end if
    this%first=>null()
    this%curr=>null()
  end subroutine destroy

end module interactions
