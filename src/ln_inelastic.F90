!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Greens functions.         !
!! Copyright (C) 2012 - 2026                                                !
!!                                                                          !
!! This file is part of libNEGF: a library for                              !
!! Non Equilibrium Green's Functions calculations                           !
!!                                                                          !
!! Developers: Alessandro Pecchia, Daniele Soccodato                        !
!! Former Contributors: Gabriele Penazzi, Luca Latessa, Aldo Di Carlo       !
!!                                                                          !
!! libNEGF is free software: you can redistribute and/or modify it          !
!! under the terms of the GNU Lesser General Public License as published    !
!! by the Free Software Foundation, either version 3 of the License, or     !
!! (at your option) any later version.                                      !
!!                                                                          !
!!  You should have received a copy of the GNU Lesser General Public        !
!!  License along with libNEGF.  If not, see                                !
!!  <http://www.gnu.org/licenses/>.                                         !
!!--------------------------------------------------------------------------!

module ln_inelastic
  use interactions
  use ln_cache
  implicit none

  public :: TInelastic

  type, abstract, extends(TInteraction) :: TInelastic

    !> sigma_r and sigma_n
    class(TMatrixCache), pointer :: sigma_r => null()
    class(TMatrixCache), pointer :: sigma_n => null()
    !> Gr and Gn are pointer alias stored in negf.
    class(TMatrixCache), pointer :: G_r => null()
    class(TMatrixCache), pointer :: G_n => null()

    !> Compute tri-diagonal block matrices
    logical :: tTridiagonal = .false.

    contains

    procedure, non_overridable :: set_Gr_pointer
    procedure, non_overridable :: set_Gn_pointer
    procedure(abst_init_Sigma_r), deferred :: init_Sigma_r
    procedure(abst_init_Sigma_n), deferred :: init_Sigma_n

  end type TInelastic

  contains

  !>  Init Sigma_r : necessary for inelastic
  subroutine abst_init_Sigma_r(this, spin)
    class(TInelastic) :: this
    integer, intent(in), optional :: spin
  end subroutine abst_init_Sigma_r

  !>  Init Sigma_n : necessary for inelastic
  subroutine abst_init_Sigma_n(this, spin)
    class(TInelastic) :: this
    integer, intent(in), optional :: spin
  end subroutine abst_init_Sigma_n


  !> Set the Gr pointer
  subroutine set_Gr_pointer(this, Gr)
    class(TInelastic), intent(inout) :: this
    class(TMatrixCache), pointer :: Gr
    this%G_r => Gr
  end subroutine set_Gr_pointer

  !> Set the Gn pointer
  subroutine set_Gn_pointer(this, Gn)
    class(TInelastic), intent(inout) :: this
    class(TMatrixCache), pointer :: Gn
    this%G_n => Gn
  end subroutine set_Gn_pointer

end module ln_inelastic
