module ln_inelastic
  use interactions
  use ln_cache
  implicit none

  public :: TInelastic

  type, abstract, extends(TInteraction) :: TInelastic

    !> sigma_r and sigma_n
    class(TMatrixCache), allocatable :: sigma_r
    class(TMatrixCache), allocatable :: sigma_n
    !> Gr and Gn are pointer alias stored in negf.
    class(TMatrixCache), pointer :: G_r => null()
    class(TMatrixCache), pointer :: G_n => null()

    contains

    procedure, non_overridable :: set_Gr_pointer
    procedure, non_overridable :: set_Gn_pointer

  end type TInelastic

  contains


  !> Set the Gr pointe
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
