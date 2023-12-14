module scba
  use ln_precision, only : dp
  use interactions
  use ln_elastic, only : TElastic
  use ln_inelastic, only : TInelastic
  use mat_def, only : z_CSR, z_DNS, create, destroy
  use sparsekit_drv, only : trace, clone

  type TScbaDriver
    !> Keep track of SCBA iteration
    integer :: scba_iter = 0

    !> SCBA Tolerance (Exact meaning may depend on model)
    real(dp) :: scba_tol = 1.0d-7

    !> SCBA error set for output
    real(dp) :: scba_err = 0.0_dp

    !> Holds the value of the previous Gn
    type(z_CSR), allocatable :: Mat_old

    !> Holds the value of the previous current
    real(dp), dimension(:), allocatable :: J_old

    !> internal status
    logical :: converged = .false.

    !> internal status
    logical :: do_write = .false.

    contains
    procedure :: init => scba_init
    procedure :: destroy => scba_destroy
    procedure :: set_scba_iter
    procedure :: scba_error
    procedure :: is_converged
    procedure :: check_Mat_convergence
    procedure :: check_J_convergence
  end type TScbaDriver

  type, extends(TScbaDriver) :: TScbaDriverElastic
    contains
    procedure :: set_scba_iter => set_scba_iter_elastic
  end type TScbaDriverElastic

  type, extends(TScbaDriver) :: TScbaDriverInelastic
    contains
    procedure :: set_scba_iter => set_scba_iter_inelastic
  end type TScbaDriverInelastic

  contains

  !> Initialize the SCBA Driver
  subroutine scba_init(this, tol, dowrite)
    class(TScbaDriver) :: this
    real(dp), intent(in) :: tol
    logical, intent(in) :: dowrite
    this%scba_tol = tol
    this%do_write = dowrite
    allocate(this%Mat_old)
  end subroutine scba_init

  !> Release internal space
  subroutine scba_destroy(this)
    class(TScbaDriver) :: this
    integer :: ii
    if (allocated(this%Mat_old%nzval)) call destroy(this%Mat_old)
    if (allocated(this%Mat_old)) deallocate(this%Mat_old)
    if (allocated(this%J_old)) deallocate(this%J_old)
    this%scba_iter = 0
    this%scba_err = 0.0_dp
    this%converged = .false.
  end subroutine scba_destroy


  !> Sets the scba_iteration in all interactions
  subroutine set_scba_iter(this, iter, interactList)
    class(TScbaDriver) :: this
    integer, intent(in) :: iter
    type(TInteractionList) :: interactList

    type(TInteractionNode), pointer :: it
    this%scba_iter = iter
    it => interactList%first
    do while (associated(it))
      call it%inter%set_scba_iter(iter)
      it => it%next
    end do
  end subroutine set_scba_iter

  !> get the actual error atteined
  function scba_error(this) result(error)
    class(TScbaDriver) :: this
    real(dp) :: error

    error = this%scba_err
  end function scba_error

  !> Check convergence on Gn matrix blocks
  subroutine check_Mat_convergence(this, Mat)
    class(TScbaDriver) :: this
    type(z_CSR), intent(in) :: Mat

    if (.not.allocated(this%Mat_old)) then
      stop 'ERROR: TScbaDriver must be initialized first'
    end if
    if (.not.allocated(this%Mat_old%nzval)) then
      call clone(Mat, this%Mat_old)
      return
    end if

    this%scba_err = maxval(abs(Mat%nzval-this%Mat_old%nzval))
    if (this%scba_err < this%scba_tol) then
      this%converged = .true.
      if (this%do_write) then
         write(*,*) "SCBA loop converged in",this%scba_iter,&
                  & " iterations with error", this%scba_err
      end if
    else
      this%converged = .false.
      if (this%do_write) then
         write(*,*) "SCBA iteration",this%scba_iter," error", this%scba_err
      end if
    end if
    call destroy(this%Mat_old)
    call clone(Mat, this%Mat_old)
  end subroutine check_Mat_convergence

  !> Check convergence on the 1st layer current
  subroutine check_J_convergence(this, J)
    class(TScbaDriver) :: this
    real(dp), dimension(:), intent(in) :: J

    if (.not.allocated(this%J_old)) then
      allocate(this%J_old, source=J)
      return
    end if

    this%scba_err = maxval(abs(this%J_old - J))
    if (this%scba_err < this%scba_tol) then
      this%converged = .true.
      if (this%do_write) then
         write(*,*) "SCBA loop converged in",this%scba_iter,&
                  & " iterations with error", this%scba_err
      end if
    else
      this%converged = .false.
      if (this%do_write) then
         write(*,*) "SCBA iteration",this%scba_iter," error", this%scba_err
      end if
    end if
    this%J_old = J
  end subroutine check_J_convergence

  function is_converged(this) result (conv)
    class(TScbaDriver) :: this
    logical :: conv

    conv = this%converged
  end function is_converged


  !> Sets the scba_iteration in all elastic interactions
  subroutine set_scba_iter_elastic(this, iter, interactList)
    class(TScbaDriverElastic) :: this
    integer, intent(in) :: iter
    type(TInteractionList) :: interactList

    type(TInteractionNode), pointer :: it
    this%scba_iter = iter
    it => interactList%first
    do while (associated(it))
      select type(pInter => it%inter)
      class is (TElastic)
        call pInter%set_scba_iter(iter)
      end select
      it => it%next
    end do
  end subroutine set_scba_iter_elastic

  !> Sets the scba_iteration in all inelastic interactions
  subroutine set_scba_iter_inelastic(this, iter, interactList)
    class(TScbaDriverInelastic) :: this
    integer, intent(in) :: iter
    type(TInteractionList) :: interactList

    type(TInteractionNode), pointer :: it
    this%scba_iter = iter
    it => interactList%first
    do while (associated(it))
      select type(pInter => it%inter)
      class is (TInelastic)
        call pInter%set_scba_iter(iter)
      end select
      it => it%next
    end do
  end subroutine set_scba_iter_inelastic


end module scba

