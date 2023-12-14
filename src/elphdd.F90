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

!> Diagonal elastic dephasing model

module elphdd

  use ln_precision, only : dp
  use interactions, only : TInteraction
  use ln_elastic, only : TElastic
  use ln_allocation, only : log_allocate, log_deallocate
  use ln_structure, only : TStruct_info
  use mat_def, only : z_dns, create

  implicit none
  private

  public :: ElPhonDephD, ElPhonDephD_create
  public :: ElPhonDephD_init

  type, extends(TElastic) :: ElPhonDephD

    !> Coupling squared per each orbital , dimension energy^2
    real(dp), allocatable, dimension(:) :: coupling
    !> Local diagonal representation of retarded self energy
    complex(dp), allocatable, dimension(:) :: sigma_r
    !> Local diagonal representation of lesser self energy
    complex(dp), allocatable, dimension(:) :: sigma_n
    !> Number of vibrational modes
    integer :: nummodes

  contains

    procedure :: add_sigma_r
    procedure :: add_sigma_n
    procedure :: get_sigma_n_blk
    procedure :: get_sigma_n_mat
    procedure :: set_Gr
    procedure :: set_Gn
    procedure :: compute_Sigma_r
    procedure :: compute_Sigma_n
    procedure :: destroy_Sigma_r
    procedure :: destroy_Sigma_n

  end type ElPhonDephD

contains

  subroutine ElPhonDephD_create(this)
    class(TInteraction), allocatable :: this
    allocate(ElPhonDephD::this)
  end subroutine ElPhonDephD_create
  !>
  ! Factory for el-ph dephasing diagonal model
  ! @param struct : contact/device partitioning infos
  ! @param coupling: coupling (energy units)
  ! @param niter: fixed number of scba iterations
  ! @param tol: scba tolerance
  subroutine ElPhonDephD_init(this, struct, coupling, niter)
    type(ElPhonDephD) :: this
    type(TStruct_info), intent(in) :: struct
    real(dp), dimension(:), intent(in) :: coupling
    integer, intent(in) :: niter

    this%descriptor = &
        & "Electron-Phonon dephasing model in fully diagonal model"

    call log_allocate(this%coupling, size(coupling))
    call log_allocate(this%sigma_r, size(coupling))
    call log_allocate(this%sigma_n, size(coupling))
    this%struct = struct
    this%coupling = coupling * coupling
    this%scba_niter = niter
    this%sigma_r = (0.0_dp, 0.0_dp)
    this%sigma_n = (0.0_dp, 0.0_dp)
    this%nummodes = 1  !Single mode (Actually n localized modes, but we
                       !treat them all together)
    this%wq = 0.0_dp   ! Zero energy mode

  end subroutine ElPhonDephD_init


  !> This interface should append
  !  the retarded self energy to ESH
  subroutine add_sigma_r(this, esh, en_index, k_index, spin)
    class(ElPhonDephD) :: this
    type(z_dns), dimension(:,:), intent(inout) :: esh
    integer, intent(in), optional :: en_index
    integer, intent(in), optional  :: k_index
    integer, intent(in), optional  :: spin

    integer :: n, nbl, ii, pl_start, pl_end

    if (this%scba_iter .eq. 0) return
    nbl = this%struct%num_PLs
    do n=1,nbl
      pl_start=this%struct%mat_PL_start(n)
      pl_end=this%struct%mat_PL_end(n)
      do ii = 1, pl_end - pl_start + 1
          ESH(n,n)%val(ii,ii) = ESH(n,n)%val(ii,ii) - &
            this%sigma_r(pl_start + ii - 1)
      end do
    end do

  end subroutine add_sigma_r

  !--------------------------------------------------------------------------
  !> This interface should append
  !  sigma_n to a passed self energy, sigma
  subroutine add_sigma_n(this, sigma, en_index, k_index, spin)
    class(ElPhonDephD) :: this
    type(z_dns), dimension(:,:), intent(inout) :: sigma
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin

    integer :: n, nbl, ii, pl_start, pl_end

    if (this%scba_iter .eq. 0) return

    nbl = this%struct%num_PLs
    do n=1,nbl
      pl_start=this%struct%mat_PL_start(n)
      pl_end=this%struct%mat_PL_end(n)
      do ii = 1, pl_end - pl_start + 1
          sigma(n,n)%val(ii,ii) = sigma(n,n)%val(ii,ii) + &
            this%sigma_n(pl_start + ii - 1)
      end do
    end do

  end subroutine add_sigma_n


  !> Returns the lesser (n) Self Energy in block format
  !
  subroutine get_sigma_n_blk(this, blk_sigma_n, en_index, k_index, spin)
    class(ElPhonDephD) :: this
    type(z_dns), dimension(:,:), intent(inout) :: blk_sigma_n
    integer, intent(in), optional  :: en_index
    integer, intent(in), optional  :: k_index
    integer, intent(in), optional  :: spin

    integer :: n, nbl, ii, nrow, pl_start, pl_end
    nbl = this%struct%num_PLs

    if (this%scba_iter .eq. 0) return

    do n = 1, nbl
      pl_start=this%struct%mat_PL_start(n)
      pl_end=this%struct%mat_PL_end(n)
      nrow = pl_end - pl_start + 1
      if (.not.allocated(blk_sigma_n(n,n)%val)) then
        call create(blk_sigma_n(n,n), nrow, nrow)
      end if
      blk_sigma_n(n,n)%val = (0.0_dp, 0.0_dp)
      do ii = 1, pl_end - pl_start + 1
         blk_sigma_n(n,n)%val(ii,ii) = this%sigma_n(pl_start+ii-1)
      end do
    enddo

  end subroutine get_sigma_n_blk

  subroutine get_sigma_n_mat(this, sigma_n, ii, jj, en_index, k_index, spin)
    class(ElPhonDephD) :: this
    type(z_dns), intent(inout) :: sigma_n
    integer, intent(in) :: ii
    integer, intent(in) :: jj
    integer, intent(in), optional  :: en_index
    integer, intent(in), optional  :: k_index
    integer, intent(in), optional  :: spin

    integer :: n, nbl, nrow, pl_start, pl_end

    nbl = this%struct%num_PLs
    if (this%scba_iter .eq. 0) return
    if (ii .ne. jj .or. ii.gt.nbl) return

    pl_start=this%struct%mat_PL_start(ii)
    pl_end=this%struct%mat_PL_end(ii)
    nrow = pl_end - pl_start + 1
    if (.not.allocated(sigma_n%val)) then
      call create(sigma_n, nrow, nrow)
    end if
    sigma_n%val = (0.0_dp, 0.0_dp)
    do n = 1, nrow
       sigma_n%val(n,n) = this%sigma_n(pl_start+n-1)
    end do

  end subroutine get_sigma_n_mat

  !> Give the Gr at given energy point to the interaction
  subroutine set_Gr(this, Gr, en_index, k_index, spin)
    class(ElPhonDephD) :: this
    type(z_dns), dimension(:,:), intent(in) :: Gr
    integer, intent(in), optional  :: en_index
    integer, intent(in), optional  :: k_index
    integer, intent(in), optional  :: spin

    integer :: n, npl, ii, pl_start, pl_end

    ! Do not update if maximum iterations have been reached
    if (this%scba_iter > this%scba_niter) return
    ! This model does not keep track of Gr at different energies because
    ! the model is local in energy. We directly calculate sigma_r
    npl = this%struct%num_PLs
    do n=1,npl
      pl_start=this%struct%mat_PL_start(n)
      pl_end=this%struct%mat_PL_end(n)
      do ii = 1, pl_end - pl_start + 1
          this%sigma_r(pl_start + ii - 1) = Gr(n,n)%val(ii,ii) * &
            & this%coupling(pl_start + ii - 1)
      end do
    end do

  end subroutine set_Gr

  !> Give the Gn at given energy point to the interaction
  subroutine set_Gn(this, Gn, en_index, k_index, spin)
    class(ElPhonDephD) :: this
    type(z_dns), dimension(:,:), intent(in) :: Gn
    integer, intent(in), optional  :: en_index
    integer, intent(in), optional  :: k_index
    integer, intent(in), optional  :: spin

    integer :: n, npl, ii, pl_start, pl_end

    ! Do not update if maximum iterations have been reached
    if (this%scba_iter > this%scba_niter) return
    ! This model does not keep track of Gr at different energies because
    ! the model is local in energy. We directly calculate sigma_r
    npl = this%struct%num_PLs
    do n=1,npl
      pl_start=this%struct%mat_PL_start(n)
      pl_end=this%struct%mat_PL_end(n)
      do ii = 1, pl_end - pl_start + 1
          this%sigma_n(pl_start + ii - 1) = Gn(n,n)%val(ii,ii) * &
            & this%coupling(pl_start + ii - 1)
      end do
    end do

  end subroutine set_Gn

  !>  Compute Sigma_r : dummy
  subroutine compute_Sigma_r(this, en_index, k_index, spin)
    class(ElPhonDephD) :: this
    integer, intent(in), optional  :: en_index
    integer, intent(in), optional  :: k_index
    integer, intent(in), optional  :: spin
  end subroutine compute_Sigma_r

  !>  Compute Sigma_n : dummy
  subroutine compute_Sigma_n(this, en_index, k_index, spin)
    class(ElPhonDephD) :: this
    integer, intent(in), optional  :: en_index
    integer, intent(in), optional  :: k_index
    integer, intent(in), optional  :: spin
  end subroutine compute_Sigma_n

  !> Destroy Sigma_r
  subroutine destroy_Sigma_r(this)
    class(ElPhonDephD) :: this
    call log_deallocate(this%sigma_r)      
  end subroutine destroy_Sigma_r

  !> Destroy Sigma_n
  subroutine destroy_Sigma_n(this)
    class(ElPhonDephD) :: this
    call log_deallocate(this%sigma_n)      
  end subroutine destroy_Sigma_n

end module elphdd
