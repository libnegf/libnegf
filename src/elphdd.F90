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
  use interactions, only : interaction
  use ln_allocation, only : log_allocate, log_deallocate
  use ln_structure, only : TStruct_info
  use mat_def, only : z_dns, create
  
  implicit none
  private

  public :: ElPhonDephD, ElPhonDephD_create

  type, extends(interaction) :: ElPhonDephD

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
    procedure :: get_sigma_n
    procedure :: set_Gr
    procedure :: set_Gn

  end type ElPhonDephD

contains

  !>
  ! Factory for el-ph dephasing diagonal model
  ! @param struct : contact/device partitioning infos
  ! @param coupling: coupling (energy units) 
  ! @param niter: fixed number of scba iterations
  ! @param tol: scba tolerance
  subroutine ElPhonDephD_create(this, struct, coupling, niter, tol)
    
    type(ElPhonDephD), intent(inout) :: this
    type(TStruct_info), intent(in) :: struct
    real(dp), dimension(:), intent(in) :: coupling
    integer, intent(in) :: niter
    real(dp), intent(in) :: tol

    this%descriptor = &
        & "Electron-Phonon dephasing model in fully diagonal model"

    call log_allocate(this%coupling, size(coupling))
    call log_allocate(this%sigma_r, size(coupling))
    call log_allocate(this%sigma_n, size(coupling))
    this%struct = struct
    this%coupling = coupling * coupling
    this%scba_niter = niter
    this%scba_tol = tol
    this%sigma_r = (0.0d0, 0.0d0)
    this%sigma_n = (0.0d0, 0.0d0)
    this%nummodes = 1  !Single 0eV mode (Actually n localized modes, but we 
                       !treat them all contemporary)

  end subroutine ElPhonDephD_create


  !> This interface should append
  !  the retarded self energy to ESH
  subroutine add_sigma_r(this, esh)
    class(ElPhonDephD) :: this
    type(z_dns), dimension(:,:), intent(inout) :: esh

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
  

  !> Returns the lesser (n) Self Energy in block format
  !  
  subroutine get_sigma_n(this, blk_sigma_n, en_index)
    class(ElPhonDephD) :: this
    type(z_dns), dimension(:,:), intent(inout) :: blk_sigma_n
    integer, intent(in) :: en_index

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

  end subroutine get_sigma_n

  !> Give the Gr at given energy point to the interaction
  subroutine set_Gr(this, Gr, en_index)
    class(ElPhonDephD) :: this
    type(z_dns), dimension(:,:), intent(in) :: Gr
    integer :: en_index

    integer :: n, npl, ii, pl_start, pl_end
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
  subroutine set_Gn(this, Gn, en_index)
    class(ElPhonDephD) :: this
    type(z_dns), dimension(:,:), intent(in) :: Gn
    integer :: en_index

    integer :: n, npl, ii, pl_start, pl_end
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


end module elphdd
