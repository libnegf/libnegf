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


module elph

  use ln_precision, only : dp
  use globals
  use ln_allocation
  use mat_def, only : create, destroy, z_DNS

  implicit none
  private

  public :: Telph
  public :: init_elph_1, destroy_elph, init_elph_2


  !> This type contains information describing different electron phonon 
  !! models: input parameters, temporary data and output quantities
  !! 
  !! Note: I don't use explicitely interfaces and different data types
  !! for different models because then we'll need all the different
  !! containers in negf container (cannot use virtualization)
  !! It may be managed in a more elegant way
  type Telph
    !> Describe the model implemented. Currently supported:
    !! 0 : dummy model, no electron-phonon interactions
    !! 1 : electron phonon dephasing limit (as in Datta, Cresti etc.)
    !!     Assumes elastic scattering and fully local (diagonal) model
    !!     Coupling is diagonal (Local Deformation Potential) and 
    !!     Self energies are diagonal as well
    !! 2 : semi-local electron-phonon dephasing
    !!     Similar to 1, but the oscillator is considered local on
    !!     more than a contiguos basis function per oscillator position Ri
    !!     It is a atom-block generalization of 1 for LCAO
    !!     Coupling is diagonal (per oscillator site, per orbital)
    !!     Self energy is an array of atomic block
    !!     An additional descriptor with the number of orbitals per atom is 
    !!     needed.
    !!     Note: the modes do not need to be local on a single atom, but
    !!     you need the orbitals on a given local phonon site to be contiguous
    
    !! Common parameters
    integer :: model = 0
    !> SCBA option: number of iterations
    !! 0 corresponds to no iterations (self energy is not calculated)
    integer :: scba_niter = 0
    
    !> Keep track of SCBA iteration 
    integer :: scba_iter = 0
    !! Model specific
    !! -----------------------------------------------------------------------
    !! Model 1
    !!
    !> Diagonal coupling. Used in local coupling models (1)
    !! Note: it is stored directly as squared value as we always use it  
    !! that way (units energy^2)
    real(dp), allocatable, dimension(:) :: coupling_array
    !> Diagonal elelents of retarded self energy. Used only in model (1)
    complex(dp), allocatable, dimension(:) :: diag_sigma_r
    !> Diagonal elelents of lesser (n) self energy. Used only in model (1)
    complex(dp), allocatable, dimension(:) :: diag_sigma_n
    !! -----------------------------------------------------------------------
    !! Model 2
    !! -----------------------------------------------------------------------
    !> Block self energies, used in model (2)
    type(z_DNS), allocatable, dimension(:) :: atmblk_sigma_n
    type(z_DNS), allocatable, dimension(:) :: atmblk_sigma_r
    !> for model 2 we put the coupling in an atom block array, to avoid
    !  conversion all the time
    type(z_DNS), allocatable, dimension(:) :: atmcoupling
    !> An array specifying how many orbital per atom (assumed contiguous)
    !  used in model (2)
    integer, allocatable, dimension(:) :: orbsperatm
    !> From orbsperatom, aa work array containing starting orbital index
    !  for each atom, to speed up some patchworking
    integer, allocatable, dimension(:) :: atmorbstart
    !> For each atom, determine in which PL it sits, to accelerate block-sparse
    !  assignment
    integer, allocatable, dimension(:) :: atmpl
    !!-------------------------------------------------------------------------

    !> Number of active modes
    integer :: nummodes
     


    integer :: numselmodes
    logical, dimension(:), pointer :: selmodes => null()
    real(dp), dimension(:), pointer :: Wq => null()
    real(dp), dimension(:), pointer :: Nq => null()
 
    real(dp), dimension(:), pointer :: Mq => null()

    real(dp) :: self_range(2) !Energy interval for Sigma< 
    real(dp) :: Erange(2)     !Integration interval

    integer :: scba_iterations
    logical :: Selfene_Gr
    logical :: Selfene_Gless
    logical :: Selfene_Hilb
    logical :: memory
    logical :: check
  end type Telph

contains

  !>
  ! Initialize the el-ph structure when model = 1 (elastic model)
  ! @param elph: electron-phonon container
  ! @param coupling: coupling (energy units) 
  ! @param niter: fixed number of scba iterations
  subroutine init_elph_1(elph, coupling, niter)
    Type(Telph), intent(inout) :: elph
    real(dp), dimension(:), allocatable, intent(in) :: coupling
    integer :: niter

    elph%model = 1
    elph%coupling_array = coupling * coupling
    elph%scba_niter = niter
    call log_allocate(elph%diag_sigma_r, size(coupling))
    call log_allocate(elph%diag_sigma_n, size(coupling))
    elph%diag_sigma_r = 0.d0
    elph%diag_sigma_n = 0.d0
    elph%nummodes = 1  !Single 0eV mode (Actually n localized modes, but we 
                       !treat them all contemporary)

  end subroutine init_elph_1

  !>
  ! Initialize the el-ph structure when model = 2 (elastic quasi-local model)
  ! @param elph: electron-phonon container
  ! @param coupling: coupling per orbital (energy units) 
  ! @param orbsperatm: number of orbitals per each atom
  ! @param niter: fixed number of scba iterations
  ! @param pl_start: PL partitioning: used to accelerate block-sparse assignment
  !                  note: needs to contain NPL+1 elements, the last one 
  !                  is norbs+1, as in dftb
  subroutine init_elph_2(elph, coupling, orbsperatm, niter, pl_start)
    Type(Telph), intent(inout) :: elph
    real(dp), dimension(:), allocatable, intent(in) :: coupling
    integer, dimension(:), allocatable, intent(in) :: orbsperatm
    integer, dimension(:), pointer, intent(in) :: pl_start
    integer :: niter, ii, jj, natm, ierr

    !Check input size
    if (size(coupling).ne.sum(orbsperatm)) then
      stop 'Error: coupling and orbsperatom not compatible'
    end if
    elph%model = 2
    elph%scba_niter = niter
    elph%orbsperatm = orbsperatm
    natm = size(orbsperatm)
    call log_allocate(elph%atmorbstart, natm)
    call log_allocate(elph%atmpl, natm)
    elph%atmorbstart(1) = 1
    do ii = 2,natm
      elph%atmorbstart(ii) = sum(elph%orbsperatm(1:ii-1)) + 1
    enddo
    allocate(elph%atmblk_sigma_r(natm),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate atmblk_sigma_r'
    allocate(elph%atmblk_sigma_n(natm),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate atmblk_sigma_n'
    allocate(elph%atmcoupling(natm),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate atmcoupling'
    do ii = 1,natm
      call create(elph%atmblk_sigma_r(ii),orbsperatm(ii),orbsperatm(ii))
      elph%atmblk_sigma_r(ii)%val = 0.d0
      call create(elph%atmblk_sigma_n(ii),orbsperatm(ii),orbsperatm(ii))
      elph%atmblk_sigma_n(ii)%val = 0.d0
      call create(elph%atmcoupling(ii),orbsperatm(ii),orbsperatm(ii))
      elph%atmcoupling(ii)%val = 0.d0
    end do
    elph%nummodes = 1  !Single 0eV mode (Actually n localized modes, but we 
                       !treat them all contemporary)
    ! Assign coupling
    do ii = 1,natm
      do jj = 1,elph%orbsperatm(ii)
      elph%atmcoupling(ii)%val(jj,jj) = coupling(jj + elph%atmorbstart(jj) - 1)
      end do
    end do
    ! Determine atmpl
    elph%atmpl = 0
    do ii = 1,natm
      do jj = 1, size(pl_start) - 1
        if (elph%atmorbstart(ii).ge.pl_start(jj).and. &
            elph%atmorbstart(ii).lt.pl_start(jj + 1)) then
          elph%atmpl(ii) = jj
        end if
      end do
    end do
    ! Check that they are all assigned
    do ii = 1,natm
      if (elph%atmpl(ii).eq.0) then
        write(*,*) elph%atmpl
        stop 'atmpl not correctly set'
      end if
    end do

  end subroutine init_elph_2


  !>
  ! Destroy elph structure when model = 1 (elastic model)
  subroutine destroy_elph_1(elph)
    Type(Telph) :: elph

    elph%model = 0
    call log_deallocate(elph%coupling_array)
    call log_deallocate(elph%diag_sigma_r)
    call log_deallocate(elph%diag_sigma_n)

  end subroutine destroy_elph_1

  !>
  ! Destroy elph structure when model = 1 (elastic model)
  subroutine destroy_elph_2(elph)
    Type(Telph) :: elph

    integer :: ii, ierr

    elph%model = 0
    do ii=1,size(elph%orbsperatm)
      call destroy(elph%atmblk_sigma_r(ii))
      call destroy(elph%atmblk_sigma_n(ii))
      call destroy(elph%atmcoupling(ii))
    end do
    deallocate(elph%atmcoupling, stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not deallocate atmblk_sigma_r'
    deallocate(elph%atmblk_sigma_n, stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not deallocate atmblk_sigma_n'
    call log_deallocate(elph%orbsperatm)
    call log_deallocate(elph%atmorbstart)
    call log_deallocate(elph%atmpl)

  end subroutine destroy_elph_2


  !>
  ! el-ph destruction interface
  subroutine destroy_elph(elph)
     Type(Telph) :: elph

     if (elph%model .eq. 0) then
       return
     else if (elph%model .eq. 1) then
       call destroy_elph_1(elph)
     else if (elph%model .eq. 2) then
       call destroy_elph_2(elph)
     else
       write(*,*) 'Warning, not implemented'
     end if

  end subroutine destroy_elph



  subroutine init_elph(elph,nummodes)
    Type(Telph) :: elph
    integer, optional :: nummodes

    if (present(nummodes)) then
       elph%nummodes = nummodes
       elph%numselmodes = nummodes
    else
       elph%nummodes = 0
       elph%numselmodes = 0
    end if
    elph%scba_iterations = 0 ! starts from 0  
    elph%scba_iter = 0       ! initialize at 0

    if (elph%nummodes .eq. 0) then
      
      elph%Selfene_Gr = .false.
      elph%Selfene_Gless =.false.
      elph%Selfene_Hilb = .false.
     
    else

      call log_allocatep(elph%selmodes, elph%nummodes)
      call log_allocatep(elph%Wq, elph%nummodes)
      call log_allocatep(elph%Nq, elph%nummodes)
      call log_allocatep(elph%Mq, elph%nummodes)
 
      elph%Wq = 0.0_dp/27.2114_dp  !0.050_dp  ! 50 meV phonon
      elph%Nq = 0.0_dp              ! should be bose-einstain
      elph%Mq = 0.0_dp/27.2114_dp    ! 100 meV phonon coupling
      elph%selmodes = .true.
 
      elph%Selfene_Gr = .true.
      elph%Selfene_Gless = .true.
      elph%Selfene_Hilb = .true.
      
      elph%memory = .true.
      elph%check = .false. 
    end if

  end subroutine init_elph



end module elph


!*******
!*******

!subroutine pippo(n)
!   integer, intent(in) :: n
!
!   integer, dimension(:), allocatable :: A
!   real(dp), dimension(:), pointer :: P
!   double DEPRECATED
!
!   allocate(A(n))

!   call log_allocate(A,n)     ! safe allocation in libNEGF

!   operazioni su A
!   call sub(A)

!   call log_deallocate(A)

!end subroutine pippo


!subroutine sub(B)
!   integer, dimension(:) :: B

!   integer :: i

!   do i = 1, size(B)
!      B(i) = 5
!   end do

! end subroutine sub


!subroutine sub2(B,n)
!   integer :: B(*)
!   integer :: n

!   integer :: i

!   do i = 1, n
!      B(i) = 5
!   end do

! end subroutine sub2

