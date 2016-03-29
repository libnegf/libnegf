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
  use mat_def, only : create, destroy, z_DNS, z_CSR, z_COO
  use sparsekit_drv, only : prealloc_mult, prealloc_sum, coo2csr, extract, csr2dns

  implicit none
  private

  public :: Telph
  public :: init_elph_1, destroy_elph, init_elph_2, init_elph_3


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
    !! 3 : as 2, but I will set the coupling as csr matrix including 
    !!     overlap M->MS/2+SM/2 and treat each atomic oscillator mode separately
    !!     This is needed to verify whether neglecting overlap is ok
    
    !! Common parameters
    integer :: model = 0
    !> SCBA option: number of iterations
    !! 0 corresponds to no iterations (self energy is not calculated)
    integer :: scba_niter = 0
    
    !> Keep track of SCBA iteration 
    integer :: scba_iter = 0

    !> SCBA Tolerance (Exact meaning may depend on model)
    real(dp) :: scba_tol = 1.0d-7

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

    !> Generic CSR self energies per each mode, now used in model (3)
    type(z_CSR), allocatable, dimension(:) :: csr_sigma_n
    type(z_CSR), allocatable, dimension(:) :: csr_sigma_r
    !> for model 3, CSR couplings
    type(z_CSR), allocatable, dimension(:) :: csr_couplings


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
    !  assignment. Used in model 2,3
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
    elph%scba_iter = 0
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
      elph%atmcoupling(ii)%val(jj,jj) = coupling(jj + elph%atmorbstart(ii) - 1)
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
  ! Initialize the el-ph structure when model = 3 (semilocal S masked)
  ! @param elph: electron-phonon container
  ! @param coupling: coupling per orbital (energy units) 
  ! @param orbsperatm: number of orbitals per each atom
  ! @param niter: fixed number of scba iterations
  ! @param pl_start: PL partitioning: used to accelerate block-sparse assignment
  !                  note: needs to contain NPL+1 elements, the last one 
  !                  is norbs+1, as in dftb
  subroutine init_elph_3(elph, coupling, orbsperatm, niter, pl_start, over)
    Type(Telph), intent(inout) :: elph
    real(dp), dimension(:), allocatable, intent(in) :: coupling
    integer, dimension(:), allocatable, intent(in) :: orbsperatm
    integer, dimension(:), pointer, intent(in) :: pl_start
    type(z_CSR), pointer, intent(in) :: over

    integer :: niter, ii, jj, natm, ierr, norbs, offset, iimode
    type(z_COO) :: tmp_coo
    type(z_CSR) :: work1, work2, work3, over_device

    !Check input size
    if (size(coupling).ne.sum(orbsperatm)) then
      stop 'Error: coupling and orbsperatom not compatible'
    end if
    elph%model = 3
    elph%scba_niter = niter
    elph%scba_iter = 0
    elph%orbsperatm = orbsperatm
    natm = size(orbsperatm)
    norbs = size(coupling)
    elph%nummodes = natm ! One localized oscillator per atom

    call extract(over, 1, norbs, 1, norbs, over_device)

    !! Check how many modes we really need, 
    !! in some cases we can have all zero couplings
    do ii = 1,natm
      offset = sum(orbsperatm(1:ii-1))
      !! if all couplings are zero, just skip it.
      if (all(coupling(offset+1:offset+orbsperatm(ii)) .eq. 0.0d0)) then
        elph%nummodes = elph%nummodes - 1
      end if
    end do

    allocate(elph%csr_sigma_r(elph%nummodes), stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate csr_sigma_r'
    allocate(elph%csr_sigma_n(elph%nummodes), stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate csr_sigma_n'
    allocate(elph%csr_couplings(elph%nummodes), stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate csr_couplings'
    
    iimode = 0
    do ii = 1, natm
     offset = sum(orbsperatm(1:ii-1))
      !! if all couplings are zero, just skip it.
      if (all(coupling(offset+1:offset+orbsperatm(ii)) .eq. 0.0d0)) then
        cycle
      end if
      iimode = iimode + 1
      !! I assemble directly the very sparse CSR
      call create(work1, norbs, norbs, orbsperatm(ii))

         !   write(*,*) "orbsperatom", orbsperatm(ii), ii, offset, norbs, natm
      work1%rowpnt(1:offset) = 1
      work1%rowpnt(offset+1:norbs+1) = orbsperatm(ii) + 1
      do jj = 1,orbsperatm(ii)
        work1%nzval(jj) = coupling(jj + offset)
        work1%rowpnt(jj + offset) = jj
        work1%colind = jj + offset
      end do
      !! M -> M*S/2 + S*M/2
      call prealloc_mult(work1, over_device, (0.5d0, 0.0d0), work2)
      call prealloc_mult(over_device, work1, (0.5d0, 0.0d0), work3)
      call destroy(work1)
      call prealloc_sum(work2, work3, elph%csr_couplings(iimode))
      call destroy(work2)
      call destroy(work3)

    end do
    call destroy(over_device)

    ! Determine atmpl
!!$    elph%atmpl = 0
!!$    do ii = 1,natm
!!$      do jj = 1, size(pl_start) - 1
!!$        if (elph%atmorbstart(ii).ge.pl_start(jj).and. &
!!$            elph%atmorbstart(ii).lt.pl_start(jj + 1)) then
!!$          elph%atmpl(ii) = jj
!!$        end if
!!$      end do
!!$    end do
!!$    ! Check that they are all assigned
!!$    do ii = 1,natm
!!$      if (elph%atmpl(ii).eq.0) then
!!$        write(*,*) elph%atmpl
!!$        stop 'atmpl not correctly set'
!!$      end if
!!$    end do

    end subroutine init_elph_3


  !>
  ! Destroy elph structure when model = 1 (elastic model)
  subroutine destroy_elph_3(elph)
    Type(Telph) :: elph

    integer :: ii, ierr

    elph%model = 0
! These should be cleaned up after every energy point
    deallocate(elph%csr_couplings, stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not deallocate csr_couplings'
    deallocate(elph%csr_sigma_r, stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not deallocate csr_sigma_r'
    deallocate(elph%csr_sigma_n, stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not deallocate csr_sigma_n'
    call log_deallocate(elph%orbsperatm)

  end subroutine destroy_elph_3


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
  ! Destroy elph structure when model = 2 (block diagonal elastic model)
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
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not deallocate atmcoupling'
    deallocate(elph%atmblk_sigma_n, stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not deallocate atmblk_sigma_n'
    deallocate(elph%atmblk_sigma_r, stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not deallocate atmblk_sigma_r'
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
     else if (elph%model .eq. 3) then
       call destroy_elph_3(elph)
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
