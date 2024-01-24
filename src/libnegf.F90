!! libNEGF: a general library for Non-Equilibrium Green's functions.        !
!! Copyright (C) 2012                                                       !
!!                                                                          !
!! This file is part of libNEGF: a library for                              !
!! Non Equilibrium Green's Function calculation                             !
!!                                                                          !
!! Developers: Alessandro Pecchia, Gabriele Penazzi                         !
!! Former Contributors: Luca Latessa, Aldo Di Carlo                        !
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


module libnegf

 use ln_precision
 use ln_constants
 use ln_allocation
 use lib_param
 use ln_cache
 use globals, only : LST
 use mpi_globals, only : id, id0, numprocs, negf_cart_init, check_cart_comm, &
      & globals_mpi_init => negf_mpi_init
 use input_output
 use ln_structure
 use rcm_module
 use mat_def
 use ln_extract
 use sparsekit_drv
 use integrations
 use iso_c_binding
 use system_calls
 use elph, only : interaction_models
 use interactions, only : get_max_wq
#:if defined("MPI")
 use libmpifx_module, only : mpifx_comm
#:endif
 use clock
 implicit none
 private

 public :: log_deallocatep
 public :: r_CSR, z_CSR, r_DNS, z_DNS, create, destroy   !from matdef
 public :: HAR, eovh, pi, kb, units, set_drop, DELTA_SQ, DELTA_W, DELTA_MINGO ! from ln_constants
 public :: convertCurrent, convertHeatCurrent, convertHeatConductance ! from ln_constants
 public :: Tnegf
 public :: set_bp_dephasing
 public :: set_elph_dephasing, set_elph_block_dephasing, set_elph_s_dephasing
 public :: set_elph_polaroptical, set_elph_nonpolaroptical, destroy_interactions
 public :: interaction_models
 public :: set_clock, write_clock
 public :: writeMemInfo, writePeakInfo
 public :: dns2csr, csr2dns, nzdrop

 public :: id, id0
#:if defined("MPI")
 public :: negf_mpi_init
 public :: negf_cart_init !from mpi_globals
 public :: set_cartesian_bare_comms
#:endif
 public :: set_kpoints
 public :: set_mpi_bare_comm

 !Input and work flow procedures
 public :: lnParams
 public :: init_negf, destroy_negf
 public :: init_contacts, init_structure, init_basis
 public :: get_params, set_params, set_scratch, set_outpath, create_scratch, set_scba_tolerances
 public :: init_ldos, set_ldos_intervals, set_ldos_indexes, set_tun_indexes

 public :: create_HS, destroy_HS, set_H, set_S, set_S_id, read_HS, pass_HS, copy_HS
 public :: set_readOldDMsgf, set_readOldTsgf, set_computation
 public :: set_convfactor, set_fictcont
 public :: read_negf_in
 public :: negf_version
 public :: destroy_contact_matrices ! cleanup matrices in Tnegf container (H,S)
 public :: destroy_surface_green_cache ! Clean surface green cache (useful for memory cache)
 public :: destroy_DM ! cleanup matrices in Tnegf container (rho,rhoE)
 private :: block_partition ! chop structure into PLs (CAREFUL!!!)
                            ! H need to be already ordered properly
 public :: negf_partition_info  !write down partition info
 public :: find_cblocks        ! Find interacting contact block
 public :: set_ref_cont, print_tnegf
 public :: associate_transmission, associate_current, associate_ldos
 public :: associate_lead_currents
 public :: get_energies, pass_DM, get_DM, get_currents

 public :: compute_density_dft      ! high-level wrapping
                                    ! Extract HM and SM
                                    ! run DM calculation
 public :: compute_density_efa      ! high-level wrapping
                                    ! Extract HM and SM
                                    ! run DM calculation
 public :: compute_density_quasiEq  ! Run DM calculation with
                                    ! quasi-equilibr. approximation
 public :: compute_current          ! high-level wrapping routines
                                    ! Extract HM and SM
                                    ! run total current calculation
 public :: compute_meir_wingreen    ! Meir-Wingreen contact currents
 public :: compute_layer_current    ! high-level wrapping routines
                                    ! layer currents calculation
 public :: compute_dephasing_transmission ! high-level wrapping routines
                                          ! Extract HM and SM
                                          ! run total current calculation
 public :: write_tunneling_and_dos  ! Print tunneling and dot to file
                                    ! Note: for debug purpose. I/O should be managed
                                    ! by calling program
 public :: compute_ldos             ! wrapping to compute ldos
 public :: return_dos_mat           ! return pointer to dos_proj matrix

 public :: compute_phonon_current   ! High-level wrapping to
                                    ! compute phonon transmission
                                    ! and heat currents
 public :: thermal_conductance

 public :: reorder, sort, swap            ! not used
 public :: printcsr   ! debugging routines
 public :: printcsrij   ! debugging routines
 public :: getel   ! debugging routines

 integer, parameter :: VBT=70
 integer, parameter :: MAXNUMPLs = 10000
 integer, parameter, public :: READ_SGF = 0
 integer, parameter, public :: COMP_SGF = 1
 integer, parameter, public :: COMPSAVE_SGF = 2
  !-----------------------------------------------------------------------------
  !> Contains all the general parameters to be passed as input to library
  !! which are compatible with iso_c_binding representations
  !! 1-1 correspondance to input parameter in type Tnegf
  type, bind(c) :: lnparams
   !! General
   !> verbosity, > 100 is maximum
   integer(c_int) :: verbose
   !> Managing SGF readwrite for DM: 0: Read 1: compute 2: comp & save
   integer(c_int)  :: readOldDM_SGFs
   !> Managing SGF readwrite for Tunn: 0: Read 1: compute 2: comp & save
   integer(c_int)  :: readOldT_SGFs
   !> SGF cache destination: 0 for disk, 1 for memory, 2 for a dummy cache (no save performed)
   integer(c_int) :: SGFcache
   !> Spin component (for io)
   integer(c_int)  :: spin
   !> k-point index (for io)
   integer(c_int) :: ikpoint
   !> Spin degeneracy
   real(c_double) :: g_spin
   !> Imaginary delta
   real(c_double) :: delta
   !> delta model for phonon GF
   integer(c_int) :: deltaModel
   !> Maximal energy in Mingo delta model
   real(c_double) :: wmax
   !> Additional delta to force more broadening in the DOS
   real(c_double) :: dos_delta
   !> Energy conversion factor
   real(c_double) :: eneconv
   !> Weight for k-point integration, output integrals are scaled accordingly
   real(c_double) :: kwght
   !> Conduction band edge
   real(c_double) :: ec
   !> Valence band edge
   real(c_double) :: ev
   !> Safe guard energy below Ec
   real(c_double) :: deltaec
   !> Safe guard energy above Ev
   real(c_double) :: deltaev
   !! Real axis integral
   !> Minimum energy for real axis (current integration, DOS, tunneling)
   real(c_double) :: emin
   !> Maximum energy for real axis
   real(c_double) :: emax
   !> Energy step for real axis
   real(c_double) :: estep
   !> Energy step for coarse integrations
   real(c_double) :: estep_coarse
   !! Contacts info
   !> Electron electrochemical potential
   real(c_double) :: mu_n(MAXNCONT)
   !> Hole electrochemical potential
   real(c_double) :: mu_p(MAXNCONT)
   !> Electron electrochemical Potential (for dft)
   real(c_double) :: mu(MAXNCONT)
   !> Contact DOS for WBA
   real(c_double) :: contact_dos(MAXNCONT)
   !> Logical value: is the contact WB?
   logical(c_bool)  :: fictcont(MAXNCONT)
   !> Electronic temperature for each contact (Density Matrix)
   real(c_double) :: kbT_dm(MAXNCONT)
   !> Electronic temperature for each contact (Transmission)
   real(c_double) :: kbT_t(MAXNCONT)
   !> SCBA tolerance for inelastic loop
   real(c_double) :: scba_inelastic_tol
   !> SCBA tolerance for elastic loop
   real(c_double) :: scba_elastic_tol
   !! Contour integral
   !> Number of points for n
   integer(c_int) :: np_n(2)
   !> Number of points for p
   integer(c_int) :: np_p(2)
   !> Number of real axis points
   integer(c_int) :: np_real
   !> ! Number of kT extending integrations
   integer(c_int) :: n_kt
   !> Number of poles
   integer(c_int) :: n_poles
   !! Emitter and collector for transmission or Meir-Wingreen
   !! (only emitter in this case)
   !> Emitter contact list (or lead for integration in MW)
   integer(c_int) :: ni(MAXNCONT)
   !> Collector contact list
   integer(c_int) :: nf(MAXNCONT)
   !> Should I calculate the density ("D") or the energy weighted density ("E")?
   character(kind=c_char, len=1) :: dore  ! Density or En.Density
   !> Reference contact is set to maximum or minimum Fermi level
   integer(c_int) :: min_or_max
 end type lnparams
  !-----------------------------------------------------------------------------

contains

  !--------------------------------------------------------------------
  !>  Init libNEGF
  !!  General initializations of libNEGF are currently done via files.
  !!  "negf.in" contains all parameters information
  !!  all needed info about the structure and matrix format.
  !!  H and S are also read-in from files
  !!  Some parameters are still hard-coded and need to be set from api
  !--------------------------------------------------------------------
  subroutine init_negf(negf)
    type(Tnegf) :: negf

    real(dp) :: kpoints(3,1), kweights(1)
    integer :: local_kindex(1)

    call set_defaults(negf)
    negf%form%formatted = .true.
    negf%form%type = "PETSc"
    negf%form%fmt = "F"

    ! Allocate zero contacts by default. The actual number of contacts
    ! can be set calling init_contacts again.
    call init_contacts(negf, 0)

    ! Initialize a default gamma point
    kpoints(:,1) = (/0.0_dp, 0.0_dp, 0.0_dp /)
    kweights(1) = 1.0_dp
    local_kindex(1) = 1
    call set_kpoints(negf, kpoints, kweights, local_kindex, 0)

  end subroutine init_negf

   !!===================================================================
   !! INPUT Routines
   !!===================================================================

   !!-------------------------------------------------------------------
   !! Passing & reading H,S
   !!-------------------------------------------------------------------

   !>--------------------------------------------------------------------
   !!  Read H and S from file
   !!  @param[in] negf: libnegf container instance
   !!  @param[in] real_path (string): filename for real part of target matrix in PETSC format
   !!  @param[in] imag_path (string): filename for imaginary part of target matrix in PETSC format.
   !!  @param[in] target_matrix (integer): controlo flag which specify if the Hamiltonian
   !!      or the Overlap is parsed. 0 for the Hamiltonian, 1 for the Overlap
   !!      formatted (logical, optional): true (default) if file is formatted.
   !!
   !!   Note: up to now both the real and imaginary part files must have the same
   !!     indexes and number of non zero elements, even if zero values appear
   subroutine read_HS(negf, real_path, imag_path, target_matrix, formatted)
     type(Tnegf), intent(inout) :: negf
     character(len=*), intent(in) :: real_path, imag_path
     integer, intent(in) :: target_matrix
     logical, intent(in), optional :: formatted

     logical :: formatted_opt
     logical :: doesexist
     character(11) :: fmtstring
     type(format) :: fmt

     if (present(formatted)) then
       formatted_opt = formatted
     else
       formatted_opt = .true.
     endif

     if(formatted_opt) then
       fmtstring = 'formatted'
       fmt%formatted = .true.
     else
       fmtstring = 'unformatted'
       fmt%formatted = .false.
     endif

     fmt%type = 'PETSc'  !! only PETSc now supported here. Maybe need to add support for UPT?
     fmt%fmt = 'F'       !! only Full matrix supported. We could also support upper or lower

     inquire(file=trim(real_path), exist= doesexist)
     inquire(file=trim(imag_path), exist= doesexist)
     if (.not.doesexist) then
       write(*,*) "libNEGF error. Matrix files not found"
       stop
     endif

     open(401, file=real_path, form=trim(fmtstring))
     open(402, file=imag_path, form=trim(fmtstring))

     call create_HS(negf, 1)
     if (target_matrix.eq.0) then
       if (.not.associated(negf%HS(1)%H)) allocate(negf%HS(1)%H)
       call read_H(401,402,negf%HS(1)%H,fmt)
       negf%H => negf%HS(1)%H
     else if (target_matrix.eq.1) then
       if (.not.associated(negf%HS(1)%S)) allocate(negf%HS(1)%S)
       call read_H(401,402,negf%HS(1)%S,fmt)
       negf%S => negf%HS(1)%S
     else
       write(*,*) "libNEGF error. Wrong target_matrix: must be 0 (H) or 1 (S)"
       stop
     endif
     close(401)
     close(402)

   end subroutine read_HS

  !-------------------------------------------------------------------
  !> Create the HS container for k-dependent H(k) and S(k)
  !
  subroutine create_HS(negf, nHS)
    type(Tnegf) :: negf
    integer, intent(in) :: nHS

    if (.not.allocated(negf%HS)) then
       !print*,'create HS container with',nHS,'  elements'
       allocate(negf%HS(nHS))
    else
       if (size(negf%HS) .ne. nHS) then
         print*, "creating HS container with size",nHS
         print*, "HS container already present with size",size(negf%HS)
         error stop 'ERROR: HS container already created with different size'
       end if
    end if

  end subroutine create_HS

  !-------------------------------------------------------------------
  !> Create the DM container for k-dependent rho(k) and rhoE(k)
  !
  !subroutine create_DM(negf, nDM)
  !  type(Tnegf) :: negf
  !  integer, intent(in) :: nDM
  !
  !  allocate(negf%DM(nDM))
  !end subroutine create_DM


  !-------------------------------------------------------------------
  !> Pass H from memory in CSR format
  !! @param[in] negf: libnegf container instance
  !! @param[in] nrow: number of rows
  !! @param[in] nzval: number of non zero values
  !! @param[in] colind: column indexes
  !! @param[in] rowpnt: row pointers
  !! @param[in] iKS: k-index of H
  subroutine set_H(negf, nrow, nzval, colind, rowpnt, iKS)
    type(Tnegf) :: negf
    integer :: nrow
    complex(dp) :: nzval(*)
    integer :: colind(*)
    integer :: rowpnt(*)
    integer, optional :: iKS

    integer :: nnz, i, base, ii

    if (present(iKS)) then
      ii = iKS
    else
      ii = 1
    end if

    if (ii > size(negf%HS)) then
       stop "Error: set_H with index > allocated array. Call create_HS with correct size"
    else
      if (.not.associated(negf%HS(ii)%H)) then 
        allocate(negf%HS(ii)%H)
      else
        call destroy(negf%HS(ii)%H)
      endif
    end if

    base = 0
    if (rowpnt(1) == 0) base = 1

    nnz = rowpnt(nrow+1)-rowpnt(1)

    call create(negf%HS(ii)%H,nrow,nrow,nnz)

    do i = 1, nnz
      negf%HS(ii)%H%nzval(i) = nzval(i)
      negf%HS(ii)%H%colind(i) = colind(i) + base
    enddo
    do i = 1,nrow+1
      negf%HS(ii)%H%rowpnt(i) = rowpnt(i) + base
    enddo
    negf%HS(ii)%internalHS=.true.
    negf%H => negf%HS(ii)%H

  end subroutine set_H

  !-------------------------------------------------------------------
  !> Pass S from memory in CSR format
  !! @param[in] negf: libnegf container instance
  !! @param[in] nrow: number of rows
  !! @param[in] nzval: number of non zero values
  !! @param[in] colind: column indexes
  !! @param[in] rowpnt: row pointers
  !! @param[in] iKS: k-index of S
  subroutine set_S(negf, nrow, nzval, colind, rowpnt, iKS)
    type(Tnegf) :: negf
    integer :: nrow
    complex(dp) :: nzval(*)
    integer :: colind(*)
    integer :: rowpnt(*)
    integer, optional :: iKS

    integer :: nnz, i, base, ii

    if (present(iKS)) then
      ii = iKS
    else
      ii = 1
    end if

    if (ii > size(negf%HS)) then
       stop "Error: set_S with index > allocated array. Call create_HS with correct size"
    else
      if (.not.associated(negf%HS(ii)%S)) then 
        allocate(negf%HS(ii)%S)
      else
        call destroy(negf%HS(ii)%S)
      endif
    end if

    base = 0
    if (rowpnt(1) == 0) base = 1

    nnz = rowpnt(nrow+1)-rowpnt(1)

    call create(negf%HS(ii)%S,nrow,nrow,nnz)

    do i = 1, nnz
      negf%HS(ii)%S%nzval(i) = nzval(i)
      negf%HS(ii)%S%colind(i) = colind(i) + base
    enddo
    do i = 1,nrow+1
      negf%HS(ii)%S%rowpnt(i) = rowpnt(i) + base
    enddo
    negf%HS(ii)%internalHS=.true.
    negf%S => negf%HS(ii)%S

  end subroutine set_S

  !-------------------------------------------------------------------
  !> Set S as identity
  !! @param[in] negf: libnegf container instance
  !! @param[in] nrow: number of rows
  !! @param[in] iKS: k-index of S
  subroutine set_S_id(negf, nrow, iKS)
    type(Tnegf) :: negf
    integer, intent(in) :: nrow
    integer, intent(in), optional :: iKS

    integer :: ii
    if (present(iKS)) then
      ii = iKS
    else
      ii = 1
    end if

    if (ii > size(negf%HS)) then
       stop "Error: set_S_id with index > allocated array. Call create_HS with correct size"
    else
      if (.not.associated(negf%HS(ii)%S)) allocate(negf%HS(ii)%S)
    end if

    call create_id(negf%HS(ii)%S, nrow)
    negf%HS(ii)%isSid = .true.
    negf%S => negf%HS(ii)%S

  end subroutine set_S_id

  !------------------------------------------------------------------
  !> Assign H,S pointers to externally allocated matrices
  !! @param [in]  negf: libnegf container instance
  !! @param [in] H: target z_CSR hamiltonian
  !! @param [in] S: target z_CSR overlap (optional, default to identity)
  !! @param[in] iKS: k-index of H and S
  subroutine pass_HS(negf,H,S,iKS)
    type(Tnegf) :: negf
    type(z_CSR), pointer, intent(in) :: H
    type(z_CSR), pointer, intent(in), optional :: S
    integer, intent(in), optional :: iKS

    integer :: ii
    if (present(iKS)) then
      ii = iKS
    else
      ii = 1
    end if

    if (ii > size(negf%HS)) then
      print*,'Passing HS index',ii,' HS container with size',size(negf%HS)
      stop "Error: pass_HS with index > allocated array. Call create_HS with correct size"
    endif

    negf%HS(ii)%H => H
    if (present(S)) then
       negf%HS(ii)%S => S
    else
       negf%HS(ii)%isSid=.true.
       allocate(negf%HS(ii)%S)
       call create_id(negf%HS(ii)%S,negf%HS(ii)%H%nrow)
    endif
    negf%HS(ii)%internalHS = .false.
    negf%H => negf%HS(ii)%H
    negf%S => negf%HS(ii)%S

  end subroutine pass_HS

  ! -----------------------------------------------------
  !  Allocate and copy H,S
  ! -----------------------------------------------------
  subroutine copy_HS(negf,H,S,iKS)
    type(Tnegf) :: negf
    type(z_CSR), intent(in) :: H
    type(z_CSR), intent(in), optional :: S
    integer, intent(in), optional :: iKS

    integer :: ii
    if (present(iKS)) then
      ii = iKS
    else
      ii = 1
    end if

    if (ii > size(negf%HS)) then
       stop "Error: copy_HS with index > allocated array. Call create_HS with correct size"
    else
      if (.not.associated(negf%HS(ii)%H)) allocate(negf%HS(ii)%H)
      if (.not.associated(negf%HS(ii)%S)) allocate(negf%HS(ii)%S)
    end if

    call create(negf%HS(ii)%H,H%nrow,H%ncol,H%nnz)
    negf%HS(ii)%H%nzval = H%nzval
    negf%HS(ii)%H%colind = H%colind
    negf%HS(ii)%H%rowpnt = H%rowpnt
    negf%HS(ii)%H%sorted = H%sorted

    if (present(S)) then
       negf%HS(ii)%isSid=.false.
       call create(negf%HS(ii)%S,S%nrow,S%ncol,S%nnz)
       negf%HS(ii)%S%nzval = S%nzval
       negf%HS(ii)%S%colind = S%colind
       negf%HS(ii)%S%rowpnt = S%rowpnt
       negf%HS(ii)%S%sorted = S%sorted
    else
       negf%HS(ii)%isSid=.true.
       call create_id(negf%HS(ii)%S,negf%HS(ii)%H%nrow)
    endif

    negf%HS(ii)%internalHS = .true.
    negf%H => negf%HS(ii)%H
    negf%S => negf%HS(ii)%S

  end subroutine copy_HS

  !!-------------------------------------------------------------------
  !! Setting structure and partitioning
  !!-------------------------------------------------------------------
  !------------------------------------------------------------------
  !> Set device/contact and PLs partition information
  !! @param [in] negf: libnegf container instance
  !! @param [in] ncont: number of contacts
  !! @param [in] contend: on which hamiltonian index where are contacts
  !!               ending? Array with size ncont
  !! @param [in] surfend: on which index is Device surface ending before
  !!              the corresponding contact (would be beginning contact-1)
  !! @param [in] npl: number of principal layers
  !! @param [in] plend: indexes where each principal layer ends
  !! @param [in] cblk: array with index of interacting blocks for each
  !!               contact(fortran indexing. If cblk is not known, use
  !!               find_cblocks
  !!
  !! If nbl = 0 the code will try to guess an automatic partitioning and
  !! plend, cblk will be ignored.
  !!
  !! Example: device goes from 1 to 60. Contacts from 61 to 80 and to
  !! 81 to 100. Only 1 PL:
  !!  ncont = 2
  !!  contend = [80, 100]
  !!  surfend = [60, 80]
  !!  npl = 1
  !!  plend = [60]
  subroutine init_structure(negf, ncont, surfstart, surfend, contend, npl, plend, cblk)
     type(Tnegf) :: negf
     integer, intent(in) :: ncont
     integer, intent(in) :: surfstart(:), surfend(:), contend(:)
     integer, intent(in) :: npl
     integer, intent(in) :: plend(:)
     integer, intent(inout), allocatable :: cblk(:)

     integer, allocatable :: plend_tmp(:), cblk_tmp(:)
     integer :: npl_tmp


     ! Make sure we called init_contacts in a consistent way.
     if (size(negf%cont) .ne. ncont) then
       write(*, *) 'size(negf%cont)=',size(negf%cont),'<->  ncont=',ncont
       stop "Error in set_structure: ncont not compatible with previous initialization."
     end if
     ! More sanity checks.
     if (size(surfstart) .ne. ncont) then
       write(*, *) 'size(surfstart)=',size(surfstart),'<->  ncont=',ncont
       stop "Error in set_structure: surfend and ncont mismatch"
     end if
     if (size(surfend) .ne. ncont) then
       write(*, *) 'size(surfend)=',size(surfend),'<->  ncont=',ncont
       stop "Error in set_structure: surfend and ncont mismatch"
     end if
     if (size(contend) .ne. ncont) then
       write(*, *) 'size(contend)=',size(contend),'<->  ncont=',ncont
       stop "Error in set_structure: contend and ncont mismatch"
     end if
     if (npl .ne. 0 .and. size(plend) .ne. npl) then
       write(*, *) 'size(plend)=',size(plend),'<->  npl=',npl
       stop "Error in set_structure: plend and npl mismatch"
     end if

     if (npl .eq. 0) then
       if (.not.allocated(negf%HS)) then
         stop "Error in init_structure: invoking block_partition but H not created"
         if (.not.associated(negf%HS(1)%H)) then
           stop "Error in init_structure: invoking block_partition but H not created"
         end if
       end if
       ! supposedly performs an internal block partitioning but it is not reliable.
       call log_allocate(plend_tmp, MAXNUMPLs)
       call block_partition(negf%HS(1)%H, surfend(1), contend, surfend, ncont, npl_tmp, plend_tmp)
       call find_cblocks(negf%HS(1)%H, ncont, npl_tmp, plend_tmp, surfstart, contend, cblk)
       call create_Tstruct(ncont, npl_tmp, plend_tmp, surfstart, surfend, contend, cblk, negf%str)
       call log_deallocate(plend_tmp)
     else
       call create_Tstruct(ncont, npl, plend, surfstart, surfend, contend, cblk, negf%str)
     end if

  end subroutine init_structure



  !> Initialize basis
  subroutine init_basis(negf, coords, nCentral, matrixIndices, latticeVects)
    type(Tnegf) :: negf
    real(dp), intent(in) :: coords(:,:)
    integer, intent(in) :: nCentral
    integer, intent(in) :: matrixIndices(:)
    real(dp), intent(in), optional :: latticeVects(:,:)

    if (present(latticeVects)) then
       call create_TBasis(negf%basis, coords, nCentral, lattVecs=latticeVects, &
             & basisToMatrix=matrixIndices)
    else
       call create_TBasis(negf%basis, coords, nCentral, basisToMatrix=matrixIndices)
    end if

  end subroutine init_basis

  !> Initialize contanct data
  subroutine init_contacts(negf, ncont)
    type(Tnegf) :: negf
    integer, intent(in) :: ncont

    integer :: ii

    ! Make sure that the number of contacts is compatible with naming formatting.
    if (ncont .gt. 99) then
      stop "Too many contacts. Cannot assign default names."
    end if
    ! Deallocate existing contacts if any, then allocate.
    if (allocated(negf%cont)) then
      deallocate(negf%cont)
    end if
    allocate(negf%cont(ncont))

    ! Initialize the structure members to sensible defaults.
    do ii = 1, ncont
      ! Whether the contacts are ficticious and DOS to be used if the contact is
      ! ficticious.
      negf%cont(ii)%FictCont = .false.
      negf%cont(ii)%contact_DOS = 0.0_dp
      ! Electrochemical potentials.
      negf%cont(ii)%mu = 0.0_dp
      negf%cont(ii)%mu_n = 0.0_dp
      negf%cont(ii)%mu_p = 0.0_dp
      ! Electronic temperature for the density matrix calculation.
      negf%cont(ii)%kbT_dm = 0.0_dp
       ! Electronic temperature for the transmission calculation.
      negf%cont(ii)%kbT_t = 0.0_dp
      ! Initialize the names to a default ContactXX, where XX is an index.
      write (negf%cont(ii)%name , "(A7, I2.2)") "Contact", ii
    end do

  end subroutine init_contacts

  !> subroutine used to setup kpoints
  !  k-point sampling must be expressed in reduced coordinates, i.e.
  !  either [-0.5..+0.5]x[-0.5..+0.5] (Gamma-centered) or [0..1]x[0..1] (I quadrant)
  !  kpoints(:)  kweights(:)  are global
  !  local_kindex(:) is a local array storing the local indices
  !  kSamplingType: 0 = Gamma-centered, no inversion
  !                 1 = Shifted in the I quadrant (0..1)x(0..1), no inversion
  subroutine set_kpoints(negf, kpoints, kweights, local_kindex, kSamplingType)
    type(Tnegf) :: negf
    real(dp), intent(in) :: kpoints(:,:)
    real(dp), intent(in) :: kweights(:)
    integer, intent(in) :: local_kindex(:)
    integer, intent(in) :: kSamplingType

    integer :: ii
    real(dp) :: shift(3)

    if (size(kpoints,2) /= size(kweights)) then
       STOP 'Error: size of kpoints do not match'
    end if
    if (allocated(negf%kpoints)) then
       call log_deallocate(negf%kpoints)
    end if
    call log_allocate(negf%kpoints,3,size(kweights))
    if (kSamplingType == 0) then
      negf%kpoints = kpoints
    else if (kSamplingType == 1) then
      ! Try to guess if the system is 2D or 3D
      ! If all k-components are 0 along x or y
      shift = [-0.5_dp, -0.5_dp, 0.0_dp]
      if (all(kpoints(1,:)==0.0_dp)) then
         shift(1) = 0.0_dp
      end if
      if (all(kpoints(2,:)==0.0_dp)) then
         shift(2) = 0.0_dp
      end if
      do ii = 1, size(kweights)
        negf%kpoints(:,ii) = kpoints(:,ii) + shift
      end do
    else
      stop "kSamplingType must be either 0 or 1"
    end if
    if (allocated(negf%kweights)) then
       call log_deallocate(negf%kweights)
    end if
    call log_allocate(negf%kweights,size(kweights))
    negf%kweights = kweights
    if (allocated(negf%local_k_index)) then
       call log_deallocate(negf%local_k_index)
    end if
    call log_allocate(negf%local_k_index,size(local_kindex))
    negf%local_k_index = local_kindex

    if (id0) then
      write(*,*) 'k-points used in NEGF:'
      do ii = 1, size(kweights)
        write(*,*) negf%kpoints(:,ii), negf%kweights(ii)
      end do
    end if
  end subroutine set_kpoints


  !!-------------------------------------------------------------------
  !! Get/Set parameters container
  !!-------------------------------------------------------------------

  !> Set paramters from libnegf. Useful to get defaults
  !!  or to only set some values
  subroutine get_params(negf, params)
    type(Tnegf) :: negf
    type(lnParams), intent(out) :: params

    integer :: nn

    params%verbose = negf%verbose
    params%readOldDM_SGFs = negf%readOldDM_SGFs
    params%readOldT_SGFs = negf%readOldT_SGFs
    params%g_spin = negf%g_spin
    params%delta = negf%delta
    params%deltaModel = negf%deltaModel
    params%wmax = negf%wmax
    params%dos_delta = negf%dos_delta
    nn = size(negf%cont)
    params%mu_n(1:nn) = negf%cont(1:nn)%mu_n
    params%mu_p(1:nn) = negf%cont(1:nn)%mu_p
    params%mu(1:nn) = negf%cont(1:nn)%mu
    params%contact_dos(1:nn) = negf%cont(1:nn)%contact_dos
    params%FictCont(1:nn) = negf%cont(1:nn)%FictCont
    params%kbT_dm(1:nn) = negf%cont(1:nn)%kbT_dm
    params%kbT_t(1:nn) = negf%cont(1:nn)%kbT_t
    params%mu_n(nn+1:MAXNCONT) = 0.0_dp
    params%mu_p(nn+1:MAXNCONT) = 0.0_dp
    params%mu(nn+1:MAXNCONT) = 0.0_dp
    params%contact_dos(nn+1:MAXNCONT) = 0.0_dp
    params%FictCont(nn+1:MAXNCONT) = .false.
    params%kbT_dm(nn+1:MAXNCONT) = 0.0_dp
    params%kbT_t(nn+1:MAXNCONT) = 0.0_dp
    params%scba_inelastic_tol = negf%scba_inelastic_tol
    params%scba_elastic_tol = negf%scba_elastic_tol
    if (nn == 0) then
      params%mu_n(1) = negf%mu_n
      params%mu_p(1) = negf%mu_p
      params%mu(1) = negf%mu
      params%kbT_dm(1) = negf%kbT
    end if
    params%Np_n = negf%Np_n
    params%Np_real = negf%Np_real
    params%n_kt = negf%n_kt
    params%n_poles = negf%n_poles
    params%Ec = negf%Ec
    params%Ev = negf%Ev
    params%DeltaEc = negf%DeltaEc
    params%DeltaEv = negf%DeltaEv
    params%Emin = negf%Emin
    params%Emax = negf%Emax
    params%Estep = negf%Estep
    nn = size(negf%ni)
    params%ni=0; params%ni(1:nn) = negf%ni(1:nn)
    params%nf=0; params%nf(1:nn) = negf%nf(1:nn)
    params%eneconv = negf%eneconv
    params%spin = negf%spin
    params%kwght = negf%kwght
    params%ikpoint = negf%ikpoint
    params%DorE = negf%DorE
    params%min_or_max = negf%min_or_max
    params%SGFcache = get_surface_green_cache_type(negf)

  end subroutine get_params

  function get_surface_green_cache_type(negf) result(idx)
    type(Tnegf) :: negf
    integer :: idx

    select type (sgf => negf%surface_green_cache)
    type is (TMatrixCacheDisk)
      idx = 0
    type is (TMatrixCacheMem)
      idx = 1
    type is (TMatrixCacheDummy)
      idx = 2
    class default
      idx = 2
    end select
  end function get_surface_green_cache_type

  !> Assign parameters to libnegf
  subroutine set_params(negf, params)
    type(Tnegf) :: negf
    type(lnParams), intent(in) :: params

    integer :: nn, tmp

    negf%verbose = params%verbose
    negf%readOldDM_SGFs = params%readOldDM_SGFs
    negf%readOldT_SGFs = params%readOldT_SGFs
    negf%g_spin = params%g_spin
    negf%delta = params%delta
    negf%deltaModel = params%deltaModel
    negf%wmax = params%wmax
    negf%dos_delta = params%dos_delta
    nn = size(negf%cont)
    negf%cont(1:nn)%mu_n        = params%mu_n(1:nn)
    negf%cont(1:nn)%mu_p        = params%mu_p(1:nn)
    negf%cont(1:nn)%mu          = params%mu(1:nn)
    negf%cont(1:nn)%contact_dos = params%contact_dos(1:nn)
    negf%cont(1:nn)%FictCont    = params%FictCont(1:nn)
    negf%cont(1:nn)%kbT_dm      = params%kbT_dm(1:nn)
    negf%cont(1:nn)%kbT_t       = params%kbT_t(1:nn)
    negf%scba_inelastic_tol = params%scba_inelastic_tol
    negf%scba_elastic_tol   = params%scba_elastic_tol
    if (nn == 0) then
      negf%mu   = params%mu(1)
      negf%mu_n = params%mu_n(1)
      negf%mu_p = params%mu_p(1)
      negf%kbT = params%kbT_dm(1)
    end if
    negf%Np_n = params%Np_n
    negf%Np_p = params%Np_p
    negf%Np_real = params%Np_real
    negf%n_kt = params%n_kt
    negf%n_poles = params%n_poles
    negf%Ec = params%Ec
    negf%Ev = params%Ev
    negf%DeltaEc = params%DeltaEc
    negf%DeltaEv = params%DeltaEv
    negf%Emin = params%Emin
    negf%Emax = params%Emax
    negf%Estep = params%Estep
    negf%Estep_coarse = params%Estep_coarse
    if (allocated(negf%ni)) deallocate(negf%ni)
    nn = count(params%ni .ne. 0)
    allocate(negf%ni(nn))
    negf%ni(1:nn) = params%ni(1:nn)
    if (allocated(negf%nf)) deallocate(negf%nf)
    nn = count(params%nf .ne. 0)
    allocate(negf%nf(nn))
    negf%nf(1:nn) = params%nf(1:nn)
    negf%eneconv = params%eneconv
    negf%spin = params%spin
    negf%kwght = params%kwght
    negf%ikpoint = params%ikpoint
    negf%DorE = params%DorE
    negf%min_or_max = params%min_or_max

    !! Some internal variables in libnegf are set internally
    !! after parameters are available
    call set_ref_cont(negf)

    ! The surface green cache is created if there's none. Otherwise
    ! it is reset only if the type is changed. If the energy integral
    ! is changed and the cache type is not, it must be forcibly destroyed
    ! using destroy_surface_green_cache
    tmp = get_surface_green_cache_type(negf)
    select case(params%SGFcache)
    case(0)
      if (tmp .ne. 0) then
        if (allocated(negf%surface_green_cache)) then
           call negf%surface_green_cache%destroy()
           deallocate(negf%surface_green_cache)
        end if
        negf%surface_green_cache = TMatrixCacheDisk(scratch_path=negf%scratch_path)
      end if
    case(1)
      if (tmp .ne. 1) then
        if (allocated(negf%surface_green_cache)) then
           call negf%surface_green_cache%destroy()
           deallocate(negf%surface_green_cache)
        end if
        negf%surface_green_cache = TMatrixCacheMem(tagname='SurfaceGF')
      end if
    case(2)
      if (tmp .ne. 2) then
        if (allocated(negf%surface_green_cache)) then
           call negf%surface_green_cache%destroy()
           deallocate(negf%surface_green_cache)
        end if
        negf%surface_green_cache = TMatrixCacheDummy()
      end if
    end select

  end subroutine set_params

  !----------------------------------------------------------
  subroutine set_scba_tolerances(negf, elastic_tol, inelastic_tol)
    type(Tnegf) :: negf
    real(dp), intent(in) :: elastic_tol
    real(dp), intent(in) :: inelastic_tol
    negf%scba_elastic_tol = elastic_tol
    negf%scba_inelastic_tol = inelastic_tol
  end subroutine set_scba_tolerances
  !----------------------------------------------------------
  subroutine set_scratch(negf, scratchpath)
    type(Tnegf) :: negf
    character(*), intent(in) :: scratchpath

    if (len(scratchpath)>LST) then
       print*, "ERROR: scratch string longer than",LST
       stop
    end if
    !negf%scratch_path = trim(scratchpath)//'/GS/'
    negf%scratch_path = trim(scratchpath)//'/'

    ! Update the cache object if needed.
    select type (sgf => negf%surface_green_cache)
      type is (TMatrixCacheDisk)
        sgf%scratch_path = negf%scratch_path
    end select

  end subroutine set_scratch
  !----------------------------------------------------------

  !----------------------------------------------------------
  subroutine set_outpath(negf, outpath)
    type(Tnegf) :: negf
    character(LST), intent(in) :: outpath

    negf%out_path = trim(outpath)//'/'

  end subroutine set_outpath
  !----------------------------------------------------------

  subroutine create_scratch(negf)
    type(Tnegf) :: negf

    call create_directory(trim(negf%scratch_path))

  end subroutine create_scratch


  !--------------------------------------------------------------------
  ! dos_proj methods: you can set N index intervals OR N separate index
  ! arrays. You have to initialize the data by indicating the number of
  ! ldos interval (nldos) and then you can either set the start/end
  ! indexes for intervals OR append one by one explicit arrays
  !
  !--------------------------------------------------------------------

  !> Initialize the ldos info
  !! @param [in] negf: libnegf container instance
  !! @param [in] nldos: number of intervals
  subroutine init_ldos(negf,nldos)
    type(Tnegf) :: negf
    integer, intent(in) :: nldos

    integer :: i

    if (allocated(negf%dos_proj)) then
      do i=1, size(negf%dos_proj)
        if (allocated(negf%dos_proj(i)%indexes)) then
          call log_deallocate(negf%dos_proj(i)%indexes)
        end if
      end do
      deallocate(negf%dos_proj)
    end if
    allocate(negf%dos_proj(nldos))
    negf%ndos_proj = nldos

  end subroutine init_ldos

  !> Destroy the dos_proj container
  subroutine destroy_ldos(ldos)
    type(intarray), dimension(:), allocatable :: ldos

    integer :: err, i
    do i=1, size(ldos)
      if (allocated(ldos(i)%indexes)) then
        call log_deallocate(ldos(i)%indexes)
      end if
    end do

    deallocate(ldos)

  end subroutine destroy_ldos


  !> Set ldos intervals
  !! @param [in] negf: libnegf container instance
  !! @param [in] istart(nldos) array with first interval index
  !! @param [in] iend(nldos) array with first interval index
  subroutine set_ldos_intervals(negf, nldos, istart, iend)
    type(Tnegf) :: negf
    integer, intent(in) :: nldos
    integer, intent(in) :: istart(*), iend(*)

    integer :: ii, jj

    do ii = 1, negf%ndos_proj
      call log_allocate(negf%dos_proj(ii)%indexes,iend(ii)-istart(ii)+1)
      do jj = 1, iend(ii) - istart(ii) + 1
        negf%dos_proj(ii)%indexes(jj) = istart(ii) + jj - 1
      end do
    end do

  end subroutine set_ldos_intervals

  !> Set ldos indexes arrays for a given ldos
  !!
  !! @param [in] negf: libnegf container instance
  !! @param [in] ildos: index of ldos
  !! @param [in] idx: 1D array with indexes
  subroutine set_ldos_indexes(negf, ildos, idx)
    type(Tnegf) :: negf
    integer, intent(in) :: ildos
    integer, intent(in) :: idx(:)

    integer :: ii, jj

    if (.not.allocated(negf%dos_proj(ildos)%indexes)) then
       call log_allocate(negf%dos_proj(ildos)%indexes, size(idx))
    end if
    negf%dos_proj(ildos)%indexes = idx

  end subroutine set_ldos_indexes
  ! -------------------------------------------------------------------
  !> Set tunneling projection indexes array
  subroutine set_tun_indexes(negf, idx)
    type(Tnegf) :: negf
    integer, intent(in) :: idx(:)

    integer :: ii, jj

    if (.not.allocated(negf%tun_proj%indexes)) then
      call log_allocate(negf%tun_proj%indexes, size(idx))
    else
      if (size(negf%tun_proj%indexes) /= size(idx)) then
         write(*,*) 'ERROR in set_tun_indexes size mismatch'
         return
      end if
    end if
    negf%tun_proj%indexes = idx

  end subroutine set_tun_indexes
  ! -------------------------------------------------------------------

#:if defined("MPI")

  ! Purpose: setup negf communicators
  ! globalComm
  ! cartComm
  ! energyComm
  ! kComm
  !
  subroutine negf_mpi_init(negf, cartComm, energyComm, kComm)
    type(Tnegf) :: negf
    type(mpifx_comm), intent(in), optional :: cartComm
    type(mpifx_comm), intent(in), optional :: energyComm
    type(mpifx_comm), intent(in), optional :: kComm

    integer :: err

    if (present(cartComm)) then
       if (.not.present(energyComm) .and. .not.present(kComm)) then
          stop "ERROR: cartesian communicator also requires energy and k comm"
       end if
       call check_cart_comm(cartComm, err)
       if (err /= 0) then
         stop "ERROR: pass non cartesian communicator to negf_mpi_init"
       end if
       negf%globalComm = cartComm
       negf%cartComm = cartComm
       negf%energyComm = energyComm
       negf%kComm = kComm
       id0 = (negf%cartComm%rank == 0)
    else
       if (.not.present(energyComm)) then
          stop "ERROR: negf_mpi_init needs at lest the energy communicator"
       end if
       negf%globalComm = energyComm
       negf%energyComm = energyComm
       if (present(kComm)) then
          negf%kComm = kComm
       else
          negf%kComm%id = 0
       end if
       negf%cartComm%id = 0
       id0 = (negf%energyComm%rank == 0)
    end if
    numprocs = negf%energyComm%size
    id = negf%energyComm%rank

  end subroutine negf_mpi_init

  !> Set a global mpifx communicator from a bare communicator.
  !!
  !! @param [in] negf: libnegf container instance
  !! @param [in] mpicomm: an mpi communicator
  subroutine set_mpi_bare_comm(negf, mpicomm)
    type(Tnegf), intent(inout) :: negf
    integer, intent(in) :: mpicomm

    call negf%globalComm%init(mpicomm)
    call negf%energyComm%init(mpicomm)
    call globals_mpi_init(negf%energyComm)

  end subroutine

  subroutine set_cartesian_bare_comms(negf, mpicomm, nk, cartComm, kComm)
    type(Tnegf), intent(inout) :: negf
    integer, intent(in) :: mpicomm
    integer, intent(in) :: nk
    integer, intent(out) :: cartComm
    integer, intent(out) :: kComm

    call negf%globalComm%init(mpicomm)

    call negf_cart_init(negf%globalComm, nk, negf%cartComm, negf%energyComm, negf%kComm, cartComm, kComm)
    call negf_mpi_init(negf, negf%cartComm, negf%energyComm, negf%kComm)

  end subroutine set_cartesian_bare_comms
#:else
  !> Dummy method for the C-interface, when mpi implementation is missing.
  !!
  !! @param [in] negf: libnegf container instance
  !! @param [in] mpicomm: unused.
  subroutine set_mpi_bare_comm(negf, mpicomm)
    type(Tnegf), intent(inout) :: negf
    integer, intent(in) :: mpicomm

  end subroutine
#:endif
  ! -------------------------------------------------------------------
  subroutine set_convfactor(negf, eneconv)
    type(Tnegf) :: negf
    real(dp) :: eneconv

    negf%eneconv=eneconv

  end subroutine set_convfactor

  ! -------------------------------------------------------------------
  subroutine set_fictcont(negf,cont,dos)
    type(Tnegf) :: negf
    integer :: cont
    real(dp) :: DOS

    negf%cont(cont)%FictCont = .true.
    negf%cont(cont)%contact_DOS = DOS

  end subroutine set_fictcont

  ! -------------------------------------------------------------------

  subroutine set_computation(negf,DorE)
    type(Tnegf) :: negf
    character(1) :: DorE           !Density or En.Density

    negf%DorE=DorE
  end subroutine set_computation

  ! -------------------------------------------------------------------
  subroutine set_readOldDMsgf(negf,flag)
    type(Tnegf) :: negf
    integer :: flag !between 0:2

    negf%ReadOldDM_SGFs = flag
  end subroutine set_readOldDMsgf

  ! -------------------------------------------------------------------
  subroutine set_readOldTsgf(negf,flag)
    type(Tnegf) :: negf
    integer :: flag ! between 0:2

    negf%ReadOldT_SGFs = flag
  end subroutine set_readOldTsgf

  !--------------------------------------------------------------------
  !> Initialize and set parameters from input file negf.in
  subroutine read_negf_in(negf)
    type(Tnegf) :: negf

    integer :: ncont, nbl, ii, jj, ist, iend
    integer, dimension(:), allocatable :: PL_end, cont_end, surf_end
    integer, dimension(:), allocatable :: surf_start, cblk
    character(32) :: tmp
    character(LST) :: file_re_H, file_im_H, file_re_S, file_im_S

    open(101, file="negf.in", form='formatted')

    read(101,*) tmp, file_re_H
    read(101,*) tmp, file_im_H

    read(101,*) tmp, file_re_S
    read(101,*) tmp, file_im_S

    if (trim(file_re_H).ne.'memory') then
      call read_HS(negf, file_re_H, file_im_H, 0)
    end if
    if (trim(file_re_S).ne.'memory') then
      if (trim(file_re_S).eq.'identity') then
         call set_S_id(negf, negf%H%nrow, 1)
      end if
    else
      call read_HS(negf, file_re_S, file_im_S, 1)
    end if

    !! A dummy descriptor must be present at the beginning of each line
    !! Used for debug
    read(101,*) tmp, ncont
    call log_allocate(cblk, ncont)
    call log_allocate(surf_start,ncont)
    call log_allocate(surf_end,ncont)
    call log_allocate(cont_end,ncont)

    read(101,*) tmp, nbl
    if (nbl .gt. 0) then
       call log_allocate(PL_end,nbl)
       read(101,*) tmp, PL_end(1:nbl)
    end if
    read(101,*) tmp, surf_start(1:ncont)
    read(101,*) tmp, surf_end(1:ncont)
    read(101,*) tmp, cont_end(1:ncont)

    call find_cblocks(negf%HS(1)%H, ncont, nbl, PL_end, surf_start, cont_end, cblk)
    call init_structure(negf, ncont, surf_start, surf_end, cont_end, nbl, PL_end, cblk)

    call log_deallocate(PL_end)
    call log_deallocate(cont_end)
    call log_deallocate(surf_end)
    call log_deallocate(surf_start)
    call log_deallocate(cblk)

    read(101,*)  tmp, negf%Ec, negf%Ev
    read(101,*)  tmp, negf%DeltaEc, negf%DeltaEv
    read(101,*)  tmp, negf%Emin, negf%Emax, negf%Estep
    if (ncont.gt.0) then
      read(101,*) tmp, negf%cont(1:ncont)%kbT_dm
    else
      read(101,*) tmp, negf%kbT
    endif
    read(101,*)  tmp, negf%kwght
    read(101,*)  tmp, negf%Np_n(1:2)
    read(101,*)  tmp, negf%Np_p(1:2)
    read(101,*)  tmp, negf%Np_real
    read(101,*)  tmp, negf%n_kt
    read(101,*)  tmp, negf%n_poles
    read(101,*)  tmp, negf%g_spin
    read(101,*)  tmp, negf%delta
    read(101,*)  tmp, negf%ndos_proj
    if (allocated(negf%dos_proj)) then
      deallocate(negf%dos_proj)   !DAR
    endif
    allocate(negf%dos_proj(negf%ndos_proj))
    do ii = 1, negf%ndos_proj
      read(101,*) tmp, ist, iend
      call log_allocate(negf%dos_proj(ii)%indexes, iend-ist+1)
      do jj = 1, iend-ist+1
        negf%dos_proj(ii)%indexes(jj) = ist + jj - 1
      end do
    end do
    if (ncont.gt.0) then
      read(101,*) tmp, negf%cont(1:ncont)%mu_n  ! Will be the Electrochemical potential
      read(101,*) tmp, negf%cont(1:ncont)%mu_p  ! hole potentials
    else
      read(101,*) tmp, negf%mu_n                ! Will be the Electrochemical potential
      read(101,*) tmp, negf%mu_p                ! hole potentials
    end if
    !! Internally a different mu is used for dft-like integrations
    !! we define it as equal to mu_n in negf.in
    negf%cont(1:ncont)%mu = negf%cont(1:ncont)%mu_n

    close(101)

  end subroutine read_negf_in

  !--------------------------------------------------------------------
  subroutine negf_version(negf)
    type(Tnegf) :: negf
    character(3), parameter :: GITVER= "${_GITREVISION}$"
    character(10),parameter :: DATE= "${_COMPDATE}$"

    write(*,'(a21,a20,2x,a10)') '(libnegf) version: 1.',TRIM(GITVER), &
                                         TRIM(DATE)

  end subroutine negf_version

!--------------------------------------------------------------------
  subroutine negf_partition_info(negf)
      type(Tnegf) :: negf

      integer :: i

      write(*,*) "(LibNEGF) Partitioning:"
      write(*,*) "Number of blocks: ",negf%str%num_Pls
      !write(*,*) negf%str%mat_PL_end(:)
      write(*,*) "Contact interactions:",negf%str%cblk(:)

      open(1001,file='blocks.dat')
        write(1001,*) 1
        do i = 1, negf%str%num_Pls
           write(1001,*)  negf%str%mat_PL_end(i)
        enddo
      close(1001)

  end subroutine negf_partition_info

  !--------------------------------------------------------------------
  !> Destroy all the info defined in initialization.
  !! To run at the very end of libnegf usage
  subroutine destroy_negf(negf)
    type(Tnegf) :: negf

    call destroy_HS(negf)
    call kill_Tstruct(negf%str)
    call destroy_TBasis(negf%basis)
    if (allocated(negf%dos_proj)) then
       call destroy_ldos(negf%dos_proj)
    end if
    if (allocated(negf%tun_proj%indexes)) then
       call log_deallocate(negf%tun_proj%indexes)
    end if
    if (allocated(negf%en_grid)) then
       deallocate(negf%en_grid)
    end if
    if (allocated(negf%tunn_mat)) then
       call log_deallocate(negf%tunn_mat)
    end if
    if (allocated(negf%curr_mat)) then
       call log_deallocate(negf%curr_mat)
    end if
    if (allocated(negf%ldos_mat)) then
         call log_deallocate(negf%ldos_mat)
    end if
    if (allocated(negf%currents)) then
      call log_deallocate(negf%currents)
    end if
    if (allocated(negf%kpoints)) then
      call log_deallocate(negf%kpoints)
    end if
    if (allocated(negf%kweights)) then
      call log_deallocate(negf%kweights)
    end if
    if (allocated(negf%local_k_index)) then
      call log_deallocate(negf%local_k_index)
    end if

    call destroy_interactions(negf)

    call destroy_DM(negf)
    call destroy_contact_matrices(negf)
    call destroy_surface_green_cache(negf)
    call destroy_cache_space(negf)

  end subroutine destroy_negf

  !> Destroy the surface green cache.
  subroutine destroy_surface_green_cache(negf)
    type(Tnegf) :: negf

    if (allocated(negf%surface_green_cache)) then
          call negf%surface_green_cache%destroy()
    end if

  end subroutine destroy_surface_green_cache

  !--------------------------------------------------------------------
  !> Copy the energy axis on all processors (for output, plot, debug)
  !! @param [in] negf: negf container
  !! @param [out] energies: energy values, it can eb allocated internally
  subroutine get_energies(negf, energies)
    type(Tnegf), intent(in) :: negf
    complex(dp), allocatable :: energies(:)

    if (.not.allocated(energies)) then
      allocate(energies(size(negf%en_grid)))
    end if
    energies = negf%en_grid(:)%Ec

  end subroutine get_energies

  !--------------------------------------------------------------------
  !>
  !! Associate an input pointer with the internal pointer of
  !! transmissions. Return NULL if internal pointer is not
  !! associated
  subroutine associate_transmission(negf, tr_pointer)
    type(TNegf), pointer, intent(in)  :: negf
    real(dp), dimension(:,:), pointer, intent(inout) :: tr_pointer

    if (allocated(negf%tunn_mat)) then
      tr_pointer => negf%tunn_mat
    else
      tr_pointer => NULL()
    end if

  end subroutine associate_transmission

  !--------------------------------------------------------------------
  !>
  !!  Associate an input pointer with the internal pointer of
  !! dos_proj
  subroutine associate_ldos(negf, ldos_pointer)
    type(TNegf), pointer, intent(in)  :: negf
    real(dp), dimension(:,:), pointer, intent(inout) :: ldos_pointer

    if (allocated(negf%ldos_mat)) then
      ldos_pointer => negf%ldos_mat
    else
      ldos_pointer => NULL()
    end if

  end subroutine associate_ldos

  !--------------------------------------------------------------------
  !>
  !!  Associate an input pointer with the internal pointer of
  !! currents
  subroutine associate_current(negf, curr_pointer)
    type(TNegf), pointer, intent(in)  :: negf
    real(dp), dimension(:,:), pointer, intent(inout) :: curr_pointer

    if (allocated(negf%curr_mat)) then
      curr_pointer => negf%curr_mat
    else
      curr_pointer => NULL()
    end if

  end subroutine associate_current


  !--------------------------------------------------------------------
  subroutine associate_lead_currents(negf, curr)
    type(TNegf), pointer, intent(in)  :: negf
    real(dp), dimension(:), pointer, intent(inout) :: curr

    if (allocated(negf%currents)) then
      curr => negf%currents
    else
      curr => null()
    end if

  end subroutine associate_lead_currents


  !--------------------------------------------------------------------
  !>
  !! Get currents by copy.
  !! @param [in] negf: negf container
  !! @param [out] currents: current values, it can eb allocated internally
  subroutine get_currents(negf, currents)
    type(TNegf), intent(in)  :: negf
    real(dp), intent(out) :: currents(:)

    currents = negf%currents(:)
  end subroutine get_currents

  !> Get density matrix CSR sparse arrays by copy
  !! @param [in] negf: negf container
  !! @param [out] nzval: number of non zero values
  !! @param [out] nrow: number of rows
  !! @param [out] rowpnt (int array): row pointer indexes
  !! @param [out] colind (int array): column indexes array
  !! @param [out] nzval (complex array): non zero values
  subroutine get_DM(negf, nnz, nrow, rowpnt, colind, nzval)
    type(TNegf), intent(in)  :: negf
    integer, intent(out) :: nnz, nrow
    integer, intent(out) :: rowpnt(:), colind(:)
    real(dp), intent(out) :: nzval(:)

    nnz = negf%rho%nnz
    nrow = negf%rho%nrow
    rowpnt = negf%rho%rowpnt
    colind = negf%rho%colind
    nzval = real(negf%rho%nzval)
  end subroutine get_DM

  !--------------------------------------------------------------------
  subroutine create_DM(negf)
    type(Tnegf) :: negf

    if (negf%internalDM) then
      if (.not.associated(negf%rho)) allocate(negf%rho)
      if (.not.associated(negf%rho_eps)) allocate(negf%rho_eps)
    end if

  end subroutine create_DM

  ! -----------------------------------------------------
  !  Pass an externally allocated density matrix
  ! -----------------------------------------------------
  subroutine pass_DM(negf,rho, rhoE)
    type(Tnegf) :: negf
    type(z_CSR), pointer, intent(in), optional :: rho
    type(z_CSR), pointer, intent(in), optional :: rhoE

    if (present(rho)) then
       negf%rho => rho
       if(allocated(negf%rho%nzval)) then
          call destroy(negf%rho)
       endif
    endif

    if (present(rhoE)) then
       negf%rho_eps => rhoE
       if (allocated(negf%rho_eps%nzval)) then
          call destroy(negf%rho_eps)
       endif
    end if

    negf%internalDM = .false.

  end subroutine pass_DM

  !--------------------------------------------------------------------
  subroutine destroy_HS(negf)
    type(Tnegf) :: negf

    integer :: ii

    if (.not.allocated(negf%HS)) return

    do ii = 1, size(negf%HS)
      if (negf%HS(ii)%internalHS) then
        if (associated(negf%HS(ii)%H)) then
          if (allocated(negf%HS(ii)%H%nzval)) then
             !print*,'(destroy) deallocate negf%H',%LOC(negf%H%nzval)
             call destroy(negf%HS(ii)%H)
          end if
          deallocate(negf%HS(ii)%H)
          nullify(negf%HS(ii)%H)
        endif

        if (associated(negf%HS(ii)%S)) then
          if (allocated(negf%HS(ii)%S%nzval)) then
             !print*,'(destroy) deallocate negf%S',%LOC(negf%S%nzval)
             call destroy(negf%HS(ii)%S)
          end if
          deallocate(negf%HS(ii)%S)
          nullify(negf%HS(ii)%S)
        endif
      endif
    end do
    deallocate(negf%HS)

  end subroutine destroy_HS

  !--------------------------------------------------------------------
  subroutine destroy_DM(negf)
    type(Tnegf) :: negf

    integer :: ii

    if (.not.negf%internalDM) return

    if (associated(negf%rho)) then
      if (allocated(negf%rho%nzval)) then
         call destroy(negf%rho)
      end if
      deallocate(negf%rho)
      nullify(negf%rho)
    end if

    if (associated(negf%rho_eps)) then
      if (allocated(negf%rho_eps%nzval)) then
         call destroy(negf%rho_eps)
      end if
      deallocate(negf%rho_eps)
      nullify(negf%rho_eps)
    endif

   ! if (.not.allocated(negf%DM)) return

   ! do ii = 1, size(negf%DM)
   !   if (negf%DM(ii)%internalDM) then
   !     if (associated(negf%DM(ii)%rho)) then
   !       if (allocated(negf%DM(ii)%rho%nzval)) then
   !          !print*,'(destroy) deallocate negf%rho',%LOC(negf%rho%nzval)
   !          call destroy(negf%DM(ii)%rho)
   !       end if
   !       deallocate(negf%DM(ii)%rho)
   !       nullify(negf%DM(ii)%rho)
   !     endif

   !     if (associated(negf%DM(ii)%rho_eps)) then
   !       if (allocated(negf%DM(ii)%rho_eps%nzval)) then
   !          !print*,'(destroy) deallocate negf%rho_eps',%LOC(negf%rho_eps%nzval)
   !          call destroy(negf%DM(ii)%rho_eps)
   !       end if
   !       deallocate(negf%DM(ii)%rho_eps)
   !       nullify(negf%DM(ii)%rho_eps)
   !     endif
   !   endif
   ! end do
   ! deallocate(negf%DM)

  end subroutine destroy_DM

  !-------------------------------------------------------------------------------
  ! Compact collection of calls to extract device/contact H and S
  ! and compute density matrix using contour + real axis integration
  ! Should be used for dftt calculations
  !
  ! NOTE: the returned DensMat and EnMat are masked with S
  ! Matrix Structure:  CSR
  !                    %nrow=%ncol=(Full squared Hamiltonian size)
  !                    %nnz = Only non-zero elements of the blocks
  !
  !                    +-----+--+--+--+
  !                    !  D  !C1!C2!C3!  masked with the S matrix
  !                    !     !  !  !  !
  !                    +-----+--+--+--+
  !                    ! C1  !0 !0 !0 !  The lower part of DensMat
  !                    +-----+--+--+--+  is filled with 0.
  !                    ! C2  !0 !0 !0 !
  !                    +-----+--+--+--+  negf%outer=0,1,2 is used
  !                    ! C3  !0 !0 !0 !  in order to compute Ci
  !                    +-----+--+--+--+
  !-------------------------------------------------------------------------------
  subroutine compute_density_dft(negf)
    type(Tnegf) :: negf

    !Temporary hack
    negf%H => negf%HS(1)%H
    negf%S => negf%HS(1)%S

    call extract_cont(negf)

    !! Did anyone passed externally allocated DM? If not, create it
    call create_DM(negf)

    ! Reference contact for contour/real axis separation
    call set_ref_cont(negf)

    if (negf%Np_n(1)+negf%Np_n(2)+negf%n_poles.gt.0) then
      call contour_int_def(negf)
      call contour_int(negf)
    endif

    if (negf%Np_real.gt.0) then
      call real_axis_int_def(negf)
      call real_axis_int(negf)
    endif

    call destroy_contact_matrices(negf)

  end subroutine compute_density_dft


  !-------------------------------------------------------------------------------
  ! Compact collection of calls to extract device/contact H and S
  ! and compute density matrix
  !
  ! It has been used to interface libnegf to TiberCAD
  ! Computes density for CB or VB semiconductor
  !-------------------------------------------------------------------------------
  subroutine compute_density_efa(negf, q, particle)

    type(Tnegf) :: negf
    real(dp), dimension(:) :: q
    integer :: particle  ! +1 for electrons, -1 for holes
    complex(dp), dimension(:), allocatable :: q_tmp

    integer :: k

    if (particle /= +1 .and. particle /= -1) then
       write(*,*) "libNEGF error. In compute_density_efa, unknown particle"
       stop
    endif

    call extract_cont(negf)

    call create_DM(negf)
    ! it is not zeroed out at end of integration, so we have to do it here
    ! (otherwise k-integration from tibercad will not work)
    if (allocated(negf%rho%nzval)) then
      call destroy(negf%rho)
    end if

    ! Reference contact for contour/real axis separation
    call set_ref_cont(negf)

    ! Contour integral for equilibrium calculations
    if (particle == 1) then
      if (negf%str%num_conts > 0) then
        negf%muref = negf%cont(1)%mu_n
      else
        negf%muref = negf%mu_n
      endif

      if (negf%Np_n(1)+negf%Np_n(2)+negf%n_poles.gt.0) then
         call contour_int_n_def(negf)
         call contour_int(negf)
      else
         ! HACKING: THIS WAY COMPUTES DM FOR ALL CONTACTS
         negf%refcont = negf%str%num_conts+1
      endif
    else ! particle == -1
      if (negf%str%num_conts > 0) then
        negf%muref = negf%cont(1)%mu_p
      else
        negf%muref = negf%mu_p
      endif

      if (negf%Np_p(1)+negf%Np_p(2)+negf%n_poles.gt.0) then
         call contour_int_p_def(negf)
         call contour_int(negf)
      else
         ! HACKING: THIS WAY COMPUTES DM FOR ALL CONTACTS
         negf%refcont = negf%str%num_conts+1
      endif
    endif

    ! Real axis integral for non-equilibrium calculations
    if (negf%Np_real.gt.0) then
       if (particle == 1) then
          negf%particle = 1
          call real_axis_int_n_def(negf)
       else
          negf%particle = -1
          call real_axis_int_p_def(negf)
       endif
       ! we use contour_int here because it integrates Gr, while
       ! real_axis_int integrates Gn
       !call contour_int(negf)
       call real_axis_int(negf)
    endif

    ! We need not to include S:
    ! rho(r) = sum_ij ui(r) Pij uj(r)
    ! On the mesh nodes:
    ! rho(rk) = sum_ii Pii ui(rk)^2
    if (negf%rho%nrow.gt.0) then
       call log_allocate(q_tmp, negf%rho%nrow)

       call getdiag(negf%rho, q_tmp)

       do k = 1, size(q)
          q(k) = real(q_tmp(k))
       enddo

       call log_deallocate(q_tmp)
    else
       q = 0.0_dp
    endif

    call destroy_contact_matrices(negf)

  end subroutine compute_density_efa

  subroutine compute_density_quasiEq(negf, q, particle, &
                                     Ec, Ev, mu_n, mu_p)
    !In/Out
    type(Tnegf) :: negf
    real(dp), dimension(:) :: q, Ec, Ev, mu_n, mu_p
    integer :: particle  ! +1 for electrons, -1 for holes

    if (particle /= +1 .and. particle /= -1) then
       write(*,*) "libNEGF error. In compute_density_quasiEq, unknown particle"
       stop
    endif


    call extract_cont(negf)
    q = 0.0_dp

    if (particle == 1) then
      call quasiEq_int_n(negf, mu_n, Ec, q)
    else ! particle == -1
      call quasiEq_int_p(negf, mu_p, Ev, q)
    endif
   
    call destroy_contact_matrices(negf)

  end subroutine compute_density_quasiEq
  !-------------------------------------------------------------------------------
  subroutine compute_ldos(negf)
    type(Tnegf) :: negf

    call extract_cont(negf)
    call tunneling_int_def(negf)
    call ldos_int(negf)
    call destroy_contact_matrices(negf)

  end subroutine compute_ldos

  !-------------------------------------------------------------------------------
  subroutine compute_current(negf)
    type(Tnegf) :: negf

    integer :: fu
    open(newunit=fu, file='H.dat')
    call zprint_csrcoo(fu,negf%HS(1)%H,'r')
    close(fu)

    if ( negf%interactList%counter > 0 ) then
       if (get_max_wq(negf%interactList) == 0.0_dp) then
          call compute_dephasing_transmission(negf)
       end if
    else
       call compute_landauer(negf);
    endif

  end subroutine compute_current

  !-------------------------------------------------------------------------------
  !> Calculate current, tunneling and, if specified, density of states using
  !! Landauer formula. DOS is calculated during the T(E) loop
  !! @param negf input/output container
  subroutine compute_landauer(negf)

    type(Tnegf) :: negf

    call extract_cont(negf)
    call tunneling_int_def(negf)
    
    call tunneling_and_dos(negf)

    call electron_current(negf)

    call destroy_contact_matrices(negf)

  end subroutine compute_landauer

  !-------------------------------------------------------------------------------
  !> Calculate current and tunnelling for elastic el-ph dephasing models.
  !> Since the "real" landauer-like formula is not implemented yet, we use
  !> a dirty trick only valid for 2 contacts.
  !! @param negf input/output container
  subroutine compute_dephasing_transmission(negf)

    type(Tnegf) :: negf
    real(dp), allocatable, dimension(:) :: occupations

    if (negf%str%num_conts .ne. 2) then
      error stop "Effective transmission is only supported for 2 electrodes"
    end if

    call tunneling_int_def(negf)
    ! Dirty trick. Set the contact population to 1 on the final contact and
    ! 1 on the initial one.
    allocate(occupations(2))
    occupations(negf%ni(1)) = 1.0_dp
    occupations(negf%nf(1)) = 0.0_dp

    call meir_wingreen(negf, fixed_occupations=occupations)
    ! Assign the current matrix values to the transmission.
    negf%tunn_mat = negf%curr_mat

    call electron_current(negf)

  end subroutine compute_dephasing_transmission

  !-------------------------------------------------------------------------------
  !> Calculate current, tunneling and, if specified, density of states using
  !! Landauer formula. DOS is calculated during the T(E) loop
  !! @param negf input/output container
  subroutine compute_meir_wingreen(negf)

    type(Tnegf) :: negf

    call tunneling_int_def(negf)
    call meir_wingreen(negf)
    call destroy_contact_matrices(negf)

  end subroutine compute_meir_wingreen

  !-------------------------------------------------------------------------------
  !> Calculate current, tunneling and, if specified, density of states using
  subroutine compute_layer_current(negf)

    type(Tnegf) :: negf

    call tunneling_int_def(negf)
    call layer_current(negf)
    call destroy_contact_matrices(negf)

  end subroutine compute_layer_current

  ! --------------------------------------------------------------------------------
  ! GP Left in MPI version for debug purpose only. This will write a separate
  ! file for every ID, which is not possible on all architectures
  subroutine write_current(negf)
    type(Tnegf) :: negf

    integer :: i1
    logical :: lex
    character(6) :: idstr

    write(idstr,'(i6.6)') id
    inquire(file=trim(negf%out_path)//'current_'//idstr//'.dat',EXIST=lex)

    if (lex) then
       open(101,file=trim(negf%out_path)//'current_'//idstr//'.dat',position='APPEND')
    else
       open(101,file=trim(negf%out_path)//'current_'//idstr//'.dat')
    endif

    do i1=1,size(negf%currents)

       write(101,'(1x,a,i3,i3,a,i3,a,ES14.5,a,ES14.5,a)') 'contacts:',negf%ni(i1),negf%nf(i1), &
            '; k-point:',negf%ikpoint,'; current:', negf%currents(i1),' A'

    end do

    close(101)

  end subroutine write_current
  !-------------------------------------------------------------------------------

  !---- RETURN THE DOS MATRIX ---------------------------------------------------------
  subroutine return_dos_mat(negf, esteps, npoints, ldos)
    type(Tnegf), intent(in) :: negf
    integer, intent(in) :: esteps, npoints
    real(dp), dimension(:,:) :: ldos

    integer :: i, j

    if (allocated(negf%ldos_mat) .and. &
       & (esteps .eq. size(negf%ldos_mat,1)) .and. &
       & (npoints .eq. size(negf%ldos_mat,2))) then

       ldos = negf%ldos_mat
    end if

  end subroutine return_dos_mat


  !---- SAVE TUNNELING AND DOS ON FILES -----------------------------------------------
  ! GP Left in MPI version for debug purpose only. This will write a separate
  ! file for every ID, which is not possible on all architectures
  subroutine write_tunneling_and_dos(negf)

    type(Tnegf), intent(in) :: negf

    integer :: Nstep, i, i1, idos_proj, size_ni, iu
    character(6) :: ofKP, idstr
    real(dp) :: E

    if (allocated(negf%tunn_mat)) then

      Nstep = size(negf%tunn_mat,1)
      size_ni = size(negf%tunn_mat,2)

      write(ofKP,'(i6.6)') negf%ikpoint
      write(idstr,'(i6.6)') id

      open(newunit=iu,file=trim(negf%out_path)//'tunneling_'//ofKP//'_'//idstr//'.dat')

      !negf%eneconv=1.0_dp

      do i = 1,Nstep

        E=(negf%Emin+negf%Estep*(i-1))

        WRITE(iu,'(E17.8,20(E17.8))') E*negf%eneconv, &
            (negf%tunn_mat(i,i1), i1=1,size_ni)

      enddo

      close(iu)

    endif

    if (allocated(negf%ldos_mat) .and. negf%ndos_proj.gt.0) then

        Nstep = size(negf%ldos_mat,1)

        write(ofKP,'(i6.6)') negf%ikpoint
        write(idstr,'(i6.6)') id

        open(newunit=iu,file=trim(negf%out_path)//'localDOS_'//ofKP//'_'//idstr//'.dat')

        do i = 1,Nstep

          E=(negf%Emin+negf%Estep*(i-1))
          WRITE(iu,'((E17.8))',advance='NO') E*negf%eneconv
          do idos_proj = 1, negf%nDOS_proj
            WRITE(iu,'((E17.8))',advance='NO') negf%ldos_mat(i,idos_proj)/negf%eneconv
          end do
          write(iu,*)

        end do

        close(iu)

    endif

  end subroutine write_tunneling_and_dos

  !-------------------------------------------------------------------------------
  subroutine compute_phonon_current(negf)

    type(Tnegf) :: negf


    call extract_cont(negf)

    call tunneling_int_def(negf)

    call phonon_tunneling(negf)

    call phonon_current(negf)

    !!GP Locally writing energy dependent data is not meaningful in the MPI
    !implementation, because the gathering is done externally.
    ! An implementation node by node is still active, for debugging purposes
    !call write_tunneling_and_dos(negf)

    call destroy_contact_matrices(negf)

  end subroutine compute_phonon_current

  !---------------------------------------------------------------------------
  ! Sets the Reference contact for non-eq calculations
  !
  ! The behaviour depends on how negf%min_or_max has been set.
  !
  ! min_or_max = 0 : refcont is chosen at the minimum   mu
  ! min_or_max = 1 : refcont is chosen at the maximum   mu
  ! min_or_max = * : refcont is set to 0 (no referece)
  !
  subroutine set_ref_cont(negf)

    type(TNegf) :: negf

    integer :: nc_vec(1), ncont

    ncont = negf%str%num_conts

    if (ncont > 0) then
      if (negf%min_or_max .eq. 0) then
         negf%muref = minval(negf%cont(1:ncont)%mu)
         nc_vec = minloc(negf%cont(1:ncont)%mu)
      elseif (negf%min_or_max .eq. 1) then
         negf%muref = maxval(negf%cont(1:ncont)%mu)
         nc_vec = maxloc(negf%cont(1:ncont)%mu)
      else
         nc_vec = ncont + 1
         negf%muref = maxval(negf%cont(1:ncont)%mu)
      endif
      negf%refcont = nc_vec(1)
    else
      negf%muref = negf%mu
      negf%refcont = 0
    endif

  end subroutine set_ref_cont

  !> Print TNegf state, for debug
  subroutine print_tnegf(negf)
    type(TNegf) :: negf

    call print_all_vars(negf,6)
  end subroutine print_tnegf

  !////////////////////////////////////////////////////////////////////////
  ! RCM algorithm for reordering.
  !
  ! Actually not used because it is not suitable for contacted structures
  !////////////////////////////////////////////////////////////////////////
  subroutine reorder(mat)
    type(z_CSR) :: mat


    type(z_CSR) :: P, Tmp

    integer, dimension(:), allocatable :: perm
    integer :: i, nrow

    nrow=mat%nrow

    call log_allocate(perm,nrow)

    call genrcm(nrow, mat%nnz, mat%rowpnt, mat%colind, perm)

    call create(P,nrow,nrow,nrow)

    do i=1,nrow
       P%nzval(i)=1
       P%colind(i)=perm(i)
       P%rowpnt(i)=i
    enddo
    P%rowpnt(nrow+1)=nrow+1


    call create(Tmp,nrow,nrow,mat%nnz)

    call zamub_st(P,mat,Tmp)

    call ztransp_st(P)

    call zamub_st(Tmp,P,mat)

    call destroy(P,Tmp)

    call log_deallocate(perm)

  end subroutine reorder

  !>  Authomatic Block partitioning. The matrix must be already sorted.
  subroutine block_partition(mat,nrow,cont_end,surf_end,ncont,nbl,blks)

    !> The matrix to be partitioned.
    type(z_CSR), intent(in) :: mat
    !> The number of row to partition.
    integer, intent(in) :: nrow
    !> The indices indicating the end of the contact.
    integer, dimension(:), intent(in) :: cont_end
    !> The indices indicating the end of the scattering region surface
    !> (last orbitals before corresponding contact.)
    integer, dimension(:), intent(in) :: surf_end
    !> The number of contacts.
    integer, intent(in) :: ncont
    !> The number of blocks.
    integer, intent(out) :: nbl
    !> The array with the end index for each block.
    integer, dimension(:), intent(inout) :: blks

    integer :: j, k, i
    integer :: i1, i2

    integer :: rn, rnold, tmax, rmax, maxmax
    integer :: dbuff, minsize, minv, maxv

    minsize = 0
     do i1 = 1, ncont
        maxv = 0
        minv = 400000000
        do k = surf_end(i1)+1, cont_end(i1)
           do i = mat%rowpnt(k), mat%rowpnt(k+1)-1
            if (mat%colind(i).le.nrow .and.  mat%colind(i).lt.minv) minv = mat%colind(i)
            if (mat%colind(i).le.nrow .and.  mat%colind(i).gt.maxv) maxv = mat%colind(i)
           end do
        end do
        if (maxv-minv+1 .gt. minsize) minsize = maxv - minv + 1
    end do

    ! The current algorithm does not work when the minimum block
    ! size is 1. We fix the minimum possible size to 2 as temporary fix.
    minsize = max(minsize, 2)

    ! Find maximal stancil of the matrix and on which row
    !  ( Xx     )
    !  ( xXxx   )
    !  (  xXxxx )  <- maxmax = 3 ; rmax = 3
    !  (   xXx  )
    maxmax = 0
    do j=1,nrow
       tmax = 0
       do i = mat%rowpnt(j), mat%rowpnt(j+1) - 1
           if ( mat%colind(i).le.nrow .and. (mat%colind(i)-j) .gt. tmax) then
                tmax = mat%colind(i)-j
           endif
       enddo

       if(tmax .gt. maxmax) then
          maxmax = tmax
          rmax = j
       endif

       dbuff = maxmax        ! dbuff should be linked to maxmax
       minsize = max((dbuff+1)/2,minsize)
    enddo

    ! Define central block
    rn = rmax - maxmax/2 - dbuff

    if(rn-dbuff.ge.0) then

       blks(1) = rn-1  ! fine del blocco precedente

       nbl = 1

       do

          do j = rn, minsize, -1

             rnold = rn
             i1 = mat%rowpnt(j-minsize+1)
             i2 = mat%rowpnt(j+1) - 1

             !k = maxval(mat%colind(i1:i2))
             k = 0
             do i = i1, i2
                if ( mat%colind(i).le.nrow .and. (mat%colind(i)) .gt. k) then
                   k = mat%colind(i)
                endif
             enddo

             if(k.lt.rn) then
                rn = j
                nbl = nbl + 1
                if (nbl.gt.MAXNUMPLs) call errormsg()
                blks(nbl) = j-1 ! fine del blocco precedente
                exit
             endif
          enddo

          if(rn.le.minsize .or. rnold.eq.rn) then
             exit
          endif

       enddo

       rn = rmax - maxmax/2 - dbuff

    else
       nbl= 0
       rn = 1

    endif

    do

       do j = rn, nrow-minsize+1, 1

          rnold = rn
          i1 = mat%rowpnt(j)
          i2 = mat%rowpnt(j+minsize) - 1

          !k = minval(mat%colind(i1:i2))
          k = nrow
          do i = i1, i2
             if ( mat%colind(i).le.nrow .and. (mat%colind(i)) .lt. k) then
                k = mat%colind(i)
             endif
          enddo


          if(k.gt.rn) then
             rn = j
             nbl = nbl + 1
             if (nbl.gt.MAXNUMPLs) call errormsg()
             blks(nbl) = j-1 ! fine del blocco
             exit
          endif
       enddo

       if(nrow-rn.le.minsize .or. rnold.eq.rn) then
          exit
       endif
    enddo

    nbl = nbl + 1
    if (nbl.gt.MAXNUMPLs) call errormsg()
    blks(nbl) = nrow

    ! Sorting blocks

    do i = 1, nbl
       do j = i+1, nbl
          if(blks(j).lt.blks(i)) then
             k = blks(i)
             blks(i) = blks(j)
             blks(j) = k
          endif
       enddo
    enddo


  end subroutine block_partition

!----------------------------------------------------------------------------
  subroutine errormsg()

     write(*,*) "ERROR: Maximum number of PLs exceeded"
     write(*,*) "increase the value of MAXNUMPLS in libnegf.F90"
     write(*,*) "and recompile the library"

     STOP !should rise an exception

  end subroutine errormsg

!----------------------------------------------------------------------------
  subroutine find_cblocks(mat ,ncont, nbl, PL_end, surf_start, cont_end, cblk)
    type(z_CSR), intent(in) :: mat
    integer, intent(in) :: ncont
    integer, intent(in) :: nbl
    integer, dimension(:), intent(in) :: PL_end
    integer, dimension(:), intent(in) :: surf_start
    integer, dimension(:), intent(in) :: cont_end
    integer, dimension(:), allocatable, intent(out) :: cblk

    integer :: j1,k,i,min,max
    integer, dimension(:), allocatable :: PL_start

    call log_allocate(PL_start,nbl)
    call log_allocate(cblk,ncont)

    PL_start(1) = 1

    do i = 2, nbl
       PL_start(i) = PL_end(i-1) + 1
    enddo


    do j1 = 1, ncont

       max = 0
       min = 400000000

       do k = surf_start(j1), cont_end(j1)

          do i = mat%rowpnt(k), mat%rowpnt(k+1)-1

             if (mat%colind(i).le.PL_end(nbl) .and.  mat%colind(i).lt.min) min = mat%colind(i)
             if (mat%colind(i).le.PL_end(nbl) .and.  mat%colind(i).gt.max) max = mat%colind(i)

          end do

       end do

       do k = 1, nbl

          if( max .le. PL_end(k) ) then
             cblk(j1) = k

             if( min .ge. PL_start(k) ) then
                exit
             else
                write(*,*) "(LibNEGF) Partitioning:"
                write(*,*) "Number of blocks: ",nbl
                write(*,*) "PL_end: ",PL_end(1:nbl)
                write(*,*) "Contact interaction: ",cblk(j1)
                write(*,'(a,i3,a)') " ERROR: contact",j1," interacting with more than one block"
                write(*,*) "min ",min,"max ",max
                stop
             end if

          end if

       end do

    end do

    call log_deallocate(PL_start)

  end subroutine find_cblocks

!----------------------------------------------------------------------------

  Subroutine sort(blks, Ipt)
    ! *
    ! ***********************************
    ! * Sort Array X(:) in ascendent order.
    ! * If present Ipt, a pointer with the
    ! * changes is returned in Ipt.
    ! ***********************************

    Integer, Intent (inout) :: blks(:)
    Integer, Intent (out), Optional :: Ipt(:)

    Integer :: Rtmp
    Integer :: i, j

    If (Present(Ipt)) Then
       Forall (i=1:Size(blks)) Ipt(i) = i

       Do i = 2, Size(blks)
          Rtmp = blks(i)
          Do j = i-1, 1, -1
             If (Rtmp < blks(j)) Then
                blks(j+1) = blks(j)
                call Swap(blks, j, j+1)
             Else
                Exit
             End If
          End Do
          blks(j+1) = Rtmp
       End Do
    Else
       Do i = 2, Size(blks)
          Rtmp = blks(i)
          Do j = i-1, 1, -1
             If (Rtmp < blks(j)) Then
                blks(j+1) = blks(j)
             Else
                Exit
             End If
          End Do
          blks(j+1) = Rtmp
       End Do
    End If

    Return
  End Subroutine sort

  ! ***********************************
  ! * Swaps elements I and J of array X(:).
  ! ***********************************
  Subroutine Swap(X, i, j)

    Integer, Intent (inout) :: X(:)
    Integer, Intent (in) :: i, j

    Integer :: Itmp

    Itmp = X(I)
    X(I) = X(J)
    X(J) = Itmp

    Return
  End Subroutine Swap

  !------------------------------------------
  subroutine printcsr(id, mat)
    integer :: id
    type(z_csr) :: mat

    call zprint_csrcoo(id, mat, 'c')
  end subroutine printcsr

  subroutine printcsrij(id, mat, i, j)
    integer :: id
    type(z_csr) :: mat
    integer :: i, j

    write(id,*) i,j,getelement(i,j,mat)

  end subroutine printcsrij

  function getel(mat,i,j) result(el)
    type(z_csr) :: mat
    integer :: i, j
    complex(dp) :: el

    el = getelement(i,j,mat)

  end function getel
  
  subroutine aggregate_vec(v_in, thres, v_out, start_idx, end_idx)
      real(dp), dimension(:), intent(in) :: v_in
      real(dp), intent(in) :: thres
      real(dp), dimension(:), allocatable, intent(out) :: v_out
      integer, dimension(:), allocatable, intent(out) :: start_idx, end_idx

      real(dp) :: avg
      integer :: i, rs, re

      allocate(start_idx(0))
      allocate(end_idx(0))
      allocate(v_out(0))

      start_idx = [start_idx, 1]
      do i = 1, size(v_in)-1
         if (abs(v_in(i+1) - v_in(i)) > thres) then
            start_idx = [start_idx, i+1]
         endif
      end do

      do i = 1, size(v_in)-1
         if (abs(v_in(i+1) - v_in(i)) > thres) then
                 end_idx = [end_idx, i]
         endif
      end do
      end_idx = [end_idx, size(v_in)]

      do i = 1, size(start_idx)
         rs = start_idx(i)
         re = end_idx(i)
         avg = sum(v_in(rs:re))/real(size(v_in(rs:re)), dp)
         v_out = [v_out, avg]
      end do

   end subroutine aggregate_vec

end module libnegf
