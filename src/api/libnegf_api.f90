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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Fortran 77 style subroutines for communication with LIBNEGF
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! The order you should invoke the subroutines during initialisation is the
!! following:
!!
!!   negf_init_session
!!   negf_init                (after negf_init_session)
!!                             See in libnegf.F90
!!
!! After initialisation you can call the following methods:
!!
!!   negf_set_verbosity
!!   negf_get_version
!!   negf_set_iteration       (set scf iteration, default 1)
!!   negf_set_outer           (outer blocks of DM, default 2)
!!   negf_set_kpoint          (integer value for file handling)
!!   negf_set_reference       (sets ref. contact, default mu_max)
!!   negf_set_output          (set output path)
!!   negf_set_scratch         (set scratch path)
!!   negf_set_writetunn       (writes tunneling files)
!!   negf_set_writeldos       (writes ldos files)
!!   negf_write_partition     (writes H partitioning)
!!   negf_density             (computes density matrix)
!!   negf_current             (computes tunneling and current)
!!
!! In order to destroy a UPT instance call
!!
!!   negf_destruct_libnegf      (clean the data container)
!!   negf_destruct_session      (after upt_destruct_libnegf)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!>
!! Returns the size of the handler array
!! @param[out]  handlerSize Contains the size of the handler array on exit.
subroutine negf_gethandlersize(handlerSize) bind(C)
  use iso_c_binding, only : c_int  ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  implicit none
  integer(c_int), intent(out) :: handlerSize  ! if:var:out

  handlerSize = DAC_handlerSize

end subroutine negf_gethandlersize

!>
!! Initialises a new LIBNEGF instance
!! @param [out]  handler  Contains the handler for the new instance on return
subroutine negf_init_session(handler) bind(C)
  use iso_c_binding, only : c_int  ! if:mod:use
  use libnegfAPICommon ! if:mod:use
  implicit none
  integer(c_int) :: handler(DAC_handlerSize) ! if:var:out

  !type(TNEGF), pointer :: pNEGF
  type(NEGFpointers) :: LIB

  IF( size(transfer(LIB, handler)) > size(handler) ) then
     write(*,*) size(transfer(LIB, handler)), size(handler)
     stop 'Handler size mismatch'
  ENDIF

  NULLIFY(LIB%pNEGF)
  ALLOCATE(LIB%pNEGF)

  handler(:) = 0
  handler = transfer(LIB, handler, size(handler))

  ! call here a specific initialization method

end subroutine negf_init_session

!>
!! Get library version. Reads svn versions=> broken with git
!! @param  handler  Contains the handler for the new instance on return
subroutine negf_get_version(handler) bind(C)
  use iso_c_binding, only : c_int  ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf  ! if:mod:use
  implicit none
  integer(c_int) :: handler(DAC_handlerSize)  ! if:var:inout

  !type(TNEGF), pointer :: pNEGF
  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)

  call negf_version(LIB%pNEGF)

end subroutine negf_get_version

!!* Initialises a new LIBNEGF instance
!!* @param  handler  Contains the handler for the new instance on return
subroutine negf_init(handler) bind(C)
  use iso_c_binding, only : c_int  ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  implicit none
  integer(c_int), intent(inout) :: handler(DAC_handlerSize)  ! if:var:inout

  !type(TNEGF), pointer :: pNEGF
  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)

  call init_negf(LIB%pNEGF)

end subroutine negf_init

!!* Pass structure
subroutine negf_init_structure(handler, ncont, contend, surfend, npl, plend, cblk) bind(c)
  use iso_c_binding, only : c_int  ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  implicit none
  integer(c_int), intent(inout) :: handler(DAC_handlerSize)  ! if:var:inout
  integer(c_int), intent(in), value :: ncont ! if:var:inout
  integer(c_int), intent(in), value :: npl ! if:var:inout
  integer(c_int), intent(in) :: surfend(*) ! if:var:inout
  integer(c_int), intent(in) :: contend(*) ! if:var:inout
  integer(c_int), intent(in) :: plend(*)   ! if:var:inout
  integer(c_int), intent(in) :: cblk(*) ! if:var:inout

  integer, allocatable :: surfend_al(:), contend_al(:), plend_al(:), cblk_al(:)
  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)

  allocate(surfend_al(ncont))
  allocate(contend_al(ncont))
  allocate(cblk_al(ncont))
  allocate(plend_al(npl))

  surfend_al(1:ncont) = surfend(1:ncont)
  contend_al(1:ncont) = contend(1:ncont)
  cblk_al(1:ncont) = cblk(1:ncont)
  plend_al(1:npl) = plend(1:npl)

  call init_structure(LIB%pNEGF, ncont, contend_al, surfend_al, npl, plend_al, cblk_al)

end subroutine negf_init_structure

!!* Pass contacts
subroutine negf_init_contacts(handler, ncont) bind(c)
  use iso_c_binding, only : c_int  ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  implicit none
  integer(c_int), intent(inout) :: handler(DAC_handlerSize)  ! if:var:inout
  integer(c_int), intent(in), value :: ncont ! if:var:inout
  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)

  call init_contacts(LIB%pNEGF, ncont)

end subroutine negf_init_contacts

!!* Passing parameters
subroutine negf_set_params(handler, params) bind(c)
  use iso_c_binding, only : c_int  ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  implicit none
  integer(c_int), intent(inout) :: handler(DAC_handlerSize)  ! if:var:inout
  type(lnparams), intent(in) :: params ! if:var:inout

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
  call set_params(LIB%pNEGF, params)
end subroutine negf_set_params

!!* Getting parameters
subroutine negf_get_params(handler, params) bind(c)
  use iso_c_binding, only : c_int  ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  implicit none
  integer(c_int), intent(inout) :: handler(DAC_handlerSize)  ! if:var:inout
  type(lnparams), intent(inout) :: params ! if:var:inout

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
  call get_params(LIB%pNEGF, params)
end subroutine negf_get_params

!!* Passing Hamiltonian from memory
!!* @param  handler  Contains the handler for the new instance on return
subroutine negf_set_h(handler, nrow, A, JA, IA)
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  use ln_precision      ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer, intent(in), value :: nrow     ! if:var:in
  complex(dp) :: A(*) ! if:var:in
  integer :: JA(*)    ! if:var:in
  integer :: IA(*)    ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
  call set_H(LIB%pNEGF,nrow, A, JA, IA)
end subroutine negf_set_h

!!* Passing Overlap from memory
!!* @param  handler  Contains the handler for the new instance on return
subroutine negf_set_s(handler, nrow, A, JA, IA)
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  use ln_precision      ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer, intent(in), value :: nrow     ! if:var:out
  complex(dp) :: A(*) ! if:var:out
  integer :: JA(*)    ! if:var:out
  integer :: IA(*)    ! if:var:out

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
  call set_S(LIB%pNEGF,nrow, A, JA, IA)
end subroutine negf_set_s

!> Same as negf_set_s, but pass separately real
!! and imaginary part (for ctypes c-python interface)
subroutine negf_set_s_cp(handler, nrow, a_re, a_im, ja, ia) bind(c)
  use iso_c_binding, only : c_int, c_double
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  use ln_precision      ! if:mod:use
  implicit none
  integer(c_int) :: handler(DAC_handlerSize)  ! if:var:in
  integer(c_int), intent(in), value :: nrow     ! if:var:out
  real(c_double) :: a_re(*) ! if:var:out
  real(c_double) :: a_im(*) ! if:var:out
  integer(c_int) :: ja(*)    ! if:var:out
  integer(c_int) :: ia(*)    ! if:var:out

  type(NEGFpointers) :: LIB
  integer :: nnz, ii
  complex(dp), allocatable :: A(:)

  LIB = transfer(handler, LIB)
  nnz = IA(nrow+1) - IA(1)
  allocate(A(nnz))
  do ii = 1, nnz
    A(ii) = cmplx(A_re(ii), A_im(ii), kind=kind(1.0d0))
  end do
  call set_S(LIB%pNEGF,nrow, A, JA, IA)
end subroutine negf_set_s_cp

!> Same as negf_set_h, but pass separately real
!! and imaginary part (for ctypes c-python interface)
subroutine negf_set_h_cp(handler, nrow, a_re, a_im, ja, ia) bind(c)
  use iso_c_binding, only : c_int, c_double
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  use ln_precision      ! if:mod:use
  implicit none
  integer(c_int) :: handler(DAC_handlerSize)  ! if:var:in
  integer(c_int), intent(in), value :: nrow     ! if:var:out
  real(c_double) :: a_re(*) ! if:var:out
  real(c_double) :: a_im(*) ! if:var:out
  integer(c_int) :: ja(*)    ! if:var:out
  integer(c_int) :: ia(*)    ! if:var:out

  type(NEGFpointers) :: LIB
  integer :: nnz, ii
  complex(dp), allocatable :: A(:)

  LIB = transfer(handler, LIB)
  nnz = ia(nrow+1) - ia(1)
  allocate(A(nnz))
  do ii = 1,nnz
    A(ii) = cmplx(a_re(ii), a_im(ii), kind=kind(1.0d0))
  end do
  call set_H(LIB%pNEGF,nrow, A, JA, IA)
end subroutine negf_set_h_cp

!!* Passing Overlap from memory
!!* @param  handler  Contains the handler for the new instance on return
subroutine negf_set_s_id(handler, nrow) bind(C)
  use iso_c_binding, only : c_int
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  use ln_precision      ! if:mod:use
  implicit none
  integer(c_int) :: handler(DAC_handlerSize)  ! if:var:in
  integer(c_int), intent(in), value:: nrow     ! if:var:out

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)

  call set_S_id(LIB%pNEGF,nrow)

end subroutine negf_set_s_id

subroutine negf_print_mat(handler) bind(C)
  use iso_c_binding, only : c_int  ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  implicit none
  integer(c_int), intent(in) :: handler(DAC_handlerSize)  ! if:var:in

  !type(TNEGF), pointer :: pNEGF
  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)

  !print*,"(libNEGF) write "//trim(LIB%pNEGF%scratch_path)//"H.dat"
  open(666, file=trim(LIB%pNEGF%scratch_path)//'H.dat')
  call printcsr(666,LIB%pNEGF%H)
  close(666)

  !print*,"(libNEGF) write "//trim(LIB%pNEGF%scratch_path)//"S.dat"
  open(666, file=trim(LIB%pNEGF%scratch_path)//'S.dat')
  call printcsr(666,LIB%pNEGF%S)
  close(666)

end subroutine negf_print_mat

!!* Fill parameters from input file negf.in
!!* @param  handler  Contains the handler for the new instance on return
subroutine negf_read_input(handler) bind(C)
  use iso_c_binding, only : c_int  ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  implicit none
  integer(c_int), intent(in) :: handler(DAC_handlerSize)  ! if:var:in

  !type(TNEGF), pointer :: pNEGF
  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
  call read_negf_in(LIB%pNEGF)

end subroutine negf_read_input

!!* Fill parameters from input file negf.in
!!* @param  handler  Contains the handler for the new instance on return
  !TODO: We can (and maybe should) wrap string in input to accept c-like
subroutine negf_read_hs(handler, real_path, imag_path, target_matrix) bind(C)
  use iso_c_binding, only : c_int, c_char, c_null_char  ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  use globals           ! if:mod:use
  implicit none
  integer(c_int), intent(in) :: handler(DAC_handlerSize)  ! if:var:in
  character(kind=c_char), intent(in) :: real_path(*) ! if:var:in
  character(kind=c_char), intent(in) :: imag_path(*) ! if:var:in
  integer(c_int), intent(in), value :: target_matrix !if:var:in

  character(len=LST) :: freal_path ! if:var:in
  character(len=LST) :: fimag_path ! if:var:in
  integer :: nn, n_char
  !type(TNEGF), pointer :: pNEGF
  type(NEGFpointers) :: LIB

  !Padding and converting to fortran string
  n_char = 1
  do while(real_path(n_char).ne.c_null_char)
    n_char = n_char + 1
  end do
  freal_path = " "
  do nn=1,n_char-1
    freal_path(nn:nn) = real_path(nn)
  end do
  n_char = 1
  do while(imag_path(n_char).ne.c_null_char)
    n_char = n_char + 1
  end do
  fimag_path = " "
  do nn=1,n_char-1
    fimag_path(nn:nn) = imag_path(nn)
  end do

  LIB = transfer(handler, LIB)

  call read_HS(LIB%pNEGF, freal_path, fimag_path, target_matrix)

end subroutine negf_read_hs


!!* Destroys a certain LIBNEGF instance
!!* @param  handler  Handler for the instance to destroy
subroutine negf_destruct_session(handler) bind(C)
  use iso_c_binding, only : c_int  ! if:mod:use
  use libnegfAPICommon             ! if:mod:use
  implicit none
  integer(c_int), intent(in) :: handler(DAC_handlerSize)  ! if:var:in

  type(NEGFpointers) :: LIB
  integer :: err

  LIB = transfer(handler, LIB)

  if (associated(LIB%pNEGF)) deallocate(LIB%pNEGF, stat= err)

  if (err.ne.0) write(*,*) '(negf_destructsession) Deallocation error'

end subroutine negf_destruct_session


!!* Clean the data containers of a given LIBNEGF instance.
!!* @param handler Number for the LIBNEGF instance to destroy.
subroutine negf_destruct_libnegf(handler) bind(C)
  use iso_c_binding, only : c_int  ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf  ! if:mod:use
  implicit none
  integer(c_int), intent(in) :: handler(DAC_handlerSize)  ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
  call destroy_negf(LIB%pNEGF)

end subroutine negf_destruct_libnegf

!!* Destructs a given LIBNEGF instance.
!!* @param handler Number for the LIBNEGF instance to destroy.
subroutine negf_set_verbosity(handler,verbose_lev) bind(C)
  use iso_c_binding, only : c_int  ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  implicit none
  integer(c_int), intent(in) :: handler(DAC_handlerSize)  ! if:var:in
  integer(c_int), intent(in) :: verbose_lev               ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)

  LIB%pNEGF%verbose = verbose_lev

end subroutine negf_set_verbosity

!>
!! Solve the Landauer problem: calculate transmission and
!! density of states according to previously specified parameters
!! @param[in]  handler: handler Number for the LIBNEGF instance
subroutine negf_solve_landauer(handler) bind(C)
  use iso_c_binding, only : c_int  ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf   ! if:mod:use
  implicit none
  integer(c_int), intent(in) :: handler(DAC_handlerSize)  ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
  call compute_current(LIB%pNEGF)
end subroutine negf_solve_landauer

!>
!! Calculate the density matrix for the dft problem
!! @param[in]  handler: handler Number for the LIBNEGF instance
subroutine negf_solve_density_dft(handler) bind(C)
  use iso_c_binding, only : c_int  ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf   ! if:mod:use
  implicit none
  integer(c_int), intent(in) :: handler(DAC_handlerSize)  ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
  call compute_density_dft(LIB%pNEGF)
end subroutine negf_solve_density_dft


!>
!! Get current value for a specific couple of leads
!!  @param[in] handler: handler Number for the LIBNEGF instance
!! @param [in] lead_pair: specifies which leads are considered
!!             for retrieving current (as ordered in ni, nf)
!!  @param[in] unitoOfH: units modifer for energy (write"unknown"
!!                   for default)
!!  @param[in] unitOfJ: units modifier for current (write"unknown"
!!                  for default)
!!  @param[out] current: current value
subroutine negf_get_current(handler, leadPair, unitOfH, unitOfJ, current)
  use libnegfAPICommon  ! if:mod:use
  use libnegf   ! if:mod:use
  use ln_constants ! if:mod:use
  use globals  ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: leadPair      !if:var:in
  real(dp) :: current       !if:var:inout
  character(SST) :: unitOfH !if:var:in
  character(SST) :: unitOfJ !if:var:in

  type(NEGFpointers) :: LIB
  type(unit) :: unitH, unitJ

  unitH%name=trim(unitOfH)
  unitJ%name=trim(unitOfJ)
  LIB = transfer(handler, LIB)
  current = LIB%pNEGF%currents(leadPair) ! just take first value (2 contacts)
  ! units conversion.
  current = current * convertCurrent(unitH, unitJ)
end subroutine negf_get_current

!> Copy the energy axis on all processors (for output, plot, debug)
!! Uses a fixed size array interface
!! @param [in] handler:  handler Number for the LIBNEGF instance
!! @param [out] npoints (int):
!! @param [out] re_en (array): real part of energy values
!! @param [out] im_en (array): imaginary part of energy values
!! @param [in] copy (int): 0 if you want only to fill npoints, any value
!!               if you want to perform the actual copy
subroutine negf_get_energies(handler, npoints, re_en, im_en, copy) bind(c)
  use iso_c_binding, only : c_int, c_double, c_double_complex  ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf   ! if:mod:use
  implicit none
  integer(c_int) :: handler(DAC_handlerSize)  ! if:var:in
  integer(c_int), intent(out) ::npoints ! if:var:in
  real(c_double), intent(out) :: re_en(*)  ! if:var:in
  real(c_double), intent(out) :: im_en(*)  ! if:var:in
  integer(c_int), intent(in), value :: copy ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
  if (copy.eq.0) then
    npoints = size(LIB%pNEGF%en_grid)
  else
    npoints = size(LIB%pNEGF%en_grid)
    re_en(1:npoints) = real(LIB%pNEGF%en_grid(1:npoints)%Ec)
    im_en(1:npoints) = aimag(LIB%pNEGF%en_grid(1:npoints)%Ec)
  end if
end subroutine negf_get_energies

!> Copy the energy axis on all processors (for output, plot, debug)
!! Uses a fixed size array interface
!! @param [in] handler:  handler Number for the LIBNEGF instance
!! @param [out] nnz (int): number of non zero values
!! @param [out] nrow (int): number of rows
!! @param [out] rowpnt (int array): row pointer indexes, size nrow+1
!! @param [out] colind (int array): column indexes array, size nnz
!! @param [out] re_nzval (double array): non zero values (real part), size nnz
!! @param [out] im_nzval (double array): non zero values (imag part), size nnz
!! @param [in] copy (int): 0 if you want only to fill npoints, any value
!!               if you want to perform the actual copy
subroutine negf_get_dm(handler, nnz, nrow, rowpnt, colind, re_nzval, im_nzval, copy) bind(c)
  use iso_c_binding, only : c_int, c_double, c_double_complex  ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf   ! if:mod:use
  implicit none
  integer(c_int) :: handler(DAC_handlerSize)  ! if:var:in
  integer(c_int), intent(out) ::nnz ! if:var:in
  integer(c_int), intent(out) ::nrow ! if:var:in
  integer(c_int), intent(out) :: rowpnt(*)  ! if:var:in
  integer(c_int), intent(out) :: colind(*)  ! if:var:in
  real(c_double), intent(out) :: re_nzval(*)  ! if:var:in
  real(c_double), intent(out) :: im_nzval(*)  ! if:var:in
  integer(c_int), intent(in), value :: copy ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
  nnz = LIB%pNEGF%rho%nnz
  nrow = LIB%pNEGF%rho%nrow
  if (.not.(copy.eq.0)) then
    rowpnt(1:nrow+1) = LIB%pNEGF%rho%rowpnt(1:nrow+1)
    colind(1:nnz) = LIB%pNEGF%rho%colind(1:nnz)
    re_nzval(1:nnz) = real(LIB%pNEGF%rho%nzval(1:nnz))
    im_nzval(1:nnz) = aimag(LIB%pNEGF%rho%nzval(1:nnz))
  end if
end subroutine negf_get_dm

!> Copy the currents values
!! Uses a fixed size array interface
!! @param [in] handler:  handler Number for the LIBNEGF instance
!! @param [out] npairs (int): number of active lead pairs (# current values)
!! @param [out] currents (array): real part of energy values
!! @param [in] copy (int): 0 if you want only to fill npoints, any value
!!               if you want to perform the actual copy
subroutine negf_get_currents(handler, npairs, currents, copy) bind(c)
  use iso_c_binding, only : c_int, c_double, c_double_complex  ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf   ! if:mod:use
  implicit none
  integer(c_int) :: handler(DAC_handlerSize)  ! if:var:in
  integer(c_int), intent(out) :: npairs ! if:var:in
  real(c_double), intent(out) :: currents(*)  ! if:var:in
  integer(c_int), intent(in), value :: copy ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
  if (copy.eq.0) then
    npairs = size(LIB%pNEGF%currents)
  else
    npairs = size(LIB%pNEGF%currents)
    currents(1:npairs) = LIB%pNEGF%currents(1:npairs)
  end if
end subroutine negf_get_currents

!> Pass pointer to transmission output to a compatible C pointer
!!  @param[in]  handler:  handler Number for the LIBNEGF instance
!!  @param[out] tr_shape: shape of transmission n-array (in fortran)
!!  @param[out] tr_pointer: C pointer to data
subroutine negf_associate_transmission(handler, tr_shape, tr_pointer) bind(c)
  use iso_c_binding, only : c_int, c_double, c_loc, c_ptr   ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf   ! if:mod:use
  implicit none
  integer(c_int) :: handler(DAC_handlerSize)  ! if:var:in
  integer(c_int), intent(out) :: tr_shape(2) ! if:var:out
  type(c_ptr), intent(out) :: tr_pointer ! if:var:out

  type(NEGFpointers) :: LIB
  real(c_double), dimension(:,:), pointer :: f_p

  LIB = transfer(handler, LIB)
  call associate_transmission(LIB%pNEGF, f_p)
  tr_pointer = c_loc(f_p(1,1))
  tr_shape = shape(f_p)

end subroutine negf_associate_transmission

!> Pass pointer to transmission output to a compatible C pointer
!!  @param[in]  handler:  handler Number for the LIBNEGF instance
!!  @param[out] ldos_shape: shape of ldos n-array (in fortran)
!!  @param[out] ldos_pointer: C pointer to data
subroutine negf_associate_ldos(handler, ldos_shape, ldos_pointer) bind(c)
  use iso_c_binding, only : c_int, c_double, c_loc, c_ptr   ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf   ! if:mod:use
  implicit none
  integer(c_int) :: handler(DAC_handlerSize)  ! if:var:in
  integer(c_int), intent(out) :: ldos_shape(2) ! if:var:out
  type(c_ptr), intent(out) :: ldos_pointer ! if:var:out

  type(NEGFpointers) :: LIB
  real(c_double), dimension(:,:), pointer :: f_p

  LIB = transfer(handler, LIB)
  call associate_ldos(LIB%pNEGF, f_p)
  ldos_pointer = c_loc(f_p(1,1))
  ldos_shape = shape(f_p)

end subroutine negf_associate_ldos

!!> Set ldos intervals
!! @param [in] handler: handler Number for the LIBNEGF instance
!! @param [in] istart(nldos) array with first interval index
!! @param [in] iend(nldos) array with first interval index
subroutine negf_set_ldos_intervals(handler, nldos, istart, iend) bind(c)
  use iso_c_binding, only : c_int ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf   ! if:mod:use
  implicit none
  integer(c_int) :: handler(DAC_handlerSize)  ! if:var:in
  integer(c_int), intent(in), value :: nldos ! if:var:in
  integer(c_int), intent(in) :: istart(*)  ! if:var:in
  integer(c_int), intent(in) :: iend(*)  ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
  call set_ldos_intervals(LIB%pNEGF, nldos, istart, iend)
end subroutine negf_set_ldos_intervals

!> Initialize the ldos container
!! @param [in] handler: handler Number for the LIBNEGF instance
!! @param [in] nldos: number of intervals
subroutine negf_init_ldos(handler, nldos) bind(c)
  use iso_c_binding, only : c_int ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf   ! if:mod:use
  implicit none
  integer(c_int) :: handler(DAC_handlerSize)  ! if:var:in
  integer(c_int), intent(in), value :: nldos ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
  call init_ldos(LIB%pNEGF, nldos)
end subroutine negf_init_ldos

!> Set ldos indexes arrays for a given ldos
!!
!! @param [in] handler: handler Number for the LIBNEGF instance
!! @param [in] ildos: index of ldos (fortran indexing)
!! @param [in] idx_size: size of index array
!! @param [in] idx: array with indexes
subroutine negf_set_ldos_indexes(handler, ildos, idx_size, idx) bind(c)
  use iso_c_binding, only : c_int ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf   ! if:mod:use
  implicit none
  integer(c_int) :: handler(DAC_handlerSize)  ! if:var:in
  integer(c_int), intent(in), value :: ildos ! if:var:in
  integer(c_int), intent(in) ::idx_size ! if:var:in
  integer(c_int), intent(in) :: idx(*)  ! if:var:in

  type(NEGFpointers) :: LIB
  integer(c_int), allocatable :: idx_tmp(:)

  LIB = transfer(handler, LIB)
  allocate(idx_tmp(idx_size))
  idx_tmp(1:idx_size) = idx(1:idx_size)
  call set_ldos_indexes(LIB%pNEGF, ildos, idx_tmp)
end subroutine negf_set_ldos_indexes

!>
!!  Write tunneling and density of states (if any) to file
!!  @param
!!*        handler:  handler Number for the LIBNEGF instance
!!*        path: string specifying the output file
!!  NOTE: when running parallel code every node will write a
!!  separate bunch of energy points. You should implement I/O
!!  OUT of the library and use this routine for testing/debugging
!!  ONLY
subroutine negf_write_tunneling_and_dos(handler) bind(C)
  use iso_c_binding, only : c_int   ! if:mod:use
  use libnegfAPICommon ! if:mod:use
  use libnegf          ! if:mode:use
  integer(c_int) :: handler(DAC_handlerSize)  ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
  call write_tunneling_and_dos(LIB%pNEGF)

end subroutine negf_write_tunneling_and_dos


!!* Sets iteration in self-consistent loops
!!* @param handler Number for the LIBNEGF instance to destroy.
!! SAME, IS IT REALLY NEEDED??
subroutine negf_set_output(handler, out_path)
  use libnegfAPICommon    ! if:mod:use
  use globals             ! if:mod:use
  use libnegf             ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  character(LST) :: out_path(1)    ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)

  LIB%pNEGF%out_path = trim( out_path(1) ) // '/'

end subroutine negf_set_output


!!* Sets iteration in self-consistent loops
!!* @param handler Number for the LIBNEGF instance to destroy.
subroutine negf_set_scratch(handler, scratch_path)
  use libnegfAPICommon    ! if:mod:use
  use globals             ! if:mod:use
  use libnegf             ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  character(LST) :: scratch_path(1)    ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)

  LIB%pNEGF%scratch_path = trim( scratch_path(1) ) // '/'

end subroutine negf_set_scratch

!!* Instructs whether the device-contact part of the Density-Mat
!!* needs to be computed (Relevant for charges in non-orthogonal bases)
!!* outer = 0   No calculations
!!* outer = 1   Computes upper diagonal blocks
!!* outer = 2   Computes both upper and lower blocks (needed by dftb+)
!!* @param handler Number for the LIBNEGF instance to destroy.
subroutine negf_set_outer(handler, outer)
  use libnegfAPICommon    ! if:mod:use
  use libnegf             ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: outer              ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)

  LIB%pNEGF%outer = outer

end subroutine negf_set_outer


subroutine negf_set_kpoint(handler, kpoint)
  use libnegfAPICommon    ! if:mod:use
  use libnegf             ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: kpoint             ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)

  LIB%pNEGF%kpoint = kpoint

end subroutine negf_set_kpoint


subroutine negf_set_reference(handler, minmax)
  use libnegfAPICommon    ! if:mod:use
  use libnegf             ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: minmax             ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)

  LIB%pNEGF%min_or_max = minmax

end subroutine negf_set_reference


subroutine negf_write_partition(handler)
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)

  call negf_partition_info(LIB%pNEGF)

end subroutine negf_write_partition

!>
!!* Compute current for a given LIBNEGF instance.
!!* @param [in] handler Number for the LIBNEGF instance to destroy.
subroutine negf_current(handler, current, unitOfH, unitOfJ)
  use libnegfAPICommon  ! if:mod:use
  use libnegf   ! if:mod:use
  use ln_constants ! if:mod:use
  use globals  ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  real(dp) :: current       !if:var:inout
  character(SST) :: unitOfH !if:var:in
  character(SST) :: unitOfJ !if:var:in

  type(NEGFpointers) :: LIB
  type(unit) :: unitH, unitJ

  unitH%name=trim(unitOfH)
  unitJ%name=trim(unitOfJ)

  LIB = transfer(handler, LIB)

  call compute_current(LIB%pNEGF)

  current = LIB%pNEGF%currents(1) ! just take first value (2 contacts)

  ! units conversion.
  current = current * convertCurrent(unitH, unitJ)

  call write_tunneling_and_dos(LIB%pNEGF)

end subroutine negf_current

!> Print TNegf container for debug
subroutine negf_print_tnegf(handler) bind(c)
  use iso_c_binding, only : c_int   ! if:mod:use
  use libnegfAPICommon  ! if:mod:use
  use libnegf   ! if:mod:use
  implicit none
  integer(c_int) :: handler(DAC_handlerSize)  ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
  call print_tnegf(LIB%pNEGF)

end subroutine negf_print_tnegf
