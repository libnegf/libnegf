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
subroutine negf_gethandlersize(handlerSize)
  use libnegfAPICommon  ! if:mod:use
  implicit none
  integer :: handlerSize  ! if:var:out

  handlerSize = DAC_handlerSize

end subroutine negf_gethandlersize

!>
!! Initialises a new LIBNEGF instance
!! @param [out]  handler  Contains the handler for the new instance on return
subroutine negf_init_session(handler)
  use libnegfAPICommon  ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:out
  
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
subroutine negf_get_version(handler)
  use libnegfAPICommon  ! if:mod:use
  use libnegf  ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:inout
  
  !type(TNEGF), pointer :: pNEGF
  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
 
  call negf_version(LIB%pNEGF)

end subroutine negf_get_version

!!* Initialises a new LIBNEGF instance
!!* @param  handler  Contains the handler for the new instance on return
subroutine negf_init(handler)
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:inout

  !type(TNEGF), pointer :: pNEGF
  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
 
  call init_negf(LIB%pNEGF)
  
end subroutine negf_init

!!* Passing Hamiltonian from memory 
!!* @param  handler  Contains the handler for the new instance on return
subroutine negf_set_h(handler, nrow, A, JA, IA)
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  use ln_precision      ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:inout
  integer :: nrow     ! if:var:out
  complex(dp) :: A(*) ! if:var:out 
  integer :: JA(*)    ! if:var:out
  integer :: IA(*)    ! if:var:out

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
  integer :: handler(DAC_handlerSize)  ! if:var:inout
  integer :: nrow     ! if:var:out
  complex(dp) :: A(*) ! if:var:out 
  integer :: JA(*)    ! if:var:out
  integer :: IA(*)    ! if:var:out

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
 
  call set_S(LIB%pNEGF,nrow, A, JA, IA)
  
end subroutine negf_set_s 

!!* Passing Overlap from memory                        
!!* @param  handler  Contains the handler for the new instance on return
subroutine negf_set_s_id(handler, nrow)
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  use ln_precision      ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:inout
  integer :: nrow     ! if:var:out

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
 
  call set_S_id(LIB%pNEGF,nrow)
  
end subroutine negf_set_s_id


subroutine negf_print_mat(handler)
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:inout

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
subroutine negf_read_input(handler)
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:inout

  !type(TNEGF), pointer :: pNEGF
  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
  !print*,'(libNEGF) reading negf.in ...' 
  call read_negf_in(LIB%pNEGF)
  
end subroutine negf_read_input

!!* Fill parameters from input file negf.in
!!* @param  handler  Contains the handler for the new instance on return
subroutine negf_read_hs(handler, real_path, imag_path, target_matrix)
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  use globals           ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:inout
  character(LST) :: real_path ! if:var:in
  character(LST) :: imag_path ! if:var:in
  integer :: target_matrix !if:var:in

  !type(TNEGF), pointer :: pNEGF
  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
 
  call read_HS(LIB%pNEGF, real_path, imag_path, target_matrix)
  
end subroutine negf_read_hs


!!* Destroys a certain LIBNEGF instance
!!* @param  handler  Handler for the instance to destroy
subroutine negf_destruct_session(handler)
  use libnegfAPICommon                 ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in

  type(NEGFpointers) :: LIB
  integer :: err

  LIB = transfer(handler, LIB)
 
  if (associated(LIB%pNEGF)) deallocate(LIB%pNEGF, stat= err)

  if (err.ne.0) write(*,*) '(negf_destructsession) Deallocation error'

end subroutine negf_destruct_session


!!* Clean the data containers of a given LIBNEGF instance.
!!* @param handler Number for the LIBNEGF instance to destroy.
subroutine negf_destruct_libnegf(handler)
  use libnegfAPICommon  ! if:mod:use
  use libnegf  ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
  call destroy_negf(LIB%pNEGF)

end subroutine negf_destruct_libnegf

!!* Destructs a given LIBNEGF instance.
!!* @param handler Number for the LIBNEGF instance to destroy.
subroutine negf_set_verbosity(handler,verbose_lev)
  use libnegfAPICommon  ! if:mod:use   
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: verbose_lev               ! if:var:in

  type(NEGFpointers) :: LIB
  
  LIB = transfer(handler, LIB) 
   
  LIB%pNEGF%verbose = verbose_lev

end subroutine negf_set_verbosity

!>
!! Solve the Landauer problem: calculate transmission and 
!! density of states according to previously specified parameters
!! @param[in]  handler: handler Number for the LIBNEGF instance
subroutine negf_solve_landauer(handler)
  use libnegfAPICommon  ! if:mod:use 
  use libnegf   ! if:mod:use 
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB) 
  call compute_current(LIB%pNEGF)
end subroutine negf_solve_landauer

!>
!! Get the transmission and density of states real energy
!! array size
!! This information can be used to correctly allocate
!! return arrays.
!! Note: in mpi run every node could have a different 
!! energy interval
!! @param [in]  handler: handler Number for the LIBNEGF instance
!! @param [out] arraySize: size of energy range array
subroutine negf_get_energygrid_size(handler, arraySize)
  use libnegfAPICommon  ! if:mod:use 
  use libnegf   ! if:mod:use 
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: arraySize !if:var:out

  integer :: work
  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB) 
  arraySize = size(LIB%pNEGF%en_grid)

end subroutine negf_get_energygrid_size

!>
!! Get the transmission between a specific lead pair 
!! @param [in]  handler: handler Number for the LIBNEGF instance
!! @param [in] lead_pair: specifies which leads are considered 
!!             for retrieving current (as ordered in ni, nf)
!! @param [in] nsteps: number of energy points, Energies and transmission
!!             should be allocated according to it
!! @param [inout] transmission: array where the transmission is 
!!             copied (note: must be already allocated)
!! TODO: do we really need fixed
subroutine negf_get_transmission(handler, lead_pair, nsteps, energies, transmission)
  use libnegfAPICommon  ! if:mod:use 
  use libnegf   ! if:mod:use 
  use ln_constants !if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: lead_pair ! if:var:in
  integer :: nsteps ! if:var:in
  real(dp), dimension(nsteps) :: energies  ! if:var:inout
  real(dp), dimension(nsteps) :: transmission  ! if:var:inout

  type(NEGFpointers) :: LIB

   LIB = transfer(handler, LIB) 
   energies(:) = real(LIB%pNEGF%en_grid(:)%Ec)
   transmission(:) = LIB%pNEGF%tunn_mat(:,lead_pair)

end subroutine negf_get_transmission

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

!>
!!  Write tunneling and density of states (if any) to file
!!  @param 
!!*        handler:  handler Number for the LIBNEGF instance
!!*        path: string specifying the output file
!!  NOTE: when running parallel code every node will write a 
!!  separate bunch of energy points. You should implement I/O 
!!  OUT of the library and use this routine for testing/debugging
subroutine negf_write_tunneling_and_dos(handler)
  use libnegfAPICommon ! if:mod:use
  use libnegf          ! if:mode:use
  integer :: handler(DAC_handlerSize)  ! if:var:in
  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
  call write_tunneling_and_dos(LIB%pNEGF)
  
end subroutine negf_write_tunneling_and_dos

!!* Compute charge Density given LIBNEGF instance.
!!* @param handler Number for the LIBNEGF instance to destroy.
subroutine negf_density_efa(handler,ndofs,density,particle)
  use libnegfAPICommon  ! if:mod:use 
  use ln_precision      !if:mod:use
  use libnegf           ! if:mod:use 
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: ndofs                     ! if:var:in 
  real(dp) :: density(ndofs)           ! if:var:out
  !! particle: +1 for electrons, -1 for holes
  integer :: particle                  ! if:var:in

  type(NEGFpointers) :: LIB
  
  LIB = transfer(handler, LIB) 

  call compute_density_efa(LIB%pNEGF, density, particle)

end subroutine negf_density_efa

!!* Compute charge Density given LIBNEGF instance.
!!* @param handler Number for the LIBNEGF instance to destroy.
subroutine negf_density_dft(handler,ndofs,density)
  use libnegfAPICommon  ! if:mod:use 
  use ln_precision      !if:mod:use
  use libnegf           ! if:mod:use 
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: ndofs                     ! if:var:in 
  real(dp) :: density(ndofs)           ! if:var:in

  type(NEGFpointers) :: LIB
  
  LIB = transfer(handler, LIB) 

  call compute_density_dft(LIB%pNEGF)

  ! TO DO: OUTPUT 
  !density = LIB%pNEGF%

end subroutine negf_density_dft

!!* Sets iteration in self-consistent loops
!!* @param handler Number for the LIBNEGF instance to destroy.
subroutine negf_set_iteration(handler, iter)
  use libnegfAPICommon  ! if:mod:use 
  use libnegf           ! if:mod:use 
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: iter                      ! if:var:in

  type(NEGFpointers) :: LIB
  
  LIB = transfer(handler, LIB) 

  LIB%pNEGF%iteration = iter

end subroutine negf_set_iteration

!!* Sets iteration in self-consistent loops
!!* @param handler Number for the LIBNEGF instance to destroy.
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

  LIB%pNEGF%minmax = minmax

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
