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
!! Fortran 77 style subroutines for communication with UPTIGHT LIB
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


!!* Returns the size of the handler array
!!* @param  handlerSize  Contains the size of the handler array on exit.
subroutine negf_gethandlersize(handlerSize)
  use libnegfAPICommon  ! if:mod:use
  implicit none
  integer :: handlerSize  ! if:var:out

  handlerSize = DAC_handlerSize

end subroutine negf_gethandlersize

!!* Initialises a new LIBNEGF instance
!!* @param  handler  Contains the handler for the new instance on return
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

!!* Get library version. Reads svn versions=> broken with git  
!!* @param  handler  Contains the handler for the new instance on return
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
subroutine negf_read_hs(handler)
  use libnegfAPICommon  ! if:mod:use
  use libnegf           ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:inout

  !type(TNEGF), pointer :: pNEGF
  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
 
  call read_HS(LIB%pNEGF)
  
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
  use libnegfAPICommon  ! if:mod:use  use negf_param  
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: verbose_lev               ! if:var:in

  type(NEGFpointers) :: LIB
  
  LIB = transfer(handler, LIB) 
   
  LIB%pNEGF%verbose = verbose_lev

end subroutine negf_set_verbosity

!!* Compute current for a given LIBNEGF instance.
!!* @param handler Number for the LIBNEGF instance to destroy.
subroutine negf_current(handler, current, unitOfH, unitOfJ)
  use libnegfAPICommon  ! if:mod:use  use negf_param 
  use libnegf   ! if:mod:use 
  use ln_constants ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  real(dp) :: current       !if:var:inout
  character(SST) :: unitOfH(1) !if:var:in
  character(SST) :: unitOfJ(1) !if:var:in

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

!!* Compute charge Density given LIBNEGF instance.
!!* @param handler Number for the LIBNEGF instance to destroy.
subroutine negf_density_efa(handler,ndofs,density)
  use libnegfAPICommon  ! if:mod:use  use negf_param 
  use ln_precision      !if:mod:use
  use libnegf           ! if:mod:use 
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: ndofs                     ! if:var:in 
  real(dp) :: density(ndofs)           ! if:var:in

  type(NEGFpointers) :: LIB
  
  LIB = transfer(handler, LIB) 

  call compute_density_efa(LIB%pNEGF, density)

end subroutine negf_density_efa

!!* Compute charge Density given LIBNEGF instance.
!!* @param handler Number for the LIBNEGF instance to destroy.
subroutine negf_density_dft(handler,ndofs,density)
  use libnegfAPICommon  ! if:mod:use  use negf_param 
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
  use libnegfAPICommon  ! if:mod:use  use negf_param 
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
  use libnegfAPICommon    ! if:mod:use  use negf_param 
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
  use libnegfAPICommon    ! if:mod:use  use negf_param 
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
  use libnegfAPICommon    ! if:mod:use  use negf_param 
  use libnegf             ! if:mod:use 
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: outer              ! if:var:in

  type(NEGFpointers) :: LIB
  
  LIB = transfer(handler, LIB) 

  LIB%pNEGF%outer = outer

end subroutine negf_set_outer


subroutine negf_set_kpoint(handler, kpoint)
  use libnegfAPICommon    ! if:mod:use  use negf_param 
  use libnegf             ! if:mod:use 
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: kpoint             ! if:var:in
  
  type(NEGFpointers) :: LIB
  
  LIB = transfer(handler, LIB) 

  LIB%pNEGF%kpoint = kpoint

end subroutine negf_set_kpoint


subroutine negf_set_reference(handler, minmax)
  use libnegfAPICommon    ! if:mod:use  use negf_param 
  use libnegf             ! if:mod:use 
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: minmax             ! if:var:in  
  
  ! minmax = 0 takes minimum; minmax = 1 takes maximum mu
  type(NEGFpointers) :: LIB
  
  LIB = transfer(handler, LIB) 

  LIB%pNEGF%minmax = minmax

end subroutine negf_set_reference

subroutine negf_set_writetunn(handler, flag)
  use libnegfAPICommon  ! if:mod:use  use negf_param  
  use libnegf           ! if:mod:use 
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: flag               ! if:var:in

  type(NEGFpointers) :: LIB
  
  LIB = transfer(handler, LIB) 
 
  LIB%pNEGF%writeTunn = .true.
  if(flag.eq.0) LIB%pNEGF%writeTunn = .false.
 
end subroutine negf_set_writetunn


subroutine negf_set_writeldos(handler, flag)
  use libnegfAPICommon  ! if:mod:use  use negf_param  
  use libnegf           ! if:mod:use 
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: flag               ! if:var:in

  type(NEGFpointers) :: LIB
  
  LIB = transfer(handler, LIB) 
 
  LIB%pNEGF%writeLDOS = .true.
  if(flag.eq.0) LIB%pNEGF%writeLDOS = .false.
 
end subroutine negf_set_writeldos

subroutine negf_write_partition(handler)
  use libnegfAPICommon  ! if:mod:use  use negf_param  
  use libnegf           ! if:mod:use 
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  
  type(NEGFpointers) :: LIB
  
  LIB = transfer(handler, LIB) 

  call negf_partition_info(LIB%pNEGF)

end subroutine negf_write_partition

