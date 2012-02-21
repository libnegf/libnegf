!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Fortran 77 style subroutines for communication with UPTIGHT LIB
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! The order you should invoke the subroutines during initialisation is the
!! following:
!!
!!   negf_initsession
!!   negf_fillbasicparameters (after upt_initsession)
!!   negf_inituptight         (after upt_fillbasicparameters)
!!
!! After initialisation you can call the following methods:
!!
!!   negf_createhamiltonian    
!!   negf_lanczos              (after upt_createhamiltonian)
!!   negf_getzcsrhamiltonian   (after upt_createhamiltonian)
!!
!! In order to destroy a UPT instance call
!!
!!   negf_destructuptight
!!   negf_destructsession      (after upt_destructuptight)
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

subroutine negf_getversion(handler)
  use libnegfAPICommon  ! if:mod:use
  use libnegf  ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:inout
  
  !type(TNEGF), pointer :: pNEGF
  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
 
  call negf_version(LIB%pNEGF)

end subroutine negf_getversion

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

!!* Initialises a new LIBNEGF instance
!!* @param  handler  Contains the handler for the new instance on return
!!* @param  infile 
subroutine negf_fillparameters(handler, infile)
  use libnegfAPICommon  ! if:mod:use
  use libnegf, only : negf_version
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: infile                    ! if:var:in

  type(NEGFpointers) :: LIB

  LIB = transfer(handler, LIB)
 
end subroutine negf_fillparameters

!!* Destroys a certain LIBNEGF instance
!!* @param  handler  Handler for the instance to destroy
subroutine negf_destruct_session(handler)
  use libnegfAPICommon  ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in

  type(NEGFpointers) :: LIB
  integer :: err

  LIB = transfer(handler, LIB)
 
  if (associated(LIB%pNEGF)) deallocate(LIB%pNEGF, stat= err)

  if (err.ne.0) write(*,*) '(negf_destructsession) Deallocation error'

end subroutine negf_destruct_session


!!* Destructs a given LIBNEGF instance.
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


subroutine negf_set_verbosity(handler,verbose_lev)
  use libnegfAPICommon  ! if:mod:use  use negf_param  
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: verbose_lev               ! if:var:in

  type(NEGFpointers) :: LIB
  
  LIB = transfer(handler, LIB) 
   
  LIB%pNEGF%verbose = verbose_lev

end subroutine negf_set_verbosity


subroutine negf_current(handler)
  use libnegfAPICommon  ! if:mod:use  use negf_param 
  use libnegf  ! if:mod:use 
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: verbose_lev               ! if:var:in

  type(NEGFpointers) :: LIB
  
  LIB = transfer(handler, LIB) 

  call extract_compute_current(LIB%pNEGF)

end subroutine negf_current
