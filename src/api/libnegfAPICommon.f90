!!* Contains the type definitions and constants needed by the API routines.
module libnegfAPICommon

  use lib_param
  implicit none
  private

  public :: DAC_handlerSize, TNegf   !, NEGFPointers

  !!* Contains a pointer to a TUPTIn and an OUPT instance
  !type NEGFPointers
  !   type(TNegf), pointer :: plibNEGF
  !end type NEGFPointers
  
  ! Size handler 4 bytes * 4 = 16 bytes
  integer, parameter :: DAC_handlerSize = 4  

  
end module libnegfAPICommon



