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


!!* Contains the type definitions and constants needed by the API routines.
module libnegfAPICommon

  use lib_param
  implicit none
  private

  public :: DAC_handlerSize, NEGFPointers
  public :: TNegf 
  public :: convert_c_string

  !!* Contains a pointer to a TUPTIn and an OUPT instance
  type NEGFPointers
     type(TNegf), pointer :: pNEGF
  end type NEGFPointers
  
  ! Size handler 4 bytes * 4 = 16 bytes
  integer, parameter :: DAC_handlerSize = 4  


  contains

    !!* Convert a c_char array to a fortran string
    subroutine convert_c_string(c_str, f_str)
      use iso_c_binding, only : c_char, c_null_char  ! if:mod:use
      use globals             ! if:mod:use
      implicit none
      character(kind=c_char), intent(in) :: c_str(*) ! if:var:in
      character(len=*), intent(inout)    :: f_str    ! if:var:inout

      integer :: nn, n_char

      !Padding and converting to fortran string
      n_char = 1
      do while(c_str(n_char).ne.c_null_char)
        n_char = n_char + 1
      end do
  
      ! need this, otherwise it will result padded with char(0)
      f_str = " "

      if (len(f_str) .lt. n_char) then
        n_char = len(f_str)
      end if

      do nn=1,n_char-1
        f_str(nn:nn) = c_str(nn)
      end do
      
    end subroutine convert_c_string
  
end module libnegfAPICommon



