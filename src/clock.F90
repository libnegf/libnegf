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


module clock

  implicit none
  private
  
 integer, PARAMETER :: NMAXCLKS=5
 integer, public, save :: t1(NMAXCLKS),t2(NMAXCLKS),cr,cm,nclks=0,cpos=0
 logical, public, save :: sus(NMAXCLKS)=.false.

 public :: message_clock, set_clock, write_clock
 
contains

 subroutine message_clock(message)

   character(len=*), intent(in) :: message
   integer :: l_mess
   character(2) :: str_mess, str_dots

   l_mess=len(message)
   write(str_mess,'(I2)') l_mess
   write(str_dots,'(I2)') 54-l_mess

   if (nclks.gt.0 .and. cpos.gt.0) then
       write(6,*)
       sus(nclks) = .true.
   endif

   write(6,FMT='(A'//str_mess//','//str_dots//'("."))',ADVANCE='NO') message 

   flush(6)

   cpos=54

   call set_clock

 end subroutine message_clock

 subroutine set_clock

   if (nclks.lt.NMAXCLKS) then
      nclks = nclks + 1
      call SYSTEM_CLOCK(t1(nclks),cr,cm) 
   endif

 end subroutine set_clock


 subroutine write_clock(message)
   character(len=*), intent(in), optional :: message

   integer :: l_mess
   character(2) :: str_mess, str_dots

   if (nclks.gt.0) then
      call SYSTEM_CLOCK(t2(nclks),cr,cm) 
      if (sus(nclks)) then
         write(6,FMT='(54("."))',ADVANCE='NO')  
         sus(nclks)=.false.
      endif   
      if (present(message)) then
        l_mess=len(message)
        write(str_mess,'(I2)') l_mess
        write(str_dots,'(I2)') 54-l_mess
        write(6,FMT='(A'//str_mess//','//str_dots//'("."))',ADVANCE='NO') message 
      end if
      write(6,*) (t2(nclks)-t1(nclks))*1.0/cr,"sec"
      nclks=nclks-1
      cpos=0
   endif

 end subroutine write_clock


end module clock
