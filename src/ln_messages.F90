module ln_messages
  use mpi_globals
  private

  public :: error_msg

  contains

  subroutine error_msg(str_msg, info)
    character(len=*) :: str_msg
    integer, optional :: info

    if (id0) then
       write(*,*)
       if (present(info)) then
         write(*,*) str_msg, info
       else
         write(*,*) str_msg
       end if  
    end if
    error stop 1
  end subroutine error_msg

end module ln_messages
