module system_calls
 use iso_c_binding
 implicit none
 private

 public :: create_directory
 public :: remove_directory
 public :: remove_file
 
 logical, parameter :: verbose = .false.

 integer, parameter :: INT_EEXIST = 1
 integer, parameter :: INT_ENODIR = 2
 integer, parameter :: INT_EPERM = 3

 interface
   subroutine makedir(dirname, error) bind(C, name="makedir")
     use iso_c_binding, only: c_char, c_int 
     character(kind=C_CHAR), intent(in) :: dirname(*)
     integer(C_INT) :: error    
   end subroutine 
 end interface

 interface
   subroutine removedir(dirname, error) bind(C, name="removedir")
     use iso_c_binding, only: c_char, c_int 
     character(kind=C_CHAR), intent(in) :: dirname(*)
     integer(C_INT) :: error    
   end subroutine 
 end interface

 interface
   subroutine removefile(filename, error) bind(C, name="removefile")
     use iso_c_binding, only: c_char, c_int 
     character(kind=C_CHAR), intent(in) :: filename(*)
     integer(C_INT) :: error    
   end subroutine 
 end interface


 contains

 subroutine create_directory(dirname, error)
   character(len=*), intent(in) :: dirname
   integer, intent(out), optional :: error 

   integer(c_int) :: err
   integer :: ios


   !open(9834, file=trim(dirname)//'/test', iostat=ios)
   !if (ios == 0) then 
   !  if (verbose) print*,'folder '//trim(dirname)//' already exists'
   !  close(9834)
   !  call remove_file(trim(dirname)//'/test')
   !  return
   !end if   

   call makedir(f_to_c_string(dirname), err) 

   if (err == 0) then
     if (verbose) print*,'Folder "'//trim(dirname)//'" created'
   else
     if (present(error)) then
       select case (error) 
       case(INT_EEXIST)
       case(INT_ENODIR)      
         print*,'Folder name error "'//trim(dirname)//'"' 
       case(INT_EPERM)   
         print*,'Permission denied in mkdir "'//trim(dirname)//'"'  
       end select
     endif  
   end if

 end subroutine create_directory 


 subroutine remove_directory(dirname, error)
   character(len=*), intent(in) :: dirname 
   integer, intent(out), optional :: error 
 
   integer(c_int) :: err

   call removedir(f_to_c_string(dirname), err) 
 
   if (err == 0) then
     if (verbose) print*,'Folder "'//trim(dirname)//'" removed'
   else
     print*,'error: could not remove folder. errno=',err
     if (present(error)) then
       error = err
     else    
       stop
     endif  
   end if
 
 end subroutine remove_directory 


 subroutine remove_file(filename, error)
   character(len=*), intent(in) :: filename 
   integer, intent(out), optional :: error 

   integer(c_int) :: err
   logical :: ex

   inquire(file=trim(filename), exist=ex)
   if (.not.ex) return

   call removefile(f_to_c_string(filename), err) 

   if (err == 0) then
     if (verbose) print*,'"'//trim(filename)//'" removed'
   else
     print*,'error: could not remove "'//trim(filename)//'"'
     if (present(error)) then
       error = err
     else    
       stop
     endif  
   end if

 end subroutine remove_file


 !Convert a fortran string into a C character array (NULL terminated)
 pure function f_to_c_string (f_string) result (c_string)
   use, intrinsic :: iso_c_binding, only: c_char, c_null_char
   implicit none   

   character(len=*), intent(in) :: f_string
   character(len=1, kind=c_char), dimension(len_trim(f_string)+1) :: c_string

   integer :: n, i
   n= len_trim(f_string)
   do i = 1, n
     c_string(i) = f_string(i:i)
   end do
   c_string(n+1) = c_null_char

 end function f_to_c_string


end module system_calls
