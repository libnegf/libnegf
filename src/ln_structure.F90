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


module ln_structure

implicit none
private

public ::  TStruct_info, create_TStruct, print_TStruct, kill_TStruct
public :: neighbor_map

type TStruct_Info
   integer, dimension(:), Pointer :: mat_PL_start => NULL() !ind(..)
   integer, dimension(:), Pointer :: mat_PL_end => NULL()    !ind(..)

   integer, dimension(:), Pointer :: mat_B_start => NULL()  ! starting bound    
   integer, dimension(:), Pointer :: mat_C_start => NULL()  ! starting real contact
   integer, dimension(:), Pointer :: mat_C_end => NULL()  ! end real contact
   integer, dimension(:), Pointer :: cont_dim => NULL()   ! total contact dim

   integer, dimension(:), Pointer :: cblk => NULL()        !contact int block
   integer :: central_dim
   integer :: total_dim
   integer :: num_PLs
   integer :: num_conts
   integer :: active_cont
end type TStruct_Info


type neighbor_map
   integer, DIMENSION(:), ALLOCATABLE :: nn
end type neighbor_map



contains
  
  ! INIT TStructure_Info:
  ! This subroutine requires the definition of
  ! ncont:                   number of contacts
  ! nbl:                     number of PLs
  ! PL_end(nbl)              array of size nbl with the end of each block
  ! cont_end(ncont)          array containing the end of each contact
  ! surf_end(ncont)          array containing the end of each contact surface
  ! cblk(ncont)              array containing the PL-contact position
  !
  subroutine create_TStruct(ncont, nbl, PL_end, cont_end, surf_end, cblk, str)
    integer, intent(in) :: ncont
    integer, intent(in) :: nbl
    integer, dimension(:) :: PL_end, cont_end, surf_end, cblk   
    type(TStruct_Info), intent(inout) :: str

    
    integer :: i

    str%num_conts = ncont
    str%num_PLs = nbl
    str%active_cont = 1

    !  -------------------------------- - - -
    !  | S       |  PL1   |    PL2    |
    !  -------------------------------- - - -
    !  B_start   C_start              C_end 
    !
    !  B_end = C_start - 1

   
    if(ncont.gt.0) then
       allocate(str%mat_B_start(ncont))
       allocate(str%mat_C_start(ncont))
       allocate(str%mat_C_end(ncont))
       allocate(str%cblk(ncont))
       allocate(str%cont_dim(ncont))
    endif

    if(ncont.gt.0) then
       str%mat_B_start(1) =  PL_end(nbl)+1
       str%mat_C_start(1) =  surf_end(1)+1 
       str%mat_C_end(1)   =  cont_end(1)    
       str%cblk(1) = cblk(1)
       str%cont_dim(1) =  str%mat_C_end(1) - str%mat_B_start(1) + 1
    endif
       

    do i=2,ncont
       str%mat_B_start(i) = cont_end(i-1) + 1  
       str%mat_C_start(i) = surf_end(i) + 1 
       str%mat_C_end(i)   = cont_end(i)
       str%cblk(i) = cblk(i)
       str%cont_dim(i) =  str%mat_C_end(i) - str%mat_B_start(i) + 1
    enddo

    allocate(str%mat_PL_start(nbl+1))
    allocate(str%mat_PL_end(nbl))
    
    str%mat_PL_start(1) = 1
    str%mat_PL_end(1) = PL_end(1)

    do i=2,nbl
       str%mat_PL_start(i) = PL_end(i-1)+1
       str%mat_PL_end(i) = PL_end(i)
    enddo
    str%mat_PL_start(nbl+1)= PL_end(nbl)+1

    str%central_dim = PL_end(nbl)
    
    if(ncont.gt.0) then
       str%total_dim = cont_end(ncont)
    else
       str%total_dim = str%central_dim
    endif

  end subroutine create_TStruct
  
  ! --------------------------------------------------------
  subroutine kill_TStruct(str)
    type(TStruct_Info) :: str
    
    if(associated(str%mat_B_start)) deallocate(str%mat_B_start) 
    if(associated(str%mat_C_start)) deallocate(str%mat_C_start)
    if(associated(str%mat_C_end))   deallocate(str%mat_C_end)
    if(associated(str%cblk))  deallocate(str%cblk)

    if(associated(str%mat_PL_start)) deallocate(str%mat_PL_start)
    if(associated(str%mat_PL_end)) deallocate(str%mat_PL_end)
    if(associated(str%cont_dim)) deallocate(str%cont_dim)
  end subroutine kill_TStruct

  ! --------------------------------------------------------  
  subroutine print_Tstruct(str)
    type(Tstruct_Info) :: str
    integer :: i

    write(*,*) 'Hamiltonian Structure:'
    write(*,*) 'num contacts:',str%num_conts 
    write(*,*) 'num layers:',str%num_Pls
    
    do i=1,str%num_PLs
       write(*,*) 'PL:',i, str%mat_PL_start(i), str%mat_PL_end(i)
    enddo

    write(*,*) 'central dim:',str%central_dim

    do i=1,str%num_conts
      write(*,*) 'S: ',i, str%mat_B_start(i), str%mat_C_start(i)-1      
      write(*,*) 'C: ',i, str%mat_C_start(i), str%mat_C_end(i)
      write(*,*) 'cont dim:',str%mat_C_end(i)-str%mat_C_start(i)+1
    enddo  
   
    write(*,*) 'total dim:',str%total_dim

  end subroutine print_Tstruct

end module ln_structure
