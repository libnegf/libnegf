module structure

implicit none
private

public ::  TStruct_info, create_TStruct, print_TStruct, kill_TStruct


type TStruct_Info
   integer, dimension(:), Pointer :: mat_PL_start => NULL() !ind(..)
   integer, dimension(:), Pointer :: mat_PL_end => NULL()    !ind(..)

   integer, dimension(:), Pointer :: mat_B_start => NULL()  ! starting bound    
   integer, dimension(:), Pointer :: mat_C_start => NULL()  ! starting real contact
   integer, dimension(:), Pointer :: mat_C_end => NULL()  ! end real contact

   integer, dimension(:), Pointer :: cblk => NULL()        !contact int block
   integer :: central_dim
   integer :: total_dim
   integer :: num_PLs
   integer :: num_conts
   integer :: active_cont
end type TStruct_Info


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
  subroutine create_TStruct(ncont,nbl, PL_end, cont_end, surf_end, cblk, str)
    integer, intent(in) :: ncont
    integer, intent(in) :: nbl
    integer :: PL_end(*), cont_end(*), surf_end(*), cblk(*)   
    type(TStruct_Info), intent(out) :: str

    
    integer :: i
    
    str%num_conts = ncont
    str%num_PLs = nbl
    str%active_cont = 1

    if(ncont.gt.0) then
       allocate(str%mat_B_start(ncont))
       allocate(str%mat_C_start(ncont))
       allocate(str%mat_C_end(ncont))
       allocate(str%cblk(ncont))
    endif
    
    if(ncont.gt.0) then
       str%mat_B_start(1) =  PL_end(nbl)+1
       str%mat_C_start(1) =  surf_end(1)+1 
       str%mat_C_end(1)   =  cont_end(1)    
       str%cblk(1) = cblk(1)
    endif
       

    do i=2,ncont
       str%mat_B_start(i) = cont_end(i-1) + 1  
       str%mat_C_start(i) = surf_end(i) + 1 
       str%mat_C_end(i)   = cont_end(i)
       str%cblk(i) = cblk(i)
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

end module structure
