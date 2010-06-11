module structure

implicit none
private

public ::  TStruct_info, create_TStruct, kill_TStruct


type TStruct_Info
   integer, dimension(:), Pointer :: mat_PL_start  !ind(..)
   integer, dimension(:), Pointer :: mat_PL_end    !ind(..)

   integer, dimension(:), Pointer :: mat_B_start  ! starting bound    
   integer, dimension(:), Pointer :: mat_C_start  ! starting real contact
   integer, dimension(:), Pointer :: mat_C_end    ! end real contact

   integer, dimension(:), Pointer :: cblk          !contact int block
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
    str$num_PLs = nbl
    str%active_cont = 1

    allocate(str%mat_B_start(ncont))
    allocate(str%mat_C_start(ncont))
    allocate(str%mat_C_end(ncont))
    allocate(str%cblk(ncont))
    
    str%mat_B_start(1) =  PL_end(nbl)+1   
    str%mat_C_start(1) =  surf_end(1)+1 
    str%mat_C_end(1)   =  cont_end(1)    

    do i=2,ncont
       str%mat_B_start(i) = cont_end(i-1) + 1  
       str%mat_C_start(i) = surf_end(i-1) + 1 
       str%mat_C_end(i)   = cont_end(i-1)
       str%cblk(i) = cblk(i)
    enddo

    allocate(str%mat_PL_start(nbl))
    allocate(str%mat_PL_end(nbl))
    
    str%mat_PL_start(1) = 1
    str%mat_PL_end(1) = PL_end(1)

    do i=2,nbl
       str%mat_PL_start(i) = PL_end(i-1)+1
       str%mat_PL_end(i) = PL_end(i)
    enddo

    str%central_dim = PL_end(nbl)
    str%total_dim = cont_end(ncont)
    
  end subroutine create_TStruct
  
  ! --------------------------------------------------------
  subroutine kill_TStruct(str)
    type(TStruct_Info) :: str
    
    deallocate(str%mat_B_start)
    deallocate(str%mat_C_start)
    deallocate(str%mat_C_end)
    deallocate(str%cblk)
    deallocate(str%mat_PL_start)
    deallocate(str%mat_PL_end)
  end subroutine kill_TStruct
  
  

end module structure
