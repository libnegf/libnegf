module structure

implicit none
private

public ::  TStruct_info, create_TStruct, kill_TStruct

! these should be removed and use TStruct_Info
public :: ind, nbl, atmblk, indblk, cstartblk, cblk, cindblk
! --------------------------------------------

type TStruct_Info
   integer, dimension(:), Pointer :: PL_start      !iatm(1)
   integer, dimension(:), Pointer :: PL_end        !iatm(2)
   integer, dimension(:), Pointer :: mat_PL_start  !ind(..)
   integer, dimension(:), Pointer :: mat_PL_end    !ind(..)

   integer, dimension(:), Pointer :: mat_S_start  !    
   integer, dimension(:), Pointer :: mat_C_start  !
   integer, dimension(:), Pointer :: mat_C_end    !

   integer, dimension(:), Pointer :: cblk          !contact int block
   integer :: central_dim
   integer :: num_PLs
   integer :: num_conts
   integer :: active_cont
end type TStruct_Info


 integer, DIMENSION(:), ALLOCATABLE, save :: ind   
 ! iAtomStart(N) starting from 1
 
 integer, SAVE :: nbl
 integer, DIMENSION(:), ALLOCATABLE, SAVE :: atmblk
 !Indice dei vari blocchi nella Hamiltoniana 
 integer, DIMENSION(:), ALLOCATABLE, SAVE :: indblk
 !Indice di partenza dei blocchi relativi ai contatti
 integer, DIMENSION(:), ALLOCATABLE, SAVE :: cstartblk
 !Posizione dei contatti come numero di blocco
 integer, DIMENSION(:), ALLOCATABLE, SAVE :: cblk
 !Indice di posizione dei contatti nella rispettiva matrice
 !(contatti ordinati come negli altri array e indice che parte da 1)
 INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: cindblk


contains
  
  subroutine create_TStruct(ncont,iatm,iatc,str)
    integer, intent(in) :: ncont
    integer, intent(in) :: iatm(2)
    integer, intent(in) :: iatc(3,*)
    type(TStruct_Info), intent(out) :: str

    
    integer :: i
    
    str%num_conts = ncont
    str%active_cont = 1
    !allocate(str%cont_start(ncont))
    !allocate(str%cont_end(ncont))
    allocate(str%mat_C_start(ncont))
    allocate(str%mat_C_end(ncont))
    allocate(str%cblk(ncont))
    
    do i=1,ncont
       !str%cont_start(i) = iatc(3,i)
       !str%cont_end(i) = iatc(2,i)
       !str%mat_C_start(i) = ind(str%cont_start(i))+1
       !str%mat_C_end(i) = ind(str%cont_end(i)+1)
       str%cblk(i) = cblk(i)
    enddo
    
    str%num_PLs = nbl
    allocate(str%PL_start(nbl))
    allocate(str%PL_end(nbl))
    allocate(str%mat_PL_start(nbl))
    allocate(str%mat_PL_end(nbl))
    
    do i=1,nbl
       str%PL_start(i) = atmblk(i)
       str%mat_PL_start(i) = ind(str%PL_start(i))+1
    enddo
    do i=1,nbl-1
       str%PL_end(i) = atmblk(i+1)-1
       str%mat_PL_end(i) = ind(str%PL_end(i)+1)
    enddo
    str%PL_end(nbl) = iatm(2)
    str%mat_PL_end(nbl) = ind(iatm(2)+1)   
    
  end subroutine create_TStruct
  
  ! --------------------------------------------------------
  subroutine kill_TStruct(str)
    type(TStruct_Info) :: str
    
    !deallocate(str%cont_start)
    !deallocate(str%cont_end)
    deallocate(str%mat_C_start)
    deallocate(str%mat_C_end)
    deallocate(str%cblk)
    deallocate(str%PL_start)
    deallocate(str%PL_end)
    deallocate(str%mat_PL_start)
    deallocate(str%mat_PL_end)
  end subroutine kill_TStruct
  
  

end module structure
