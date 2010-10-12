Module mat_def
  use precision  
  use allocation
  implicit none
  private

public :: z_CSR,z_CSC,z_MSR,z_COO,z_EXT_COO,z_DNS
public :: r_CSR,r_CSC,r_MSR,r_COO,r_DNS, z_vec, z_RGM

public :: create, init, recreate, destroy, create_id
public :: print_mat, read_mat, writemem

interface create
   module procedure zcreate_CSR
   module procedure rcreate_CSR
   module procedure zcreate_CSC
   module procedure rcreate_CSC
   module procedure zcreate_MSR
   module procedure rcreate_MSR
   module procedure zcreate_COO
   module procedure rcreate_COO
   module procedure zcreate_EXT_COO
   module procedure zcreate_DNS
   module procedure rcreate_DNS
end interface

interface create_id
   module procedure zcreate_id_CSR
end interface

interface recreate
   module procedure zrecreate_CSR
   module procedure rrecreate_CSR
end interface

interface init
   module procedure zinit_CSR
end interface



interface destroy
   module procedure zdestroy_CSR
   module procedure rdestroy_CSR
   module procedure zdestroy_CSC
   module procedure rdestroy_CSC
   module procedure zdestroy_MSR
   module procedure rdestroy_MSR
   module procedure zdestroy_COO
   module procedure rdestroy_COO
   module procedure zdestroy_EXT_COO
   module procedure zdestroy_DNS
   module procedure rdestroy_DNS
end interface

interface print_mat
   module procedure zPrint_CSR
   module procedure rPrint_CSR
   module procedure zPrint_CSC
   module procedure rPrint_CSC
   module procedure zPrint_MSR
   module procedure rPrint_MSR
   module procedure zPrint_COO
   module procedure rPrint_COO
end interface

interface read_mat
   module procedure zRead_CSR  
end interface

interface writemem
   module procedure zwriteMem_CSR
end interface

Type z_CSR
  integer :: nnz
  integer :: nrow 
  integer :: ncol
  complex(kind=dp), DIMENSION(:), ALLOCATABLE :: nzval  
  integer, DIMENSION(:), ALLOCATABLE :: colind
  integer, DIMENSION(:), ALLOCATABLE :: rowpnt
end Type z_CSR

!CSC Complex Structure definition (Compressed Sparse Column format)

Type z_CSC
  integer :: nnz
  integer :: nrow 
  integer :: ncol
  complex(kind=dp), DIMENSION(:), ALLOCATABLE :: nzval 
  integer, DIMENSION(:), ALLOCATABLE :: rowind 
  integer, DIMENSION(:), ALLOCATABLE :: colpnt 
end Type z_CSC

!MSR Complex Structure definition (Modified Sparse Row Format)

Type z_MSR
  integer :: nnz
  integer :: nrow
  integer :: ncol
  complex(kind=dp), DIMENSION(:), ALLOCATABLE :: nzval 
  integer, DIMENSION(:), ALLOCATABLE :: index 
end Type z_MSR

Type z_COO
  integer :: nnz
  integer :: nrow
  integer :: ncol
  integer, DIMENSION(:), ALLOCATABLE :: index_i 
  integer, DIMENSION(:), ALLOCATABLE :: index_j 
  complex(kind=dp), DIMENSION(:), ALLOCATABLE :: nzval  
end Type z_COO

Type z_EXT_COO
  integer :: nnz
  integer :: nrow
  integer :: ncol
  integer, DIMENSION(:), ALLOCATABLE :: index_i 
  integer, DIMENSION(:), ALLOCATABLE :: index_j 
  complex(kind=dp), DIMENSION(:), ALLOCATABLE :: nzval
  logical, DIMENSION(:), ALLOCATABLE :: first
end Type z_EXT_COO 

Type z_DNS
  integer :: nrow
  integer :: ncol
  complex(kind=dp), DIMENSION(:,:), ALLOCATABLE :: val  
end Type z_DNS 

Type z_vec
   integer :: len
   complex(dp), dimension(:), allocatable :: val 
   integer, dimension(:), allocatable :: ind 
end Type z_vec

Type z_RGM
   integer :: nrow
   type(z_vec), dimension(:), allocatable :: row
end Type z_RGM


 
!CSR Complex Structure definition (Compressed Sparse Row Format)
Type r_CSR
  integer :: nnz
  integer :: nrow
  integer :: ncol
  real(kind=dp), DIMENSION(:), ALLOCATABLE :: nzval  
  integer, DIMENSION(:), ALLOCATABLE :: colind 
  integer, DIMENSION(:), ALLOCATABLE :: rowpnt 
end Type r_CSR

!CSC Complex Structure definition (Compressed Sparse Column format)

Type r_CSC
  integer :: nnz
  integer :: nrow 
  integer :: ncol
  real(kind=dp), DIMENSION(:), ALLOCATABLE :: nzval 
  integer, DIMENSION(:), ALLOCATABLE :: rowind 
  integer, DIMENSION(:), ALLOCATABLE :: colpnt 
end Type r_CSC

!MSR Real Structure definition (Modified Sparse Row Format)

Type r_MSR
  integer :: nnz
  integer :: ncol
  integer :: nrow
  real(kind=dp), DIMENSION(:), ALLOCATABLE :: nzval 
  integer, DIMENSION(:), ALLOCATABLE :: index 
end Type r_MSR

Type r_COO
  integer :: nnz
  integer :: nrow
  integer :: ncol
  integer, DIMENSION(:), ALLOCATABLE :: index_i 
  integer, DIMENSION(:), ALLOCATABLE :: index_j 
  real(kind=dp), DIMENSION(:), ALLOCATABLE :: nzval  
end Type r_COO 

Type r_DNS
  integer :: nrow
  integer :: ncol
  real(kind=dp), DIMENSION(:,:), ALLOCATABLE :: val  
end Type r_DNS 
! *******************************************************************
contains

!Utilities per z_CSR format

subroutine zcreate_CSR(mat,nrow,ncol,nnz)
  type(z_CSR) :: mat
  integer :: nrow, ncol, nnz, ierr

  mat%nnz=nnz
  mat%nrow=nrow
  mat%ncol=ncol

  if(nnz.ne.0) then
     call log_allocate(mat%nzval,nnz)
     call log_allocate(mat%colind,nnz)
  endif

  call log_allocate(mat%rowpnt,nrow+1)
    
end subroutine zcreate_CSR
! ------------------------------------------------------------------

subroutine zinit_CSR(mat)
  type(z_CSR) :: mat

  integer :: k

  if(mat%nnz.ne.mat%nrow) STOP 'cannot initialize matrix (nnz != nrow)' 

  do k=1,Mat%nrow
     Mat%nzval(k)=0.d0
     Mat%colind(k)=k
     Mat%rowpnt(k)=k
  enddo
  Mat%rowpnt(mat%nrow+1)=mat%nrow+1
  
end subroutine zinit_CSR
! ------------------------------------------------------------------
subroutine zcreate_id_CSR(mat,nrow)
  type(z_CSR) :: mat
  integer nrow, i

  call zcreate_CSR(mat,nrow,nrow,nrow)
  do i=1,nrow
     mat%nzval(i)=1.d0
     mat%rowpnt(i)=i
     mat%colind(i)=i
  enddo
  mat%rowpnt(nrow+1)=nrow+1

end subroutine zcreate_id_CSR


! ------------------------------------------------------------------
subroutine zrecreate_CSR(mat)
  type(z_CSR) :: mat
  integer :: nrow, ncol, nnz, ierr
  integer :: k

  nnz=mat%nnz
  nrow=mat%nrow
  ncol=mat%ncol

  call zdestroy_CSR(mat)
  call zcreate_CSR(mat,ncol,ncol,ncol)
    
   do k=1,ncol
     mat%nzval(k)=0.d0
     mat%colind(k)=k
     mat%rowpnt(k)=k
  enddo
  mat%rowpnt(ncol+1)=ncol+1

end subroutine zrecreate_CSR
! ------------------------------------------------------------------
subroutine zdestroy_CSR(mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat8)

  type(z_CSR) :: mat1
  type(z_CSR), optional :: mat2,mat3,mat4,mat5,mat6,mat7,mat8
  
  mat1%nnz=0
  mat1%nrow=0
  mat1%ncol=0

  if (allocated(mat1%nzval)) then
     call log_deallocate(mat1%nzval)
     call log_deallocate(mat1%colind)
  endif
  call log_deallocate(mat1%rowpnt)

  if (present(mat2)) then
     mat2%nnz=0
     mat2%nrow=0
     mat2%ncol=0
     if (allocated(mat2%nzval)) then
        call log_deallocate(mat2%nzval)
        call log_deallocate(mat2%colind)
     endif
     call log_deallocate(mat2%rowpnt)
  else
     return
  endif

  if (present(mat3)) then
     mat3%nnz=0
     mat3%nrow=0
     mat3%ncol=0
     if (allocated(mat3%nzval)) then
        call log_deallocate(mat3%nzval)
        call log_deallocate(mat3%colind)
     endif
     call log_deallocate(mat3%rowpnt)
  else
     return
  end if

  if (present(mat4)) then
     mat4%nnz=0
     mat4%nrow=0
     mat4%ncol=0
     if (allocated(mat4%nzval)) then
        call log_deallocate(mat4%nzval)
        call log_deallocate(mat4%colind)
     endif
     call log_deallocate(mat4%rowpnt)
  else
     return
  endif

  if (present(mat5)) then
     mat5%nnz=0
     mat5%nrow=0
     mat5%ncol=0
     if (allocated(mat5%nzval)) then
        call log_deallocate(mat5%nzval)
        call log_deallocate(mat5%colind)
     endif
     call log_deallocate(mat5%rowpnt)
  else
     return
  endif

  if (present(mat6)) then
     mat6%nnz=0
     mat6%nrow=0
     mat6%ncol=0
     if (allocated(mat6%nzval)) then
        call log_deallocate(mat6%nzval)
        call log_deallocate(mat6%colind)
     endif
     call log_deallocate(mat6%rowpnt)
  else
     return
  endif

  if (present(mat7)) then
     mat7%nnz=0
     mat7%nrow=0
     mat7%ncol=0
     if (allocated(mat7%nzval)) then
        call log_deallocate(mat7%nzval)
        call log_deallocate(mat7%colind)
     endif
     call log_deallocate(mat7%rowpnt)
  else
     return
  endif

  if (present(mat8)) then
     mat8%nnz=0
     mat8%nrow=0
     mat8%ncol=0
     if (allocated(mat8%nzval)) then
        call log_deallocate(mat8%nzval)
        call log_deallocate(mat8%colind)
     endif
     call log_deallocate(mat8%rowpnt)
  else
     return
  endif


end subroutine zdestroy_CSR
! ------------------------------------------------------------------
subroutine zPrint_CSR(id,sp,fmt)

  type(z_CSR) :: sp
  integer :: id,i
  logical :: fmt
  
  if(fmt) then
     write(id,*) 'Nrow: ',sp%nrow
     write(id,*) 'Ncol: ',sp%ncol 
     write(id,*) 'Nnz: ',sp%nnz
     
     write(id,*) 'Nzval array'
     do i=1,sp%nnz
        write(id,*) sp%nzval(i)
     enddo
     write(id,*) 'Colind array'
     do i=1,sp%nnz
        write(id,*) sp%colind(i)
     enddo
     write(id,*) 'Rowpnt array'
     if(sp%nrow.gt.0) then
     do i=1,sp%nrow+1   
        write(id,*) sp%rowpnt(i)
     enddo
     endif
  else
     write(id) sp%nrow
     write(id) sp%ncol 
     write(id) sp%nnz
     
     !write(id) 'Nzval array'
     do i=1,sp%nnz
        write(id) sp%nzval(i)
     enddo
     !write(id) 'Colind array'
     do i=1,sp%nnz
        write(id) sp%colind(i)
     enddo
     !write(id) 'Rowpnt array'
     if(sp%nrow.gt.0) then
     do i=1,sp%nrow+1
        write(id) sp%rowpnt(i)
     enddo
     endif
  endif

end subroutine zPrint_CSR
! ------------------------------------------------------------------
subroutine zRead_CSR(id,sp,fmt)

  type(z_CSR) :: sp
  integer :: id,i
  character(20) :: tmpst
  logical :: fmt

  if(fmt) then
     read(id,*) tmpst,sp%nrow
     read(id,*) tmpst,sp%ncol 
     read(id,*) tmpst,sp%nnz
     
     read(id,*) tmpst
     do i=1,sp%nnz
        read(id,*) sp%nzval(i)
     enddo
     read(id,*) tmpst
     do i=1,sp%nnz
        read(id,*) sp%colind(i)
     enddo
     read(id,*) tmpst
     
     if(sp%nrow.gt.0) then
     do i=1,sp%nrow+1
        read(id,*) sp%rowpnt(i)
     enddo
     endif
  else
     read(id) sp%nrow
     read(id) sp%ncol 
     read(id) sp%nnz
     
     do i=1,sp%nnz
        read(id) sp%nzval(i)
     enddo
     
     do i=1,sp%nnz
        read(id) sp%colind(i)
     enddo
     
     if(sp%nrow.gt.0) then
     do i=1,sp%nrow+1
        read(id) sp%rowpnt(i)
     enddo
     endif
  endif

end subroutine zRead_CSR
! ------------------------------------------------------------------
!Utilities per z_CSC format

subroutine zcreate_CSC(mat,nrow,ncol,nnz)

  type(z_CSC) :: mat
  integer :: nrow, ncol, nnz, ierr

  if(nrow.eq.0.or.ncol.eq.0) STOP 'ERROR: nrow or ncol = 0'

  mat%nnz=nnz
  mat%nrow=nrow
  mat%ncol=ncol
  if(nnz.ne.0) then
     call log_allocate(mat%nzval,nnz)
     call log_allocate(mat%rowind,nnz)
  endif

  call log_allocate(mat%colpnt,ncol+1)

end subroutine zcreate_CSC
! ------------------------------------------------------------------

subroutine zdestroy_CSC(mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat8)

  type(z_CSC) :: mat1
  type(z_CSC), optional :: mat2,mat3,mat4,mat5,mat6,mat7,mat8
  
  mat1%nnz=0
  mat1%nrow=0
  mat1%ncol=0

  if (allocated(mat1%nzval)) then
     call log_deallocate(mat1%nzval)
     call log_deallocate(mat1%rowind)
  endif
  call log_deallocate(mat1%colpnt)

  if (present(mat2)) then
     mat2%nnz=0
     mat2%nrow=0
     mat2%ncol=0
     if (allocated(mat2%nzval)) then
        call log_deallocate(mat2%nzval)
        call log_deallocate(mat2%rowind)
     endif
     call log_deallocate(mat2%colpnt)
  else
     return
  endif

  if (present(mat3)) then
     mat3%nnz=0
     mat3%nrow=0
     mat3%ncol=0
     if (allocated(mat3%nzval)) then
        call log_deallocate(mat3%nzval)
        call log_deallocate(mat3%rowind)
     endif
     call log_deallocate(mat3%colpnt)
  else
     return
  endif

  if (present(mat4)) then
     mat4%nnz=0
     mat4%nrow=0
     mat4%ncol=0
     if (allocated(mat4%nzval)) then
        call log_deallocate(mat4%nzval)
        call log_deallocate(mat4%rowind)
     endif
     call log_deallocate(mat4%colpnt)
  else
    return
  endif

  if (present(mat5)) then
     mat5%nnz=0
     mat5%nrow=0
     mat5%ncol=0
     if (allocated(mat5%nzval)) then
        call log_deallocate(mat5%nzval)
        call log_deallocate(mat5%rowind)
     endif
     call log_deallocate(mat5%colpnt)
  else
     return
  endif

  if (present(mat6)) then
     mat6%nnz=0
     mat6%nrow=0
     mat6%ncol=0
     if (allocated(mat6%nzval)) then
        call log_deallocate(mat6%nzval)
        call log_deallocate(mat6%rowind)
     endif
     call log_deallocate(mat6%colpnt)
  else
     return
  endif

  if (present(mat7)) then
     mat7%nnz=0
     mat7%nrow=0
     mat7%ncol=0
     if (allocated(mat7%nzval)) then
        call log_deallocate(mat7%nzval)
        call log_deallocate(mat7%rowind)
     endif
     call log_deallocate(mat7%colpnt)
  else
     return
  endif

  if (present(mat8)) then
     mat8%nnz=0
     mat8%nrow=0
     mat8%ncol=0
     if (allocated(mat8%nzval)) then
        call log_deallocate(mat8%nzval)
        call log_deallocate(mat8%rowind)
     endif
     call log_deallocate(mat8%colpnt)
  else
     return
  endif


end subroutine zdestroy_CSC
! ------------------------------------------------------------------

subroutine zPrint_CSC(sp)  !Identical to Print_CSR

  type(z_CSC) :: sp
  integer :: i
  
  do i=1,sp%nnz
    write(*,*) sp%nzval(i)
  enddo  

end subroutine zprint_CSC
! ------------------------------------------------------------------
!Utilities per MSR format

subroutine zcreate_MSR(mat,nrow,ncol,nnz)

  type(z_MSR) :: mat
  integer :: nrow, ncol, nnz, ierr

  if(nrow.eq.0.or.ncol.eq.0) STOP 'ERROR: nrow or ncol = 0'

  mat%nnz=nnz
  mat%nrow=nrow
  mat%ncol=ncol
  if(nnz.ne.0) then
     call log_allocate(mat%nzval,nnz+1)
     call log_allocate(mat%index,nnz+1)
  endif
end subroutine zcreate_MSR
! ------------------------------------------------------------------
subroutine zdestroy_MSR(mat)

  type(z_MSR) :: mat
  
  mat%nnz=0
  mat%nrow=0
  mat%ncol=0

  if (allocated(mat%nzval)) then 
     call log_deallocate(mat%nzval)
     call log_deallocate(mat%index)
  endif

end subroutine zdestroy_MSR
! ------------------------------------------------------------------

subroutine zprint_MSR(sp)

  type(z_MSR) :: sp
  integer :: i
  
  do i=1,sp%nnz
    write(*,*) sp%nzval(i)
  enddo  

end subroutine zprint_MSR
! ------------------------------------------------------------------


subroutine zcreate_COO(mat,nrow,ncol,nnz)

  type(z_COO) :: mat
  integer :: nrow, ncol, nnz, ierr

  if(nrow.eq.0.or.ncol.eq.0) STOP 'ERROR: nrow or ncol = 0'

  mat%nnz=nnz
  mat%nrow=nrow
  mat%ncol=ncol

  if(nnz.ne.0) then
     call log_allocate(mat%nzval,nnz)
     call log_allocate(mat%index_i,nnz)
     call log_allocate(mat%index_j,nnz)
  endif

end subroutine zcreate_COO
! ------------------------------------------------------------------
subroutine zdestroy_COO(mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat8)

  type(z_COO) :: mat1
  type(z_COO), optional :: mat2,mat3,mat4,mat5,mat6,mat7,mat8
  
  mat1%nnz=0
  mat1%nrow=0
  mat1%ncol=0

  if (allocated(mat1%nzval)) then
     call log_deallocate(mat1%nzval)
     call log_deallocate(mat1%index_i)
     call log_deallocate(mat1%index_j)
  end if

  if (present(mat2)) then
     mat2%nnz=0
     mat2%nrow=0
     mat2%ncol=0
     if (allocated(mat2%nzval)) then
        call log_deallocate(mat2%nzval)
        call log_deallocate(mat2%index_i)
        call log_deallocate(mat2%index_j)
     end if
  else
     return
  endif

  if (present(mat3)) then
     mat3%nnz=0
     mat3%nrow=0
     mat3%ncol=0
     if (allocated(mat3%nzval)) then
        call log_deallocate(mat3%nzval)
        call log_deallocate(mat3%index_i)
        call log_deallocate(mat3%index_j)
     end if
  else
     return
  endif

  if (present(mat4)) then
     mat4%nnz=0
     mat4%nrow=0
     mat4%ncol=0
     if (allocated(mat4%nzval)) then
        call log_deallocate(mat4%nzval)
        call log_deallocate(mat4%index_i)
        call log_deallocate(mat4%index_j)
     end if
  else
     return
  endif

  if (present(mat5)) then
     mat5%nnz=0
     mat5%nrow=0
     mat5%ncol=0
     if (allocated(mat5%nzval)) then
        call log_deallocate(mat5%nzval)
        call log_deallocate(mat5%index_i)
        call log_deallocate(mat5%index_j)
     end if
  else
     return
  endif

  if (present(mat6)) then
     mat6%nnz=0
     mat6%nrow=0
     mat6%ncol=0
     if (allocated(mat6%nzval)) then
        call log_deallocate(mat6%nzval)
        call log_deallocate(mat6%index_i)
        call log_deallocate(mat6%index_j)
     end if
  else
     return
  endif

  if (present(mat7)) then
     mat7%nnz=0
     mat7%nrow=0
     mat7%ncol=0
     if (allocated(mat7%nzval)) then
        call log_deallocate(mat7%nzval)
        call log_deallocate(mat7%index_i)
        call log_deallocate(mat7%index_j)
     endif
  else
     return
  end if

  if (present(mat8)) then
     mat8%nnz=0
     mat8%nrow=0
     mat8%ncol=0
     if (allocated(mat8%nzval)) then
        call log_deallocate(mat8%nzval)
        call log_deallocate(mat8%index_i)
        call log_deallocate(mat8%index_j)
     endif
  else
     return
  end if


end subroutine zdestroy_COO
! ------------------------------------------------------------------

subroutine zprint_COO(id,sp,fmt)

  type(z_COO) :: sp
  integer :: id,i
  logical :: fmt
  
  if(fmt) then
     write(id,*) sp%nrow,'Nrow'
     write(id,*) sp%ncol,'Ncol'
     write(id,*) sp%nnz,'Nnz'
  
     do i=1,sp%nnz  
        write(id,*) sp%index_i(i),sp%index_j(i),sp%nzval(i)
     enddo

  else

     write(id) sp%nrow
     write(id) sp%ncol 
     write(id) sp%nnz
  
     do i=1,sp%nnz  
        write(id) sp%index_i(i),sp%index_j(i),sp%nzval(i)
     enddo

  end if

end subroutine zprint_COO
! ------------------------------------------------------------------
! ------------------------------------------------------------------


subroutine zcreate_EXT_COO(mat,nrow,ncol,nnz)

  type(z_EXT_COO) :: mat
  integer :: nrow, ncol, nnz, ierr

  if(nrow.eq.0.or.ncol.eq.0) STOP 'ERROR: nrow or ncol = 0'

  mat%nnz=nnz
  mat%nrow=nrow
  mat%ncol=ncol

  if(nnz.ne.0) then
     call log_allocate(mat%nzval,nnz)
     call log_allocate(mat%index_i,nnz)
     call log_allocate(mat%index_j,nnz)
     call log_allocate(mat%first,nnz)
  endif

  mat%first(:) = .false.
  
end subroutine zcreate_EXT_COO
! ------------------------------------------------------------------

subroutine zdestroy_EXT_COO(mat)

  type(z_EXT_COO) :: mat
  
  mat%nnz=0
  mat%nrow=0
  mat%ncol=0
  if (allocated(mat%nzval)) then 
     call log_deallocate(mat%nzval)
     call log_deallocate(mat%index_i)
     call log_deallocate(mat%index_j)
     call log_deallocate(mat%first)
  endif

end subroutine zdestroy_EXT_COO
! ------------------------------------------------------------------

subroutine zcreate_DNS(mat,nrow,ncol)

  type(z_DNS) :: mat
  integer :: nrow, ncol, ierr

  if(nrow.eq.0.or.ncol.eq.0) STOP 'ERROR: nrow or ncol = 0'

  mat%ncol=ncol
  mat%nrow=nrow
  call log_allocate(mat%val,nrow,ncol)

end subroutine zcreate_DNS
! ------------------------------------------------------------------

subroutine zdestroy_DNS(mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat8)

  type(z_DNS) :: mat1
  type(z_DNS), optional :: mat2,mat3,mat4,mat5,mat6,mat7,mat8
  
  mat1%nrow=0
  mat1%ncol=0

  if (allocated(mat1%val)) then
     call log_deallocate(mat1%val)
  end if

  if (present(mat2)) then
     mat2%nrow=0
     mat2%ncol=0
     if (allocated(mat2%val)) then
        call log_deallocate(mat2%val)
     end if
  else
     return
  endif

  if (present(mat3)) then
     mat3%nrow=0
     mat3%ncol=0
     if (allocated(mat3%val)) then
        call log_deallocate(mat3%val)
     end if
  else
     return
  endif

  if (present(mat4)) then
     mat4%nrow=0
     mat4%ncol=0
     if (allocated(mat4%val)) then
        call log_deallocate(mat4%val)
     endif
  else
     return
  endif

  if (present(mat5)) then
     mat5%nrow=0
     mat5%ncol=0
     if (allocated(mat5%val)) then
        call log_deallocate(mat5%val)
     end if
  else
     return
  endif

  if (present(mat6)) then
     mat6%nrow=0
     mat6%ncol=0
     if (allocated(mat6%val)) then
        call log_deallocate(mat6%val)
     end if
  else
     return
  endif

  if (present(mat7)) then
     mat7%nrow=0
     mat7%ncol=0
     if (allocated(mat7%val)) then
        call log_deallocate(mat7%val)
     endif
  else
     return
  end if

  if (present(mat8)) then
     mat8%nrow=0
     mat8%ncol=0
     if (allocated(mat8%val)) then
        call log_deallocate(mat8%val)
     endif
  else
     return
  end if

end subroutine zdestroy_DNS


! ------------------------------------------------------------------
! ==================================================================

!Utilities per r_CSR format

subroutine rcreate_CSR(mat,nrow,ncol,nnz)

  type(r_CSR) :: mat
  integer :: nrow, ncol, nnz, ierr

  if(nrow.eq.0.or.ncol.eq.0) STOP 'ERROR: nrow or ncol = 0'

  mat%nnz=nnz
  mat%nrow=nrow
  mat%ncol=ncol
  if(nnz.ne.0) then
     call log_allocate(mat%nzval,nnz)
     call log_allocate(mat%colind,nnz)
  endif

  call log_allocate(mat%rowpnt,nrow+1)

end subroutine rcreate_CSR
! ------------------------------------------------------------------
subroutine rrecreate_CSR(mat)
  type(r_CSR) :: mat
  integer :: nrow, ncol, nnz, ierr
  integer :: k

  nnz=mat%nnz
  nrow=mat%nrow
  ncol=mat%ncol

  call rdestroy_CSR(mat)
  call rcreate_CSR(mat,ncol,ncol,ncol)
    
   do k=1,ncol
     mat%nzval(k)=0.d0
     mat%colind(k)=k
     mat%rowpnt(k)=k
  enddo
  mat%rowpnt(ncol+1)=ncol+1

end subroutine rrecreate_CSR
! ------------------------------------------------------------------

subroutine rdestroy_CSR(mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat8)

  type(r_CSR) :: mat1
  type(r_CSR), optional :: mat2,mat3,mat4,mat5,mat6,mat7,mat8
  
  mat1%nnz=0
  mat1%nrow=0
  mat1%ncol=0

  if (allocated(mat1%nzval)) then
     call log_deallocate(mat1%nzval)
     call log_deallocate(mat1%colind)
  endif
  call log_deallocate(mat1%rowpnt)

  if (present(mat2)) then
     mat2%nnz=0
     mat2%nrow=0
     mat2%ncol=0
     if (allocated(mat2%nzval)) then
        call log_deallocate(mat2%nzval)
        call log_deallocate(mat2%colind)
     endif
     call log_deallocate(mat2%rowpnt)
  else
     return
  endif

  if (present(mat3)) then
     mat3%nnz=0
     mat3%nrow=0
     mat3%ncol=0
     if (allocated(mat3%nzval)) then
        call log_deallocate(mat3%nzval)
        call log_deallocate(mat3%colind)
     endif
     call log_deallocate(mat3%rowpnt)
  else
     return
  end if

  if (present(mat4)) then
     mat4%nnz=0
     mat4%nrow=0
     mat4%ncol=0
     if (allocated(mat4%nzval)) then
        call log_deallocate(mat4%nzval)
        call log_deallocate(mat4%colind)
     endif
     call log_deallocate(mat4%rowpnt)
  else
     return
  endif

  if (present(mat5)) then
     mat5%nnz=0
     mat5%nrow=0
     mat5%ncol=0
     if (allocated(mat5%nzval)) then
        call log_deallocate(mat5%nzval)
        call log_deallocate(mat5%colind)
     endif
     call log_deallocate(mat5%rowpnt)
  else
     return
  endif

  if (present(mat6)) then
     mat6%nnz=0
     mat6%nrow=0
     mat6%ncol=0
     if (allocated(mat6%nzval)) then
        call log_deallocate(mat6%nzval)
        call log_deallocate(mat6%colind)
     endif
     call log_deallocate(mat6%rowpnt)
  else
     return
  endif

  if (present(mat7)) then
     mat7%nnz=0
     mat7%nrow=0
     mat7%ncol=0
     if (allocated(mat7%nzval)) then
        call log_deallocate(mat7%nzval)
        call log_deallocate(mat7%colind)
     endif
     call log_deallocate(mat7%rowpnt)
  else
     return
  endif

  if (present(mat8)) then
     mat8%nnz=0
     mat8%nrow=0
     mat8%ncol=0
     if (allocated(mat8%nzval)) then
        call log_deallocate(mat8%nzval)
        call log_deallocate(mat8%colind)
     endif
     call log_deallocate(mat7%rowpnt)
  else
     return
  endif


end subroutine rdestroy_CSR

! ------------------------------------------------------------------

subroutine rPrint_CSR(id,sp)

  type(r_CSR) :: sp
  integer :: i,id
  
  write(id,*) 'Nrow is ',sp%nrow
  write(id,*) 'Nnz is ',sp%nnz

  write(id,*) 'Nzval array'
  do i=1,sp%nnz
    write(id,*) sp%nzval(i)
  enddo 
  write(id,*) 'Colind array'
  do i=1,sp%nnz
    write(id,*) sp%colind(i)
  enddo  
  write(id,*) 'Rowpnt array'
  do i=1,sp%nrow+1
    write(id,*) sp%rowpnt(i)
  enddo 

end subroutine rPrint_CSR
! ------------------------------------------------------------------
!Utilities per z_CSC format

subroutine rcreate_CSC(mat,nrow,ncol,nnz)

  type(r_CSC) :: mat
  integer :: nrow, ncol, nnz, ierr
  
  if(nrow.eq.0.or.ncol.eq.0) STOP 'ERROR: nrow or ncol = 0'
   
  mat%nnz=nnz
  mat%nrow=nrow
  mat%ncol=ncol
  if (nnz.ne.0) then
     call log_allocate(mat%nzval,nnz)
     call log_allocate(mat%rowind,nnz)
  endif
  call log_allocate(mat%colpnt,ncol+1)

end subroutine rcreate_CSC
! ------------------------------------------------------------------
subroutine rdestroy_CSC(mat)

  type(r_CSC) :: mat
  
  mat%nnz=0
  mat%nrow=0
  mat%ncol=0
  if (allocated(mat%nzval)) then 
     call log_deallocate(mat%nzval)
     call log_deallocate(mat%rowind)
  endif  
  call log_deallocate(mat%colpnt)

end subroutine rdestroy_CSC
! ------------------------------------------------------------------

subroutine rPrint_CSC(sp)  !Identical to Print_CSR

  type(r_CSC) :: sp
  integer :: i
  
  do i=1,sp%nnz
    write(*,*) sp%nzval(i)
  enddo  

end subroutine rprint_CSC
! ------------------------------------------------------------------
!Utilities per MSR format

subroutine rcreate_MSR(mat,nrow,ncol,nnz)

  type(r_MSR) :: mat
  integer :: nrow, ncol, nnz, ierr

  if(nrow.eq.0.or.ncol.eq.0) STOP 'ERROR: nrow or ncol = 0'

  mat%nnz=nnz
  mat%nrow=nrow
  mat%ncol=ncol
  if (nnz.ne.0) then
     call log_allocate(mat%nzval,nnz+1)
     call log_allocate(mat%index,nnz+1)
  endif
end subroutine rcreate_MSR
! ------------------------------------------------------------------
subroutine rdestroy_MSR(mat)

  type(r_MSR) :: mat
  
  mat%nnz=0
  mat%nrow=0
  if (allocated(mat%nzval)) then 
     call log_deallocate(mat%nzval)
     call log_deallocate(mat%index)
  endif
end subroutine rdestroy_MSR
! ------------------------------------------------------------------

subroutine rprint_MSR(sp)

  type(r_MSR) :: sp
  integer :: i
  
  do i=1,sp%nnz
    write(*,*) sp%nzval(i)
  enddo  

end subroutine rprint_MSR
! ------------------------------------------------------------------


subroutine rcreate_COO(mat,nrow,ncol,nnz)

  type(r_COO) :: mat
  integer :: nrow, ncol, nnz, ierr

  if(nrow.eq.0.or.ncol.eq.0) STOP 'ERROR: nrow or ncol = 0'

  mat%nnz=nnz
  mat%nrow=nrow
  mat%ncol=ncol
  if (nnz.ne.0) then
     call log_allocate(mat%nzval,nnz)
     call log_allocate(mat%index_i,nnz)
     call log_allocate(mat%index_j,nnz)
  endif
end subroutine rcreate_COO
! ------------------------------------------------------------------

subroutine rdestroy_COO(mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat8)

  type(r_COO) :: mat1
  type(r_COO), optional :: mat2,mat3,mat4,mat5,mat6,mat7,mat8
  
  mat1%nnz=0
  mat1%nrow=0
  mat1%ncol=0

  if (allocated(mat1%nzval)) then
     call log_deallocate(mat1%nzval)
     call log_deallocate(mat1%index_i)
     call log_deallocate(mat1%index_j)
  end if

  if (present(mat2)) then
     mat2%nnz=0
     mat2%nrow=0
     mat2%ncol=0
     if (allocated(mat2%nzval)) then
        call log_deallocate(mat2%nzval)
        call log_deallocate(mat2%index_i)
        call log_deallocate(mat2%index_j)
     end if
  else
     return
  endif

  if (present(mat3)) then
     mat3%nnz=0
     mat3%nrow=0
     mat3%ncol=0
     if (allocated(mat3%nzval)) then
        call log_deallocate(mat3%nzval)
        call log_deallocate(mat3%index_i)
        call log_deallocate(mat3%index_j)
     endif
  else
     return
  endif

  if (present(mat4)) then
     mat4%nnz=0
     mat4%nrow=0
     mat4%ncol=0
     if (allocated(mat4%nzval)) then
        call log_deallocate(mat4%nzval)
        call log_deallocate(mat4%index_i)
        call log_deallocate(mat4%index_j)
     end if
  else
     return
  endif

  if (present(mat5)) then
     mat5%nnz=0
     mat5%nrow=0
     mat5%ncol=0
     if (allocated(mat5%nzval)) then
        call log_deallocate(mat5%nzval)
        call log_deallocate(mat5%index_i)
        call log_deallocate(mat5%index_j)
     end if
  else
     return
  endif

  if (present(mat6)) then
     mat6%nnz=0
     mat6%nrow=0
     mat6%ncol=0
     if (allocated(mat6%nzval)) then
        call log_deallocate(mat6%nzval)
        call log_deallocate(mat6%index_i)
        call log_deallocate(mat6%index_j)
     end if
  else
     return
  endif

  if (present(mat7)) then
     mat7%nnz=0
     mat7%nrow=0
     mat7%ncol=0
     if (allocated(mat7%nzval)) then
        call log_deallocate(mat7%nzval)
        call log_deallocate(mat7%index_i)
        call log_deallocate(mat7%index_j)
     endif
  else
     return
  end if

  if (present(mat8)) then
     mat8%nnz=0
     mat8%nrow=0
     mat8%ncol=0
     if (allocated(mat8%nzval)) then
        call log_deallocate(mat8%nzval)
        call log_deallocate(mat8%index_i)
        call log_deallocate(mat8%index_j)
     endif
  else
     return
  end if


end subroutine rdestroy_COO

! ------------------------------------------------------------------

subroutine rcreate_DNS(mat,nrow,ncol)

 type(r_DNS) :: mat
 integer :: nrow, ncol, ierr
  

  if(nrow.eq.0.or.ncol.eq.0) STOP 'ERROR: nrow or ncol = 0'

  mat%ncol=ncol
  mat%nrow=nrow
  call log_allocate(mat%val,nrow,ncol)

end subroutine rcreate_DNS
! ------------------------------------------------------------------

subroutine rdestroy_DNS(mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat8)
  
  type(r_DNS) :: mat1
  type(r_DNS), optional :: mat2,mat3,mat4,mat5,mat6,mat7,mat8
  
  mat1%nrow=0
  mat1%ncol=0

  if (allocated(mat1%val)) then
     call log_deallocate(mat1%val)
  end if

  if (present(mat2)) then
     mat2%nrow=0
     mat2%ncol=0
     if (allocated(mat2%val)) then
        call log_deallocate(mat2%val)
     end if
  else
     return
  endif

  if (present(mat3)) then
     mat3%nrow=0
     mat3%ncol=0
     if (allocated(mat3%val)) then
        call log_deallocate(mat3%val)
     endif
  else
     return
  endif

  if (present(mat4)) then
     mat4%nrow=0
     mat4%ncol=0
     if (allocated(mat4%val)) then
        call log_deallocate(mat4%val)
     end if
  else
     return
  endif

  if (present(mat5)) then
     mat5%nrow=0
     mat5%ncol=0
     if (allocated(mat5%val)) then
        call log_deallocate(mat5%val)
     end if
  else
     return
  endif

  if (present(mat6)) then
     mat6%nrow=0
     mat6%ncol=0
     if (allocated(mat6%val)) then
        call log_deallocate(mat6%val)
     end if
  else
     return
  endif

  if (present(mat7)) then
     mat7%nrow=0
     mat7%ncol=0
     if (allocated(mat7%val)) then
        call log_deallocate(mat7%val)
     endif
  else
     return
  end if

  if (present(mat8)) then
     mat8%nrow=0
     mat8%ncol=0
     if (allocated(mat8%val)) then
        call log_deallocate(mat8%val)
     endif
  else
     return
  end if

end subroutine rdestroy_DNS
! -----------------------------------------------------------------

subroutine rprint_COO(id,sp)

  type(r_COO) :: sp
  integer :: id,i
  
  do i=1,sp%nnz
    write(id,*) sp%index_i(i),sp%index_j(i),sp%nzval(i)
  enddo  

end subroutine rprint_COO
! ------------------------------------------------------------------

subroutine zwriteMem_CSR(id,mat)

  type(z_CSR) :: mat
  integer :: id,dec
  integer(8) :: mem
  character(3) :: str

  mem=mat%nnz*(2*dp+4)+(mat%nrow+2)*4
  call memstr(mem,dec,str)
  write(id,*) 'Memory Used=',mem/dec,str

end subroutine zwriteMem_CSR

! -----------------------------------------------------------------




end module mat_def
