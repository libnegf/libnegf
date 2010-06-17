module load
 
  use mat_def
  implicit none
  private

  public :: load_HS


contains

  subroutine load_HS(H,S,fmt)
    type(z_CSR) :: H,S  
    integer :: nrow, ncol, nnz
    logical :: fmt

    ! --------------------------------------------
    if (fmt) then
       open(11, file='H.dat', FORM='FORMATTED')
       read(11,*) nrow
       read(11,*) ncol  
       read(11,*) nnz
       call create(H,nrow,ncol,nnz)
       
       rewind 11
       call read_mat(11,H,.true.)
       close(11)
    else
       open(11, file='H.dat', FORM='UNFORMATTED')
       read(11) nrow
       read(11) ncol  
       read(11) nnz
       call create(H,nrow,ncol,nnz)
       
       rewind 11
       call read_mat(11,H,.false.)
       close(11)
    endif


    if (fmt) then
       open(11, file='S.dat', FORM='FORMATTED')
       read(11,*) nrow
       read(11,*) ncol  
       read(11,*) nnz
       call create(S,nrow,ncol,nnz)
       
       rewind 11
       call read_mat(11,S,.true.)
       close(11)  
    else
       open(11, file='S.dat', FORM='UNFORMATTED')       
       read(11) nrow
       read(11) ncol  
       read(11) nnz
       call create(S,nrow,ncol,nnz)
       
       rewind 11
       call read_mat(11,S,.false.)
       close(11)       
    endif


  end subroutine load_HS


end module load
