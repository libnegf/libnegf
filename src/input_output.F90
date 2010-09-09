module input_output
  use precision
  use allocation
  use mat_def
  use sparsekit_drv
  implicit none
  private

  public :: read_H, read_S 


  contains

subroutine read_H(idR,idI, zmat)
  Integer :: idR, idI 
  Type(z_CSR) :: zmat
  Type(r_COO) :: mat, mat1
  Type(z_COO) :: mat2
  Integer i, count
  Character(10) :: tmp1, tmp2, tmp3, tmp4


  read (idR,*) tmp1, tmp2, tmp3, mat%nrow, mat%ncol
  read (idR,*) tmp1, tmp4, tmp3, mat%nnz
  read (idR,*)
  read (idR,*)
  write(*,*) 'The number of rows (Hreal) is', mat%nrow
  write(*,*) 'The number of columns (Hreal) is', mat%ncol
  write(*,*) 'The number of non zero elements (Hreal) is', mat%nnz

  call create(mat,mat%nrow, mat%ncol, mat%nnz)
  call create(mat1,mat%nrow, mat%ncol, mat%nnz)
  
  do i=1,mat%nnz
     read(idR,*) mat%index_j(i), mat%index_i(i), mat%nzval(i)
  enddo
  !write(*,*) mat%index_j, mat%index_i, mat%nzval
  
  
  read (idI,*) tmp1, tmp2, tmp3, mat1%nrow, mat1%ncol
  read (idI,*) tmp1, tmp4, tmp3, mat1%nnz
  read (idI,*)
  read (idI,*)
  write(*,*) 'The number of rows (Himm) is', mat1%nrow
  write(*,*) 'The number of columns (Himm) is', mat1%ncol
  write(*,*) 'The number of non zero elements (Himm) is', mat1%nnz


  if(mat1%nrow.ne.mat%nrow) stop 'nrow Error'  
  if(mat%nnz.ne.mat1%nnz) stop 'nnz Error' 
  if(mat%nnz.eq.0) stop 'nnz Error'   

  do i=1,mat%nnz
     read(idI,*) mat1%index_j(i), mat1%index_i(i), mat1%nzval(i)
     if(mat1%index_j(i).ne.mat%index_j(i)) stop 'Index Error'
     if(mat1%index_i(i).ne.mat%index_i(i)) stop 'Index Error'     
  enddo
  !write(*,*) mat%index_j, mat%index_i, mat1%nzval
  
  !write(*,*) 'Hcomplex is'

  count = 0
  do i=1,mat%nnz  
     if(  abs(mat%nzval(i)).gt.EPS .or. abs(mat1%nzval(i)).gt.EPS  ) then
        count = count + 1
     endif
  enddo
   
  call create(mat2,mat%nrow,mat%ncol,count)

  count = 0
  do i=1,mat%nnz  
     if(  abs(mat%nzval(i)).gt.EPS .or. abs(mat1%nzval(i)).gt.EPS  ) then
       count = count + 1
       mat2%index_i(count) = mat%index_i(i)
       mat2%index_j(count) = mat%index_j(i)
       mat2%nzval(count)=(1.0_dp,0.0_dp)*mat%nzval(i)+(0.0_dp,1.0_dp)*mat1%nzval(i)        
     endif
  enddo
  !write(*,*) mat2%nzval                                    

  call destroy(mat,mat1)

!conversione COO-CSR 

  CALL create(zmat,mat2%nrow,mat2%ncol,mat2%nnz)

  CALL coo2csr(mat2,zmat)
  

  call destroy(mat2)               !deallocation Hcomplex

  write(*,*) 'matrix H read'

end subroutine read_H

! -----------------------------------------------------------------

subroutine read_S(id, mat)

  Type(r_COO) :: mat
  Integer i
  Integer :: id
  Character tmp1, tmp2, tmp3, tmp4

  read (id,*) tmp1, tmp2, tmp3, mat%nrow, mat%ncol
  read (id,*) tmp1, tmp4, tmp3, mat%nnz
  read (id,*)
  read (id,*)
  write(*,*) 'The number of rows (Hreal) is', mat%nrow
  write(*,*) 'The number of columns (Hreal) is', mat%ncol
  write(*,*) 'The number of non zero elements (Hreal) is', mat%nnz

  call log_allocate(mat%index_j,mat%nnz)        !allocation column index
  call log_allocate(mat%index_i,mat%nnz)        !allocation row index
  call log_allocate(mat%nzval,mat%nnz)          !allocation Sreal

  do i=1,mat%nnz
     read(id,*) mat%index_j(i), mat%index_i(i), mat%nzval(i)
  enddo
  write(*,*) mat%index_j, mat%index_i, mat%nzval

  call log_deallocate(mat%index_j)              !deallocation column index
  call log_deallocate(mat%index_i)              !deallocation row index
  call log_deallocate(mat%nzval)                !deallocation Sreal 

  !Inserire conversione COO-CSR

end subroutine read_S



end module input_output
