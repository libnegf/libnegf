module input_output
  use ln_precision
  use ln_allocation
  use mat_def
  use sparsekit_drv
  implicit none
  private

  public :: read_H, format

  type format
     character(5) :: type ! 'PETSc','UPT'  
     logical :: formatted ! formatted/unformatted
     character(1) :: fmt  ! U,L,F
  end type format

  contains

subroutine read_H(idR,idI, zmat,fmt)
  Integer :: idR, idI 
  Type(z_CSR) :: zmat
  Type(format) :: fmt

  !locals
  Type(r_COO) :: matr, mati
  Type(z_COO) :: mat2
  Integer :: i, count, k, nnz
  Character(10) :: tmp1, tmp2, tmp3, tmp4

  if (trim(fmt%type).eq.'PETSc') then

     read (idR,*) tmp1, tmp2, tmp3, matr%nrow, matr%ncol
     read (idR,*) tmp1, tmp4, tmp3, nnz
     read (idR,*)
     read (idR,*)
     !write(*,*) 'The number of rows (Hreal) is', matr%nrow
     !write(*,*) 'The number of columns (Hreal) is', matr%ncol
     !write(*,*) 'The number of non zero elements (Hreal) is', nnz

  elseif (trim(fmt%type).eq.'UPT') then

     if (fmt%formatted) then
        read(idR,*) matr%nrow, nnz
     else
        read(idR) matr%nrow, nnz
     endif
     !write(*,*) 'The number of rows (Hreal) is', matr%nrow
     !write(*,*) 'The number of non zero elements (Hreal) is', nnz

  end if

  select case(fmt%fmt)
  case('F')
     call create(matr,matr%nrow, matr%nrow, nnz)
     call create(mati,matr%nrow, matr%nrow, nnz)
  case('L','U')
     call create(matr,matr%nrow, matr%nrow, 2*nnz-matr%nrow)
     call create(mati,matr%nrow, matr%nrow, 2*nnz-matr%nrow) 
  end select

  if (fmt%formatted) then
     if (trim(fmt%type).eq.'PETSc') then  

        k = 0
        do i=1,nnz
           k=k+1
           read(idR,*) matr%index_j(k), matr%index_i(k), matr%nzval(k)
      
           select case(fmt%fmt)
           case('U','L') 
              if(matr%index_j(k).ne. matr%index_i(k)) then
                 k = k + 1
                 matr%index_i(k) =  matr%index_j(k-1)
                 matr%index_j(k) =  matr%index_i(k-1)
                 matr%nzval(k) = matr%nzval(k-1) 
              endif
           case('F')
           end select
        enddo

     elseif (trim(fmt%type).eq.'UPT') then
        print*,'Read Formatted UPT'
        k = 0
        do i=1,nnz
           k = k + 1           
           read(idR,*) matr%index_i(k), matr%index_j(k), matr%nzval(k)
           
           select case(fmt%fmt)
           case('U','L') 
              if(matr%index_j(k).ne. matr%index_i(k)) then
                 k = k + 1
                 matr%index_i(k) =  matr%index_j(k-1)
                 matr%index_j(k) =  matr%index_i(k-1)
                 matr%nzval(k) = matr%nzval(k-1)
              endif
           case('F')
           end select
        enddo     

     endif

  else
     if (trim(fmt%type).eq.'PETSc') then  

        k = 0
        do i=1,nnz

           k = k + 1
           read(idR) matr%index_j(k), matr%index_i(k), matr%nzval(k)
          
           select case(fmt%fmt)
           case('U','L') 
              if(matr%index_j(k).ne. matr%index_i(k)) then
                 k = k + 1
                 matr%index_i(k) =  matr%index_j(k-1)
                 matr%index_j(k) =  matr%index_i(k-1)
                 matr%nzval(k) = matr%nzval(k-1)
              endif
           case('F')
           end select
        enddo

     elseif (trim(fmt%type).eq.'UPT') then
        k = 0
        do i=1,nnz
           k = k + 1
           read(idR) matr%index_i(k), matr%index_j(k), matr%nzval(k)
           
           select case(fmt%fmt)
           case('U','L') 
              if(matr%index_j(k).ne. matr%index_i(k)) then
                 k = k + 1
                 matr%index_i(k) =  matr%index_j(k-1)
                 matr%index_j(k) =  matr%index_i(k-1)
                 matr%nzval(k) = matr%nzval(k-1)
              endif
           case('F')
           end select
        enddo     

     endif

  endif


  !write(*,*) matr%index_j, matr%index_i, matr%nzval
  ! =============================================================================
  !  Imaginary part
  ! =============================================================================
  if (trim(fmt%type).eq.'PETSc') then

     read (idI,*) tmp1, tmp2, tmp3, mati%nrow, mati%ncol
     read (idI,*) tmp1, tmp4, tmp3, nnz
     read (idI,*)
     read (idI,*)
     !write(*,*) 'The number of rows (Himm) is', mati%nrow
     !write(*,*) 'The number of columns (Himm) is', mati%ncol
     !write(*,*) 'The number of non zero elements (Himm) is', nnz

  elseif (trim(fmt%type).eq.'UPT') then

     if (fmt%formatted) then
        read(idI,*) matr%nrow, nnz
     else
        read(idI) matr%nrow, nnz
     endif
     !write(*,*) 'The number of rows (Himm) is', matr%nrow
     !write(*,*) 'The number of non zero elements (Hreal) is', nnz

  end if
  
  if(mati%nrow.ne.matr%nrow) stop 'nrow Error'  
  if(matr%nnz.ne.mati%nnz) stop 'nnz Error' 
  if(matr%nnz.eq.0) stop 'nnz Error'   

  k = 0
  do i=1,nnz

     k = k + 1
     if (fmt%formatted) then   
        if (trim(fmt%type).eq.'PETSc') then    
           read(idI,*) mati%index_j(k), mati%index_i(k), mati%nzval(k)
        elseif (trim(fmt%type).eq.'UPT') then
           read(idI,*) mati%index_i(k), mati%index_j(k), mati%nzval(k)
        endif
     else
        if (trim(fmt%type).eq.'PETSc') then
           read(idI) mati%index_j(k), mati%index_i(k), mati%nzval(k)
        elseif (trim(fmt%type).eq.'UPT') then
           read(idI) mati%index_i(k), mati%index_j(k), mati%nzval(k)
        endif
     endif

     if(mati%index_j(k).ne.matr%index_j(k)) stop 'Index Error'
     if(mati%index_i(k).ne.matr%index_i(k)) stop 'Index Error'

     select case(fmt%fmt)
     case('U','L') 
        if(mati%index_j(k).ne. mati%index_i(k)) then
           k = k + 1
           mati%index_i(k) =  mati%index_j(k-1)
           mati%index_j(k) =  mati%index_i(k-1)
           mati%nzval(k) = - mati%nzval(k-1)
        endif
     case('F')
     end select
  enddo
  !write(*,*) mat%index_j, mat%index_i, mati%nzval

  !write(*,*) 'Hcomplex is'
  
  count = k
  !do i = 1, k  
  !   if(  abs(matr%nzval(i)).gt.EPS .or. abs(mati%nzval(i)).gt.EPS  ) then
  !      count = count + 1
  !   endif
  !enddo
   
  call create(mat2,matr%nrow,matr%nrow,count)

  count = 0
  do i = 1, k  
  !   if(  abs(matr%nzval(i)).gt.EPS .or. abs(mati%nzval(i)).gt.EPS  ) then
       count = count + 1
       mat2%index_i(count) = matr%index_i(i)
       mat2%index_j(count) = matr%index_j(i)
       mat2%nzval(count)= cmplx(matr%nzval(i), mati%nzval(i), dp)        
  !   endif
  enddo
  

  call destroy(matr,mati)

  !COO-CSR conversion 
  CALL create(zmat,mat2%nrow,mat2%ncol,mat2%nnz)

  zmat%nzval=(0.d0,0.d0)

  CALL coo2csr(mat2,zmat)

  call destroy(mat2)               !deallocation Hcomplex

  !write(*,*) '(readH) matrix H read'


end subroutine read_H

! -----------------------------------------------------------------




end module input_output
