!----------------------------------------------------------------------
! libNEGF  
!
! Source codes for sparse matrix operations and algebra 
!----------------------------------------------------------------------
! Source codes originally coded in Sparsekit
!           Y. Saad, Sep. 21 1989                                      
!----------------------------------------------------------------------
! 
! Re-written in modern Fortran in 2022
! A. Pecchia 
! D. Soccodato
!
!----------------------------------------------------------------------

module skit_formats
  use ln_precision
  implicit none
  private

  public :: dnscsr
  public :: csrdns
  public :: submat
  
  interface dnscsr
    module procedure :: rdnscsr  
    module procedure :: zdnscsr  
  end interface dnscsr

  interface csrdns
    module procedure :: rcsrdns  
    module procedure :: zcsrdns    
  end interface csrdns

  interface submat
    module procedure :: rsubmat    
    module procedure :: zsubmat    
  end interface submat
        
  contains
  !-----------------------------------------------------------------------
  ! Dense   to    Compressed Row Sparse 
  !----------------------------------------------------------------------- 
  !
  ! converts a densely stored matrix into a row orientied
  ! compactly sparse matrix. ( reverse of csrdns )
  ! Note: this routine does not check whether an element 
  ! is small. It considers that a(i,j) is zero if it is exactly
  ! equal to zero: see test below.
  !-----------------------------------------------------------------------
  ! on entry:
  !---------
  !
  ! nrow  = row-dimension of a
  ! ncol  = column dimension of a
  ! nzmax = maximum number of nonzero elements allowed. This
  !         should be set to be the lengths of the arrays a and ja.
  ! dns   = input nrow x ncol (dense) matrix.
  ! ndns  = first dimension of dns. 
  !
  ! on return:
  !---------- 
  ! 
  ! a, ja, ia = value, column, pointer  arrays for output matrix 
  !
  ! ierr  = integer error indicator: 
  !         ierr .eq. 0 means normal retur
  !         ierr .eq. i means that the the code stopped while
  !         processing row number i, because there was no space left in
  !         a, and ja (as defined by parameter nzmax).
  !----------------------------------------------------------------------- 
  subroutine rdnscsr(nrow,ncol,nzmax,dns,ndns,a,ja,ia,ierr)
    integer, intent(in) :: nrow, ncol, nzmax    
    integer, intent(in) :: ndns    
    real(dp), intent(in) :: dns(ndns,*)
    real(dp), intent(inout) :: a(:)
    integer, intent(inout) :: ia(:),ja(:)
    integer, intent(out) :: ierr  

    integer :: i, j, next

    ierr = 0
    next = 1
    ia(1) = 1
    do i = 1, nrow
      do j = 1, ncol 
         if (abs(dns(i,j)) .eq. 0.0_dp) cycle
         if (next .gt. nzmax) then
            ierr = i
            return
         endif
         ja(next) = j
         a(next) = dns(i,j)
         next = next+1
       end do    
       ia(i+1) = next
     end do
  end subroutine rdnscsr
  !----------------------------------------------------------------------- 
  subroutine zdnscsr(nrow,ncol,nzmax,dns,ndns,a,ja,ia,ierr)
    integer, intent(in) :: nrow, ncol, nzmax    
    integer, intent(in) :: ndns    
    complex(dp), intent(in) :: dns(ndns,*)
    complex(dp), intent(inout) :: a(:)
    integer, intent(inout) :: ia(:),ja(:)
    integer, intent(out) :: ierr  

    integer :: i, j, next

    ierr = 0
    next = 1
    ia(1) = 1
    do i = 1, nrow
      do j = 1, ncol 
         if (abs(dns(i,j)) .eq. 0.0_dp) cycle
         if (next .gt. nzmax) then
            ierr = i
            return
         endif
         ja(next) = j
         a(next) = dns(i,j)
         next = next+1
       end do    
       ia(i+1) = next
     end do
  end subroutine zdnscsr
  
  !-----------------------------------------------------------------------
  ! Compressed Sparse Row    to    Dense 
  !-----------------------------------------------------------------------
  !
  ! converts a row-stored sparse matrix into a densely stored one
  !
  ! On entry:
  !---------- 
  !
  ! nrow  = row-dimension of a
  ! ncol  = column dimension of a
  ! a, 
  ! ja, 
  ! ia    = input matrix in compressed sparse row format. 
  !         (a=value array, ja=column array, ia=pointer array)
  ! dns   = array where to store dense matrix
  ! ndns  = first dimension of array dns 
  !
  ! on return: 
  !----------- 
  ! dns   = the sparse matrix a, ja, ia has been stored in dns(ndns,*)
  ! 
  ! ierr  = integer error indicator. 
  !         ierr .eq. 0  means normal return
  !         ierr .eq. i  means that the code has stopped when processing
  !         row number i, because it found a column number .gt. ncol.
  ! 
  !----------------------------------------------------------------------- 
  subroutine rcsrdns(nrow, ncol,a, ja, ia, dns, ndns, ierr) 
    integer, intent(in) :: nrow, ncol
    real(dp), intent(in) :: a(:)
    integer, intent(in) :: ja(:), ia(:)
    real(dp), intent(inout) :: dns(ndns,*)
    integer, intent(out) :: ierr

    integer :: i, j, k

    ierr = 0
    do i=1, nrow
       do j=1, ncol
          dns(i,j) = 0.0_dp
       end do
    end do

    do i = 1, nrow
      do k = ia(i), ia(i+1)-1
         j = ja(k) 
         if (j .gt. ncol) then
            ierr = i
            return
         endif
         dns(i,j) = a(k)
      end do    
    end do

  end subroutine rcsrdns

  !----------------------------------------------------------------------- 
  subroutine zcsrdns(nrow, ncol,a, ja, ia, dns, ndns, ierr) 
    integer, intent(in) :: nrow, ncol
    complex(dp), intent(in) :: a(:)
    integer, intent(in) :: ja(:), ia(:)
    complex(dp), intent(inout) :: dns(ndns,*)
    integer, intent(out) :: ierr

    integer :: i, j, k

    ierr = 0
    do i=1, nrow
       do j=1, ncol
          dns(i,j) = 0.0_dp
       end do
    end do

    do i = 1, nrow
      do k = ia(i), ia(i+1)-1
         j = ja(k) 
         if (j .gt. ncol) then
            ierr = i
            return
         endif
         dns(i,j) = a(k)
      end do    
    end do

  end subroutine zcsrdns


  !-----------------------------------------------------------------------
  ! extracts the submatrix A(i1:i2,j1:j2) and puts the result in 
  ! matrix ao,iao,jao
  !---- In place: ao,jao,iao may be the same as a,ja,ia.
  !-------------- 
  ! on input
  !---------
  ! n = row dimension of the matrix 
  ! i1,i2 = two integers with i2 .ge. i1 indicating the range of rows to be
  !          extracted. 
  ! j1,j2 = two integers with j2 .ge. j1 indicating the range of columns 
  !         to be extracted.
  !         * There is no checking whether the input values for i1, i2, j1,
  !           j2 are between 1 and n. 
  ! a,
  ! ja,
  ! ia    = matrix in compressed sparse row format. 
  !
  ! job = job indicator: if job .ne. 1 then the real values in a are NOT
  !         extracted, only the column indices (i.e. data structure) are.
  !         otherwise values as well as column indices are extracted...
  !         
  ! on output
  !-------------- 
  ! nr  = number of rows of submatrix 
  ! nc  = number of columns of submatrix 
  !   * if either of nr or nc is nonpositive the code will quit.
  !
  ! ao,
  ! jao,iao = extracted matrix in general sparse format with jao containing
  ! the column indices,and iao being the pointer to the beginning 
  ! of the row,in arrays a,ja.
  subroutine rsubmat(n, job, i1, i2, j1, j2, a, ja, ia, nr, nc, ao, jao, iao)  
    integer, intent(in) :: n, job, i1, i2, j1, j2 
    integer, intent(in) :: ia(:),ja(:)
    real(dp), intent(in) :: a(:)
    integer, intent(out) :: nr, nc 
    real(dp), intent(inout) :: ao(:)
    integer, intent(inout) :: jao(:),iao(:)

    integer :: i, j, k, k1, k2, klen
    nr = i2-i1+1
    nc = j2-j1+1
    if ( nr .le. 0 .or. nc .le. 0) return
    klen = 0
    ! 
    ! simple procedure. proceeds row-wise...
    ! 
    do i = 1, nr
       ii = i1+i-1
       k1 = ia(ii)
       k2 = ia(ii+1)-1
       iao(i) = klen+1
       do k = k1, k2
          j = ja(k)
          if (j .ge. j1 .and. j .le. j2) then
             klen = klen+1
             if (job .eq. 1) ao(klen) = a(k)
             jao(klen) = j - j1+1
          endif
       end do
    end do 
    iao(nr+1) = klen+1
      
  end subroutine rsubmat

  subroutine zsubmat(n, job, i1, i2, j1, j2, a, ja, ia, nr, nc, ao, jao, iao)  
    integer, intent(in) :: n, job, i1, i2, j1, j2 
    integer, intent(in) :: ia(:),ja(:)
    complex(dp), intent(in) :: a(:)
    integer, intent(out) :: nr, nc 
    complex(dp), intent(inout) :: ao(:)
    integer, intent(inout) :: jao(:),iao(:)

    integer :: i, j, k, k1, k2, klen
    nr = i2-i1+1
    nc = j2-j1+1
    if ( nr .le. 0 .or. nc .le. 0) return
    klen = 0
    ! 
    ! simple procedure. proceeds row-wise...
    ! 
    do i = 1, nr
       ii = i1+i-1
       k1 = ia(ii)
       k2 = ia(ii+1)-1
       iao(i) = klen+1
       do k = k1, k2
          j = ja(k)
          if (j .ge. j1 .and. j .le. j2) then
             klen = klen+1
             if (job .eq. 1) ao(klen) = a(k)
             jao(klen) = j - j1+1
          endif
       end do
    end do 
    iao(nr+1) = klen+1
      
  end subroutine zsubmat

end module skit_formats

! FORMAT OPERATIONS
! -----------------------
!    subroutine zcoocsr(nrow, nnz, aa, ia, ja, bb, jb, ib)
!    
!    subroutine zcsrcoo(nrow, job, nnz, aa, ja, ia, innz, bb, ib, jb, ierr)
!    
!    subroutine zcsrcsc(nrow, job, ipos, aa, colind, rowpnt, bb, rowind, colpnt)
!
!    subroutine zcsrcsc2(nrow, ncol, job, ipos, aa, colind, rowpnt, bb, rowind, colpnt)
!    

! UNARY OPERATIONS
! -----------------------
!    function getelm(ii, jj, aa, ja, ia, iadd, sorted) result(res)
!
!    function zgetelm(ii, jj, aa, ja, ia, iadd, sorted) result(res)
!       
!    subroutine zcsort(nrow, aa, ja, ia, iw, vals)
!
!    subroutine ztransp(nrow,ncol,aa,ja,ia,iwk,ierr)
!    
!    subroutine bandwidth(ncol, ja, ia, ml, mu, a_bw, bndav)
!    
!    subroutine getdia(nrow, ncol, job, aa, ja, ia, len, dd, idiag, ioff)


