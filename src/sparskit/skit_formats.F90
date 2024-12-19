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
  public :: coocsr
  public :: csrcoo
  public :: csrcsc
  public :: csrcsc2
  public :: submat
  
  interface dnscsr
    module procedure :: ddnscsr  
    module procedure :: zdnscsr  
  end interface dnscsr

  interface csrdns
    module procedure :: dcsrdns  
    module procedure :: zcsrdns    
  end interface csrdns
  
  interface coocsr
    module procedure :: dcoocsr  
    module procedure :: zcoocsr    
  end interface coocsr

  interface csrcoo
    module procedure :: dcsrcoo  
    module procedure :: zcsrcoo    
  end interface csrcoo

  interface csrcsc
    module procedure :: dcsrcsc
    module procedure :: zcsrcsc    
  end interface csrcsc

  interface csrcsc2
    module procedure :: dcsrcsc2
    module procedure :: zcsrcsc2   
  end interface csrcsc2

  interface submat
    module procedure :: dsubmat    
    module procedure :: zsubmat    
  end interface submat
        
  contains
  !-----------------------------------------------------------------------
  ! dnscsr:  Dense   to    Compressed Row Sparse 
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
  subroutine ddnscsr(nrow,ncol,nzmax,dns,a,ja,ia,ierr)
    integer, intent(in) :: nrow, ncol, nzmax    
    real(dp), intent(in) :: dns(:,:)
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
  end subroutine ddnscsr
  !----------------------------------------------------------------------- 
  subroutine zdnscsr(nrow,ncol,nzmax,dns,a,ja,ia,ierr)
    integer, intent(in) :: nrow, ncol, nzmax    
    complex(dp), intent(in) :: dns(:,:)
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
  ! csrdns:  Compressed Sparse Row    to    Dense 
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
  subroutine dcsrdns(nrow, ncol,a, ja, ia, dns, ierr) 
    integer, intent(in) :: nrow, ncol
    real(dp), intent(in) :: a(:)
    integer, intent(in) :: ja(:), ia(:)
    real(dp), intent(inout) :: dns(:,:)
    integer, intent(out) :: ierr

    integer :: i, j, k

    ierr = 0
    dns = 0.0_dp

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

  end subroutine dcsrdns

  !----------------------------------------------------------------------- 
  subroutine zcsrdns(nrow, ncol,a, ja, ia, dns, ierr) 
    integer, intent(in) :: nrow, ncol
    complex(dp), intent(in) :: a(:)
    integer, intent(in) :: ja(:), ia(:)
    complex(dp), intent(inout) :: dns(:,:)
    integer, intent(out) :: ierr

    integer :: i, j, k

    ierr = 0
    dns = (0.0_dp, 0.0_dp)

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
  ! coocsr:   Coordinate     to   Compressed Sparse Row 
  !----------------------------------------------------------------------- 
  ! converts a matrix that is stored in coordinate format
  !  a, ir, jc into a row general sparse ao, jao, iao format.
  !
  ! on entry:
  !--------- 
  ! nrow  = dimension of the matrix 
  ! nnz = number of nonzero elements in matrix
  ! a,
  ! ir, 
  ! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
  !         nonzero elements of the matrix with a(k) = actual real value of
  !     the elements, ir(k) = its row number and jc(k) = its column 
  !   number. The order of the elements is arbitrary. 
  !
  ! on return:
  !----------- 
  ! ir  is destroyed
  !
  ! ao, jao, iao = matrix in general sparse matrix format with ao 
  !   continung the real values, jao containing the column indices, 
  ! and iao being the pointer to the beginning of the row, 
  ! in arrays ao, jao.
  !
  ! Notes:
  !------ This routine is NOT in place.  See coicsr
  !
  !------------------------------------------------------------------------
  subroutine dcoocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
    integer, intent(in) :: nrow, nnz   
    real(dp), intent(in) :: a(:) 
    integer, intent(in) :: ir(:), jc(:)
    real(dp), intent(inout) :: ao(:)
    integer, intent(inout) :: jao(:), iao(:)

    integer :: i, j, k, k0, iad
    real(dp) :: x

    iao = 0
    do k = 1, nnz
         iao(ir(k)) = iao(ir(k))+1
    end do 
      
    k = 1
    do j = 1,nrow+1
       k0 = iao(j)
       iao(j) = k
       k = k+k0
    end do
    do k = 1, nnz
       i = ir(k)
       j = jc(k)
       x = a(k)
       iad = iao(i)
       ao(iad) =  x
       jao(iad) = j
       iao(i) = iad+1
    end do
    do j = nrow,1,-1
       iao(j+1) = iao(j)
    end do
    iao(1) = 1

  end subroutine dcoocsr

  subroutine zcoocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
    integer, intent(in) :: nrow, nnz   
    complex(dp), intent(in) :: a(:) 
    integer, intent(in) :: ir(:), jc(:)
    complex(dp), intent(inout) :: ao(:)
    integer, intent(inout) :: jao(:), iao(:)

    integer :: i, j, k, k0, iad
    complex(dp) :: x

    iao = 0
    do k = 1, nnz
         iao(ir(k)) = iao(ir(k))+1
    end do 
      
    k = 1
    do j = 1,nrow+1
       k0 = iao(j)
       iao(j) = k
       k = k+k0
    end do
    do k = 1, nnz
       i = ir(k)
       j = jc(k)
       x = a(k)
       iad = iao(i)
       ao(iad) =  x
       jao(iad) = j
       iao(i) = iad+1
    end do
    do j = nrow,1,-1
       iao(j+1) = iao(j)
    end do
    iao(1) = 1
  end subroutine zcoocsr

  !----------------------------------------------------------------------- 
  ! csrcoo:  Compressed Sparse Row      to      Coordinate 
  !----------------------------------------------------------------------- 
  ! converts a matrix that is stored in coordinate format
  !  a, ir, jc into a row general sparse ao, jao, iao format.
  !
  ! on entry: 
  !---------
  ! nrow  = dimension of the matrix.
  ! job   = integer serving as a job indicator. 
  !         if job = 1 fill in only the array ir, ignore jc, and ao.
  !         if job = 2 fill in ir, and jc but not ao 
  !         if job = 3 fill in everything.
  !         The reason why these options are provided is that on return 
  !         ao and jc are the same as a, ja. So when job = 3, a and ja are
  !         simply copied into ao, jc.  When job=2, only jc and ir are
  !         returned. With job=1 only the array ir is returned. Moreover,
  !         the algorithm is in place:
  !      call csrcoo (nrow,1,nzmax,a,ja,ia,nnz,a,ia,ja,ierr) 
  !         will write the output matrix in coordinate format on a, ja,ia.
  !
  ! a,
  ! ja,
  ! ia    = matrix in compressed sparse row format.
  ! nzmax = length of space available in ao, ir, jc.
  !         the code will stop immediatly if the number of
  !         nonzero elements found in input matrix exceeds nzmax.
  ! 
  ! on return:
  !----------- 
  ! ao, ir, jc = matrix in coordinate format.
  !
  ! nnz        = number of nonzero elements in matrix.
  ! ierr       = integer error indicator.
  !         ierr .eq. 0 means normal retur
  !         ierr .eq. 1 means that the the code stopped 
  !         because there was no space in ao, ir, jc 
  !         (according to the value of  nzmax).
  ! 
  ! NOTES: 1)This routine is PARTIALLY in place: csrcoo can be called with 
  !         ao being the same array as as a, and jc the same array as ja. 
  !         but ir CANNOT be the same as ia. 
  !         2) note the order in the output arrays, 
  !------------------------------------------------------------------------
  subroutine dcsrcoo(nrow,job,nzmax,a,ja,ia,nnz,ao,ir,jc,ierr)
    integer, intent(in) :: nrow, job, nzmax   
    real(dp), intent(in) :: a(:) 
    integer, intent(in) :: ia(:), ja(:)
    real(dp), intent(inout) :: ao(:)
    integer, intent(inout) :: ir(:), jc(:)
    integer, intent(out) :: nnz, ierr
    
    integer :: i, k1, k2, k

    ierr = 0
    nnz = ia(nrow+1)-1
    if (nnz .gt. nzmax) then
       ierr = 1
       return
    endif

    !  copy backward to allow for in-place processing. 
    do i=nrow,1,-1
       k1 = ia(i+1)-1
       k2 = ia(i)
       do k=k1,k2,-1
          ir(k) = i
       end do
    end do

    if (job >= 3) then
      do k=1,nnz
        ao(k) = a(k)
      end do 
    end if  
    if (job >= 2) then
      do k=1,nnz
        jc(k) = ja(k)
      end do  
    end if  
  end subroutine dcsrcoo

  subroutine zcsrcoo(nrow,job,nzmax,a,ja,ia,nnz,ao,ir,jc,ierr)
    integer, intent(in) :: nrow, job, nzmax   
    complex(dp), intent(in) :: a(:) 
    integer, intent(in) :: ia(:), ja(:)
    complex(dp), intent(inout) :: ao(:)
    integer, intent(inout) :: ir(:), jc(:)
    integer, intent(out) :: nnz, ierr
    
    integer :: i, k1, k2, k

    ierr = 0
    nnz = ia(nrow+1)-1
    if (nnz .gt. nzmax) then
       ierr = 1
       return
    endif

    !  copy backward to allow for in-place processing. 
    do i=nrow,1,-1
       k1 = ia(i+1)-1
       k2 = ia(i)
       do k=k1,k2,-1
          ir(k) = i
       end do
    end do

    if (job >= 3) then
      do k=1,nnz
        ao(k) = a(k)
      end do 
    end if  
    if (job >= 2) then
      do k=1,nnz
        jc(k) = ja(k)
      end do  
    end if  
  end subroutine zcsrcoo

  !-----------------------------------------------------------------------
  ! Compressed Sparse Row     to      Compressed Sparse Column
  !
  ! (transposition operation)   Not in place. 
  !----------------------------------------------------------------------- 
  ! -- not in place --
  ! this subroutine transposes a matrix stored in a, ja, ia format.
  ! ---------------
  ! on entry:
  !----------
  ! n = dimension of A.
  ! job = integer to indicate whether to fill the values (job.eq.1) of the
  !         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
  !
  ! ipos  = starting position in ao, jao of the transposed matrix.
  !         the iao array takes this into account (thus iao(1) is set to ipos.)
  !         Note: this may be useful if one needs to append the data structure
  !         of the transpose to that of A. In this case use for example
  !                call csrcsc (n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2)) 
  !   for any other normal usage, enter ipos=1.
  ! a = real array of length nnz (nnz=number of nonzero elements in input 
  !         matrix) containing the nonzero elements.
  ! ja  = integer array of length nnz containing the column positions
  !     of the corresponding elements in a.
  ! ia  = integer of size n+1. ia(k) contains the position in a, ja of
  !   the beginning of the k-th row.
  !
  ! on return:
  ! ---------- 
  ! output arguments:
  ! ao  = real array of size nzz containing the "a" part of the transpose
  ! jao = integer array of size nnz containing the column indices.
  ! iao = integer array of size n+1 containing the "ia" index array of
  !   the transpose. 
  !
  !-----------------------------------------------------------------------     
  subroutine dcsrcsc(n,job,ipos,a,ja,ia,ao,jao,iao)
    integer, intent(in) :: n, job, ipos
    real(dp), intent(in) :: a(:)
    integer, intent(in) :: ja(:), ia(:)
    real(dp), intent(inout) :: ao(:)
    integer, intent(inout) :: jao(:), iao(:)

    call dcsrcsc2(n,n,job,ipos,a,ja,ia,ao,jao,iao)
  end subroutine dcsrcsc
  
  subroutine zcsrcsc(n,job,ipos,a,ja,ia,ao,jao,iao)
    integer, intent(in) :: n, job, ipos
    complex(dp), intent(in) :: a(:)
    integer, intent(in) :: ja(:), ia(:)
    complex(dp), intent(inout) :: ao(:)
    integer, intent(inout) :: jao(:), iao(:)

    call zcsrcsc2(n,n,job,ipos,a,ja,ia,ao,jao,iao)
  end subroutine zcsrcsc

  !-----------------------------------------------------------------------
  ! Compressed Sparse Row     to      Compressed Sparse Column
  !
  ! (transposition operation)   Not in place. 
  !----------------------------------------------------------------------- 
  ! Rectangular version.  n is number of rows of CSR matrix,
  !                       n2 (input) is number of columns of CSC matrix.
  !----------------------------------------------------------------------- 
  ! -- not in place --
  ! this subroutine transposes a matrix stored in a, ja, ia format.
  ! ---------------
  ! on entry:
  !----------
  ! n = number of rows of CSR matrix.
  ! n2    = number of columns of CSC matrix.
  ! job = integer to indicate whether to fill the values (job.eq.1) of the
  !         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
  !
  ! ipos  = starting position in ao, jao of the transposed matrix.
  !         the iao array takes this into account (thus iao(1) is set to ipos.)
  !         Note: this may be useful if one needs to append the data structure
  !         of the transpose to that of A. In this case use for example
  !                call csrcsc2 (n,n,1,ia(n+1),a,ja,ia,a,ja,ia(n+2)) 
  !   for any other normal usage, enter ipos=1.
  ! a = real array of length nnz (nnz=number of nonzero elements in input 
  !         matrix) containing the nonzero elements.
  ! ja  = integer array of length nnz containing the column positions
  !     of the corresponding elements in a.
  ! ia  = integer of size n+1. ia(k) contains the position in a, ja of
  !   the beginning of the k-th row.
  !
  ! on return:
  ! ---------- 
  ! output arguments:
  ! ao  = real array of size nzz containing the "a" part of the transpose
  ! jao = integer array of size nnz containing the column indices.
  ! iao = integer array of size n+1 containing the "ia" index array of
  !   the transpose. 
  !
  !----------------------------------------------------------------------- 
  subroutine dcsrcsc2(n,n2,job,ipos,a,ja,ia,ao,jao,iao)
    integer, intent(in) :: n, n2, job, ipos
    real(dp), intent(in) :: a(:)
    integer, intent(in) :: ja(:), ia(:)
    real(dp), intent(inout) :: ao(:)
    integer, dimension(:), intent(inout) :: jao, iao

    integer :: i, j, k, next

    !----- compute lengths of rows of transp(A) ----------------
    do i = 1, n2+1
       iao(i) = 0
    end do
    do i = 1, n
       do k = ia(i), ia(i+1)-1 
          j = ja(k)+1
          iao(j) = iao(j)+1
       end do
    end do
    !----- compute pointers from lengths -----------------------
    iao(1) = ipos 
    do i= 1, n2
      iao(i+1) = iao(i) + iao(i+1)
    end do
    !---- now do the actual copying ---------------------------- 
    do i = 1, n
       do k= ia(i), ia(i+1)-1 
          j = ja(k) 
          next = iao(j)
          if (job .eq. 1)  ao(next) = a(k)
          jao(next) = i
          iao(j) = next+1
       end do
    end do
    !---------- reshift iao and leave ---------------------- 
    do i = n2, 1, -1
       iao(i+1) = iao(i)
    end do
    iao(1) = ipos
  end subroutine dcsrcsc2

  subroutine zcsrcsc2(n,n2,job,ipos,a,ja,ia,ao,jao,iao)
    integer, intent(in) :: n, n2, job, ipos
    complex(dp), intent(in) :: a(:)
    integer, intent(in) :: ja(:), ia(:)
    complex(dp), intent(inout) :: ao(:)
    integer, intent(inout) :: jao(:), iao(:)

    integer :: i, j, k, next

    !----- compute lengths of rows of transp(A) ----------------
    do i = 1, n2+1
       iao(i) = 0
    end do
    do i = 1, n
       do k = ia(i), ia(i+1)-1 
          j = ja(k)+1
          iao(j) = iao(j)+1
       end do
    end do
    !----- compute pointers from lengths -----------------------
    iao(1) = ipos 
    do i= 1, n2
      iao(i+1) = iao(i) + iao(i+1)
    end do
    !---- now do the actual copying ---------------------------- 
    do i = 1, n
       do k= ia(i), ia(i+1)-1 
          j = ja(k) 
          next = iao(j)
          if (job .eq. 1)  ao(next) = a(k)
          jao(next) = i
          iao(j) = next+1
       end do
    end do
    !---------- reshift iao and leave ---------------------- 
    do i = n2, 1, -1
       iao(i+1) = iao(i)
    end do
    iao(1) = ipos
  end subroutine zcsrcsc2


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
  subroutine dsubmat(n, job, i1, i2, j1, j2, a, ja, ia, nr, nc, ao, jao, iao)  
    integer, intent(in) :: n, job, i1, i2, j1, j2 
    integer, intent(in) :: ia(:),ja(:)
    real(dp), intent(in) :: a(:)
    integer, intent(out) :: nr, nc 
    real(dp), intent(inout) :: ao(:)
    integer, intent(inout) :: jao(:),iao(:)

    integer :: i, j, k, k1, k2, klen, ii
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
      
  end subroutine dsubmat

  subroutine zsubmat(n, job, i1, i2, j1, j2, a, ja, ia, nr, nc, ao, jao, iao)  
    integer, intent(in) :: n, job, i1, i2, j1, j2 
    integer, intent(in) :: ia(:),ja(:)
    complex(dp), intent(in) :: a(:)
    integer, intent(out) :: nr, nc 
    complex(dp), intent(inout) :: ao(:)
    integer, intent(inout) :: jao(:),iao(:)

    integer :: i, j, k, k1, k2, klen, ii
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

