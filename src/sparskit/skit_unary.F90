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
module skit_unary
  use ln_precision
  implicit none
  private

  public :: getelm
  public :: csort
  public :: transp
  public :: bandwidth
  public :: getdia
  public :: amask 

  interface getelm
    module procedure :: rgetelm    
    module procedure :: zgetelm    
  end interface getelm
  
  interface csort
    module procedure :: rcsort  
    module procedure :: zcsort    
  end interface csort
  
  interface vperm
    module procedure :: ivperm  
    module procedure :: rvperm   
    module procedure :: zvperm   
  end interface vperm

  interface transp 
    module procedure :: rtransp  
    module procedure :: ztransp   
  end interface transp

  interface getdia 
    module procedure :: rgetdia
    module procedure :: zgetdia  
  end interface getdia

  interface amask 
    module procedure :: ramask 
    module procedure :: zamask
  end interface amask 

  contains

  !-----------------------------------------------------------------------
  ! getelm: 
  !-----------------------------------------------------------------------
  !     this function returns the element a(i,j) of a matrix a, 
  !     for any pair (i,j).  the matrix is assumed to be stored 
  !     in compressed sparse row (csr) format. getelm performs a
  !     binary search in the case where it is known that the elements 
  !     are sorted so that the column indices are in increasing order. 
  !     also returns (in iadd) the address of the element a(i,j) in 
  !     arrays a and ja when the search is successsful (zero if not).
  !----- 
  !     first contributed by noel nachtigal (mit). 
  !     recoded jan. 20, 1991, by y. saad [in particular
  !     added handling of the non-sorted case + the iadd output] 
  !-----------------------------------------------------------------------
  ! on entry: 
  !---------- 
  !     i      = the row index of the element sought (input).
  !     j      = the column index of the element sought (input).
  !     a      = the matrix a in compressed sparse row format (input).
  !     ja     = the array of column indices (input).
  !     ia     = the array of pointers to the rows' data (input).
  !     sorted = logical indicating whether the matrix is knonw to 
  !              have its column indices sorted in increasing order 
  !              (sorted=.true.) or not (sorted=.false.).
  !              (input). 
  ! on return:
  !----------- 
  !     getelm = value of a(i,j). 
  !     iadd   = address of element a(i,j) in arrays a, ja if found,
  !              zero if not found. (output) 
  !
  !     note: the inputs i and j are not checked for validity. 
  !-----------------------------------------------------------------------
  !     noel m. nachtigal october 28, 1990 -- youcef saad jan 20, 1991.
  !----------------------------------------------------------------------- 
  function rgetelm(i,j,a,ja,ia,iadd,sorted) result(getelm) 
    integer, intent(in) :: i,j
    integer, intent(in) :: ia(:), ja(:) 
    integer, intent(inout) :: iadd
    real(dp), intent(in) :: a(:)
    logical, intent(in) :: sorted 

    integer :: ibeg, iend, imid, k
    real(dp) :: getelm
      
    iadd = 0 
    getelm = 0.0_dp
    ibeg = ia(i)
    iend = ia(i+1)-1
      
    if (.not. sorted) then 
       do k = ibeg, iend
          if (ja(k) .eq.  j) then
             iadd = k 
             if (iadd .ne. 0) then
                getelm = a(iadd)
                return 
             end if 
          endif
       end do
    else
       ! binary search
       do
          imid = ( ibeg + iend ) / 2
          if (ja(imid).eq.j) then
             iadd = imid 
             if (iadd .ne. 0) then
                getelm = a(iadd)
                return 
             end if 
          endif
          if (ibeg .ge. iend) then
             if (iadd .ne. 0) then
                getelm = a(iadd)
                return 
             end if 
          end if
          if (ja(imid).gt.j) then
             iend = imid -1
          else 
             ibeg = imid +1
          endif
       end do 
    endif
  end function rgetelm

  function zgetelm(i,j,a,ja,ia,iadd,sorted) result(getelm) 
    integer, intent(in) :: i,j
    integer, intent(in) :: ia(:), ja(:) 
    integer, intent(inout) :: iadd
    complex(dp), intent(in) :: a(:)
    logical, intent(in) :: sorted 

    complex(dp) :: getelm
    integer :: ibeg, iend, imid, k
      
    iadd = 0 
    getelm = (0.0_dp, 0.0_dp)
    ibeg = ia(i)
    iend = ia(i+1)-1
      
    if (.not. sorted) then 
       do k = ibeg, iend
          if (ja(k) .eq.  j) then
             iadd = k 
             if (iadd .ne. 0) then
                getelm = a(iadd)
                return 
             end if 
          endif
       end do
    else
       ! binary search
       do
          imid = ( ibeg + iend ) / 2
          if (ja(imid).eq.j) then
             iadd = imid 
             if (iadd .ne. 0) then
                getelm = a(iadd)
                return 
             end if 
          endif
          if (ibeg .ge. iend) then
             if (iadd .ne. 0) then
                getelm = a(iadd)
                return 
             end if 
          end if
          if (ja(imid).gt.j) then
             iend = imid -1
          else 
             ibeg = imid +1
          endif
       end do 
    endif
  end function zgetelm

  !-----------------------------------------------------------------------
  ! getdia:
  !-----------------------------------------------------------------------
  ! this subroutine extracts a given diagonal from a matrix stored in csr 
  ! format. the output matrix may be transformed with the diagonal removed
  ! from it if desired (as indicated by job.) 
  !----------------------------------------------------------------------- 
  ! our definition of a diagonal of matrix is a vector of length nrow
  ! (always) which contains the elements in rows 1 to nrow of
  ! the matrix that are contained in the diagonal offset by ioff
  ! with respect to the main diagonal. if the diagonal element
  ! falls outside the matrix then it is defined as a zero entry.
  ! thus the proper definition of diag(*) with offset ioff is 
  !
  !     diag(i) = a(i,ioff+i) i=1,2,...,nrow
  !     with elements falling outside the matrix being defined as zero.
  ! 
  !----------------------------------------------------------------------- 
  ! 
  ! on entry:
  !---------- 
  !
  ! nrow  = integer. the row dimension of the matrix a.
  ! ncol  = integer. the column dimension of the matrix a.
  ! job   = integer. job indicator.  if job = 0 then
  !         the matrix a, ja, ia, is not altered on return.
  !         if job.ne.0  then getdia will remove the entries
  !         collected in diag from the original matrix.
  !         this is done in place.
  !
  ! a,ja,
  !    ia = matrix stored in compressed sparse row a,ja,ia,format
  ! ioff  = integer,containing the offset of the wanted diagonal
  !   the diagonal extracted is the one corresponding to the
  !   entries a(i,j) with j-i = ioff.
  !   thus ioff = 0 means the main diagonal
  !
  ! on return:
  !----------- 
  ! len   = number of nonzero elements found in diag.
  !         (len .le. min(nrow,ncol-ioff)-max(1,1-ioff) + 1 )
  !
  ! diag  = double precision array of length nrow containing the wanted diagonal.
  !   diag contains the diagonal (a(i,j),j-i = ioff ) as defined 
  !         above. 
  !
  ! idiag = integer array of  length len, containing the poisitions 
  !         in the original arrays a and ja of the diagonal elements
  !         collected in diag. a zero entry in idiag(i) means that 
  !         there was no entry found in row i belonging to the diagonal.
  !         
  ! a, ja,
  !    ia = if job .ne. 0 the matrix is unchanged. otherwise the nonzero
  !         diagonal entries collected in diag are removed from the 
  !         matrix and therefore the arrays a, ja, ia will change.
  !   (the matrix a, ja, ia will contain len fewer elements) 
  ! 
  !----------------------------------------------------------------------c
  subroutine rgetdia(nrow,ncol,job,a,ja,ia,len,diag,idiag,ioff)
    integer, intent(in) :: nrow, ncol, job   
    real(dp), intent(inout) :: a(:)
    integer, intent(inout) :: ia(:), ja(:)
    integer, intent(in) :: ioff 
    integer, intent(out) :: len
    real(dp), intent(inout) :: diag(:)
    integer, intent(inout) :: idiag(:)
  
    integer :: istart, max, iend, i, kold, k, kdiag, ko
    
    istart = max(0,-ioff)
    iend = min(nrow,ncol-ioff)
    len = 0
    do i = 1, nrow
       idiag(i) = 0
       diag(i) = 0.0d0 
    end do

    ! extract diagonal elements
    do i = istart+1, iend
       do k= ia(i),ia(i+1) -1
          if (ja(k)-i .eq. ioff) then
             diag(i)= a(k)
             idiag(i) = k
             len = len+1
             exit
          endif
       end do 
    end do   
    
    if (job .eq. 0 .or. len .eq.0) return

    !  remove diagonal elements and rewind structure
    ko = 0
    do i = 1, nrow 
       kold = ko
       kdiag = idiag(i) 
       do k= ia(i), ia(i+1)-1 
          if (k .ne. kdiag) then
             ko = ko+1
             a(ko) = a(k)
             ja(ko) = ja(k)
          endif
       end do 
       ia(i) = kold+1
    end do   

    !  redefine ia(nrow+1)
    ia(nrow+1) = ko+1
  end subroutine rgetdia
  !----------------------------------------------------------------------

  subroutine zgetdia(nrow,ncol,job,a,ja,ia,len,diag,idiag,ioff)
    integer, intent(in) :: nrow, ncol, job   
    complex(dp), intent(inout) :: a(:)
    integer, intent(inout) :: ia(:), ja(:)
    integer, intent(in) :: ioff 
    integer, intent(out) :: len
    complex(dp), intent(inout) :: diag(:)
    integer, intent(inout) :: idiag(:)
  
    integer :: istart, max, iend, i, kold, k, kdiag, ko
    
    istart = max(0,-ioff)
    iend = min(nrow,ncol-ioff)
    len = 0
    do i = 1, nrow
       idiag(i) = 0
       diag(i) = 0.0d0 
    end do

    ! extract diagonal elements
    do i = istart+1, iend
       do k= ia(i),ia(i+1) -1
          if (ja(k)-i .eq. ioff) then
             diag(i)= a(k)
             idiag(i) = k
             len = len+1
             exit 
          endif
       end do 
    end do   
    
    if (job .eq. 0 .or. len .eq.0) return

    !  remove diagonal elements and rewind structure
    ko = 0
    do i = 1, nrow 
       kold = ko
       kdiag = idiag(i) 
       do k= ia(i), ia(i+1)-1 
          if (k .ne. kdiag) then
             ko = ko+1
             a(ko) = a(k)
             ja(ko) = ja(k)
          endif
       end do 
       ia(i) = kold+1
    end do   

    !  redefine ia(nrow+1)
    ia(nrow+1) = ko+1
  end subroutine zgetdia
      
  !------------------------------------------------------------------------
  ! In-place transposition routine.
  !------------------------------------------------------------------------
  ! this subroutine transposes a matrix stored in compressed sparse row 
  ! format. the transposition is done in place in that the arrays a,ja,ia
  ! of the transpose are overwritten onto the original arrays.
  !------------------------------------------------------------------------
  ! on entry:
  !--------- 
  ! nrow  = integer. The row dimension of A.
  ! ncol  = integer. The column dimension of A.
  ! a = real array of size nnz (number of nonzero elements in A).
  !         containing the nonzero elements 
  ! ja  = integer array of length nnz containing the column positions
  !     of the corresponding elements in a.
  ! ia  = integer of size n+1, where n = max(nrow,ncol). On entry
  !         ia(k) contains the position in a,ja of  the beginning of 
  !         the k-th row.
  !
  ! iwk = integer work array of same length as ja.
  !
  ! on return:
  !----------
  !
  ! ncol  = actual row dimension of the transpose of the input matrix.
  !         Note that this may be .le. the input value for ncol, in
  !         case some of the last columns of the input matrix are zero
  !         columns. In the case where the actual number of rows found
  !         in transp(A) exceeds the input value of ncol, transp will
  !         return without completing the transposition. see ierr.
  ! a,
  ! ja,
  ! ia  = contains the transposed matrix in compressed sparse
  !         row format. The row dimension of a, ja, ia is now ncol.
  !
  ! ierr  = integer. error message. If the number of rows for the
  !         transposed matrix exceeds the input value of ncol,
  !         then ierr is  set to that number and transp quits.
  !         Otherwise ierr is set to 0 (normal return).
  !
  ! Note: 
  !----- 1) If you do not need the transposition to be done in place
  !         it is preferrable to use the conversion routine csrcsc 
  !         (see conversion routines in formats).
  !      2) the entries of the output matrix are not sorted (the column
  !         indices in each are not in increasing order) use csrcsc
  !         if you want them sorted.
  !----------------------------------------------------------------------
  subroutine rtransp(nrow,ncol,a,ja,ia,iwk,ierr)
    integer, intent(in) :: nrow
    real(dp), intent(inout) :: a(:) 
    integer, intent(inout) :: ia(:), ja(:), iwk(:)
    integer, intent(out) :: ierr, ncol

    real(dp) :: t, t1
    integer :: k, jcol, nnz, i, j, l, init, inext
    logical :: found

    ierr = 0
    nnz = ia(nrow+1)-1
    jcol = 0
    do k = 1, nnz
       jcol = max(jcol,ja(k))
    end do
    if (jcol .gt. ncol) then
       ierr = jcol
       return
    endif
    ! convert to coordinate format. use iwk for row indices
    ncol = jcol
    do i = 1, nrow
       do k = ia(i), ia(i+1)-1
          iwk(k) = i
       end do
    end do   
    !   find pointer array for transpose. 
    do i = 1, ncol+1
       ia(i) = 0
    end do
    do k = 1, nnz
       i = ja(k)
       ia(i+1) = ia(i+1)+1
    end do 
    ia(1) = 1 
    do i = 1, ncol
       ia(i+1) = ia(i) + ia(i+1)
    end do 
     
    !  loop for a cycle in chasing process. 
    init = 1
    k = 0
    mainloop: do 
       t = a(init)
       i = ja(init)
       j = iwk(init)
       iwk(init) = -1
       do 
          k = k+1     
          ! current row number is i.  determine  where to go. 
          l = ia(i)
          !  save the chased element. 
          t1 = a(l)
          inext = ja(l)
          ! then occupy its location.
          a(l)  = t
          ja(l) = j
          ! update pointer information for next element to be put in row i. 
          ia(i) = l+1
          ! determine  next element to be chased
          if (iwk(l) .lt. 0) exit
          
          t = t1
          i = inext
          j = iwk(l)
          iwk(l) = -1
          if (k .ge. nnz) exit mainloop
       end do

       do 
          init = init + 1
          if (init .gt. nnz) exit mainloop
          if (iwk(init) .ge. 0) exit
       end do
       ! restart chasing --  
    end do mainloop

    do i = ncol,1,-1 
       ia(i+1) = ia(i)
    end do
    ia(1) = 1

  end subroutine rtransp

  subroutine ztransp(nrow,ncol,a,ja,ia,iwk,ierr)
    integer, intent(in) :: nrow
    complex(dp), intent(inout) :: a(:) 
    integer, intent(inout) :: ia(:), ja(:), iwk(:)
    integer, intent(out) :: ierr, ncol

    complex(dp) :: t, t1
    integer :: k, jcol, nnz, i, j, l, init, inext
    logical :: found

    ierr = 0
    nnz = ia(nrow+1)-1
    jcol = 0
    do k = 1, nnz
       jcol = max(jcol,ja(k))
    end do
    if (jcol .gt. ncol) then
       ierr = jcol
       return
    endif
    ! convert to coordinate format. use iwk for row indices
    ncol = jcol
    do i = 1, nrow
       do k = ia(i), ia(i+1)-1
          iwk(k) = i
       end do
    end do   
    !   find pointer array for transpose. 
    do i = 1, ncol+1
       ia(i) = 0
    end do
    do k = 1, nnz
       i = ja(k)
       ia(i+1) = ia(i+1)+1
    end do 
    ia(1) = 1 
    do i = 1, ncol
       ia(i+1) = ia(i) + ia(i+1)
    end do 
     
    !  loop for a cycle in chasing process. 
    init = 1
    k = 0
    mainloop: do 
       t = a(init)
       i = ja(init)
       j = iwk(init)
       iwk(init) = -1
       do 
          k = k+1     
          ! current row number is i.  determine  where to go. 
          l = ia(i)
          !  save the chased element. 
          t1 = a(l)
          inext = ja(l)
          ! then occupy its location.
          a(l)  = t
          ja(l) = j
          ! update pointer information for next element to be put in row i. 
          ia(i) = l+1
          ! determine  next element to be chased
          if (iwk(l) .lt. 0) exit
          
          t = t1
          i = inext
          j = iwk(l)
          iwk(l) = -1
          if (k .ge. nnz) exit mainloop
       end do

       do 
          init = init + 1
          if (init .gt. nnz) exit mainloop
          if (iwk(init) .ge. 0) exit
       end do
       ! restart chasing --  
    end do mainloop

    do i = ncol,1,-1 
       ia(i+1) = ia(i)
    end do
    ia(1) = 1

  end subroutine ztransp



  !-----------------------------------------------------------------------
  ! This routine sorts the elements of  a matrix (stored in Compressed
  ! Sparse Row Format) in increasing order of their column indices within 
  ! each row. It uses a form of bucket sort with a cost of O(nnz) where
  ! nnz = number of nonzero elements. 
  ! requires an integer work array of length 2*nnz.  
  !-----------------------------------------------------------------------
  ! on entry:
  !--------- 
  ! n     = the row dimension of the matrix
  ! a     = the matrix A in compressed sparse row format.
  ! ja    = the array of column indices of the elements in array a.
  ! ia    = the array of pointers to the rows. 
  ! iwork = integer work array of length max ( n+1, 2*nnz ) 
  !         where nnz = (ia(n+1)-ia(1))  ) .
  ! values= logical indicating whether or not the real values a(*) must 
  !         also be permuted. if (.not. values) then the array a is not
  !         touched by csort and can be a dummy array. 
  ! 
  ! on return:
  !----------
  ! the matrix stored in the structure a, ja, ia is permuted in such a
  ! way that the column indices are in increasing order within each row.
  ! iwork(1:nnz) contains the permutation used  to rearrange the elements.
  !----------------------------------------------------------------------- 
  subroutine rcsort(nrow,a,ja,ia,iwork,values)
    integer, intent(in) :: nrow
    real(dp), intent(inout) :: a(:) 
    integer, intent(inout), dimension(:) :: ja, ia, iwork 
    logical, intent(in) :: values

    integer i, k, j, ifirst, nnz, next, ko, irow  
    
    iwork = 0
    do i = 1, nrow
      do k = ia(i), ia(i+1)-1 
         j = ja(k)+1
         iwork(j) = iwork(j)+1
      end do
    end do
    ! compute pointers from lengths. 
    iwork(1) = 1
    do i = 1, nrow
       iwork(i+1) = iwork(i) + iwork(i+1)
    end do

    ! get the positions of the nonzero elements in order of columns.
    ifirst = ia(1) 
    nnz = ia(nrow+1)-ifirst
    do i = 1, nrow
       do k = ia(i),ia(i+1)-1 
          j = ja(k) 
          next = iwork(j) 
          iwork(nnz+next) = k
          iwork(j) = next+1
       end do
    end do
    ! convert to coordinate format
    do i=1, nrow
       do k=ia(i), ia(i+1)-1 
          iwork(k) = i
       end do
    end do

    ! loop to find permutation: for each element find the correct 
    ! position in (sorted) arrays a, ja. Record this in iwork. 
 
    do k = 1, nnz
       ko = iwork(nnz+k) 
       irow = iwork(ko)
       next = ia(irow)

       ! the current element should go in next position in row. iwork
       ! records this position. 
       iwork(ko) = next
       ia(irow)  = next+1
    end do

    ! perform an in-place permutation of the  arrays.

    call vperm(nnz, ja, iwork)                  !Commento: era inizialmente call vperm(nnz, ja(ifirst), iwork)
    if (values) call vperm(nnz, a, iwork)       !Commento: era inizialmente call vperm(nnz, a(ifirst), iwork)

    ! reshift the pointers of the original matrix back.
    do i = nrow, 1, -1
       ia(i+1) = ia(i)
    end do
    ia(1) = ifirst 
  end subroutine rcsort
  ! -------- Cmplx version -------------------------------------------------
  subroutine zcsort(nrow,a,ja,ia,iwork,values) 
    integer, intent(in) :: nrow
    complex(dp), intent(inout) :: a(:) 
    integer, intent(inout) :: ja(:), ia(:), iwork(:) 
    logical, intent(in) :: values

    integer i, k, j, ifirst, nnz, next, ko, irow  
      
    iwork = 0
    do i = 1, nrow
      do k = ia(i), ia(i+1)-1 
         j = ja(k)+1
         iwork(j) = iwork(j)+1
      end do
    end do
    ! compute pointers from lengths. 
    iwork(1) = 1
    do i = 1, nrow
       iwork(i+1) = iwork(i) + iwork(i+1)
    end do

    ! get the positions of the nonzero elements in order of columns.
    ifirst = ia(1) 
    nnz = ia(nrow+1)-ifirst
    do i = 1, nrow
       do k = ia(i),ia(i+1)-1 
          j = ja(k) 
          next = iwork(j) 
          iwork(nnz+next) = k
          iwork(j) = next+1
       end do
    end do
    ! convert to coordinate format
    do i=1, nrow
       do k=ia(i), ia(i+1)-1 
          iwork(k) = i
       end do
    end do

    ! loop to find permutation: for each element find the correct 
    ! position in (sorted) arrays a, ja. Record this in iwork. 
 
    do k = 1, nnz
       ko = iwork(nnz+k) 
       irow = iwork(ko)
       next = ia(irow)

       ! the current element should go in next position in row. iwork
       ! records this position. 
       iwork(ko) = next
       ia(irow)  = next+1
    end do

    ! perform an in-place permutation of the  arrays.

    call vperm(nnz, ja, iwork) 
    if (values) call vperm(nnz, a, iwork) 

    ! reshift the pointers of the original matrix back.
    do i = nrow, 1, -1
       ia(i+1) = ia(i)
    end do
    ia(1) = ifirst 

  end subroutine zcsort

  !-----------------------------------------------------------------------
  ! this subroutine performs an in-place permutation of an integer vector 
  ! ix according to the permutation array perm(*), i.e., on return, 
  ! the vector x satisfies,
  !
  ! ix(perm(j)) :== ix(j), j=1,2,.., n
  !
  !-----------------------------------------------------------------------
  ! on entry:
  !---------
  ! n   = length of vector x.
  ! perm  = integer array of length n containing the permutation  array.
  ! ix  = input vector
  !
  ! on return:
  !---------- 
  ! ix  = vector x permuted according to ix(perm(*)) :=  ix(*)
  !----------------------------------------------------------------------c
  subroutine ivperm(n, ix, perm) 
    integer, intent(in) :: n
    integer, intent(inout) :: ix(n)
    integer, intent(inout) :: perm(n)
      
    integer :: ii, k, init
    integer :: tmp, tmp1, next

    init      = 1
    tmp = ix(init)  
    ii  = perm(init)
    perm(init)= -perm(init)
    k = 0

    mainloop: do
       do
          k = k+1
          tmp1    = ix(ii) 
          ix(ii)     = tmp
          next    = perm(ii) 
          if (next .lt. 0 ) exit
          if (k .gt. n) exit mainloop
          tmp       = tmp1
          perm(ii)  = - perm(ii)
          ii        = next 
       end do

       ! reinitilaize cycle --
       do 
          init = init+1
          if (init .gt. n) exit mainloop
          if (perm(init) .ge. 0) exit
       end do
       tmp = ix(init)
       ii  = perm(init)
       perm(init)=-perm(init)
    end do mainloop 
    
    do k = 1, n
       perm(k) = -perm(k)
    end do 
 
  end subroutine ivperm
  ! ---- Real version ----------------------------------------------------
  subroutine rvperm(n, x, perm) 
    integer, intent(in) :: n
    real(dp), intent(inout) :: x(n)
    integer, intent(inout) :: perm(n)
      
    integer :: ii, k, init, next
    real(dp) :: tmp, tmp1

    init      = 1
    tmp = x(init)  
    ii  = perm(init)
    perm(init)= -perm(init)
    k = 0

    mainloop: do
       do
          k = k+1
          tmp1    = x(ii) 
          x(ii)     = tmp
          next    = perm(ii) 
          if (next .lt. 0 ) exit
          if (k .gt. n) exit mainloop
          tmp       = tmp1
          perm(ii)  = - perm(ii)
          ii        = next 
       end do

       ! reinitilaize cycle --
       do 
          init = init+1
          if (init .gt. n) exit mainloop
          if (perm(init) .ge. 0) exit
       end do
       tmp = x(init)
       ii  = perm(init)
       perm(init)=-perm(init)
    end do mainloop 
    
    do k = 1, n
       perm(k) = -perm(k)
    end do 
  end subroutine rvperm
  ! ---- Cmplx version ----------------------------------------------------
  subroutine zvperm(n, x, perm) 
    integer, intent(in) :: n
    complex(dp), intent(inout) :: x(n)
    integer, intent(inout) :: perm(n)
      
    integer :: ii, k, init, next
    complex(dp) :: tmp, tmp1

    init      = 1
    tmp = x(init)  
    ii  = perm(init)
    perm(init)= -perm(init)
    k = 0

    mainloop: do
       do
          k = k+1
          tmp1    = x(ii) 
          x(ii)     = tmp
          next    = perm(ii) 
          if (next .lt. 0 ) exit
          if (k .gt. n) exit mainloop
          tmp       = tmp1
          perm(ii)  = - perm(ii)
          ii        = next 
       end do

       ! reinitilaize cycle --
       do 
          init = init+1
          if (init .gt. n) exit mainloop
          if (perm(init) .ge. 0) exit
       end do
       tmp = x(init)
       ii  = perm(init)
       perm(init)=-perm(init)
    end do mainloop 
    
    do k = 1, n
       perm(k) = -perm(k)
    end do 
  end subroutine zvperm
    
    
  !-----------------------------------------------------------------------
  ! this routine computes the lower, upper, maximum, and average 
  ! bandwidths.     revised -- July 12, 2001  -- bug fix -- YS. 
  !-----------------------------------------------------------------------
  ! On Entry:
  !----------
  ! n     = integer. column dimension of matrix
  ! a     = real array containing the nonzero elements of the matrix
  !         the elements are stored by columns in order
  !         (i.e. column i comes before column i+1, but the elements
  !         within each column can be disordered).
  ! ja    = integer array containing the row indices of elements in a
  ! ia    = integer array containing of length n+1 containing the
  !         pointers to the beginning of the columns in arrays a and ja.
  !         It is assumed that ia(*) = 1 and ia(n+1) = nzz+1.
  !
  ! on return
  !----------
  ! ml    = lower bandwidth as defined by
  !        ml = max(i-j | all  a(i,j).ne. 0)
  ! mu    = upper bandwidth as defined by
  !        mu = max ( j-i | all  a(i,j).ne. 0 )
  ! iband =  maximum bandwidth as defined by
  !         iband = Max (  Max [ j | a(i,j) .ne. 0 ] - 
  !                        Min [ j | a(i,j) .ne. 0 ] )
  ! bndav = Average Bandwidth          
  !-----------------------------------------------------------------------
  subroutine bandwidth(n, ja, ia, ml, mu, iband, bndav)
    integer, intent(in) :: n    
    integer, intent(in) :: ja(:), ia(:)
    integer, intent(out) :: ml, mu, iband
    real(dp), intent(out) :: bndav

    integer :: max 
    integer :: j0, j1, jminc, jmaxc, i
    
    ml = -n
    mu = -n
    bndav = 0.0_dp
    iband = 0 
    do i = 1, n
      j0 = ia(i)
      j1 = ia(i+1) - 1
      jminc = ja(j0)
      jmaxc = ja(j1)
      ml = max(ml,i-jminc)
      mu = max(mu,jmaxc-i)
      iband = max(iband,jmaxc-jminc+1)
      bndav = bndav+real( jmaxc-jminc+1, dp)
    end do
    bndav = bndav/real(n,dp)

  end subroutine bandwidth


!-----------------------------------------------------------------------
! This subroutine builds a sparse matrix from an input matrix by 
! extracting only elements in positions defined by the mask jmask, imask
!-----------------------------------------------------------------------
! On entry:
!---------
! nrow  = integer. row dimension of input matrix 
! ncol  = integer. Column dimension of input matrix.
!
! a,
! ja,
! ia  = matrix in Compressed Sparse Row format
!
! jmask,
! imask = matrix defining mask (pattern only) stored in compressed
!         sparse row format.
!
! nzmax = length of arrays c and jc. see ierr.
! 
! On return:
!-----------
!
! a, ja, ia and jmask, imask are unchanged.
!
! c
! jc, 
! ic  = the output matrix in Compressed Sparse Row format.
! 
! ierr  = integer. serving as error message.c
!         ierr = 1  means normal return
!         ierr .gt. 1 means that amask stopped when processing
!         row number ierr, because there was not enough space in
!         c, jc according to the value of nzmax.
!
! work arrays:
!------------- 
! iw  = logical work array of length ncol.
!
! note: 
!------ the  algorithm is in place: c, jc, ic can be the same as 
! a, ja, ia in which cas the code will overwrite the matrix c
! on a, ja, ia
!
!-----------------------------------------------------------------------

  subroutine ramask(nrow,ncol,a,ja,ia,jmask,imask,c,jc,ic,iw,nzmax,ierr)
      real(dp) :: a(:),c(:) 
      integer, intent(in) :: nrow, ncol, nzmax
      integer, intent(in) :: ia(:),ja(:), jmask(:)
      integer, intent(inout) :: jc(:), ic(:), imask(:) 
      integer, intent(out) :: ierr
      logical :: iw(:)
 
      integer :: ii, j, k, len, k1, k2

      ierr = 0
      len = 0

      iw = .false.

!     unpack the mask for row ii in iw
      do ii=1, nrow
!     save pointer in order to be able to do things in place
         do k=imask(ii), imask(ii+1)-1
            iw(jmask(k)) = .true.
         end do
!     add umasked elemnts of row ii
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         ic(ii) = len+1
         do k = k1,k2 
            j = ja(k)
            if (iw(j)) then
               len = len+1
               if (len .gt. nzmax) then
                  ierr = ii
                  return
               endif
               jc(len) = j
               c(len) = a(k)
            endif
         end do       
     
         do k = imask(ii), imask(ii+1)-1
            iw(jmask(k)) = .false.
         end do
      end do 
      ic(nrow+1)=len+1

   end subroutine ramask

  subroutine zamask(nrow,ncol,a,ja,ia,jmask,imask,c,jc,ic,iw,nzmax,ierr)
      complex(dp) :: a(:),c(:) 
      integer, intent(in) :: nrow, ncol, nzmax
      integer, intent(in) :: ia(:),ja(:), jmask(:)
      integer, intent(inout) :: jc(:), ic(:), imask(:) 
      integer, intent(out) :: ierr
      logical :: iw(:)
 
      integer :: ii, j, k, len, k1, k2

      ierr = 0
      len = 0

      iw = .false.

!     unpack the mask for row ii in iw
      do ii=1, nrow
!     save pointer in order to be able to do things in place
         do k=imask(ii), imask(ii+1)-1
            iw(jmask(k)) = .true.
         end do
!     add umasked elemnts of row ii
         k1 = ia(ii)
         k2 = ia(ii+1)-1
         ic(ii) = len+1
         do k = k1,k2 
            j = ja(k)
            if (iw(j)) then
               len = len+1
               if (len .gt. nzmax) then
                  ierr = ii
                  return
               endif
               jc(len) = j
               c(len) = a(k)
            endif
         end do       
     
         do k = imask(ii), imask(ii+1)-1
            iw(jmask(k)) = .false.
         end do
      end do 
      ic(nrow+1)=len+1

   end subroutine zamask

end module skit_unary

! UNARY OPERATIONS
!       
!    subroutine zcsort(nrow, aa, ja, ia, iw, vals)
!
!    subroutine ztransp(nrow,ncol,aa,ja,ia,iwk,ierr)
!    

