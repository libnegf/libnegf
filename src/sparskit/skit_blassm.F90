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

module skit_blassm
  use ln_precision
  implicit none
  private

  public :: amub, amubs, aplb, aplb1, cplsamub
  public :: as1pls2b, copymat

  interface amub
    module procedure :: ramub    
    module procedure :: zamub    
  end interface amub
  
  interface amubs
    !module procedure :: ramub    
    module procedure :: zamubs    
  end interface amubs

  interface aplb
    module procedure :: raplb    
    module procedure :: zaplb    
  end interface aplb 

  interface aplb1
    module procedure :: raplb1    
    module procedure :: zaplb1    
  end interface aplb1 

  interface cplsamub
    !module procedure :: rcplsamub   
    module procedure :: zcplsamub   
  end interface cplsamub 

  interface as1pls2b
    !module procedure :: ras1pls2b     
    module procedure :: zas1pls2b 
  end interface as1pls2b

  interface copymat 
    !module procedure :: rcopymat     
    module procedure :: zcopymat 
  end interface copymat

  contains
  !-----------------------------------------------------------------------
  ! performs the matrix by matrix product C = A B 
  !-----------------------------------------------------------------------
  ! on entry:
  ! ---------
  ! nrow  = integer. The row dimension of A = row dimension of C
  ! ncol  = integer. The column dimension of B = column dimension of C
  ! job   = integer. Job indicator. When job = 0, only the structure
  !                  (i.e. the arrays jc, ic) is computed and the
  !                  real values are ignored.
  !
  ! a,
  ! ja,
  ! ia   = Matrix A in compressed sparse row format.
  ! 
  ! b, 
  ! jb, 
  ! ib    =  Matrix B in compressed sparse row format.
  !
  ! nzmax = integer. The  length of the arrays c and jc.
  !         amub will stop if the result matrix C  has a number 
  !         of elements that exceeds exceeds nzmax. See ierr.
  ! 
  ! on return:
  !----------
  ! c, 
  ! jc, 
  ! ic    = resulting matrix C in compressed sparse row sparse format.
  !           
  ! ierr  = integer. serving as error message. 
  !         ierr = 0 means normal return,
  !         ierr .gt. 0 means that amub stopped while computing the
  !         i-th row  of C with i=ierr, because the number 
  !         of elements in C exceeds nzmax.
  !
  ! work arrays:
  !------------
  ! iw    = integer work array of length equal to the number of
  !         columns in A.
  ! Note: 
  !-------
  !   The row dimension of B is not needed. However there is no checking 
  !   on the condition that ncol(A) = nrow(B). 
  !
  !----------------------------------------------------------------------- 
  subroutine ramub(nrow,ncol,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr) 
    integer, intent(in) :: nrow, ncol, nzmax, job 
    real(dp), intent(in) :: a(:), b(:)
    integer, intent(in) :: ja(:), jb(:), ia(:), ib(:), iw(:)
    real(dp), intent(inout) :: c(:)
    integer, intent(inout) :: jc(:), ic(:)
    integer, intent(out) :: ierr 

    integer :: len, ii, jj, kk, ka, kb, jpos, jcol
    real(dp) :: scal 
    logical :: values

    values = (job .ne. 0) 
    len = 0
    ic(1) = 1 
    ierr = 0
    iw = 0

    do ii = 1, nrow 
      do ka = ia(ii), ia(ii+1)-1 
        if (values) then
           scal = a(ka)
        end if   
        jj = ja(ka)
        do kb = ib(jj), ib(jj+1)-1
           jcol = jb(kb)
           jpos = iw(jcol)
           if (jpos .eq. 0) then
              len = len + 1
              if (len .gt. nzmax) then
                 ierr = ii
                 return
              endif
              jc(len) = jcol
              iw(jcol)= len
              if (values) then
                 c(len)  = scal*b(kb)
              endif
           else
              if (values) then
                 c(jpos) = c(jpos) + scal*b(kb)
              endif
           endif
         end do
       end do
       do kk = ic(ii), len
          iw(jc(kk)) = 0
       end do
       ic(ii+1) = len+1
     end do
  end subroutine ramub
  
  subroutine zamub(nrow,ncol,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr) 
    integer, intent(in) :: nrow, ncol, nzmax, job 
    complex(dp), intent(in) :: a(:), b(:)
    integer, intent(in) :: ja(:), jb(:), ia(:), ib(:), iw(:)
    complex(dp), intent(inout) :: c(:)
    integer, intent(inout) :: jc(:), ic(:)
    integer, intent(out) :: ierr 

    integer :: len, ii, jj, kk, ka, kb, jpos, jcol
    complex(dp) :: scal 
    logical :: values

    values = (job .ne. 0) 
    len = 0
    ic(1) = 1 
    ierr = 0
    iw = 0

    do ii = 1, nrow 
      do ka = ia(ii), ia(ii+1)-1 
        if (values) then
           scal = a(ka)
        end if   
        jj = ja(ka)
        do kb = ib(jj), ib(jj+1)-1
           jcol = jb(kb)
           jpos = iw(jcol)
           if (jpos .eq. 0) then
              len = len + 1
              if (len .gt. nzmax) then
                 ierr = ii
                 return
              endif
              jc(len) = jcol
              iw(jcol)= len
              if (values) then
                 c(len)  = scal*b(kb)
              endif
           else
              if (values) then
                 c(jpos) = c(jpos) + scal*b(kb)
              endif
           endif
        end do
      end do
      do kk = ic(ii), len
         iw(jc(kk)) = 0
      end do
      ic(ii+1) = len+1
    end do
  end subroutine zamub
       
  !-----------------------------------------------------------------------
  ! performs the matrix by matrix product C = s*A B 
  !-----------------------------------------------------------------------
  ! All the rest is like zamub
  subroutine zamubs(nrow,ncol,job,a,ja,ia,s,b,jb,ib,c,jc,ic,nzmax,iw,ierr)            !Commento: è identica a zamub
    integer, intent(in) :: nrow, ncol, nzmax, job 
    complex(dp), intent(in) :: a(:), b(:), s
    integer, intent(in) :: ja(:), ia(:), jb(:), ib(:), iw(:)
    complex(dp), intent(inout) :: c(:)
    integer, intent(inout) :: jc(:), ic(:)
    integer, intent(out) :: ierr 

    integer :: len, ii, jj, kk, ka, kb, jpos, jcol
    complex(dp) :: scal 
    logical :: values

    values = (job .ne. 0) 
    len = 0
    ic(1) = 1 
    ierr = 0
    iw = 0

    do ii = 1, nrow 
      do ka = ia(ii), ia(ii+1)-1 
        if (values) then
           scal = a(ka)*s
        end if   
        jj = ja(ka)
        do kb = ib(jj), ib(jj+1)-1
           jcol = jb(kb)
           jpos = iw(jcol)
           if (jpos .eq. 0) then
              len = len + 1
              if (len .gt. nzmax) then
                 ierr = ii
                 return
              endif
              jc(len) = jcol
              iw(jcol)= len
              if (values) then
                 c(len)  = scal*b(kb)
              endif
           else
              if (values) then
                 c(jpos) = c(jpos) + scal*b(kb)
              endif
           endif
        end do
      end do
      do kk = ic(ii), len
         iw(jc(kk)) = 0
      end do
      ic(ii+1) = len+1
    end do
  end subroutine zamubs

  !-----------------------------------------------------------------------
  ! performs the matrix sum  C = A+B. 
  !-----------------------------------------------------------------------
  ! on entry:
  ! ---------
  ! nrow  = integer. The row dimension of A and B
  ! ncol  = integer. The column dimension of A and B.
  ! job   = integer. Job indicator. When job = 0, only the structure
  !                  (i.e. the arrays jc, ic) is computed and the
  !                  real values are ignored.
  !
  ! a,
  ! ja,
  ! ia   = Matrix A in compressed sparse row format.
  ! 
  ! b, 
  ! jb, 
  ! ib  =  Matrix B in compressed sparse row format.
  !
  ! nzmax = integer. The  length of the arrays c and jc.
  !         amub will stop if the result matrix C  has a number 
  !         of elements that exceeds exceeds nzmax. See ierr.
  ! 
  ! on return:
  !----------
  ! c, 
  ! jc, 
  ! ic  = resulting matrix C in compressed sparse row sparse format.
  !     
  ! ierr  = integer. serving as error message. 
  !         ierr = 0 means normal return,
  !         ierr .gt. 0 means that amub stopped while computing the
  !         i-th row  of C with i=ierr, because the number 
  !         of elements in C exceeds nzmax.
  !
  ! work arrays:
  !------------
  ! iw  = integer work array of length equal to the number of
  !         columns in A.
  !-----------------------------------------------------------------------
  subroutine raplb(nrow,ncol,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
    integer, intent(in) :: nrow, ncol, nzmax, job 
    real(dp), intent(in) :: a(:), b(:)
    integer, intent(in) :: ja(:), ia(:), jb(:), ib(:), iw(:)
    real(dp), intent(inout) :: c(:)
    integer, intent(inout) :: jc(:), ic(:)
    integer, intent(out) :: ierr 
    
    integer :: ii, jj, kk, ka, kb, jcol, jpos   
    logical :: values
    
    values = (job .ne. 0) 
    ierr = 0
    len = 0
    ic(1) = 1 
    iw = 0
      
    do ii=1, nrow
       do ka=ia(ii), ia(ii+1)-1 
          len = len+1
          jcol    = ja(ka)
          if (len .gt. nzmax) then
             ierr = ii
             return   
          end if      
          jc(len) = jcol 
          if (values) c(len)  = a(ka) 
          iw(jcol)= len
       end do    
         
       do kb=ib(ii),ib(ii+1)-1
          jcol = jb(kb)
          jpos = iw(jcol)
          if (jpos .eq. 0) then
             len = len+1
             if (len .gt. nzmax) then
                ierr = ii
                return   
             end if      
             jc(len) = jcol
             if (values) c(len)  = b(kb)
             iw(jcol)= len
          else
             if (values) c(jpos) = c(jpos) + b(kb)
          endif
       end do
       do kk=ic(ii), len
          iw(jc(kk)) = 0
       end do 
       ic(ii+1) = len+1
    end do
  end subroutine raplb
  
  subroutine zaplb(nrow,ncol,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
    integer, intent(in) :: nrow, ncol, nzmax, job 
    complex(dp), intent(in) :: a(:), b(:)
    integer, intent(in) :: ja(:), ia(:), jb(:), ib(:), iw(:)
    complex(dp), intent(inout) :: c(:)
    integer, intent(inout) :: jc(:), ic(:)
    integer, intent(out) :: ierr 
    
    integer :: ii, jj, kk, ka, kb, jcol, jpos   
    logical :: values
    
    values = (job .ne. 0) 
    ierr = 0
    len = 0
    ic(1) = 1 
    iw = 0
      
    do ii=1, nrow
       do ka=ia(ii), ia(ii+1)-1 
          len = len+1
          jcol    = ja(ka)
          if (len .gt. nzmax) then
             ierr = ii
             return   
          end if      
          jc(len) = jcol 
          if (values) c(len)  = a(ka) 
          iw(jcol)= len
       end do    
         
       do kb=ib(ii),ib(ii+1)-1
          jcol = jb(kb)
          jpos = iw(jcol)
          if (jpos .eq. 0) then
             len = len+1
             if (len .gt. nzmax) then
                ierr = ii
                return   
             end if      
             jc(len) = jcol
             if (values) c(len)  = b(kb)
             iw(jcol)= len
          else
             if (values) c(jpos) = c(jpos) + b(kb)
          endif
       end do
       do kk=ic(ii), len
          iw(jc(kk)) = 0
       end do 
       ic(ii+1) = len+1
    end do
  end subroutine zaplb
 
  !-----------------------------------------------------------------------
  ! performs the matrix sum  C = A+B for matrices in sorted CSR format.
  ! the difference with aplb  is that the resulting matrix is such that
  ! the elements of each row are sorted with increasing column indices in
  ! each row, provided the original matrices are sorted in the same way. 
  !-----------------------------------------------------------------------
  ! on entry:
  ! ---------
  ! nrow  = integer. The row dimension of A and B
  ! ncol  = integer. The column dimension of A and B.
  ! job   = integer. Job indicator. When job = 0, only the structure
  !                  (i.e. the arrays jc, ic) is computed and the
  !                  real values are ignored.
  !
  ! a,
  ! ja,
  ! ia   = Matrix A in compressed sparse row format with entries sorted
  ! 
  ! b, 
  ! jb, 
  ! ib  =  Matrix B in compressed sparse row format with entries sorted
  !        ascendly in each row   
  !
  ! nzmax = integer. The  length of the arrays c and jc.
  !         amub will stop if the result matrix C  has a number 
  !         of elements that exceeds exceeds nzmax. See ierr.
  ! 
  ! on return:
  !----------
  ! c, 
  ! jc, 
  ! ic  = resulting matrix C in compressed sparse row sparse format
  !         with entries sorted ascendly in each row. 
  !     
  ! ierr  = integer. serving as error message. 
  !         ierr = 0 means normal return,
  !         ierr .gt. 0 means that amub stopped while computing the
  !         i-th row  of C with i=ierr, because the number 
  !         of elements in C exceeds nzmax.
  !
  ! Notes: 
  !-------
  !     this will not work if any of the two input matrices is not sorted
  !-----------------------------------------------------------------------
  subroutine raplb1(nrow,ncol,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,ierr)
    integer, intent(in) :: nrow, ncol, nzmax, job 
    real(dp), intent(in) :: a(:), b(:)
    integer, intent(in) :: ja(:), ia(:), jb(:), ib(:), iw(:)
    real(dp), intent(inout) :: c(:)
    integer, intent(inout) :: jc(:), ic(:)
    integer, intent(out) :: ierr 
  
    integer :: ii, j1, j2, ka, kb, kc, kamax, kbmax    
    logical :: values
    
    values = (job .ne. 0) 
    ierr = 0
    kc = 1
    ic(1) = kc 

    do ii=1, nrow
       ka = ia(ii)
       kb = ib(ii)
       kamax = ia(ii+1)-1
       kbmax = ib(ii+1)-1 
       do while (ka .le. kamax .or. kb .le. kbmax)
          if (ka .le. kamax) then
             j1 = ja(ka)
          else
             j1 = ncol+1
          endif
          if (kb .le. kbmax) then 
             j2 = jb(kb)         
          else 
             j2 = ncol+1
          endif
       
          if (j1 .eq. j2) then 
             jc(kc) = j1
             if (values) c(kc) = a(ka)+b(kb)
             ka = ka+1
             kb = kb+1
             kc = kc+1
          else if (j1 .lt. j2) then
             jc(kc) = j1
             if (values) c(kc) = a(ka)
             ka = ka+1
             kc = kc+1
          else if (j1 .gt. j2) then
             jc(kc) = j2
             if (values) c(kc) = b(kb)
             kb = kb+1
             kc = kc+1
          endif
          if (kc .gt. nzmax) then
             ierr = ii
             return
          end if   
       end do
       ic(i+1) = kc
     end do 
  end subroutine raplb1
  
  subroutine zaplb1(nrow,ncol,job,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,ierr)
    integer, intent(in) :: nrow, ncol, nzmax, job 
    complex(dp), intent(in) :: a(:), b(:)
    integer, intent(in) :: ja(:), ia(:), jb(:), ib(:), iw(:)
    complex(dp), intent(inout) :: c(:)
    integer, intent(inout) :: jc(:), ic(:)
    integer, intent(out) :: ierr 
  
    integer :: ii, j1, j2, ka, kb, kc, kamax, kbmax    
    logical :: values
    
    values = (job .ne. 0) 
    ierr = 0
    kc = 1
    ic(1) = kc 

    do ii=1, nrow
       ka = ia(ii)
       kb = ib(ii)
       kamax = ia(ii+1)-1
       kbmax = ib(ii+1)-1 
       do while (ka .le. kamax .or. kb .le. kbmax)
          if (ka .le. kamax) then
             j1 = ja(ka)
          else
             j1 = ncol+1
          endif
          if (kb .le. kbmax) then 
             j2 = jb(kb)         
          else 
             j2 = ncol+1
          endif
       
          if (j1 .eq. j2) then 
             jc(kc) = j1
             if (values) c(kc) = a(ka)+b(kb)
             ka = ka+1
             kb = kb+1
             kc = kc+1
          else if (j1 .lt. j2) then
             jc(kc) = j1
             if (values) c(kc) = a(ka)
             ka = ka+1
             kc = kc+1
          else if (j1 .gt. j2) then
             jc(kc) = j2
             if (values) c(kc) = b(kb)
             kb = kb+1
             kc = kc+1
          endif
          if (kc .gt. nzmax) then
             ierr = ii
             return
          end if   
       end do
       ic(i+1) = kc
    end do 
  end subroutine zaplb1
      
  !-----------------------------------------------------------------------
  ! performs the matrix sum  C = A+sB. 
  ! Nota: inserita il 22/02/2006 in correzione all'originale
  !-----------------------------------------------------------------------
  ! on entry:
  ! ---------
  ! nrow  = integer. The row dimension of A and B
  ! ncol  = integer. The column dimension of A and B.
  ! job   = integer. Job indicator. When job = 0, only the structure
  !                  (i.e. the arrays jc, ic) is computed and the
  !                  real values are ignored.
  !
  ! a,
  ! ja,
  ! ia   = Matrix A in compressed sparse row format.
  !
  ! s    = Scalar factor for B
  ! 
  ! b, 
  ! jb, 
  ! ib  =  Matrix B in compressed sparse row format.
  !
  ! nzmax = integer. The  length of the arrays c and jc.
  !         amub will stop if the result matrix C  has a number 
  !         of elements that exceeds exceeds nzmax. See ierr.
  ! 
  ! on return:
  !----------
  ! c, 
  ! jc, 
  ! ic  = resulting matrix C in compressed sparse row sparse format.
  !     
  ! ierr  = integer. serving as error message. 
  !         ierr = 0 means normal return,
  !         ierr .gt. 0 means that amub stopped while computing the
  !         i-th row  of C with i=ierr, because the number 
  !         of elements in C exceeds nzmax.
  !
  ! work arrays:
  !------------
  ! iw  = integer work array of length equal to the number of
  !         columns in A.
  !
  !-----------------------------------------------------------------------
  subroutine raplsb(nrow,ncol,job,a,ja,ia,s,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
    integer, intent(in) :: nrow, ncol, nzmax, job 
    real(dp), intent(in) :: a(:), b(:), s
    integer, intent(in) :: ja(:), ia(:), jb(:), ib(:), iw(:)
    real(dp), intent(inout) :: c(:)
    integer, intent(inout) :: jc(:), ic(:)
    integer, intent(out) :: ierr 

    integer :: ii, kk, ka, kb, len, jcol, jpos
    logical :: values
     
    values = (job .ne. 0) 
    ierr = 0
    len = 0
    ic(1) = 1 
    iw = 0

    do ii = 1, nrow
       do ka = ia(ii), ia(ii+1)-1 
          len = len+1
          jcol    = ja(ka)
          if (len .gt. nzmax) then
             ierr = ii
             return
          end if   
          jc(len) = jcol 
          if (values) c(len)  = a(ka) 
          iw(jcol)= len
       end do

       do kb = ib(ii), ib(ii+1)-1
          jcol = jb(kb)
          jpos = iw(jcol)
          if (jpos .eq. 0) then
             len = len+1
             if (len .gt. nzmax) then
                ierr = ii
                return
             end if      
             jc(len) = jcol
             if (values) c(len)  = s*b(kb)
             iw(jcol)= len
          else
             if (values) c(jpos) = c(jpos) + s*b(kb)
          endif
       end do
       do kk = ic(ii), len
          iw(jc(kk)) = 0
       end do
       ic(ii+1) = len+1
    end do 
  end subroutine raplsb
  
  subroutine zaplsb(nrow,ncol,job,a,ja,ia,s,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
    integer, intent(in) :: nrow, ncol, nzmax, job 
    complex(dp), intent(in) :: a(:), b(:), s
    integer, intent(in) :: ja(:), ia(:), jb(:), ib(:), iw(:)
    complex(dp), intent(inout) :: c(:)
    integer, intent(inout) :: jc(:), ic(:)
    integer, intent(out) :: ierr 

    integer :: ii, kk, ka, kb, len, jcol, jpos
    logical :: values
     
    values = (job .ne. 0) 
    ierr = 0
    len = 0
    ic(1) = 1 
    iw = 0

    do ii = 1, nrow
       do ka = ia(ii), ia(ii+1)-1 
          len = len+1
          jcol    = ja(ka)
          if (len .gt. nzmax) then
             ierr = ii
             return
          end if   
          jc(len) = jcol 
          if (values) c(len)  = a(ka) 
          iw(jcol)= len
       end do

       do kb = ib(ii), ib(ii+1)-1
          jcol = jb(kb)
          jpos = iw(jcol)
          if (jpos .eq. 0) then
             len = len+1
             if (len .gt. nzmax) then
                ierr = ii
                return
             end if      
             jc(len) = jcol
             if (values) c(len)  = s*b(kb)
             iw(jcol)= len
          else
             if (values) c(jpos) = c(jpos) + s*b(kb)
          endif
       end do
       do kk = ic(ii), len
          iw(jc(kk)) = 0
       end do
       ic(ii+1) = len+1
    end do 
  end subroutine zaplsb

  !-----------------------------------------------------------------------
  ! performs the matrix by matrix product C = C + s A B 
  !-----------------------------------------------------------------------
  ! on entry:
  ! ---------
  ! nrow  = integer. The row dimension of A = row dimension of C
  ! ncol  = integer. The column dimension of B = column dimension of C
  ! job   = integer. Job indicator. When job = 0, only the structure
  !                  (i.e. the arrays jc, ic) is computed and the
  !                  real values are ignored.
  !
  ! a,
  ! ja,
  ! ia   = Matrix A in compressed sparse row format.
  ! 
  ! b, 
  ! jb, 
  ! ib    =  Matrix B in compressed sparse row format.
  !
  ! nzmax = integer. The  length of the arrays c and jc.
  !         amub will stop if the result matrix C  has a number 
  !         of elements that exceeds exceeds nzmax. See ierr.
  ! 
  ! on return:
  !----------
  ! c, 
  ! jc, 
  ! ic    = resulting matrix C in compressed sparse row sparse format.
  !           
  ! ierr  = integer. serving as error message. 
  !         ierr = 0 means normal return,
  !         ierr .gt. 0 means that amub stopped while computing the
  !         i-th row  of C with i=ierr, because the number 
  !         of elements in C exceeds nzmax.
  !
  ! work arrays:
  !------------
  ! iw    = integer work array of length equal to the number of
  !         columns in A.
  ! Note: 
  !-------
  !   The row dimension of B is not needed. However there is no checking 
  !   on the condition that ncol(A) = nrow(B). 
  !
  !----------------------------------------------------------------------- 
  subroutine zcplsamub(nrow,ncol,job,a,ja,ia,s,b,jb,ib,c,jc,ic,nzmax,iw,ierr) 
    integer, intent(in) :: nrow, ncol, nzmax, job 
    complex(dp), intent(in) :: a(:), b(:), s
    integer, intent(in) :: ja(:), ia(:), jb(:), ib(:), iw(:)
    complex(dp), intent(inout) :: c(:)
    integer, intent(inout) :: jc(:), ic(:)
    integer, intent(out) :: ierr 
      
    integer :: ii, jj, kk, ka, kb, jpos, jcol, len
    complex(dp) :: scal 
    logical :: values
    
    values = (job .ne. 0) 
    len = 0
    ic(1) = 1 
    ierr = 0
    iw = 0
      
    do ii = 1, nrow 
       do ka=ia(ii), ia(ii+1)-1 
          if (values) scal = a(ka)
          jj = ja(ka)
          do kb=ib(jj),ib(jj+1)-1
             jcol = jb(kb)
             jpos = iw(jcol)
             if (jpos .eq. 0) then
                len = len+1
                if (len .gt. nzmax) then
                   ierr = ii
                   return
                endif
                jc(len) = jcol
                iw(jcol)= len
                if (values) c(len)  = c(len) + alpha*scal*b(kb)
             else
                if (values) c(jpos) = c(jpos) + alpha*scal*b(kb)
             endif
          end do
       end do 
       do kk = ic(ii), len
          iw(jc(kk)) = 0
       end do
       ic(ii+1) = len+1
    end do
  end subroutine zcplsamub 

  !-----------------------------------------------------------------------
  ! performs the matrix sum  C = s1*A+s2*B. 
  !-----------------------------------------------------------------------
  ! Upgrade 29/3 
  ! 
  !
  ! on entry:
  ! ---------
  ! nrow  = integer. The row dimension of A and B
  ! ncol  = integer. The column dimension of A and B.
  ! job   = integer. Job indicator. When job = 0, only the structure
  !                  (i.e. the arrays jc, ic) is computed and the
  !                  real values are ignored.
  !
  ! a,
  ! ja,
  ! ia   = Matrix A in compressed sparse row format.
  ! 
  ! s1   = complex scalar
  ! s2   = complex scalar
  !
  ! b, 
  ! jb, 
  ! ib  =  Matrix B in compressed sparse row format.
  !
  ! nzmax = integer. The  length of the arrays c and jc.
  !         amub will stop if the result matrix C  has a number 
  !         of elements that exceeds exceeds nzmax. See ierr.
  ! 
  ! on return:
  !----------
  ! c, 
  ! jc, 
  ! ic  = resulting matrix C in compressed sparse row sparse format.
  !     
  ! ierr  = integer. serving as error message. 
  !         ierr = 0 means normal return,
  !         ierr .gt. 0 means that amub stopped while computing the
  !         i-th row  of C with i=ierr, because the number 
  !         of elements in C exceeds nzmax.
  !
  ! work arrays:
  !------------
  ! iw  = integer work array of length equal to the number of
  !         columns in A.
  !
  subroutine zas1pls2b(nrow,ncol,job,a,ja,ia,s1,s2,b,jb,ib,c,jc,ic,nzmax,iw,ierr)
    integer, intent(in) :: nrow, ncol, nzmax, job 
    complex(dp), intent(in) :: a(:), b(:), s1, s2
    integer, intent(in) :: ja(:), ia(:), jb(:), ib(:), iw(:)
    complex(dp), intent(inout) :: c(:)
    integer, intent(inout) :: jc(:), ic(:)
    integer, intent(out) :: ierr 

    integer :: ii, kk, ka, kb, len, jcol, jpos
    logical :: values
     
    values = (job .ne. 0) 
    ierr = 0
    len = 0
    ic(1) = 1 
    iw = 0

    do ii = 1, nrow
       do ka = ia(ii), ia(ii+1)-1 
          len = len+1
          jcol    = ja(ka)
          if (len .gt. nzmax) then
             ierr = ii
             return
          end if   
          jc(len) = jcol 
          if (values) c(len)  = s1 * a(ka) 
          iw(jcol)= len
       end do

       do kb = ib(ii), ib(ii+1)-1
          jcol = jb(kb)
          jpos = iw(jcol)
          if (jpos .eq. 0) then
             len = len+1
             if (len .gt. nzmax) then
                ierr = ii
                return
             end if      
             jc(len) = jcol
             if (values) c(len)  = s2 * b(kb)
             iw(jcol)= len
          else
             if (values) c(jpos) = c(jpos) + s2 * b(kb)
          endif
       end do
       do kk = ic(ii), len
          iw(jc(kk)) = 0
       end do
       ic(ii+1) = len+1

  end subroutine zas1pls2b

  !----------------------------------------------------------------------
  ! copies the matrix a, ja, ia, into the matrix ao, jao, iao. 
  !----------------------------------------------------------------------
  ! on entry:
  !---------
  ! nrow  = row dimension of the matrix 
  ! a,
  ! ja,
  ! ia    = input matrix in compressed sparse row format. 
  ! ipos  = integer. indicates the position in the array ao, jao
  !         where the first element should be copied. Thus 
  !         iao(1) = ipos on return. 
  ! job   = job indicator. if (job .ne. 1) the values are not copies 
  !         (i.e., pattern only is copied in the form of arrays ja, ia).
  !
  ! on return:
  !----------
  ! ao,
  ! jao,
  ! iao   = output matrix containing the same data as a, ja, ia.
  !-----------------------------------------------------------------------
  subroutine zcopymat(nrow,a,ja,ia,ao,jao,iao,ipos,job) 
    integer, intent(in) :: nrow, ipos, job     
    complex(dp), intent(in) :: a(:)
    integer, intent(in) :: ja(:), ia(:)
    complex(dp), intent(inout) :: ao(:)
    integer, intent(in) :: jao(:), iao(:)

    integer kst, i, k 
      
    kst = ipos - ia(1) 
    do i = 1, nrow+1
       iao(i) = ia(i) + kst
    end do
    do k = ia(1), ia(nrow+1)-1
       jao(kst+k)= ja(k)
    end do
    if (job .ne. 0) then
      do k = ia(1), ia(nrow+1)-1
         ao(kst+k) = a(k)
      end do
    end if  
  end subroutine zcopymat

end module skit_blassm
