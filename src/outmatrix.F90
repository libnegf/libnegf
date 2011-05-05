module outmatrix

  use ln_precision
  
  implicit none
  private

  public :: printmat, printmat_f, outmat, outmat_c


  interface printmat_f
    module procedure printmat_fc,printmat_fr
  end interface
  
  contains
!--------------------------------------------------

    subroutine printmat(gs,n,k)

      integer :: n,k,n1,k1
      complex(kind=dp) :: gs(n,n)
      integer :: i1,i2,length
      character :: kk(2), line(2)
      character, allocatable :: STR(:)

      n1=n
      k1=k

      if (k1.gt.99) then 
        k1=99 
      endif
      if (k1.gt.n1) then
        k1=n1
      endif

      write (UNIT=line,FMT='(I2)') k1
      read (UNIT=line,FMT='(A2)') kk

      length= LEN('('//kk//'(ES10.2,ES9.2))')
      allocate(STR(length))
      STR='('//kk//'(ES10.2,ES9.2))'

      write (UNIT=*,FMT=STR) ((gs(i1,i2),i2=1,k1),i1=1,k1)  
      print *," "

      deallocate(STR)

    end subroutine printmat
    !--------------------------------------------------
    subroutine printmat_fc(iu,FRM,gs,n,k)

      integer :: n,k,n1,k1
      complex(kind=dp) :: gs(n,n)
      integer :: i1,i2,iu,length
      character :: kk(3), line(3),FRM(*)
      character, allocatable :: STR(:)

      n1=n
      k1=k

      if (k1.gt.n1) then
        k1=n1
      endif
      if (k1.gt.150) then 
        k1=150 
      endif
      if (k1.gt.0.and.k1.le.99) then
        write (UNIT=line,FMT='(I2)') k1
        read (UNIT=line,FMT='(A2)') kk
      endif
      if (k1.gt.99.and.k1.le.150) then 
        write (UNIT=line,FMT='(I3)') k1
        read (UNIT=line,FMT='(A3)') kk
      endif

      length= LEN('('//kk//'(ES10.2,ES9.2))')
      allocate(STR(length))
      STR='('//kk//'(ES10.2,ES9.2))'

      write (iu,STR) ((gs(i1,i2),i2=1,k1),i1=1,k1)  
      write (iu,*) " "

      deallocate(STR)

    end subroutine printmat_fc
    !--------------------------------------------------
    subroutine printmat_fr(iu,FRM,gs,n,k)

      integer :: n,k,n1,k1
      real(kind=dp) :: gs(n,n)
      integer :: i1,i2,iu,length
      character :: kk(3), line(3),FRM(*)
      character, allocatable :: STR(:)

      k1=k
      n1=n

      if (k1.gt.n1) then
        k1=n1
      endif
      if (k1.gt.150) then
        k1=150
      endif

      if (k1.gt.0.and.k1.le.99) then
        write (UNIT=line,FMT='(I2)') k1
        read (UNIT=line,FMT='(A2)') kk
      endif
      if (k1.gt.99.and.k1.le.150) then 
        write (UNIT=line,FMT='(I3)') k1
        read (UNIT=line,FMT='(A3)') kk
      endif

      length= LEN('('//kk//'(ES10.2)')
      allocate(STR(length))
      STR='('//kk//'(ES10.2))'

      !write (6,'(/,a,a,a)') 'format:',STR,'!' 
      write (iu,FMT=STR) ((gs(i1,i2),i2=1,k1),i1=1,k1)  
      write (iu,*) " "

      deallocate(STR)

    end subroutine printmat_fr

    !--------------------------------------------------
    subroutine printvect(A,n,lo,hi)

      integer :: n,lo,hi
      complex(kind=dp) :: A(n)
      integer :: i

      if (hi.gt.n.and.lo.lt.1) return

      do i=lo,hi;  print '(I3,ES17.9)', i,abs(A(i)); 
      enddo
      print *," ";
      
    end subroutine printvect
    !--------------------------------------------------

    subroutine  outmat(ndim,mat,lunit,format)
      !
      integer i,j,k,l,lunit
      integer ndim
      character(12), OPTIONAL :: format
      !  number of eigenstates, filled levels and half occupied states

      real*8 mat(ndim,ndim), nel

      !format='U'

      if (.not.PRESENT(format)) then
        format='FORMATTED'
      endif
      !   write(*,'(a1,a,a1)') '%',format,'%'

      if (format(1:1).eq.'F') then

        write (lunit,'(a,i8,a,i8)') '# matrix dimension = ', ndim, ' x ', ndim
        write (lunit,'(a)') '#'

        do i = 1, ndim
          do j = 1, i
            if (dabs(mat(i,j))>EPS) write(lunit,10) i, j, mat(i,j)
          enddo
        enddo

      else

        write (lunit) '# matrix dimension = ', ndim, ' x ', ndim
        write (lunit) '#'

        do i = 1, ndim
          do j = 1, i
            if (dabs(mat(i,j))>EPS) write(lunit) i, j, mat(i,j)
          enddo
        enddo

      endif

      ! format specification
      !
10    format(2i8,f20.10)

    END subroutine outmat
    !--------------------------------------------------
    
    subroutine  outmat_c(lunit,fmt,A,n,ndim)

      implicit none

      integer :: i,j
      integer :: ndim,n,lunit
      logical :: fmt
      complex(kind=dp) :: A(n,n)

      if (fmt) then
        do i = 1, ndim 
          do j = 1, ndim
            if(abs(dble(A(i,j))).gt.EPS.or.abs(dimag(A(i,j))).gt.EPS) then
              write(lunit,'(2i8,(f20.10,f20.10))') i, j, A(i,j)
            endif
          enddo
        enddo
      else
        do i = 1, ndim 
          do j = 1, ndim
            if(abs(dble(A(i,j))).gt.EPS.or.abs(dimag(A(i,j))).gt.EPS) then
              write(lunit) i, j, A(i,j)
            endif
          enddo
        enddo
      endif

    end subroutine outmat_c

  
  end module outmatrix
