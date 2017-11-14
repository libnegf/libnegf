!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Green's functions.        !
!! Copyright (C) 2012                                                       !
!!                                                                          ! 
!! This file is part of libNEGF: a library for                              !
!! Non Equilibrium Green's Function calculation                             ! 
!!                                                                          !
!! Developers: Alessandro Pecchia, Gabriele Penazzi                         !
!! Former Conctributors: Luca Latessa, Aldo Di Carlo                        !
!!                                                                          !
!! libNEGF is free software: you can redistribute it and/or modify          !
!! it under the terms of the GNU Lesse General Public License as published  !
!! by the Free Software Foundation, either version 3 of the License, or     !
!! (at your option) any later version.                                      !
!!                                                                          !
!!  You should have received a copy of the GNU Lesser General Public        !
!!  License along with libNEGF.  If not, see                                !
!!  <http://www.gnu.org/licenses/>.                                         !  
!!--------------------------------------------------------------------------!


module outmatrix

  use ln_precision
  
  implicit none
  private

  public :: printmat, printmat_f, outmat, outmat_c, inmat_c
  public :: direct_out_c, direct_in_c

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
      character(20) :: STR

      n1=n
      k1=k

      if (k1.gt.99) then 
        k1=99 
      endif
      if (k1.gt.n1) then
        k1=n1
      endif

      write(STR,'("(",I2,"(ES10.2,ES9.2))")') k1

      write (UNIT=*,FMT=STR) ((gs(i1,i2),i2=1,k1),i1=1,k1)  
      print *," "

    end subroutine printmat
    !--------------------------------------------------
    subroutine printmat_fc(iu,FRM,gs,n,k)

      integer :: n,k,n1,k1
      complex(kind=dp) :: gs(n,n)
      integer :: i1,i2,iu,length
      character :: FRM(*)
      character(20) :: STR

      n1=n
      k1=k

      if (k1.gt.n1) then
        k1=n1
      endif
      if (k1.gt.150) then 
        k1=150 
      endif
      
      if (k1.gt.0.and.k1.le.99) then
        write(STR,'("(",I2,"(ES10.2,ES9.2))")') k1
      endif
      if (k1.gt.99.and.k1.le.150) then 
        write(STR,'("(",I3,"(ES10.2,ES9.2))")') k1
      endif

      write (iu,STR) ((gs(i1,i2),i2=1,k1),i1=1,k1)  
      write (iu,*) " "

    end subroutine printmat_fc
    !--------------------------------------------------
    subroutine printmat_fr(iu,FRM,gs,n,k)

      integer :: n,k,n1,k1
      real(kind=dp) :: gs(n,n)
      integer :: i1,i2,iu,length
      character :: FRM(*)
      character(20) :: STR

      k1=k
      n1=n

      if (k1.gt.n1) then
        k1=n1
      endif
      if (k1.gt.150) then
        k1=150
      endif

      if (k1.gt.0.and.k1.le.99) then
        write(STR,'("(",I2,"(ES10.2))")') k1
      endif
      if (k1.gt.99.and.k1.le.150) then 
        write(STR,'("(",I3,"(ES10.2))")') k1
      endif

      !write (6,'(/,a,a,a)') 'format:',STR,'!' 
      write (iu,FMT=STR) ((gs(i1,i2),i2=1,k1),i1=1,k1)  
      write (iu,*) " "

    end subroutine printmat_fr

    !--------------------------------------------------
    !subroutine printvect(A,n,lo,hi)
    ! 
    !  integer :: n,lo,hi
    !  complex(kind=dp) :: A(n)
    !  integer :: i
    !
    !  if (hi.gt.n.and.lo.lt.1) return
    !
    !  do i=lo,hi;  print '(I3,ES17.9)', i,abs(A(i)); 
    !  enddo
    !  print *," ";
    !  
    !end subroutine printvect
    !--------------------------------------------------

    subroutine  outmat(ndim,mat,lunit,format)
      !
      integer i,j,lunit
      integer ndim
      character(12), OPTIONAL :: format
      !  number of eigenstates, filled levels and half occupied states

      real(dp) mat(ndim,ndim)

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
    
    subroutine  outmat_c(lunit,fmt,A,n,ndim,acc_in)

      implicit none

      integer :: i,j
      integer :: ndim,n,lunit
      logical :: fmt
      complex(kind=dp) :: A(n,n)
      real(dp), optional :: acc_in
      real(dp) :: acc 

      acc = EPS
      if (present(acc_in)) acc = acc_in  

      if (fmt) then
        do j = 1, ndim 
          do i = 1, ndim
            if(abs(dble(A(i,j))).gt.acc.or.abs(aimag(A(i,j))).gt.acc) then
              write(lunit,'(2i8,(f20.10,f20.10))') i, j, A(i,j)
            endif
          enddo
        enddo
      else
        do j = 1, ndim 
          do i = 1, ndim
            if(abs(dble(A(i,j))).gt.acc.or.abs(aimag(A(i,j))).gt.acc) then
                    write(lunit) i, j, A(i,j)
            endif
          enddo
        enddo
      endif

    end subroutine outmat_c

    subroutine  inmat_c(lunit,fmt,A,n,ndim)

      implicit none

      integer :: i,j
      integer :: ndim,n,lunit
      logical :: fmt
      complex(kind=dp) :: A(n,n)
      complex(kind=dp) :: mat_el 

      if (fmt) then
        do
           read (lunit,*,end=100) i,j,mat_el
           A(i,j)=mat_el 
         enddo
      else
         do
           read (lunit,end=100) i,j,mat_el
           A(i,j)=mat_el 
         enddo
      endif

100   return

    end subroutine inmat_c

    subroutine direct_out_c(lunit,A,n)

       implicit none

       integer :: lunit, n
       complex(dp) :: A(n,n)

       integer :: i,j, k

       do j = 1, n
         do i = 1, n
        
            k = (j-1)*n+i

            write(lunit,rec = k)  A(i,j)

         end do
       end do 
         
    end subroutine direct_out_c


    subroutine direct_in_c(lunit,A,n)

       implicit none

       integer :: lunit, n
       complex(dp) :: A(n,n)

       integer :: i,j, k

       do j = 1, n
         do i = 1, n
        
            k = (j-1)*n+i

            read(lunit,rec = k)  A(i,j)

         end do
       end do 
         
    end subroutine direct_in_c

  
  end module outmatrix
