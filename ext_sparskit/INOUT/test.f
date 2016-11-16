c----end-of-dump--------------------------------------------------------
c-----+--------------------------------------------------------------------+
      program test
      
      integer k,nn
      complex*16 a(30)
      character ptrfmt*16,indfmt*16,valfmt*40

      ifmt=104
      nnz=7     

      do k=1,nnz
         a(k) = (3.42342d0,-23.34503d0)
c          a(k) = -3.1233d0
      enddo
       if (ifmt .ge. 100) then
         ihead = ifmt/100
         ifmt = ifmt-100*ihead
         len = ihead+ifmt+3
         nperli = 80/len
c     
         if (len .le. 9 ) then
            write(valfmt,102) nperli,len,ifmt
         elseif (ifmt .le. 9) then
            write(valfmt,103) nperli,len,ifmt
         else 
            write(valfmt,104) nperli,len,ifmt
         endif
c     
 102   format(1h(,i2,6h(1h(,F,i1,1h.,i1,6h,1h,,F,i1,1h.,i1,5h,2h) ,2h)))
 103   format(1h(,i2,6h(1h(,F,i2,1h.,i1,6h,1h,,F,i2,1h.,i1,5h,2h) ,2h)))
 104   format(1h(,i2,6h(1h(,F,i2,1h.,i2,6h,1h,,F,i2,1h.,i2,5h,2h) ,2h))) 
C     
      else
 
         len = ifmt + 7
         lent= 2*len+3
         nperli = 80/lent
         if (nperli .le. 9) then
            if (len .le. 9 ) then
               write(valfmt,105) nperli,len,ifmt,len,ifmt
            elseif (ifmt .le. 9) then
               write(valfmt,106) nperli,len,ifmt,len,ifmt
            else 
               write(valfmt,107) nperli,len,ifmt,len,ifmt
            endif
         else 
            if (len .le. 9 ) then
               write(valfmt,108) nperli,len,ifmt,len,ifmt
            elseif (ifmt .le. 9) then
               write(valfmt,109) nperli,len,ifmt,len,ifmt
            else 
               write(valfmt,110) nperli,len,ifmt,len,ifmt
            endif
         endif
      endif
c
 105  format(1h(,i1,6h(1h(,D,i1,1h.,i1,6h,1h,,D,i1,1h.,i1,5h,2h) ,2h)))
 106  format(1h(,i1,6h(1h(,D,i2,1h.,i1,6h,1h,,D,i2,1h.,i1,5h,2h) ,2h)))
 107  format(1h(,i1,6h(1h(,D,i2,1h.,i2,6h,1h,,D,i2,1h.,i2,5h,2h) ,2h)))
 108  format(1h(,i2,6h(1h(,D,i1,1h.,i1,6h,1h,,D,i1,1h.,i1,5h,2h) ,2h)))
 109  format(1h(,i2,6h(1h(,D,i2,1h.,i1,6h,1h,,D,i2,1h.,i1,5h,2h) ,2h)))
 110  format(1h(,i2,6h(1h(,D,i2,1h.,i2,6h,1h,,D,i2,1h.,i2,5h,2h) ,2h)))
c
      iounit=30
      open(iounit,file='test.dat')
      k=mod(nnz, nperli)
      nn=nnz-k     
c 
      write(iounit,valfmt,err=1000) ( a(i), i = 1, nn )
      write(valfmt,ix) k,len,ifmt,len,ifmt
      write(iounit,valfmt,err=1000) ( a(i), i = 1, k )
      close(iounit)
     
c      open(30,file='test.dat')
c      read(30,*) (a(k),k=1,3)
c      write(*,103) (a(k),k=1,3)
c      close(30)
 1000  continue

      a=0.d0
      do k=1,nnz
         write(*,*) a(k) 
      enddo     
      open(iounit,file='test.dat')
      read(iounit,*) (a(i), i=1,nnz)
      close(iounit)
     
      do k=1,nnz
         write(*,*) a(k) 
      enddo     

c102   format(4(1h(,D9.2,D10.2,1h),1x) )
c103   format(1h(,e22.14,1h,,e22.14,1h))
      end program
c----end-of-dump--------------------------------------------------------
