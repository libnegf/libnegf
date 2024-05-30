module blockpartition
     
    if (npl .eq. 0) then
       if (.not.allocated(negf%HS)) then
         stop "Error in init_structure: invoking block_partition but H not created"
         if (.not.associated(negf%HS(1)%H)) then
           stop "Error in init_structure: invoking block_partition but H not created"
         end if
       end if
       ! supposedly performs an internal block partitioning but it is not reliable.
       call log_allocate(plend_tmp, MAXNUMPLs)
       call block_partition(negf%HS(1)%H, surfend(1), contend, surfend, ncont, npl_tmp, plend_tmp)
       call find_cblocks(negf%HS(1)%H, ncont, npl_tmp, plend_tmp, surfstart, contend, cblk)
       call create_Tstruct(ncont, npl_tmp, plend_tmp, surfstart, surfend, contend, cblk, negf%str)
       call log_deallocate(plend_tmp)
     end if
     
  !////////////////////////////////////////////////////////////////////////
  ! RCM algorithm for reordering.
  !
  ! Actually not used because it is not suitable for contacted structures
  !////////////////////////////////////////////////////////////////////////
  subroutine reorder(mat)
    type(z_CSR) :: mat


    type(z_CSR) :: P, Tmp

    integer, dimension(:), allocatable :: perm
    integer :: i, nrow

    nrow=mat%nrow

    call log_allocate(perm,nrow)

    call genrcm(nrow, mat%nnz, mat%rowpnt, mat%colind, perm)

    call create(P,nrow,nrow,nrow)

    do i=1,nrow
       P%nzval(i)=1
       P%colind(i)=perm(i)
       P%rowpnt(i)=i
    enddo
    P%rowpnt(nrow+1)=nrow+1


    call create(Tmp,nrow,nrow,mat%nnz)

    call zamub_st(P,mat,Tmp)

    call ztransp_st(P)

    call zamub_st(Tmp,P,mat)

    call destroy(P,Tmp)

    call log_deallocate(perm)

  end subroutine reorder

  !>  Authomatic Block partitioning. The matrix must be already sorted.
  subroutine block_partition(mat,nrow,cont_end,surf_end,ncont,nbl,blks)

    !> The matrix to be partitioned.
    type(z_CSR), intent(in) :: mat
    !> The number of row to partition.
    integer, intent(in) :: nrow
    !> The indices indicating the end of the contact.
    integer, dimension(:), intent(in) :: cont_end
    !> The indices indicating the end of the scattering region surface
    !> (last orbitals before corresponding contact.)
    integer, dimension(:), intent(in) :: surf_end
    !> The number of contacts.
    integer, intent(in) :: ncont
    !> The number of blocks.
    integer, intent(out) :: nbl
    !> The array with the end index for each block.
    integer, dimension(:), intent(inout) :: blks

    integer :: j, k, i
    integer :: i1, i2

    integer :: rn, rnold, tmax, rmax, maxmax
    integer :: dbuff, minsize, minv, maxv

    minsize = 0
     do i1 = 1, ncont
        maxv = 0
        minv = 400000000
        do k = surf_end(i1)+1, cont_end(i1)
           do i = mat%rowpnt(k), mat%rowpnt(k+1)-1
            if (mat%colind(i).le.nrow .and.  mat%colind(i).lt.minv) minv = mat%colind(i)
            if (mat%colind(i).le.nrow .and.  mat%colind(i).gt.maxv) maxv = mat%colind(i)
           end do
        end do
        if (maxv-minv+1 .gt. minsize) minsize = maxv - minv + 1
    end do

    ! The current algorithm does not work when the minimum block
    ! size is 1. We fix the minimum possible size to 2 as temporary fix.
    minsize = max(minsize, 2)

    ! Find maximal stancil of the matrix and on which row
    !  ( Xx     )
    !  ( xXxx   )
    !  (  xXxxx )  <- maxmax = 3 ; rmax = 3
    !  (   xXx  )
    maxmax = 0
    do j=1,nrow
       tmax = 0
       do i = mat%rowpnt(j), mat%rowpnt(j+1) - 1
           if ( mat%colind(i).le.nrow .and. (mat%colind(i)-j) .gt. tmax) then
                tmax = mat%colind(i)-j
           endif
       enddo

       if(tmax .gt. maxmax) then
          maxmax = tmax
          rmax = j
       endif

       dbuff = maxmax        ! dbuff should be linked to maxmax
       minsize = max((dbuff+1)/2,minsize)
    enddo

    ! Define central block
    rn = rmax - maxmax/2 - dbuff

    if(rn-dbuff.ge.0) then

       blks(1) = rn-1  ! fine del blocco precedente

       nbl = 1

       do

          do j = rn, minsize, -1

             rnold = rn
             i1 = mat%rowpnt(j-minsize+1)
             i2 = mat%rowpnt(j+1) - 1

             !k = maxval(mat%colind(i1:i2))
             k = 0
             do i = i1, i2
                if ( mat%colind(i).le.nrow .and. (mat%colind(i)) .gt. k) then
                   k = mat%colind(i)
                endif
             enddo

             if(k.lt.rn) then
                rn = j
                nbl = nbl + 1
                if (nbl.gt.MAXNUMPLs) call errormsg()
                blks(nbl) = j-1 ! fine del blocco precedente
                exit
             endif
          enddo

          if(rn.le.minsize .or. rnold.eq.rn) then
             exit
          endif

       enddo

       rn = rmax - maxmax/2 - dbuff

    else
       nbl= 0
       rn = 1

    endif

    do

       do j = rn, nrow-minsize+1, 1

          rnold = rn
          i1 = mat%rowpnt(j)
          i2 = mat%rowpnt(j+minsize) - 1

          !k = minval(mat%colind(i1:i2))
          k = nrow
          do i = i1, i2
             if ( mat%colind(i).le.nrow .and. (mat%colind(i)) .lt. k) then
                k = mat%colind(i)
             endif
          enddo


          if(k.gt.rn) then
             rn = j
             nbl = nbl + 1
             if (nbl.gt.MAXNUMPLs) call errormsg()
             blks(nbl) = j-1 ! fine del blocco
             exit
          endif
       enddo

       if(nrow-rn.le.minsize .or. rnold.eq.rn) then
          exit
       endif
    enddo

    nbl = nbl + 1
    if (nbl.gt.MAXNUMPLs) call errormsg()
    blks(nbl) = nrow

    ! Sorting blocks

    do i = 1, nbl
       do j = i+1, nbl
          if(blks(j).lt.blks(i)) then
             k = blks(i)
             blks(i) = blks(j)
             blks(j) = k
          endif
       enddo
    enddo


  end subroutine block_partition

!----------------------------------------------------------------------------
  subroutine errormsg()

     write(*,*) "ERROR: Maximum number of PLs exceeded"
     write(*,*) "increase the value of MAXNUMPLS in libnegf.F90"
     write(*,*) "and recompile the library"

     STOP !should rise an exception

  end subroutine errormsg

!----------------------------------------------------------------------------
  subroutine find_cblocks(mat ,ncont, nbl, PL_end, surf_start, cont_end, cblk)
    type(z_CSR), intent(in) :: mat
    integer, intent(in) :: ncont
    integer, intent(in) :: nbl
    integer, dimension(:), intent(in) :: PL_end
    integer, dimension(:), intent(in) :: surf_start
    integer, dimension(:), intent(in) :: cont_end
    integer, dimension(:), allocatable, intent(out) :: cblk

    integer :: j1,k,i,min,max
    integer, dimension(:), allocatable :: PL_start

    call log_allocate(PL_start,nbl)
    call log_allocate(cblk,ncont)

    PL_start(1) = 1

    do i = 2, nbl
       PL_start(i) = PL_end(i-1) + 1
    enddo


    do j1 = 1, ncont

       max = 0
       min = 400000000

       do k = surf_start(j1), cont_end(j1)

          do i = mat%rowpnt(k), mat%rowpnt(k+1)-1

             if (mat%colind(i).le.PL_end(nbl) .and.  mat%colind(i).lt.min) min = mat%colind(i)
             if (mat%colind(i).le.PL_end(nbl) .and.  mat%colind(i).gt.max) max = mat%colind(i)

          end do

       end do

       do k = 1, nbl

          if( max .le. PL_end(k) ) then
             cblk(j1) = k

             if( min .ge. PL_start(k) ) then
                exit
             else
                write(*,*) "(LibNEGF) Partitioning:"
                write(*,*) "Number of blocks: ",nbl
                write(*,*) "PL_end: ",PL_end(1:nbl)
                write(*,*) "Contact interaction: ",cblk(j1)
                write(*,'(a,i3,a)') " ERROR: contact",j1," interacting with more than one block"
                write(*,*) "min ",min,"max ",max
                stop
             end if

          end if

       end do

    end do

    call log_deallocate(PL_start)

  end subroutine find_cblocks

!----------------------------------------------------------------------------

  Subroutine sort(blks, Ipt)
    ! *
    ! ***********************************
    ! * Sort Array X(:) in ascendent order.
    ! * If present Ipt, a pointer with the
    ! * changes is returned in Ipt.
    ! ***********************************

    Integer, Intent (inout) :: blks(:)
    Integer, Intent (out), Optional :: Ipt(:)

    Integer :: Rtmp
    Integer :: i, j

    If (Present(Ipt)) Then
       Forall (i=1:Size(blks)) Ipt(i) = i

       Do i = 2, Size(blks)
          Rtmp = blks(i)
          Do j = i-1, 1, -1
             If (Rtmp < blks(j)) Then
                blks(j+1) = blks(j)
                call Swap(blks, j, j+1)
             Else
                Exit
             End If
          End Do
          blks(j+1) = Rtmp
       End Do
    Else
       Do i = 2, Size(blks)
          Rtmp = blks(i)
          Do j = i-1, 1, -1
             If (Rtmp < blks(j)) Then
                blks(j+1) = blks(j)
             Else
                Exit
             End If
          End Do
          blks(j+1) = Rtmp
       End Do
    End If

    Return
  End Subroutine sort

  ! ***********************************
  ! * Swaps elements I and J of array X(:).
  ! ***********************************
  Subroutine Swap(X, i, j)

    Integer, Intent (inout) :: X(:)
    Integer, Intent (in) :: i, j

    Integer :: Itmp

    Itmp = X(I)
    X(I) = X(J)
    X(J) = Itmp

    Return
  End Subroutine Swap


  subroutine aggregate_vec(v_in, thres, v_out, start_idx, end_idx)
      real(dp), dimension(:), intent(in) :: v_in
      real(dp), intent(in) :: thres
      real(dp), dimension(:), allocatable, intent(out) :: v_out
      integer, dimension(:), allocatable, intent(out) :: start_idx, end_idx

      real(dp) :: avg
      integer :: i, rs, re

      allocate(start_idx(0))
      allocate(end_idx(0))
      allocate(v_out(0))

      start_idx = [start_idx, 1]
      do i = 1, size(v_in)-1
         if (abs(v_in(i+1) - v_in(i)) > thres) then
            start_idx = [start_idx, i+1]
         endif
      end do

      do i = 1, size(v_in)-1
         if (abs(v_in(i+1) - v_in(i)) > thres) then
                 end_idx = [end_idx, i]
         endif
      end do
      end_idx = [end_idx, size(v_in)]

      do i = 1, size(start_idx)
         rs = start_idx(i)
         re = end_idx(i)
         avg = sum(v_in(rs:re))/real(size(v_in(rs:re)), dp)
         v_out = [v_out, avg]
      end do

   end subroutine aggregate_vec

!--------------------------------------------------------------------
  subroutine negf_partition_info(negf)
      type(Tnegf) :: negf

      integer :: i

      write(*,*) "(LibNEGF) Partitioning:"
      write(*,*) "Number of blocks: ",negf%str%num_Pls
      !write(*,*) negf%str%mat_PL_end(:)
      write(*,*) "Contact interactions:",negf%str%cblk(:)

      open(1001,file='blocks.dat')
        write(1001,*) 1
        do i = 1, negf%str%num_Pls
           write(1001,*)  negf%str%mat_PL_end(i)
        enddo
      close(1001)

  end subroutine negf_partition_info

end module blockpartition

