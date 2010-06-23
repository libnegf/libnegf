module population
  use precision
  use mat_def

  implicit none
  private
  public :: mulliken

  ! -------------------------------------------------------------
contains

  subroutine mulliken(DensMat,S,qmulli)
    type(z_CSR) :: S, DensMat
    real(dp), dimension(:) :: qmulli


    integer :: ii, ka, jj, kb, jcol, nrow
    real(dp) :: dd, qtot

    qmulli=0.d0
    nrow = size(qmulli)

    do ii=1, DensMat%nrow
       do ka=DensMat%rowpnt(ii), DensMat%rowpnt(ii+1)-1 
          dd = DensMat%nzval(ka)
          jj = DensMat%colind(ka)
          
          do kb=S%rowpnt(jj),S%rowpnt(jj+1)-1
             jcol = S%colind(kb)
             if (jcol .eq. ii .and. jcol.le.nrow) then
                qmulli(jcol) = qmulli(jcol) + real(dd*S%nzval(kb))
             endif
          enddo
       enddo
    enddo
    
    open(11,file='qmulli.dat')
    qtot = 0.d0
    do ii = 1, nrow
       write(11,*) ii,qmulli(ii)
       qtot = qtot+qmulli(ii)
    enddo
    close(11)
    
    write(*,*) 'qtot=',qtot
  
  end subroutine mulliken

  ! -------------------------------------------------------------

end module population
