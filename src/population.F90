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


module population
  use ln_precision
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
    real(dp) :: qtot
    complex(dp) :: dd

    qmulli=0.0_dp
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
    qtot = 0.0_dp
    do ii = 1, nrow
       write(11,*) ii,qmulli(ii)
       qtot = qtot+qmulli(ii)
    enddo
    close(11)
    
    write(*,*) 'qtot=',qtot
  
  end subroutine mulliken

  ! -------------------------------------------------------------

end module population
