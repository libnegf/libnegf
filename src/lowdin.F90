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


module lowdin
  
  use ln_precision

  implicit none
  private

  public :: ssqr, lowdin_trans


contains

  subroutine ssqr(OP,S,Shalf,dim)

    implicit none

    complex(kind=dp), PARAMETER :: alfa= (1.d0,0.d0)
    complex(kind=dp), PARAMETER :: beta= (0.d0,0.d0)

    integer :: dim
    complex(kind=dp), DIMENSION(dim,dim) :: S,Shalf
    character(1) :: OP

    integer :: LWORK, info, i,j
    real(kind=dp), DIMENSION(dim) :: W
    real(kind=dp), DIMENSION(3*dim) :: RWork
    complex(kind=dp), DIMENSION(2*dim) :: Work
    complex(kind=dp), DIMENSION(dim,dim) :: U, D

    LWORK=2*dim

    ! The Lowdin transformation is obtained by
    ! digonalizing S:   S'=U^ S U
    ! taking all eigenvalues squared and going back S = U (S')^(1/2) U^

    ! First part consists of digonalizing real-symmetric matrix S
    ! call LAPACK routine

    U(:,:)=S(:,:)

    call ZHEEV( 'V', 'U', dim, U, dim, W, Work, LWORK, RWORK, info)

    if (info.ne.0) then

       write(*,*) 'Diagonalization of S failed'

    else   

       do i=1,dim
          do j=1,dim

             if (OP.eq.'+') then
                !---- Compute S^(1/2) ---------------------
                D(j,i)= U(i,j)*dsqrt(W(j)) 
             else
                !---- Compute S^(-1/2) --------------------
                D(j,i)= U(i,j)/dsqrt(W(j))
             endif

          enddo
       enddo

       call zgemm('N','N', dim,dim,dim, alfa, U, dim, D, dim, beta, Shalf, dim)

    endif

  end subroutine ssqr
  !----------------------------------------------------------------------
  !----------------------------------------------------------------------
  subroutine ssqr_iterative(OP,S,Shalf,dim)

    implicit none

    real(kind=dp), PARAMETER :: alfa= 1.d0
    real(kind=dp), PARAMETER :: beta= 0.d0

    integer :: dim
    real(kind=dp), DIMENSION(dim,dim) :: S,Shalf
    character(1) :: OP

    integer :: i
    real(dp) :: r
    
    Shalf = 0.d0

    do i=1,dim
       S(i,i) = S(i,i) - 1.d0
       Shalf(i,i) = 1.d0
    enddo
    
    if (OP.eq.'+') r = 0.5d0
    if (OP.eq.'-') r = -0.5d0    

 
    !Shalf = Shalf + (r k) S

  end subroutine ssqr_iterative
  !----------------------------------------------------------------------
  ! Implement Lowdin transform: HL=S^(OP 1/2) H S^(OP 1/2) --------------
  ! OP is a character, can be '+' or '-' 

  subroutine lowdin_trans(OP,H,S,HL,dim)

    implicit none
    complex(kind=dp), PARAMETER :: alfa= (1.d0,0.d0)
    complex(kind=dp), PARAMETER :: beta= (0.d0,0.d0)

    integer :: dim
    complex(kind=dp), DIMENSION(dim,dim) :: S,H,HL
    character(1) :: OP

    complex(kind=dp), DIMENSION(dim,dim) :: Sh,TMP

    call ssqr(OP, S, Sh, dim)

    call zgemm('N','N', dim,dim,dim, alfa, H, dim, Sh, dim, beta, TMP, dim)
    call zgemm('N','N', dim,dim,dim, alfa, Sh, dim,TMP, dim, beta, HL, dim)

  end subroutine lowdin_trans


end module lowdin
