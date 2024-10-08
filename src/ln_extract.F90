!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Greens functions.         !
!! Copyright (C) 2012 - 2026                                                !
!!                                                                          !
!! This file is part of libNEGF: a library for                              !
!! Non Equilibrium Green's Functions calculations                           !
!!                                                                          !
!! Developers: Alessandro Pecchia, Daniele Soccodato                        !
!! Former Contributors: Gabriele Penazzi, Luca Latessa, Aldo Di Carlo       !
!!                                                                          !
!! libNEGF is free software: you can redistribute and/or modify it          !
!! under the terms of the GNU Lesser General Public License as published    !
!! by the Free Software Foundation, either version 3 of the License, or     !
!! (at your option) any later version.                                      !
!!                                                                          !
!!  You should have received a copy of the GNU Lesser General Public        !
!!  License along with libNEGF.  If not, see                                !
!!  <http://www.gnu.org/licenses/>.                                         !
!!--------------------------------------------------------------------------!


module ln_extract

  use ln_precision
  use ln_allocation
  use ln_constants
  use mat_def
  use sparsekit_drv
  use ln_structure, only: TStruct_Info
  use lib_param

  implicit none

  private

  public :: extract_cont
  public :: destroy_contact_matrices

contains

  !--------------------------------------------------------------------------------------------

  subroutine extract_cont(negf)

    Type(Tnegf) :: negf

    Integer :: i, ncont, i1, i2, j1, j2
    Integer :: cstart(MAXNCONT),cend(MAXNCONT)
    Integer :: ncdim(MAXNCONT), surfdim(MAXNCONT)

    ncont =  negf%str%num_conts
    do i=1,ncont
       !print*, '(Cont)',i,negf%str%mat_C_start(i),negf%str%mat_B_start(i),negf%str%mat_C_end(i)
       cstart(i) = negf%str%mat_B_start(i)
       cend(i)   = negf%str%mat_C_end(i)
       ncdim(i)  = cend(i)-cstart(i)+1
       surfdim(i) = negf%str%mat_C_start(i) - negf%str%mat_B_start(i)
    enddo

    do i=1,ncont
       !print*, '(ext_HC)',i,cstart(i),cend(i),surfdim(i),ncdim(i)
       call extract(negf%H,cstart(i),cend(i),cstart(i),cend(i),negf%cont(i)%HC)
       !print*, '(ext_SC)',i,cstart(i),cend(i),surfdim(i),ncdim(i)
       call extract(negf%S,cstart(i),cend(i),cstart(i),cend(i),negf%cont(i)%SC)
    enddo

    do i=1,ncont
       !print*, '(int) extract central-contact',i
       i1 = negf%str%mat_PL_start( negf%str%cblk(i) )
       i2 = negf%str%mat_PL_end( negf%str%cblk(i) )
       j1 = cstart(i);
       j2 = j1+(ncdim(i)+surfdim(i))/2-1 !Note this is Surf+1PL
       !print*, 'block HMC:',i1,i2,j1,j2
       call extract(negf%H,i1,i2,j1,j2,negf%cont(i)%HMC)
       !print*, 'block SMC:',i1,i2,j1,j2
       call extract(negf%S,i1,i2,j1,j2,negf%cont(i)%SMC)
    enddo

  end subroutine extract_cont

  !--------------------------------------------------------------------
  !> Destroy matrices created runtime in libnegf
  subroutine destroy_contact_matrices(negf)
    type(Tnegf) :: negf
    integer :: i

    if (allocated(negf%cont)) then
       do i = 1, size(negf%cont)
          if (allocated(negf%cont(i)%HC%val)) call destroy(negf%cont(i)%HC)
          if (allocated(negf%cont(i)%SC%val)) call destroy(negf%cont(i)%SC)
          if (allocated(negf%cont(i)%HMC%val)) call destroy(negf%cont(i)%HMC)
          if (allocated(negf%cont(i)%SMC%val)) call destroy(negf%cont(i)%SMC)
       enddo
    end if

  end subroutine destroy_contact_matrices

end module ln_extract
