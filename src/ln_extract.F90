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


 module negf_ln_extract
  
  use negf_ln_precision
  use negf_ln_allocation
  use negf_ln_constants
  use negf_mat_def
  use negf_sparsekit_drv
  use negf_ln_structure, only: TStruct_Info
  use negf_lib_param

  implicit none

  private

  public :: extract_cont 

contains

  !--------------------------------------------------------------------------------------------

  subroutine extract_cont(negf)

    Type(Tnegf) :: negf

    Integer :: i, ncont, i1, i2, j1, j2
    Integer :: cstart(MAXNCONT),cend(MAXNCONT)
    Integer :: ncdim(MAXNCONT), surfdim(MAXNCONT)

    ncont =  negf%str%num_conts
    do i=1,ncont
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

end  module negf_ln_extract
