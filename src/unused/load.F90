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


module load
 
  use mat_def
  implicit none
  private

  public :: load_HS


contains

  subroutine load_HS(H,S,fmt)
    type(z_CSR) :: H,S  
    integer :: nrow, ncol, nnz
    logical :: fmt

    ! --------------------------------------------
    if (fmt) then
       open(11, file='H.dat', FORM='FORMATTED')
       read(11,*) nrow
       read(11,*) ncol  
       read(11,*) nnz
       call create(H,nrow,ncol,nnz)
       
       rewind 11
       call read_mat(11,H,.true.)
       close(11)
    else
       open(11, file='H.dat', FORM='UNFORMATTED')
       read(11) nrow
       read(11) ncol  
       read(11) nnz
       call create(H,nrow,ncol,nnz)
       
       rewind 11
       call read_mat(11,H,.false.)
       close(11)
    endif


    if (fmt) then
       open(11, file='S.dat', FORM='FORMATTED')
       read(11,*) nrow
       read(11,*) ncol  
       read(11,*) nnz
       call create(S,nrow,ncol,nnz)
       
       rewind 11
       call read_mat(11,S,.true.)
       close(11)  
    else
       open(11, file='S.dat', FORM='UNFORMATTED')       
       read(11) nrow
       read(11) ncol  
       read(11) nnz
       call create(S,nrow,ncol,nnz)
       
       rewind 11
       call read_mat(11,S,.false.)
       close(11)       
    endif


  end subroutine load_HS


end module load
