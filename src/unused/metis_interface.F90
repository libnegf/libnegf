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


 module negf_metis_interface

  use negf_mat_def
  
  implicit none
  private

  public ::  METISpartition


contains

  subroutine METISpartition(mat,nbl,part)
    type(z_CSR) :: mat
    
    integer :: nbl
    integer :: n
    
    integer :: volume
    integer :: numflag
    integer :: wghts
    integer :: wgtflag

    integer, dimension(:) :: part


    integer, dimension(:), allocatable :: options
    integer, dimension(:), allocatable :: vwgt  
    integer, dimension(:), allocatable :: vsize
    

    external METIS_PartGraphVKway

    numflag = 1
    wghts = 0
    wgtflag = 0
    n = mat%nrow

    call log_allocate(vwgt, 0)
    call log_allocate(vsize, 0)
    call log_allocate(options, 5)
    options(1) = 0


    !call METIS_PartGraphVKway(n, mat%rowpnt, mat%colind, vwgt, vsize, wgtflag, &
    !                          numflag, nbl, options, volume, part)

    call METIS_PartGraphKway(n, mat%rowpnt, mat%colind, vwgt, vsize, wgtflag, &
                              numflag, nbl, options, volume, part)    

    call log_deallocate(vwgt)
    call log_deallocate(vsize)
    call log_deallocate(options)
    

  end subroutine METISpartition


  !----------------------------------------------------------

end  module negf_metis_interface
