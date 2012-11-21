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


program testreorder

  use ln_precision
  use ln_constants
  use mat_def
  use ln_allocation
  use ln_structure
  use input_output
  use sparsekit_drv
  use inversions
  use iterative
  use libnegf
  use lib_param


  implicit none

  Type(Tnegf), target :: negf
  Type(Tnegf), pointer :: pnegf

  integer :: iofile
  integer :: nbl, i
  integer, dimension(:), allocatable :: part, PL
  integer, dimension(2) :: tmp

  pnegf => negf

  print*,'(main) testreorder'

  negf%file_re_H='H_real.dat'
  negf%file_im_H='H_imm.dat'
  negf%file_struct='driver'
  negf%verbose = 10
  !negf%eneconv = HAR  ! to convert Kb 
  negf%isSid = .true.  

  negf%form%formatted = .true.
  negf%form%type = 'UPT'
  negf%form%fmt = 'U'

  call init_negf(pnegf)

  iofile = 101

  !open(iofile,file='H_orig.dat')
  !call zprint_csrcoo(iofile,pnegf%H,'r')
  !close(iofile)
print*, pnegf%H%nrow
print*, maxval(pnegf%H%colind)

print*, 'before reorder'

  call reorder(pnegf%H)

print*, 'reorder done'

print*, pnegf%H%nrow
print*, maxval(pnegf%H%colind)

print*, 'save reordered H'

  open(iofile,file='H_reord.dat')
  call zprint_csrcoo(iofile,pnegf%H,'r')
  close(iofile)

print*, 'before partition'

  call log_allocate(part,pnegf%H%nrow+1)

  call partition2(pnegf%H,nbl,part)

  call log_allocate(PL,nbl)

  PL = part(1:nbl)

  call log_deallocate(part)

  call sort(PL)

  open(iofile,file='H_part.dat')

  do i=1,nbl
     write(iofile,*) PL(i)
  enddo

  close(iofile)

print*, 'partition done'


print*, 'before create T'
  call kill_TStruct(pnegf%str)
  call create_Tstruct(0, nbl, PL, tmp, tmp, tmp, pnegf%str)
  call log_deallocate(PL)

print*,'H%nrow',pnegf%H%nrow
print*,'H%ncol',pnegf%H%ncol
print*,'H%nnz',pnegf%H%nnz,size(pnegf%H%nzval)
print*,'S%nrow',pnegf%S%nrow
print*,'S%ncol',pnegf%S%ncol
print*,'S%nnz',pnegf%S%nnz,size(pnegf%S%nzval)


print*, 'compute dos'

 ! call compute_dos(pnegf)


  call destroy_negf(pnegf)




  call writepeakinfo(6)
  call writememinfo(6)


end program testreorder
