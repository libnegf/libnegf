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


program testneqint

  use negf_ln_precision
  use negf_ln_constants
  use negf_ln_allocation
  use negf_libnegf
  use negf_mat_def
  use negf_sparsekit_drv
  use negf_lib_param
  use negf_ln_extract

  implicit none

  Type(Tnegf), target :: negf
  Type(Tnegf), pointer :: pnegf

  pnegf => negf

  print*,'(main) testneqint'

  negf%file_re_H='H_real2.dat'
  negf%file_im_H='H_imm2.dat'
  negf%file_struct='driver'
  negf%verbose = 10
  negf%eneconv = HAR !to convert kb
  negf%isSid = .true.
  negf%form%formatted = .true.
  negf%form%type = 'PETSc'
  negf%form%fmt = 'F'
  negf%ReadOldSGF = 2

  print*,'(main) init'

  call init_negf(pnegf)

  print*,'(main) extract device'

  call extract_device(pnegf)

  print*,'(main) extract contacts'

  call extract_cont(pnegf)

  print*,'(main) contour over real axis'

  call real_axis_int(pnegf)

  write(*,*) 'trace=', real(trace(pnegf%rho))*negf%spin

  call destroy(pnegf%rho)

  print*,'(main) destroy negf'

  call destroy_negf(pnegf)

  call writepeakinfo(6)
  call writememinfo(6)

end program testneqint  
