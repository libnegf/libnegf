program Testint

  use precision
  use constants
  use mat_def
  use allocation
  use structure
  use input_output
  use sparsekit_drv
  use inversions
  use iterative
  use libnegf
  use lib_param

  implicit none

  Type(Tnegf), target :: negf
  Type(Tnegf), pointer :: pnegf


  pnegf => negf

  print*,'(main) testdos'

  negf%file_re_H='H_real.dat'
  negf%file_im_H='H_imm.dat'
  negf%file_struct='driver'
  negf%verbose = 10
  negf%eneconv = HAR  ! to convert Kb 
  negf%isSid = .true.  

  print*,'(main) init'

  call init_negf(pnegf)

  print*,'(main) computing dos'

  call compute_dos(pnegf)

  print*,'(main) destroy negf'

  call destroy_negf(pnegf)

  call writememinfo(6)


end program Testint


! negf:0 XXXXXXXXXXXX  negf:END
! pointer 
