program Testdos

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


  pnegf => negf

  print*,'(main) testdos'

  negf%file_re_H='H_real.dat'
  negf%file_im_H='H_imm.dat'
  negf%file_re_S='S_real.dat'
  negf%file_im_S='S_imm.dat'
  negf%file_struct='driver'
  negf%verbose = 10
  negf%eneconv = 1.d0 ! HAR  ! to convert Kb 
  negf%isSid = .false.
  negf%form%formatted = .false.
  negf%form%type = 'UPT'
  negf%form%fmt = 'F'

  print*,'(main) init'

  call init_negf(pnegf)

  print*,'(main) computing dos'

  call compute_dos(pnegf)

  print*,'(main) destroy negf'

  call destroy_negf(pnegf)

  call writepeakinfo(6)
  call writememinfo(6)


end program Testdos


! negf:0 XXXXXXXXXXXX  negf:END
! pointer 
