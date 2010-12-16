program testcontactdos

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
  use extract

  implicit none

  Type(Tnegf), target :: negf
  Type(Tnegf), pointer :: pnegf


  pnegf => negf

  print*,'(main) testcontactdos'

  negf%file_re_H='H_real2.dat'
  negf%file_im_H='H_imm2.dat'
!  negf%file_re_S='S_real.dat'
!  negf%file_im_S='S_imm.dat'
  negf%file_struct='driver'
  negf%verbose = 10
  negf%eneconv = 1.d0 ! HAR  ! to convert Kb 
  negf%isSid = .true.
  negf%form%formatted = .true.
  negf%form%type = 'PETSc'
  negf%form%fmt = 'F'
  negf%ReadOldSGF = 2

  print*,'(main) init'

  call init_negf(pnegf)

  print*, '(main) extract device'

  call extract_device(pnegf)

  print*, '(main) extract contacts'

  call extract_cont(pnegf)
  

  print*,'(main) computing dos'

  call compute_dos(pnegf)

  print*,'(main) destroy negf'

  call destroy_negf(pnegf)

  call writepeakinfo(6)
  call writememinfo(6)


end program testcontactdos


! negf:0 XXXXXXXXXXXX  negf:END
! pointer 
