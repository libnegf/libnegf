program Testint

  use ln_precision
  use ln_constants
  use ln_allocation
  use libnegf
  use mat_def
  use sparsekit_drv
  use lib_param
  use ln_extract

  implicit none

  Type(Tnegf), target :: negf
  Type(Tnegf), pointer :: pnegf


  pnegf => negf

  print*,'(main) testint'

  negf%file_re_H='H_real2.dat'
  negf%file_im_H='H_imm2.dat'
!  negf%file_re_S='S_real2.dat'
!  negf%file_im_S='S_imm2.dat'
  negf%file_struct='driver'
  negf%verbose = 80
  negf%eneconv = HAR  ! to convert Kb 
  negf%isSid = .true.
  negf%form%formatted = .true.
  negf%form%type = 'PETSc'
  negf%form%fmt = 'F'
  negf%ReadOldSGF = 2
  negf%DorE = 'D'

  print*,'(main) init'

  call init_negf(pnegf)

  print*, '(main) extract device'

  call extract_device(pnegf)

  print*, '(main) extract contacts'

  call extract_cont(pnegf)

  print*,'(main) contour electrons'

  call contour_int(pnegf)

  write(*,*) 'trace=', real(trace(pnegf%rho))

  !call destroy(pnegf%rho)

  print*,'(main) real axis integration'

  call real_axis_int(pnegf)

  write(*,*) 'trace=', real(trace(pnegf%rho))

!  print*,'(main) contour holes'

!  call contour_int_p(pnegf)

!  write(*,*) 'trace=', real(trace(pnegf%rho))*negf%spin

  print*,'(main) destroy negf'

  call destroy_negf(pnegf)

  call writepeakinfo(6)
  call writememinfo(6)


end program Testint


! negf:0 XXXXXXXXXXXX  negf:END
! pointer 
