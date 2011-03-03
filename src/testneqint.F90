program testneqint

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
