program testcomplexband

  use ln_precision
  use ln_constants
  use mat_def
  use ln_allocation
  use ln_structure
  use input_output
  use libnegf
  use lib_param
  use clock

  implicit none

  Type(Tnegf), target :: negf
  Type(Tnegf), pointer :: pnegf

  Type(z_CSR), target :: H,S

  pnegf => negf

  print*,'(main) testcomplexband'
  call set_defaults(negf)


  negf%file_re_H='H_real.dat'
  negf%file_im_H='H_imm.dat'
!  negf%file_re_S='S_real.dat'
!  negf%file_im_S='S_imm.dat'
  negf%file_struct='driver'
  negf%verbose = 10
  negf%eneconv = HAR !1.d0 HAR  ! to convert Kb 
  negf%isSid = .true.
  negf%form%formatted = .true.
  negf%form%type = 'UPT'
  negf%form%fmt = 'U'
  negf%ReadOldSGF = 2

  pnegf%H => H
  pnegf%S => S

  print*,'(main) init'

  call init_negf(pnegf)

  print*,'(main) computing cmplx bands'

  call comp_complex_bands(pnegf)

  print*,'(main) destroy negf'

  call destroy_negf(pnegf)

  call writepeakinfo(6)
  call writememinfo(6)

end program testcomplexband
