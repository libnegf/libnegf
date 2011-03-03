module ln_precision

  integer, parameter :: sp = selected_real_kind(6,30)
  integer, parameter :: dp = selected_real_kind(14,100)

  real(dp), parameter :: EPS=1.d-10

  interface
    real(8) function DLAMCH(C)
      character C
    end function DLAMCH
  end interface
  
contains 

  function get_machine_prec() result(racc)
    real(dp) :: racc
    racc=DLAMCH('Precision')
  end function get_machine_prec

end module ln_precision
