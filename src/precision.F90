module precision

  integer, parameter :: sp = selected_real_kind(6,30)
  integer, parameter :: dp = selected_real_kind(14,100)

  real(kind=dp), parameter :: EPS=1.d-10

  interface
    real(dp) function DLAMCH(C)
      character C
    end function DLAMCH
  end interface
  
contains 

  function get_machine_prec() result(racc)
    real(dp) :: racc
    racc=DLAMCH('Precision')
  end function get_machine_prec

end module precision
