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


program Testiterativedns

  use negf_ln_precision
  use negf_ln_constants
  use negf_mat_def
  use negf_ln_allocation
  use negf_ln_structure
  use negf_input_output
  use negf_sparsekit_drv
  use inversions
  use negf_iterative
  use negf_libnegf
  use negf_lib_param
  use negf_clock

  implicit none

  Type(Tnegf), target :: negf
  Type(Tnegf), pointer :: pnegf
  Integer :: outer, N, k, i
  Complex(dp) :: Ec
  Type(z_CSR) :: SelfEneR(1), Gr
  Type(z_CSR) :: Tlc(1), Tcl(1), gsurfR(1)


  pnegf => negf

  print*,'(main) testiterativedns'

  negf%file_re_H='H_real.dat'
  negf%file_im_H='H_imm.dat'
  negf%file_re_S='S_real.dat'
  negf%file_im_S='S_imm.dat'
  negf%file_struct='driver'
  negf%verbose = 10
  negf%eneconv = 1.d0  !HAR  ! to convert Kb
  negf%isSid = .false.
  negf%form%formatted = .false.
  negf%form%type = 'UPT'
  negf%form%fmt = 'F'


  outer = 0

  print*,'(main) init'

  call init_negf(pnegf)

  print*,'(main) computing dos'


    N = nint((negf%Emax-negf%Emin)/negf%Estep)

    print*,'N=',N
    print*,'delta=',negf%delta

print*,'(main) create SelfEner, Tlc, Tcl, gsurfR'

    call create(SelfEneR(1),0,0,0)
    call create(Tlc(1),0,0,0)
    call create(Tcl(1),0,0,0)
    call create(gsurfR(1),0,0,0)

print*,'(main) open dos.dat'

    open(101,file='dos.dat')

    do i=1,N

       Ec=(negf%Emin+i*negf%Estep)+negf%delta*(0.d0,1.d0)

       call message_clock('Compute Green`s function')

       call calculate_Gr(negf%H,negf%S,Ec,SelfEneR,Tlc,Tcl,gsurfR,Gr,negf%str,outer)

       call write_clock

       negf%dos = -aimag( trace(Gr) )/pi

       write(101,*) real(Ec), negf%dos
       print*, real(Ec), negf%dos, Gr%nnz

       call destroy(Gr)

    enddo

    call destroy(SelfEneR(1))
    call destroy(Tcl(1))
    call destroy(Tlc(1))
    call destroy(gsurfR(1))

    close(101)


  print*,'(main) destroy negf'

  call destroy_negf(pnegf)

  call writepeakinfo(6)
  call writememinfo(6)


end program Testiterativedns


! negf:0 XXXXXXXXXXXX  negf:END
! pointer
