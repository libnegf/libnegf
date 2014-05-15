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


program Testdos

  use libnegf
  use lib_param
  use integrations

  implicit none

  Type(Tnegf), target :: negf
  Type(Tnegf), pointer :: pnegf

  pnegf => negf

  print*,'libnegf hello world'

  print*,'(main) init'

  call init_negf(pnegf)


  print*,'Ready to do stuff'

  call read_HS(pnegf)
  call read_negf_in(pnegf)

  call negf_partition_info(pnegf)

  !print*,'Try a basic call'

  !call compute_density_dft(pnegf)

  print*,'(main) destroy negf'

 call destroy_negf(pnegf)

  print*,'done'



end program Testdos


! negf:0 XXXXXXXXXXXX  negf:END
! pointer 
