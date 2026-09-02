!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Greens functions.         !
!! Copyright (C) 2012 - 2026                                                !
!!                                                                          !
!! This file is part of libNEGF: a library for                              !
!! Non Equilibrium Green's Functions calculations                           !
!!                                                                          !
!! Developers: Alessandro Pecchia, Daniele Soccodato                        !
!! Former Contributors: Gabriele Penazzi, Luca Latessa, Aldo Di Carlo       !
!!                                                                          !
!! libNEGF is free software: you can redistribute and/or modify it          !
!! under the terms of the GNU Lesser General Public License as published    !
!! by the Free Software Foundation, either version 3 of the License, or     !
!! (at your option) any later version.                                      !
!!                                                                          !
!!  You should have received a copy of the GNU Lesser General Public        !
!!  License along with libNEGF.  If not, see                                !
!!  <http://www.gnu.org/licenses/>.                                         !
!!--------------------------------------------------------------------------!


module ln_math
  use ln_precision

  public :: cross_product
  public :: volume

contains

  subroutine cross_product(vec1, vec2, vec3)
    real(dp), dimension(3), intent(in) :: vec1, vec2
    real(dp), dimension(3), intent(out) :: vec3

    vec3 = 0.0_dp

    vec3(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
    vec3(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
    vec3(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)

  end subroutine cross_product

  function volume(vec1, vec2, vec3) result (vol)
    real(dp), dimension(3), intent(in) :: vec1, vec2, vec3
    real(dp) :: vol

    real(dp), dimension(3) :: cross_prod

    call cross_product(vec2, vec3, cross_prod)
    vol = dot_product(vec1, cross_prod)

  end function volume

end module ln_math
