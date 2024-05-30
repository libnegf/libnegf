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

module cublas_interface
  use iso_c_binding
  implicit none
  private

  public :: cublasSetMathMode

  integer, parameter, public ::  CUBLAS_TENSOR_OP_MATH = 1  ! DEPRECATED
  integer, parameter, public ::  CUBLAS_PEDANTIC_MATH = 2
  integer, parameter, public ::  CUBLAS_TF32_TENSOR_OP_MATH = 3
  integer, parameter, public ::  CUBLAS_DEFAULT_MATH = 0


  interface
  integer(4) function cublasSetMathMode(handle, mode) &
        & bind(c, name='cublasSetMathMode')
    use cublas_v2
    type(cublasHandle), value :: handle
    integer, intent(in), value :: mode
  end function cublasSetMathMode
  end interface
end module cublas_interface
