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
