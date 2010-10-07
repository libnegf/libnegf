program readHS

  use mat_def
  use precision

  implicit none

  Type(z_CSR) :: zmat
  Type(r_COO) :: mat, mat1
  Type(z_coo) :: mat2

  Integer :: id = 10

  open(id, file='H_real.dat', form='formatted')

  call read_H(id,zmat)

 ! open(id, file='S_real.dat', form='formatted')

 !call read_S(id, mat)

  close(id)

end program readHS

