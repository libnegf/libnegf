program Testdos

  use precision
  use constants
  use mat_def
  use allocation
  use input_output
  use sparsekit_drv
  use inversions

  implicit none

  Type(z_CSR) :: zmat, Id, ESH
  Type(z_DNS) :: Inv
  real(dp) :: delta
  Integer :: i, k,  N
  real(dp) :: Emin, Emax, Estep, dos
  Complex(dp) :: E
  !Complex(dp), dimension(:,:), allocatable :: Inv

  open(101, file='H_real.dat', form='formatted')
  open(102, file='H_imm.dat', form='formatted')   !open imaginary part of H

  call read_H(101,102,zmat)

  close(101)
  close(102)

  Emin = -2.d0
  Emax = 2.d0
  Estep = 1.d-3

  delta = 3.d-3

print*, Emin,Emax

  N = nint((Emax-Emin)/Estep)
 
print*, 'N=',N

  call create(Id,zmat%nrow)
  call create(Inv,zmat%nrow,zmat%nrow)

  open(101,file='dos.dat')

  do i=1,N
     
     E=(Emin+i*Estep)+delta*(0.d0,1.d0)

     !call zcreate_id_CSR(Id, zmat%nrow)

     call prealloc_sum(Id,zmat,E,(-1.d0,0.d0),ESH)


     call zINV_PARDISO(ESH,zmat%nrow,Inv%val)
     

     dos = 0.d0
     do k=1,zmat%nrow
        dos = dos - aimag(Inv%val(k,k))/pi
     enddo

     write(101,*) real(E), dos
     write(*,*) real(E), dos

     call destroy(ESH)
     
  enddo
  
  close(101)

  call destroy(Id)
  call destroy(Inv)


end program Testdos
