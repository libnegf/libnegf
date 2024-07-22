program test
  use ln_precision
  use mat_def
  use lib_param, only : cublasHandle
  use libnegf, only : cublasInitialize, cublasFinalize
  use cudautils
  implicit none

  !type TArray
  !  type(TMat), allocatable :: Arr(:,:)
  !end type TArray

  integer, parameter :: Npl = 5 
  integer, parameter :: Nr = 2000, Nc = 3000
  type(z_DNS), allocatable :: A(:,:)
  type(z_DNS) :: work1, work2, myMat, myMat_dag 
  complex(dp) :: alpha, beta
  real(dp), allocatable :: R(:,:)
  integer :: istat
  type(cublasHandle) :: hcublas
 
  call cublasInitialize(hcublas)

  allocate(A(Npl,Npl))
  print*,'create A11'
  call createAll(A(1,1), Nr, Nr)
  print*,'create A12'
  call createAll(A(1,2), Nr, Nc)

  allocate(R(Nr,Nr))
  call random_number(R)
  A(1,1)%val = R
  call random_number(R)
  A(1,1)%val = A(1,1)%val + (0.0_dp,1.0_dp)*R
  deallocate(R)
  
  allocate(R(Nr,Nc))
  call random_number(R)
  A(1,2)%val = R
  call random_number(R)
  A(1,2)%val = A(1,2)%val + (0.0_dp,1.0_dp)*R
  deallocate(R)

  !print*,'createGPU A11'
  !call createGPU(A(1,1))
  print*,'copyGPU A11'
  call copyToGPU(A(1,1))
  !print*,'createGPU A12'
  !call createGPU(A(1,2))
  print*,'copyGPU A12'
  call copyToGPU(A(1,2))

  print*,'create work'
  call createAll(work1, Nr, Nc)
  !print*,'createGPU work'
  !call createGPU(work1)
  print*,'matmul'
  alpha = cmplx(1.0,0.0,dp)
  beta = cmplx(0.0,0.0,dp)
  call matmul_gpu(hcublas,alpha,A(1,1),A(1,2),beta,work1)
  
  print*,'copy Back'
  call copyFromGPU(work1)

  call create(work2, Nr, Nc)

  work2%val = matmul(A(1,1)%val, A(1,2)%val)

  if (any(abs(work1%val-work2%val) > 1e-10)) then
    print*,'TEST failed'
  else    
    print*,'TEST Ok'
  end if

  !print*, work%val(1,1), work%val(2,1), work%val(3,1)

  print*,'finalize'
  !call deleteGPU(work1)
  !call deleteGPU(A(1,2))
 
  call destroyAll(work1)
  call destroyAll(A(1,2))
  call destroy(work2)

  print*, '' 
  print*, 'TEST dagger' 
  call createAll(work1, Nr, Nr)
  call dagger_gpu(hcublas, A(1,1), work1)
  call copyFromGPU(work1)

  call create(work2, Nr, Nr)
  work2%val = conjg(transpose(A(1,1)%val))

  if (any(abs(work1%val-work2%val) > 1e-10)) then
    print*,'TEST failed'
  else    
    print*,'TEST Ok'
  end if

  print*,'finalize'
 
  call destroyAll(work1)
  call destroyAll(A(1,1))
  call destroy(work2)

  !call create(myMat, 4, 4)
  !myMat%val(:,:) = (0.0_dp, 0.0_dp)
  !myMat%val(1,1) = (1.0_dp,1.0_dp)
  !myMat%val(1,2) = (0.0_dp,-1.0_dp)
  !myMat%val(1,3) = (0.0_dp,4.0_dp)
  !myMat%val(1,4) = (3.0_dp,0.0_dp)
  !print*, 'A:'
  !call print_myMat(myMat%val)
  !print*, ''

  !call copyToGPU(myMat)
  !call create(myMat_dag, 4, 4)
  !call createGPU(myMat_dag)
  
  !call copyFromGPU(myMat_dag)
  !call print_myMat(myMat_dag%val)

  !call deleteGPU(myMat)
  !call deleteGPU(myMat_dag)
  !call destroy(myMat)
  !call destroy(myMat_dag)

  call cublasFinalize(hcublas)

contains
    subroutine print_myMat(mat)
        complex(dp), intent(in) :: mat(:,:)
        integer :: nrow, i

        nrow = size(mat,1)
 
        do i = 1, nrow
           print*, mat(i,:)
           print*, '' 
        end do
    end subroutine print_myMat
end program test

    
    
