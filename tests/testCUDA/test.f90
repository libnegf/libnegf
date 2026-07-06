program test
  use ln_precision
  use ln_allocation
  use mat_def
  use lib_param, only : cublasHandle, cusolverDnHandle
  use libnegf, only : cublasInitialize, cublasFinalize
  use cudautils
  implicit none

  !type TArray
  !  type(TMat), allocatable :: Arr(:,:)
  !end type TArray

  integer, parameter :: Npl = 5 
  integer, parameter :: Nr = 5000, Nc = 300
  type(z_DNS), allocatable :: A(:,:)
  type(z_DNS) :: work1, work2, work3
  complex(dp) :: alpha, beta
  real(dp), allocatable :: R(:,:)
  real(dp) :: tr1, tr2
  integer :: istat, ii, ndevs
  integer(8) :: totalmem, freemem
  type(cublasHandle) :: hblas
  type(cusolverDnHandle) :: hsolver
 
  call getDeviceCount(ndevs)
  print*,'found ',ndevs,' GPUs'
  call getDevice(ii)
  print*,'current device number: ',ii
  
  call getDevMemInfo(freemem, totalmem)
  print*, 'Total device memory: ',totalmem
  print*, ' Free device memory: ',freemem 
  
  call cublasInitialize(hblas)
  call cusolverInitialize(hsolver)

  allocate(A(Npl,Npl))
  print*,'create A11'
  call createAll(A(1,1), Nr, Nr)
  print*,'create A12'
  call createAll(A(1,2), Nr, Nc)
  print*,'create A22'
  call createAll(A(2,2), Nr, Nr)

  allocate(R(Nr,Nr))
  call random_number(R)
  A(1,1)%val = R
  A(2,2)%val = R
  call random_number(R)
  A(1,1)%val = A(1,1)%val + (0.0_dp,1.0_dp)*R
  A(2,2)%val = A(2,2)%val + (0.0_dp,1.0_dp)*R
  deallocate(R)
  print*,'maxval abs(A11)',maxval(abs(A(1,1)%val))


  allocate(R(Nr,Nc))
  call random_number(R)
  A(1,2)%val = R
  call random_number(R)
  A(1,2)%val = A(1,2)%val + (0.0_dp,1.0_dp)*R
  deallocate(R)

  !/////////////////////////////////////////////////////////////////////
  istat = 0

  print*, '' 
  print*,'1.1 TEST matmul GPU: C = A*B'
  print*,'copy A11 to GPU'
  call copyToGPU(A(1,1))
  print*,'copyGPU A12'
  call copyToGPU(A(1,2))

  print*,'create work'
  call createAll(work1, Nr, Nc)
  alpha = cmplx(1.0,0.0,dp)
  beta = cmplx(0.0,0.0,dp)
  call matmul_gpu(hblas,alpha,A(1,1),A(1,2),beta,work1)
  
  print*,'copy Back'
  call copyFromGPU(work1)

  print*,'matmul CPU'
  call create(work2, Nr, Nc)
  work2%val = matmul(A(1,1)%val, A(1,2)%val)

  print*,'TEST check'
  if (any(abs(work1%val-work2%val) > 1e-10)) then
    print*,'TEST failed'
    istat = 1001
  else    
    print*,'TEST Ok'
  end if
  print*,'-------------------------'


  call createAll(A(2,1), Nc, Nr)
  A(2,1)%val = conjg(transpose(A(1,2)%val))
  call destroyAll(A(1,2))
  call copyToGPU(A(2,1))

  print*,'1.1 TEST matmul GPU: C = A*B^H'
  call matmul_gpu(hblas,alpha,A(1,1),A(2,1),beta,work1,'dag_2nd')
  
  print*,'copy Back'
  call copyFromGPU(work1)

  print*,'TEST check'
  if (any(abs(work1%val-work2%val) > 1e-10)) then
    print*,'TEST failed'
    istat = 10021 
  else    
    print*,'TEST Ok'
  end if

  print*,'1.2 TEST matmul GPU: C = C + A*B^H'
  beta = cmplx(1.0,0.0,dp)
  call matmul_gpu(hblas,alpha,A(1,1),A(2,1),beta,work1,'dag_2nd')
  
  print*,'copy Back'
  call copyFromGPU(work1)

  print*,'TEST check'
  if (any(abs(work1%val-2.0*work2%val) > 1e-10)) then
    print*,'TEST failed'
    istat = 10022 
  else    
    print*,'TEST Ok'
  end if

  call destroyAll(A(2,1))
  call destroyAll(work1)
  call destroy(work2)
  print*,'-------------------------'
  
  !/////////////////////////////////////////////////////////////////////

  print*,'3. TEST trace'
  call trace_gpu(hblas, A(1,1), tr1)
  print*,'TEST check'
  tr2 = real(0.0,dp)
  do ii = 1, Nr
    tr2 = tr2 + real(A(1,1)%val(ii,ii))
  end do
  if (abs(tr1-tr2) > 1e-10) then
    print*,'TEST failed', tr1, tr2
    istat = 1003 
  else    
    print*,'TEST Ok'
  end if
  print*,'-------------------------'

  print*
  print*,'4. TEST hardcopy mat '
  print*,'copy A22 to GPU'
  call copyToGPU(A(2,2))
  
  print*,'create work'
  call createAll(work1, Nr, Nr)

  !/////////////////////////////////////////////////////////////////////
  call copy_mat_gpu(hblas, A(2,2), work1)
  print*,'copy Back'
  call copyFromGPU(work1)
  print*,'TEST check'
  if (any(abs(work1%val-A(2,2)%val) > 1e-13)) then
    print*,'TEST failed'
    istat = 1004
  else    
    print*,'TEST Ok'
  end if
  print*,'-------------------------'

  !/////////////////////////////////////////////////////////////////////
  print*, '' 
  print*, '5. TEST dagger GPU' 
  call dagger_gpu(hblas, A(1,1), work1)
  print*,'copy Back'
  call copyFromGPU(work1)

  print*, 'dagger CPU' 
  call create(work2, Nr, Nr)
  work2%val = conjg(transpose(A(1,1)%val))

  print*,'TEST check'
  if (any(abs(work1%val-work2%val) > 1e-10)) then
    print*,'TEST failed'
    istat = 1005
  else    
    print*,'TEST Ok'
  end if
  print*,'finalize'
  call destroy(work2)
  print*,'-------------------------'

  !/////////////////////////////////////////////////////////////////////
  !print*
  !print*,'6. TEST kernelsum GPU'
  !alpha = cmplx(1.0,0.0,dp)
  !beta = cmplx(1.0,0.0,dp)
  !call kernelsum_gpu(work1, alpha, A(1,1), beta, A(2,2))
  !print*,'copy Back'
  !call copyFromGPU(work1)
  !print*,'TEST check'
  !if (any(abs(work1%val-A(1,1)%val-A(2,2)%val) > 1e-10)) then
  !  print*,'TEST failed'
  !  istat = 1006
  !else    
  !  print*,'TEST Ok'
  !end if
  !print*,'-------------------------'

  !/////////////////////////////////////////////////////////////////////
  print*
  print*,'7. TEST matsum GPU C = A + B'
  alpha = cmplx(1.0,0.0,dp)
  beta = cmplx(1.0,0.0,dp)
  call matsum_gpu(hblas, alpha, A(1,1), beta, A(2,2), work1) 
  print*,'copy Back'
  call copyFromGPU(work1)

  print*,'TEST check'
  if (any(abs(work1%val-A(1,1)%val-A(2,2)%val) > 1e-10)) then
    print*,'TEST failed'
    istat = 1007
  else    
    print*,'TEST Ok'
  end if
  print*,'-------------------------'
  
  !/////////////////////////////////////////////////////////////////////
  print*
  print*,'7b. TEST matsum in-place GPU C = C + B'
  call matsum_gpu(hblas, alpha, work1, beta, A(2,2), work1) 
  print*,'copy Back'
  call copyFromGPU(work1)
  print*,'TEST check'
  if (any(abs(work1%val-A(1,1)%val-2.0_dp*A(2,2)%val) > 1e-10)) then
    print*,'TEST failed'
    istat = 1007
  else    
    print*,'TEST Ok'
  end if
  call destroyAll(A(2,2))
  print*,'-------------------------'


  !/////////////////////////////////////////////////////////////////////
  print*, '8. TEST spectral as C = i[A - A^H]' 
  alpha = cmplx(0.0,1.0,dp)
  beta = cmplx(0.0,-1.0,dp)
  call matsum_gpu(hblas, alpha, A(1,1), beta, A(1,1), work1,'dag_2nd') 
  print*,'copy Back'
  call copyFromGPU(work1)
  print*,'TEST check'
  if (any(abs(A(1,1)%val-conjg(transpose(A(1,1)%val))+alpha*work1%val) > 1e-10)) then
    print*,'TEST failed'
    istat = 1008
  else    
    print*,'TEST Ok'
  end if
  print*,'finalize'
  print*,'-------------------------'


  print*
  print*,'9. TEST spectral GPU'
  call spectral_gpu(hblas, A(1,1), work1) 

  print*,'copy Back'
  call copyFromGPU(work1)

  print*,'TEST check'
  beta = cmplx(0.0,1.0,dp)
  if (any(abs(A(1,1)%val-conjg(transpose(A(1,1)%val))+beta*work1%val) > 1e-10)) then
    print*,'TEST failed'
    istat = 1009
  else    
    print*,'TEST Ok'
  end if

  print*,'finalize'
  call destroyAll(work1)
  call destroyAll(A(1,1))
  print*,'-------------------------'
 
  !/////////////////////////////////////////////////////////////////////
  print*, '' 
  print*, '10. TEST inverse GPU' 
  print*,'create A11'
  call createAll(A(1,1), Nr, Nr)
  allocate(R(Nr,Nr))
  call random_number(R)
  A(1,1)%val = R
  call random_number(R)
  A(1,1)%val = A(1,1)%val + (0.0_dp,1.0_dp)*R
  deallocate(R)
  print*,'maxval abs(A11)',maxval(abs(A(1,1)%val))

  call copyToGPU(A(1,1))
  call createAll(work1, Nr, Nr)

  call inverse_gpu(hblas, hsolver, A(1,1), work1, istat)
  if (istat > 0) then
     print*, 'Error in inverse:',istat
     istat = 1010
  end if

  print*,'TEST A*A^-1 = Id'
  call createAll(work2, Nr, Nr)
  alpha = cmplx(1.0,0.0,dp)
  beta = cmplx(0.0,0.0,dp)
  call matmul_gpu(hblas,alpha,A(1,1),work1,beta,work2)
  print*,'copy Back'
  call copyFromGPU(work2)

  ! Set an identity matrix
  work1%val = (0.0_dp, 0.0_dp)
  do ii = 1, Nr
     work1%val(ii,ii) = (1.0_dp, 0.0_dp)
  end do  
 
  print*,'check'
  if (any(abs(work1%val-work2%val) > 1e-10)) then
    print*,'TEST failed'
    istat = 1010
  else    
    print*,'TEST Ok'
  end if

  print*,'finalize'
  call destroyAll(work1)
  call destroyAll(A(1,1))
  call destroyAll(work2)
  print*,'-------------------------'
  print*,"Memory summary:"
  call writePeakInfo(6)
  call writeMemInfo(6)
  
  call cublasFinalize(hblas)
  call cusolverFinalize(hsolver)

  if (istat /= 0) then
    error stop "Error found in one of the HIP tests"
  end if

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

    
    
