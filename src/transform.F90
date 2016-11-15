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



module transform

   use ln_precision
   use mpi_globals
   use libmpifx_module, only : mpifx_comm
   use, intrinsic :: iso_c_binding
   implicit none   
   include 'fftw3.f03'
   private

   public hilbert, fourier, hilbert_shift   
   public convolve_retarded   
       

 contains

  subroutine Fourier(G)
  

   COMPLEX(dp), DIMENSION(:) :: G   

   !TYPE (DFTI_DESCRIPTOR), POINTER :: my_desc_1
   INTEGER :: N, m, stat, i
   REAL(dp) :: s

   N = size(G)

   !stat = DftiCreateDescriptor (my_desc_1, DFTI_DOUBLE, DFTI_COMPLEX, 1, N)
  
   !stat = DftiCommitDescriptor (my_desc_1)

   !stat = DftiComputeForward (my_desc_1, G)

   !stat = DftiFreeDescriptor (my_desc_1)
  
   G = G/N   

  end subroutine Fourier


  subroutine Hilbert(G)
  

   COMPLEX(dp), DIMENSION(:) :: G   

   !TYPE (DFTI_DESCRIPTOR), POINTER :: my_desc_1
   INTEGER :: N, m, stat, i
   REAL(dp) :: s

   N = size(G)

   !stat = DftiCreateDescriptor (my_desc_1, DFTI_DOUBLE, DFTI_COMPLEX, 1, N)
  
   !stat = DftiCommitDescriptor (my_desc_1)

   !stat = DftiComputeForward (my_desc_1, G)

   !stat = DftiFreeDescriptor (my_desc_1)
     
   G(1) = 0.0_dp
   do i = 0 , N/2-1
     G(i+1) = -(0.0_dp,1.0_dp) * G(i+1) / N
   end do
   G(N/2+1) = 0.0_dp
   do i = N/2+1 , N-1
     G(i+1) = (0.0_dp,1.0_dp) * G(i+1) / N
   end do



   !stat = DftiCreateDescriptor(my_desc_1, DFTI_DOUBLE, DFTI_COMPLEX, 1, N)
 
   !stat = DftiCommitDescriptor(my_desc_1)

   !stat = DftiComputeBackWard(my_desc_1, G)

   !stat = DftiFreeDescriptor (my_desc_1)

  end subroutine Hilbert


  subroutine Hilbert_shift(G,Wq)
  

   COMPLEX(dp), DIMENSION(:) :: G   
   REAL(dp) :: Wq
    

   !TYPE (DFTI_DESCRIPTOR), POINTER :: my_desc_1
   INTEGER :: N, m, stat, i
   REAL(dp) :: s

   N = size(G)

   !stat = DftiCreateDescriptor (my_desc_1, DFTI_DOUBLE, DFTI_COMPLEX, 1, N)
  
   !stat = DftiCommitDescriptor (my_desc_1)

   !stat = DftiComputeForward (my_desc_1, G)

   !stat = DftiFreeDescriptor (my_desc_1)
     


   G(1) = 0

   do i = 0 , N/2-1

     s = i * 1.0_dp

     G(i+1) = 2.0_dp * G(i+1) * sin(Wq*s) / N

   enddo

   G(N/2+1) = 0

   do i = N/2 + 1, N-1

     s = (i - N) * 1.0_dp

     G(i+1) = - 2.0_dp * G(i+1) * sin(Wq*s) / N

   end do
   

   !stat = DftiCreateDescriptor(my_desc_1, DFTI_DOUBLE, DFTI_COMPLEX, 1, N)
 
   !stat = DftiCommitDescriptor(my_desc_1)

   !stat = DftiComputeBackWard(my_desc_1, G)

   !stat = DftiFreeDescriptor (my_desc_1)

   
  end subroutine Hilbert_shift

  !****************************************************************************
  ! SUBROUTINE USED TO COMPUTE THE FFT 
  !
  !
  !****************************************************************************
  subroutine convolve_retarded(mympi,L,Ar,Br,C)
    type(mpifx_comm), intent(in) :: mympi
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:), allocatable, target :: Ar
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:), allocatable, target :: Br
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:), allocatable, target :: C
    integer, intent(in) :: L 

    ! local variables
    type(C_PTR) :: planF, planB
    complex(C_DOUBLE_COMPLEX), allocatable :: data1(:), data2(:), data3(:)
    complex(C_DOUBLE_COMPLEX), allocatable :: buffer(:)
    complex(C_DOUBLE_COMPLEX), pointer :: pAr(:,:), pBr(:,:), pC(:,:)
    integer :: local_L, i, j, k, m, nrow, ncol, err, numproc
 
    nrow = size(Ar,1)
    ncol = size(Ar,2)
    local_L = size(Ar,3)
    numproc = mympi%size

    allocate(buffer(2*Local_L), stat=err)
    allocate(data1(2*L), stat=err)
    allocate(data2(2*L), stat=err)
    allocate(data3(2*L), stat=err)
    if (err/=0) stop 'fftw_serial: allocation error'
 
    !   create plan for in-place forward/backward DFT 
    planF = fftw_plan_dft_1d(2*L, data1, data1, FFTW_FORWARD, FFTW_MEASURE)
    planB = fftw_plan_dft_1d(2*L, data3, data3, FFTW_BACKWARD, FFTW_MEASURE)
 
    ! RESHAPE in place !
    call c_f_pointer(C_LOC(Ar), pAr, [nrow*ncol,local_L])
    call c_f_pointer(C_LOC(Br), pBr, [nrow*ncol,local_L])
    call c_f_pointer(C_LOC(C), pC, [nrow*ncol,local_L])
 
 
    ! walk the matrix like a vector, chunking by numprocs 
    do j = 1, nrow*ncol, numproc
      do i = 0, numproc-1
         if (j+i > nrow*ncol) exit
         do m = 1, local_L
            buffer(2*m-1)=  pAr(j+i, m)
            buffer(2*m)  = -conjg(pAr(j+i, m))
         end do
         call mpifx_gather(mympi,buffer,data3,i)
         do m = 1, local_L
            data1(m) = data3(2*m-1)
            data1(2*L-m) = data3(2*m)
         end do
         do m = 1, local_L
            buffer(2*m-1)=  pBr(j+i, m)
            buffer(2*m)  = -conjg(pBr(j+i, m))
         end do
         call mpifx_gather(mympi,buffer,data3,i)
         do m = 1, local_L
            data2(m) = data3(2*m-1)
            data2(2*L-m) = data3(2*m)
         end do
      end do
 
      call fftw_execute_dft(planF, data1, data1)
      call fftw_execute_dft(planF, data2, data2)
 
      ! no dot_product: complex arrays
      do k = 1, 2*L
         data3(k) = data1(k)*data2(k)
      end do
 
      call fftw_execute_dft(planB, data3, data3)
 
      data3 = data3/(2*L)
 
      do i = 0, numproc-1
        if (j+i > nrow*ncol) exit
        call mpifx_scatter(mympi,data3,buffer,i)
        pC(j+i,1:local_L) = buffer(1:local_L)
      end do
    end do
 
    call fftw_destroy_plan(planF)
    call fftw_destroy_plan(planB)
    deallocate(data1,data2,data3)
    deallocate(buffer)

 end subroutine convolve_retarded

end module transform       
