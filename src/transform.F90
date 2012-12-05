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
   use MKL_DFTI
   implicit none   
   private

   public hilbert, fourier, hilbert_shift   
       

contains

subroutine Fourier(G)
  

   COMPLEX(dp), DIMENSION(:) :: G   

   TYPE (DFTI_DESCRIPTOR), POINTER :: my_desc_1
   INTEGER :: N, m, stat, i
   REAL(dp) :: s

   N = size(G)

   stat = DftiCreateDescriptor (my_desc_1, DFTI_DOUBLE, DFTI_COMPLEX, 1, N)
  
   stat = DftiCommitDescriptor (my_desc_1)

   stat = DftiComputeForward (my_desc_1, G)

   stat = DftiFreeDescriptor (my_desc_1)
  
   G = G/N   

end subroutine Fourier


subroutine Hilbert(G)
  

   COMPLEX(dp), DIMENSION(:) :: G   

   TYPE (DFTI_DESCRIPTOR), POINTER :: my_desc_1
   INTEGER :: N, m, stat, i
   REAL(dp) :: s

   N = size(G)

   stat = DftiCreateDescriptor (my_desc_1, DFTI_DOUBLE, DFTI_COMPLEX, 1, N)
  
   stat = DftiCommitDescriptor (my_desc_1)

   stat = DftiComputeForward (my_desc_1, G)

   stat = DftiFreeDescriptor (my_desc_1)
     
   G(1) = 0.0_dp
   do i = 0 , N/2-1
     G(i+1) = -(0.0_dp,1.0_dp) * G(i+1) / N
   end do
   G(N/2+1) = 0.0_dp
   do i = N/2+1 , N-1
     G(i+1) = (0.0_dp,1.0_dp) * G(i+1) / N
   end do



   stat = DftiCreateDescriptor(my_desc_1, DFTI_DOUBLE, DFTI_COMPLEX, 1, N)
 
   stat = DftiCommitDescriptor(my_desc_1)

   stat = DftiComputeBackWard(my_desc_1, G)

   stat = DftiFreeDescriptor (my_desc_1)

end subroutine Hilbert


subroutine Hilbert_shift(G,Wq)
  

   COMPLEX(dp), DIMENSION(:) :: G   
   REAL(dp) :: Wq
    

   TYPE (DFTI_DESCRIPTOR), POINTER :: my_desc_1
   INTEGER :: N, m, stat, i
   REAL(dp) :: s

   N = size(G)

   stat = DftiCreateDescriptor (my_desc_1, DFTI_DOUBLE, DFTI_COMPLEX, 1, N)
  
   stat = DftiCommitDescriptor (my_desc_1)

   stat = DftiComputeForward (my_desc_1, G)

   stat = DftiFreeDescriptor (my_desc_1)
     


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
   

   stat = DftiCreateDescriptor(my_desc_1, DFTI_DOUBLE, DFTI_COMPLEX, 1, N)
 
   stat = DftiCommitDescriptor(my_desc_1)

   stat = DftiComputeBackWard(my_desc_1, G)

   stat = DftiFreeDescriptor (my_desc_1)

   
end subroutine Hilbert_shift


end module transform       
