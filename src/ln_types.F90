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

module ln_types
  use ln_precision, only : dp, sp
  use mat_def, only : z_CSR
  implicit none
  private

  public ::  TEnGrid
  public :: intArray, complArray, realArray
  public :: TMatrixArrayHS, TMatrixArrayDM

  type intArray
    integer, dimension(:), allocatable :: indexes
  end type intArray

  type complArray
    complex(dp), dimension(:), allocatable :: array
  end type complArray

  type realArray
    real(dp), dimension(:), allocatable :: array
  end type realArray

  ! Array used to store the k-dependent matrices
  type TMatrixArrayHS
    type(z_CSR), pointer :: H => null()    ! Points to externally allocated H
    type(z_CSR), pointer :: S => null()
    integer    :: iS = 1                   ! spin-index (global)
    logical    :: isSid = .false.          ! True if overlap S == Id
    logical    :: internalHS = .true.      ! tells HS are internally allocated
  end type TMatrixArrayHS

  ! Array used to store the k-dependent matrices
  type TMatrixArrayDM
    type(z_CSR), pointer :: rho => null()  ! Holding output Matrix
    integer    :: iS = 1                   ! spin-index (global)
    logical    :: internalDM = .true.      ! tells DM is internally allocated
  end type TMatrixArrayDM

  !! Structure used to define energy points for the integration
  !! For every point we define
  !!     path (1,2 or 3): the energy point belongs to a real axis
  !!     integration (1), a complex plane integration (2) or a
  !!     pole summation (3)
  !!     pt_path: relative point number within a single path
  !!     pt: absolute point number along the whole integration path
  !!     cpu: cpu assigned to the calculation of the given energy point
  !!     Ec: energy value
  !!     wght: a weight used in final summation to evaluate integrals
  type TEnGrid
     integer :: path
     integer :: pt_path
     integer :: pt
     integer :: pt_cpu
     integer :: cpu
     complex(dp) :: Ec
     complex(dp) :: wght
  end type TEnGrid

  !-----------------------------------------------------------------------------

end module ln_types

