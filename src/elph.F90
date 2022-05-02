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

module elph

  use ln_precision, only : dp
  use globals
  use ln_allocation
  use mat_def, only : create, destroy, z_DNS, z_CSR, z_COO
  use sparsekit_drv, only : prealloc_mult, prealloc_sum, coo2csr, extract, csr2dns

  implicit none
  private

  public :: Telph
  public :: interaction_models

  !> Describe the model implemented. Currently supported:
  !! 0 : dummy model, no electron-phonon interactions
  !! 1 : electron phonon dephasing limit (as in Datta, Cresti etc.)
  !!     Assumes elastic scattering and fully local (diagonal) model
  !!     Coupling is diagonal (Local Deformation Potential) and
  !!     Self energies are diagonal as well
  !! 2 : semi-local electron-phonon dephasing
  !!     Similar to 1, but the oscillator is considered local on
  !!     more than a contiguos basis function per oscillator position Ri
  !!     It is a atom-block generalization of 1 for LCAO
  !!     Coupling is diagonal (per oscillator site, per orbital)
  !!     Self energy is an array of atomic block
  !!     An additional descriptor with the number of orbitals per atom is
  !!     needed.
  !!     Note: the modes do not need to be local on a single atom, but
  !!     you need the orbitals on a given local phonon site to be contiguous
  !! 3 : as 2, but I will set the coupling as csr matrix including
  !!     overlap M->MS/2+SM/2 and treat each atomic oscillator mode separately
  !!     This is needed to verify whether neglecting overlap is ok
  type TModelsEnum
    integer :: dummy = 0
    integer :: dephdiagonal = 1
    integer :: dephatomblock = 2
    integer :: dephoverlap = 3
    integer :: polaroptical = 4
    integer :: nonpolaroptical = 5
    integer :: acousticinel = 6
    integer :: matrixcoupling = 7
  end type TModelsEnum

  type(TModelsEnum), parameter :: interaction_models = TModelsEnum()

  !> This type contains information describing different electron phonon
  !! models: input parameters, temporary data and output quantities
  !!
  !! Note: I don't use explicitely interfaces and different data types
  !! for different models because then we'll need all the different
  !! containers in negf container (cannot use virtualization)
  !! It may be managed in a more elegant way
  type Telph
    !! Common parameters
    integer :: model = interaction_models%dummy

    !> An array specifying how many orbital per atom (assumed contiguous)
    !  used in model semilocal
    integer, allocatable, dimension(:) :: orbsperatm

    logical :: Selfene_Gr
    logical :: Selfene_Gless
    logical :: Selfene_Hilb
    logical :: memory
    logical :: check

  end type Telph

end module elph
