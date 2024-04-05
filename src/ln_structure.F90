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


module ln_structure
use ln_precision
implicit none
private

public :: TStruct_info, create_TStruct, kill_TStruct
public :: TBasisCenters, create_TBasis, destroy_TBasis
public :: TNeighbourMap

type TStruct_Info
   integer, dimension(:), allocatable :: mat_PL_start  !ind(..)
   integer, dimension(:), allocatable :: mat_PL_end    !ind(..)

   integer, dimension(:), allocatable :: mat_B_start  ! starting bound
   integer, dimension(:), allocatable :: mat_C_start  ! starting real contact
   integer, dimension(:), allocatable :: mat_C_end    ! end real contact
   integer, dimension(:), allocatable :: cont_dim     ! total contact dim

   integer, dimension(:), allocatable :: cblk         !contact int block
   integer :: central_dim = 0
   integer :: total_dim = 0
   integer :: num_PLs = 0
   integer :: num_conts = 0
   integer :: active_cont = 0
   contains
   procedure :: print_Tstruct
end type TStruct_Info


!> Type to store geometry information
!> It is needed for elph coupling but also local currents

type TBasisCenters
  real(dp), dimension(:,:), allocatable :: x
  integer :: nCentralAtoms
  integer, dimension(:), allocatable :: matrixToBasis
  real(dp) :: lattVecs(3,3) = 0.0_dp
  integer :: transportDirection = 3
end type TBasisCenters

!> type to store neighbour maps between basis centers
type TNeighbourMap
   integer, DIMENSION(:), ALLOCATABLE :: nn
end type TNeighbourMap

contains

  ! INIT TStructure_Info:
  ! This subroutine requires the definition of
  ! ncont:                   number of contacts
  ! nbl:                     number of PLs
  ! PL_end(nbl)              array of size nbl with the end of each block
  ! surf_start(ncont)        array containing the start of each contact surface
  ! surf_end(ncont)          array containing the end of each contact surface
  ! cont_end(ncont)          array containing the end of each contact
  ! cblk(ncont)              array containing the PL-contact position
  !
  subroutine create_TStruct(ncont, nbl, PL_end, surf_start, surf_end, cont_end, cblk, str)
    integer, intent(in) :: ncont
    integer, intent(in) :: nbl
    integer, dimension(:), intent(in) :: PL_end, surf_start, surf_end, cont_end, cblk
    type(TStruct_Info), intent(inout) :: str


    integer :: i

    str%num_conts = ncont
    str%num_PLs = nbl
    str%active_cont = 1

    !  -------------------------------- - - -
    !  | S       |  PL1   |    PL2    |
    !  -------------------------------- - - -
    !  B_start   C_start              C_end
    !
    !  B_end = C_start - 1


    if(ncont.gt.0) then
      if (allocated(str%mat_B_start)) then
        deallocate(str%mat_B_start)
        deallocate(str%mat_C_start)
        deallocate(str%mat_C_end)
        deallocate(str%cont_dim)
      end if
      allocate(str%mat_B_start(ncont))
      allocate(str%mat_C_start(ncont))
      allocate(str%mat_C_end(ncont))
      allocate(str%cont_dim(ncont))
    endif
    if (allocated(str%cblk)) then
      deallocate(str%cblk)
    end if
    allocate(str%cblk(ncont))

    do i=1,ncont
       str%mat_B_start(i) = surf_start(i)
       str%mat_C_start(i) = surf_end(i) + 1
       str%mat_C_end(i)   = cont_end(i)
       str%cblk(i) = cblk(i)
       str%cont_dim(i) =  str%mat_C_end(i) - str%mat_B_start(i) + 1
    enddo

    if (allocated(str%mat_PL_start)) then
      deallocate(str%mat_PL_start)
      deallocate(str%mat_PL_end)
    end if
    allocate(str%mat_PL_start(nbl+1))
    allocate(str%mat_PL_end(nbl))

    str%mat_PL_start(1) = 1
    str%mat_PL_end(1) = PL_end(1)

    do i=2,nbl
       str%mat_PL_start(i) = PL_end(i-1)+1
       str%mat_PL_end(i) = PL_end(i)
    enddo
    str%mat_PL_start(nbl+1)= PL_end(nbl)+1

    str%central_dim = PL_end(nbl)

    if(ncont.gt.0) then
       str%total_dim = cont_end(ncont)
    else
       str%total_dim = str%central_dim
    endif

  end subroutine create_TStruct

  ! --------------------------------------------------------
  subroutine kill_TStruct(str)
    type(TStruct_Info) :: str

    if(allocated(str%mat_B_start)) deallocate(str%mat_B_start)
    if(allocated(str%mat_C_start)) deallocate(str%mat_C_start)
    if(allocated(str%mat_C_end))   deallocate(str%mat_C_end)
    if(allocated(str%cblk))  deallocate(str%cblk)

    if(allocated(str%mat_PL_start)) deallocate(str%mat_PL_start)
    if(allocated(str%mat_PL_end)) deallocate(str%mat_PL_end)
    if(allocated(str%cont_dim)) deallocate(str%cont_dim)
  end subroutine kill_TStruct

  ! --------------------------------------------------------
  subroutine print_Tstruct(str,io)
    class(Tstruct_Info) :: str
    integer, intent(in) :: io

    integer :: i

    write(io,*) 'Hamiltonian Structure:'
    write(io,*) 'num contacts:',str%num_conts
    write(io,*) 'num layers:',str%num_Pls

    do i=1,str%num_PLs
       write(io,*) 'PL:',i, str%mat_PL_start(i), str%mat_PL_end(i)
    enddo

    write(io,*) 'central dim:',str%central_dim

    do i=1,str%num_conts
      write(io,*) 'S: ',i, str%mat_B_start(i), str%mat_C_start(i)-1
      write(io,*) 'C: ',i, str%mat_C_start(i), str%mat_C_end(i)
      write(io,*) 'cont dim:',str%mat_C_end(i)-str%mat_C_start(i)+1
    enddo
    write(io,*) 'cblk: ',str%cblk
    write(io,*) 'total dim:',str%total_dim

  end subroutine print_Tstruct

  !--------------------------------------------------------------------
  ! Initialize TBasisCenters
  subroutine create_TBasis(this, coord, nCentral, lattVecs, basisToMatrix, matrixToBasis, transportDirection)
    type(TBasisCenters) :: this
    real(dp), intent(in) :: coord(:,:)
    integer, intent(in) :: nCentral
    ! lattice vectors stored column-wise, e.g., a1 = latVecs(:,1) 
    real(dp), intent(in), optional :: lattVecs(3,3)
    ! basis index to matrix index
    integer, intent(in), optional :: basisToMatrix(:)
    ! matrix index to basis index
    integer, intent(in), optional :: matrixToBasis(:)
    ! transport direction for coordinate systems different from NEGF
    integer, intent(in), optional :: transportDirection

    integer :: ii, jj, nn

    if (present(basisToMatrix) .and. present(matrixToBasis)) then
      stop 'Internal error in create_TBasis: two exclusive optional arrays'
    end if

    allocate(this%x(size(coord,1), size(coord,2)))
    this%x = coord
    this%nCentralAtoms = nCentral

    if (present(lattVecs)) then
       this%lattVecs = lattVecs
    end if

    if (present(matrixToBasis)) then
      allocate(this%matrixToBasis(size(matrixToBasis)))
      this%matrixToBasis = matrixToBasis
    end if
    !  basisToMatrix(1) = ind(1) = 1 !sp3
    !  basisToMatrix(2) = ind(2) = 5 !sp3
    !  basisToMatrix(3) = ind(3) = 9 !...
    ! => matrixToBasis(ind(ii):ind(ii+1)-1) = ii
    if (present(basisToMatrix)) then
      nn = size(basisToMatrix)
      allocate(this%matrixToBasis(basisToMatrix(nn)-1))
      do ii = 1, nn-1
        this%matrixToBasis(basisToMatrix(ii):basisToMatrix(ii+1)-1) = ii
      end do
    end if

    if (present(transportDirection)) then
      this%transportDirection = transportDirection
    end if

  end subroutine create_TBasis

  !--------------------------------------------------------------------
  subroutine destroy_TBasis(this)
    type(TBasisCenters) :: this

    if (allocated(this%x)) deallocate(this%x)
    if (allocated(this%matrixToBasis)) deallocate(this%matrixToBasis)

  end subroutine destroy_TBasis

end module ln_structure
