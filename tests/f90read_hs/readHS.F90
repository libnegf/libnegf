! Module read the real-space HS as written from DFTB+
! The real space Hamiltonian is unfolded following the logic:
! loop  iAt1 = 1, nAtoms
!    loop iAt2 = nNeigbours(iAt1) but keeping ONLY the blocks of H
!    such that the folded image of folded(iAt2) >= iAt1
!
! So only the lower-diagonal part of H is written to file   
!
module readHS
  use constants, only : dp    
  implicit none

  ! All public data structures to be filled when reading from file
  real(dp), allocatable :: H0(:), S(:)
  integer, allocatable :: nNeighbours(:)
  integer, allocatable :: iNeighbour(:,:)
  integer, allocatable :: iAtomStart(:)
  integer, allocatable :: iPair(:,:)
  integer, allocatable :: img2CentCell(:)
  real(dp), allocatable :: cellVec(:,:)
  integer, allocatable :: iIndexVec(:)
  integer :: nAtoms
  
  type intArray
    integer, allocatable :: iCell(:,:) 
  end type intArray
  
  type(intArray), allocatable :: cells(:)

  type TOrbitals
    integer, allocatable :: nOrbAtom(:)
    integer :: mOrb 
  end type TOrbitals

  type(TOrbitals) :: orb

  contains

  subroutine read_dftb_hs()    

    integer :: fu1, fu2
    integer :: int1, int2, int3, int4, int5, iVec(3)
    integer :: iAt1, iAt2, iNeigh, maxneig, matSize=0
    integer :: nOrbs1, nOrbs2, nNeigh, last
    character(50) :: strForm 

    open(newunit=fu1, file='hamreal1.dat')
    open(newunit=fu2, file='overreal.dat')

    read(fu1, *) !read # header line  
    read(fu2, *) !read # header line
    read(fu1, *) nAtoms
    read(fu2, *) int1
    !write(*,*) 'Number of atoms: ',nAtoms
    if (int1 /= nAtoms) then
      error stop "Number of atoms do not match"
    end if
    allocate(nNeighbours(nAtoms))
    allocate(iAtomStart(0:nAtoms))
    allocate(cells(nAtoms))
    allocate(img2CentCell(nAtoms))
    if (allocated(orb%nOrbAtom)) then
       if (size(orb%nOrbAtom) .ne. nAtoms) then
          error stop "Incompatible nAtoms between shifs and H,S"
       end if
    else
       allocate(orb%nOrbAtom(nAtoms))
    end if 
    read(fu1, *) !read # header line  
    read(fu2, *) !read # header line
    ! Read the list of atoms with number of neighbours and orbitals
    iAtomStart(0) = 1
    do iAt1 = 1, nAtoms
      read(fu1, *) int1, nNeigh, nOrbs1
      read(fu2, *) int2, int3, int4
      if (int3 .ne. nNeigh) then
         write(strForm,"(I0,2x,I0)") int3, nNeigh
         error stop "Different neighbours, "//trim(strForm)
      end if       
      if (int4 .ne. nOrbs1) then
         write(strForm,"(I0,2x,I0)") int4, nOrbs1
         error stop "Different orbitals, "//trim(strForm)
      end if       
      matSize = matSize + nOrbs1 * nOrbs1 * nNeigh 
      nNeighbours(iAt1) = nNeigh - 1
      orb%nOrbAtom(iAt1) = nOrbs1
      allocate(cells(iAt1)%iCell(3,0:nNeighbours(iAt1)))
      iAtomStart(iAt1) = iAtomStart(iAt1-1) + orb%nOrbAtom(iAt1)
      ! This is a moked array defined as the identity
      ! iNeighbours already maps into the folded atoms
      img2CentCell(iAt1) = iAt1
    end do
    orb%mOrb = maxval(orb%nOrbAtom) 
    maxneig = maxval(nNeighbours)
    !write(*,*) 'Max number of orbitals: ',orb%mOrb
    !write(*,*) 'Max number of neighbours: ',maxneig
    !write(*,*) 'Size of allocated sparse: ',matSize
    allocate(H0(matSize))
    allocate(S(matSize))
    allocate(iNeighbour(0:maxneig,nAtoms))
    allocate(iPair(0:maxneig,nAtoms))
    ! In order to fill iPair we need the actual list of neighbours
    ! We need to read the matrix once
    iPair = 0
    last = 0  
    do iAt1 = 1, nAtoms
      do iNeigh = 0, nNeighbours(iAt1)
        read(fu1, *) !read # description line
        read(fu2, *) !read # description line
        read(fu1, *) int1, int2, iAt2, iVec(:)
        read(fu2, *) int3, int4, int5, cells(iAt1)%iCell(:,iNeigh)
        !print*, iAt1, iNeigh
        !print*, int1, int2, iAt2, int5
        if (int1 .ne. iAt1) then
           write(strForm,"(I0,2x,I0)") int1, iAt1   
           error stop "Wrong atom index, "//trim(strForm)
        end if   
        if (int2 .ne. iNeigh) then
           write(strForm,"(I0,2x,I0)") int2, iNeigh   
           error stop "Wrong neighbour index, "//trim(strForm)
        end if   
        if (int5 .ne. iAt2) then
           write(strForm,"(I0,2x,I0)") int5, iAt2   
           error stop "Different neighbour atoms, "//trim(strForm)
        end if   
        if (any(iVec .ne. cells(iAt1)%iCell(:,iNeigh))) then
           write(strForm,"(I0,2x,I0)") iAt1, iNeigh  
           error stop "Different cells for atom, neigh "//trim(strForm)   
        end if

        !print*,'neighbour ',iAt2
        iNeighbour(iNeigh,iAt1) = iAt2
        nOrbs1 = orb%nOrbAtom(iAt1)
        nOrbs2 = orb%nOrbAtom(iAt2)
        ! ind = iPair(iNeigh,iAt1) + 1  => ind = iPair(0,1) + 1 == 1 =>  
        ! sparse(ind:ind+nOrb1*nOrb2-1) 
        iPair(iNeigh, iAt1) = last
        read(fu1, *)  !read # Matrix   
        read(fu2, *)  !read # Matrix 
        !write(strForm, "(A,I0,A)") "(", nOrbs2, "ES24.15)"
        read(fu1, *) H0(last+1:last+nOrbs1*nOrbs2)
        read(fu2, *) S(last+1:last+nOrbs1*nOrbs2)
        last = last + nOrbs1 * nOrbs2  
      end do
    end do

    close(fu1)
    close(fu2)
    
    call create_cellvec(cells, cellVec, iIndexVec)

  end subroutine read_dftb_hs

  ! converts cells into the arrays cellVec/iIndexVec used in the folding routines
  subroutine create_cellvec(cells, cellVec, iIndexVec)
    type(intArray), intent(in) :: cells(:)
    real(dp), intent(out), allocatable :: cellVec(:,:)
    integer, intent(out), allocatable :: iIndexVec(:)

    integer :: iAt1, iAt2, iNeigh, ii, jj, kk, iVec, mt
    integer :: imax(3), imin(3), tmpVec(3)
    integer :: nAtoms, nNeigh, nCells
    integer, allocatable :: iCellVec(:,:)
    logical :: found
    character(7) :: frm

    nAtoms = size(cells)
    ! Find minimal/maximal indices in all directions
    imax = 0
    imin = 0
    do iAt1 = 1, nAtoms
      do ii = 1, 3
        !print*,iAt1, ii
        !write(frm,"(A,I0,A)") "(",size(cells(iAt1)%iCell,2),"I3)"
        !write(*,frm) cells(iAt1)%iCell(ii,:)
        mt = maxval(cells(iAt1)%iCell(ii,:))
        !print*,'maxval:',ii,mt
        if (mt > imax(ii)) imax(ii) = mt   
        mt = minval(cells(iAt1)%iCell(ii,:))
        !print*,'minval:',ii,mt
        if (mt < imin(ii)) imin(ii) = mt
      end do  
    end do

    !print*,'imin:',imin
    !print*,'imax:',imax

    ! compute the total number of cells 
    nCells = 1
    do ii = 1, 3
      nCells = nCells * (imax(ii)-imin(ii)+1)
    end do 

    !print*,'nCells',nCells

    allocate(iCellVec(3,nCells))
    ! populate cellVec with all cell combinations
    iVec = 0
    do ii = imin(1), imax(1)
      do jj = imin(2), imax(2)
        do kk = imin(3), imax(3)
           iVec = iVec + 1
           icellVec(:, iVec) = [ii, jj, kk]  
        end do
      end do
    end do  
    
    ! Sort the cell Vectors by distance to 0
    ! i.e. (0,0,0), (-1,0,0), (1,0,0), ...
    do ii = 1, nCells
      do jj = ii+1, nCells
         if (mynorm(iCellVec(:,ii)) .gt. mynorm(iCellVec(:,jj))) then
            call swap(iCellVec(:,ii), iCellVec(:,jj))
         end if
      end do
    end do   
    !print*,'Sorted cells:'
    !do ii = 1, nCells
    !  print*, iCellVec(:,ii)
    !end do   

    ! the index vector is the size of the supercell >> nAtoms 
    allocate( iIndexVec(nAtoms*size(icellVec,2)) )
    if (allocated(img2CentCell)) deallocate(img2CentCell)
    allocate( img2CentCell(nAtoms*size(icellVec,2)) )

    ! we need to scan cells in order to find the right cell index
    ! we also need to correct iNeighbour and img2CentCell    
    loop:do iAt1 = 1, nAtoms
      nNeigh = size(cells(iAt1)%iCell,2)-1
      do iNeigh = 0, nNeigh
         !scan cellVec to find the index
         tmpVec = cells(iAt1)%iCell(:,iNeigh)
         found = .false.
         do ii = 1, nCells
            if (all( icellVec(:,ii).eq.tmpVec(:) )) then
               ! unfolding atoms into periodic replicas
               ! Originally iNeighbour contains the folded atom which is 
               ! here placed to the correct cell and then iNeighbour is redefined.
               iAt2 = iNeighbour(iNeigh,iAt1) + (ii-1)*Natoms
               iIndexVec(iAt2) = ii
               img2CentCell(iAt2) = iNeighbour(iNeigh, iAt1)
               iNeighbour(iNeigh, iAt1) = iAt2
               found = .true.   
            end if
         end do
         if (.not.found) then
            print*, tmpVec   
            error stop "Cell index not found"
         end if   
      end do
    end do loop


    allocate(cellVec(3,size(iCellVec,2)))
    cellVec = 1.0_dp * iCellVec
    deallocate(iCellVec)

    contains
    integer function mynorm(vec)
      integer :: vec(3)
      mynorm = dot_product(vec,vec)
    end function mynorm 

    subroutine swap(vec1, vec2)
      integer :: vec1(3), vec2(3)
      tmpVec = vec1
      vec1 = vec2
      vec2 = tmpVec
    end subroutine swap

  end subroutine create_cellvec


  !> Writes a sparse matrix to a file.
  subroutine writeSparse(fname, sparse, iNeighbour, nNeighbourSK, iAtomStart, iPair, img2CentCell,&
      & iCellVec, cellVec)

    !> Name of the file to write the matrix to.
    character(len=*), intent(in) :: fname

    !> Sparse matrix.
    real(dp), intent(in) :: sparse(:)

    !> Neighbour list index.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours.
    integer, intent(in) :: nNeighbourSK(:)

    !> Offset array in the square matrix.
    integer, intent(in) :: iAtomStart(:)

    !> Pair indexing array.
    integer, intent(in) :: iPair(0:,:)

    !> Mapping of the atoms to the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Index of the cell translation vectors for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Cell translation vectors.
    real(dp), intent(in) :: cellVec(:,:)

    integer :: fd, nAtom
    integer :: iAt1, iAt2, iAt2f, iNeigh, iOrig, nOrb1, nOrb2
    character(20) :: strForm

    !if (.not. tIoProc) then
    !  return
    !end if

    nAtom = size(nNeighbourSK)

    open(newunit=fd, file=fname, form="formatted", status="replace")
    write(fd, "(A1,A10)") "#", "NATOM"
    write(fd, "(1X,I10)") nAtom
    write(fd, "(A1,A10,A10,A10)") "#", "IATOM", "NNEIGH", "NORB"
    do iAt1 = 1, nAtom
      write(fd, "(1X,I10,I10,I10)") iAt1, nNeighbourSK(iAt1) + 1, iAtomStart(iAt1+1)&
          & - iAtomStart(iAt1)
    end do

    do iAt1 = 1, nAtom
      nOrb1 = iAtomStart(iAt1+1) - iAtomStart(iAt1)
      do iNeigh = 0, nNeighbourSK(iAt1)
        iOrig = iPair(iNeigh,iAt1) + 1
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        nOrb2 = iAtomStart(iAt2f+1) - iAtomStart(iAt2f)
        write(strForm, "(A,I0,A)") "(", nOrb2, "ES24.15)"
        write(fd, "(A1,A10,A10,A10,3A10)") "#", "IATOM1", "INEIGH", "IATOM2F", "ICELL(1)",&
            & "ICELL(2)", "ICELL(3)"
        write(fd, "(1X,I10,I10,I10,3I10)") iAt1, iNeigh, iAt2f, int(cellVec(:,iCellVec(iAt2)))
        write(fd, "(A1,A)") "#", " MATRIX"
        write(fd, strForm) sparse(iOrig:iOrig+nOrb1*nOrb2-1)
      end do
    end do
    close(fd)

  end subroutine writeSparse
   

  !> Writes Diagonal of a sparse matrix to a file.
  subroutine writeSparseDiagonal(fname, sparse, iNeighbour, nNeighbourSK, iAtomStart, iPair, img2CentCell,&
      & iCellVec, cellVec)

    !> Name of the file to write the matrix to.
    character(len=*), intent(in) :: fname

    !> Sparse matrix.
    real(dp), intent(in) :: sparse(:)

    !> Neighbour list index.
    integer, intent(in) :: iNeighbour(0:,:)

    !> Number of neighbours.
    integer, intent(in) :: nNeighbourSK(:)

    !> Offset array in the square matrix.
    integer, intent(in) :: iAtomStart(:)

    !> Pair indexing array.
    integer, intent(in) :: iPair(0:,:)

    !> Mapping of the atoms to the central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Index of the cell translation vectors for each atom.
    integer, intent(in) :: iCellVec(:)

    !> Cell translation vectors.
    real(dp), intent(in) :: cellVec(:,:)

    integer :: fd, nAtom
    integer :: iAt1, iAt2, iAt2f, ii, iOrig, nOrb1, nOrb2
    character(20) :: strForm
    real(dp) :: tmp(9,9)

    !if (.not. tIoProc) then
    !  return
    !end if

    nAtom = size(nNeighbourSK)

    open(newunit=fd, file=fname, form="formatted", status="replace")

    do iAt1 = 1, nAtom
      nOrb1 = iAtomStart(iAt1+1) - iAtomStart(iAt1)
      iOrig = iPair(0,iAt1)
      !iAt2 = iNeighbour(0, iAt1)
      !iAt2f = img2CentCell(iAt2)
      !nOrb2 = iAtomStart(iAt2f+1) - iAtomStart(iAt2f)
      !if (nOrb1 /= nOrb2) stop
      tmp(1:nOrb1,1:nOrb1) = reshape( &
            & sparse(iOrig+1:iOrig+nOrb1*nOrb1),(/nOrb1,nOrb1/) )
      do ii = 1, nOrb1
        write(fd, strForm) tmp(ii,ii)
      end do  
    end do
    close(fd)

  end subroutine writeSparseDiagonal


  ! Application of a shift H_ij = H_ij + (S_ij V_j + V_i S_ij)/2
  subroutine apply_shifts(ham, over, shift)
    real(dp), intent(inout) :: ham(:)
    real(dp), intent(in) :: over(:)
    real(dp), intent(in) :: shift(:)

    integer :: nAtom, iAt1, iAt2, iAt2f, iBlk1, iBlk2 
    integer :: nOrb1, nOrb2, iOrig1, iNeigh, ii
    real(dp) :: tmpH(orb%mOrb,orb%mOrb), tmpS(orb%mOrb,orb%mOrb)
    real(dp) :: shiftBlk1(orb%mOrb,orb%mOrb)
    real(dp) :: shiftBlk2(orb%mOrb,orb%mOrb)
  
    nAtom = size(orb%nOrbAtom)
    shiftBlk1 = 0.0_dp
    shiftBlk2 = 0.0_dp

    do iAt1 = 1, nAtom
      nOrb1 = orb%nOrbAtom(iAt1)
      do iNeigh = 0, nNeighbours(iAt1)
        iAt2 = iNeighbour(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        nOrb2 = orb%nOrbAtom(iAt2f)
        iOrig1 = iPair(iNeigh, iAt1)
        iBlk1 = iAtomStart(iAt1)-1
        iBlk2 = iAtomStart(iAt2f)-1
        do ii = 1, nOrb1
          shiftBlk1(ii,ii) = shift(iBlk1+ii)
        end do
        do ii = 1, nOrb2
          shiftBlk2(ii,ii) = shift(iBlk2+ii)
        end do

        tmpS(1:nOrb2,1:nOrb1) = reshape( &
            & over(iOrig1+1:iOrig1+nOrb2*nOrb1),(/nOrb2,nOrb1/) )

        !  H_ij = 1/2 S_ik (delta_kj V_j) + 1/2 (V_i delta_ik) S_kj
        !  1/2 (S21 V1 + V2 S21)
        tmpH(1:nOrb2,1:nOrb1) = 0.5_dp * ( &
            & matmul(tmpS(1:nOrb2,1:nOrb1), shiftBlk1(1:nOrb1,1:nOrb1)) + &
            & matmul(shiftBlk2(1:nOrb2,1:nOrb2), tmpS(1:nOrb2,1:nOrb1)) )

        ham(iOrig1+1:iOrig1+nOrb2*nOrb1) = &
            & ham(iOrig1+1:iOrig1+nOrb2*nOrb1) + &
            & reshape(tmpH(1:nOrb2, 1:nOrb1), (/nOrb2*nOrb1/))

      end do
    end do

  end subroutine apply_shifts

  !> Read the Hamiltonian potential shifts from file
  subroutine readShifts(fShifts, nAtom, nSpin, shiftPerL)

    !> filename where shifts are stored
    character(*), intent(in) :: fShifts

    !> number of atoms and spin blocks
    integer, intent(in) :: nAtom, nSpin

    !> potential shifts (shell,atom,spin) charge/mag is used
    real(dp), allocatable, intent(inout) :: shiftPerL(:,:,:)

    integer :: fdH, nAtomSt, nSpinSt, mOrbSt, mShellSt, ii, jj


    open(newunit=fdH, file=fShifts, form="formatted")
    read(fdH, *) nAtomSt, mShellSt, mOrbSt, nSpinSt

    if (nAtom .ne. nAtomSt) then
       error stop "Incompatible nSpin between shifts and H,S"
    end if   
    if (nSpin .ne. nSpinSt) then
       error stop "Incompatible nSpin between shifts and H,S"
    end if   
    if (allocated(orb%nOrbAtom)) then
       if (size(orb%nOrbAtom) .ne. nAtomSt) then
          error stop "Incompatible nAtoms between shifs and H,S"
       end if
    else
       allocate(orb%nOrbAtom(nAtomSt))
    end if 

    allocate(shiftPerL(mShellSt,nAtomSt,nSpinSt))
    shiftPerL(:,:,:) = 0.0_dp

    do ii = 1, nAtomSt
      read(fdH, *) orb%nOrbAtom(ii), (shiftPerL(:,ii,jj), jj = 1, nSpinSt)
    end do

    close(fdH)

  end subroutine readShifts

end module readHS
