program test_readHS
  use constants
  use ln_constants, only : eovh, HAR
  use readHS
  use matconv
  use libnegf
  use lib_param
  use integrations
  use libmpifx_module
  implicit none

  integer, parameter :: norbs = 4
  real(dp), parameter :: eV2Hartree = HAR
  Type(Tnegf), target :: negf
  Type(Tnegf), pointer :: pnegf
  Type(lnParams) :: params
  integer, allocatable :: surfstart(:), surfend(:), contend(:), plend(:), cblk(:)
  integer, allocatable :: iCellVec(:,:)
  real(dp), allocatable :: mu(:), kt(:), tunn_ref(:)
  real(dp) :: Ef, current, bias
  real(dp), allocatable :: kPoints(:,:), tunnMat(:,:)
  integer :: ierr, impierr, nPLs, ii, nSteps
  type(mpifx_Comm) :: globalComm, cartComm, kComm, enComm
  type(z_CSR), target :: csrHam, csrOvr
  type(z_CSR), pointer :: pcsrHam, pcsrOvr

  call mpifx_init(impierr);
  call globalComm%init();

  pnegf => negf

  write(*,*) 'Initializing libNEGF'
  call init_negf(pnegf)

  write(*,*) 'Setup MPI communicators (cartesian grid)'
  call negf_cart_init(globalComm, 1, cartComm, enComm, kComm)
  call negf_mpi_init(pnegf, cartComm, enComm, kComm)

  write(*,*) 'Import Hamiltonian'
  call read_dftb_hs()

  !call writeSparse("hamreal1_w.dat", H0, iNeighbour, nNeighbours, iAtomStart, iPair, img2CentCell,&
  !    & iIndexVec, cellVec)

  allocate(kPoints(3,1))
  kPoints = 0.0_dp

  ! CONTACT DEFINITION
  !  ---------------------------------- - - -
  !  | S        ||  PL1    |   PL2    |
  !  ---------------------------------- - - -
  !  surfstart  surfend               contend

  ! Si 2x2
  !surfend   = [320*norbs, 384*norbs]
  !surfstart = [320*norbs+1, 384*norbs+1]
  !contend   = [384*norbs, 448*norbs]
  !plend = [32*norbs, 64*norbs, 96*norbs, 128*norbs, 160*norbs, 192*norbs, 224*norbs, &
  !      & 256*norbs, 288*norbs, 320*norbs]

  surfend   = [80*norbs,   120*norbs]
  surfstart = [80*norbs+1, 120*norbs+1]
  contend   = [120*norbs,  160*norbs]
  plend = [40*norbs, 80*norbs]
  nPLs = 2
  cblk = [2, 1]
  
  ! Note: number of contacts should be set first
  write(*,*) 'Set contact and device structures'
  call init_contacts(pnegf, 2)
  call init_structure(pnegf, 2, surfstart, surfend, contend, nPLs, plend, cblk)

  ! Setting other parameters
  Ef = -5.7_dp  ! eV
  bias = 0.05_dp           ! V
  mu = [(Ef-bias)/eV2Hartree, (Ef+bias)/eV2Hartree]
  kt = [1.0d-5, 1.0d-5]
  
  ! Here we set the parameters, only the ones different from default
  call get_params(pnegf, params)
  params%verbose = 30
  params%Emin = -6.0_dp/eV2Hartree
  params%Emax = -5.5_dp/eV2Hartree
  params%Estep = 0.02_dp/eV2Hartree
  params%mu(1:2) = mu
  params%kbT_t(1:2) = kt
    
  params%spin = 1
  params%ikpoint = 1
  params%kwght = 1.0_dp
  call set_params(pnegf, params)
  
  write(*,*) 'create csr Hamiltonian'
  call init(csrHam, iAtomStart, iNeighbour, nNeighbours, img2CentCell, orb)
  call init(csrOvr, csrHam)

  write(*,*) 'fold H0 to csr'
  call foldToCSR(csrHam, H0, kPoints(:,1), iAtomStart, iPair, iNeighbour, nNeighbours,&
      & img2CentCell, iIndexVec, cellVec, orb)
  write(*,*) 'fold S to csr'
  call foldToCSR(csrOvr, S, kPoints(:,1), iAtomStart, iPair, iNeighbour, nNeighbours,&
      & img2CentCell, iIndexVec, cellVec, orb)

  write(*,*) 'create HS container of size 1'
  call create_HS(negf, 1)

  pcsrHam => csrHam
  pcsrOvr => csrOvr

  write(*,*) 'pass HS to Negf'
  call pass_HS(pnegf, pcsrHam, pcsrOvr)

  write(*,*) 'Compute current'
  call compute_current(pnegf)
 
  ! GATHER MPI partial results on node 0
  if (allocated(pnegf%tunn_mat)) then
    tunnMat = pnegf%tunn_mat
  else
    error stop "Error Transmission not created"
  end if
  
  current = pnegf%currents(1)*eovh

  call mpifx_reduceip(enComm, tunnMat, MPI_SUM)
  call mpifx_reduceip(enComm, current, MPI_SUM)

  ! Checking against references
  tunn_ref = [5.999357, 5.998753, 5.997118, 5.990532, 5.888102, 2.030200, 2.002432, 2.000657, &
           &  2.000262, 2.000128, 2.000070, 2.000041, 2.000026, 2.000017, 2.000011, 2.000007, &
           &  2.000005, 2.000003, 2.000001, 2.000000, 1.999999, 1.999999, 1.999998, 1.999997, &
           &  1.999997, 1.999997]
  
  ! I/O transmission and current
  ierr = 0
  if (enComm%lead) then
     nSteps = int((params%Emax - params%Emin)/params%Estep) + 1  
     do ii = 1, nSteps
       write(*,*) params%Emin + (ii-1)*params%Estep, tunnMat(ii,1), TunnMat(ii,1)-tunn_ref(ii)
     end do
     write(*,*) 'Current ',current
     write(*,*) 'Reference Current ',1.5496255e-5
  
     if (any(abs(tunnMat(:,1) - tunn_ref(:)) > 1e-5)) then
        write(*,*) "Tunneling reference not met"
        ierr = 1
     end if  
     
     if (abs(current - 1.549625501099260E-005)> 1e-5) then
        write(*,*) "Current reference not met"
        ierr = 2
     end if  
  end if

  call mpifx_barrier(enComm, impierr)

  write(*,*) 'Destroy negf'
  call destroy_negf(pnegf)
  
  !call writePeakInfo(6)
  !call writeMemInfo(6)
  
  call mpifx_finalize();
 
  if (ierr /= 0) then 
     error stop "Errors found"   
  end if   
  write(*,*) 'Done'
  
  contains
  !> utility to sum up partial results over SK communicator
  subroutine add_ks_results(kscomm, mat, matSKRes)

    type(mpifx_comm), intent(in) :: kscomm

    !> sum total
    real(dp), allocatable, intent(inout) :: mat(:,:)

    !> k-resolved sum
    real(dp), allocatable, intent(inout)  :: matSKRes(:,:,:)

    if (allocated(mat)) then
      call mpifx_reduceip(kscomm, mat, MPI_SUM)
    endif

    if (allocated(matSKRes)) then
      call mpifx_reduceip(kscomm, matSKRes, MPI_SUM)
    endif

  end subroutine add_ks_results

end program test_readHS
