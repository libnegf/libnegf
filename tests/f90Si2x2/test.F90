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

  integer, parameter :: norbs = 9 
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

  if (globalComm%lead) write(*,*) 'Initializing libNEGF'
  call init_negf(pnegf)

  if (globalComm%lead)  write(*,*) 'Setup MPI communicators (cartesian grid)'
  call negf_cart_init(globalComm, 1, cartComm, enComm, kComm)
  call negf_mpi_init(pnegf, cartComm, enComm, kComm)

  if (globalComm%lead)  write(*,*) 'Import Hamiltonian'
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
  surfend   = [320*norbs, 384*norbs]
  surfstart = [320*norbs+1, 384*norbs+1]
  contend   = [384*norbs, 448*norbs]
  plend = [32*norbs, 64*norbs, 96*norbs, 128*norbs, 160*norbs, 192*norbs, 224*norbs, &
        & 256*norbs, 288*norbs, 320*norbs]
  nPLs = 10 
  cblk = [10, 1]
  
  ! Note: number of contacts should be set first
  if (globalComm%lead)  write(*,*) 'Set contact and device structures'
  call init_contacts(pnegf, 2)
  call init_structure(pnegf, 2, surfstart, surfend, contend, nPLs, plend, cblk)

  ! Setting other parameters
  Ef = -3.7_dp  ! eV
  bias = 0.05_dp           ! V
  mu = [(Ef-bias)/eV2Hartree, (Ef+bias)/eV2Hartree]
  kt = [1.0d-5, 1.0d-5]
  
  ! Here we set the parameters, only the ones different from default
  call get_params(pnegf, params)
  params%verbose = 100 
  params%Emin = -4.560_dp/eV2Hartree
  params%Emax = -4.245_dp/eV2Hartree
  params%Estep = 0.005_dp/eV2Hartree
  params%mu(1:2) = mu
  params%kbT_t(1:2) = kt
    
  params%spin = 1
  params%ikpoint = 1
  params%kwght = 1.0_dp
  call set_params(pnegf, params)
  
  if (globalComm%lead)  write(*,*) 'create csr Hamiltonian'
  call init(csrHam, iAtomStart, iNeighbour, nNeighbours, img2CentCell, orb)
  call init(csrOvr, csrHam)

  if (globalComm%lead)  write(*,*) 'fold H0 to csr'
  call foldToCSR(csrHam, H0, kPoints(:,1), iAtomStart, iPair, iNeighbour, nNeighbours,&
      & img2CentCell, iIndexVec, cellVec, orb)
  if (globalComm%lead)  write(*,*) 'fold S to csr'
  call foldToCSR(csrOvr, S, kPoints(:,1), iAtomStart, iPair, iNeighbour, nNeighbours,&
      & img2CentCell, iIndexVec, cellVec, orb)

  if (globalComm%lead)  write(*,*) 'create HS container of size 1'
  call create_HS(negf, 1)

  pcsrHam => csrHam
  pcsrOvr => csrOvr

  if (globalComm%lead)  write(*,*) 'pass HS to Negf'
  call pass_HS(pnegf, pcsrHam, pcsrOvr)

  if (globalComm%lead)  write(*,*) 'Compute current'
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

  
  ! I/O transmission and current
  ierr = 0
  if (globalComm%lead) then
     open(111,file="transmission.dat")   
     nSteps = int((params%Emax - params%Emin)/params%Estep) + 1  
     do ii = 1, nSteps
       write(111,*) params%Emin + (ii-1)*params%Estep, tunnMat(ii,1)
     end do
     close(111)
     allocate(tunn_ref(nSteps))
     open(112,file="transmission_ref.dat")   
     do ii = 1, nSteps
       read(112,*) bias, tunn_ref(ii)
     end do
     close(112)
     write(*,*) 'Current ',current
     !write(*,*) 'Reference Current ',
  
     if (any(abs(tunnMat(1:nSteps,1) - tunn_ref(1:nSteps)) > 1e-5)) then
        write(*,*) maxval(abs(tunnMat(1:nSteps,1) - tunn_ref(1:nSteps)))
        write(*,*) "Tunneling reference not met"
        ierr = 1
     end if  
     
     !if (abs(current - 1.549625501099260E-005)> 1e-5) then
     !   write(*,*) "Current reference not met"
      !  ierr = 0
     !end if  
  end if

  call mpifx_barrier(enComm, impierr)

  if (globalComm%lead)   write(*,*) 'Destroy negf'
  call destroy_negf(pnegf)
  
  !call writePeakInfo(6)
  !call writeMemInfo(6)
  
  call mpifx_finalize();
 
  if (ierr /= 0) then 
     error stop "Errors found"   
  end if   
  if (globalComm%lead) write(*,*) 'Done'
  
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
