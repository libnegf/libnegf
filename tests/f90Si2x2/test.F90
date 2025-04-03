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
  integer, allocatable :: iCellVec(:,:), iKProcs(:)
  real(dp), allocatable :: mu(:), kt(:), tunn_ref(:)
  real(dp) :: Ef, current, bias
  real(dp), allocatable :: kPoints(:,:), tunnMat(:,:), tunnMatSK(:,:)
  integer :: ierr, impierr, nPLs, ii, jj, nSteps, nK, nn, nGroups, iKS
  type(mpifx_Comm) :: globalComm, cartComm, kComm, enComm
  type(z_CSR), target :: csrHam, csrOvr
  type(z_CSR), pointer :: pcsrHam, pcsrOvr

  call mpifx_init(impierr);
  call globalComm%init();
  nGroups = 1 

  pnegf => negf

  if (globalComm%lead) write(*,*) 'Initializing libNEGF'
  call init_negf(pnegf)

  if (globalComm%lead)  write(*,*) 'Setup MPI communicators (cartesian grid)'
  call negf_cart_init(globalComm, nGroups, cartComm, enComm, kComm)
  call negf_mpi_init(pnegf, cartComm, enComm, kComm)

  if (globalComm%lead)  write(*,*) 'Import Hamiltonian'
  call read_dftb_hs()

  !call writeSparse("hamreal1_w.dat", H0, iNeighbour, nNeighbours, iAtomStart, iPair, &
  !    & img2CentCell, iIndexVec, cellVec)

  ! Generation of 2-d k-mesh assuming simple orthogonal lattice vectors
  ! k-mesh is simply defined on a square-lattice [1,0]x[0,1]
  ! Monkhorst-Pack nxn is simply a grid with n points at 1/n+1
  ! Reduction by symmetry is possible reducing the BZ by 2: [1/2,0]x[0,1]
  ! k = 2*pi*[1/a, 1/b]  R = [n*a,m*b]
  ! exp(ik*R) = exp[i*2*pi*(n+m)]
  ! NOTE: Transport is ASSUMED along z
  ! +-------+
  ! | o | o |
  ! +---+---+      
  ! | o | o |
  ! +-------+
  ! 0       1 
  nn = 2 
  nk = nn*nn/2 ! Reduced by symmetry
  
  allocate(kPoints(3,nk))
  kPoints = 0.0_dp
  nk = 0
  do ii = 1, nn, 2
    do jj = 1, 2*nn, 2
       nk = nk + 1
       kPoints(:,nk) = 1.0_dp/(2*nn)*[ii, jj, 0]
    end do   
  end do
  
  ! Assign k-points to groups
  if (mod(nk,nGroups) /= 0) then
    error stop "nk not commensurate with nGroups"
  end if
  allocate(iKProcs(nk/nGroups))

  jj = 0
  do ii = 1, nk
    if (mod(ii, nGroups) == kComm%rank) then
        jj = jj + 1  
        iKProcs(jj)  = ii
    end if
  end do  

  if (globalComm%lead) then
     write(*,*) 'k-points'
     do ii = 1, nk
       write(*,'(I0,A,3f14.6,A)') ii,'[',kPoints(:,ii),']'
     end do
  end if


  ! CONTACT DEFINITION
  !  ---------------------------------- - - -
  !  | S        ||  PL1    |   PL2    |
  !  ---------------------------------- - - -
  !  surfstart  surfend               contend

  ! ----------------------------------------------------------------------------------
  ! Si 2x2 structure definition
  ! ----------------------------------------------------------------------------------
  surfend   = [320*norbs, 384*norbs]
  surfstart = [320*norbs+1, 384*norbs+1]
  contend   = [384*norbs, 448*norbs]
  plend = [32*norbs, 64*norbs, 96*norbs, 128*norbs, 160*norbs, 192*norbs, 224*norbs, &
        & 256*norbs, 288*norbs, 320*norbs]
  nPLs = 10 
  cblk = [10, 1]
  
  ! Note: number of contacts should be set first
  ! ----------------------------------------------------------------------------------
  if (globalComm%lead)  write(*,*) 'Set contact and device structures'
  call init_contacts(pnegf, 2)
  call init_structure(pnegf, 2, surfstart, surfend, contend, nPLs, plend, cblk)

  ! Setting transmission parameters
  ! ----------------------------------------------------------------------------------
  Ef = -4.4_dp  ! eV
  bias = 0.05_dp           ! V
  mu = [(Ef-bias)/eV2Hartree, (Ef+bias)/eV2Hartree]
  kt = [1.0d-5, 1.0d-5]
  
  ! Here we set the parameters, only the ones different from default
  ! ----------------------------------------------------------------------------------
  call get_params(pnegf, params)
  params%verbose = 100 
  params%Emin = -4.560_dp/eV2Hartree
  params%Emax = -4.485_dp/eV2Hartree
  params%Estep = 0.005_dp/eV2Hartree
  ! nStep = (-4.485 + 4.560)/0.005 + 1 = 16
  params%delta = 1.e-5_dp  !Already in Hartree
  params%mu(1:2) = mu
  params%kbT_t(1:2) = kt
  ! setting spin and k-index & kwght
  params%spin = 1
  params%ikpoint = 1
  params%kwght = 1.0_dp/nk
  call set_params(pnegf, params)

  ! ----------------------------------------------------------------------------------
  ! Prepare Hamiltonian and Overlap CSR sparsity map
  ! ----------------------------------------------------------------------------------
  if (globalComm%lead)  write(*,*) 'create csr Hamiltonian'
  call init(csrHam, iAtomStart, iNeighbour, nNeighbours, img2CentCell, orb)
  call init(csrOvr, csrHam)
  pcsrHam => csrHam
  pcsrOvr => csrOvr
  current = 0.0_dp

  ! Loop over k-points
  ! ----------------------------------------------------------------------------------
  kloop:do ii = 1, size(iKProcs) 

    iKS = iKProcs(ii)

    if (globalComm%lead)  write(*,*) 'fold H0 to csr, k=',iKS
    call foldToCSR(csrHam, H0, kPoints(:,iKS), iAtomStart, iPair, iNeighbour, nNeighbours,&
      & img2CentCell, iIndexVec, cellVec, orb)
    if (globalComm%lead)  write(*,*) 'fold S to csr, k=',iKS
    call foldToCSR(csrOvr, S, kPoints(:,iKS), iAtomStart, iPair, iNeighbour, nNeighbours,&
      & img2CentCell, iIndexVec, cellVec, orb)

    if (globalComm%lead)  write(*,*) 'create HS container of size 1'
    call create_HS(negf, 1)

    if (globalComm%lead)  write(*,*) 'pass HS to Negf'
    call pass_HS(pnegf, pcsrHam, pcsrOvr)

    if (globalComm%lead)  write(*,*) 'Compute current'
    call compute_current(pnegf)
 
    if (.not.allocated(pnegf%tunn_mat)) then
      error stop "Error Transmission not created"
    end if
    if (.not.allocated(tunnMat)) then
      allocate(tunnMat(size(pnegf%tunn_mat,1),1))
      tunnMat = 0.0_dp
    end if 
    tunnMat = tunnMat + pnegf%tunn_mat
    
    if (.not.allocated(tunnMatSK)) then
       allocate(tunnMatSK(size(tunnMat,1),nk))
    end if
    tunnMatSK(:,iKS) = pnegf%tunn_mat(:,1)

    current = current + pnegf%currents(1)*eovh

  end do kloop
    
  ! GATHER MPI partial results on enComm%rank == 0
  call mpifx_reduceip(enComm, tunnMat, MPI_SUM)
  call mpifx_reduceip(enComm, tunnMatSK, MPI_SUM)
  call mpifx_reduceip(enComm, current, MPI_SUM)
    
  ! GATHER MPI partial results on kComm%rank == 0
  if (enComm%lead) then
    call mpifx_reduceip(kComm, tunnMat, MPI_SUM)
    call mpifx_reduceip(kComm, tunnMatSK, MPI_SUM)
    call mpifx_reduceip(kComm, current, MPI_SUM)
  end if

  ! ----------------------------------------------------------------------------------
  ! I/O transmission and current
  ! ----------------------------------------------------------------------------------
  ierr = 0
  if (globalComm%lead) then
     open(111,file="transmission.dat")   
     nSteps = size(tunnMat,1)
     !nint((params%Emax - params%Emin)/params%Estep) + 1 
     do ii = 1, nSteps
       write(111,*) (params%Emin + (ii-1)*params%Estep)*eV2Hartree, tunnMat(ii,1)
     end do
     close(111)
     open(113,file="transmission_kpoints.dat")   
     do ii = 1, nSteps
        write(113,advance='no',fmt='(f20.6)') (params%Emin + (ii-1)*params%Estep)*eV2Hartree
        do jj = 1, nK
          write(113,advance='no',fmt='(es20.8)') tunnMatSK(ii,jj)
        end do 
        write(113,*) 
     end do
     close(113)
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
  ! ----------------------------------------------------------------------------------

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
