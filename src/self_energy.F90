module self_energy
  use ln_precision, only: lp => dp
  use mat_def, only : x_DNS => z_DNS, create, destroy, assignment(=)
  use mpi, MPI_XX_COMPLEX => MPI_DOUBLE_COMPLEX
  !use ln_precision, only: lp => sp
  !use mat_def, only : x_DNS => c_DNS, create, destroy, assignment(=)
  !use mpi, MPI_XX_COMPLEX => MPI_COMPLEX
  use mpi_globals
  implicit none
  private

  public :: selfenergy

  type, public :: TMatPointer
    type(x_DNS), pointer :: pMat
  end type TMatPointer

  contains

  subroutine selfenergy(comm2d, GG, fac_minus, fac_plus, iEhbaromega, &
                               &  izr, izc, KK, NK, NKloc, NE, NEloc, kindices_map, Sigma)
    integer, intent(in) :: comm2d
    type(TMatPointer), intent(in) :: GG(0:,0:)
    complex(lp), intent(in) :: fac_minus
    complex(lp), intent(in) :: fac_plus
    integer, intent(in) :: izr(:), izc(:)
    real(lp), intent(in) :: KK(0:,0:,0:)
    integer, intent(in) :: iEhbaromega
    integer, intent(in) :: NK
    integer, intent(in) :: NKloc
    integer, intent(in) :: NE
    integer, intent(in) :: NEloc
    integer, intent(in) :: kindices_map(0:,0:)
    type(TMatPointer), intent(in) :: Sigma(0:,0:)

    ! locals
    integer :: Np, Mp, iQ, iK, iQ2, iQglo, iQglo2, iKglo, mu, nu
    integer :: iEminus, iEplus, iE, iEglo, iin
    integer :: myid, commsize, ndims = 2
    integer :: dims(2), coords(2), coordsH(2)
    type(x_DNS), target :: rbuff1, rbuff2, sbuffH, rbuffH
    type(x_DNS), pointer :: pbuff1, pbuff2, pGG, pSigma
    real(lp), allocatable :: KKbuf(:)
    integer :: ierr
    integer :: rqE1,rqE2,rqE3,rqE4, rqH1, rqH2
    integer :: statusE(MPI_STATUS_SIZE,4), statusH(MPI_STATUS_SIZE,2)

    integer :: msource, mdest, psource, pdest
    integer :: hsource, hdest
    integer :: ndiff

    dims(1) = NK/NKloc
    dims(2) = NE/NEloc

    call MPI_comm_size(comm2d, commsize, ierr)
    call MPI_comm_rank(comm2d, myid, ierr)
    call MPI_cart_coords(comm2d, myid, ndims, coords, ierr)
    call MPI_barrier(comm2d, ierr)

    Np = size(GG(0,0)%pMat%val,1)
    Mp = size(GG(0,0)%pMat%val,2)
    call create(rbuff1, Np, Mp)
    call create(rbuff2, Np, Mp)
    call create(sbuffH, Np, Mp)
    call create(rbuffH, Np, Mp)

    ! Sigma_ij(iK, iE) = Sum_iQ   KK(|z_i-z_j|, iK, iQ) *
    !               * (fac_minus * GG_ij(iQ, E-wq) + fac_plus * GG_ij(iQ, E+wq))
    qloop:do iQ = 0, NKloc-1
      iQglo = kindices_map(iQ,coords(1))-1
      eloop:do iE = 0, NEloc-1
        iEglo = iE + coords(2)*NEloc
        !////////////////////////////////////////////////////////////////////////////////////
        !// Communications of G(k, E-hwq)
        !////////////////////////////////////////////////////////////////////////////////////
        !//  012345 012345 012345 012345    iE < NEloc = 6
        !//  -------------------------------------------------------------------------------
        !//  000000 000011 111111 112222
        !//  012345 678901 234567 890123    iEglo = 1 + 6 = 7
        !// |oooooo|oooooo|oooooo|oooooo|
        !//  ^-----|-E----|--^              ihbarOmega = 7; ndiff = 1 iMinus = 0
        !//     ^--|----E-|-----^                           ndiff = 1 iMinus = 3
        !////////////////////////////////////////////////////////////////////////////////////

        !// pbuff1 points to G(k,E-wq)
        !// checks if iE-iEhbaromega is on the same processor => no communication
        !// if iEglo<iEhbaromega => the processor is on the lower end of the energy grid
        !//                      => communication is local or forced to be by truncation of G
        !//                         such that G(E<0) = G(0), e.g., G(E-hw) = G(Elow)
        iEminus = iE - iEhbaromega
        if (iEminus >= 0 .or. iEglo < iEhbaromega) then
          if (iEglo < iEhbaromega) iEminus = 0
          pbuff1 => GG(iEminus, iQ)%pMat
          mdest = MPI_PROC_NULL
          msource = MPI_PROC_NULL
        end if

        ! MPI communications
        if (dims(2) > 1) then
          if ( iEhbaromega >= NEloc ) then
            ndiff = iEhbaromega / NEloc
          else
            ndiff = 1
          end if

          call MPI_Cart_shift( comm2d, 1, ndiff, msource, mdest, ierr)

          if (mdest /= MPI_PROC_NULL .and. iE < iEhbaromega) then
             iEminus = mod(iE + (ndiff+1)*NEloc - iEhbaromega, NEloc)
             pGG => GG(iEminus, iQ)%pMat
             call MPI_Isend(pGG%val, Mp*Np, MPI_XX_COMPLEX, mdest, 41, comm2d, rqE1, ierr)
          end if

          if (msource /= MPI_PROC_NULL .and. iE < iEhbaromega) then
             call MPI_Irecv(rbuff1%val, Mp*Np, MPI_XX_COMPLEX, msource, 41, comm2d, rqE2, ierr)
             pbuff1 => rbuff1
          end if
        end if

        !////////////////////////////////////////////////////////////////////////////////////
        !// Communications of G(k, E+hwq)
        !////////////////////////////////////////////////////////////////////////////////////
        !// checks if iE+iEhbaromega is on the same processor => no communication
        !// pbuff2 points to G(k,E+wq)
        iEplus = iE + iEhbaromega
        if (iEplus < NEloc .or. iEglo+iEhbaromega > NE) then
          if (iEglo+iEhbaromega > NE) iEplus = NEloc-1
          pbuff2 => GG(iEplus, iQ)%pMat
          pdest = MPI_PROC_NULL
          psource = MPI_PROC_NULL
        end if

        if (dims(2) > 1) then
          if ( iEhbaromega >= NEloc ) then
            ndiff = iEhbaromega / NEloc
          else
            ndiff = 1
          end if

          call MPI_Cart_shift( comm2d, 1, -ndiff, psource, pdest, ierr)

          if (pdest /= MPI_PROC_NULL .and. iE >= NEloc - iEhbaromega) then
             iEplus = mod(iEglo + iEhbaromega, NEloc)
             pGG => GG(iEplus, iQ)%pMat
             call MPI_Isend(pGG%val, Mp*Np, MPI_XX_COMPLEX, pdest, 42, comm2d, rqE3, ierr)
          end if

          if (psource /= MPI_PROC_NULL .and. iE >= NEloc - iEhbaromega) then
             call MPI_Irecv(rbuff2%val, Mp*Np, MPI_XX_COMPLEX, psource, 42, comm2d, rqE4, ierr)
             pbuff2 => rbuff2
          end if
        end if

        ! Wait for complete communications
        if (dims(2) > 1) then
          if(mdest /= MPI_PROC_NULL .and. iE<iEhbaromega) call MPI_Wait(rqE1, statusE(:,1), ierr);
          if(msource /= MPI_PROC_NULL  .and. iE<iEhbaromega) call MPI_Wait(rqE2, statusE(:,2), ierr);

          if(pdest /= MPI_PROC_NULL .and. iE >= NEloc-iEhbaromega) call MPI_Wait(rqE3, statusE(:,3), ierr);
          if(psource /= MPI_PROC_NULL .and. iE >= NEloc-iEhbaromega) call MPI_Wait(rqE4, statusE(:,4), ierr);
        end if

        !//   Compute:
        !//   sbuffH_ij(iQ) = (fac_minus * GG_ij(iQ, E-wq) + fac_plus * GG_ij(iQ, E+wq))
        !$OMP PARALLEL DO private (nu)
        do nu = 1, Mp
          sbuffH%val(:,nu) = (fac_minus*pbuff1%val(:,nu) + fac_plus*pbuff2%val(:,nu))
        end do
        !$OMP END PARALLEL DO

        !// The Outer loop performs the summation over iQ
        !// This Inner loop updates Sigma(iK):
        !// Sigma_ij(iK, iE) = Sum_iQ   KK(|z_i-z_j|, iK, iQ) * sbuffH_ij(iQ)
        do iK = 0, NKloc-1
          iKglo = kindices_map(iK,coords(1))-1
          pSigma => Sigma(iE, iK)%pMat

          !$OMP PARALLEL private (nu, KKbuf)
          allocate(KKbuf(Np))
          !$OMP DO
          do nu = 1, Mp
            KKbuf(:) = KK(abs(izr(1:Np)-izc(nu)), iKglo, iQglo)
            pSigma%val(:,nu) = pSigma%val(:,nu) + KKbuf(:) * sbuffH%val(:,nu)
          end do
          !$OMP END DO
          deallocate(KKbuf)
          !$OMP END PARALLEL
        end do


        !////////////////////////////////////////////////////////////////////////////////////
        !// Communications over the k-grid (Note k-grid is periodic)
        !// Suppose iQglo = 0, 1, 2, 3, 4, 5
        !// dims(1) == 3 =>
        !// cpu#0: iQ = 0, 1; coords(1)=0 => iQglo = 0, 1
        !// cpu#1: iQ = 0, 1; coords(1)=1 => iQglo = 2, 3
        !// cpu#2: iQ = 0, 1; coords(1)=2 => iQglo = 4, 5
        !//
        !// ii = 1: #0->#1, #1->#2,  #2->#0
        !// cpu#0: (recv sbuffH(iQglo=4,5) from cpu#2  coordsH(1)=2) => iQglo2 = 4,5  ()
        !// cpu#1: (recv sbuffH(iQglo=0,1) from cpu#0  coordsH(1)=0) => iQglo2 = 0,1  ()
        !// cpu#2: (recv sbuffH(iQglo=2,3) from cpu#1  coordsH(1)=1) => iQglo2 = 2,3  ()
        !// ii = 2: #0->#2, #1->#0,  #2->#1
        !// cpu#0: (recv sbuffH(iQglo=2,3) from cpu#1  coordsH(1)=1) => iQglo2 = 2,3  ()
        !// cpu#1: (recv sbuffH(iQglo=4,5) from cpu#2  coordsH(1)=2) => iQglo2 = 4,5  ()
        !// cpu#2: (recv sbuffH(iQglo=0,1) from cpu#0  coordsH(1)=0) => iQglo2 = 0,1  ()
        !////////////////////////////////////////////////////////////////////////////////////
        if (dims(1) > 1) then
          do iin = 1, dims(1)-1
            call MPI_Cart_shift(comm2d, 0, iin, hsource, hdest, ierr)
            call MPI_Cart_coords(comm2d, hsource, ndims, coordsH, ierr)

            iQglo2 = kindices_map(iQ,coordsH(1))-1

            call MPI_Isend(sbuffH%val, Mp*Np, MPI_XX_COMPLEX, hdest, 43, comm2d, rqH1, ierr)
            call MPI_Irecv(rbuffH%val, Mp*Np, MPI_XX_COMPLEX, hsource, 43, comm2d, rqH2, ierr)

            call MPI_Wait(rqH2, statusH(:,2), ierr)

            ! Sigma_ij(iK, iE) = Sum_iQ   KK(|z_i-z_j|, iK, iQ2) * rbuffH_ij(iQ2)
            do iK = 0, NKloc-1
              iKglo = kindices_map(iK,coords(1))-1
              pSigma => Sigma(iE, iK)%pMat
              !$OMP PARALLEL private (nu, KKbuf)
              allocate(KKbuf(Np))
              !$OMP DO
              do nu = 1, Mp
                KKbuf(:) = KK(abs(izr(1:Np)-izc(nu)), iKglo, iQglo2)
                pSigma%val(:,nu) = pSigma%val(:,nu) + KKbuf(:) * rbuffH%val(:,nu)
              end do
              !$OMP END DO
              deallocate(KKbuf)
              !$OMP END PARALLEL
            end do

            call MPI_Wait(rqH1, statusH(:,1), ierr)
          end do
        end if

      end do eloop
    end do qloop

    call destroy(rbuff1)
    call destroy(rbuff2)
    call destroy(sbuffH)
    call destroy(rbuffH)


  end subroutine selfenergy


end module self_energy

