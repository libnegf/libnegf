!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Greens functions.         !
!! Copyright (C) 2012 - 2026                                                !
!!                                                                          !
!! This file is part of libNEGF: a library for                              !
!! Non Equilibrium Green's Functions calculations                           !
!!                                                                          !
!! Developers: Alessandro Pecchia, Daniele Soccodato                        !
!! * This module was originally conceived by Sebastian Achilles (JFZ)
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

module self_energy
  use ln_precision, only: lp => dp
  use mat_def, only : x_DNS => z_DNS, create, destroy, assignment(=)
  use mpi_f08, MPI_XX_COMPLEX => MPI_DOUBLE_COMPLEX
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
    type(MPI_Comm), intent(in) :: comm2d
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
    integer :: Np, Mp, iQ, iK, iQglo, iQglo2, iKglo, nu
    integer :: iEminus, iEplus, iE, iEglo, iin
    integer :: myid, commsize, ndims = 2
    integer :: dims(2), coords(2), coordsH(2)
    type(x_DNS), target :: rbuff1, rbuff2, sbuffK, rbuffK
    type(x_DNS), pointer :: pbuff1, pbuff2, pGG1, pGG2, pSigma
    real(lp), allocatable :: KKbuf(:)
    integer :: ierr
    type(MPI_Request) :: rqE1,rqE2,rqE3,rqE4, rqK1, rqK2
    type(MPI_Status) :: statE(4), statK(2)

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
    call create(sbuffK, Np, Mp)
    call create(rbuffK, Np, Mp)

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
        !//       ^|------|E-----|--^                       ndiff = 2 iMinus = 5
        !//^|------|E-----|-^                               ndiff = 1 iMinus = 0
        !////////////////////////////////////////////////////////////////////////////////////
        nullify(pbuff1)
        nullify(pbuff2)
        !// pbuff1 => G(k,E-wq)
        !// checks if iE-iEhbaromega is on the same processor => no communication
        !//          for E<0 take G(E<0) = 0
        mdest = MPI_PROC_NULL
        msource = MPI_PROC_NULL
        iEminus = iE - iEhbaromega
        if (iEminus >= 0) then
          pbuff1 => GG(iEminus, iQ)%pMat
        else !(iEminus < 0) 
          rbuff1%val=(0.0_lp, 0.0_lp)  
          pbuff1 => rbuff1
        end if   

        ! -------------------------------------------------------------------------------
        ! MPI communications. 
        ! Communications take place if (iE < iEhbaromega .and. iEglo >= iEhbaromega)
        ! The latter is necessary to avoid 1-to-many communications
        ! If we take G(0) for all G(E<0) it is possible that rank 0 has to send
        ! G's to more than 1 rank, whenever ihbaromega > NEloc
        ! Example behaviour for ihbaromega=7 and NEloc=6, for iE=0:
        ! iEglo=0; ndiff=0 => no communication 
        ! iEglo=6; ihbaromega=7; NEloc=6 => ndiff=1; iEminus=0 
        ! => rank=0 SEND TO rank=1 GG(0,iQ)
        ! iEglo=12; ihbaromega=7; ndiff=12/6 - 5/6 = 2 => iEminus=(18-7)%6=5
        ! => rank=0 SEND TO rank=2 GG(5,iQ)   
        ! We have a 1 TO many communication
        ! -------------------------------------------------------------------------------
        if (dims(2) > 1) then
          if (iE < iEhbaromega) then
            if (iEglo - iEhbaromega < 0 ) then
              rbuff1%val=(0.0_lp, 0.0_lp)  
              pbuff1 => rbuff1
              ! takes care of the negative side where -1 turns into 0
              ! so +1 is added at the end if (mod(iEhbaromega-iEglo, NEloc) /= 0)
              !  ______________________
              !  333333222222221111111100000000111111112222222233333333
              ! --ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo----------
              !   345670123456701234567012345670123456701234567012345670
              !        ^-------------------^
              ! -----------------------------------------------------------------------------
              !rank=           0 iE           4 iEglo           4 iEminus         -16     -16
              !ndiff = (ihwq-iEglo)/NEloc + iEglo/NEloc = 16/8 + 0 = 2 
              !rank=           2 iE           4 iEglo          20 iEminus         -16      0
              !ndiff = iEglo/NEloc - (iEglo-ihwq)/NEloc = 20/8 - (0) = 2 
              ! -----------------------------------------------------------------------------
              !rank=           0 iE           5 iEglo           5 iEminus         -15     -15
              !ndiff = (ihwq-iEglo)/NEloc + iEglo/NEloc + 1 = 15/8 + 0 + 1 = 2 
              !rank=           2 iE           5 iEglo          21 iEminus         -15      1
              !ndiff = iEglo/NEloc - (iEglo-ihwq)/NEloc = 21/8 - 1/8  = 2 
              ! -----------------------------------------------------------------------------
              ndiff = (iEhbaromega-iEglo)/NEloc + iEglo/NEloc
              if (mod(iEhbaromega-iEglo, NEloc) /= 0) then
                ndiff = ndiff + 1
              end if
            else
              ! get processor distance between iEglo and iEglo-ihbarOmega
              ndiff = iEglo/NEloc - (iEglo-iEhbaromega)/NEloc
            end if
           
            iEminus = mod(iE + (ndiff+1)*NEloc - iEhbaromega, NEloc)

            call MPI_Cart_shift( comm2d, 1, ndiff, msource, mdest, ierr)
     
            if (mdest /= MPI_PROC_NULL) then 
               pGG1 => GG(iEminus, iQ)%pMat
               call MPI_Isend(pGG1%val, Mp*Np, MPI_XX_COMPLEX, mdest, 41, comm2d, rqE1, ierr)
            end if
            
            if (msource /= MPI_PROC_NULL) then
               call MPI_Irecv(rbuff1%val, Mp*Np, MPI_XX_COMPLEX, msource, 41, comm2d, rqE2, ierr)
               pbuff1 => rbuff1
            end if

            ! Wait for complete communications
            if (mdest/=MPI_PROC_NULL) call MPI_Wait(rqE1, statE(1), ierr);
            if (msource/=MPI_PROC_NULL) call MPI_Wait(rqE2, statE(2), ierr);
                 
          end if
        end if

        !////////////////////////////////////////////////////////////////////////////////////
        !// Communications of G(k, E+hwq)
        !////////////////////////////////////////////////////////////////////////////////////
        !// checks if iE+iEhbaromega is on the same processor => no communication
        !// pbuff2 => G(k,E+wq)
        pdest = MPI_PROC_NULL
        psource = MPI_PROC_NULL
        iEplus = iE + iEhbaromega
        if (iEplus < NEloc) then
          pbuff2 => GG(iEplus, iQ)%pMat
        else
          rbuff2%val=(0.0_lp, 0.0_lp)  
          pbuff2 => rbuff2
        end if

        if (dims(2) > 1) then
          if (iE + iEhbaromega >= NEloc) then
            if (iEglo + iEhbaromega >= NE ) then
              rbuff2%val=(0.0_lp, 0.0_lp)  
              pbuff2 => rbuff2
            end if
      
            ! get processor distance between iEglo+ihbarOmega and iEglo
            ndiff = (iEglo+iEhbaromega)/NEloc - iEglo/NEloc
            iEplus = mod(iEglo + iEhbaromega, NEloc)
            
            call MPI_Cart_shift( comm2d, 1, -ndiff, psource, pdest, ierr)
   
            if (pdest /= MPI_PROC_NULL) then
               pGG2 => GG(iEplus, iQ)%pMat
               call MPI_Isend(pGG2%val, Mp*Np, MPI_XX_COMPLEX, pdest, 42, comm2d, rqE3, ierr)
            end if
          
            if (psource /= MPI_PROC_NULL) then
               call MPI_Irecv(rbuff2%val, Mp*Np, MPI_XX_COMPLEX, psource, 42, comm2d, rqE4, ierr)
               pbuff2 => rbuff2
            end if

            if(pdest/=MPI_PROC_NULL) call MPI_Wait(rqE3, statE(3), ierr);
            if(psource/=MPI_PROC_NULL) call MPI_Wait(rqE4, statE(4), ierr);
               
          end if
        end if
        
        !//   Compute:
        !//   sbuffK_ij(iQ) = (fac_minus * GG_ij(iQ, E-wq) + fac_plus * GG_ij(iQ, E+wq))

        !$OMP PARALLEL DO private (nu)
        do nu = 1, Mp
          sbuffK%val(:,nu) = (fac_minus*pbuff1%val(:,nu) + fac_plus*pbuff2%val(:,nu))
        end do
        !$OMP END PARALLEL DO

        !// The Outer loop performs the summation over iQ
        !// This Inner loop updates Sigma(iK):
        !// Sigma_ij(iK, iE) = Sum_iQ   KK(|z_i-z_j|, iK, iQ) * sbuffK_ij(iQ)
        do iK = 0, NKloc-1
          iKglo = kindices_map(iK,coords(1))-1
          pSigma => Sigma(iE, iK)%pMat

          !$OMP PARALLEL private (nu, KKbuf)
          allocate(KKbuf(Np))
          !$OMP DO
          do nu = 1, Mp
            KKbuf(:) = KK(abs(izr(1:Np)-izc(nu)), iKglo, iQglo)
            pSigma%val(:,nu) = pSigma%val(:,nu) + KKbuf(:) * sbuffK%val(:,nu)
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
        !// cpu#0: (recv sbuffK(iQglo=4,5) from cpu#2  coordsH(1)=2) => iQglo2 = 4,5  ()
        !// cpu#1: (recv sbuffK(iQglo=0,1) from cpu#0  coordsH(1)=0) => iQglo2 = 0,1  ()
        !// cpu#2: (recv sbuffK(iQglo=2,3) from cpu#1  coordsH(1)=1) => iQglo2 = 2,3  ()
        !// ii = 2: #0->#2, #1->#0,  #2->#1
        !// cpu#0: (recv sbuffK(iQglo=2,3) from cpu#1  coordsH(1)=1) => iQglo2 = 2,3  ()
        !// cpu#1: (recv sbuffK(iQglo=4,5) from cpu#2  coordsH(1)=2) => iQglo2 = 4,5  ()
        !// cpu#2: (recv sbuffK(iQglo=0,1) from cpu#0  coordsH(1)=0) => iQglo2 = 0,1  ()
        !////////////////////////////////////////////////////////////////////////////////////
        if (dims(1) > 1) then
          do iin = 1, dims(1)-1
            ! Compute the communication points applying a shift in dimension 0 (k-grid)
            ! NOTE: this means a communication of 1->2,3,4,5... 2->3,4,5,6...
            !       maybe a nearest-neighbours communication 1->2->3->4 ... 
            !       is more efficient ? => sbuff = rbuff and send to neighbour
            ! First shift used to compute iQglo2:
            call MPI_Cart_shift(comm2d, 0, iin, hsource, hdest, ierr)
            call MPI_Cart_coords(comm2d, hsource, ndims, coordsH, ierr)
            iQglo2 = kindices_map(iQ,coordsH(1))-1

            ! Actual shift used in send/recv:
            call MPI_Cart_shift(comm2d, 0, 1, hsource, hdest, ierr)

            call MPI_Isend(sbuffK%val, Mp*Np, MPI_XX_COMPLEX, hdest, 43, comm2d, rqK1, ierr)
            call MPI_Irecv(rbuffK%val, Mp*Np, MPI_XX_COMPLEX, hsource, 43, comm2d, rqK2, ierr)

            call MPI_Wait(rqK2, statK(2), ierr)

            ! Sigma_ij(iK, iE) = Sum_iQ   KK(|z_i-z_j|, iK, iQ2) * rbuffK_ij(iQ2)
            do iK = 0, NKloc-1
              iKglo = kindices_map(iK,coords(1))-1
              pSigma => Sigma(iE, iK)%pMat
              !$OMP PARALLEL private (nu, KKbuf)
              allocate(KKbuf(Np))
              !$OMP DO
              do nu = 1, Mp
                KKbuf(:) = KK(abs(izr(1:Np)-izc(nu)), iKglo, iQglo2)
                pSigma%val(:,nu) = pSigma%val(:,nu) + KKbuf(:) * rbuffK%val(:,nu)
              end do
              !$OMP END DO
              deallocate(KKbuf)
              !$OMP END PARALLEL
            end do

            call MPI_Wait(rqK1, statK(1), ierr)

            !$OMP PARALLEL DO
            do nu = 1, Mp
              sbuffK%val(:,nu) = rbuffK%val(:,nu)
            end do
            !$OMP END PARALLEL DO

          end do
        end if

      end do eloop
    end do qloop

    call destroy(rbuff1)
    call destroy(rbuff2)
    call destroy(sbuffK)
    call destroy(rbuffK)


  end subroutine selfenergy


end module self_energy

