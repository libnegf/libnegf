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

module iterative_cpu
  use ln_precision
  use ln_constants, only : pi, minusOne
  use ln_allocation
  use mat_def
  use sparsekit_drv
  use inversions
  use ln_enums 
  use ln_structure, only : TStruct_Info
  use lib_param, only : MAXNCONT, Tnegf, intarray
  use interactions, only : TInteraction, TInteractionList, TInteractionNode

  implicit none
  private
  public :: calculate_gsmr_blocks
  public :: calculate_Gr_tridiag_blocks
  public :: calculate_Gn_tridiag_blocks
  public :: calculate_single_transmission_2_contacts
  public :: calculate_single_transmission_N_contacts
  public :: check_convergence_trid

contains

  !***********************************************************************
  !
  !  g_small right (gsmr) calculation - write on memory
  !
  !***********************************************************************

  subroutine calculate_gsmr_blocks(negf,ESH,sbl,ebl,gsmr,keepall)

    !***********************************************************************
    !Input:
    !ESH: sparse matrices array ESH(nbl,nbl)
    !
    !nbl (number of layers)  global needed   (indblk not anymore)
    !
    !Output:
    !sparse matrices array global variable gsmr(nbl) is available in memory
    !single blocks are allocated internally, array Gr(nbl,nbl)
    !must be allocated externally
    !***********************************************************************

    implicit none

    !In/Out
    type(Tnegf), intent(in) :: negf
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl,ebl                 ! start block, end block
    type(z_DNS), dimension(:), intent(inout) :: gsmr
    logical, intent(in), optional :: keepall

    !Work
    !type(z_DNS), dimension(:,:), allocatable :: INV
    type(z_DNS) :: work1, work2
    integer :: nrow, M, N
    integer :: i, nbl
    logical :: keep

    if (sbl.lt.ebl) return

    nbl = size(ESH,1)

    if (nbl.eq.1) return

    keep = .true.
    if (present(keepall)) then
       keep = keepall
    end if

    nrow=ESH(sbl,sbl)%nrow

    call create(gsmr(sbl),nrow,nrow)

    call compGreen(gsmr(sbl),ESH(sbl,sbl),nrow)


    do i=sbl-1,ebl,-1

       call prealloc_mult(ESH(i,i+1),gsmr(i+1),minusOne,work1)

       if (.not.keep) then
          call destroy(gsmr(i+1))
       end if

       call prealloc_mult(work1,ESH(i+1,i),work2)

       call destroy(work1)

       call prealloc_sum(ESH(i,i),work2,work1)

       call destroy(work2)

       call create(gsmr(i),work1%nrow,work1%nrow)

       call compGreen(gsmr(i),work1,work1%nrow)

       call destroy(work1)

    end do


  end subroutine calculate_gsmr_blocks


  !***********************************************************************
  !
  !  Diagonal, Subdiagonal, Superdiagonal blocks of Green Retarded
  !  Gr(nbl,nbl) - writing on memory
  !
  !***********************************************************************

  subroutine calculate_Gr_tridiag_blocks(negf,ESH,gsmr,Gr,sbl,ebl)

    !***********************************************************************
    !Input:
    !ESH: dense matrices array ESH(nbl,nbl)
    !sbl, ebl : block indexes
    ! If only sbl is specified, it calculates Gr(sbl, sbl)
    ! If sbl > ebl, it calculates Gr(ebl:sbl, ebl:sbl), Gr(ebl:sbl + 1, ebl:sbl),
    !    Gr(ebl:sbl, ebl:sbl + 1) (need gsml)
    ! If sbl < ebl, it calculates Gr(sbl:ebl, sbl:ebl), Gr(sbl:ebl - 1, sbl:ebl),
    !    Gr(sbl:ebl, sbl:ebl - 1) (need gsmr)
    !
    !
    !Output:
    !sparse matrices array global variable Gr(nbl,nbl) is available in
    !memory - single blocks are allocated internally, array Gr(nbl,nbl)
    !must be allocated externally
    !***********************************************************************

    implicit none

    !In/Out
    type(Tnegf), intent(in) :: negf
    type(z_DNS), dimension(:), intent(in) :: gsmr
    type(z_DNS), dimension(:,:), intent(inout) :: Gr
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl
    integer, intent(in), optional :: ebl

    !Work
    integer :: i,nrow,nbl
    type(z_DNS) :: work1, work2, work3

    nbl = size(ESH,1)

    if (sbl.gt.nbl) return
    if (sbl.lt.1) return

    if (.not.present(ebl)) then
       if (nbl.eq.1) then
          nrow = ESH(sbl,sbl)%nrow
          call create(Gr(sbl,sbl),nrow,nrow)
          call compGreen(Gr(sbl,sbl),ESH(sbl,sbl),nrow)
       else
          nrow = ESH(sbl,sbl)%nrow
          call create(work1,nrow,nrow)
          work1%val = ESH(sbl,sbl)%val
          if (sbl+1.le.nbl) then
             call prealloc_mult(ESH(sbl,sbl+1),gsmr(sbl+1),work2)
             call prealloc_mult(work2,ESH(sbl+1,sbl),work3)
             call destroy(work2)
             call prealloc_sum(work1,work3,minusOne,work2)
             call destroy(work3)
             work1%val = work2%val
             call destroy(work2)
          endif
          if (sbl-1.ge.1) then
             error stop "Error: Gr_tridiag requires gsml"
          endif

          call create(Gr(sbl,sbl),nrow,nrow)
          call compGreen(Gr(sbl,sbl),work1,nrow)
          call destroy(work1)
       endif
       return
    endif


    !***
    !Diagonal, Subdiagonal and Superdiagonal blocks
    !***
    if ((ebl.ge.sbl).and.(ebl.gt.1).and.(sbl.gt.1)) THEN
       do i=sbl,ebl,1
          call prealloc_mult(gsmr(i),ESH(i,i-1),work1)
          call prealloc_mult(work1,Gr(i-1,i-1),minusOne,Gr(i,i-1))
          call destroy(work1)

          call prealloc_mult(ESH(i-1,i),gsmr(i),work2)
          call prealloc_mult(Gr(i-1,i-1),work2,minusOne,Gr(i-1,i))

          call prealloc_mult(Gr(i,i-1),work2,minusOne,work1)
          call destroy(work2)

          call prealloc_sum(gsmr(i),work1,Gr(i,i))
          call destroy(work1)
       end do
    ELSE
       error stop "Error: Gr_tridiag requires gsml"
    endif

  end subroutine calculate_Gr_tridiag_blocks

  ! Implements a new algorithm based on an iterative scheme to solve
  ! [ES - H - SigmaR] Gr = Sigma< Ga
  ! The subroutine uses the gsmr(:) computed before and makes an iteration
  ! upward to build gsmn and then downward to build the 3-diagonal blocks
  subroutine calculate_Gn_tridiag_blocks(negf,ESH,SelfEneR,frm,ref,struct,gsmr,Gr,Gn)
    type(TNegf), intent(in) :: negf
    type(z_DNS), dimension(:,:), intent(in) :: ESH, Gr
    type(z_DNS), dimension(:), intent(in) :: SelfEneR, gsmr
    real(dp), dimension(:), intent(in) :: frm
    integer, intent(in) :: ref
    type(Tstruct_info), intent(in) :: struct
    type(z_DNS), dimension(:,:), intent(inout) :: Gn

    !Work
    type(z_DNS), dimension(:,:), allocatable :: Sigma_n
    type(z_DNS) :: work1, Ga, Gam
    complex(dp) :: frmdiff
    integer :: i, j, fp
    integer :: nbl, ncont, cb

    ncont = struct%num_conts
    nbl = struct%num_PLs

    !build Sigma_n from SelfEneR
    call allocate_blk_dns(Sigma_n, nbl)
    call init_tridiag_blk(Sigma_n, ESH)

    call add_sigma_n(negf, Sigma_n)

    ! Add contact self-energies
    do j=1,ncont
      frmdiff = frm(j) - frm(ref)
      if (j.NE.ref .AND. ABS(frmdiff).GT.EPS) THEN
        cb=struct%cblk(j) ! block corresponding to contact j
        call zspectral(SelfEneR(j),SelfEneR(j),0,Gam)
        Sigma_n(cb,cb)%val = Sigma_n(cb,cb)%val + frmdiff*Gam%val
        call destroy(Gam)
      endif
    end do

    call calculate_sigma_n()

    call zdagger(Gr(1,1), Ga)
    call prealloc_mult(Sigma_n(1,1), Ga, work1)
    call prealloc_mult(Gr(1,1), work1, Gn(1,1))
    call destroy(Ga, work1)

    if (nbl .eq. 1) then
       call destroy_tridiag_blk(Sigma_n)
       call deallocate_blk_dns(Sigma_n)
       return
    end if
    !Explicit formulae:
    !Gn(i+1,i) = gsmr(i+1)*[Sigma(i+1,i)Ga(i,i) + Sigma(i+1,i+1)Ga(i+1,i) - Tr(i+1,i)Gn(i,i)]
    !Gn(i,i+1) = [Gr(i,i)Sigma(i,i+1) + Gr(i,i+1)Sigma(i+1,i+1) - Gn(i,i)Ta(i,i+1)] * gsma(i+1)
    !Use Hermitian property of Gn:
    !Gn(i,i+1) = Gn(i+1,i)^dag
    !Gn(i+1,i+1) = gsmr(i+1) * [Sigma(i+1,i)Ga(i,i+1) + Sigma(i+1,i+1)Ga(i+1,i+1) - Tr(i+1,i)Gn(i,i+1)]
    !Implementation exploits cumulative sum of prealloc_mult, C = C + A*B

    do i = 1, nbl-1

        call zdagger(Gr(i,i), Ga)
        call prealloc_mult(Sigma_n(i+1,i), Ga, work1)
        call destroy(Ga)

        call zdagger(Gr(i,i+1), Ga)
        call prealloc_mult(Sigma_n(i+1,i+1), Ga, work1)
        call destroy(Ga)

        call prealloc_mult(ESH(i+1,i), Gn(i,i), minusOne, work1)

        call prealloc_mult(gsmr(i+1), work1, Gn(i+1,i))
        call destroy(work1)

        call destroy(Gn(i,i+1))
        call zdagger(Gn(i+1,i), Gn(i,i+1))

        call zdagger(Gr(i+1,i), Ga)
        call prealloc_mult(Sigma_n(i+1,i), Ga, work1)
        call destroy(Ga)

        call zdagger(Gr(i+1,i+1), Ga)
        call prealloc_mult(Sigma_n(i+1,i+1), Ga, work1)
        call destroy(Ga)

        call prealloc_mult(ESH(i+1,i), Gn(i,i+1), minusOne, work1)

        call prealloc_mult(gsmr(i+1), work1, Gn(i+1,i+1))
        call destroy(work1)

    end do

    call destroy_tridiag_blk(Sigma_n)
    call deallocate_blk_dns(Sigma_n)

    contains
    ! Recursive calculation of Sigma_n:
    ! gns(i+1) = gsmr(i+1) Sigma(i+1,i+1) gsmr(i+1)^dag
    ! Sigma(i,i) = Sigma(i,i) + Tr(i,i+1) gns(i+1) Ta(i+1,i)
    !                         - Tr(i,i+1) gsmr(i+1) Sigma(i+1,i)
    !                         - Sigma(i,i+1) gsmr^dag(i+1) Ta(i+1,i)]
    !
    subroutine calculate_sigma_n()
      !Work
      type(z_DNS) :: work, gns, gsmrDag, ESHdag

      ! if nbl = 1 => Sigma_n(1,1) is ready
      if (nbl.eq.1) return
      !g^n(nbl) = gsmr(nbl) Sigma(nbl,nbl) gsma(nbl)
      call zdagger(gsmr(nbl),gsmrDag)
      call prealloc_mult(gsmr(nbl), Sigma_n(nbl,nbl), work1)
      call prealloc_mult(work1, gsmrDag, gns)
      call destroy(gsmrDag, work1)

      do i = nbl-1, 1, -1
        !work1 = Tr(i,i+1) gns(i+1) Ta(i+1,i)
        ! Tr(i,i+1) = ESH(i,i+1);  Ta(i+1,i) = ESH(i,i+1)^dag
        call zdagger(ESH(i,i+1), ESHdag)
        call prealloc_mult(ESH(i,i+1), gns, work)
        call prealloc_mult(work, ESHdag, Sigma_n(i,i))
        call destroy(work, gns)

        !work2 = Sigma(i,i+1) gsmr^dag(i+1) Ta(i+1,i)
        call zdagger(gsmr(i+1), gsmrDag)
        call prealloc_mult(Sigma_n(i,i+1), gsmrDag, work)
        call prealloc_mult(work, ESHdag, minusone, Sigma_n(i,i))
        call destroy(work, gsmrDag, ESHdag)

        !work3 = ESH(i,i+1) gsmr(i+1) Sigma(i+1,i)
        call prealloc_mult(ESH(i,i+1), gsmr(i+1), work)
        call prealloc_mult(work, Sigma_n(i+1,i), minusone, Sigma_n(i,i))
        call destroy(work)

        if (i > 1) then
          !gns(i) = gsmr(i) * Sigma_n(i,i) * gsmr^dag(i)
          call zdagger(gsmr(i), gsmrDag)
          call prealloc_mult(gsmr(i), Sigma_n(i,i), work)
          call prealloc_mult(work, gsmrDag, gns)
          call destroy(work, gsmrDag)
        end if

      end do

    end subroutine calculate_sigma_n

  end subroutine calculate_Gn_tridiag_blocks

  !--------------------------------------------------------------------------
  ! Add all Self-enrgies to Sigma_n
  subroutine add_sigma_n(negf, sigma_n)
    class(TNegf) :: negf
    type(z_DNS) :: sigma_n(:,:)

    type(TInteractionNode), pointer :: it
    it => negf%interactList%first

    do while (associated(it))
      call it%inter%add_sigma_n(sigma_n, negf%iEloc, negf%iKloc, negf%spin)
      it => it%next
    end do
  end subroutine add_sigma_n
  !************************************************************************
  !
  ! Subroutine for transmission calculation
  !
  !************************************************************************
  ! NOTE:
  !
  !  This subroutine was hacked quickly to obain effiecient tunneling calcs
  !  Useful only when there are 2 contacts
  !                ===================
  !************************************************************************
  subroutine calculate_single_transmission_2_contacts(negf,ni,nf,ESH,SelfEneR,cblk,tun_proj,Gr,TUN)
    type(Tnegf), intent(in) :: negf
    integer, intent(in) :: ni,nf
    type(z_DNS), intent(in) :: SelfEneR(MAXNCONT)
    type(z_DNS), intent(in) :: ESH(:,:)
    integer, intent(in) :: cblk(:)
    type(intarray), intent(in) :: tun_proj
    type(z_DNS), dimension(:,:), intent(in) :: Gr
    real(dp), intent(out) :: TUN

    !Work variables
    Integer :: ct1, bl1
    logical, dimension(:), allocatable :: tun_mask
    Type(z_DNS) :: work1, work2, GAM1_dns, GA, TRS, AA
    Complex(dp), parameter ::    j = (0.0_dp,1.0_dp)  ! CMPX unity

    if (size(cblk).gt.2) then
       write(*,*) "ERROR: calculate_single_transmission_2_contacts is valid only for 2 contacts"
       TUN = 0.0_dp
       return
    endif

    !Arrange contacts in way that order between first and second is always the
    !same (always ct1 < ct2)

    if (cblk(ni).lt.cblk(nf)) then
       ct1=ni;
    else
       ct1=nf;
    endif

    bl1=cblk(ct1);

    call zdagger(Gr(bl1,bl1),GA)

    ! Computes the Gamma matrices
    call zspectral(SelfEneR(ct1),SelfEneR(ct1),0,GAM1_dns)

    ! Work to compute transmission matrix (Gamma G Gamma G)
    call prealloc_mult(GAM1_dns,Gr(bl1,bl1),work1)

    call prealloc_mult(work1,GAM1_dns,work2)

    call destroy(work1)

    call prealloc_mult(work2,GA,work1)

    call destroy(work2)

    call create(AA,GA%nrow,GA%ncol)

    AA%val = j * (Gr(bl1,bl1)%val-GA%val)

    call destroy(GA)

    call prealloc_mult(GAM1_dns,AA,work2)

    call destroy(GAM1_dns,AA)

    call create(TRS,work1%nrow,work1%ncol)

    TRS%val = work2%val - work1%val

    call get_tun_mask(ESH, bl1, tun_proj, tun_mask)

    TUN = abs( real(trace(TRS, tun_mask)) )

    call log_deallocate(tun_mask)

    call destroy(TRS,work1,work2)

  end subroutine calculate_single_transmission_2_contacts

  !************************************************************************
  !
  ! Subroutine for transmission calculation (GENERIC FOR N CONTACTS)
  !
  !************************************************************************
  subroutine calculate_single_transmission_N_contacts(negf,ni,nf,ESH,SelfEneR,cblk,tun_proj,gsmr,Gr,TUN)
    type(Tnegf), intent(in) :: negf
    integer, intent(in) :: ni,nf
    type(z_DNS), intent(in) :: SelfEneR(MAXNCONT)
    type(z_DNS), intent(in) :: ESH(:,:)
    type(z_DNS), dimension(:),intent(in) :: gsmr
    type(z_DNS), dimension(:,:),intent(in) :: Gr
    integer, intent(in) :: cblk(:)
    type(intarray), intent(in) :: tun_proj
    real(dp), intent(out) :: TUN

    !Work variables
    Integer :: ct1, ct2, bl1, bl2, i, nbl
    logical, dimension(:), allocatable :: tun_mask
    Type(z_DNS) :: work1, work2, GAM1_dns, GAM2_dns, GA, TRS
    Real(kind=dp) :: max

    !Arrange contacts in way that order between first and second is always the
    !same (always ct1 < ct2)

    if (cblk(ni).lt.cblk(nf)) then
       ct1=ni;ct2=nf;
    else
       ct1=nf;ct2=ni;
    endif

    bl1=cblk(ct1); bl2=cblk(ct2);
    nbl = size(cblk)
    ! in this way nt1 < nt2 by construction

    if ( nbl.gt.1 .and. (bl2-bl1).gt.1) then

       ! Compute column-blocks of Gr(i,bl1) up to i=bl2
       ! Gr(i,bl1) = -gr(i) T(i,i-1) Gr(i-1,bl1)
       do i = bl1+1, bl2
          !Checks whether previous block is non null.
          !If so next block is also null => TUN = 0
          max=maxval(abs(Gr(i-1,bl1)%val))

          if (max.lt.EPS) then
             TUN = EPS*EPS !for log plots
             !Destroy also the block adjecent to diagonal since
             !this is not deallocated anymore in calling subroutine
             if (i.gt.(bl1+1)) call destroy(Gr(i-1,bl1))
             return
          endif

          !Checks whether block has been created, if not do it
          if (.not.allocated(Gr(i,bl1)%val)) then

             call prealloc_mult(gsmr(i),ESH(i,i-1),(-1.0_dp, 0.0_dp),work1)

             call prealloc_mult(work1,Gr(i-1,bl1),Gr(i,bl1))

             call destroy(work1)

          endif

          ! avoid destroying blocks closer to diagonal
          if (i.gt.(bl1+2)) call destroy(Gr(i-1,bl1))

       end do

    endif
    ! Computes the Gamma matrices
    call zspectral(SelfEneR(ct1),SelfEneR(ct1),0,GAM1_dns)
    call zspectral(SelfEneR(ct2),SelfEneR(ct2),0,GAM2_dns)

         write(*,*) '~-~-~-~-~-~-~-~-~-~-~-~-~-~-'
         write(*,*) 'N_conts: TRS= GAM2 * Gr(bl2,bl1)* GAM1 * Gr(bl2,bl1)^+'
         write(*,*) 'N_conts: sum_GAM1_dns=', sum(ABS(GAM1_dns%val))
         write(*,*) 'N_conts: sum_Gr(',bl2,bl1,')=', sum(ABS(Gr(bl2,bl1)%val))
         write(*,*) ''

    ! Work to compute transmission matrix (Gamma2 Gr Gamma1 Ga)
    call prealloc_mult(GAM2_dns,Gr(bl2,bl1),work1)
         write(*,*) 'N_conts: sum_GAM2_dns=', sum(ABS(GAM2_dns%val))

    call destroy(GAM2_dns)

    call prealloc_mult(work1,GAM1_dns,work2)

         write(*,*) 'N_conts: sum_work1=', sum(ABS(work1%val))
         write(*,*) 'N_conts: sum_work2', sum(ABS(work2%val))
    call destroy(work1)

    call destroy(GAM1_dns)

    call zdagger(Gr(bl2,bl1),GA)

    if (bl2.gt.bl1+1) call destroy( Gr(bl2,bl1) )

    call prealloc_mult(work2,GA,TRS)

    call destroy(work2)

    call destroy(GA)

    call get_tun_mask(ESH, bl2, tun_proj, tun_mask)

    TUN = abs( real(trace(TRS, tun_mask)) )

    call log_deallocate(tun_mask)

    call destroy(TRS)

  end subroutine calculate_single_transmission_N_contacts

  ! Based on projection indices build a logical mask just on contact block
  subroutine get_tun_mask(ESH,nbl,tun_proj,tun_mask)
    Type(z_DNS), intent(in) :: ESH(:,:)
    integer, intent(in) :: nbl
    type(intarray), intent(in) :: tun_proj
    logical, intent(out), allocatable :: tun_mask(:)

    integer :: ii, istart, iend, ind

    call log_allocate(tun_mask, ESH(nbl,nbl)%nrow)

    if (allocated(tun_proj%indexes)) then
       tun_mask = .false.

       ! set the start/end indices of nbl
       ! NB: istart has offset -1 to avoid +/-1 operations
       istart = 0
       do ii = 1, nbl-1
          istart = istart + ESH(ii,ii)%nrow
       end do
       iend = istart + ESH(nbl,nbl)%nrow + 1

       ! select the indices in tun_proj
       do ii = 1, size(tun_proj%indexes)
          ind = tun_proj%indexes(ii)
          if (ind > istart .and. ind < iend) then
             tun_mask(ind - istart) = .true.
          end if
       end do
    else
       tun_mask = .true.
    end if

  end subroutine get_tun_mask

  subroutine check_convergence_trid(negf,T,nome,gpu)
    type(Tnegf), intent(in) :: negf
    type(z_DNS), dimension(:,:), intent(in) :: T
    character(3), intent(in) :: nome
    logical, intent(in) :: gpu

    integer :: nbl, i
    real(dp) :: summ

    nbl = size(T,1)

    if (gpu) then
       write(*,*) '~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-'
    else
       write(*,*) '~-~-~-~-',nome,' check convergence CPU: ~-~-~-~-'
       write(*,*) '       ',nome,'(',1,1,')=', sum(ABS(T(1,1)%val))

       do i= 2,nbl
          write(*,*) '       ',nome,'(',i,i,')=', sum(ABS(T(i,i)%val))

          write(*,*) '       ',nome,'(',i,i-1,')=', sum(ABS(T(i,i-1)%val))

          write(*,*) '       ',nome,'(',i-1,i,')=', sum(ABS(T(i-1,i)%val))
       end do
       write(*,*) '~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-'
    endif
  end subroutine check_convergence_trid

  !---------------------------------------------------

  subroutine allocate_blk_dns(blkM,nbl)
    type(z_DNS), dimension(:,:), allocatable :: blkM
    integer :: nbl, ierr

    allocate(blkM(nbl,nbl),stat=ierr)
    if (ierr.ne.0) error stop 'ALLOCATION ERROR: could not allocate block-Matrix'

  end subroutine allocate_blk_dns

  !---------------------------------------------------

  subroutine deallocate_blk_dns(blkM)
    type(z_DNS), dimension(:,:), allocatable :: blkM
    integer :: ierr

    deallocate(blkM,stat=ierr)
    if (ierr.ne.0) error stop 'DEALLOCATION ERROR: could not deallocate block-Matrix'

  end subroutine deallocate_blk_dns

  !---------------------------------------------------

  subroutine init_tridiag_blk(Matrix,S)
    type(z_DNS), dimension(:,:) :: Matrix,S

    integer :: nbl, j

    nbl = SIZE(Matrix,1)

    call create(Matrix(1,1),S(1,1)%nrow,S(1,1)%ncol)
    Matrix(1,1)%val=(0.0_dp,0.0_dp)
    do j=2,nbl-1
       call create(Matrix(j-1,j),S(j-1,j)%nrow,S(j-1,j)%ncol)
       Matrix(j-1,j)%val=(0.0_dp,0.0_dp)
       call create(Matrix(j,j),S(j,j)%nrow,S(j,j)%ncol)
       Matrix(j,j)%val=(0.0_dp,0.0_dp)
       call create(Matrix(j,j-1),S(j,j-1)%nrow,S(j,j-1)%ncol)
       Matrix(j,j-1)%val=(0.0_dp,0.0_dp)
    end do
    if (nbl.gt.1) then
       call create(Matrix(nbl,nbl),S(nbl,nbl)%nrow,S(nbl,nbl)%ncol)
       Matrix(nbl,nbl)%val=(0.0_dp,0.0_dp)
       call create(Matrix(nbl-1,nbl),S(nbl-1,nbl)%nrow,S(nbl-1,nbl)%ncol)
       Matrix(nbl-1,nbl)%val=(0.0_dp,0.0_dp)
       call create(Matrix(nbl,nbl-1),S(nbl,nbl-1)%nrow,S(nbl,nbl-1)%ncol)
       Matrix(nbl,nbl-1)%val=(0.0_dp,0.0_dp)
    endif

  end subroutine init_tridiag_blk

  !---------------------------------------------------

  subroutine destroy_tridiag_blk(M)
    type(z_DNS), dimension(:,:), allocatable :: M

    integer :: i, nbl

    if (.not.allocated(M)) return

    nbl=size(M,1)

    do i=1,nbl
      if (allocated(M(i,i)%val)) then
        call destroy(M(i,i))
      end if
    end do
    do i=2,nbl
      if (allocated(M(i-1,i)%val)) then
        call destroy(M(i-1,i))
      end if
      if (allocated(M(i,i-1)%val)) then
        call destroy(M(i,i-1))
      end if
    end do

  end subroutine destroy_tridiag_blk


end module iterative_cpu

