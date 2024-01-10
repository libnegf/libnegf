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


module iterative

  use ln_precision
  use ln_constants, only : pi, i_unit => j, minusOne
  use ln_allocation
  use mat_def
  use sparsekit_drv
  use inversions
  use elph
  use ln_structure, only : TStruct_Info
  use lib_param, only : MAXNCONT, Tnegf, intarray, TInteractionList
  use mpi_globals, only : id, numprocs, id0
  use outmatrix, only : outmat_c, inmat_c, direct_out_c, direct_in_c
  use clock
  use ln_cache
  use ln_elastic
  use ln_inelastic
  !use transform

  implicit none
  private

  public :: calculate_transmissions
  public :: calculate_transmissions_and_dos

  public :: calculate_Gr
  public :: calculate_Gn_neq_components
  public :: calculate_elastic_scba
  public :: calculate_gsmr_blocks
  !public :: calculate_gsml_blocks
  public :: calculate_Gr_tridiag_blocks
  !public :: calculate_Gr_column_blocks
  public :: calculate_Gn_tridiag_blocks
  public :: calculate_Gr_outer
  public :: calculate_Gn_outer

  public :: iterative_meir_wingreen
  public :: iterative_layer_current
  public :: transmission_BP_corrected

  public :: destroy_all_blk

  logical, parameter :: debug=.false.
  ! These are here temporarily ...
  type(z_DNS), dimension(:), allocatable :: gsmr
  !type(z_DNS), dimension(:), allocatable :: gsml
  type(z_DNS), dimension(:,:), allocatable :: Gr
  type(z_DNS), dimension(:,:), allocatable :: ESH
  type(z_DNS), dimension(:,:), allocatable :: Gn


CONTAINS

  !****************************************************************************
  !
  ! Driver for computing Equilibrium Green Retarded (Outer blocks included)
  ! writing on memory
  !
  !****************************************************************************

  subroutine calculate_Gr(negf,E,SelfEneR,Tlc,Tcl,gsurfR,Grout,outer)

    !****************************************************************************
    !
    !Input
    !negf:    negf data container
    !E:        Energy point
    !SelfEneR: matrices array containing contacts Self Energy
    !Tlc:      matrices array containing contacts-device interaction blocks (ES-H)
    !Tcl:      matrices array containing device-contacts interaction blocks (ES-H)
    !gsurfR:   matrices array containing contacts surface green
    !outer:    optional parameter (0,1,2).
    !
    !Output:
    !Grout: Retarded Green's function (Device + Contacts overlap regions -> effective conductor)
    !   outer = 0  no outer parts are computed
    !   outer = 1  only D/C part is computed
    !   outer = 2  D/C and C/D parts are computed
    !              (needed for K-points calculations)
    !
    !*****************************************************************************

    type(Tnegf), intent(inout) :: negf
    complex(dp), intent(in) :: E
    type(z_DNS), dimension(:), intent(in) :: SelfEneR
    type(z_DNS), dimension(:), intent(inout) :: Tlc, Tcl, gsurfR
    type(z_CSR), intent(out) :: Grout
    integer, intent(in) :: outer

    !Work
    type(z_CSR) :: ESH_tot, Ain
    integer :: i,ierr, nbl, ncont,ii,n

    nbl = negf%str%num_PLs
    ncont = negf%str%num_conts

    ! Take CSR H,S and build ES-H in dense blocks
    call prealloc_sum(negf%H,negf%S,minusOne,E,ESH_tot)

    call allocate_blk_dns(ESH,nbl)

    call zcsr2blk_sod(ESH_tot, ESH, negf%str%mat_PL_start)

    call destroy(ESH_tot)

    associate(cblk=>negf%str%cblk)
    do i=1,ncont
      ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
    end do
    end associate

    !! Add interaction self energy contribution, if any
    call add_sigma_r(negf, ESH)

    call allocate_gsm(gsmr,nbl)
    call calculate_gsmr_blocks(ESH,nbl,2)

    call allocate_blk_dns(Gr,nbl)

    call calculate_Gr_tridiag_blocks(ESH,1)
    call calculate_Gr_tridiag_blocks(ESH,2,nbl)

    call destroy_tridiag_blk(ESH)
    call deallocate_blk_dns(ESH)

    call destroy_gsm(gsmr)
    call deallocate_gsm(gsmr)

    !! Deliver Gr to interaction models if any
    call set_Gr(negf, Gr)

    call blk2csr(Gr,negf%str,negf%S,Grout)

    SELECT CASE (outer)
    CASE(0)
    CASE(1)
      call calculate_Gr_outer(Tlc,Tcl,gsurfR,negf%str,.FALSE.,Grout)
    CASE(2)
      call calculate_Gr_outer(Tlc,Tcl,gsurfR,negf%str,.TRUE.,Grout)
    end SELECT

    call destroy_blk(Gr)
    DEALLOCATE(Gr)

  end subroutine calculate_Gr

  !****************************************************************************
  !
  ! Driver for computing G_n contributions due to all contacts MINUS reference:
  ! Reference is necessary when splitting into contour + real-axis integration
  !
  !   Sum   [f_j(E)-f_r(E)] Gr Gam_j Ga
  !   j!=r
  !
  ! NOTE: Setting reference to ncont+1 removes the reference part
  !****************************************************************************

  subroutine calculate_Gn_neq_components(negf,E,SelfEneR,Tlc,Tcl,gsurfR,frm,Glout,outblocks,Grout)

    !****************************************************************************
    !
    !Input
    !negf:    negf data container
    !E:        Energy point
    !SelfEneR: matrices array containing contacts Self Energy
    !Tlc:      matrices array containing contacts-device interaction blocks (ES-H)
    !Tcl:      matrices array containing device-contacts interaction blocks (ES-H)
    !gsurfR:   matrices array containing contacts surface green
    !outblocks: optional parameter (0,1,2).
    !
    !Output:
    !Gn: NE GF (Device + Contacts overlap regions -> effective conductor)
    !   outblocks = 0  no outer parts are computed
    !   outblocks = 1  only D/C part is computed
    !   outblocks = 2  D/C and C/D parts are computed
    !              (needed for K-points calculations)
    !
    !Optional: Grout (used with Gn to calculate Gp)
    !   outblokcs -> same options as above
    !*****************************************************************************

    implicit none

    !In/Out
    type(Tnegf), intent(inout) :: negf
    type(z_DNS), dimension(:), intent(in)  :: SelfEneR, gsurfR, Tlc, Tcl
    real(dp), intent(in)  :: E
    real(dp), dimension(:), intent(in)  :: frm
    type(z_CSR), intent(inout), optional  :: Glout, Grout
    integer, intent(in), optional  :: outblocks

    !Work
    integer :: ref
    complex(dp) :: Ec
    integer :: i,ierr,ncont,nbl, lbl, rbl
    integer, dimension(:), allocatable :: Gr_columns
    type(z_CSR) :: ESH_tot, Gl
    logical :: mask(MAXNCONT)


    nbl = negf%str%num_PLs
    ncont = negf%str%num_conts
    ref = negf%refcont
    associate (cblk=>negf%str%cblk, indblk=>negf%str%mat_PL_start)

    Ec = cmplx(E,0.0_dp,dp)

    ! Take CSR H,S and build ES-H in dense blocks
    call prealloc_sum(negf%H,negf%S,minusOne,Ec,ESH_tot)

    call allocate_blk_dns(ESH,nbl)

    call zcsr2blk_sod(ESH_tot,ESH, negf%str%mat_PL_start)

    call destroy(ESH_tot)

    do i=1,ncont
      ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val - SelfEneR(i)%val
    end do

    ! Add interaction self energies (if any)
    call add_sigma_r(negf, ESH)

    call allocate_gsm(gsmr,nbl)
    call calculate_gsmr_blocks(ESH,nbl,2)

    call allocate_blk_dns(Gr,nbl)

    ! compute Gr(1,1)
    call calculate_Gr_tridiag_blocks(ESH,1)
    ! compute Gr(n,n), Gr(n-1,n), Gr(n, n-1);  n = 2 .. nbl
    call calculate_Gr_tridiag_blocks(ESH,2,nbl)

    !Passing Gr to interactions
    call set_Gr(negf, Gr)

    !Computing device G_n
    call allocate_blk_dns(Gn,nbl)
    call init_tridiag_blk(Gn,ESH)

    call calculate_Gn_tridiag_blocks(negf,ESH,SelfEneR,frm,ref,negf%str,Gn)

    call destroy_gsm(gsmr)
    call deallocate_gsm(gsmr)

    !Passing G^n to interactions
    call set_Gn(negf, Gn)

    if (present(Glout)) then
      call blk2csr(Gn,negf%str,negf%S,Glout)
    end if

    if (present(Grout)) then
      call blk2csr(Gr,negf%str,negf%S,Grout)
    end if
    end associate



    !Computing the 'outer' blocks (device/contact overlapping elements)
    if (present(Glout)) then
      SELECT CASE (outblocks)
      CASE(0)
      CASE(1)
        call calculate_Gn_outer(Tlc,gsurfR,negf%str,frm,ref,.false.,Glout)
      CASE(2)
        call calculate_Gn_outer(Tlc,gsurfR,negf%str,frm,ref,.true.,Glout)
      end SELECT
    end if

    if (present(Grout)) then
      SELECT CASE (outblocks)
      CASE(0)
      CASE(1)
        call calculate_Gr_outer(Tlc,Tcl,gsurfR,negf%str,.false.,Grout)
      CASE(2)
        call calculate_Gr_outer(Tlc,Tcl,gsurfR,negf%str,.true.,Grout)
      end SELECT
    end if

    call destroy_all_blk(negf)

  end subroutine calculate_Gn_neq_components


  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  !---------------------------------------------------------------------
  subroutine calculate_elastic_scba(negf,E,SelfEneR,Tlc,Tcl,gsurfR,frm,scba_niter, &
            &  scba_tol, scba_error)
    type(Tnegf), intent(inout) :: negf
    type(z_DNS), dimension(:), intent(in)  :: SelfEneR, gsurfR, Tlc, Tcl
    real(dp), intent(in)  :: E
    real(dp), dimension(:), intent(in)  :: frm
    integer, intent(in) :: scba_niter
    real(dp), intent(in) :: scba_tol
    real(dp), intent(inout) :: scba_error

    logical :: tDestroyGn, tDestroyESH, tDestroyGr
    Type(z_CSR) :: csrGn
    integer :: scba_iter, outer = 0

    tDestroyGn = negf%tDestroyGn
    tDestroyGr = negf%tDestroyGr
    tDestroyESH = negf%tDestroyESH

    call negf%scbaDriverElastic%init(tol = scba_tol, dowrite = .false.)
    scba_iter = 0

    do while (.not.negf%scbaDriverElastic%is_converged() .and. scba_iter <= scba_niter)
      call negf%scbaDriverElastic%set_scba_iter(scba_iter, negf%interactList)

      negf%tDestroyGn = .true.
      negf%tDestroyGr = .true.
      negf%tDestroyESH = .true.

      call destroy_all_blk(negf)

      if (.not.tDestroyGn) negf%tDestroyGn = .false.
      if (.not.tDestroyGr) negf%tDestroyGr = .false.
      if (.not.tDestroyESH) negf%tDestroyESH = .false.

      call calculate_Gn_neq_components(negf,E,SelfEneR,Tlc,Tcl,gsurfR,frm,csrGn,outer)

      call negf%scbaDriverElastic%check_Mat_convergence(csrGn)
      call destroy(csrGn)

      scba_iter = scba_iter + 1

    enddo

    scba_error = negf%scbaDriverElastic%scba_error()

    call negf%scbaDriverElastic%destroy()

  end subroutine calculate_elastic_scba

  !---------------------------------------------------------------------
  !>
  !  Iterative algorithm implementing Meir Wingreen formula for a given
  !  electrode
  !  Note: self consistent born approximation is not accounted for here
  !  It is assumed that the el-ph container already includes the
  !  desired values. SCBA loop should be run outside
  !  Many operations from calls_neq_ph are repeated here, as it is
  !  assumed that A and Gn are not available at the time of the call
  !
  !
  !  I_i = Tr[\Sigma_{i}^{n} A - \Gamma_{i} G^{n}] =
  !      = Tr[\Gamma_{i}( f_{i} A - G^{n} )]
  !
  ! The subroutine assumes that Gr and Gn are available
  ! set negf%tDestroyGr,Gn = .false. when computing Gn
  !---------------------------------------------------------------------
  subroutine iterative_meir_wingreen(negf,E,SelfEneR,frm,curr_mat)
    type(Tnegf), intent(inout) :: negf
    real(dp) :: E
    type(z_DNS), dimension(:) :: SelfEneR
    real(dp), dimension(:) :: frm
    real(dp), dimension(:) :: curr_mat

    !Work
    complex(dp) :: Ec, tmp
    integer :: ii, ncont, lead, lead_blk
    type(z_DNS) :: work1, Gam, A

    ncont = negf%str%num_conts

    do ii = 1, ncont
      lead = ii
      lead_blk = negf%str%cblk(lead)
      call zspectral(SelfEneR(lead),SelfEneR(lead), 0, Gam)
      call prealloc_mult(Gam, Gn(lead_blk, lead_blk), work1)

      call zspectral(Gr(lead_blk, lead_blk), Gr(lead_blk, lead_blk), 0, A)
      tmp = cmplx(-frm(lead),0.0_dp, dp)
      call prealloc_mult(Gam, A, tmp, work1)
      curr_mat(ii) = -real(trace(work1))

      call destroy(work1)
      call destroy(Gam)
      call destroy(A)
    end do

    call destroy_all_blk(negf)

  end subroutine iterative_meir_wingreen

  !---------------------------------------------------------------------
  !  Calculate the layer current per unit energy
  !
  !    I_LL'(E) = Tr[(ES-H)_LL' * Gn_L'L(E)-(ES-H)_L'L * Gn_LL'(E)]
  !
  ! The subroutine assumes that Gn(ii, ii+1) blocks are available
  ! set negf%tDestroyGn = .false. when computing Gn
  !
  ! Note : ESH(ii,ii+1) are not changed by diagonal Sigma^r
  !        With non-diagonal sigma^r ESH must be recomputed
  !
  subroutine iterative_layer_current(negf,E,curr_mat,ldos_mat)
    type(Tnegf), intent(inout) :: negf
    real(dp), intent(in) :: E
    real(dp), dimension(:), intent(inout) :: curr_mat
    real(dp), dimension(:), intent(inout) :: ldos_mat

    integer :: nbl, ii
    type(z_DNS) :: work1
    type(z_CSR) :: ESH_tot
    complex(dp) :: Ec

    nbl = negf%str%num_PLs

    if (size(curr_mat) .ne. nbl-1) then
       stop 'ERROR: curr_mat with wrong size in iterative_layer_current'
    end if

    Ec=cmplx(E,0.0_dp,dp)
    call prealloc_sum(negf%H, negf%S, minusOne, Ec, ESH_tot)
    call allocate_blk_dns(ESH,nbl)
    call zcsr2blk_sod(ESH_tot, ESH, negf%str%mat_PL_start)
    call destroy(ESH_tot)

    do ii = 1, nbl-1
      call prealloc_mult(ESH(ii,ii+1),Gn(ii+1,ii),work1)
      call prealloc_mult(Gn(ii,ii+1),ESH(ii+1,ii),minusOne, work1)
      curr_mat(ii) = real(i_unit*trace(work1))
      call destroy(work1)
      ldos_mat(ii) = trace(Gn(ii,ii))
    end do
    ldos_mat(nbl) = trace(Gn(nbl,nbl))

    call destroy_all_blk(negf)

  end subroutine iterative_layer_current


  !--------------------------------------------------------------------------
  ! Add all Self-enrgies to the Hamiltonian
  subroutine add_sigma_r(negf, ESH)
    class(TNegf) :: negf
    type(z_DNS) :: ESH(:,:)

    type(TInteractionNode), pointer :: it
    it => negf%interactList%first
    do while (associated(it))
      call it%inter%add_sigma_r(ESH, negf%iE, negf%iKpoint, negf%spin)
      it => it%next
    end do

  end subroutine add_sigma_r

  !--------------------------------------------------------------------------
  ! Add all Self-enrgies to Sigma_n
  subroutine add_sigma_n(negf, sigma_n)
    class(TNegf) :: negf
    type(z_DNS) :: sigma_n(:,:)

    type(TInteractionNode), pointer :: it
    it => negf%interactList%first

    do while (associated(it))
      call it%inter%add_sigma_n(sigma_n, negf%iE, negf%iKpoint, negf%spin)
      it => it%next
    end do
  end subroutine add_sigma_n


  !--------------------------------------------------------------------------
  !> Store Gr
  !>
  subroutine cache_Gr(negf, Gr, en_index, k_index, spin, tridiagonal)
    class(TNegf) :: negf
    type(z_dns), dimension(:,:), intent(in) :: Gr
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin
    logical, intent(in), optional :: tridiagonal

    type(TMatLabel) :: label
    integer :: ii, jj, nbl
    logical :: tTridiag

    if (.not.associated(negf%G_r)) then
       allocate(TMatrixCacheMem::negf%G_r)
       select type(p => negf%G_r)
       type is(TMatrixCacheMem)
          p%tagname='G_r'
       end select
    end if

    label%kpoint = 0
    label%energy_point = 0
    label%spin = 0
    tTridiag = .true.

    if (present(k_index)) then
      label%kpoint = k_index
    end if
    if (present(k_index)) then
      label%energy_point = en_index
    end if
    if (present(spin)) then
      label%spin = spin
    end if
    if (present(tridiagonal)) then
      tTridiag = tridiagonal
    end if

    ! store tri- diagonal blocks
    do ii = 1, size(Gr,1)
      label%row_block = ii
      label%col_block = ii
      call negf%G_r%add(Gr(ii,ii), label)
      if (tTridiag .and. ii < size(Gr,1)) then
         label%col_block = ii + 1
         call negf%G_r%add(Gr(ii,ii+1), label)
         label%row_block = ii + 1
         label%col_block = ii
         call negf%G_r%add(Gr(ii+1,ii), label)
      end if
    end do

  end subroutine cache_Gr

  !--------------------------------------------------------------------------
  !> store Gn
  !>
  subroutine cache_Gn(negf, Gn, en_index, k_index, spin, tridiagonal)
    class(TNegf) :: negf
    type(z_dns), dimension(:,:), intent(in) :: Gn
    integer, intent(in), optional :: en_index
    integer, intent(in), optional :: k_index
    integer, intent(in), optional :: spin
    logical, intent(in), optional :: tridiagonal

    type(TMatLabel) :: label
    integer :: ii, jj, nbl
    logical :: tTridiag

    if (.not.associated(negf%G_n)) then
       allocate(TMatrixCacheMem::negf%G_n)
       select type(p => negf%G_n)
       type is(TMatrixCacheMem)
          p%tagname='G_n'
       end select
    end if

    label%kpoint = 0
    label%energy_point = 0
    label%spin = 0
    tTridiag = .true.

    if (present(k_index)) then
      label%kpoint = k_index
    end if
    if (present(k_index)) then
      label%energy_point = en_index
    end if
    if (present(spin)) then
      label%spin = spin
    end if
    if (present(tridiagonal)) then
      tTridiag = tridiagonal
    end if

    ! Store diagonal and upper diagonal
    do ii = 1, size(Gn,1)
      label%row_block = ii
      label%col_block = ii
      call negf%G_n%add(Gn(ii,ii), label)
      if (tTridiag .and. ii < size(Gn,1)) then
         label%col_block = ii + 1
         call negf%G_n%add(Gn(ii,ii+1), label)
         label%row_block = ii + 1
         label%col_block = ii
         call negf%G_n%add(Gn(ii+1,ii), label)
      end if
    end do

  end subroutine cache_Gn

  ! Provides Gr to the interaction models.
  ! In some case the self/energies are computed
  subroutine set_Gr(negf, Gr)
    type(TNegf) :: negf
    type(z_DNS), intent(in) :: Gr(:,:)

    integer :: iE, iK, iSpin
    type(TInteractionNode), pointer :: it
    it => negf%interactList%first
    iE=negf%iE
    iK=negf%iKpoint
    iSpin=negf%spin

    do while (associated(it))
      select type(pInter => it%inter)
      class is (Telastic)
        call it%inter%set_Gr(Gr, iE, iK, iSpin)
      class is (TInelastic)
        ! cache Gr in negf container and pass the pointer
        call cache_Gr(negf, Gr, iE, iK, iSpin, pInter%tTridiagonal)
        call pInter%set_Gr_pointer(negf%G_r)
      end select
      it => it%next
    end do

  end subroutine set_Gr

  ! Provides Gn to the interaction models.
  ! In some case the self/energies are computed
  subroutine set_Gn(negf, Gn)
    type(TNegf) :: negf
    type(z_DNS), intent(in) :: Gn(:,:)

    integer :: iE, iK, iSpin
    type(TInteractionNode), pointer :: it
    it => negf%interactList%first
    iE=negf%iE
    iK=negf%iKpoint
    iSpin=negf%spin

    do while (associated(it))
      select type(pInter => it%inter)
      class is (Telastic)
        call it%inter%set_Gn(Gn, iE, iK, iSpin)
      class is (TInelastic)
        ! cache Gn in negf container and pass the pointer
        call cache_Gn(negf, Gn, iE, iK, iSpin, pInter%tTridiagonal)
        call pInter%set_Gn_pointer(negf%G_n)
      end select
      it => it%next
    end do

  end subroutine set_Gn

  !------------------------------------------------------------------------------!
  ! Transmission_BP_corrected
  !
  ! The routine implements the current-probes, i.e. it assumes that the BPs are
  ! elastic dephasing sources that fullfill current conservation at any energy
  ! For a nice discussion see papers by D. Ryndyk and
  ! M. Kilgour, D. Segal, The J. of Chem Phys 144, 124107 (2016)
  !
  ! LIMITATIONS:
  !
  ! This routine is limited to 2 contacts
  ! It uses dense matrices and assumes there is only 1 PL in the device
  !------------------------------------------------------------------------------!

  subroutine transmission_BP_corrected(negf,SelfEneR,tun_mat)

    type(Tnegf), intent(inout) :: negf
    type(z_DNS), dimension(:) :: SelfEneR
    real(dp), dimension(:) :: tun_mat

    integer :: nn, mm, cont, ni, nf
    integer :: NumOrbs, ncont
    integer, allocatable ::  pls(:), ple(:)
    real(dp), allocatable  :: Trans(:,:), W(:,:), GG(:,:), R(:)
    type(z_DNS) :: Gam1, Gam2, GreenR, GreenA, Tmp1, Tmp2


    NumOrbs=negf%str%central_dim
    ncont = negf%str%num_conts

    allocate(Trans(NumOrbs+ncont,NumOrbs+ncont))
    Trans=0.0_dp

    call create(Gam1, NumOrbs, NumOrbs)
    call create(Gam2, NumOrbs, NumOrbs)

    call blk2dns(Gr,negf%str,GreenR)
    call zdagger(GreenR, GreenA)

    allocate(pls(ncont))
    allocate(ple(ncont))

    do nn = 1, ncont
      pls(nn) = negf%str%mat_PL_start(negf%str%cblk(nn))
      ple(nn) = negf%str%mat_PL_end(negf%str%cblk(nn))
    end do

    do nn = 1, NumOrbs+ncont

      if (nn > NumOrbs) then
         cont = nn - NumOrbs
         call zspectral(SelfEneR(cont),SelfEneR(cont),0,Tmp1)
         Gam1%val(pls(cont):ple(cont),pls(cont):ple(cont)) = Tmp1%val
         call destroy(Tmp1)
      else
         Gam1%val = 0.0_dp
         Gam1%val(nn,nn)=negf%bp_deph%coupling(nn)
      end if

      do mm = 1, NumOrbs+ncont

        if (mm > NumOrbs) then
           cont = mm - NumOrbs
           call zspectral(SelfEneR(cont),SelfEneR(cont),0,Tmp2)
           Gam2%val(pls(cont):ple(cont),pls(cont):ple(cont)) = Tmp2%val
           call destroy(Tmp2)
        else
           Gam2%val = 0.0_dp
           Gam2%val(mm,mm)=negf%bp_deph%coupling(mm)
        end if

        ! Compute coherent transmission: Tr[Gam1 Gr Gam2 Ga]
        ! The inefficient quadruple loop has been substituted with
        ! M*M multiplications exploting OMP parallelism
        call prealloc_mult(Gam1,GreenR,Tmp1)
        call prealloc_mult(Tmp1,Gam2,Tmp2)
        call destroy(Tmp1)
        call prealloc_mult(Tmp2,GreenA,Tmp1)
        call destroy(Tmp2)

        Trans(nn,mm) = real(trace(Tmp1))

        call destroy(Tmp1)

      end do
    end do

    call destroy(Gam1, Gam2, GreenR, GreenA)
    deallocate(pls, ple)

    allocate(W(NumOrbs,NumOrbs))
    allocate(R(NumOrbs+ncont))

    do nn = 1, NumOrbs+ncont
       R(nn) = 1.0_dp
       do mm= 1, NumOrbs+ncont
          if (mm.ne.nn) then
             R(nn) = R(nn)-Trans(mm,nn)
          end if
       end do
    end do

    do nn = 1, NumOrbs
      do mm = 1, NumOrbs
         W(nn,mm) = -Trans(nn,mm)
      end do
      W(nn,nn) = 1.0_dp - R(nn)
    end do

    deallocate(R)
    allocate(GG(NumOrbs,NumOrbs))

    call inverse(GG,W,NumOrbs)

    deallocate(W)

    allocate(R(NumOrbs))

    do nn = 1, size(negf%ni)
      ni = negf%ni(nn)
      nf = negf%nf(nn)
      R =  matmul(Trans(NumOrbs+ni,:),GG)
      tun_mat(nn) = Trans(NumOrbs+ni,NumOrbs+nf) + dot_product(R, Trans(:,NumOrbs+nf))
    end do

    deallocate(Trans,GG)
    deallocate(R)

  end subroutine transmission_BP_corrected

  !------------------------------------------------------------------------------!
  !DAR end
  !------------------------------------------------------------------------------!

  !**********************************************************************
  subroutine init_tridiag_blk(Matrix,S)
    type(z_DNS), dimension(:,:) :: Matrix,S

    integer :: nbl, j

    nbl = SIZE(Matrix,1)

    call create(Matrix(1,1),S(1,1)%nrow,S(1,1)%ncol)
    Matrix(1,1)%val=(0.0_dp,0.0_dp)
    if (nbl.gt.1) then
       do j=2,nbl
         call create(Matrix(j-1,j),S(j-1,j)%nrow,S(j-1,j)%ncol)
         Matrix(j-1,j)%val=(0.0_dp,0.0_dp)
         call create(Matrix(j,j),S(j,j)%nrow,S(j,j)%ncol)
         Matrix(j,j)%val=(0.0_dp,0.0_dp)
         call create(Matrix(j,j-1),S(j,j-1)%nrow,S(j,j-1)%ncol)
         Matrix(j,j-1)%val=(0.0_dp,0.0_dp)
       end do
    endif

  end subroutine init_tridiag_blk

  !---------------------------------------------------------------------
  subroutine destroy_all_blk(negf)
    type(Tnegf), intent(in) :: negf

    if (negf%tDestroyESH) then
      call destroy_tridiag_blk(ESH)
      if (allocated(ESH)) deallocate(ESH)
    end if
    if (negf%tDestroyGr) then
      call destroy_blk(Gr)
      if (allocated(Gr)) deallocate(Gr)
    end if
    if (negf%tDestroyGn) then
       call destroy_blk(Gn)
       if (allocated(Gn)) deallocate(Gn)
    end if
  end subroutine destroy_all_blk

  !**********************************************************************
  subroutine destroy_gsm(gsm)
    type(z_DNS), dimension(:) :: gsm
    integer :: i, i1, nbl

    nbl=size(gsm,1)

    do i=1,nbl
      if (allocated(gsm(i)%val)) call destroy(gsm(i))
    end do

  end subroutine destroy_gsm

  !**********************************************************************
  subroutine destroy_blk(M)
    type(z_DNS), dimension(:,:), allocatable :: M
    integer :: i, i1, nbl

    if (.not.allocated(M)) return

    nbl=size(M,1)

    do i=1,nbl
      do i1=1,nbl
        if (ALLOCATED(M(i1,i)%val)) THEN
          call destroy(M(i1,i))
        end if
      end do
    end do

  end subroutine destroy_blk

  !**********************************************************************
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

  !***********************************************************************
  !
  !  g_small right (gsmr) calculation - write on memory
  !
  !***********************************************************************

  subroutine calculate_gsmr_blocks(ESH,sbl,ebl,keepall)

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
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl,ebl                 ! start block, end block
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
  !  g_small left (gsml) calculation - write on memory
  !
  !***********************************************************************

 ! subroutine calculate_gsml_blocks(ESH,sbl,ebl)

 !   !***********************************************************************
 !   !Input:
 !   !ESH: sparse matrices array ESH(nbl,nbl)
 !   !
 !   !nbl (number of layers), indblk(nbl+1) global variables needed
 !   !
 !   !Output:
 !   !sparse matrices array global variable gsml(nbl) is available in memory
 !   !single blocks are allocated internally, array Gr(nbl,nbl)
 !   !must be allocated externally
 !   !***********************************************************************


 !   implicit none

 !   !In/Out
 !   type(z_DNS), dimension(:,:), intent(in) :: ESH
 !   integer, intent(in) :: sbl,ebl                       ! start block, end block

 !   !Work
 !   type(z_DNS) :: work1, work2
 !   integer :: nrow
 !   integer :: i, nbl
 !   !type(z_DNS) :: INV(sbl,sbl)

 !   if (sbl.gt.ebl) return

 !   nbl = size(ESH,1)

 !   if (nbl.eq.1) return

 !   nrow=ESH(sbl,sbl)%nrow

 !   call create(gsml(sbl),nrow,nrow)

 !   call compGreen(gsml(sbl),ESH(sbl,sbl),nrow)


 !   do i=sbl+1,ebl

 !     nrow=ESH(i,i)%nrow

 !     call prealloc_mult(ESH(i,i-1),gsml(i-1),(-1.0_dp, 0.0_dp),work1)

 !     call prealloc_mult(work1,ESH(i-1,i),work2)

 !     call destroy(work1)

 !     call prealloc_sum(ESH(i,i),work2,work1)

 !     call destroy(work2)

 !     call create(gsml(i),work1%nrow,work1%nrow)

 !     call compGreen(gsml(i),work1,work1%nrow)

 !     call destroy(work1)

 !   end do

 !   if (debug) then
 !     WRITE(*,*) '********************'
 !     WRITE(*,*) 'calculate_gsml done'
 !     WRITE(*,*) '********************'
 !   endif

 ! end subroutine calculate_gsml_blocks





  !***********************************************************************
  !
  !  Diagonal, Subdiagonal, Superdiagonal blocks of Green Retarded
  !  Gr(nbl,nbl) - writing on memory
  !
  !***********************************************************************

  subroutine calculate_Gr_tridiag_blocks(ESH,sbl,ebl)

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
    type(z_DNS), dimension(:,:) :: ESH
    integer :: sbl
    integer, optional :: ebl

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
          stop "Error: Gr_tridiag requires gsml"
          !call prealloc_mult(ESH(sbl,sbl-1),gsml(sbl-1),work2)
          !call prealloc_mult(work2,ESH(sbl-1,sbl),work3)
          !call destroy(work2)
          !call prealloc_sum(work1,work3,(-1.0_dp, 0.0_dp),work2)
          !call destroy(work3)
          !work1%val = work2%val
          !call destroy(work2)
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
      do i=sbl,ebl,-1
        stop "Error: Gr_tridiag requires gsml"
        !call prealloc_mult(gsml(i),ESH(i,i+1),work1)
        !call prealloc_mult(work1,Gr(i+1,i+1),(-1.0_dp,0.0_dp),Gr(i,i+1))
        !call destroy(work1)

        !call prealloc_mult(ESH(i+1,i),gsml(i),work2)
        !call prealloc_mult(Gr(i+1,i+1),work2,(-1.0_dp, 0.0_dp),Gr(i+1,i))

        !call prealloc_mult(Gr(i,i+1),work2,(-1.0_dp,0.0_dp),work1)
        !call destroy(work2)

        !call prealloc_sum(gsml(i),work1,Gr(i,i))
        !call destroy(work1)
      end do
    endif

  end subroutine calculate_Gr_tridiag_blocks

  !**************************************************************************
  !
  !  Calculate Green Retarded column "n" - writing on memory
  !
  !**************************************************************************
 ! subroutine calculate_Gr_column_blocks(ESH,n,indblk)

 !   !***********************************************************************
 !   !Input:
 !   !ESH: sparse matrices array ESH(nbl,nbl)
 !   !n: n umber of column to be calculated
 !   !
 !   !global variables needed: nbl (number of layers), indblk(nbl+1),
 !   !Gr diagonal, subadiagonal and superdiagonal, gsmr(:) for
 !   !downgoing and gsml(:) for upgoing
 !   !
 !   !Output:
 !   !sparse matrices array global variable Gr(:,n) is available in
 !   !memory - single blocks are allocated internally, array Gr(nbl,nbl)
 !   !must be allocated externally
 !   !***********************************************************************

 !   implicit none

 !   !In/Out
 !   type(z_DNS), dimension(:,:), intent(in) :: ESH
 !   integer, intent(in) :: n
 !   integer, dimension(:), intent(in) :: indblk

 !   !Work
 !   integer :: i,nrow,ncol,nbl
 !   type(z_DNS) :: work1
 !   real(dp) :: max

 !   nbl = size(ESH,1)

 !   if (n.GT.nbl) THEN
 !     STOP 'Error in calculate_Grcol : n is greater than nbl'
 !   endif

 !   !***************************************
 !   !  Downgoing (j>=n+2 && n<nbl-1)
 !   !
 !   !   G_j,n = -gR_jj T_j,j-1 G_j-1,n
 !   !
 !   !***************************************
 !   if (n.LT.(nbl-1)) THEN

 !     do i=n+2,nbl

 !       max=MAXVAL(ABS(Gr(i-1,n)%val))
 !       if (max.GT.EPS) THEN
 !         call prealloc_mult(gsmr(i),ESH(i,i-1),(-1.0_dp, 0.0_dp),work1)
 !         call prealloc_mult(work1,Gr(i-1,n),Gr(i,n))
 !         call destroy(work1)
 !       else
 !         ! WHEN BLOCK IS SMALLER THAN EPS IT IS NOT CREATED
 !         exit
 !       end if

 !     end do

 !   endif
 !   !*************************************
 !   !   Up-going (j<=n-2 && n>2)
 !   !
 !   !   G_j,n = -gL_jj T_j,j+1 G_j+1,n
 !   !
 !   !*************************************

 !   if (n.GT.2) THEN

 !     do i=n-2,1,(-1)

 !       max=MAXVAL(ABS(Gr(i+1,n)%val))

 !       if (max.GT.EPS) THEN
 !         call prealloc_mult(gsml(i),ESH(i,i+1),(-1.0_dp, 0.0_dp),work1)
 !         call prealloc_mult(work1,Gr(i+1,n),Gr(i,n))
 !         call destroy(work1)
 !       else
 !         ! WHEN BLOCK IS SMALLER THAN EPS IT IS NOT CREATED
 !         exit
 !       endif

 !     end do

 !   endif

 !   if (debug) then
 !     WRITE(*,*) '******************************'
 !     WRITE(*,*) 'calculate_Grcol done column',n
 !     WRITE(*,*) '******************************'
 !   endif

 ! end subroutine calculate_Gr_column_blocks

  !****************************************************************************
  !
  !  Calculate Green Retarded - writing on memory (optimized on mask)
  !
  !****************************************************************************

  subroutine Gr_blk2csr(P,nbl,indblk,A)

    !****************************************************************************
    !Input:
    !P: CSR matrix containing masking pattern
    !
    !global variable needed: nbl, indblk(nbl+1), Gr(:,:)
    !
    !Output:
    !A: sparse matrix containing Green Retarded of device (allocated internally)
    !****************************************************************************

    implicit none

    !In/Out
    integer :: nbl
    integer, dimension(:), pointer :: indblk
    type(z_CSR) :: A, P, GrCsr

    !Work
    integer :: i, j, i1, ix, iy, x, y, col, oldx

    !create A with same pattern of P
    call create(A,P%nrow,P%ncol,P%nnz)
    A%rowpnt(:)=P%rowpnt(:)
    A%colind(:)=P%colind(:)
    A%nzval = (0.0_dp,0.0_dp)

    !If only one block is present, concatenation is not needed
    !and it's implemented in a more trivial way
    if (nbl.EQ.1) THEN

      call create(GrCsr,Gr(1,1)%nrow,Gr(1,1)%ncol,Gr(1,1)%nrow*Gr(1,1)%ncol)
      call dns2csr(Gr(1,1),GrCsr)

      call mask(GrCsr,P,A)
      call destroy(GrCsr)

    ELSE

      !Cycle upon all rows
      x = 1
      do i = 1, A%nrow
        !Choose which block (row) we're dealing with
        oldx = x

        !Check if row is in same block of previous or in next block. Not needed
        !(and not allowed not to exceed indblk index boundaries) if we're in the last block
        if (oldx.EQ.nbl) THEN
          x = oldx
        ELSE
          do ix = oldx, oldx+1
            if ( (i.GE.indblk(ix)).AND.(i.LT.indblk(ix+1)) ) x = ix
          end do
        endif

        !Offset: i1 is the index for separate blocks
        i1 = i - indblk(x) + 1
        !Cycle upon columns
        do j = A%rowpnt(i), A%rowpnt(i+1) -1
          !Choose which block column we're dealing with
          y = 0
          if (x.EQ.1) THEN
            if ( (A%colind(j).GE.indblk(x)).AND.(A%colind(j).LT.indblk(x + 1)) ) then
              y = 1
            ELSEif ( (A%colind(j).GE.indblk(x + 1)).AND.(A%colind(j).LT.indblk(x + 2)) ) then
              y = 2
            endif
          elseif (x.eq.nbl) then
            if ( (A%colind(j).GE.indblk(x)).AND.(A%colind(j).LT.indblk(x + 1)) ) then
              y = nbl
            ELSEif ( (A%colind(j).GE.indblk(x - 1)).AND.(A%colind(j).LT.indblk(x)) ) then
              y = nbl - 1
            endif
          ELSE
            do iy = x-1, x+1
              if ( (A%colind(j).GE.indblk(iy)).AND.(A%colind(j).LT.indblk(iy + 1)) ) y = iy
            end do
          endif
          if (y.eq.0) then
            write(*,*)
            write(*,*) 'ERROR in Gr_blk2csr: probably wrong PL size',x
            write(*,*) 'row',i,A%colind(j)
            write(*,*) 'block indeces:',indblk(1:nbl)
            stop
          endif

          col = A%colind(j) - indblk(y) + 1

          A%nzval(j) = Gr(x,y)%val(i1,col)

        end do

      end do

    endif

    !if (debug) call writePeakInfo(6)
    if (debug) then
      WRITE(*,*) '**********************'
      WRITE(*,*) 'calculate_GreenR done'
      WRITE(*,*) '**********************'
    endif

  end subroutine Gr_blk2csr


  ! Implements a new algorithm based on an iterative scheme to solve
  ! [ES - H - SigmaR] Gr = Sigma< Ga
  ! The subroutine uses the gsmr(:) computed before and makes an iteration
  ! upward to build gsmn and then downward to build the 3-diagonal blocks
  subroutine calculate_Gn_tridiag_blocks(negf,ESH,SelfEneR,frm,ref,struct,Gn)
    type(TNegf), intent(in) :: negf
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    type(z_DNS), dimension(:), intent(in) :: SelfEneR
    real(dp), dimension(:), intent(in) :: frm
    integer, intent(in) :: ref
    type(Tstruct_info), intent(in) :: struct
    type(z_DNS), dimension(:,:), intent(inout) :: Gn

    !Work
    type(z_DNS), dimension(:,:), allocatable :: Sigma_n
    type(z_DNS) :: work1, Ga, Gam
    complex(dp) :: frmdiff
    integer :: i, j
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
      call prealloc_mult(gsmr(nbl), Sigma_n(nbl,nbl), work)
      call prealloc_mult(work, gsmrDag, gns)
      call destroy(gsmrDag, work)

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
        call prealloc_mult(work, ESHdag, minusOne, Sigma_n(i,i))
        call destroy(work, gsmrDag, ESHdag)

        !work3 = ESH(i,i+1) gsmr(i+1) Sigma(i+1,i)
        call prealloc_mult(ESH(i,i+1), gsmr(i+1), work)
        call prealloc_mult(work, Sigma_n(i+1,i), minusOne, Sigma_n(i,i))
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

  ! blk-sparse to dense converision
  subroutine blk2dns(G,str,Gdns)
    type(z_DNS), dimension(:,:), intent(in) :: G
    type(Tstruct_info), intent(in) :: str
    type(z_DNS) :: Gdns

    integer :: n, m, ii, jj, pl_start1, pl_end1, pl_start2, pl_end2, nbl, nrows

    nbl = str%num_PLs
    nrows = str%central_dim
    call create(Gdns,nrows,nrows)

    do n = 1, nbl
      do m= 1, nbl
         pl_start1 = str%mat_PL_start(n)
         pl_end1 = str%mat_PL_end(n)
         pl_start2 = str%mat_PL_start(m)
         pl_end2 = str%mat_PL_end(m)

         do ii = 1, pl_end1 - pl_start1 + 1
           do jj = 1, pl_end2 - pl_start2 + 1
               Gdns%val(pl_start1 + ii - 1,pl_start2 + jj - 1) = G(n,m)%val(ii,jj)
           end do
         end do
      end do
    end do

  end subroutine blk2dns


  !Concatenation for every contact in G_n. Performs a sum on elements, not a replacement
  !Similar to calculate_GreenR2, except for sum of elements
  !Note: to backup old version zconcat calls (and Glsub deallocations) must be
  !      uncommented and all this part removed
  !If only one block is present, concatenation is not needed and it's implemented in a
  !more trivial way
  subroutine blk2csr(G,struct,P,Gcsr)

    type(z_DNS), dimension(:,:) :: G
    type(Tstruct_info), intent(in) :: struct
    type(z_CSR) :: Gcsr
    type(z_CSR) :: P, G_sp

    integer :: nbl, oldx, row, col, iy, ix, x, y, ii, jj, nrows

    nbl = struct%num_PLs
    nrows = struct%mat_PL_end(nbl)

    !create Gcsr with same pattern of P
    call create(Gcsr,P%nrow,P%ncol,P%nnz)
    Gcsr%rowpnt = P%rowpnt
    Gcsr%colind = P%colind
    Gcsr%nzval = (0.0_dp, 0.0_dp)

    associate(indblk=>struct%mat_PL_start)
    !Cycle upon all rows
    x = 1
    do ii = 1, nrows
      !Search block x containing row ii
      oldx = x
      if (oldx.EQ.nbl) THEN
        x = oldx
      ELSE
        do ix = oldx, oldx+1
          if ( (ii.GE.indblk(ix)).AND.(ii.LT.indblk(ix+1)) ) x = ix
        end do
      endif

      !Offset: row is the index for separate blocks
      row = ii - indblk(x) + 1

      !Cycle upon columns of Gcsr (which has been ALREADY MASKED by S)
      do jj = Gcsr%rowpnt(ii), Gcsr%rowpnt(ii+1) -1
        if (Gcsr%colind(jj).gt.nrows) CYCLE
        !Choose which block column we're dealing with
        y = 0
        if (x.eq.1) then
          if ( (Gcsr%colind(jj).GE.indblk(x)).AND.(Gcsr%colind(jj).LT.indblk(x + 1)) ) then
            y = 1
          ELSEif ( (Gcsr%colind(jj).GE.indblk(x + 1)).AND.(Gcsr%colind(jj).LT.indblk(x + 2)) ) then
            y = 2
          endif
        elseif (x.eq.nbl) then
          if ( (Gcsr%colind(jj).GE.indblk(x)).AND.(Gcsr%colind(jj).LT.indblk(x + 1)) ) then
            y = nbl
          ELSEif ( (Gcsr%colind(jj).GE.indblk(x - 1)).AND.(Gcsr%colind(jj).LT.indblk(x)) ) then
            y = nbl - 1
          endif
        else
          do iy = x-1, x+1
            if ( (Gcsr%colind(jj).GE.indblk(iy)).AND.(Gcsr%colind(jj).LT.indblk(iy + 1)) ) y = iy
          end do
        endif

        if (y.EQ.0) THEN
          write(*,*)
          write(*,*) 'ERROR in blk2csr: probably wrong PL size', x
          write(*,*) 'row',ii,nrows,Gcsr%colind(jj)
          write(*,*) 'block indeces:',indblk(1:nbl)
          STOP
        endif

        col = Gcsr%colind(jj) - indblk(y) + 1

        if (allocated(G(x,y)%val)) THEN
          Gcsr%nzval(jj) = Gcsr%nzval(jj) + G(x,y)%val(row,col)
        endif

      end do

    end do
    end associate

  end subroutine blk2csr


  !****************************************************************************
  !
  !  Calculate Green Retarded in the
  !  contacts regions, where overlap with device orbitals is non-zero
  !  (only the upper part, needed in both K=0 and K points calculations,
  !   or both upper and lower parts)
  !  writing on memory
  !
  !****************************************************************************

  subroutine calculate_Gr_outer(Tlc,Tcl,gsurfR,struct,lower,Aout)

    !****************************************************************************
    !Input:
    !Tlc: sparse arraymatrices array containing contact-device interacting blocks
    !     ESH-H
    !Tlc: sparse arraymatrices array containing device-contact interacting blocks
    !     ESH-H
    !gsurfR: sparse matrices array containing contacts surface green
    !lower: if .true., also lower parts are calculated and concatenated
    !
    !global variables needed: nbl, indblk(nbl+1), ncont, cblk(ncont),
    !cindblk(ncont), Gr in the diagonal block related to the interaction layer
    !
    !Output:
    !Aout: sparse matrix containing density matrix in the region
    !      corresponding to non-zero overlap
    !
    !****************************************************************************

    implicit none

    !In/Out
    type(z_DNS), dimension(:), intent(in) :: Tlc,Tcl,gsurfR
    logical, intent(in) :: lower
    type(Tstruct_info), intent(in) :: struct
    type(z_CSR), intent(inout) :: Aout


    !Work
    type(z_DNS) :: work1, Grcl, Grlc
    type(z_CSR) :: GrCSR, TCSR
    integer :: i,cb,nrow_tot,i1,j1
    integer :: ncont, nbl

    ncont = struct%num_conts
    nbl = struct%num_PLs
    nrow_tot = struct%total_dim

    if (.not.allocated(Aout%nzval)) THEN
      call create(Aout,nrow_tot,nrow_tot,0)
      Aout%rowpnt(:)=1
    endif

    do i=1,ncont

      !Numero di blocco del contatto
      cb=struct%cblk(i)
      call prealloc_mult(Gr(cb,cb),Tlc(i),minusOne,work1)
      call prealloc_mult(work1,gsurfR(i),Grlc)

      call destroy(work1)

      j1=nzdrop(Grlc,EPS)
      call create(GrCSR,Grlc%nrow,Grlc%ncol,j1)
      call dns2csr(Grlc,GrCSR)
      call destroy(Grlc)
      j1=nzdrop(Tlc(i),EPS)
      call create(TCSR,Tlc(i)%nrow,Tlc(i)%ncol,j1)
      call dns2csr(Tlc(i),TCSR)
      call zmask_realloc(GrCSR,TCSR)
      call destroy(TCSR)

      !Concatenazione di Asub nella posizione corrispondente
      i1=struct%mat_PL_start(cb)
      j1=struct%mat_B_start(i)

      call concat(Aout,GrCSR,i1,j1)

      call destroy(GrCSR)

      if (lower) THEN

        call prealloc_mult(gsurfR(i),Tcl(i),minusOne, work1)
        call prealloc_mult(work1, Gr(cb,cb), Grcl)

        call destroy(work1)

        j1=nzdrop(Grcl,EPS)
        call create(GrCSR,Grcl%nrow,Grcl%ncol,j1)
        call dns2csr(Grcl,GrCSR)
        call destroy(Grcl)
        j1=nzdrop(Tcl(i),EPS)
        call create(TCSR,Tcl(i)%nrow,Tcl(i)%ncol,j1)
        call dns2csr(Tcl(i),TCSR)
        call zmask_realloc(GrCSR,TCSR)
        call destroy(TCSR)

        i1 = struct%mat_B_start(i)-struct%central_dim+struct%mat_PL_start(nbl+1)-1
        j1 = struct%mat_PL_start(cb)

        call concat(Aout,GrCSR,i1,j1)

        call destroy(GrCSR)

      endif

    end do

    if (debug) then
      WRITE(*,*) '********************'
      WRITE(*,*) 'Outer_GreenR done'
      WRITE(*,*) '********************'
    endif

  end subroutine calculate_Gr_outer


  subroutine calculate_Gn_outer(Tlc,gsurfR,struct,frm,ref,lower,Gn_out)
    type(z_DNS), dimension(:), intent(in) :: Tlc, gsurfR
    real(dp), dimension(:), intent(in) :: frm
    type(Tstruct_info), intent(in) :: struct
    integer, intent(in) :: ref
    logical, intent(in) :: lower
    type(z_CSR) :: Gn_out

    !Work
    type(z_DNS) :: gsurfA, work1, work2, work3, Gn_lc
    type(z_CSR) :: GnCSR, TCSR
    integer :: k,cbk,i1,j1,nrow_tot
    integer :: ncont, nbl
    complex(dp) :: frmdiff

    ncont = struct%num_conts
    nbl = struct%num_PLs

    if (.not.allocated(Gn_out%nzval)) THEN
      nrow_tot = struct%total_dim
      call create(Gn_out,nrow_tot,nrow_tot,0)
      Gn_out%rowpnt(:)=1
    endif

    do k=1,ncont
      ! Sigma_n(cbk,ck) = Sigma_r(cbk,ck) = 0  => Tlc^r = Tlc^a = Tlc = T(cbk,ck)
      ! Gn(cbk,ck) =  -Gn(cbk,cbk)*Tlc*gA(cb) -(fk-fr)*Gr(cbk,cbk)*Tlc*j(gsurfR-gsurfA)]

      cbk=struct%cblk(k)
      call zdagger(gsurfR(k),gsurfA)
      call prealloc_mult(Tlc(k),gsurfA,work2)
      call destroy(gsurfA)

      call prealloc_mult(Gn(cbk,cbk),work2,minusOne,work3)
      call destroy(work2)

      !Checks that Fermi levels are sufficiently different and contact is not reference
      if ((ABS(frm(k)-frm(ref)).GT.EPS).AND.(k.NE.ref)) THEN

        frmdiff = cmplx(frm(ref)-frm(k),0.0_dp,dp)
        call zspectral(gsurfR(k),gsurfR(k),0,work1)

        call prealloc_mult(Tlc(k),work1,work2)
        call destroy(work1)

        call prealloc_mult(Gr(cbk,cbk),work2,frmdiff,work1)
        call destroy(work2)
      else
        call create(work1, Gr(cbk,cbk)%nrow, Tlc(k)%ncol)
        work1%val=(0.0_dp, 0.0_dp)
      end if

      call prealloc_sum(work3,work1,Gn_lc)
      call destroy(work1)
      call destroy(work3)

      call mask(Gn_lc,Tlc(k))
      i1=nzdrop(Gn_lc,EPS)

      if (i1.gt.0) THEN
        call create(GnCSR,Gn_lc%nrow,Gn_lc%ncol,i1)
        call dns2csr(Gn_lc,GnCSR)
        call destroy(Gn_lc)

        ! GnCSR is concatenated (added) to Gn_out
        i1=struct%mat_PL_start(cbk)
        j1=struct%mat_B_start(k)-struct%central_dim+struct%mat_PL_start(nbl+1)-1
        call concat(Gn_out,GnCSR,i1,j1)

        ! lower Gn(c,cb) outer part is computed via Gn(c,cb) = Gn(cb,c)+
        if (lower) THEN
          call zdagger(GnCSR,TCSR)
          call concat(Gn_out,TCSR,j1,i1)
          call destroy(TCSR)
        endif

        call destroy(GnCSR)
      else
        call destroy(Gn_lc)
      end if

    end do


  end subroutine calculate_Gn_outer

  !---------------------------------------------------

  subroutine calculate_transmissions(H,S,Ec,SelfEneR,ni,nf,str,tun_proj,tun_mat)
    Type(z_CSR) :: H
    Type(z_CSR) :: S
    Complex(dp) :: Ec
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR
    Integer :: ni(:)
    Integer :: nf(:)
    Type(z_CSR) :: ESH_tot
    Type(TStruct_Info) :: str
    type(intarray), intent(in) :: tun_proj
    Real(dp), Dimension(:) :: tun_mat

    ! Local variables
    Real(dp) :: tun
    Integer :: nbl,ncont,ibl
    Integer :: i, ierr, icpl, nit, nft, nt, nt1

    nbl = str%num_PLs
    ncont = str%num_conts

    if (ncont == 1) then
      tun_mat = 0.0_dp
      return
    end if

    !Calculation of ES-H and brak into blocks
    call prealloc_sum(H,S,minusOne,Ec,ESH_tot)

    call allocate_blk_dns(ESH,nbl)

    call zcsr2blk_sod(ESH_tot,ESH,str%mat_PL_start)
    call destroy(ESH_tot)

    !Inclusion of the contact Self-Energies to the relevant blocks
    do i=1,ncont
      ibl = str%cblk(i)
      ESH(ibl,ibl)%val = ESH(ibl,ibl)%val-SelfEneR(i)%val
    end do

    nit=ni(1)
    nft=nf(1)
    ! find the contact with smaller block index
    if (str%cblk(nit).lt.str%cblk(nft)) then
      nt = str%cblk(nit)
    else
      nt = str%cblk(nft)
    endif

    ! Fall here when there are 2 contacts for fast transmission
    if (ncont == 2 .and. size(ni) == 1 .and. nt == 1) then
      call allocate_gsm(gsmr,nbl)
      call calculate_gsmr_blocks(ESH,nbl,2,.false.)
      call allocate_blk_dns(Gr,nbl)
      call calculate_Gr_tridiag_blocks(ESH,1)
      call calculate_single_transmission_2_contacts(nit,nft,ESH,SelfEneR,str%cblk,tun_proj,tun)
      tun_mat(1) = tun

    else
      ! MULTITERMINAL case
      call allocate_gsm(gsmr,nbl)
      call calculate_gsmr_blocks(ESH,nbl,2)

      do icpl = 1, size(ni)

        !Computation of transmission(s) between contacts ni(:) -> nf(:)
        nit=ni(icpl)
        nft=nf(icpl)

        ! find the largest contact block between the two terminals
        if (str%cblk(nit).gt.str%cblk(nft)) then
          nt1 = str%cblk(nit)
        else
          nt1 = str%cblk(nft)
        endif

        if (icpl == 1) then
          ! Iterative calculation of Gr down to nt
          nt = nt1
          call allocate_blk_dns(Gr,nbl)
          call calculate_Gr_tridiag_blocks(ESH,1)
          if (nt.gt.1) then
            call calculate_Gr_tridiag_blocks(ESH,2,nt)
          end if
        else
          ! When more contacts are present sometimes we can re-use previous GF
          ! if nt1 > nt extend the Gr calculation
          if (nt1 .gt. nt) then
            call calculate_Gr_tridiag_blocks(ESH,nt+1,nt1)
            nt = nt1
          endif
        end if

        call calculate_single_transmission_N_contacts(nit,nft,ESH,SelfEneR,str%cblk,tun_proj,tun)

        tun_mat(icpl) = tun

      end do
    end if

    !Deallocate energy-dependent matrices
    call destroy_gsm(gsmr)
    call deallocate_gsm(gsmr)

    call destroy_tridiag_blk(Gr)
    call deallocate_blk_dns(Gr)

    call destroy_tridiag_blk(ESH)
    call deallocate_blk_dns(ESH)

  end subroutine calculate_transmissions

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

  subroutine calculate_single_transmission_2_contacts(ni,nf,ESH,SelfEneR,cblk,tun_proj,TUN)
    integer, intent(in) :: ni,nf
    type(z_DNS), intent(in) :: SelfEneR(MAXNCONT)
    type(z_DNS), intent(in) :: ESH(:,:)
    integer, intent(in) :: cblk(:)
    type(intarray), intent(in) :: tun_proj
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
  subroutine calculate_single_transmission_N_contacts(ni,nf,ESH,SelfEneR,cblk,tun_proj,TUN)
    integer, intent(in) :: ni,nf
    type(z_DNS), intent(in) :: SelfEneR(MAXNCONT)
    type(z_DNS), intent(in) :: ESH(:,:)
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

          call prealloc_mult(gsmr(i),ESH(i,i-1),minusOne,work1)

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

    ! Work to compute transmission matrix (Gamma2 Gr Gamma1 Ga)
    call prealloc_mult(GAM2_dns,Gr(bl2,bl1),work1)

    call destroy(GAM2_dns)

    call prealloc_mult(work1,GAM1_dns,work2)

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

  !---------------------------------------------------!
  !Subroutine for transmission and dos calculation    !
  !---------------------------------------------------!

  subroutine calculate_transmissions_and_dos(H,S,Ec,SelfEneR,Gs,ni,nf,str,tun_proj,TUN_MAT,dos_proj,LEDOS)
    Type(z_CSR), intent(in) :: H
    Type(z_CSR), intent(in) :: S
    Complex(dp), intent(in) :: Ec
    Type(z_DNS), Dimension(MAXNCONT), intent(in) :: SelfEneR, Gs
    Integer, intent(in) :: ni(:)
    Integer, intent(in) :: nf(:)
    Type(TStruct_Info), intent(in) :: str
    type(intarray), intent(in) :: tun_proj
    type(intarray), dimension(:), intent(in) :: dos_proj
    Real(dp), Dimension(:), intent(inout) :: TUN_MAT
    Real(dp), Dimension(:), intent(inout) :: LEDOS

    ! Local variables
    Type(z_CSR) :: ESH_tot, GrCSR
    Type(z_DNS), Dimension(:,:), allocatable :: ESH
    Type(r_CSR) :: Grm                          ! Green Retarded nella molecola
    real(dp), dimension(:), allocatable :: diag
    Real(dp) :: tun
    Complex(dp) :: zc
    Integer :: nbl,ncont, ierr
    Integer :: nit, nft, icpl
    Integer :: iLDOS, i2, i
    Character(1) :: Im


    nbl = str%num_PLs
    ncont = str%num_conts
    Im = 'I'

    !Calculation of ES-H and brak into blocks
    call prealloc_sum(H,S,minusOne,Ec,ESH_tot)

    call allocate_blk_dns(ESH,nbl)
    call zcsr2blk_sod(ESH_tot,ESH,str%mat_PL_start)
    call destroy(ESH_tot)

    !Inclusion of the contact Self-Energies to the relevant blocks
    do i=1,ncont
      ESH(str%cblk(i),str%cblk(i))%val = ESH(str%cblk(i),str%cblk(i))%val-SelfEneR(i)%val
    end do

    call allocate_gsm(gsmr,nbl)
    call calculate_gsmr_blocks(ESH,nbl,2)

    call allocate_blk_dns(Gr,nbl)
    ! call create
    call calculate_Gr_tridiag_blocks(ESH,1)
    call calculate_Gr_tridiag_blocks(ESH,2,nbl)

    !Computation of transmission(s) between contacts ni(:) -> nf(:)
    do icpl=1,size(ni)

      nit=ni(icpl)
      nft=nf(icpl)

      select case(ncont)
      case(1)
        tun = 0.0_dp
      case(2)
        call calculate_single_transmission_2_contacts(nit,nft,ESH,SelfEneR,str%cblk,tun_proj,tun)
      case default
        call calculate_single_transmission_N_contacts(nit,nft,ESH,SelfEneR,str%cblk,tun_proj,tun)
      end select

      TUN_MAT(icpl) = tun

    end do

    !Deallocate energy-dependent matrices
    call destroy_gsm(gsmr)
    call deallocate_gsm(gsmr)
    call destroy_tridiag_blk(ESH)
    call deallocate_blk_dns(ESH)

    ! Destroy only off-diagonal blocks
    do i=2,nbl
      call destroy(Gr(i-1,i))
      call destroy(Gr(i,i-1))
    end do

    call create(Grm,str%mat_PL_start(nbl+1)-1,str%mat_PL_start(nbl+1)-1,0)

    Grm%rowpnt(:)=1

    do i=1,nbl
      call create(GrCSR,Gr(i,i)%nrow,Gr(i,i)%ncol,Gr(i,i)%nrow*Gr(i,i)%ncol)
      call dns2csr(Gr(i,i),GrCSR)
      !Concatena direttamente la parte immaginaria per il calcolo della doS
      zc=minusOne/pi

      call concat(Grm,zc,GrCSR,Im,str%mat_PL_start(i),str%mat_PL_start(i))
      call destroy(Gr(i,i))
      call destroy(GrCSR)
    end do

    call deallocate_blk_dns(Gr)

    !Compute LDOS on the specified intervals
    if (size(dos_proj).gt.0) then
      call log_allocate(diag, Grm%nrow)
      call getdiag(Grm,diag)
      if (size(dos_proj) == size(diag) .and. dos_proj(1)%indexes(1) == 0) then
          LEDOS(:) = LEDOS(:) + diag(:)
      else
        do iLDOS=1,size(dos_proj)
          do i = 1, size(dos_proj(iLDOS)%indexes)
            i2 = dos_proj(iLDOS)%indexes(i)
            if (i2 .le. str%central_dim) then
              LEDOS(iLDOS) = LEDOS(iLDOS) + diag(i2)
            end if
          end do
        end do
      endif
      call log_deallocate(diag)
    endif

    call destroy(Grm)

  end subroutine calculate_transmissions_and_dos


  !---------------------------------------------------


  subroutine allocate_gsm(gsm,nbl)
    type(z_DNS), dimension(:), allocatable :: gsm
    integer :: nbl, ierr

    if (.not.allocated(gsm)) then
      allocate(gsm(nbl),stat=ierr)
    end if
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate gsm'

  end subroutine allocate_gsm

  !---------------------------------------------------

  subroutine allocate_blk_dns(blkM,nbl)
    type(z_DNS), dimension(:,:), allocatable :: blkM
    integer :: nbl, ierr

    allocate(blkM(nbl,nbl),stat=ierr)
    if (ierr.ne.0) stop 'ALLOCATION ERROR: could not allocate block-Matrix'

  end subroutine allocate_blk_dns

  !---------------------------------------------------

  subroutine deallocate_gsm(gsm)
    type(z_DNS), dimension(:), allocatable :: gsm
    integer :: ierr

    deallocate(gsm,stat=ierr)
    if (ierr.ne.0) stop 'DEALLOCATION ERROR: could not deallocate gsmr'

  end subroutine deallocate_gsm

  !---------------------------------------------------

  subroutine deallocate_blk_dns(blkM)
    type(z_DNS), dimension(:,:), allocatable :: blkM
    integer :: ierr

    deallocate(blkM,stat=ierr)
    if (ierr.ne.0) stop 'DEALLOCATION ERROR: could not deallocate block-Matrix'

  end subroutine deallocate_blk_dns

end module iterative


  !****************************************************************************
  !
  ! Calculate G_p=iG> contributions due to el-ph
  ! Writing on memory
  !
  !****************************************************************************
!!$  subroutine calculate_Gp_ph(negf,ESH,iter,Gp)
!!$
!!$    type(Tnegf) :: negf
!!$    type(z_DNS), dimension(:,:) :: ESH, Gp
!!$    integer :: iter
!!$
!!$    Type(z_DNS), dimension(:,:), allocatable :: Sigma_ph_p
!!$    Type(z_DNS) :: Ga, work1, work2
!!$    integer :: n, k, nbl, nrow, ierr
!!$
!!$    nbl = negf%str%num_PLs
!!$    ALLOCATE(Sigma_ph_p(nbl,nbl),stat=ierr)
!!$    if (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_n'
!!$
!!$
!!$    do n = 1, nbl
!!$
!!$      nrow = ESH(n,n)%nrow
!!$
!!$      call create(Sigma_ph_p(n,n), nrow, nrow)
!!$
!!$      Sigma_ph_p(n,n)%val = (0.0_dp, 0.0_dp)
!!$      if (iter .gt. 0) then
!!$         call read_blkmat(Sigma_ph_p(n,n),negf%scratch_path,'Sigma_ph_p_',n,n,negf%iE)
!!$      else
!!$         call write_blkmat(Sigma_ph_p(n,n),negf%scratch_path,'Sigma_ph_p_',n,n,negf%iE)
!!$      endif
!!$
!!$    end do
!!$
!!$    do n = 1, nbl
!!$
!!$      do k = 1, nbl
!!$
!!$         if (Gr(n,k)%nrow.gt.0) then
!!$            call zdagger(Gr(n,k),Ga)
!!$            call prealloc_mult(Gr(n,k), Sigma_ph_p(k,k), work1)
!!$            call prealloc_mult(work1, Ga, work2)
!!$            Gp(n,n)%val = Gp(n,n)%val + work2%val
!!$            call destroy(work1, work2, Ga)
!!$         endif
!!$
!!$      end do
!!$
!!$    end do
!!$
!!$    do n = 1, nbl
!!$      call destroy(Sigma_ph_p(n,n))
!!$    end do
!!$
!!$    DEALLOCATE(Sigma_ph_p)
!!$
!!$  end subroutine calculate_Gp_ph
