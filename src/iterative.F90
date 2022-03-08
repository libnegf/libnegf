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
  use ln_constants, only : pi
  use ln_allocation
  use mat_def
  use sparsekit_drv
  use inversions
  use elph
  use ln_structure, only : TStruct_Info
  use lib_param, only : MAXNCONT, Tnegf, intarray
  use mpi_globals, only : id, numprocs, id0                                 !DAR
  use outmatrix, only : outmat_c, inmat_c, direct_out_c, direct_in_c
  use clock
  !use transform

  implicit none
  private

  public :: calculate_transmissions
  public :: calculate_transmissions_and_dos

  public :: calculate_Gr
  public :: calculate_Gn_neq_components

  public :: sigma_ph_n
  public :: sigma_ph_p
  public :: sigma_ph_r
  public :: sigma_ph_r_z
  public :: check_sigma_ph_r
  public :: check_Gl_Gr

  public :: rebuild_dns
  public :: calculate_gsmr_blocks
  public :: calculate_gsml_blocks
  public :: calculate_Gr_tridiag_blocks
  public :: calculate_Gr_column_blocks
  public :: calculate_Gn_tridiag_blocks
  public :: calculate_Gr_outer
  public :: calculate_Gn_outer
  !public :: Outer_A_dns
  public :: iterative_meir_wingreen
  public :: transmission_BP_corrected

  public :: complete_sigma_ph_r

  public :: create_scratch
  public :: destroy_scratch

  logical, parameter :: debug=.false.

  type(z_DNS), dimension(:), allocatable :: gsmr
  type(z_DNS), dimension(:), allocatable :: gsml
  type(z_DNS), dimension(:,:), allocatable :: Gr

  type(z_DNS), dimension(:,:,:), allocatable, SAVE :: GGn
  type(z_DNS), dimension(:,:,:), allocatable, SAVE :: GGp
  type(z_DNS), dimension(:,:,:), allocatable, SAVE :: GGr
  type(z_DNS), dimension(:,:,:), allocatable, SAVE :: Sgn
  type(z_DNS), dimension(:,:,:), allocatable, SAVE :: Sgp
  type(z_DNS), dimension(:,:,:), allocatable, SAVE :: Sgr
  logical, parameter :: memory = .true.

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
    type(z_DNS), dimension(:,:), allocatable :: ESH
    type(z_CSR) :: ESH_tot, Ain
    integer :: i,ierr, nbl, ncont,ii,n

    nbl = negf%str%num_PLs
    ncont = negf%str%num_conts

    ! Take CSR H,S and build ES-H in dense blocks
    call prealloc_sum(negf%H,negf%S,(-1.0_dp, 0.0_dp),E,ESH_tot)

    call allocate_blk_dns(ESH,nbl)

    call zcsr2blk_sod(ESH_tot, ESH, negf%str%mat_PL_start)

    call destroy(ESH_tot)

    associate(cblk=>negf%str%cblk)
    do i=1,ncont
      ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
    end do
    end associate

    !! Add interaction self energy contribution, if any
    if (allocated(negf%inter)) call negf%inter%add_sigma_r(ESH)

    call allocate_gsm(gsmr,nbl)
    call calculate_gsmr_blocks(ESH,nbl,2)

    call allocate_blk_dns(Gr,nbl)

    call calculate_Gr_tridiag_blocks(ESH,1)
    call calculate_Gr_tridiag_blocks(ESH,2,nbl)

    call destroy_ESH(ESH)
    call deallocate_blk_dns(ESH)

    call destroy_gsm(gsmr)
    call deallocate_gsm(gsmr)


    ! SAVE ON FILES/MEMORY (for elph).........................
    if (negf%elph%numselmodes.gt.0 .and. negf%elph%model .eq. -1) then
      ! save diagonal blocks of Gn = -i G<
      do i = 1, nbl
        call write_blkmat(Gr(i,i),negf%scratch_path,'G_r_',i,i,negf%iE)
      end do
    endif
    !..........................................................
    !! Deliver Gr to interaction models if any
    if (allocated(negf%inter)) call negf%inter%set_Gr(Gr, negf%iE)
    !-----------------------------------------------------------

    call blk2csr(Gr,negf%str,negf%S,Grout)

    SELECT CASE (outer)
    CASE(0)
    CASE(1)
      call calculate_Gr_outer(Tlc,Tcl,gsurfR,negf%str,.FALSE.,Grout)
    CASE(2)
      call calculate_Gr_outer(Tlc,Tcl,gsurfR,negf%str,.TRUE.,Grout)
    end SELECT

    !Distruzione dell'array Gr
    call destroy_blk(Gr)
    DEALLOCATE(Gr)

  end subroutine calculate_Gr


  !****************************************************************************
  !
  ! Driver for computing G_n contributions due to all contacts but reference:
  !
  !   Sum   [f_j(E)-f_r(E)] Gr Gam_j Ga
  !   j!=r
  !
  !****************************************************************************

  subroutine calculate_Gn_neq_components(negf,E,SelfEneR,Tlc,Tcl,gsurfR,frm,Glout,outblocks)

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
    !*****************************************************************************

    implicit none

    !In/Out
    type(Tnegf), intent(inout) :: negf
    type(z_CSR), intent(inout)  :: Glout
    type(z_DNS), dimension(:), intent(in)  :: SelfEneR, gsurfR, Tlc, Tcl
    real(dp), intent(in)  :: E
    real(dp), dimension(:), intent(in)  :: frm
    integer, intent(in)  :: outblocks

    !Work
    integer :: ref, iter
    complex(dp) :: Ec
    integer :: i,ierr,ncont,nbl, lbl, rbl
    type(z_DNS), dimension(:,:), allocatable :: ESH
    type(z_DNS), dimension(:,:), allocatable :: Gn
    integer, dimension(:), allocatable :: Gr_columns
    type(z_CSR) :: ESH_tot, Gl
    logical :: mask(MAXNCONT)


    nbl = negf%str%num_PLs
    ncont = negf%str%num_conts
    ref = negf%refcont
    associate (cblk=>negf%str%cblk, indblk=>negf%str%mat_PL_start)

    Ec = cmplx(E,0.0_dp,dp)

    ! Take CSR H,S and build ES-H in dense blocks
    call prealloc_sum(negf%H,negf%S,(-1.0_dp, 0.0_dp),Ec,ESH_tot)

    call allocate_blk_dns(ESH,nbl)

    call zcsr2blk_sod(ESH_tot,ESH, negf%str%mat_PL_start)

    call destroy(ESH_tot)

    do i=1,ncont
      ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
    end do

    !! Add interaction self energy if any and initialize scba counter
    if (allocated(negf%inter)) then
      call negf%inter%add_sigma_r(ESH)
      if (allocated(negf%inter)) iter = negf%inter%scba_iter
    end if

    call allocate_gsm(gsmr,nbl)
    call allocate_gsm(gsml,nbl)

    ! Determine the leftmost and rightmost contact blocks to determine
    ! which column blocks are needed and hence which gsmr and gsml. In the case
    ! of interactions we need to be able to calculate all columns.
    if (allocated(negf%inter)) then
      rbl = 1
      lbl = nbl - 1
    else
      mask = .true.
      mask(ref) = .false.
      rbl = minval(cblk(1:ncont),mask(1:ncont))
      lbl = maxval(cblk(1:ncont),mask(1:ncont))
    endif

    call calculate_gsmr_blocks(ESH,nbl,rbl+1)
    call calculate_gsml_blocks(ESH,1,lbl-1)

    call allocate_blk_dns(Gr,nbl)

    ! -------------------------------------------------------------
    ! 1. rbl>lbl  => lbl+1=rbl-1 => compute first Gr(rbl-1,rbl-1)
    ! 2. rbl<lbl  => lbl=rbl-2 has been computed
    ! calculate_Gr does not compute if sbl>nbl or sbl<1
    call calculate_Gr_tridiag_blocks(ESH,rbl)
    call calculate_Gr_tridiag_blocks(ESH,rbl+1,nbl)
    call calculate_Gr_tridiag_blocks(ESH,rbl-1,1)
    !Passing Gr to interaction that builds Sigma_n
    if (allocated(negf%inter)) call negf%inter%set_Gr(Gr, negf%iE)

    !Computes the columns of Gr for the contacts != reference
    ! Keep track of the calculated column indices in the array Gr_columns.
    allocate(Gr_columns(ncont))
    Gr_columns = 0
    do i=1,ncont
      if (i.NE.ref) THEN
        Gr_columns(i) = cblk(i)
        call calculate_Gr_column_blocks(ESH,cblk(i),indblk)
      endif
    end do
    
    !If not interactions are present we can already destroy gsmr, gsml.
    !Otherwise they are still needed to calculate columns ont-the-fly.
    if (.not.allocated(negf%inter)) then
      call destroy_gsm(gsmr)
      call deallocate_gsm(gsmr)
      call destroy_gsm(gsml)
      call deallocate_gsm(gsml)
    end if

    !Computing device G_n
    call allocate_blk_dns(Gn,nbl)

    call init_blkmat(Gn,ESH)

    call calculate_Gn_tridiag_blocks(ESH,SelfEneR,frm,ref,negf%str,Gn)

    !Adding el-ph part: G^n = G^n + G^r Sigma^n G^a (at first call does nothing)
    !NOTE:  calculate_Gn has factor [f_i - f_ref], hence all terms will contain this factor
    if (allocated(negf%inter)) then
      call calculate_Gn_tridiag_elph_contributions(negf,ESH,iter,Gn,Gr_columns)
      call destroy_gsm(gsmr)
      call deallocate_gsm(gsmr)
      call destroy_gsm(gsml)
      call deallocate_gsm(gsml)
    end if

    !Passing G^n to interaction that builds Sigma^n
    if (allocated(negf%inter)) call negf%inter%set_Gn(Gn, negf%iE)

    call blk2csr(Gn,negf%str,negf%S,Glout)

    end associate

    !Computing the 'outer' blocks (device/contact overlapping elements)
    SELECT CASE (outblocks)
    CASE(0)
    CASE(1)
      call calculate_Gn_outer(Tlc,gsurfR,SelfEneR,negf%str,frm,ref,.false.,Glout)
    CASE(2)
      call calculate_Gn_outer(Tlc,gsurfR,SelfEneR,negf%str,frm,ref,.true.,Glout)
    end SELECT

    call destroy_blk(Gn)
    DEALLOCATE(Gn)

    call destroy_blk(Gr)
    DEALLOCATE(Gr)

    call destroy_ESH(ESH)
    DEALLOCATE(ESH)

  end subroutine calculate_Gn_neq_components

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
  !  It implements the form without the reference electrode:
  !
  !  I_i = Tr[\Sigma_{i}^{n}A - \Gamma_{i}G^{n}] =
  !      = Tr[\Gamma_{i}(f_{i}A - f_{ref}A - G^{n,l\neq ref} - G^{n,l\neq ref}_{\phi})]
  !
  !  where G^{n,l\neq ref} is the component including no el-ph
  !  If i=ref it reduces to
  !  I_i = Tr[ \Gamma_{i}(-G^{n,l\neq ref}  - G^{n,l\neq ref}_{\phi}) ]
  !---------------------------------------------------------------------
  subroutine iterative_meir_wingreen(negf,E,SelfEneR,Tlc,Tcl,gsurfR,frm,curr_mat)
    type(Tnegf), intent(inout) :: negf
    real(dp) :: E
    type(z_DNS), dimension(:) :: SelfEneR, gsurfR, Tlc, Tcl
    real(dp), dimension(:) :: frm
    real(dp), dimension(:) :: curr_mat

    !Work
    complex(dp) :: Ec, tmp
    integer :: i,ierr,ncont,nbl,lbl,ii,n
    integer :: pl_start, pl_end
    integer :: ref, iter, lead, lead_blk, ref_blk
    type(z_DNS), dimension(:,:), allocatable :: ESH, Gn
    integer, dimension(:), allocatable :: Gr_columns
    type(z_DNS) :: work1, work2, Gam, A
    type(z_CSR) :: ESH_tot, Gl

    nbl = negf%str%num_PLs
    ncont = negf%str%num_conts
    ref = negf%refcont
    ref_blk = negf%str%cblk(ref)
    associate (cblk=>negf%str%cblk)

    iter = 0
    if (allocated(negf%inter)) then
      iter = negf%inter%scba_iter
    end if

    Ec=cmplx(E,0.0_dp,dp)

    call prealloc_sum(negf%H,negf%S,(-1.0_dp, 0.0_dp),Ec,ESH_tot)

    call allocate_blk_dns(ESH,nbl)

    call zcsr2blk_sod(ESH_tot, ESH, negf%str%mat_PL_start)

    call destroy(ESH_tot)

    ! Add contact self energies
    do i=1,ncont
      ESH(cblk(i),cblk(i))%val = ESH(cblk(i),cblk(i))%val-SelfEneR(i)%val
    end do

    !! Add el-ph self energy if any
    if (allocated(negf%inter)) call negf%inter%add_sigma_r(ESH)

    call allocate_gsm(gsmr,nbl)
    call allocate_gsm(gsml,nbl)

    call calculate_gsmr_blocks(ESH,nbl,2)

    if ( ncont.gt.1 ) then
      call calculate_gsml_blocks(ESH,1,nbl-2)
    endif

    call allocate_blk_dns(Gr,nbl)

    call calculate_Gr_tridiag_blocks(ESH,1)
    call calculate_Gr_tridiag_blocks(ESH,2,nbl)

    !! Give Gr to interaction model if any
    if (allocated(negf%inter)) call negf%inter%set_Gr(Gr, negf%iE)

    !! Never calculate outer blocks
    call allocate_blk_dns(Gn, nbl)

    call init_blkmat(Gn,ESH)

    !! TEMPORARY AND INEFFICIENT:
    !! CALCULATE THE FULL Gn WHEN THE CONTACT BLOCK WHOULD BE ENOUGH
    allocate(Gr_columns(ncont))
    Gr_columns = 0
    do i=1,ncont
      if (i.NE.ref) THEN
        Gr_columns(i) = cblk(i)
        call calculate_Gr_column_blocks(ESH, cblk(i), negf%str%mat_PL_start)
      endif
    end do
    call calculate_Gn_tridiag_blocks(ESH,SelfEneR,frm,ref,negf%str,Gn)

    if (allocated(negf%inter)) then
      call calculate_Gn_tridiag_elph_contributions(negf,ESH,iter,Gn, Gr_columns)
    end if

    ! The gsmr, gsml are used to calculate columns on-the-fly in calculate_Gn_tridiag_elph_contributions, we
    ! can destroy them here.
    call destroy_gsm(gsmr)
    call deallocate_gsm(gsmr)
    call destroy_gsm(gsml)
    call deallocate_gsm(gsml)

    do i=1,size(negf%ni)
      lead = negf%ni(i)
      lead_blk = cblk(lead)
      call zspectral(SelfEneR(lead),SelfEneR(lead), 0, Gam)
      if (lead.eq.ref) then
        call prealloc_mult(Gam, Gn(lead_blk, lead_blk), (-1.0_dp, 0.0_dp), work1)
        curr_mat(i) = real(trace(work1))
        call destroy(work1)
      else
        call zspectral(Gr(lead_blk, lead_blk), Gr(lead_blk, lead_blk), 0, A)
        tmp = frm(lead)-frm(ref)
        call prealloc_sum(A, Gn(lead_blk, lead_blk), tmp, (-1.0_dp, 0.0_dp), work1)
        call destroy(A)
        call prealloc_mult(Gam, work1, work2)
        call destroy(work1)
        curr_mat(i) = real(trace(work2))
        call destroy(work2)
      endif
      call destroy(Gam)
    end do
    !Convert to output CSR format.
    call blk2csr(Gn,negf%str,negf%S,Gl)
    DEALLOCATE(Gn)

    call destroy_blk(Gr)
    DEALLOCATE(Gr)

    call destroy_ESH(ESH)
    DEALLOCATE(ESH)
      
    end associate

  end subroutine iterative_meir_wingreen

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
         Gam1%val(nn,nn)=negf%deph%bp%coupling(nn)
      end if

      do mm = 1, NumOrbs+ncont

        if (mm > NumOrbs) then
           cont = mm - NumOrbs
           call zspectral(SelfEneR(cont),SelfEneR(cont),0,Tmp2)
           Gam2%val(pls(cont):ple(cont),pls(cont):ple(cont)) = Tmp2%val
           call destroy(Tmp2)
        else
           Gam2%val = 0.0_dp
           Gam2%val(mm,mm)=negf%deph%bp%coupling(mm)
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
  subroutine init_blkmat(Matrix,S)
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

  end subroutine init_blkmat


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
    type(z_DNS), dimension(:,:) :: M
    integer :: i, i1, nbl

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
  subroutine destroy_ESH(ESH)

    integer :: i, nbl
    type(z_DNS), dimension(:,:) :: ESH

    nbl=size(ESH,1)

    do i=1,nbl
      call destroy(ESH(i,i))
    end do
    do i=2,nbl
      call destroy(ESH(i-1,i))
      call destroy(ESH(i,i-1))
    end do

  end subroutine destroy_ESH


  !***********************************************************************
  !> Update the value of el-ph Retarded  from block matrix
  !  Self Energy for the dephasing model
  !
  !***********************************************************************
  subroutine update_elph_r(negf, Gr)
    type(Tnegf), intent(inout) :: negf
    type(z_DNS), dimension(:,:), intent(in) :: Gr

    type(z_DNS) :: work1, work2, work3
    integer :: npl, nblk, n, ii, jj, norbs, indstart, indend, nnz

    npl = negf%str%num_PLs

    select case(negf%elph%model)
    case (0)
      return
    case (1)
      stop "Deprecated"
    case(2)
      stop "Deprecated"
    case(3)
      stop "Deprecated"
    case default
      stop 'Elph model not yet implemented'
    end select

  end subroutine update_elph_r

  !***********************************************************************
  !> Update the value of el-ph inscattering from block matrix
  !  Self Energy for the dephasing model
  !
  !***********************************************************************
  subroutine update_elph_n(negf, Gn)
    type(Tnegf), intent(inout) :: negf
    type(z_DNS), dimension(:,:), intent(in) :: Gn

    type(z_DNS) :: work1, work2, work3
    integer :: npl, nblk, n, ii, jj, norbs, indstart, indend, nnz

    npl = negf%str%num_PLs

    select case(negf%elph%model)
    case (0)
      return
    case (1)
      stop "Deprecated"
    case(2)
      stop "Deprecated"
    case(3)
      stop "Deprecated"
    case default
      stop 'Elph model not yet implemented'
    end select

  end subroutine update_elph_n

  !***********************************************************************
  !
  !  Reloads the retarded elph self-energy and add it to ES-H
  !
  !***********************************************************************
  subroutine add_elph_sigma_r(negf, ESH, elph)
    type(Tnegf), intent(in) :: negf
    type(Telph), intent(in) :: elph
    type(z_DNS), dimension(:,:), intent(inout) :: ESH

    type(z_DNS), dimension(:,:), allocatable :: Sigma_ph_r, sigma_blk
    integer :: n, nbl, nrow, ierr, ii, jj, nblk, norbs, indstart, indend

    nbl = negf%str%num_PLs

    ! At first loop we should add no
    if (elph%scba_iter .eq. 0 .or. elph%model .eq. 0) then
      return
    end if

    select case(elph%model)

    case(0)
      return
    case(1)
      stop "Deprecated"
    case (2)
      stop "Deprecated"
    case(3)
      stop "Deprecated"
    ! Old Alex one
    case(-1)
      nbl = negf%str%num_PLs
      allocate(Sigma_ph_r(nbl,nbl),stat=ierr)
      if (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_r'

      do n = 1, nbl
        nrow = ESH(n,n)%nrow
        call create(Sigma_ph_r(n,n), nrow, nrow)
        Sigma_ph_r(n,n)%val = (0.0_dp, 0.0_dp)
        if (elph%scba_iter .gt. 0) then
          call read_blkmat(Sigma_ph_r(n,n),negf%scratch_path,'Sigma_ph_r_',n,n,negf%iE)
        else
          call write_blkmat(Sigma_ph_r(n,n),negf%scratch_path,'Sigma_ph_r_',n,n,negf%iE)
        end if
        ESH(n,n)%val = ESH(n,n)%val - Sigma_ph_r(n,n)%val
        call destroy(Sigma_ph_r(n,n))
      end do
      deallocate(Sigma_ph_r)
    case default
      stop 'Elph model not yet implemented'
    end select

  end subroutine add_elph_sigma_r

  !***********************************************************************
  !
  !  Construct a sparse matrix starting from the sparse matrices array
  !  related to the blocks
  !
  !***********************************************************************

  subroutine rebuild_dns(Atot,A,n,indb)

    !***********************************************************************
    !Input:
    !A: sparse matrices array
    !n: A dimension (A(n,n))
    !indb: blocks indeces array (contains position of first elements of diagonal
    !      blocks, included an additional ending address "last row+1")
    !
    !Output:
    !Atot: total sparse matrix (allocated internally)
    !***********************************************************************

    implicit none

    !In/Out
    type(z_CSR) :: Atot  !Allocato internamente
    integer :: n
    integer, dimension(n+1) :: indb
    type(z_cSR), dimension(n,n) :: A

    !Work
    integer :: i,j,Atot_nrow,i1,j1

    Atot_nrow=indb(n+1)-1

    call create(Atot,Atot_nrow,Atot_nrow,0)
    Atot%rowpnt(:)=1

    do i=1,n
      do j=1,n

        if (A(i,j)%nrow.GT.0) THEN
          i1=indb(i)
          j1=indb(j)
          call concat(Atot,A(i,j),i1,j1)
        endif

      end do
    end do

  end subroutine rebuild_dns

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

      call prealloc_mult(ESH(i,i+1),gsmr(i+1),(-1.0_dp, 0.0_dp),work1)

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

  subroutine calculate_gsml_blocks(ESH,sbl,ebl)

    !***********************************************************************
    !Input:
    !ESH: sparse matrices array ESH(nbl,nbl)
    !
    !nbl (number of layers), indblk(nbl+1) global variables needed
    !
    !Output:
    !sparse matrices array global variable gsml(nbl) is available in memory
    !single blocks are allocated internally, array Gr(nbl,nbl)
    !must be allocated externally
    !***********************************************************************


    implicit none

    !In/Out
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: sbl,ebl                       ! start block, end block

    !Work
    type(z_DNS) :: work1, work2
    integer :: nrow
    integer :: i, nbl
    !type(z_DNS) :: INV(sbl,sbl)

    if (sbl.gt.ebl) return

    nbl = size(ESH,1)

    if (nbl.eq.1) return

    nrow=ESH(sbl,sbl)%nrow

    call create(gsml(sbl),nrow,nrow)

    call compGreen(gsml(sbl),ESH(sbl,sbl),nrow)


    do i=sbl+1,ebl

      nrow=ESH(i,i)%nrow

      call prealloc_mult(ESH(i,i-1),gsml(i-1),(-1.0_dp, 0.0_dp),work1)

      call prealloc_mult(work1,ESH(i-1,i),work2)

      call destroy(work1)

      call prealloc_sum(ESH(i,i),work2,work1)

      call destroy(work2)

      call create(gsml(i),work1%nrow,work1%nrow)

      call compGreen(gsml(i),work1,work1%nrow)

      call destroy(work1)

    end do

    if (debug) then
      WRITE(*,*) '********************'
      WRITE(*,*) 'calculate_gsml done'
      WRITE(*,*) '********************'
    endif

  end subroutine calculate_gsml_blocks





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
          call prealloc_sum(work1,work3,(-1.0_dp, 0.0_dp),work2)
          call destroy(work3)
          work1%val = work2%val
          call destroy(work2)
        endif
        if (sbl-1.ge.1) then
          call prealloc_mult(ESH(sbl,sbl-1),gsml(sbl-1),work2)
          call prealloc_mult(work2,ESH(sbl-1,sbl),work3)
          call destroy(work2)
          call prealloc_sum(work1,work3,(-1.0_dp, 0.0_dp),work2)
          call destroy(work3)
          work1%val = work2%val
          call destroy(work2)
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
        call prealloc_mult(work1,Gr(i-1,i-1),(-1.0_dp,0.0_dp),Gr(i,i-1))
        call destroy(work1)

        call prealloc_mult(ESH(i-1,i),gsmr(i),work2)
        call prealloc_mult(Gr(i-1,i-1),work2,(-1.0_dp, 0.0_dp),Gr(i-1,i))

        call prealloc_mult(Gr(i,i-1),work2,(-1.0_dp,0.0_dp),work1)
        call destroy(work2)

        call prealloc_sum(gsmr(i),work1,Gr(i,i))
        call destroy(work1)
      end do
    ELSE
      do i=sbl,ebl,-1
        call prealloc_mult(gsml(i),ESH(i,i+1),work1)
        call prealloc_mult(work1,Gr(i+1,i+1),(-1.0_dp,0.0_dp),Gr(i,i+1))
        call destroy(work1)

        call prealloc_mult(ESH(i+1,i),gsml(i),work2)
        call prealloc_mult(Gr(i+1,i+1),work2,(-1.0_dp, 0.0_dp),Gr(i+1,i))

        call prealloc_mult(Gr(i,i+1),work2,(-1.0_dp,0.0_dp),work1)
        call destroy(work2)

        call prealloc_sum(gsml(i),work1,Gr(i,i))
        call destroy(work1)
      end do
    endif

  end subroutine calculate_Gr_tridiag_blocks

  !**************************************************************************
  !
  !  Calculate Green Retarded column "n" - writing on memory
  !
  !**************************************************************************

  subroutine calculate_Gr_column_blocks(ESH,n,indblk)

    !***********************************************************************
    !Input:
    !ESH: sparse matrices array ESH(nbl,nbl)
    !n: n umber of column to be calculated
    !
    !global variables needed: nbl (number of layers), indblk(nbl+1),
    !Gr diagonal, subadiagonal and superdiagonal, gsmr(:) for
    !downgoing and gsml(:) for upgoing
    !
    !Output:
    !sparse matrices array global variable Gr(:,n) is available in
    !memory - single blocks are allocated internally, array Gr(nbl,nbl)
    !must be allocated externally
    !***********************************************************************

    implicit none

    !In/Out
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    integer, intent(in) :: n
    integer, dimension(:), intent(in) :: indblk

    !Work
    integer :: i,nrow,ncol,nbl
    type(z_DNS) :: work1
    real(dp) :: max

    nbl = size(ESH,1)

    if (n.GT.nbl) THEN
      STOP 'Error in calculate_Grcol : n is greater than nbl'
    endif

    !***************************************
    !  Downgoing (j>=n+2 && n<nbl-1)
    !
    !   G_j,n = -gR_jj T_j,j-1 G_j-1,n
    !
    !***************************************
    if (n.LT.(nbl-1)) THEN

      do i=n+2,nbl

        max=MAXVAL(ABS(Gr(i-1,n)%val))
        if (max.GT.EPS) THEN
          call prealloc_mult(gsmr(i),ESH(i,i-1),(-1.0_dp, 0.0_dp),work1)
          call prealloc_mult(work1,Gr(i-1,n),Gr(i,n))
          call destroy(work1)
        else
          ! WHEN BLOCK IS SMALLER THAN EPS IT IS NOT CREATED
          exit
        end if

      end do

    endif
    !*************************************
    !   Up-going (j<=n-2 && n>2)
    !
    !   G_j,n = -gL_jj T_j,j+1 G_j+1,n
    !
    !*************************************

    if (n.GT.2) THEN

      do i=n-2,1,(-1)

        max=MAXVAL(ABS(Gr(i+1,n)%val))

        if (max.GT.EPS) THEN
          call prealloc_mult(gsml(i),ESH(i,i+1),(-1.0_dp, 0.0_dp),work1)
          call prealloc_mult(work1,Gr(i+1,n),Gr(i,n))
          call destroy(work1)
        else
          ! WHEN BLOCK IS SMALLER THAN EPS IT IS NOT CREATED
          exit
        endif

      end do

    endif

    if (debug) then
      WRITE(*,*) '******************************'
      WRITE(*,*) 'calculate_Grcol done column',n
      WRITE(*,*) '******************************'
    endif

  end subroutine calculate_Gr_column_blocks

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


  !****************************************************************************
  !
  ! Calculate G_n contributions for all contacts (except reference)
  !
  !   Sum   [f_j(E)-f_r(E)] Gr Gam_j Ga
  !   j!=r
  !
  ! TRICK: in order to have the usual Gn set ref>ncont and make sure f_r = 0
  !
  !****************************************************************************

  subroutine calculate_Gn_tridiag_blocks(ESH,SelfEneR,frm,ref,struct,Gn)

    !******************************************************************************
    !Input:
    !ESH(nbl,nbl): sparse matrices array ES-H
    !SelfEneR(ncont): sparse matrices array with contacts Self Energy
    !frm(ncont): Fermi distribution value for each contact
    !ref:  reference contact
    !
    !global variables needed: nbl, indblk(nbl+1), ncont, cblk(ncont), Gr(:,:)
    !diagonal, subdiagonal, overdiagonal and in colums Gr(:,cb) where cb are
    !layers interacting with all contacts but collector
    !
    !Output:
    !Gl: sparse matrix containing G  contacts contribution
    !*******************************************************************************

    implicit none

    !In/Out
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    type(z_DNS), dimension(:,:), intent(inout) :: Gn
    type(z_DNS), dimension(:), intent(in) :: SelfEneR
    type(Tstruct_info), intent(in) :: struct
    real(dp), dimension(:), intent(in) :: frm
    integer, intent(in) :: ref

    !Work
    Type(z_DNS) :: Gam
    type(z_DNS) :: work1,Ga
    integer :: i,j,cb
    integer :: ncont, nbl
    complex(dp) :: frmdiff

    ncont = struct%num_conts
    nbl = struct%num_PLs

    !*******************************************
    ! Contact Iteration
    !*******************************************
    do j=1,ncont

      ! NOTE: this soubroutine uses prealloc_mult that performs
      ! C = C + A*B
      if (j.NE.ref .AND. ABS(frm(j)-frm(ref)).GT.EPS) THEN

        cb=struct%cblk(j) ! block corresponding to contact j

        call zspectral(SelfEneR(j),SelfEneR(j),0,Gam)

        frmdiff = cmplx(frm(j)-frm(ref),0.0_dp,dp)
        ! Computation of Gl(1,1) = Gr(1,cb) Gam(cb) Ga(cb,1)
        if (allocated(Gr(1,cb)%val)) then
          call prealloc_mult(Gr(1,cb),Gam,work1)
          call zdagger(Gr(1,cb),Ga)
          call prealloc_mult(work1,Ga,frmdiff,Gn(1,1))
          call destroy(work1, Ga)
        else
          Gn(1,1)%val=(0.0_dp,0.0_dp)
        endif

        ! Computation of all tridiagonal blocks
        do i=2,nbl

          ! Computation of Gl(i,j) = Gr(i,cb) Gam(cb) Ga(cb,j)
          ! Both Gr(i,cb) and Gr(j,cb) must be non-zero
          if (Gr(i-1,cb)%nrow.gt.0 .and. Gr(i,cb)%nrow.gt.0) then
            call prealloc_mult(Gr(i-1,cb),Gam,work1)
            call zdagger(Gr(i,cb),Ga)
            call prealloc_mult(work1,Ga,frmdiff,Gn(i-1,i))
            call destroy(work1)

            call prealloc_mult(Gr(i,cb),Gam,work1)
            call prealloc_mult(work1,Ga,frmdiff,Gn(i,i))

            call destroy(work1, Ga)

            call prealloc_mult(Gr(i,cb),Gam,work1)
            call zdagger(Gr(i-1,cb),Ga)
            call prealloc_mult(work1,Ga,frmdiff,Gn(i,i-1))

            call destroy(work1, Ga)
          else
            Gn(i  ,i)%val=(0.0_dp,0.0_dp)
            Gn(i-1,i)%val=(0.0_dp,0.0_dp)
            Gn(i,i-1)%val=(0.0_dp,0.0_dp)
          endif

        end do

        call destroy(Gam)

      endif

    end do

  end subroutine calculate_Gn_tridiag_blocks

  subroutine calculate_gns(ESH,SelfEneR,struct,gns)
    implicit none

    !In/Out
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    type(z_DNS), dimension(:,:), intent(inout) :: gns
    type(z_DNS), dimension(:), intent(in) :: SelfEneR
    type(Tstruct_info), intent(in) :: struct

    !Work
    type(z_DNS) :: work, work1, work2, work3, factors 
    type(z_DNS) :: Ga, gsmrDag, ESHdag
    integer :: i,j,cb
    integer :: nbl

    nbl = struct%num_PLs
    
    !g^n(nbl) = Gr(nbl,nbl) Sigma(nbl,nbl) Ga(nbl,nbl)
    call zdagger(Gr(nbl,nbl),Ga)
    call prealloc_mult(Gr(nbl,nbl), SelfEneR(nbl,nbl), work1)
    call prealloc_mult(work1, Ga, gns(nbl))
    call destroy(Ga, work1)

    !gns(i) = gsmr(i) * [Sigma(i,i) + ESH(i,i+1) gns(i+1) ESH^dag(i+1,i) - ESH(i,i+1) gsmr(i+1) Sigma(i+1,i) - 
    !                    Sigma(i,i+1) gsmr^dag(i+1) ESH^dag(i+1,i)] * gsmr^dag(i)
    do i = 1, nbl-1, -1
        !work1 = ESH(i,i+1) gns(i+1) ESH^dag(i+1,i)
        call zdagger(ESH(i,i+1), ESHdag)
        call prealloc_mult(ESH(i+1,i), gns(i+1), work)
        call prealloc_mult(work, ESHdag, work1)
        call destroy(work)

        !work2 = Sigma(i,i+1) gsmr^dag(i+1) ESH^dag(i+1,i)
        call zdagger(gsmr(i+1), gsmrDag)
        call prealloc_mult(SelfEneR(i,i+1), gsmrDag, work)
        call prealloc_mult(work, ESHdag, work2)
        call destroy(work, gsmrDag, ESHdag)

        !work3 = ESH(i,i+1) gsmr(i+1) Sigma(i+1,i)
        call prealloc_mult(ESH(i,i+1), gsmr(i+1), work)
        call prealloc_mult(work, SelfEneR(i+1,i), work3)
        call destroy(work)

        !sum of the four factors
        call create(factors, SelfEneR(i,i)%nrow, SelfEneR(i,i)%ncol)
        factors%val = SelfEneR(i,i)%val + work1%val - work2%val - work3%val

        !gsmr(i) * factors * gsmr^dag(i)
        call zdagger(gsmr(i), gsmrDag)
        call prealloc_mult(gsmr(i), factors, work)
        call prealloc_mult(work, gsmrDag, gns(i))

        call destroy(factors, work, gsmrDag)
    end do 

  end subroutine calculate_gns
  !****************************************************************************
  !
  ! Calculate G_n contributions for all contacts (except reference)
  ! This version computes Grcol on the fly
  !
  !****************************************************************************

  subroutine calculate_Gn_tridiag_blocks2(ESH,SelfEneR,frm,ref,struct,Gn)

    !******************************************************************************
    !Input:
    !ESH(nbl,nbl): sparse matrices array ES-H
    !SelfEneR(ncont): sparse matrices array with contacts Self Energy
    !frm(ncont): Fermi diistribution value for each contact
    !ref:  reference contact
    !
    !global variables needed: nbl, indblk(nbl+1), ncont, cblk(ncont), Gr(:,:)
    !diagonal, subdiagonal, overdiagonal and in colums Gr(:,cb) where cb are
    !layers interacting with all contacts but collector
    !
    !Output:
    !Gl: sparse matrix containing G  contacts contribution
    !*******************************************************************************

    implicit none

    !In/Out
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    type(z_DNS), dimension(:,:), intent(inout) :: Gn
    type(z_DNS), dimension(:), intent(in) :: SelfEneR
    type(Tstruct_info), intent(in) :: struct
    real(dp), dimension(:), intent(in) :: frm
    integer, intent(in) :: ref

    !Work
    Type(z_DNS) :: Gam
    type(z_DNS) :: work1,Ga
    integer :: i,j,cb
    integer :: ncont, nbl
    complex(dp) :: frmdiff

    ncont = struct%num_conts
    nbl = struct%num_PLs

    !*******************************************
    ! Contact Iteration
    !*******************************************
    do j=1,ncont

      ! NOTE: this soubroutine uses prealloc_mult that performs
      ! C = C + A*B
      if (j.NE.ref .AND. ABS(frm(j)-frm(ref)).GT.EPS) THEN

        cb=struct%cblk(j) ! block corresponding to contact j

        call zspectral(SelfEneR(j),SelfEneR(j),0,Gam)

        frmdiff = cmplx(frm(j)-frm(ref),0.0_dp,dp)
        ! Computation of Gl(1,1) = Gr(1,cb) Gam(cb) Ga(cb,1)
        if (Gr(1,cb)%nrow.gt.0) then
          call prealloc_mult(Gr(1,cb),Gam,work1)
          call zdagger(Gr(1,cb),Ga)
          call prealloc_mult(work1,Ga,frmdiff,Gn(1,1))
          call destroy(work1, Ga)
        else
          Gn(1,1)%val=(0.0_dp,0.0_dp)
        endif

        ! Computation of all tridiagonal blocks
        do i=2,nbl

          ! Computation of Gr(i,cb) assuming Gr(i-1,cb) exists
          ! Assume downgoing: i > cb
          if (Gr(i-1,cb)%nrow.GT.0) THEN
            call prealloc_mult(gsmr(i),ESH(i,i-1),(-1.0_dp, 0.0_dp),work1)
            call destroy(gsmr(i))
            call prealloc_mult(work1,Gr(i-1,cb),Gr(i,cb))
            call destroy(work1)
            if (MAXVAL(ABS(Gr(i,cb)%val)).lt.EPS) call destroy(Gr(i,cb))
          endif

          ! Computation of Gl(i,j) = Gr(i,cb) Gam(cb) Ga(cb,j)
          ! Both Gr(i,cb) and Gr(j,cb) must be non-zero
          if (Gr(i-1,cb)%nrow.gt.0 .and. Gr(i,cb)%nrow.gt.0) THEN
            call prealloc_mult(Gr(i-1,cb),Gam,work1)
            call zdagger(Gr(i,cb),Ga)
            call prealloc_mult(work1,Ga,frmdiff,Gn(i-1,i))
            call destroy(work1)

            call prealloc_mult(Gr(i,cb),Gam,work1)
            call prealloc_mult(work1,Ga,frmdiff,Gn(i,i))

            call destroy(work1, Ga)

            call prealloc_mult(Gr(i,cb),Gam,work1)
            call zdagger(Gr(i-1,cb),Ga)
            call prealloc_mult(work1,Ga,frmdiff,Gn(i,i-1))

            call destroy(work1, Ga)
          ELSE
            Gn(i  ,i)%val=(0.0_dp,0.0_dp)
            Gn(i-1,i)%val=(0.0_dp,0.0_dp)
            Gn(i,i-1)%val=(0.0_dp,0.0_dp)
          endif

          if (Gr(i-1,cb)%nrow.gt.0) call destroy(Gr(i-1,cb))

        end do

        call destroy(Gam)

      endif

    end do

  end subroutine calculate_Gn_tridiag_blocks2

  !****************************************************************************
  !
  ! Calculate G_n contributions due to elph:  G_n = G_n + Gr Sigma_ph Ga
  !
  !****************************************************************************
  subroutine calculate_Gn_tridiag_elph_contributions(negf, ESH, iter, Gn, existing_Gr_cols)

    type(Tnegf), intent(in) :: negf
    type(z_DNS), dimension(:,:), intent(in) :: ESH
    type(z_DNS), dimension(:,:), intent(inout) :: Gn
    integer, dimension(:), intent(in) :: existing_Gr_cols
    integer, intent(in) :: iter

    Type(z_DNS), dimension(:,:), allocatable :: Sigma_ph_n, sigma_blk
    Type(z_DNS) :: Ga, work1, work2, sigma_tmp
    integer :: n, k, nbl, nrow, ierr, ii, jj, norbs, nblk, indstart, indend

    nbl = negf%str%num_PLs

    !! If this is the first scba cycle, there's nothing to do
    if (negf%inter%scba_iter .eq. 0) then
      return
    endif
    allocate(Sigma_ph_n(nbl,nbl),stat=ierr)
    if (ierr.NE.0) STOP 'ALLOCATION ERROR: could not allocate Sigma_ph_n'

    ! The block sigma n is made available from el-ph model
    ! Note: the elph models could not keep a copy and calculate it
    ! on the fly. You have to rely on the local copy
    if (allocated(negf%inter)) then
      call negf%inter%get_sigma_n(Sigma_ph_n, negf%ie)
    end if


    !! Calculate the diagonal and off diagonal (if needed) blocks of Gn
    !! in the assumption of diagonal self energy
    !! Gn(n,n) = Gr(n,k)*Sigma_n(k,k)*Ga(k,n)
    !! Gn(n,n+1) = Gr(n,k)*Sigma_n(k,k)*Ga(k,n+1)
    !! Gn(n,n-1) = Gr(n,k)*Sigma_n(k,k)*Ga(k,n-1)
    !! All the rows of Gr need to be available

    do k = 1, nbl
      ! Calculate k-th column on-the-fly. The reference contact column
      ! might already be available. Check if the top and bottom of the
      ! column are available.
      if (all(existing_Gr_cols .ne. k)) then
        call calculate_Gr_column_blocks(ESH, k, negf%str%mat_PL_start)
      endif

      do n = 1, nbl

        ! Zero column blocks are not created at all.
        if (allocated(Gr(n, k)%val)) then
          ! Calculate diagonal blocks Gn(n, n)
          call zdagger(Gr(n, k), Ga)
          call prealloc_mult(Gr(n, k), Sigma_ph_n(k, k), work1)
          call prealloc_mult(work1, Ga, work2)
          Gn(n,n)%val = Gn(n,n)%val + work2%val
          call destroy(work2, Ga)

          ! Computing blocks of Gn(n, n - 1).
          if (n .lt. nbl .and. allocated(Gn(n, n + 1)%val) .and. allocated(Gr(n + 1, k)%val)) then
            call zdagger(Gr(n + 1, k), Ga)
            call prealloc_mult(work1, Ga, work2)
            Gn(n, n + 1)%val = Gn(n, n + 1)%val + work2%val
            call destroy(work2, Ga)
          endif

          ! Computing blocks of Gn(n, n - 1).
          if (n .gt. 1 .and. allocated(Gn(n, n - 1)%val) .and. allocated(Gr(n - 1, k)%val)) then
            call zdagger(Gr(n - 1, k), Ga)
            call prealloc_mult(work1, Ga, work2)
            Gn(n, n - 1)%val = Gn(n, n - 1)%val + work2%val
            call destroy(work2, Ga)
          endif
          call destroy(work1)
        end if
      end do

      ! Remove column blocks of Gr if they were calculated on the fly.
      do n = 1, nbl
        if (abs(n - k) .gt. 1 .and. all(existing_Gr_cols .ne. k)) then
          call destroy(Gr(n, k))
        end if
      end do

    end do

    deallocate(Sigma_ph_n)

  end subroutine calculate_Gn_tridiag_elph_contributions

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

  ! ******************************************************************************
  ! Computes Sigma_ph_n and save it file
  ! ******************************************************************************
  subroutine Sigma_ph_n(negf,Epnt)

    type(Tnegf) :: negf
    real(dp), dimension(:),allocatable :: Epnt

    !Local variables
    type(z_DNS) :: work1,work2
    type(z_DNS) :: Mq_mat
    type(z_DNS), dimension(:,:), allocatable :: Sigma_n, G_n_interP, G_n_interN

    real(dp), dimension(:), pointer :: Wq
    real(dp), dimension(:), pointer :: Nq
    real(dp), dimension(:), pointer :: Mq
    logical, dimension(:), pointer :: selmodes

    integer :: i, m, iE, i1,i2
    integer :: nummodes, numselmodes, nbl
    real(dp) :: E1, E2, En

    nbl = negf%str%num_PLs

    selmodes => negf%elph%selmodes
    Wq => negf%elph%Wq
    Mq => negf%elph%Mq
    Nq => negf%elph%Nq
    numselmodes = negf%elph%numselmodes
    nummodes = negf%elph%nummodes

    iE = negf%iE

    allocate(G_n_interP(nbl,nbl))
    allocate(G_n_interN(nbl,nbl))
    allocate(Sigma_n(nbl,nbl))

    associate(indblk=>negf%str%mat_PL_start)
    do i = 1, nbl
      m = indblk(i+1)-indblk(i)
      call create(Sigma_n(i,i), m, m)
      Sigma_n(i,i)%val = (0.0_dp, 0.0_dp)
      call create(G_n_interP(i,i), m, m)
      call create(G_n_interN(i,i), m, m)
    end do
    end associate

    do m = 1 , nummodes

      if (.not.selmodes(m)) cycle
      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(negf,Wq(m),Epnt,i1,i2,E1,E2)

      En = real(negf%Epnt)+Wq(m)

      call interpolation(i1, i2, E1, E2, En, negf%scratch_path, &
          'G_n_', G_n_interP)

      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(negf,-Wq(m),Epnt,i1,i2,E1,E2)

      En = real(negf%Epnt)-Wq(m)

      call interpolation(i1, i2, E1, E2, En, negf%scratch_path, &
          'G_n_', G_n_interN)


      do i = 1, nbl

        i1 =  G_n_interN(i,i)%nrow

        call create(work1, i1, i1)

        work1%val = (Nq(m)+1.0_dp)*G_n_interP(i,i)%val + Nq(m)* G_n_interN(i,i)%val

        call create_id(Mq_mat,i1,Mq(m))

        call prealloc_mult( Mq_mat, work1, work2)

        call destroy(work1)

        call prealloc_mult( work2, Mq_mat, work1)

        call destroy(work2)

        Sigma_n(i,i)%val = Sigma_n(i,i)%val + work1%val

        call destroy(work1, Mq_mat)

      end do

    end do

    !Throw away all non diagonal parts
    !do i = 1, nbl
    !   i1 = Sigma_n(i,i)%nrow
    !   do m = 1, i1
    !      Sigma_n(i,i)%val(1:m-1,m) = (0.0_dp,0.0_dp)
    !      Sigma_n(i,i)%val(m+1:i1,m) = (0.0_dp,0.0_dp)
    !   end do
    !end do
    !

    do i = 1, nbl
      call write_blkmat(Sigma_n(i,i),negf%scratch_path,'Sigma_ph_n_',i,i,iE)
      call destroy(Sigma_n(i,i))
      call destroy(G_n_interN(i,i))
      call destroy(G_n_interP(i,i))
    end do

    deallocate(Sigma_n,G_n_interN, G_n_interP)

  end subroutine Sigma_ph_n


  ! ******************************************************************************
  ! Computes Sigma_ph> and save it file
  ! ******************************************************************************
  subroutine Sigma_ph_p(negf,Epnt)

    type(Tnegf) :: negf
    real(dp), dimension(:),allocatable :: Epnt

    !Local variables
    type(z_DNS) :: work1,work2
    type(z_DNS) :: Mq_mat
    type(z_DNS), dimension(:,:), allocatable :: Sigma_p, G_p_interP, G_p_interN

    real(dp), dimension(:), pointer :: Wq
    real(dp), dimension(:), pointer :: Nq
    real(dp), dimension(:), pointer :: Mq
    logical, dimension(:), pointer :: selmodes

    integer :: i, m, iE, i1,i2
    integer :: nummodes, numselmodes, nbl
    real(dp) :: E1, E2, En

    nbl = negf%str%num_PLs

    selmodes => negf%elph%selmodes
    Wq => negf%elph%Wq
    Mq => negf%elph%Mq
    Nq => negf%elph%Nq
    numselmodes = negf%elph%numselmodes
    nummodes = negf%elph%nummodes

    iE = negf%iE

    allocate(G_p_interP(nbl,nbl))
    allocate(G_p_interN(nbl,nbl))
    allocate(Sigma_p(nbl,nbl))

    associate(indblk=>negf%str%mat_PL_start)
    do i = 1, nbl
      m = indblk(i+1)-indblk(i)
      call create(Sigma_p(i,i), m, m)
      Sigma_p(i,i)%val = (0.0_dp, 0.0_dp)
      call create(G_p_interP(i,i), m, m)
      call create(G_p_interN(i,i), m, m)
    end do
    end associate

    do m = 1 , nummodes

      if (.not.selmodes(m)) cycle
      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(negf,Wq(m),Epnt,i1,i2,E1,E2)

      En = real(negf%Epnt)+Wq(m)

      call interpolation(i1, i2, E1, E2, En, negf%scratch_path, &
          'G_p_', G_p_interP)

      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(negf,-Wq(m),Epnt,i1,i2,E1,E2)
      En = real(negf%Epnt)-Wq(m)
      call interpolation(i1, i2, E1, E2, En, negf%scratch_path, &
          'G_p_', G_p_interN)

      do i = 1, nbl

        i1 =  G_p_interN(i,i)%nrow

        call create(work1, i1, i1)

        work1%val = (Nq(m)+1.0_dp)*G_p_interN(i,i)%val + Nq(m)* G_p_interP(i,i)%val

        call create_id(Mq_mat,i1,Mq(m))

        call prealloc_mult( Mq_mat, work1, work2)

        call destroy(work1)

        call prealloc_mult ( work2, Mq_mat, work1)

        call destroy(work2)

        Sigma_p(i,i)%val = Sigma_p(i,i)%val + work1%val

        call destroy(work1, Mq_mat)

      end do

    end do

    ! throw away all non diagonal parts
    !do i = 1, nbl
    !   i1 = Sigma_p(i,i)%nrow
    !   do m = 1, i1
    !      Sigma_p(i,i)%val(1:m-1,m) = (0.0_dp,0.0_dp)
    !      Sigma_p(i,i)%val(m+1:i1,m) = (0.0_dp,0.0_dp)
    !   end do
    !end do
    !

    do i = 1, nbl
      call write_blkmat(Sigma_p(i,i),negf%scratch_path,'Sigma_ph_p_',i,i,iE)
      call destroy(Sigma_p(i,i))
      call destroy(G_p_interN(i,i))
      call destroy(G_p_interP(i,i))
    end do

    deallocate(Sigma_p,G_p_interN, G_p_interP)

  end subroutine Sigma_ph_p

  !--------------------------------------------------------------------------------
  subroutine Sigma_ph_r(negf,Epnt)

    type(Tnegf) :: negf
    real(dp), dimension(:),allocatable :: Epnt

    !Local variables
    type(z_DNS) :: work1,work2
    type(z_DNS) :: Mq_mat
    type(z_DNS), dimension(:,:), allocatable :: Sigma_r, G_r_interP, G_r_interN
    type(z_DNS), dimension(:,:), allocatable :: G_n_interP, G_n_interN

    real(dp), dimension(:), pointer :: Wq
    real(dp), dimension(:), pointer :: Nq
    real(dp), dimension(:), pointer :: Mq
    logical, dimension(:), pointer :: selmodes
    integer :: i, m, iE, i1, i2
    integer :: nummodes, numselmodes, nbl
    real(dp) :: E1, E2, En

    nbl = negf%str%num_PLs


    selmodes => negf%elph%selmodes
    Mq => negf%elph%Mq
    Wq => negf%elph%Wq
    Nq => negf%elph%Nq
    nummodes = negf%elph%nummodes
    numselmodes = negf%elph%numselmodes

    iE=negf%iE

    allocate(G_r_interP(nbl,nbl))
    allocate(G_r_interN(nbl,nbl))
    allocate(G_n_interP(nbl,nbl))
    allocate(G_n_interN(nbl,nbl))
    allocate(Sigma_r(nbl,nbl))

    associate(indblk=>negf%str%mat_PL_start)
    do i = 1, nbl
      m = indblk(i+1)-indblk(i)
      call create(Sigma_r(i,i), m, m)
      Sigma_r(i,i)%val = (0.0_dp, 0.0_dp)
      call create(G_r_interP(i,i), m, m)
      call create(G_r_interN(i,i), m, m)
      call create(G_n_interP(i,i), m, m)
      call create(G_n_interN(i,i), m, m)
      !call create(G_r(i,i), m, m)
      !call read_blkmat(G_r(i,i),negf%scratch_path,'G_r_',i,i,iE)
    end do
    end associate
    do m = 1 , nummodes

      if (.not.selmodes(m)) cycle

      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(negf,Wq(m),Epnt,i1,i2,E1,E2)
      En = real(negf%Epnt)+Wq(m)
      call interpolation(i1, i2, E1, E2, En, negf%scratch_path, &
          'G_n_', G_n_interP)
      call interpolation(i1, i2, E1, E2, En, negf%scratch_path, &
          'G_r_', G_r_interP)

      ! INTERPOLATION OF GREEN'S FUNCTIONS BETWEEN GRID POINTS
      call search_points(negf,-Wq(m),Epnt,i1,i2,E1,E2)
      En = real(negf%Epnt)-Wq(m)
      call interpolation(i1, i2, E1, E2, En, negf%scratch_path, &
          'G_n_', G_n_interN)
      call interpolation(i1, i2, E1, E2, En, negf%scratch_path, &
          'G_r_', G_r_interN)

      do i = 1, nbl

        i1 = Sigma_r(i,i)%nrow

        call create(work1,i1,i1)

        work1%val = (0.0_dp,0.0_dp)

        if (negf%elph%selfene_gr) then
          !print*,'SelfEneR_Gr'
          ! Via I: should be exact for E >> Ef_max
          work1%val = (Nq(m)+1.0_dp)*G_r_interN(i,i)%val + Nq(m)*G_r_interP(i,i)%val
          ! Via II: should be exact for E << Ef_min
          !work1%val = (Nq(m)+1.0_dp)*G_r_interP(i,i)%val + Nq(m)*G_r_interN(i,i)%val
          ! Via III: should work as a compromise
          !work1%val = Nq(m)*(G_r_interP(i,i)%val + G_r_interN(i,i)%val) + G_r(i,i)%val
        endif

        if (negf%elph%selfene_gless) then
          !print*,'SelfEneR_G<'
          ! should be 1/2 [G<- - G<+] == i/2 [Gn- - Gn+]
          work1%val = work1%val + &
              (0.0_dp,0.5_dp) * (G_n_interN(i,i)%val - G_n_interP(i,i)%val)
        endif

        if (negf%elph%selfene_gless .or. negf%elph%selfene_gr) then

          call create_id(Mq_mat,i1,Mq(m))

          call prealloc_mult(Mq_mat, work1, work2)

          call destroy(work1)

          call prealloc_mult(work2, Mq_mat, work1)

          call destroy(work2)

          Sigma_r(i,i)%val = Sigma_r(i,i)%val + work1%val

          call destroy(work1, Mq_mat)

        else

          Sigma_r(i,i)%val = (0.0_dp,0.0_dp)
          call destroy(work1)

        endif

      end do

    end do

    ! throw away all non diagonal parts
    !do i = 1, nbl
    !    i1 = Sigma_r(i,i)%nrow
    !    do m = 1, i1
    !       Sigma_r(i,i)%val(1:m-1,m) = (0.0_dp,0.0_dp)
    !       Sigma_r(i,i)%val(m+1:i1,m) = (0.0_dp,0.0_dp)
    !    end do
    ! end do
    !

    do i = 1, nbl
      !print*
      !print*,'(sigma_r) Sigma_ph_r',maxval(abs(Sigma_r(i,i)%val)), iE
      call write_blkmat(Sigma_r(i,i),negf%scratch_path,'Sigma_ph_r_',i,i,iE)
      call destroy(Sigma_r(i,i))
      call destroy(G_r_interP(i,i))
      call destroy(G_r_interN(i,i))
      call destroy(G_n_interP(i,i))
      call destroy(G_n_interN(i,i))
      !call destroy(G_r(i,i))
    end do

    deallocate(Sigma_r,G_r_interP,G_r_interN,G_n_interP,G_n_interN)
    !deallocate(G_r)

  end subroutine Sigma_ph_r

  subroutine Sigma_ph_r_z(negf,z)

    type(Tnegf) :: negf
    complex(dp), intent(in) :: z

    !Local variables
    type(z_DNS) :: work1,work2
    type(z_DNS) :: Mq_mat
    type(z_DNS), dimension(:,:), allocatable :: Sigma_r, G_r
    real(dp), dimension(:), pointer :: Wq
    real(dp), dimension(:), pointer :: Nq
    real(dp), dimension(:), pointer :: Mq
    logical, dimension(:), pointer :: selmodes
    integer :: i, m, iE,i1
    integer :: nummodes, numselmodes, nbl

    nbl = negf%str%num_PLs

    selmodes => negf%elph%selmodes
    Mq => negf%elph%Mq
    Wq => negf%elph%Wq
    Nq => negf%elph%Nq
    nummodes = negf%elph%nummodes
    numselmodes = negf%elph%numselmodes

    iE=negf%iE

    allocate(Sigma_r(nbl,nbl))
    allocate(G_r(nbl,nbl))
    
    associate(indblk => negf%str%mat_PL_start)
    do i = 1, nbl
      m = indblk(i+1)-indblk(i)
      call create(Sigma_r(i,i), m, m)
      Sigma_r(i,i)%val = (0.0_dp, 0.0_dp)
      call create(G_r(i,i), m, m)
    end do
    end associate

    do m = 1 , nummodes

      if (.not.selmodes(m)) cycle

      do i = 1, nbl

        call read_blkmat(G_r(i,i),negf%scratch_path,'G_r_',i,i,iE)

        i1 = Sigma_r(i,i)%nrow

        if (negf%elph%selfene_gr) then

          call create(work1,i1,i1)

          work1%val = (0.0_dp,0.0_dp)

          work1%val = (2.0_dp*Nq(m)+1.0_dp)*G_r(i,i)%val

          call create_id(Mq_mat,i1,Mq(m))

          call prealloc_mult(Mq_mat, work1, work2)

          call destroy(work1)

          call prealloc_mult(work2, Mq_mat, work1)

          call destroy(work2)

          Sigma_r(i,i)%val = Sigma_r(i,i)%val +  work1%val

          call destroy(work1, Mq_mat)

        else
          Sigma_r(i,i)%val = (0.0_dp,0.0_dp)
        endif

      end do

    end do

    do i = 1, nbl
      call write_blkmat(Sigma_r(i,i),negf%scratch_path,'Sigma_ph_r_',i,i,iE)
      call destroy(Sigma_r(i,i))
      call destroy(G_r(i,i))
    end do

  end subroutine Sigma_ph_r_z
  ! ----------------------------------------------------------
  subroutine check_Gl_Gr(negf)
    type(Tnegf) :: negf
    real(dp), dimension(:),allocatable :: Epnt
    integer :: ioffset

    type(z_DNS), dimension(:,:), allocatable :: G_r, G_p, G_n
    type(z_DNS) :: A, T
    integer :: nbl, n, sizebl, i_start, i_stop, psize, maxpos(2)
    real(dp) :: Wmax, maxdev, tmp, maxG

    nbl = negf%str%num_PLs

    allocate(G_n(nbl,nbl))
    allocate(G_p(nbl,nbl))
    allocate(G_r(nbl,nbl))

    maxdev = 0.0_dp
    maxG = 0.0_dp
    psize = 0

    do n = 1, nbl
      sizebl = negf%str%mat_PL_start(n+1)-negf%str%mat_PL_start(n)
      call create(G_r(n,n),sizebl,sizebl)
      call create(G_n(n,n),sizebl,sizebl)
      call create(G_p(n,n),sizebl,sizebl)
      call read_blkmat(G_r(n,n), negf%scratch_path, 'G_r_',n,n, negf%iE)
      call read_blkmat(G_n(n,n),negf%scratch_path,'G_n_',n,n,negf%iE)
      call read_blkmat(G_p(n,n),negf%scratch_path,'G_p_',n,n,negf%iE)

      call zspectral(G_r(n,n),G_r(n,n),0,A)

      call create(T,sizebl,sizebl)


      T%val = G_n(n,n)%val + G_p(n,n)%val - A%val

      tmp = maxval(abs(A%val))
      if (tmp .gt. maxG) maxG=tmp

      tmp = maxval(abs(T%val))/maxval(abs(A%val))
      if (tmp .gt. maxdev) then
        maxdev = tmp
        maxpos = maxloc(abs(T%val)) + psize
      endif

      psize = psize + sizebl

      call destroy(G_r(n,n))
      call destroy(G_n(n,n))
      call destroy(G_p(n,n))
      call destroy(A)
      call destroy(T)

    end do

    !print*,'CHECK Gn+Gp=Gr-Ga',negf%iE, maxG, maxdev

    deallocate(G_n, G_p, G_r)

  end subroutine check_Gl_Gr

  ! ----------------------------------------------------------
  subroutine check_sigma_ph_r(negf)
    type(Tnegf) :: negf

    type(z_DNS), dimension(:,:), allocatable :: Sigma_r, Sigma_p, Sigma_n
    type(z_DNS) :: Gam, T
    integer :: nbl, n, sizebl, psize, maxpos(2)
    real(dp) :: maxdev, tmp, maxG

    nbl = negf%str%num_PLs

    allocate(Sigma_n(nbl,nbl))
    allocate(Sigma_p(nbl,nbl))
    allocate(Sigma_r(nbl,nbl))

    maxdev = 0.0_dp
    maxG = 0.0_dp
    psize = 0

    do n = 1, nbl
      sizebl = negf%str%mat_PL_start(n+1)-negf%str%mat_PL_start(n)
      call create(Sigma_r(n,n),sizebl,sizebl)
      call create(Sigma_n(n,n),sizebl,sizebl)
      call create(Sigma_p(n,n),sizebl,sizebl)
      call read_blkmat(Sigma_r(n,n), negf%scratch_path, 'Sigma_ph_r_',n,n, negf%iE)
      call read_blkmat(Sigma_n(n,n),negf%scratch_path,'Sigma_ph_n_',n,n,negf%iE)
      call read_blkmat(Sigma_p(n,n),negf%scratch_path,'Sigma_ph_p_',n,n,negf%iE)

      call zspectral(Sigma_r(n,n),Sigma_r(n,n),0,Gam)

      call create(T,sizebl,sizebl)

      !print*,'CHECK Sigma_ph_n+Sigma_ph_p',maxval(abs(Sigma_n(n,n)%val+Sigma_p(n,n)%val))

      T%val = Sigma_n(n,n)%val + Sigma_p(n,n)%val - Gam%val

      tmp = maxval(abs(Gam%val))
      if (tmp .gt. maxG) maxG=tmp

      tmp = maxval(abs(T%val))/maxval(abs(Gam%val))
      if (tmp .gt. maxdev) then
        maxdev = tmp
        maxpos = maxloc(abs(T%val)) + psize
      endif

      psize = psize + sizebl

      call destroy(Sigma_r(n,n))
      call destroy(Sigma_n(n,n))
      call destroy(Sigma_p(n,n))
      call destroy(Gam)
      call destroy(T)

    end do

    !print*,'CHECK Sigma_ph_r',negf%iE, maxG, maxdev

    deallocate(Sigma_n, Sigma_p, Sigma_r)


  end subroutine check_sigma_ph_r


  ! ----------------------------------------------------------
  ! Search points for interpolations.
  ! Wq has a sign (+/- Wq)
  !
  subroutine search_points(negf, Wq, Epnt, i1, i2, E1, E2)
    use energy_mesh, only : elem
    type(TNegf) :: negf
    real(dp) :: Wq
    real(dp), dimension(:), allocatable :: Epnt
    integer, intent(out) :: i1,i2
    real(dp), intent(out) :: E1, E2

    integer :: iE, iel, istart, iend, ip
    Type(elem), pointer :: pel
    real(dp) :: En

    if (Wq.eq.0) then
      i1 = negf%iE
      i2 = negf%iE
      E1 = real(negf%Epnt)
      E2 = real(negf%Epnt)
      return
    endif

    if (allocated(Epnt)) then
      ! Remove offset such that search can work on Epnt(1..N)
      iE = negf%iE - negf%Np_n(1) - negf%Np_n(2) - negf%n_poles
      En = real(negf%Epnt) + Wq !Wq carry the right sign
      !print*
      !print*,'iE', iE, real(negf%Epnt) + Wq

      if (sign(1.0_dp,Wq) .gt. 0) then
        i2 = iE + 1
        if (i2.gt.size(Epnt)) then
          i2 = size(Epnt)
        else
          do while (Epnt(i2) .lt. En)
            i2 = i2 + 1
          end do
        endif
        i1 = i2 - 1
      else
        i1 = iE - 1
        if (i1.lt.1) then
          i1 = 1
        else
          do while (Epnt(i1) .gt. En)
            i1 = i1 - 1
          end do
        endif
        i2 = i1 + 1
      endif
      E1 = Epnt(i1)
      E2 = Epnt(i2)
      ! add back the offset to the point
      i1 = i1 + negf%Np_n(1) + negf%Np_n(2) + negf%n_poles
      i2 = i2 + negf%Np_n(1) + negf%Np_n(2) + negf%n_poles

    else

      !if (.not.allocated(negf%emesh%pactive)) STOP 'emesh not initialized'

      En = real(negf%Epnt) + Wq

      if (sign(1.0_dp,Wq) .gt. 0) then
        istart = negf%emesh%iactive
        iend = negf%emesh%maxind

        elloop1: do iel = istart, iend
          pel => negf%emesh%pactive(iel)%pelem
          do ip = 1, 3
            if (pel%pnt(ip) .gt. En) then
              exit elloop1
            endif
          end do
        end do elloop1
        i1 = pel%map(ip-1)
        i2 = pel%map(ip)
        E1 = pel%pnt(ip-1)
        E2 = pel%pnt(ip)
      else
        istart = negf%emesh%iactive
        iend = 1

        elloop2: do iel = istart, iend, -1
          pel => negf%emesh%pactive(iel)%pelem
          do ip = 3, 1, -1
            if (pel%pnt(ip) .lt. En) then
              exit elloop2
            endif
          end do
        end do elloop2
        i1 = pel%map(ip)
        i2 = pel%map(ip+1)
        E1 = pel%pnt(ip)
        E2 = pel%pnt(ip+1)
      end if
    end if
    !print*
    !print*,E1,En,E2
    !print*,'interpolate between:',i1,i2

  end subroutine search_points


  subroutine interpolation(i1,i2, E1, E2, E, path, name, G_interp)
    integer, intent(in) :: i1, i2
    real(dp) :: E1, E2, E
    type(z_DNS), dimension(:,:) :: G_interp
    CHARACTER(*) :: path
    CHARACTER(*) :: name

    !local variables
    type(z_DNS) :: work1,work2
    integer :: i

    do i = 1, size(G_interp,1)

      call create(work1,G_interp(i,i)%nrow,G_interp(i,i)%ncol)
      work1%val = (0.0_dp,0.0_dp)
      call read_blkmat(work1, path, name, i, i, i1)

      if (E1.ne.E2 .and. i1.ne.i2) then

        call create(work2,G_interp(i,i)%nrow,G_interp(i,i)%ncol)
        work2%val = (0.0_dp,0.0_dp)
        call read_blkmat(work2, path, name, i, i, i2)

        G_interp(i,i)%val = ((E-E1)*work2%val + (E2-E)*work1%val)/(E2-E1)

        call destroy(work2)

      else

        G_interp(i,i)%val = work1%val

      endif

      call destroy(work1)

    end do


  end subroutine interpolation


  !*********************************************************************
  !ADD Hilbert-transform part to Sigma_ph_r
  !Need to set back FFT transforms
  !********************************************************************
  subroutine complete_sigma_ph_r(negf, Epnt, ioffset)

    type(Tnegf) :: negf
    real(dp), dimension(:) :: Epnt
    integer :: ioffset

    ! Locals
    integer :: i,j,n,m, iE, bl, sizebl,i_start,i_stop, ierr, nummodes
    real(dp), dimension(:), pointer :: Wq
    real(dp), dimension(:), pointer :: Mq
    real(dp) :: Wmax, dE, tmp

    character(4) :: ofblki
    complex(dp), dimension(:), allocatable :: temp1
    type(z_DNS), dimension(:), allocatable :: Gn_E, Sigma_r_E
    !type(z_DNS) :: work1, work2, work3, Mq_mat
    logical, dimension(:), pointer :: selmodes

    Mq => negf%elph%Mq
    Wq => negf%elph%Wq
    selmodes => negf%elph%selmodes

    n = size(Epnt)
    nummodes = size(Wq)

    Wmax = maxval(Wq)
    dE = Epnt(2)-Epnt(1)
    ! ENERGY-INTERVAL FOR SELF-ENERGIES:
    i_start = int(aint(Wmax/dE) + 1)
    i_stop =  (n - i_start)
    i_start = i_start + 1

    ! PROBLEM:
    ! The Hilb-Transf should be resticted on the sub interval
    ! Emin+m*Wmax ... Emax-m*Wmax
    ! However the number of points should be always 2^p
    ! This condition is difficult to fulfill.
    ! Possible solution: interpolation between grids.

    ! CREATE THE ARRAY OF MATRICES (one for each E-point)
    allocate(Gn_E(n),stat=ierr)
    allocate(Sigma_r_E(n),stat=ierr)
    if (ierr.ne.0) STOP 'ERROR in allocation of Gn_E'
    call log_allocate(temp1,n)

    print*
    print*,'HILBERT TRANSFORM memory:',sizebl*sizebl*n*16
    print*,'HILBERT TRANSFORM interval:',i_start,i_stop

    ! LOOP ON blocks
    do bl = 1, negf%str%num_PLs

      sizebl =  negf%str%mat_PL_start(bl+1) - negf%str%mat_PL_start(bl)
      if (bl.le.9999) write(ofblki,'(i4.4)') bl
      if (bl.gt.9999) stop 'ERROR: too many blks (> 9999)'

      ! LOAD ALL G< FROM FILES and store in Gn_E(iE)
      do iE = 1, n
        call create(Gn_E(iE),sizebl,sizebl)
        call read_blkmat(Gn_E(iE), negf%scratch_path, 'G_n_', bl, bl,iE+ioffset)
      end do

      ! LOAD ALL Sigma_r in the right interval
      do iE = i_start, i_stop
        call create(Sigma_r_E(iE),sizebl,sizebl)
        call read_blkmat(Sigma_r_E(iE), negf%scratch_path, 'Sigma_ph_r_', bl, bl, iE+ioffset)
      end do

      do m = 1, nummodes
        if (.not.selmodes(m)) cycle

        ! COMPUTE   Mq G< Mq   (assume diagonal now)
        do iE = 1, n
          !call create_id(Mq_mat,sizebl,Mq(m))
          !call prealloc_mult(Mq_mat, Gn_E(iE), work1)
          !call prealloc_mult(work1, Mq_mat, work2)
          !call destroy(work1, Mq_mat)
          !Gn_E(iE)%val = work2%val
          !call destroy(work2)
          tmp = Mq(m)*Mq(m)
          Gn_E(iE)%val = tmp*Gn_E(iE)%val
        end do

        ! PERFORM Hilbert in Energy for each (i,j)
        do j = 1, sizebl
          do i = 1, sizebl

            !          k = (j-1)*sizebl+i
            !          do iE = 1, n
            !             if (iE.le.99999) write(ofpnt,'(i5.5)') iE+ioffset
            !             filename = 'G_n_'//trim(ofblki)//'_'//trim(ofblki)//'_'//trim(ofpnt)//'.dat'
            !             open(100,file=trim(negf%scratch_path)//trim(filename), access='DIRECT', recl=4)
            !             READ(100, rec = k)  temp1(iE)
            !             close(100)
            !          end do

            ! SETUP a vector out of all G<_ij(E)
            ! Here we could perform an efficient interpolation on a regular grid 2^p
            do iE = 1, n
              temp1(iE) = Gn_E(iE)%val(i,j)
            end do

            Wq(m) = Wq(m)*2.0_dp*pi/((n-1)*dE)

            !call Hilbert_shift(temp1, Wq(m))

            Wq(m) = Wq(m)*((n-1)*dE)/(2.0_dp*pi)

            ! UPDATE the self-energies with the Hilbert part.
            do iE = i_start, i_stop
              ! Should be  -i/2 H[ G<+ - G<-  ] = 1/2 H[ Gn+ - Gn- ]
              Sigma_r_E(iE)%val(i,j) =  Sigma_r_E(iE)%val(i,j) + (0.5_dp, 0.0)* temp1(iE)
            end do


            !          do iE = i_start+1, i_stop
            !             if (iE.le.99999) write(ofpnt,'(i5.5)') iE+ioffset
            !             filename = 'Sigma_ph_r_'//trim(ofblki)//'_'//trim(ofblki)//'_'//trim(ofpnt)//'.dat'
            !             open(100,file=trim(negf%scratch_path)//trim(filename), access='DIRECT', recl=4)
            !             READ (100, rec = k) temp2
            !             temp2 = temp2 - (0.0_dp, 0.5_dp)* temp1(iE)
            !             WRITE (100, rec = k) temp2
            !             close(100)
            !          end do

          end do ! Loop on block size
        end do ! Loop on block size

      end do !Loop on modes

      do iE = i_start, i_stop
        call write_blkmat(Sigma_r_E(iE), negf%scratch_path, 'Sigma_ph_r_', bl, bl, iE+ioffset)
      end do

      do iE = 1, n
        call destroy(Gn_E(iE))
      end do
      do iE = i_start, i_stop
        call destroy(Sigma_r_E(iE))
      end do

    end do !Loop on blocks

    call log_deallocate(temp1)
    deallocate(Gn_E)
    deallocate(Sigma_r_E)

  end subroutine complete_sigma_ph_r

  subroutine create_scratch(nbl, npoints)
    integer :: nbl, npoints

    integer :: i,j,k,err

    call destroy_scratch(nbl, npoints)

    ALLOCATE(GGn(nbl,nbl,npoints),stat=err)
    ALLOCATE(GGr(nbl,nbl,npoints),stat=err)
    ALLOCATE(GGp(nbl,nbl,npoints),stat=err)
    ALLOCATE(Sgn(nbl,nbl,npoints),stat=err)
    ALLOCATE(Sgr(nbl,nbl,npoints),stat=err)
    ALLOCATE(Sgp(nbl,nbl,npoints),stat=err)

    ! Initialize everything to 0
    do k=1,npoints
      do j=1,nbl
        do i=1,nbl
          GGn(i,j,k)%nrow=0
          GGn(i,j,k)%ncol=0
          GGr(i,j,k)%nrow=0
          GGr(i,j,k)%ncol=0
          GGp(i,j,k)%nrow=0
          GGp(i,j,k)%ncol=0
          Sgn(i,j,k)%nrow=0
          Sgn(i,j,k)%ncol=0
          Sgr(i,j,k)%nrow=0
          Sgr(i,j,k)%ncol=0
          Sgp(i,j,k)%nrow=0
          Sgp(i,j,k)%ncol=0
        end do
      end do
    end do

    if (err.ne.0) then
      STOP 'ERROR: Cannot allocate GG'
    endif
    print*,'Created memory scratch',nbl,'x',nbl,'x',npoints

  end subroutine create_scratch

  subroutine destroy_scratch(nbl, npoints)
    integer :: nbl, npoints

    integer :: i, j, iE,  err

    err = 0


    do i = 1, nbl
      do j = 1, nbl
        do iE = 1, npoints

          if (allocated(GGn)) then
            if (allocated(GGn(i,j,iE)%val)) call destroy(GGn(i,j,iE))
          endif
          if (allocated(GGr)) then
            if (allocated(GGr(i,j,iE)%val)) call destroy(GGr(i,j,iE))
          endif
          if (allocated(GGp)) then
            if (allocated(GGp(i,j,iE)%val)) call destroy(GGp(i,j,iE))
          endif
          if (allocated(Sgn)) then
            if (allocated(Sgn(i,j,iE)%val)) call destroy(Sgn(i,j,iE))
          endif
          if (allocated(Sgr)) then
            if (allocated(Sgr(i,j,iE)%val)) call destroy(Sgr(i,j,iE))
          endif
          if (allocated(Sgp)) then
            if (allocated(Sgp(i,j,iE)%val)) call destroy(Sgp(i,j,iE))
          endif

        end do
      end do
    end do

    if (allocated(GGn)) DEALLOCATE(GGn,stat=err)
    if (allocated(GGr)) DEALLOCATE(GGr,stat=err)
    if (allocated(GGp)) DEALLOCATE(GGp,stat=err)
    if (allocated(Sgn)) DEALLOCATE(Sgn,stat=err)
    if (allocated(Sgr)) DEALLOCATE(Sgr,stat=err)
    if (allocated(Sgp)) DEALLOCATE(Sgp,stat=err)

    if (err.ne.0) then
      STOP 'ERROR: Cannot deallocate GG'
    endif

  end subroutine destroy_scratch

  ! READ Matrices
  subroutine read_blkmat(Matrix, path, name, i, j, iE)
    type(z_DNS), intent(inout) :: Matrix
    CHARACTER(*), intent(in) :: path
    CHARACTER(*), intent(in) :: name
    integer, intent(in) :: i, j, iE

    CHARACTER(64) :: filename
    character(4) :: ofblki, ofblkj
    character(10) :: ofpnt
    logical :: lex
    integer :: i1,i2
    complex(dp) :: mat_el

    if (memory) then
      select case(trim(name))
      case('G_n_')
        Matrix%val = GGn(i,j,iE)%val
      case('G_r_')
        Matrix%val = GGr(i,j,iE)%val
      case('G_p_')
        Matrix%val = GGp(i,j,iE)%val
      case('Sigma_ph_r_')
        Matrix%val = Sgr(i,j,iE)%val
      case('Sigma_ph_n_')
        Matrix%val = Sgn(i,j,iE)%val
      case('Sigma_ph_p_')
        Matrix%val = Sgp(i,j,iE)%val
      case default
        stop 'internal error: read_blkmat does not correspond'
      end select
      return
    endif

    Matrix%val = (0.0_dp,0.0_dp)

    if (i.le.9999) write(ofblki,'(i4.4)') i
    if (i.gt.9999) stop 'ERROR: too many blks (> 9999)'
    if (j.le.9999) write(ofblkj,'(i4.4)') j
    if (j.gt.9999) stop 'ERROR: too many blks (> 9999)'

    if (iE.le.99999) write(ofpnt,'(i5.5)') iE

    filename = trim(name)//trim(ofblki)//'_'//trim(ofblkj)//'_'//trim(ofpnt)//'.dat'

    inquire(file=trim(path)//trim(filename),EXIST=lex)
    if (.not.lex) then
      RETURN
      !WRITE(*,*) 'ERROR: FILE '//trim(filename)//' doES NOT EXIST'
      !STOP
    endif

    open(9091,file=trim(path)//trim(filename), access='STREAM')

    call inmat_c(9091,.false.,Matrix%val,Matrix%nrow,Matrix%ncol)

    !open(9001,file=trim(path)//trim(filename), access='DIRECT', recl=4)
    !call direct_in_c(9001,Matrix%val,Matrix%nrow)

    close(9091)

  end subroutine read_blkmat

  ! WRITE Matrices
  subroutine write_blkmat(Matrix, path, name, i, j, iE)
    type(z_DNS), intent(in) :: Matrix
    CHARACTER(*), intent(in) :: path
    CHARACTER(*), intent(in) :: name
    integer, intent(in) :: i, j, iE

    CHARACTER(64) :: filename
    character(4) :: ofblki, ofblkj
    character(10) :: ofpnt
    integer :: m,n

    if (memory) then

      select case(trim(name))
      case('G_n_')
        if (.not.allocated(GGn(i,j,iE)%val)) then
          call create(GGn(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        GGn(i,j,iE)%val = Matrix%val
      case('G_r_')
        if (.not.allocated(GGr(i,j,iE)%val)) then
          call create(GGr(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        GGr(i,j,iE)%val = Matrix%val
      case('G_p_')
        if (.not.allocated(GGp(i,j,iE)%val)) then
          call create(GGp(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        GGp(i,j,iE)%val = Matrix%val
      case('Sigma_ph_r_')
        if (.not.allocated(Sgr(i,j,iE)%val)) then
          call create(Sgr(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        Sgr(i,j,iE)%val = Matrix%val
      case('Sigma_ph_n_')
        if (.not.allocated(Sgn(i,j,iE)%val)) then
          call create(Sgn(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        Sgn(i,j,iE)%val = Matrix%val
      case('Sigma_ph_p_')
        if (.not.allocated(Sgp(i,j,iE)%val)) then
          call create(Sgp(i,j,iE),Matrix%nrow,Matrix%ncol)
        endif
        Sgp(i,j,iE)%val = Matrix%val
      case default
        stop 'internal error: write_blkmat does not correspond'
      end select
      return
    endif

    if (i.le.9999) write(ofblki,'(i4.4)') i
    if (i.gt.9999) stop 'ERROR: too many blks (> 9999)'
    if (j.le.9999) write(ofblkj,'(i4.4)') j
    if (j.gt.9999) stop 'ERROR: too many blks (> 9999)'

    if (iE.le.99999) write(ofpnt,'(i5.5)') iE

    filename = trim(name)//trim(ofblki)//'_'//trim(ofblkj)//'_'//trim(ofpnt)//'.dat'

    open(9001,file=trim(path)//trim(filename), access='STREAM', status='REPLACE')

    call outmat_c(9001,.false.,Matrix%val,Matrix%nrow,Matrix%ncol) !,1.0d-36)

    !open(9001,file=trim(path)//trim(filename), status='REPLACE', access='DIRECT', recl=4)
    !call direct_out_c(9001,Matrix%val,Matrix%nrow)

    close(9001)

  end subroutine write_blkmat


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
    type(z_DNS), dimension(:) :: Tlc,Tcl,gsurfR
    logical :: lower
    type(Tstruct_info), intent(in) :: struct
    type(z_CSR) :: Aout


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
      call prealloc_mult(Gr(cb,cb),Tlc(i),(-1.0_dp, 0.0_dp),work1)
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

        call prealloc_mult(gsurfR(i),Tcl(i),(-1.0_dp, 0.0_dp), work1)
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



  !****************************************************************************
  !
  !  Calculate Gn contributions for all contacts except collector, in the
  !  outer region where contacts-device overlap is non-zero - writing on
  !  memory
  !
  !****************************************************************************

  subroutine calculate_Gn_outer(Tlc,gsurfR,SelfEneR,struct,frm,ref,lower,Glout)

    !****************************************************************************
    !Input:
    !Tlc: sparse matrices array containing contacts-device interaction blocks (ES-H)
    !gsurfR: sparse matrices array containing contacts surface green
    !SelfEneR: sparse matrices array containing contacts Self Energy
    !frm: array containing Fermi distribution values for all contacts
    !ref: reference contact
    !
    !global variables needed: nbl, indblk(nbl+1), cindblk(ncont), ncont,
    !cblk(ncont), Gr(:,:), diagonal, subdiagonal, overdiagonal and
    !in colums Gr(:,cb) where cb are layers interacting with all contacts
    !but collector
    !
    !Output:
    !Glout: sparse matrix containing G_n  contributions in the region
    !       corresponding to non-zero overlap
    !
    !****************************************************************************

    implicit none

    !In/Out
    type(z_DNS), dimension(:) :: Tlc, gsurfR, SelfEneR
    real(dp), dimension(:) :: frm
    type(Tstruct_info), intent(in) :: struct
    integer :: ref
    logical :: lower
    type(z_CSR) :: Glout

    !Work
    type(z_DNS) :: Gam, gsurfA, Ga, work1, work2, work3, Glsub
    type(z_CSR) :: GlCSR, TCSR
    integer :: j,k,cb,cbj,i1,j1,nrow_tot
    integer :: ncont, nbl
    complex(dp) :: frmdiff


    ncont = struct%num_conts
    nbl = struct%num_PLs
    nrow_tot = struct%total_dim

    !Allocazione della Glout
    !Righe totali del conduttore effettivo nrow_tot
    !nrow_tot=indblk(nbl+1)-1
    !do i=1,ncont
    !   nrow_tot=nrow_tot+ncdim(i) !gsurfR(i)%nrow
    !end do
    if (.not.allocated(Glout%nzval)) THEN
      call create(Glout,nrow_tot,nrow_tot,0)
      Glout%rowpnt(:)=1
    endif
    !***
    !Iterazione su tutti i contatti "k"
    !***
    do k=1,ncont

      !Esegue le operazioni relative al contatto solo se e` valida la condizione
      !sulle distribuzioni di Fermi e se non si tratta del contatto iniettante (ref)
      if ((ABS(frm(k)-frm(ref)).GT.EPS).AND.(k.NE.ref)) THEN

        !Calcolo della Gamma corrispondente
        cb=struct%cblk(k)
        !nota:cb indica l'indice di blocco corrispondente al contatto j-esimo

        call zspectral(SelfEneR(k),SelfEneR(k),0,Gam)

        !***
        !Calcolo del contributo sulla proria regione
        !***
        frmdiff = cmplx(frm(ref)-frm(k),0.0_dp,dp)

        call zspectral(gsurfR(k),gsurfR(k),0,work1)

        !print*, 'work2=Tlc*work1=Tlc*j(gsurfR-gsurfA)'
        call prealloc_mult(Tlc(k),work1,work2)
        call destroy(work1)

        !print *, 'work1=-Gr(cb,cb)*work2=-Gr(cb,cb)*Tlc*j(gsurfR-gsurfA)'
        call prealloc_mult(Gr(cb,cb),work2,frmdiff,work1)
        call destroy(work2)

        !print*,'work2=Tlc*gsurfA'
        call zdagger(gsurfR(k),gsurfA)
        call prealloc_mult(Tlc(k),gsurfA,work2)
        call destroy(gsurfA)

        !print*,'work3=Ga*work2=Ga*Tlc*gsurfA'
        call zdagger(Gr(cb,cb),Ga)
        call prealloc_mult(Ga,work2,work3)

        call destroy(Ga)
        call destroy(work2)

        !print*,'work2=Gam*work3=Gam*Ga*Tlc*gsurfA'
        call prealloc_mult(Gam,work3,work2)
        call destroy(work3)

        !print *,'work3=-Gr*work2=-Gr*Gam*Ga*Tlc*gsurfA'
        call prealloc_mult(Gr(cb,cb),work2,frmdiff,work3)
        call destroy(work2)

        !Contributo totale sulla propria regione
        call prealloc_sum(work3,work1,Glsub)
        call destroy(work1)
        call destroy(work3)

        call mask(Glsub,Tlc(k))
        i1=nzdrop(Glsub,EPS)

        if (i1.gt.0) THEN
          call create(GlCSR,Glsub%nrow,Glsub%ncol,i1)
          call dns2csr(Glsub,GlCSR)

          !Concatenazione di Glsub nella matrice globale Glout
          i1=struct%mat_PL_start(cb)
          j1=struct%mat_B_start(k)-struct%central_dim+struct%mat_PL_start(nbl+1)-1
          call concat(Glout,GlCSR,i1,j1)

          ! compute lower outer part using (iG<)+ = iG<
          if (lower) THEN
            call zdagger(GlCSR,TCSR)
            call concat(Glout,TCSR,j1,i1)
            call destroy(TCSR)
          endif

          call destroy(GlCSR)
        end if

        call destroy(Glsub)
        !***
        !Ciclo per il calcolo dei contributi sulle regioni degli altri contatti
        !***

        do j=1,ncont

          cbj=struct%cblk(j)
          !Esegue le operazioni del ciclo solo se il j.ne.k o se
          !il blocco colonna di Gr e` non nullo (altrimenti il contributo e` nullo)

          if ((j.NE.k).AND.(Gr(cbj,cb)%nrow.NE.0 .AND. (Gr(cbj,cb)%ncol.NE.0))) THEN

            !print*,'work1=Tlc*gsurfA'
            call zdagger(gsurfR(j),gsurfA)
            call prealloc_mult(Tlc(j),gsurfA,work1)
            call destroy(gsurfA)

            !print*,'work2=Ga*work1=Ga*Tlc*gsurfA'
            call zdagger(Gr(cbj,cb),Ga)
            call prealloc_mult(Ga,work1,work2)

            call destroy(Ga)
            call destroy(work1)

            !print*,'work1=Gam*work2=Gam*Ga*Tlc*gsurfA'
            call prealloc_mult(Gam,work2,work1)
            call destroy(work2)

            !print*,'Glsub=-Gr*work1=-Gr*Gam*Ga*Tlc*gsurfA'
            call prealloc_mult(Gr(cbj,cb),work1,frmdiff,Glsub)
            call destroy(work1)

            call mask(Glsub,Tlc(j))
            i1=nzdrop(Glsub,EPS)

            if (i1.gt.0) THEN
              call create(GlCSR,Glsub%nrow,Glsub%ncol,i1)
              call dns2csr(Glsub,GlCSR)

              !Concatenazione di Glsub nella posizione corrispondente al contatto "j"
              i1=struct%mat_PL_start(cbj)
              j1=struct%mat_B_start(j)-struct%central_dim+struct%mat_PL_start(nbl+1)-1

              call concat(Glout,GlCSR,i1,j1)

              ! compute lower outer part using (iG<)+ = iG<
              if (lower) THEN
                call zdagger(GlCSR,TCSR)
                call concat(Glout,TCSR,j1,i1)
                call destroy(TCSR)
              endif

              call destroy(GlCSR)
            endif

            call destroy(Glsub)

          endif
        end do

        call destroy(Gam)

      endif

    end do

    if (debug) then
      WRITE(*,*) '********************'
      WRITE(*,*) 'Outer_Gl done'
      WRITE(*,*) '********************'
    endif

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
    Type(z_DNS), Dimension(:,:), allocatable :: ESH
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
    call prealloc_sum(H,S,(-1.0_dp, 0.0_dp),Ec,ESH_tot)

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

    !Distruzione delle Green
    do i=2, nt
      call destroy(Gr(i,i))
      call destroy(Gr(i-1,i))
      call destroy(Gr(i,i-1))
    end do
    call destroy(Gr(1,1))

    call deallocate_blk_dns(Gr)

    !Deallocate energy-dependent matrices
    call destroy_gsm(gsmr)
    call deallocate_gsm(gsmr)

    call destroy_ESH(ESH)
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
    call prealloc_sum(H,S,(-1.0_dp, 0.0_dp),Ec,ESH_tot)

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

    !Distruzione dei blocchi fuori-diagonale
    do i=2,nbl
      call destroy(Gr(i-1,i))
      call destroy(Gr(i,i-1))
    end do
    call destroy_ESH(ESH)
    call deallocate_blk_dns(ESH)

    call create(Grm,str%mat_PL_start(nbl+1)-1,str%mat_PL_start(nbl+1)-1,0)

    Grm%rowpnt(:)=1

    do i=1,nbl
      call create(GrCSR,Gr(i,i)%nrow,Gr(i,i)%ncol,Gr(i,i)%nrow*Gr(i,i)%ncol)
      call dns2csr(Gr(i,i),GrCSR)
      !Concatena direttamente la parte immaginaria per il calcolo della doS
      zc=(-1.0_dp,0.0_dp)/pi

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

    allocate(gsm(nbl),stat=ierr)
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
