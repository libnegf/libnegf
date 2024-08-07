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


!---------------------------------------------------------------------
!    Subroutine : SelfEnergy for contacts
!---------------------------------------------------------------------
module ContSelfEnergy

 use ln_precision
 use ln_constants
 use lib_param
 use ln_structure, only : Tstruct_info
 use ln_allocation
 use mat_def
 use sparsekit_drv
 use outmatrix, only : outmat_c, inmat_c
 use inversions, only : compGreen, inverse
 use clock
 use mpi_globals
 use ln_cache, only : TMatLabel, get_string_label
#:if defined("MPI")
 use libmpifx_module, only : mpifx_reduceip
#:endif
#:if defined("GPU")
 use iso_c_binding
#:endif
 implicit none
 private

 integer, PARAMETER :: VBT=70

  public :: surface_green
  public :: SelfEnergy
  public :: SelfEnergies
  public :: compute_contacts
  public :: decimation
#:if defined("GPU")
  public :: decimation_gpu
#:endif

  interface SelfEnergy
     module procedure SelfEnergy_csr
     module procedure SelfEnergy_dns
  end interface

  interface SelfEnergies
     module procedure SelfEnergies_csr
     module procedure SelfEnergies_dns
  end interface

  interface compute_contacts
     module procedure compute_contacts_csr
     module procedure compute_contacts_dns
  end interface

#:if defined("GPU")
  interface
     integer(c_int) function cu_Cdecimation(hcublas, hcusolver, h_Go_out, h_Ao_in, h_Bo_in, h_Co_in, &
                     & n, tf32, ncyc, one, mone, zero, SGFACC) bind(C, name='cu_Cdecimation')
        use iso_c_binding
        import cublasHandle
        import cusolverDnHandle

        type(cublasHandle), value :: hcublas
        type(cusolverDnHandle), value :: hcusolver
        type(c_ptr), value :: h_Go_out
        type(c_ptr), value :: h_Ao_in
        type(c_ptr), value :: h_Bo_in
        type(c_ptr), value :: h_Co_in
        integer(c_int), value :: n
        type(c_ptr), value :: ncyc
        integer(c_int), value :: tf32
        complex(c_float_complex) :: one
        complex(c_float_complex) :: mone
        complex(c_float_complex) :: zero
        real(c_float) :: SGFACC
     end function

     integer(c_int) function cu_Zdecimation(hcublas, hcusolver, h_Go_out, h_Ao_in, h_Bo_in, h_Co_in, &
                     & n, tf32, ncyc, one, mone, zero, SGFACC) bind(C, name='cu_Zdecimation')
        use iso_c_binding
        import cublasHandle
        import cusolverDnHandle

         type(cublasHandle), value :: hcublas
        type(cusolverDnHandle), value :: hcusolver
        type(c_ptr), value :: h_Go_out
        type(c_ptr), value :: h_Ao_in
        type(c_ptr), value :: h_Bo_in
        type(c_ptr), value :: h_Co_in
        integer(c_int), value :: n
        type(c_ptr), value :: ncyc
        integer(c_int), value :: tf32
        complex(c_double_complex) :: one
        complex(c_double_complex) :: mone
        complex(c_double_complex) :: zero
        real(c_double) :: SGFACC
     end function
  end interface

  interface decimation_gpu
     module procedure decimation_gpu_sp
     module procedure decimation_gpu_dp
  end interface
#:endif
contains
  !--------------------------------------------------------------------
  ! SURFACE GREEN's FUNCTION USING THE DECIMATION ITERATION
  !--------------------------------------------------------------------
  subroutine surface_green(E,HC,SC,pnegf,ncyc,GS)
    complex(dp), intent(in) :: E
    type(z_DNS), intent(in) :: HC,SC
    type(Tnegf) :: pnegf
    integer, intent(out) :: ncyc
    type(z_DNS), intent(out) :: GS


    complex(kind=dp), DIMENSION(:,:), allocatable :: Ao,Bo,Co,Go
    type(z_DNS) :: gt

    integer :: ii,i1,n0,n1,n2,n3,n4,nd,npl,ngs
    integer :: pnt,nfc,verbose,contdim,surfdim
    integer :: flag            ! flag=0 Load contact gs
                               ! flag=1 Compute
                               ! flag=2 Compute and save
    real(kind=dp) :: dens
    character(10) :: ofpnt
    logical :: lex
    type(TMatLabel) :: label

    pnt = pnegf%iE_path    ! Step of the energy integration
    ii = pnegf%activecont
    label%kpoint = pnegf%ikpoint
    label%energy_point = pnt
    label%spin = pnegf%spin
    label%row_block = ii
    label%col_block = 0

    flag = pnegf%ReadOldSGF
    verbose = pnegf%verbose
    contdim = pnegf%str%cont_dim(ii)
    surfdim = pnegf%str%mat_C_Start(ii) - pnegf%str%mat_B_Start(ii)

    !      Surf    PL1   PL2
    ! contdim = surfdim + 2 PL =>
    ! ngs = surfdim + PL = surfdim + (contdim-surfdim)/2
    ngs = (surfdim + contdim)/2

    ncyc=0
    nfc=0

    lex = pnegf%surface_green_cache%is_cached(label)

    call get_string_label(pnt, ofpnt)

    if (.not.lex .and. flag.eq.0) then
        flag = 2
    endif

    if(flag.ge.1) then
        if (id0.and.verbose.gt.VBT) call message_clock(trim("Computing SGF "//ofpnt))
    else         !*** load from file ***
        if (id0.and.verbose.gt.VBT) call message_clock(trim("Loading SGF "//ofpnt))
    endif

    call create(GS,ngs,ngs)
    GS%val=zero

    !.......... Ficticious contact ....................
    if(pnegf%cont(ii)%FictCont) then

       dens=pi*pnegf%cont(ii)%contact_DOS
       nfc=nfc+1
       do i1 = 1,ngs
          GS%val(i1,i1)=-j*dens
       end do

    else

       !   1      n0 n1 n2 n3  n4
       !   +--------+-----+-----+
       !      Surf    PL1   PL2

       n0 = surfdim
       n1 = n0+1                  !start of the real contact mat element
       nd = contdim - n0          !dimension of the real contact
       npl = nd/2                 !dimension of half real contact (1 PL!)
       n2 = n0+npl                !end of half real contact
       n3 = n2+1                  !start of the second half real contact
       n4 = contdim               !end of second half real contact

       if(flag.ge.1) then

          call log_allocate(Ao,npl,npl)
          call log_allocate(Bo,npl,npl)
          call log_allocate(Co,npl,npl)
          call log_allocate(Go,npl,npl)

          Ao=E*SC%val(n1:n2,n1:n2)-HC%val(n1:n2,n1:n2)
          Bo=E*SC%val(n1:n2,n3:n4)-HC%val(n1:n2,n3:n4)
          Co=conjg(E)*SC%val(n1:n2,n3:n4)-HC%val(n1:n2,n3:n4)
          Co=conjg(transpose(Co))
#:if defined("GPU")
          call decimation_gpu(pnegf,Go,Ao,Bo,Co,npl,.false.,ncyc)
#:else
          call decimation(Go,Ao,Bo,Co,npl,ncyc)
#:endif
          call log_deallocate(Ao)
          call log_deallocate(Bo)
          call log_deallocate(Co)

          ! Fill up remaining bits of the surface green's function
          ! Add green's function of the bound layer.....
          if (n0.gt.0) then
             call create(gt,ngs,ngs)
             gt%val=zero
             do i1=1,ngs
                gt%val(i1,i1)=one
             enddo
             !Here we define the Green's function related to bound states.
             call inverse(gt%val(n1:n2,n1:n2),Go,npl)
             call log_deallocate(Go)

             gt%val(1:n0,1:n0)=E*SC%val(1:n0,1:n0)-HC%val(1:n0,1:n0)
             gt%val(1:n0,n1:n2)=E*SC%val(1:n0,n1:n2)-HC%val(1:n0,n1:n2)
             gt%val(n1:n2,1:n0)=E*SC%val(n1:n2,1:n0)-HC%val(n1:n2,1:n0)

             call inverse(GS%val,gt%val,ngs)
             !End finally the full Green's function of the contacts.
             call destroy(gt)
          else
             GS%val = Go
             call log_deallocate(Go)
          endif

          !print*,'GS',maxval(abs(GS%val))
          !...............................................
          !*** save in file ***
          if (flag.eq.2) then

            call pnegf%surface_green_cache%add(GS, label)

          endif

       else         !*** load from file ***

         call pnegf%surface_green_cache%retrieve(GS, label)

       endif

    end if !(Fict Contact or not)


    if (id0.and.verbose.gt.VBT) call write_clock

  end subroutine surface_green
  !---------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------
  subroutine decimation(Go,Ao,Bo,Co,n,ncyc)
    integer, intent(in) :: n
    complex(dp), DIMENSION(n,n), intent(out) :: Go
    complex(dp), DIMENSION(n,n), intent(inout) :: Ao,Bo,Co
    integer, intent(out) :: ncyc

    complex(dp), ALLOCATABLE, DIMENSION(:,:) :: Ao_s, A1, B1, C1
    complex(dp), ALLOCATABLE, DIMENSION(:,:) :: GoXCo
    complex(dp), ALLOCATABLE, DIMENSION(:,:) :: GoXBo, Self
    integer :: i1, err
    logical :: okCo = .false.

    call log_allocate(Ao_s, n, n)
    Ao_s=Ao;

    do i1 = 1, 300
      ncyc=i1

      call compGreen(Go,Ao,n)

      call log_allocate(GoXCo, n, n)
      call ZGEMM('N','N',n,n,n, one, Go, n, Co, n,  zero, GoXCo, n)

      call log_allocate(C1, n, n)
      call ZGEMM('N','N',n,n,n,  one, Co, n, GoXCo, n, zero, C1, n)

      if (maxval(abs(C1)).le.SGFACC) then
         if (okCo) then
            call log_deallocate(GoXCo)
            call log_deallocate(C1)
            exit;
         else
            okCo = .true.
         endif
      else
         okCo = .false.
      endif

      call log_allocate(Self, n, n)
      call ZGEMM('N','N',n,n,n, one, Bo, n, GoXCo, n, zero, Self, n)
      call log_deallocate(GoXCo)
      Ao_s  = Ao_s - Self
      Ao    = Ao - Self
      call log_deallocate(Self)

      call log_allocate(GoXBo, n, n)
      call ZGEMM('N','N',n,n,n, one, Go, n, Bo, n,  zero, GoXBo, n)

      call log_allocate(B1, n, n)
      call ZGEMM('N','N',n,n,n,  one, Bo, n, GoXBo, n, zero, B1, n)
      Bo = B1
      call log_deallocate(B1)

      call ZGEMM('N','N',n,n,n, -one, Co, n, GoXBo, n, one, Ao, n)

      Co = C1
      call log_deallocate(C1)
      call log_deallocate(GoXBo)

    end do

    ncyc=i1

    call compGreen(Go,Ao_s,n)
    call log_deallocate(Ao_s)

  end subroutine decimation

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
#:if defined("GPU")
  subroutine decimation_gpu_sp(negf, Go_out, Ao_in, Bo_in, Co_in, n, tf32, ncyc)
     implicit none
     type(Tnegf), intent(in) :: negf
     integer, intent(in) :: n
    complex(sp), dimension(n,n), intent(out), target :: Go_out
    complex(sp), dimension(n,n), intent(in), target :: Ao_in, Bo_in, Co_in
     logical, intent(in) :: tf32
     integer, intent(out) :: ncyc

    complex(sp) :: one = (1.0_sp, 0.0_sp)
    complex(sp) :: mone = (-1.0_sp, 0.0_sp)
    complex(sp) :: zero = (0.0_sp, 0.0_sp)
    integer :: istat, tf
    type(cublasHandle) :: hh
    type(cusolverDnHandle) :: hhsol

    if (tf32) then
       tf = 1
    else
       tf = 0
    endif
    hh = negf%hcublas
    hhsol = negf%hcusolver

    istat = cu_Cdecimation(hh, hhsol, c_loc(Go_out), c_loc(Ao_in), c_loc(Bo_in), c_loc(Co_in), n, tf, c_loc(ncyc), &
            one, mone, zero, real(SGFACC,sp)*n*n)

  end subroutine decimation_gpu_sp

  subroutine decimation_gpu_dp(negf, Go_out, Ao_in, Bo_in, Co_in, n, tf32, ncyc)
    implicit none
    type(Tnegf), intent(in) :: negf
    integer, intent(in) :: n
    complex(dp), dimension(n,n), intent(out), target :: Go_out
    complex(dp), dimension(n,n), intent(in), target :: Ao_in, Bo_in, Co_in
    logical, intent(in) :: tf32
    integer, intent(out), target :: ncyc

    complex(dp) :: one = (1.0_dp, 0.0_dp)
    complex(dp) :: mone = (-1.0_dp, 0.0_dp)
    complex(dp) :: zero = (0.0_dp, 0.0_dp)
    integer :: istat, tf
    type(cublasHandle) :: hh
    type(cusolverDnHandle) :: hhsol

    if (tf32) then
       tf = 1
    else
       tf = 0
    endif

    hh = negf%hcublas
    hhsol = negf%hcusolver

    istat = cu_Zdecimation(hh, hhsol, c_loc(Go_out), c_loc(Ao_in), c_loc(Bo_in), c_loc(Co_in), n, tf, c_loc(ncyc), &
            one, mone, zero, real(SGFACC,dp)*n*n)

  end subroutine decimation_gpu_dp

#:endif

!-------------------------------------------------------------------------------

  subroutine compute_contacts_csr(Ec,pnegf,avncyc,Tlc,Tcl,SelfEneR,GS)
    complex(dp), intent(in) :: Ec
    Type(Tnegf), intent(inout) :: pnegf
    real(dp), intent(out) :: avncyc
    Type(z_CSR), Dimension(MAXNCONT), intent(in) :: Tlc, Tcl
    Type(z_CSR), Dimension(MAXNCONT), intent(out) :: SelfEneR, GS


    Type(z_DNS) :: GS_d
    Type(z_CSR) :: TpMt

    Integer :: ncyc, ncont, i, l

    ncont = pnegf%str%num_conts
    avncyc = 0.0_dp

    error stop 'Internal error: HMC has been changed to dns format'
    ! -----------------------------------------------------------------------
    !  Calculation of contact self-energies
    ! -----------------------------------------------------------------------
    ! For the time HC and SC are dense, GS is sparse (already allocated)
    ! TM and ST are sparse, SelfEneR is allocated inside SelfEnergy
    ! -----------------------------------------------------------------------

    do i= 1,ncont
       pnegf%activecont=i

       call surface_green(Ec,pnegf%cont(i)%HC,pnegf%cont(i)%SC,pnegf,ncyc,GS_d)

       l = nzdrop(GS_d,EPS)

       call create(GS(i),GS_d%nrow,GS_d%ncol,l)

       call dns2csr(GS_d,GS(i))

       call destroy(GS_d)

       avncyc = avncyc + 1.0_dp*ncyc/ncont

       call zdagger(TpMt,Tcl(i))

       call destroy(TpMt)

       call SelfEnergy( GS(i),Tlc(i),Tcl(i),SelfEneR(i) )

    enddo

  end subroutine compute_contacts_csr

  !-------------------------------------------------------------------------------
  subroutine compute_contacts_dns(Ec,pnegf,avncyc,Tlc,Tcl,SelfEneR,GS)
    complex(dp), intent(in) :: Ec
    Type(Tnegf), intent(inout) :: pnegf
    real(dp), intent(out) :: avncyc
    Type(z_DNS), Dimension(MAXNCONT), intent(inout) :: Tlc, Tcl
    Type(z_DNS), Dimension(MAXNCONT), intent(out) :: SelfEneR, GS


    Type(z_DNS) :: TpMt

    Integer :: ncyc, ncont, i

    ncont = pnegf%str%num_conts
    avncyc = 0.0_dp

    ! -----------------------------------------------------------------------
    !  Calculation of contact self-energies
    ! -----------------------------------------------------------------------
    ! For the time HC and SC are dense, GS is sparse (already allocated)
    ! TM and ST are sparse, SelfEneR is allocated inside SelfEnergy
    ! -----------------------------------------------------------------------

    do i= 1,ncont

       pnegf%activecont=i

       call surface_green(Ec,pnegf%cont(i)%HC,pnegf%cont(i)%SC,pnegf,ncyc,GS(i))

       avncyc = avncyc + ncyc*1.0_dp/ncont

       call prealloc_sum(pnegf%cont(i)%HMC,pnegf%cont(i)%SMC,minusone,Ec,Tlc(i))

       call prealloc_sum(pnegf%cont(i)%HMC,pnegf%cont(i)%SMC,minusone,conjg(Ec),TpMt)

       call zdagger(TpMt,Tcl(i))

       call destroy(TpMt)

       call SelfEnergy( GS(i),Tlc(i),Tcl(i),SelfEneR(i) )

    enddo

  end subroutine compute_contacts_dns

  !--------------------------------------------------------------------
  subroutine SelfEnergies_csr(E,ncont,GS,Tlc,Tcl,SelfEneR)
    complex(dp) :: E
    integer :: ncont
    type(z_CSR) :: GS(MAXNCONT),Tlc(MAXNCONT),Tcl(MAXNCONT)
    !OUTPUT
    type(z_CSR) :: SelfEneR(MAXNCONT)
    integer :: i

    do i = 1,ncont
       call SelfEnergy( GS(i),Tlc(i),Tcl(i),SelfEneR(i) )
    end do

  end subroutine SelfEnergies_csr

  ! -------------------------------------------------------------
  subroutine SelfEnergy_csr(GS,Tlc,Tcl,SelfEneR)

    type(z_CSR) :: GS,Tlc,Tcl
    !OUTPUT

    type(z_CSR) :: SelfEneR

    ! locals
    type(z_CSR) :: TG


    if(Tlc%nnz.eq.0) then
       call create(SelfEneR,Tlc%nrow,Tlc%nrow,1)
       SelfEneR%nnz=0
       if (id0) write(*,*) '(SelfEnergy) WARNING: SelfEne= 0'
       return
    endif

    call prealloc_mult(Tlc,GS,TG)

    call prealloc_mult(TG,Tcl,SelfEneR)

    call destroy(TG)

  end subroutine SelfEnergy_csr


  !--------------------------------------------------------------------
  subroutine SelfEnergies_dns(E,ncont,GS,Tlc,Tcl,SelfEneR)
    complex(dp) :: E
    integer :: ncont
    type(z_DNS) :: GS(MAXNCONT),Tlc(MAXNCONT),Tcl(MAXNCONT)
    !OUTPUT
    type(z_DNS) :: SelfEneR(MAXNCONT)
    integer :: i

    do i = 1,ncont
       call SelfEnergy( GS(i),Tlc(i),Tcl(i),SelfEneR(i) )
    end do

  end subroutine SelfEnergies_dns

  ! -------------------------------------------------------------
  subroutine SelfEnergy_dns(GS,Tlc,Tcl,SelfEneR)

    type(z_DNS) :: GS,Tlc,Tcl
    !OUTPUT

    type(z_DNS) :: SelfEneR

    ! locals
    type(z_DNS) :: TG

    call prealloc_mult(Tlc,GS,TG)

    call prealloc_mult(TG,Tcl,SelfEneR)

    call destroy(TG)

  end subroutine SelfEnergy_dns


end module ContSelfEnergy





!    call create(D12,npl,npl)
!    call create(D22,npl,npl)
!
!    D12%val = Vr%val(i1:i2,j1:j2)   !D12
!    D22%val = Vr%val(j1:j2,j1:j2)   !D22
!
!    call destroy(Vr)
!
!    call create(D11,npl,npl)
!    D11%val = matmul(Z0,D12%val) + matmul(Z1,D22%val)!
!
!    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!    ! Compute surface G.F.
!    call create(L11,npl,npl)
!
!    call zinv(L11%val,D11%val,npl)
!
!    call destroy(D11)
!
!    call prealloc_mult(D12,L11,GS)
!
!    call destroy(L11,D12)
