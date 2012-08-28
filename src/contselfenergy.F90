!*************************************************************************
!  Copyright (c) 2004 by Univ. Rome 'Tor Vergata'. All rights reserved.  *  
!  Authors: A. Pecchia, L. Latessa, A. Di Carlo                          *
!                                                                        *
!  Permission is hereby granted to use or copy this program for any      *
!  purpose, provided the above notices are retained on all copies, and   *
!  a notice that the code was modifided is included with the above       *
!  copyright.                                                            *
!  Permission to redistribute the code to third parties is restricted    *
!  by the licence agreement.                                             *
!*************************************************************************
!---------------------------------------------------------------------
!    Subroutine : SelfEnergy dei contatti
!---------------------------------------------------------------------
module ContSelfEnergy

 use ln_precision
 use ln_constants
 use lib_param
 use ln_structure, only : Tstruct_info
 use ln_allocation
 use mat_def
 use sparsekit_drv
 use outmatrix, only : outmat_c
 use inversions, only : block2Green, inverse 
 use clock
 use mpi_globals
 use complexbands

 implicit none
 private
 
 integer, PARAMETER :: VBT=199

  public :: surface_green !surface_green_2 
  public :: SelfEnergy
  public :: SelfEnergies
  public :: compute_contacts

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

contains 
  !--------------------------------------------------------------------
  ! SURFACE GREEN's FUNCTION USING THE DECIMATION ITERATION
  !--------------------------------------------------------------------  
  subroutine surface_green(E,HC,SC,pnegf,pnt,avncyc,GS)
    complex(dp), intent(in) :: E
    type(z_DNS), intent(in) :: HC,SC
    type(Tnegf), pointer :: pnegf
    real(dp), intent(inout) :: avncyc  ! Average num. cycles
    integer, intent(in)     :: pnt     ! Step of the energy integration
    type(z_DNS), intent(out) :: GS


    complex(kind=dp), DIMENSION(:,:), allocatable :: Ao,Bo,Co
    type(z_DNS) :: gt
    complex(kind=dp) :: mat_el

    integer :: i,i1,i2,n0,n1,n2,n3,n4,nd,npl,ngs,nkp
    integer :: ncyc,nfc,verbose,contdim,surfdim
    integer :: flag            ! flag=0 Load contact gs
                               ! flag=1 Compute 
                               ! flag=2 Compute and save
    real(kind=dp) :: dens
    character(2) :: ofcont
    character(10) :: ofkpnt
    character(10) :: ofpnt
    character(64) :: filename
    logical :: lex

    i = pnegf%activecont
    nkp = pnegf%kpoint
    flag = pnegf%ReadOldSGF
    verbose = pnegf%verbose
    contdim = pnegf%str%mat_C_end(i) - pnegf%str%mat_C_start(i) + 1
    surfdim = pnegf%str%mat_C_Start(i) - pnegf%str%mat_B_Start(i)
    ! ngs space for surface + 1 PL
    ngs = (contdim + surfdim)/2

    avncyc=0.0
    ncyc=0
    nfc=0

    write(ofcont,'(i2.2)') i

    if (nkp.le.99)  write(ofkpnt,'(i2.2)') nkp
    if (nkp.gt.99.and.id.le.999)  write(ofkpnt,'(i3.3)') nkp
    if (nkp.gt.999.and.id.le.9999)  write(ofkpnt,'(i4.4)') nkp
    if (nkp.gt.9999) stop 'ERROR: too many k-points (> 9999)'
    if (pnt.le.999) write(ofpnt,'(i3.3)') pnt  
    if (pnt.gt.999.and.pnt.le.9999) write(ofpnt,'(i4.4)') pnt  
    if (pnt.gt.9999.and.pnt.le.99999) write(ofpnt,'(i5.5)') pnt  
    if (pnt.gt.99999) stop 'ERROR: too many contour points (> 99999)'

    
    filename = 'GS'//ofcont//'_'//trim(ofkpnt)//'_'//trim(ofpnt)//'.dat'
    inquire(file=trim(pnegf%scratch_path)//filename,EXIST=lex)

    if(.not.lex.or.flag.ge.1) then

      if (id0.and.verbose.gt.VBT) call message_clock('Computing SGF '//ofpnt)

    else         !*** load from file ***

      if (id0.and.verbose.gt.VBT) call message_clock('Loading SGF '//ofpnt) 

    endif

    call create(GS,ngs,ngs)
    GS%val=(0.D0,0.D0)

    !.......... Ficticious contact ....................
    if(pnegf%FictCont(i)) then

       dens=pi*pnegf%contact_DOS(i)
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

       do i1=1,ngs
          GS%val(i1,i1)=(1.D0,0.D0)
       end do
       
       if(flag.ge.1) then
          
          !if (flag.eq.3) then
          !   ! ok only for real energies
          !   call log_allocate(Co,n5,n5)
          !   Co=conjg(transpose(SC%val))
          !   call decimation2(E,GS,HC%val,SC%val,Co,n5,n1,n2,ncyc)
          !   call log_deallocate(Co)
          !else

          call log_allocate(Ao,npl,npl)
          call log_allocate(Bo,npl,npl)
          call log_allocate(Co,npl,npl)

          Ao=E*SC%val(n1:n2,n1:n2)-HC%val(n1:n2,n1:n2)
          Bo=E*SC%val(n1:n2,n3:n4)-HC%val(n1:n2,n3:n4)
          Co=conjg(E)*SC%val(n1:n2,n3:n4)-HC%val(n1:n2,n3:n4)
          Co=conjg(transpose(Co))

          call decimation2(E,GS,Ao,Bo,Co,npl,n1,n2,ncyc)

          !call sgf_complx(E,Ao,Bo,Co,npl,GS)

          call log_deallocate(Ao)
          call log_deallocate(Bo)
          call log_deallocate(Co)
         
          ! Fill up remaining bits of the surface green's function 
          ! Add green's function of the bound layer.....
          if (n0.gt.0) then
             call create(gt,ngs,ngs)
             gt%val=(0.D0,0.D0)
             do i1=1,ngs
                gt%val(i1,i1)=(1.D0,0.D0)
             enddo
             !Here we define the Green's function related to bound states.
             call inverse(gt%val(n1:n2,n1:n2),GS%val(n1:n2,n1:n2),npl)              
             
             gt%val(1:n0,1:n0)=E*SC%val(1:n0,1:n0)-HC%val(1:n0,1:n0)         
             gt%val(1:n0,n1:n2)=E*SC%val(1:n0,n1:n2)-HC%val(1:n0,n1:n2)
             gt%val(n1:n2,1:n0)=E*SC%val(n1:n2,1:n0)-HC%val(n1:n2,1:n0) 
             
             call inverse(GS%val,gt%val,ngs)         
             !End finally the full Green's function of the contacts.
             call destroy(gt)
          endif 
          !...............................................            
          !*** save in file ***
          if (flag.eq.2) then  
             open (66,file=trim(pnegf%scratch_path)//filename, form='UNFORMATTED')

             call outmat_c(66,.false.,GS%val,ngs,ngs)
             !note: only the first of the 2 PL is saved
             close (66)
          endif
          
       else         !*** load from file ***
                
          open (65,file=trim(pnegf%scratch_path)//filename, form='UNFORMATTED')
          do
             read (65,end=100) i1,i2,mat_el
             GS%val(i1,i2)=mat_el 
          enddo
100       close(65)
          
       endif
       
       avncyc=avncyc+1.0*ncyc
       
    end if !(Fict Contact or not)

      
    if (id0.and.verbose.gt.VBT) call write_clock
    
    !if (avncyc.gt.0.and.id0.and.verbose.gt.VBT) then
    !   write(*,*) 'Number of iterations:',avncyc
    !endif

  end subroutine surface_green
  !---------------------------------------------------------------------------------------
  subroutine savedensity(g,n,E,icon,dens)

    implicit none

    integer :: n
    complex(kind=dp) :: g(n,n)
    integer :: i,icon
    complex(kind=dp) :: E,dens

    dens=(0.D0,0.D0)
    do i=1,n
      dens = dens + g(i,i) 
    enddo
    dens = (-2.D0/pi)*dens

  endsubroutine savedensity

  !--------------------------------------------------------------------
  subroutine decimation2(E,GS,Ao,Bo,Co,n,n1,n2,ncyc)

    implicit none

    complex(dp), parameter :: alfa = (1.d0,0.d0)  ! For LAPAK
    complex(dp), parameter :: beta = (0.d0,0.d0)  ! MATRIX MULT.

    integer :: n,n1,n2,ncyc

    type(z_DNS) :: GS
    complex(dp), DIMENSION(n,n) :: Ao,Bo,Co 

    integer :: i1,err, iter
    complex(dp) :: E

    complex(dp), ALLOCATABLE, DIMENSION(:,:) :: Ao_s,A1,B1,C1,A1_s
    complex(dp), ALLOCATABLE, DIMENSION(:,:) :: inAo, inAoXCo
    complex(dp), ALLOCATABLE, DIMENSION(:,:) :: inAoXBo, Gamma

    allocate(Ao_s(n,n),stat=err);
    allocate(A1(n,n),stat=err);allocate(B1(n,n),stat=err);
    allocate(C1(n,n),stat=err);allocate(A1_s(n,n),stat=err);
    allocate(inAo(n,n),stat=err);allocate(inAoXCo(n,n),stat=err);
    allocate(Gamma(n,n),stat=err);allocate(inAoXBo(n,n),stat=err);
    if (err /= 0) STOP 'no space for allocation in decimation'

    Ao_s=Ao;

    do i1=1,300
      !write(*,*) 'decimation:',i1
      !call zinv(inAo,Ao,n)
       iter = 1
       call block2Green(inAo,Ao,n,iter)

      call ZGEMM('N','N',n,n,n, alfa, inAo, n, Bo, n,  beta, inAoXBo, n)
      call ZGEMM('N','N',n,n,n, alfa, inAo, n, Co, n,  beta, inAoXCo, n) 
      call ZGEMM('N','N',n,n,n, alfa, Bo, n, inAoXCo, n, beta, Gamma, n) 

      A1_s  = Ao_s - Gamma
      A1    = Ao - Gamma

      call ZGEMM('N','N',n,n,n, -alfa, Co, n, inAoXBo, n, alfa, A1, n) 
      call ZGEMM('N','N',n,n,n,  alfa, Bo, n, inAoXBo, n, beta, B1, n)
      call ZGEMM('N','N',n,n,n,  alfa, Co, n, inAoXCo, n, beta, C1, n) 

      if((maxval(abs(Co)).le.(1.D-9)).and.(maxval(abs(C1)).le.(1.D-9))) then
        ncyc=i1
        exit;
      endif

      Ao=A1
      Ao_s=A1_s
      Bo=B1
      Co=C1

    end do
    !write(*,*) 'decimation: final inv'
    !call zinv(Gamma,A1_s,n)
    iter = 1
    call block2green(Gamma,A1_s,n,iter)

    GS%val(n1:n2,n1:n2) = Gamma(1:n,1:n)

    deallocate(Ao_s,A1,A1_s,B1,C1)
    deallocate(inAo,inAoXCo,inAoXBo,Gamma)

    return

  end subroutine decimation2

!-------------------------------------------------------------------------------

  subroutine compute_contacts_csr(Ec,pnegf,pnt,ncyc,Tlc,Tcl,SelfEneR,GS)
    complex(dp) :: Ec
    Type(Tnegf), pointer :: pnegf
    integer, intent(in) :: pnt
    Type(z_CSR), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS

    Type(z_DNS) :: GS_d
    Type(z_CSR) :: TpMt

    Integer :: nbl, ncont, i, l
    Real(dp) :: ncyc, avncyc

    nbl = pnegf%str%num_PLs
    ncont = pnegf%str%num_conts
    avncyc = 0

    ! -----------------------------------------------------------------------
    !  Calculation of contact self-energies
    ! -----------------------------------------------------------------------
    ! For the time HC and SC are dense, GS is sparse (already allocated)
    ! TM and ST are sparse, SelfEneR is allocated inside SelfEnergy
    ! -----------------------------------------------------------------------
    
    do i= 1,ncont
       pnegf%activecont=i

       call surface_green(Ec,pnegf%HC(i),pnegf%SC(i),pnegf,pnt,ncyc,GS_d)

       l = nzdrop(GS_d,EPS)
       
       call create(GS(i),GS_d%nrow,GS_d%ncol,l)
       
       call dns2csr(GS_d,GS(i))
       
       call destroy(GS_d)

       avncyc = avncyc + ncyc

       STOP 'Internal error: HMC has been changed to dns format'
       !call prealloc_sum(pnegf%HMC(i),pnegf%SMC(i),(-1.d0, 0.d0),Ec,Tlc(i))

       !call prealloc_sum(pnegf%HMC(i),pnegf%SMC(i),(-1.d0, 0.d0),conjg(Ec),TpMt)

       call zdagger(TpMt,Tcl(i))
       
       call destroy(TpMt)

       call SelfEnergy( GS(i),Tlc(i),Tcl(i),SelfEneR(i) )

    enddo

  end subroutine compute_contacts_csr
!-------------------------------------------------------------------------------

  subroutine compute_contacts_dns(Ec,pnegf,pnt,ncyc,Tlc,Tcl,SelfEneR,GS)
    complex(dp) :: Ec
    Type(Tnegf), pointer :: pnegf
    integer, intent(in) :: pnt
    Type(z_DNS), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS

    Type(z_DNS) :: TpMt

    Integer :: nbl, ncont, l, i
    Real(dp) :: ncyc, avncyc

    nbl = pnegf%str%num_PLs
    ncont = pnegf%str%num_conts
    avncyc = 0

    ! -----------------------------------------------------------------------
    !  Calculation of contact self-energies
    ! -----------------------------------------------------------------------
    ! For the time HC and SC are dense, GS is sparse (already allocated)
    ! TM and ST are sparse, SelfEneR is allocated inside SelfEnergy
    ! -----------------------------------------------------------------------
    
    do i= 1,ncont

       pnegf%activecont=i

       call surface_green(Ec,pnegf%HC(i),pnegf%SC(i),pnegf,pnt,ncyc,GS(i))
       
       avncyc = avncyc + ncyc

       call prealloc_sum(pnegf%HMC(i),pnegf%SMC(i),(-1.d0, 0.d0),Ec,Tlc(i))

       call prealloc_sum(pnegf%HMC(i),pnegf%SMC(i),(-1.d0, 0.d0),conjg(Ec),TpMt)

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

          
  !--------------------------------------------------------------------
  ! SURFACE GREEN's FUNCTION USING THE COMPLEX BANDSTRUCTURE
  ! This method is similar to the old Umerski... but it works !!
  ! It is ok for real E
  !--------------------------------------------------------------------
  ! Create Surface G.F. using complex band bloch-vectors
  ! Only for the first contact PL 
  !    ^-1    a    a   a   a ^-1
  !   g    = Z  + Z   D   D   
  !    a      11   12  22  12
  !----------------------------------------------------------------  
  subroutine sgf_complx(Ec,Z0,Z1,Z2,npl,GS)
    complex(dp), intent(in) :: Ec
    integer :: npl
    complex(dp), DIMENSION(:,:) :: Z0,Z1,Z2 
    type(z_DNS) :: GS

    ! Locals: ........................................
    type(z_DNS) :: Vr, invD12
    complex(dp), dimension(:,:), allocatable :: TT1, TT2
    real(dp), dimension(:), allocatable :: vf
    complex(dp), dimension(:), allocatable :: kzi
    integer :: i1, i2, j1, j2
    type(TStatesSummary) :: summ

    !!if(id0.and.verbose.gt.VBT) call message_clock('Computing SGF ')
    
    !  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! STEP 1: Solve Complex Bands and sort traveling bloch states
    call create(Vr,2*npl,2*npl)
    call log_allocate(vf,2*npl)
    call log_allocate(kzi,2*npl) 

    call complex_k(real(Ec),npl,Z0,Z1,kzi,Z21=Z2,Cr=Vr%val,vf=vf)

    call sort_and_normalize(2*npl,kzi,vf,Vr%val,summ)    

    !!if(id0) call write_clock     

    call log_deallocate(kzi)
    call log_deallocate(vf)

    if(summ%prop_in.ne.summ%prop_out .or. summ%evan_in.ne.summ%evan_out ) &
               STOP 'ERROR: Asymmetry found betw. IN/OUT states'   
 
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Set-up C-Matrix (Bloch vectors ordered in columns)

    ! Extract and invert D12
    call log_allocate(TT1,npl,npl)
    call log_allocate(TT2,npl,npl)

    i1 = 1; i2 = npl
    j1 = npl+1; j2 = 2*npl
    TT2=Vr%val(i1:i2,j1:j2)  !(D22)

    ! This inverse may not exist. 
    ! Probably one should invert removing the null-subspace (?)
    call create(invD12,nPL,nPL)
    call inverse(invD12%val,TT2,nPL)

    ! Compute D22*D12^-1
    TT2 = matmul(Vr%val(j1:j2,j1:j2),invD12%val)

    call destroy(Vr)

    ! Compute inverse (like G.F.)
    TT1 = matmul(Z1,TT2) + Z0
       
    call create(GS,npl,npl)   
       
    call inverse(GS%val,TT1,npl)

    deallocate(TT1,TT2)

  end subroutine sgf_complx

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


