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

 use precision
 use constants
 use parameters, only : verbose, ncont, ncdim, mbound_end, hartree, &
                        FictCont, DOS
 use allocation
 use mat_def
 use sparsekit_drv
 use outmatrix, only : outmat_c
 use inversions, only : zinv => inverse 
 use clock
 use mpi_globals
 use complexbands

 implicit none
 private
 
 integer, PARAMETER :: VBT=90

  public :: surface_green, surface_green_2, SelfEnergies, SelfEnergy

contains
  !--------------------------------------------------------------------
  ! SURFACE GREEN's FUNCTION USING THE DECIMATION ITERATION
  !--------------------------------------------------------------------  
  subroutine surface_green(E,nc,HC,SC,flag,pnt,avncyc,GS,GS_sp)
    complex(dp), intent(in) :: E
    integer, intent(in)     :: nc
    type(z_DNS), intent(in) :: HC,SC
    real(dp), intent(inout) :: avncyc  ! Average num. cycles
    integer, intent(in)     :: pnt     ! Step of the energy integration
    integer, intent(inout)  :: flag    ! flag=0 Load contact gs
                                       ! flag=1 Compute 
                                       ! flag=2 Compute and save
  
    type(z_DNS), intent(out) :: GS
    type(z_CSR), intent(out), optional :: GS_sp

    complex(kind=dp), DIMENSION(:,:), allocatable :: Ao,Bo,Co
    type(z_DNS) :: gt
    complex(kind=dp) :: mat_el

    integer :: count1,count2,count_rate,count_max
    integer :: i,i1,i2,n0,n1,n2,n3,n4,nd,npl,ngs
    integer :: ncyc,err,nfc,ncol,nnz

    real(kind=dp) :: dens,re_mat,im_mat
    character(2) :: ofcont,ofproc
    character(3) :: ofpnt
   
    logical :: lex

    avncyc=0.0
    ncyc=0
    nfc=0

    write(ofproc,'(i2.2)') id
    write(ofcont,'(i2.2)') nc
    write(ofpnt,'(i3.3)') pnt  

    inquire(file='./GS/GS'//ofcont//'_'//ofpnt//'_'//ofproc//'.dat',EXIST=lex)

    if(flag.eq.0) then
       if(lex) then
          if (id0.and.verbose.gt.VBT) call message_clock('Loading SGF ')        
       else
          flag=2
       endif
    endif
    if(flag.ge.1) then
      if (id0.and.verbose.gt.VBT) call message_clock('Computing SGF')
    endif

    i = nc

    write(ofcont,'(i2.2)') i

    ngs = (ncdim(i)+mbound_end(i))/2

    call create(GS,ngs,ngs)
    GS%val=(0.D0,0.D0)
    
    !.......... Ficticious contact ....................
    if(FictCont(i)) then
       dens=pi*DOS(i)*hartree
       nfc=nfc+1
       do i1 = 1,ncdim(i) 
          GS%val(i1,i1)=-j*dens
       end do
       !..................................................
    else
       n0 = mbound_end(i)
       n1 = n0+1                  !start of the real contact mat element
       nd = ncdim(i)-n0           !dimension of the real contact
       npl = nd/2                 !dimension of half real contact (1 PL!)
       n2 = n0+npl                !end of half real contact 
       n3 = n2+1                  !start of the second half real contact
       n4 = ncdim(i)              !end of second half real contact

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
             call zinv(gt%val(n1:n2,n1:n2),GS%val(n1:n2,n1:n2),npl)              
             
             gt%val(1:n0,1:n0)=E*SC%val(1:n0,1:n0)-HC%val(1:n0,1:n0)         
             gt%val(1:n0,n1:n2)=E*SC%val(1:n0,n1:n2)-HC%val(1:n0,n1:n2)
             gt%val(n1:n2,1:n0)=E*SC%val(n1:n2,1:n0)-HC%val(n1:n2,1:n0) 
             
             call zinv(GS%val,gt%val,ngs)         
             !End finally the full Green's function of the contacts.
             call destroy(gt)
          endif
          !...............................................            
          !*** save in file ***
          if (flag.eq.2) then  
             open (66,file='./GS/GS'//ofcont//'_'//ofpnt//'_'//ofproc//'.dat', &
                  form='UNFORMATTED')
             call outmat_c(66,.false.,GS%val(:,:),ngs,ngs)
             !note: only the first of the 2 PL is saved
             close (66)
          endif
          
       else         !*** load from file ***
          
          
          open (65,file='./GS/GS'//ofcont//'_'//ofpnt//'_'//ofproc//'.dat', &
               form='UNFORMATTED')
          do
             read (65,end=100) i1,i2,mat_el
             GS%val(i1,i2)=mat_el 
             GS%val(i2,i1)=mat_el
          enddo
100       close(65)
          
       endif
       
       avncyc=avncyc+1.0*ncyc
       
    end if !(Fict Contact or not)

    if (present(GS_sp)) then
       call create(GS_sp,ngs,ngs,nzdrop(GS,EPS))
       call dns2csr(GS,GS_sp)
       call destroy(GS)
    endif
    
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

    integer :: i1,err
    complex(dp) :: E,Ecc,dens

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
      call zinv(inAo,Ao,n)

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
    call zinv(Gamma,A1_s,n)

    GS%val(n1:n2,n1:n2) = Gamma(1:n,1:n)

    deallocate(Ao_s,A1,A1_s,B1,C1)
    deallocate(inAo,inAoXCo,inAoXBo,Gamma)

    return

  end subroutine decimation2


  !--------------------------------------------------------------------
  subroutine SelfEnergies(E,GS,Tlc,Tcl,SelfEneR)
    complex(kind=dp) :: E
    type(z_CSR) :: GS(10),Tlc(10),Tcl(10)

    !OUTPUT
    type(z_CSR) :: SelfEneR(10)

    integer :: i 

    do i = 1,ncont
       call SelfEnergy( GS(i),Tlc(i),Tcl(i),SelfEneR(i) )
    end do

  end subroutine SelfEnergies

  ! -------------------------------------------------------------
  subroutine SelfEnergy(GS,Tlc,Tcl,SelfEneR)

    type(z_CSR) :: GS,Tlc,Tcl
    !OUTPUT

    type(z_CSR) :: SelfEneR    

    ! locals 
    integer :: i,n,err,nrow,ncol
    type(z_CSR) :: TG,TT

    
    if(Tlc%nnz.eq.0) then
       call create(SelfEneR,Tlc%nrow,Tlc%nrow,1)
       SelfEneR%nnz=0
       if (id0) write(*,*) 'SelfEne= 0'    
       return
    endif
    
    call prealloc_mult(Tlc,GS,TG)
    !if (id0) write(*,*) 'Tlc GS', TG%nrow,TG%nnz
    
    call prealloc_mult(TG,Tcl,SelfEneR)
    
    !deallocate (TG,TT)
    call destroy(TG)    
    
  end subroutine SelfEnergy

  !--------------------------------------------------------------------
  ! SURFACE GREEN's FUNCTION USING THE COMPLEX BANDSTRUCTURE
  ! This method is similar to the old Umerski... but it works !!
  ! It is ok for real E
  !--------------------------------------------------------------------
  subroutine surface_green_2(Ec,nc,HC,SC,GS)
    complex(dp), intent(in) :: Ec
    integer, intent(in)     :: nc
    type(z_DNS), intent(in) :: HC,SC
    type(z_CSR), intent(out) :: GS

    ! Locals: ........................................
    integer :: PLdim, n0,n1,n2,n3,n4,n5,nmax
    type(z_DNS) :: Vr, Z0, Z1, Z2
    type(z_DNS) :: D11,D12,D22,L11,GS_d
    real(dp), dimension(:), allocatable :: vf
    complex(dp), dimension(:), allocatable :: kzi
    type(TStatesSummary) :: summ


    if(id0.and.verbose.gt.VBT) call message_clock('Computing SGF ')

    n0 = mbound_end(nc)
    n1 = n0+1                   !start of the real contact mat element
    n4 = ncdim(nc)-n0           !dimension of the real contact
    n5 = n4/2                   !dimension of half real contact (1 PL!)
    n2 = n0+n5                  !end of half real contact 
    n3 = n2+1                   !start of the second half real contact
    nmax=ncdim(nc)              !end of second half real contact

    ! Solve Complex Bands
    call create(Vr,2*n5,2*n5)
    call create(Z0,n5,n5)
    call create(Z1,n5,n5)
    call log_allocate(vf,2*n5)
    call log_allocate(kzi,2*n5) 

    Z0%val = Ec * SC%val(n1:n2,n1:n2) - HC%val(n1:n2,n1:n2)
    Z1%val = Ec * SC%val(n1:n2,n3:nmax) - HC%val(n1:n2,n3:nmax)
    Z2%val = Ec * SC%val(n3:nmax,n1:n2) - HC%val(n3:nmax,n1:n2) 

    call complex_k(real(Ec),n5,Z0%val,Z1%val,kzi,Z21=Z2%val,Cr=Vr%val,vf=vf)

    call sort_and_normalize(2*n5,kzi,vf,Vr%val,summ)    

    if(id0) call write_clock     

    call log_deallocate(kzi)
    call log_deallocate(vf)

    if(summ%prop_in.ne.summ%prop_out .or. summ%evan_in.ne.summ%evan_out ) &
               STOP 'ERROR: Asymmetry found betw. IN/OUT states'   
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Set-up C-Matrix (Bloch vectors ordered in columns)

    call create(D12,n5,n5)
    call create(D22,n5,n5)
 
    D12%val = Vr%val(1:n5,n5+1:2*n5)        !D12
    D22%val = Vr%val(n5+1:2*n5,n5+1:2*n5)   !D22      
    
    call destroy(Vr)

    call create(D11,n5,n5)    
    D11%val = matmul(Z0%val,D12%val) + matmul(Z1%val,D22%val)

    call destroy(Z0,Z1,D22)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Compute surface G.F.
    call create(L11,n5,n5)

    call zinv(L11%val,D11%val,n5)

    call destroy(D11)
    
    call prealloc_mult(D12,L11,GS_d)
   
    call destroy(L11,D12)

    call create(GS,GS_d%nrow,GS_d%ncol,nzdrop(GS_d,EPS))
    call dns2csr(GS_d,GS)

    call destroy(GS_d)

    if(id0.and.verbose.gt.VBT) call write_clock

  end subroutine Surface_green_2

end module ContSelfEnergy
