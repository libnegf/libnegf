module libnegf

 use ln_precision
 use ln_constants
 use ln_allocation
 use fermi_dist
 use lib_param
 use mpi_globals, only : id, numprocs, id0
 use input_output
 use ln_structure
 use sparsekit_drv
 use inversions
 use iterative
 use iterative_dns
 use rcm_module
 use mat_def
 use ln_extract
 use contselfenergy
 use clock

 implicit none
 private

 public :: init_negf, destroy_matrices, destroy_negf
 public :: compute_dos, compute_contacts
 public :: contour_int_n, contour_int_p, real_axis_int, contour_int
 public :: compute_current, integrate
 public :: reorder, partition, partition2
 public :: sort, swap

 integer, PARAMETER :: VBT=70

contains
  
  subroutine init_negf(negf)
    type(Tnegf), pointer :: negf
    Integer :: ncont, nbl
    Integer, dimension(:), allocatable :: PL_end, cont_end, surf_end, cblk
    character(11) :: fmtstring
 
    if(negf%form%formatted) then
       fmtstring = 'formatted'
    else
       fmtstring = 'unformatted'
    endif

    open(101, file=negf%file_re_H, form=trim(fmtstring))
    open(102, file=negf%file_im_H, form=trim(fmtstring))   !open imaginary part of H

    call read_H(101,102,negf%H,negf%form)

    close(101)
    close(102)

    if(.not.negf%isSid) then
       open(101, file=negf%file_re_S, form=trim(fmtstring))
       open(102, file=negf%file_im_S, form=trim(fmtstring))   !open imaginary part of S

       call read_H(101,102, negf%S,negf%form)
       
       close(101)
       close(102)

    else
       ! create an Id matrix for S
       call create_id( negf%S,negf%H%nrow) 

    endif

    open(101, file=negf%file_struct, form='formatted')  

    read(101,*) ncont
    read(101,*) nbl

    call log_allocate(PL_end,nbl)
    call log_allocate(cblk,ncont)
    call log_allocate(cont_end,ncont)
    call log_allocate(surf_end,ncont)

    read(101,*) PL_end(1:nbl)
    read(101,*) cont_end(1:ncont)
    read(101,*) surf_end(1:ncont)
    read(101,*) cblk(1:ncont)

    read(101,*) negf%mu_n
    read(101,*) negf%mu_p
    read(101,*) negf%Ec
    read(101,*) negf%Ev
    read(101,*) negf%DeltaEc
    read(101,*) negf%DeltaEv
    read(101,*) negf%Emin
    read(101,*) negf%Emax
    read(101,*) negf%Estep
    read(101,*) negf%kbT
    read(101,*) negf%wght
    read(101,*) negf%Np_n(1)
    read(101,*) negf%Np_n(2)
    read(101,*) negf%Np_p(1)
    read(101,*) negf%Np_p(2)
    read(101,*) negf%Np_real
    read(101,*) negf%n_kt
    read(101,*) negf%n_poles
    read(101,*) negf%spin
    read(101,*) negf%delta
    read(101,*) negf%nLDOS
    call log_allocatep(negf%LDOS,2,negf%nLDOS)
    read(101,*) negf%LDOS
    read(101,*) negf%Efermi(1:ncont)
    read(101,*) negf%mu(1:ncont)

    close(101)

    print*, '(init_negf) LDOS:'
    print*, negf%LDOS

    print*, '(init_negf) before create T'
    call create_Tstruct(ncont, nbl, PL_end, cont_end, surf_end, cblk, negf%str)
    print*, '(init_negf) done'

    call log_deallocate(PL_end)
    call log_deallocate(cblk)
    call log_deallocate(cont_end)
    call log_deallocate(surf_end)

  end subroutine init_negf

!--------------------------------------------------------------------
  subroutine init_structure(negf,ncont,nbl,PL_end,cont_end,surf_end,cblk)
    type(Tnegf), pointer :: negf
    Integer :: ncont, nbl
    Integer, dimension(:) :: PL_end, cont_end, surf_end, cblk

    call create_Tstruct(ncont, nbl, PL_end, cont_end, surf_end, cblk, negf%str)

  end subroutine init_structure

!--------------------------------------------------------------------
  subroutine destroy_negf(negf)
    type(Tnegf), pointer :: negf   
    integer :: i

    call destroy_matrices(negf)

    call kill_Tstruct(negf%str) 

    if (associated(negf%LDOS)) call log_deallocatep(negf%LDOS)

  end subroutine destroy_negf
!--------------------------------------------------------------------

  subroutine destroy_matrices(negf)
    type(Tnegf), pointer :: negf   
    integer :: i

    if (allocated(negf%H%nzval)) then
       !print*,'(destroy) deallocate negf%H',%LOC(negf%H%nzval)
       call destroy(negf%H) 
    end if
    if (allocated(negf%S%nzval)) then
       !print*,'(destroy) deallocate negf%S',%LOC(negf%S%nzval)
       call destroy(negf%S) 
    end if
    if (allocated(negf%rho%nzval)) then
       !print*,'(destroy) deallocate negf%rho',%LOC(negf%rho%nzval)
       call destroy(negf%rho) 
    end if
    if (allocated(negf%rho_eps%nzval)) then
       !print*,'(destroy) deallocate negf%rho_eps',%LOC(negf%rho_eps%nzval)
       call destroy(negf%rho_eps) 
    end if


    if (allocated(negf%HM%nzval)) call destroy(negf%HM)
    if (allocated(negf%SM%nzval)) call destroy(negf%SM)

    do i=1,negf%str%num_conts
       if (allocated(negf%HC(i)%val)) call destroy(negf%HC(i))
       if (allocated(negf%SC(i)%val)) call destroy(negf%SC(i))
       if (allocated(negf%HMC(i)%nzval)) call destroy(negf%HMC(i))
       if (allocated(negf%SMC(i)%nzval)) call destroy(negf%SMC(i))
    enddo

  end subroutine destroy_matrices

!--------------------------------------------------------------------

  subroutine compute_dos(negf) 
    type(Tnegf), pointer :: negf

    Type(z_CSR), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    Type(z_CSR) ::  Gr

    integer :: N, k, i, i1, it
    integer :: outer, nbl, ncont

    real(dp) :: ncyc
    complex(dp) :: Ec

    outer = 1

    it = negf%iteration
    nbl = negf%str%num_PLs
    ncont = negf%str%num_conts

    N = nint((negf%Emax-negf%Emin)/negf%Estep)

    print*,'N=',N
    print*,'delta=',negf%delta


    open(101,file='dos.dat')

    do i=1,N
  
       Ec=(negf%Emin+i*negf%Estep)+negf%delta*(0.d0,1.d0)

       call compute_contacts(Ec,negf,i,ncyc,Tlc,Tcl,SelfEneR,GS)
   
       call calls_eq_mem_dns(negf%HM,negf%SM,Ec,SelfEneR,Tlc,Tcl,GS,Gr,negf%str,outer)

       do i1=1,ncont
          call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
       enddo

       negf%dos = -aimag( trace(Gr) )/pi

       write(101,*) real(Ec), negf%dos
    !   print*, real(Ec), negf%dos, Gr%nnz

       call destroy(Gr)

    enddo

    call writememinfo(6)

    close(101)   

  end subroutine compute_dos

!-------------------------------------------------------------------------------

  subroutine compute_contacts(Ec,pnegf,pnt,ncyc,Tlc,Tcl,SelfEneR,GS)

    Type(Tnegf), pointer :: pnegf
    integer, intent(in) :: pnt
    Type(z_CSR), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS


    Type(z_DNS) :: GS_d
    Type(z_CSR) :: TpMt

    Integer :: nbl, ncont, l, i1
    Real(dp) :: ncyc, avncyc
    complex(dp) :: Ec

    nbl = pnegf%str%num_PLs
    ncont = pnegf%str%num_conts
    avncyc = 0

    if(pnegf%iteration .eq. 1) then 
       pnegf%ReadOldSGF = 2
    else
       pnegf%ReadOldSGF = 0
    endif

       ! -----------------------------------------------------------------------
       !  Calculation of contact self-energies
       ! -----------------------------------------------------------------------
       ! For the time HC and SC are dense, GS is sparse (already allocated)
       ! TM and ST are sparse, SelfEneR is allocated inside SelfEnergy
       ! -----------------------------------------------------------------------

       do l=1,ncont
          pnegf%activecont=l
          call surface_green(Ec,pnegf%HC(l),pnegf%SC(l),pnegf,pnt,ncyc,GS_d)
  
          i1 = nzdrop(GS_d,EPS)
  
          call create(GS(l),GS_d%nrow,GS_d%ncol,i1)
       
          call dns2csr(GS_d,GS(l))
          call destroy(GS_d)
     
          avncyc = avncyc + ncyc

       enddo

 
       ! -- GF Calculation ----------------------------------------------------
       ! 
       !Tlc= Ec*ST - TM 
       !Tcl= (conjg(Ec)*ST - TM )
       !Array di GS sparse.
       
       !Tlc: matrici di interazione (ES-H) device-contatti (l=layer,c=contact)
       !Tcl: matrici di interazione (ES-H) contatti-device (l=layer,c=contact
       
       do i1=1,ncont
          call prealloc_sum(pnegf%HMC(i1),pnegf%SMC(i1),(-1.d0, 0.d0),Ec,Tlc(i1))

          call prealloc_sum(pnegf%HMC(i1),pnegf%SMC(i1),(-1.d0, 0.d0),conjg(Ec),TpMt)

          call zdagacsr(TpMt,Tcl(i1))

          call destroy(TpMt)

       enddo

       !if (id0.and.pnegf%verbose.gt.60) call message_clock('Compute Green`s funct ')

       call SelfEnergies(Ec,ncont,GS,Tlc,Tcl,SelfEneR)

     end subroutine compute_contacts

!-------------------------------------------------------------------------------

     subroutine compute_current(negf)


       type(Tnegf), pointer :: negf

       Type(z_CSR), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS

       Real(dp), Dimension(:), allocatable :: currents 
       Real(dp), Dimension(:), allocatable :: TUN_MAT
       Real(dp), Dimension(:,:), allocatable :: TUN_TOT_MAT  
       !Real(kind=dp), Dimension(:,:,:), allocatable :: TUN_PMAT, TUN_TOT_PMAT      

       Real(dp), Dimension(:), allocatable :: mumin_array, mumax_array
       Real(dp), Dimension(:), allocatable :: LEDOS
       Real(dp), Dimension(:,:), allocatable :: LEDOS_MAT
       Real(dp) :: avncyc, ncyc, Ec_min, mumin, mumax, telec

       Complex(dp) :: Ec

       Integer, Dimension(:), pointer :: cblk, indblk
       Integer :: i, Nstep, npid, istart, iend, i1
       Integer :: size_ni, size_nf, icpl, ncont, icont
       Integer :: ierr, nbl, nit, nft
       Integer :: iLDOS

       Character(6) :: ofKP
       Logical, Parameter :: do_LEDOS=.TRUE.

       nbl = negf%str%num_PLs
       ncont = negf%str%num_conts
       cblk => negf%str%cblk
       indblk => negf%str%mat_PL_start

       negf%ni=0; negf%nf=0 
       negf%ni(1)=1 
       negf%nf(1)=2 

       Nstep=NINT((negf%Emax-negf%Emin)/negf%Estep)
       npid = int((Nstep+1)/numprocs)

       !Get out if Emax<Emin and Nstep<0
       if (Nstep.lt.0) then
          if(id0) write(*,*) '0 tunneling points;  current = 0.0'
          return
       endif

       !Extract Contacts in main
       !Tunneling set-up
       do i=1,size(negf%ni)
          if (negf%ni(i).eq.0) then
             size_ni=i-1
             exit
          endif
       enddo

       do i=1,size(negf%nf)
          if (negf%nf(i).eq.0) then
             size_nf=i-1
             exit
          endif
       enddo

       !check size_ni .ne. size_nf
       if (size_ni.ne.size_nf) then 
          size_ni=min(size_ni,size_nf)
          size_nf=min(size_ni,size_nf)
       endif

       call log_allocate(mumin_array,size_ni)
       call log_allocate(mumax_array,size_ni)

       ! find bias window for each contact pair
       do icpl=1,size_ni
          mumin_array(icpl)=min(negf%Efermi(negf%ni(icpl))-negf%mu(negf%ni(icpl)),&
               negf%Efermi(negf%nf(icpl))-negf%mu(negf%nf(icpl)))
          mumax_array(icpl)=max(negf%Efermi(negf%ni(icpl))-negf%mu(negf%ni(icpl)),&
               negf%Efermi(negf%nf(icpl))-negf%mu(negf%nf(icpl)))
       enddo

       ncyc=0
       istart = 1
       iend = npid

       call log_allocate(TUN_MAT,size_ni)
       call log_allocate(TUN_TOT_MAT,Nstep+1,size_ni)   
       !call log_allocate(TUN_PMAT,npid,size_ni,num_channels) 
       !call log_allocate(TUN_TOT_PMAT,Nstep+1,size_ni,num_channels) 
       
       if (do_LEDOS) then
          call log_allocate(LEDOS_MAT,Nstep+1,negf%nLDOS)
          call log_allocate(LEDOS,negf%nLDOS)          
          LEDOS_MAT(:,:)=0.d0
       endif


       !Loop on energy points: tunneling 
       do i1 = istart,iend

          Ec_min = negf%Emin + id*npid*negf%Estep
          Ec = (Ec_min + negf%Estep*(i1-1))*(1.d0,0.d0) !+negf%delta*(0.d0,1.d0) 
          
          call compute_contacts(Ec+negf%delta*(0.d0,1.d0),negf,i1,ncyc,Tlc,Tcl,SelfEneR,GS)

          do icont=1,ncont
             call destroy(Tlc(icont))
             call destroy(Tcl(icont))
          enddo

          if (.not.do_LEDOS) then

             call tunneling_dns(negf%HM,negf%SM,Ec,SelfEneR,negf%ni,negf%nf,size_ni, &
                  negf%str,TUN_MAT)
             TUN_TOT_MAT(i1,:) = TUN_MAT(:) * negf%wght

          else

             LEDOS(:) = 0.d0

             call tun_and_dos(negf%HM,negf%SM,Ec,SelfEneR,Gs,negf%ni,negf%nf,negf%nLDOS, &
                  negf%LDOS,size_ni,negf%str,TUN_MAT,LEDOS)

             TUN_TOT_MAT(i1,:) = TUN_MAT(:) * negf%wght
             LEDOS_MAT(i1,:) = LEDOS(:) * negf%wght

          endif

!          TUN_TOT_MAT(i1,:) = TUN_MAT(:) * negf%wght

!          if (do_LEDOS) then
!             LEDOS_MAT(i1,:) = LEDOS(:) * negf%wght
!          endif


          do icont=1,ncont
             call destroy(SelfEneR(icont))
             call destroy(GS(icont))
          enddo

       enddo !Loop on energy i1 = istart,iend


       !---- SAVE TUNNELING AND DOS ON FILES -----------------------------------------------
       write(ofKP,'(i6.6)') negf%kpoint

       open(121,file='tunneling_'//ofKP//'.dat')
       
       do i = 1,Nstep+1

          Ec=(negf%Emin+negf%Estep*(i-1))*(1,0)

          WRITE(121,'(E17.8,20(E17.8))') REAL(Ec)*negf%eneconv, &
               (TUN_TOT_MAT(i,i1), i1=1,size_ni)

       enddo

       close(121)

       if(do_LEDOS) then

          open(126,file='LEDOS_'//ofKP//'.dat')

          do i = 1,Nstep+1

             Ec=(negf%Emin+negf%Estep*(i-1))*(1,0)

             WRITE(126,'(E17.8,10(E17.8))') REAL(Ec)*negf%eneconv, & 
                  ((LEDOS_MAT(i,iLDOS)/negf%eneconv), iLDOS=1,negf%nLDOS)        

          end do

          close(126)

       endif

       !---------------------------------------------------------------------------
       !   COMPUTATION OF CURRENTS 
       !---------------------------------------------------------------------------
print*,'compute current'  

       open(101,file='current.dat')
  
       call log_allocate(currents,size_ni)
       currents(:)=0.d0
 
       do icpl=1,size_ni
   
          mumin=mumin_array(icpl)
          mumax=mumax_array(icpl)
print*,'Energy window:',mumin,mumax
          if (id0.and.mumin.lt.mumax.and.(mumin.lt.negf%Emin.or.mumax.gt.negf%Emax)) then
             write(*,*) 'WARNING: the interval Emin..Emax is smaller than the bias window'
             write(*,*) 'mumin=',mumin,'mumax=',mumax
             write(*,*) 'emin=',negf%emin,'emax=',negf%emax    
          endif

          currents(icpl)= integrate(TUN_TOT_MAT(:,icpl),mumin,mumax,negf%kbT, &
                                       negf%Emin,negf%Emax,negf%Estep)
       
       enddo

       do i1=1,size_ni
          write(*,'(1x,a,i3,i3,a,ES14.5,a,ES14.5,a)') 'contacts:',negf%ni(i1),negf%nf(i1), &
               ';  current:', currents(i1),' A'
       enddo
 
       close(101)
  
       call log_deallocate(TUN_TOT_MAT)
 
       call log_deallocate(TUN_MAT)
       !call log_deallocate(TUN_PMAT)
       !call log_deallocate(TUN_TOT_PMAT)  
 
       call log_deallocate(mumin_array)
       call log_deallocate(mumax_array)
 
       call log_deallocate(currents)

       if(do_LEDOS) then
          call log_deallocate(LEDOS)
          call log_deallocate(LEDOS_MAT)
       endif

 
     end subroutine compute_current


     !-------------------------------------------------------------------------------


!////////////////////////////////////////////////////////////////////////
!************************************************************************
!
! Function to integrate the tunneling and get the current
!
!************************************************************************

function integrate(TUN_TOT,mumin,mumax,kT,emin,emax,estep)

  implicit none

  real(dp) :: integrate
  real(dp), intent(in) :: mumin,mumax,emin,emax,estep
  real(dp), dimension(:), intent(in) :: TUN_TOT
  real(dp), intent(in) :: kT 

  REAL(dp) :: destep,kbT,TT1,TT2,E3,E4,TT3,TT4
  REAL(dp) :: E1,E2,c1,c2,curr
  INTEGER :: i,i1,N,Nstep,imin,imax
  
  curr=0.d0
  N=0
  destep=1.0d10 
  Nstep=NINT((emax-emin)/estep);

  if (kT.lt.3.d-6) then
    kbT = 1.d-5
  else
    kbT = kT     
  endif
  ! Find initial step for integration

  imin=0
  do i=0,Nstep
     E1=emin+estep*i     
     imin=i-1
     if(E1.ge.mumin-10*kbT) then 
        exit
     endif
  enddo

  ! Find final step for integration 
  imax=0
  do i=Nstep,imin,-1    
     E1=emin+estep*i 
     imax=i+1
     if(E1.le.mumax+10*kbT) then 
        exit
     endif
  enddo

  !checks if the energy interval is appropriate
  if (imin.lt.0.or.imax.gt.Nstep.or.imin.eq.imax) then
     write(*,*) 'WARNING: Wrong energy interval for current calculation'
     write(*,*) 'Suggested interval:',(mumin-10*kbT),(mumax+10*kbT)
  endif

  !rest the min and max to the actual interval
  if (imin.lt.0) imin=0
  if (imax.gt.Nstep) imax=Nstep 

  ! performs the integration with simple trapezium rule. 
  do i=imin,imax-1
     
     E1=emin+estep*i  
     TT1=TUN_TOT(i+1)     
     E2=emin+estep*(i+1)
     TT2=TUN_TOT(i+2)
     
     ! Each step is devided into substeps in order to
     ! smooth out the Fermi function
     do while (destep.ge.2*kbT) 
        N=N+1
        destep=(E2-E1)/N
     enddo
     
     ! within each substep the tunneling is linearly 
     ! interpolated
     do i1=0,N-1
        
        E3=E1+(E2-E1)*i1/N
        E4=E3+(E2-E1)/N
        TT3=( TT2-TT1 )*i1/N + TT1
        TT4=TT3 + (TT2-TT1)/N
        
        c1=2.d0*eovh*(fermi_f(E3,mumax,KbT)-fermi_f(E3,mumin,KbT))*TT3
        c2=2.d0*eovh*(fermi_f(E4,mumax,KbT)-fermi_f(E4,mumin,KbT))*TT4
        
        curr=curr+(c1+c2)*(E4-E3)/2.d0
        
     enddo
     
  enddo
  
  integrate = curr
  
end function integrate

!-----------------------------------------------------------------------
! Contour integration for density matrix 
! DOES INCLUDE FACTOR 2 FOR SPIN !! 
!-----------------------------------------------------------------------
subroutine contour_int(negf)

  type(Tnegf), pointer :: negf 
  Type(z_CSR), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
  type(z_CSR) :: GreenR, TmpMt 

  integer :: npid, istart, iend, NumPoles
  integer :: i, i1, outer, it, ncont, nbl

  real(dp), DIMENSION(:), ALLOCATABLE :: wght,pnts   ! Gauss-quadrature points
  real(dp) :: Omega, Lambda, Rad, Centre
  real(dp) :: muref, kbT, teta, alpha
  real(dp) :: ncyc, dt, Elow

  complex(dp) :: z1,z2,z_diff, zt
  complex(dp) :: Ec, ff, Pc

  it = negf%iteration
  ncont = negf%str%num_conts
  nbl = negf%str%num_PLs
  kbT = negf%kbT

  muref = negf%mu_n

  Omega = negf%n_kt * kbT
  Lambda = 2.d0* negf%n_poles * KbT * pi
  
  if(kbT.eq.0.d0) then
     NumPoles = 0
  else
     NumPoles = negf%n_poles
  endif

  outer = 2 !Compute lower-outer part of density matrix
  !outer = 0  ! no outer part

  call create(TmpMt,negf%H%nrow,negf%H%ncol,negf%H%nrow)
  call initialize(TmpMt)
  ! -----------------------------------------------------------------------
  !  Integration loop starts here
  ! -----------------------------------------------------------------------
  ! ***********************************************************************
  ! 1. INTEGRATION OVER THE CIRCLE Pi..alpha    Np(1)
  ! ***********************************************************************
  ! NEW INTEGRATION FOR COMPLEX DENSITY:
  !----------------------------------------------------
  !   g  [ /           ]     g  [ /         it   ] 
  !  --- [ | Gr(z) dz  ] =  --- [ | iGr(t)Re  dt ]  
  !  2pi [ /           ]    2pi [ /              ]
  !----------------------------------------------------

  Elow = negf%Ec
  Centre = (Lambda**2-Elow**2+(muref-Omega)**2)/(2.d0*(muref-Omega-Elow))
  Rad = Centre - Elow

  if (negf%kbT.ne.0.d0) then        
     alpha = atan(Lambda/(muref-Centre-Omega)) 
  else
     alpha = 0.1d0*pi             
  end if

  !Setting weights for gaussian integration 
  allocate(wght(negf%Np_n(1)))
  allocate(pnts(negf%Np_n(1)))

  call gauleg(pi,alpha,pnts,wght,negf%Np_n(1))

  !Computing complex integral (Common for T>=0)

  npid = int(negf%Np_n(1)/numprocs)
  istart = id*npid+1
  if(id.ne.(numprocs-1)) then 
     iend = (id+1)*npid
  else
     iend = negf%Np_n(1)
  end if

  do i = istart,iend

     if (negf%verbose.gt.VBT) then
        write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 1: point #',i,'/',iend,'  CPU=&
             &', id
     endif

     if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Green`s funct ')

     Pc = Rad*exp(j*pnts(i))
     Ec = Centre+Pc
     zt = j * Pc * negf%spin * wght(i)/(2.d0*pi)

     call compute_contacts(Ec,negf,i,ncyc,Tlc,Tcl,SelfEneR,GS)

     call calls_eq_mem_dns(negf%HM,negf%SM,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer)

     if(negf%DorE.eq.'D') then
        call concat(TmpMt,zt,GreenR,1,1)
     endif
     if(negf%DorE.eq.'E') then
        call concat(TmpMt,zt*Ec,GreenR,1,1)
     endif

     call destroy(GreenR)

     do i1=1,ncont
        call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
     enddo

     if (id0.and.negf%verbose.gt.VBT) call write_clock

  enddo

  deallocate(wght)
  deallocate(pnts)

  ! *******************************************************************************
  ! 2. INTEGRATION OVER THE SEGMENT [Ec + j*Lambda, mu(r) + Omega+j*Lambda]
  ! (Temp /= 0) OR OVER THE CIRCLE WITH TETA FROM ZERO TO ALPHA (Temp == 0)
  ! *******************************************************************************
  ! NEW INTEGRATION FOR COMPLEX DENSITY (T=0):
  !----------------------------------------------------
  !   2  [ /           ]     1  [ /         it   ] 
  !  --- [ |  Gr(z) dz ] =   -- [ | Gr(t)Rie  dt ]  
  !  2pi [ /           ]     pi [ /              ]
  !----------------------------------------------------
  ! NEW INTEGRATION FOR COMPLEX DENSITY (T>0):
  !----------------------------------------------------
  !   g  [ /                ]      g   [ /                       ] 
  !  --- [ |  Gr(z)*f(z) dz ] =   ---  [ | Gr(z)*(z2-z1)*f(z)*dt ]  
  !  2pi [ /                ]     2*pi [ /                       ]
  !----------------------------------------------------

  allocate(wght(negf%Np_n(2)))
  allocate(pnts(negf%Np_n(2)))

  if (negf%kbT.eq.0.d0) then                        ! Circle integration T=0

     call  gauleg(alpha,0.d0,pnts,wght,negf%Np_n(2))
     
  else                                          ! Segment integration T>0
     
     z1 = muref + Omega + j*Lambda
     z2 = muref - Omega + j*Lambda
     
     z_diff = z2 - z1
     
     call  gauleg(1.d0,0.d0,pnts,wght,negf%Np_n(2))    !Setting weights for integration
     
  endif

  npid = int(negf%Np_n(2)/numprocs)

  istart = id*npid+1
  if(id.ne.(numprocs-1)) then 
     iend = (id+1)*npid
  else
     iend = negf%Np_n(2)
  end if

  do i = istart,iend

     if (negf%verbose.gt.VBT) then
        write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 2: point #',i,'/',iend,'  CPU=&
             &', id
     endif

     if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Green`s funct ')

     if (negf%kbT.eq.0.d0) then                      ! Circle integration T=0            
        
        Pc = Rad*exp(j*pnts(i))
        Ec = Centre+Pc
        dt = negf%spin*wght(i)/(2.d0*pi)
        zt = dt*Pc*j
       
     else                                        ! Segment integration T>0
        
        Ec = z1 + pnts(i)*z_diff
        ff = fermi_fc(Ec,muref,KbT)
        zt = negf%spin * z_diff * ff * wght(i) / (2.d0 *pi)

     endif

     call compute_contacts(Ec,negf,negf%Np_n(1)+i,ncyc,Tlc,Tcl,SelfEneR,GS)

     call calls_eq_mem_dns(negf%HM,negf%SM,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer) 

     if(negf%DorE.eq.'D') then
        call concat(TmpMt,zt,GreenR,1,1)
     endif
     if(negf%DorE.eq.'E') then
        call concat(TmpMt,zt*Ec,GreenR,1,1)
     endif

     call destroy(GreenR)

     do i1=1,ncont
        call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
     enddo

     if (id0.and.negf%verbose.gt.VBT) call write_clock

  enddo

  deallocate(wght)
  deallocate(pnts)

  ! *******************************************************************************
  ! 3. SUMMATION OVER THE POLES ENCLOSED IN THE CONTOUR  (NumPoles)
  ! *******************************************************************************          
  ! NEW INTEGRATION FOR COMPLEX DENSITY (T>=0):
  !---------------------------------------------------------------------
  !             [  g                 ]          
  !  2 pi j* Res[ --- *Gr(z_k)f(z_k)  ] =  j*(-kb T)* Gr(z_k)     
  !             [ 2pi                ]         
  !                                         (-kb*T) <- Residue
  !---------------------------------------------------------------------

  npid = int(NumPoles/numprocs)
  istart = id*npid+1
  if(id.ne.(numprocs-1)) then 
     iend = (id+1)*npid
  else
     iend = NumPoles
  end if

  do i = istart,iend

     if (negf%verbose.gt.VBT) then
        write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 3: point #',i,'/',iend,'  CPU=&
             &', id
     endif

     if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Green`s funct ')

     Ec = muref + j * KbT *pi* (2.d0*real(i,dp) - 1.d0)   

     zt= -j*negf%spin*KbT

     call compute_contacts(Ec,negf,negf%Np_n(1)+negf%Np_n(2)+i,ncyc,Tlc,Tcl,SelfEneR,GS)

     call calls_eq_mem_dns(negf%HM,negf%SM,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer)

     if(negf%DorE.eq.'D') then
        call concat(TmpMt,zt,GreenR,1,1)
     endif
     if(negf%DorE.eq.'E') then
        call concat(TmpMt,zt*Ec,GreenR,1,1)
     endif

     call destroy(GreenR)   

     do i1=1,ncont
        call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
     enddo

     if (id0.and.negf%verbose.gt.VBT) call write_clock

  enddo

  if(negf%DorE.eq.'D') then
     call zspectral(TmpMt,TmpMt,0,negf%rho)
  endif
  if(negf%DorE.eq.'E') then
     call zspectral(TmpMt,TmpMt,0,negf%rho_eps)
  endif


  call destroy(TmpMt)


end subroutine contour_int

!-----------------------------------------------------------------------
! Contour integration for density matrix 
! DOES INCLUDE FACTOR 2 FOR SPIN !! 
!-----------------------------------------------------------------------
  subroutine contour_int_n(negf)

    type(Tnegf), pointer :: negf 
    Type(z_CSR), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    type(z_CSR) :: GreenR, TmpMt 

    integer :: npid, istart, iend, NumPoles
    integer :: i, i1, outer, it, ncont, nbl

    real(dp), DIMENSION(:), ALLOCATABLE :: wght,pnts   ! Gauss-quadrature points
    real(dp) :: Omega, Lambda
    real(dp) :: muref, kbT
    real(dp) :: ncyc

    complex(dp) :: z1,z2,z_diff, zt
    complex(dp) :: Ec, ff

    it = negf%iteration
    ncont = negf%str%num_conts
    nbl = negf%str%num_PLs
    kbT = negf%kbT
    
    muref = negf%mu_n

    Omega = negf%n_kt * kbT
    Lambda = 2.d0* negf%n_poles * KbT * pi
    NumPoles = negf%n_poles

    outer = 1 !no contacts no outer

    call create(TmpMt,negf%H%nrow,negf%H%ncol,negf%H%nrow)
    call initialize(TmpMt)


    ! *******************************************************************************
    ! 1. INTEGRATION OVER THE SEGMENT [Ec - dEc , Ec - dEc + j*Lambda]
    ! *******************************************************************************
    ! NEW INTEGRATION FOR COMPLEX DENSITY (T>0):
    !----------------------------------------------------
    !   g  [ /                ]      g   [ /                       ] 
    !  --- [ |  Gr(z)*f(z) dz ] =   ---  [ | Gr(z)*(z2-z1)*f(z)*dt ]  
    !  2pi [ /                ]     2*pi [ /                       ]
    !----------------------------------------------------

    z1 = negf%Ec + negf%DeltaEc
    z2 = negf%Ec + negf%DeltaEc + j*Lambda

    z_diff = z2 - z1

    allocate(wght(negf%Np_n(1)))
    allocate(pnts(negf%Np_n(1)))

    call gauleg(0.d0,1.d0,pnts,wght,negf%Np_n(1))

    npid = int(negf%Np_n(1)/numprocs)
    istart = id*npid+1
    if(id.ne.(numprocs-1)) then 
       iend = (id+1)*npid
    else
       iend = negf%Np_n(1)
    end if

    do i = istart,iend

       !if (negf%verbose.gt.80) then
       !   write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 1: point #',i,'/',iend,'  CPU=&
       !        &', id
       !endif

       Ec = z1 + pnts(i) * z_diff

       ff = fermi_fc(Ec,muref,KbT)

       zt = negf%spin * z_diff * ff * wght(i) / (2.d0 *pi)

       call compute_contacts(Ec,negf,i,ncyc,Tlc,Tcl,SelfEneR,GS)

       call calls_eq_mem_dns(negf%HM,negf%SM,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer)

       if(negf%DorE.eq.'D') then
          call concat(TmpMt,zt,GreenR,1,1)
       endif
       if(negf%DorE.eq.'E') then
          call concat(TmpMt,zt*Ec,GreenR,1,1)
       endif

       call destroy(GreenR) 

       do i1=1,ncont
          call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
       enddo

    enddo

    deallocate(wght)
    deallocate(pnts)


    ! *******************************************************************************
    ! 2. INTEGRATION OVER THE SEGMENT [Ec + j*Lambda, mu(r) + Omega+j*Lambda]
    ! *******************************************************************************
    ! NEW INTEGRATION FOR COMPLEX DENSITY (T>0):
    !----------------------------------------------------
    !   g  [ /                ]      g   [ /                       ] 
    !  --- [ |  Gr(z)*f(z) dz ] =   ---  [ | Gr(z)*(z2-z1)*f(z)*dt ]  
    !  2pi [ /                ]     2*pi [ /                       ]
    !----------------------------------------------------


    allocate(wght(negf%Np_n(2)))
    allocate(pnts(negf%Np_n(2)))

    call gauleg(0.d0,1.d0,pnts,wght,negf%Np_n(2))    !Setting weights for integration

    z1 = negf%Ec + negf%DeltaEc  + j*Lambda
    z2 = muref + Omega + j*Lambda

    z_diff = z2 - z1

    npid = int(negf%Np_n(2)/numprocs)

    istart = id*npid+1
    if(id.ne.(numprocs-1)) then 
       iend = (id+1)*npid
    else
       iend = negf%Np_n(2)
    end if
  
    do i = istart,iend

       !if (negf%verbose.gt.80) then
       !   write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 2: point #',i,'/',iend,'  CPU=&
       !        &', id
       !endif

       Ec = z1 + pnts(i) * z_diff

       ff = fermi_fc(Ec,muref,KbT)

       zt = negf%spin *  z_diff * ff * wght(i) / (2.d0 *pi)

       call compute_contacts(Ec,negf,negf%Np_n(1)+i,ncyc,Tlc,Tcl,SelfEneR,GS)

       call calls_eq_mem_dns(negf%HM,negf%SM,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer) 

       if(negf%DorE.eq.'D') then
          call concat(TmpMt,zt,GreenR,1,1)
       endif
       if(negf%DorE.eq.'E') then
          call concat(TmpMt,zt*Ec,GreenR,1,1)
       endif

       call destroy(GreenR)

       do i1=1,ncont
          call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
       enddo

    enddo

    deallocate(wght)
    deallocate(pnts)



    ! *******************************************************************************
    ! 3. SUMMATION OVER THE POLES ENCLOSED IN THE CONTOUR  (NumPoles)
    ! *******************************************************************************          
    ! NEW INTEGRATION FOR COMPLEX DENSITY (T>=0):
    !---------------------------------------------------------------------
    !                  [ 1                  ]          
    ! 2 pi * g * j* Res[ -- *Gr(z_k)f(z_k)  ] =  j*(-kb T)* Gr(z_k)     
    !                  [ 2pi                ]         
    !                                              (-kb*T) <- Residue
    !---------------------------------------------------------------------

    npid = int(NumPoles/numprocs)
    istart = id*npid+1
    if(id.ne.(numprocs-1)) then 
       iend = (id+1)*npid
    else
       iend = NumPoles
    end if

    do i = istart,iend

       !if (negf%verbose.gt.80) then
       !   write(6,'(a17,i3,a1,i3,a6,i3)') 'POLES: point #',i,'/',iend,'  CPU=', id
       !endif

       Ec = muref + j * KbT *pi* (2.d0*real(i,dp) - 1.d0)   

       zt= -j * KbT * negf%spin *(1.d0,0.d0) 

       call compute_contacts(Ec,negf,negf%Np_n(1)+negf%Np_n(2)+i,ncyc,Tlc,Tcl,SelfEneR,GS)

       call calls_eq_mem_dns(negf%HM,negf%SM,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer)

       if(negf%DorE.eq.'D') then
          call concat(TmpMt,zt,GreenR,1,1)
       endif
       if(negf%DorE.eq.'E') then
          call concat(TmpMt,zt*Ec,GreenR,1,1)
       endif

       call destroy(GreenR)   

       do i1=1,ncont
          call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
       enddo

    enddo

    if(negf%DorE.eq.'D') then
       call zspectral(TmpMt,TmpMt,0,negf%rho)
    endif
    if(negf%DorE.eq.'E') then
       call zspectral(TmpMt,TmpMt,0,negf%rho_eps)
    endif


    call destroy(TmpMt)

  end subroutine contour_int_n


!--------------------------------------------------------------------------------
!-----------------------------------------------------------------------
! Contour integration for density matrix 
! DOES INCLUDE FACTOR 2 FOR SPIN !! 
!-----------------------------------------------------------------------
  subroutine contour_int_p(negf)

    type(Tnegf), pointer :: negf 
    Type(z_CSR), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    type(z_CSR) :: GreenR, TmpMt 

    integer :: npid, istart, iend, NumPoles
    integer :: i, i1, outer, it, ncont, nbl

    real(dp), DIMENSION(:), ALLOCATABLE :: wght,pnts   ! Gauss-quadrature points
    real(dp) :: Omega, Lambda
    real(dp) :: muref, kbT
    real(dp) :: ncyc

    complex(dp) :: z1,z2,z_diff, zt
    complex(dp) :: Ev, ff

    it = negf%iteration
    ncont = negf%str%num_conts
    nbl = negf%str%num_PLs
    kbT = negf%kbT

    muref = negf%mu_p

    Omega = negf%n_kt * kbT
    Lambda = 2.d0* negf%n_poles * KbT * pi
    NumPoles = negf%n_poles
    
    outer = 1 !no contacts no outer

    call create(TmpMt,negf%H%nrow,negf%H%ncol,negf%H%nrow)
    call initialize(TmpMt)

    !1. INTEGRATION OVER THE SEGMENT

    z1 = negf%Ev + negf%DeltaEv + j*Lambda
    z2 = negf%Ev + negf%DeltaEv

    z_diff = z2 - z1

    allocate(wght(negf%Np_p(1)))
    allocate(pnts(negf%Np_p(1)))

    call gauleg(0.d0,1.d0,pnts,wght,negf%Np_p(1))

    npid = int(negf%Np_p(1)/numprocs)
    istart = id*npid+1
    if(id.ne.(numprocs-1)) then 
       iend = (id+1)*npid
    else
       iend = negf%Np_p(1)
    end if


    do i = istart,iend

       !if (negf%verbose.gt.80) then
       !   write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 1: point #',i,'/',iend,'  CPU=&
       !        &', id
       !endif

       Ev = z1 + pnts(i) * z_diff

       ff = (1.d0,0.d0) - fermi_fc(Ev,muref,KbT)

       zt = negf%spin * z_diff * ff * wght(i) / (2.d0 *pi)

       call compute_contacts(Ev,negf,i,ncyc,Tlc,Tcl,SelfEneR,GS)

       call calls_eq_mem_dns(negf%HM,negf%SM,Ev,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer)

       if(negf%DorE.eq.'D') then
          call concat(TmpMt,zt,GreenR,1,1)
       endif
       if(negf%DorE.eq.'E') then
          call concat(TmpMt,zt*Ev,GreenR,1,1)
       endif

       call destroy(GreenR) 

       do i1=1,ncont
          call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
       enddo

    enddo

    deallocate(wght)
    deallocate(pnts)

    ! 2. INTEGRATION OVER THE SEGMENT 

    allocate(wght(negf%Np_p(2)))
    allocate(pnts(negf%Np_p(2)))

    call gauleg(0.d0,1.d0,pnts,wght,negf%Np_p(2))    !Setting weights for integration

    z1 = muref - Omega + j*Lambda
    z2 = negf%Ev + negf%DeltaEv  + j*Lambda

    z_diff = z2 - z1

    npid = int(negf%Np_p(2)/numprocs)

    istart = id*npid+1
    if(id.ne.(numprocs-1)) then 
       iend = (id+1)*npid
    else
       iend = negf%Np_p(2)
    end if
  
    do i = istart,iend

       !if (negf%verbose.gt.80) then
       !   write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 2: point #',i,'/',iend,'  CPU=&
       !        &', id
       !endif

       Ev = z1 + pnts(i) * z_diff

       ff = (1.d0,0.d0) - fermi_fc(Ev,muref,KbT)

       zt = z_diff * negf%spin * ff * wght(i) / (2.d0 *pi)

       call compute_contacts(Ev,negf,negf%Np_p(1)+i,ncyc,Tlc,Tcl,SelfEneR,GS)

       call calls_eq_mem_dns(negf%HM,negf%SM,Ev,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer) 

       if(negf%DorE.eq.'D') then
          call concat(TmpMt,zt,GreenR,1,1)
       endif
       if(negf%DorE.eq.'E') then
          call concat(TmpMt,zt*Ev,GreenR,1,1)
       endif

       call destroy(GreenR)

       do i1=1,ncont
          call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
       enddo

    enddo

    deallocate(wght)
    deallocate(pnts)
    
    ! SUMMATION OVER THE POLES ENCLOSED IN THE CONTOUR
    
    npid = int(NumPoles/numprocs)
    istart = id*npid+1
    if(id.ne.(numprocs-1)) then 
       iend = (id+1)*npid
    else
       iend = NumPoles
    end if

    do i = istart,iend

       !if (negf%verbose.gt.80) then
       !   write(6,'(a17,i3,a1,i3,a6,i3)') 'POLES: point #',i,'/',iend,'  CPU=', id
       !endif

       Ev =  muref + j * KbT *pi* (2.d0*real(i,dp) - 1.d0)   

       zt= j*negf%spin*KbT*(1.d0,0.d0) 

       call compute_contacts(Ev,negf,negf%Np_p(1)+negf%Np_p(2)+i,ncyc,Tlc,Tcl,SelfEneR,GS)

       call calls_eq_mem_dns(negf%HM,negf%SM,Ev,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer)

       if(negf%DorE.eq.'D') then
          call concat(TmpMt,zt,GreenR,1,1)
       endif
       if(negf%DorE.eq.'E') then
          call concat(TmpMt,zt*Ev,GreenR,1,1)
       endif

       call destroy(GreenR)  

        do i1=1,ncont
          call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
       enddo

    enddo

    if(negf%DorE.eq.'D') then
       call zspectral(TmpMt,TmpMt,0,negf%rho)
    endif
    if(negf%DorE.eq.'E') then
       call zspectral(TmpMt,TmpMt,0,negf%rho_eps)
    endif

    call destroy(TmpMt)

  end subroutine contour_int_p

!--------------------------------------------!
!--------------------------------------------!
! Non equilibrium integration over real axis !
!--------------------------------------------!
!--------------------------------------------!
!-----------------------------------------------------------------------
! Contour integration for density matrix 
! DOES INCLUDE FACTOR 2 FOR SPIN !! 
!-----------------------------------------------------------------------
  subroutine real_axis_int(negf)

    type(Tnegf), pointer :: negf 
    Type(z_CSR), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS
    type(z_CSR) :: GreenR, TmpMt 

    integer :: npid, istart, iend, NumPoles, min, nc_vec(1)
    integer :: i, i1, outer, it, ncont, nbl, j1, npT

    real(dp), DIMENSION(:), allocatable :: wght,pnts   ! Gauss-quadrature points
    real(dp), DIMENSION(:), allocatable :: frm_f
    real(dp) :: Omega, Lambda
    real(dp) :: mumin, mumax
    real(dp) :: ncyc, kbT, dt

    complex(dp) :: z1,z2,z_diff, zt
    complex(dp) :: Ec, ff

    it = negf%iteration
    ncont = negf%str%num_conts
    nbl = negf%str%num_PLs
    kbT = negf%kbT

    mumin=minval(negf%Efermi(1:ncont)-negf%mu(1:ncont))
    mumax=maxval(negf%Efermi(1:ncont)-negf%mu(1:ncont))

    if (mumax.gt.mumin) then

       nc_vec = maxloc(negf%Efermi(1:ncont)-negf%mu(1:ncont))
       min = nc_vec(1)
       
       Omega = negf%n_kt * kbT
       Lambda = 2.d0 * negf%n_poles * KbT * pi
       NumPoles = negf%n_poles
       
       outer = 2 !no contacts no outer
       
       call log_allocate(frm_f,ncont)
       
       call create(TmpMt,negf%H%nrow,negf%H%ncol,negf%H%nrow)
       call initialize(TmpMt)
       
       !Compute extended number of points due to kT
       npT=nint(negf%Np_real/(mumax-mumin))*Omega
 
       allocate(pnts(negf%Np_real+2*npT))
       allocate(wght(negf%Np_real+2*npT))

       !Setting weights for gaussian integration
       call gauleg(mumin-Omega,mumax+Omega,pnts,wght,negf%Np_real+2*npT)

       !Computing real axis integral       
       npid = int((negf%Np_real+2*npT)/numprocs)
       istart = id*npid+1
       if(id.ne.(numprocs-1)) then 
          iend = (id+1)*npid
       else
          iend = negf%Np_real+2*npT
       end if
       
       do i = istart,iend

          if (negf%verbose.gt.VBT) then
             write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL neq: point #',i,'/',iend,'  CPU=&
                  &', id
          endif

          if (id0.and.negf%verbose.gt.VBT) call message_clock('Compute Green`s funct ')

          Ec = cmplx(pnts(i),negf%delta,dp)

          dt = negf%spin * wght(i)/pi
          zt=dt*(1.d0,0.d0)

          do j1=1,ncont
             frm_f(j1)=fermi_f(real(Ec),negf%Efermi(j1)-negf%mu(j1),KbT)
          enddo

          call compute_contacts(Ec,negf,i,ncyc,Tlc,Tcl,SelfEneR,GS)

          call calls_neq_mem_dns(negf%HM,negf%SM,Ec,SelfEneR,Tlc,Tcl,GS,negf%str,frm_f,min,GreenR,outer)

          do i1=1,ncont
             call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
          enddo

          if(negf%DorE.eq.'D') then
             call concat(TmpMt,zt,GreenR,1,1)
          endif
          if(negf%DorE.eq.'E') then
             call concat(TmpMt,zt*Ec,GreenR,1,1)
          endif

          call destroy(GreenR) 

          if (id0.and.negf%verbose.gt.VBT) call write_clock

       enddo

       deallocate(wght)
       deallocate(pnts)

       if(negf%DorE.eq.'D') then
          if(allocated(negf%rho%nzval)) then
             call concat(negf%rho,TmpMt,1,1)
          else
             call clone(TmpMt,negf%rho) 
          endif
       endif
       if(negf%DorE.eq.'E') then
          if(allocated(negf%rho_eps%nzval)) then
             call concat(negf%rho_eps,GreenR,1,1)
          else
             call clone(TmpMt,negf%rho_eps) 
          endif
       endif
       
       call destroy(TmpMt)
       
       call log_deallocate(frm_f)

    end if


  end subroutine real_axis_int

!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------

  subroutine gauleg(x1,x2,x,w,n)

    real(kind=dp), PARAMETER :: ACC = 1d-15

    INTEGER n
    real(kind=dp) :: x1,x2,x(n),w(n)

    INTEGER i,k,m
    real(kind=dp) :: p1,p2,p3,pp,xl,xm,z,z1

    m=(n+1)/2

    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)

    do i=1,m

       z=cos(Pi*(i-0.25d0)/(n+0.5d0))

       do
          p1=1.d0
          p2=0.d0

          ! Legendre polynomial p1 evaluated by rec. relations:
          do k=1,n
             p3=p2
             p2=p1
             p1=((2.d0*k-1.d0)*z*p2-(k-1.d0)*p3)/k
          enddo
          ! Derivative pp using the relation of p1 and p2:
          pp=n*(z*p1-p2)/(z*z-1.d0)

          ! Newton method to refine the zeros:
          z1=z
          z=z1-p1/pp

          if(abs(z-z1).le.ACC) exit
       enddo

       ! Scale the interval to x1..x2:
       x(i)=xm-xl*z
       x(n+1-i)=xm+xl*z
       w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
       w(n+1-i)=w(i)
    enddo
    
    return
    
  end subroutine gauleg


!--------------------------------------------------------------------  
  
  subroutine reorder(mat)
    type(z_CSR) :: mat
    

    type(z_CSR) :: P, Tmp

    integer, dimension(:), allocatable :: perm
    integer :: i, nrow

    nrow=mat%nrow

    call log_allocate(perm,nrow)

    call genrcm(nrow, mat%nnz, mat%rowpnt, mat%colind, perm)

    call create(P,nrow,nrow,nrow)

    do i=1,nrow
       P%nzval(i)=1
       P%colind(i)=cmplx(perm(i),dp)
       P%rowpnt(i)=i
    enddo
    P%rowpnt(nrow+1)=nrow+1

    
    call create(Tmp,nrow,nrow,mat%nnz)

    call zamub_st(P,mat,Tmp)

    call ztransp_st(P)

    call zamub_st(Tmp,P,mat)   

    call destroy(P,Tmp)

    call log_deallocate(perm)
 
  end subroutine reorder

!----------------------------------------------------------------------------

  subroutine partition(mat,nbl,part)
    type(z_CSR) :: mat
    
    integer :: nbl
    integer :: n
    
    integer :: volume
    integer :: numflag
    integer :: wghts
    integer :: wgtflag

    integer, dimension(:) :: part


    integer, dimension(:), allocatable :: options
    integer, dimension(:), allocatable :: vwgt  
    integer, dimension(:), allocatable :: vsize
    

    external METIS_PartGraphVKway

    numflag = 1
    wghts = 0
    wgtflag = 0
    n = mat%nrow

    call log_allocate(vwgt, 0)
    call log_allocate(vsize, 0)
    call log_allocate(options, 5)
    options(1) = 0


    !call METIS_PartGraphVKway(n, mat%rowpnt, mat%colind, vwgt, vsize, wgtflag, &
    !                          numflag, nbl, options, volume, part)

    call METIS_PartGraphKway(n, mat%rowpnt, mat%colind, vwgt, vsize, wgtflag, &
                              numflag, nbl, options, volume, part)    

    call log_deallocate(vwgt)
    call log_deallocate(vsize)
    call log_deallocate(options)
    

  end subroutine partition

!----------------------------------------------------------------------

  
  subroutine partition2(mat,nbl,blks)
    type(z_CSR), intent(in) :: mat
    integer, intent(out) :: nbl
    integer, dimension(:), intent(inout) :: blks 

    integer :: j, k
    integer :: i1, i2

    integer :: rn, rnold, tmax, rmax,  nrow, maxmax
    integer :: dbuff, minsize

    nrow = mat%nrow
    
    ! Find maximal stancil of the matrix and on which row
    !  ( Xx     )
    !  ( xXxx   )  
    !  (  xXxxx )  <- maxmax = 3 ; rmax = 3
    !  (   xXx  )
    maxmax = 0
    do j=1,nrow

       i1 = mat%rowpnt(j)
       i2 = mat%rowpnt(j+1) - 1
       tmax = maxval(mat%colind(i1:i2)) - j

       if(tmax .gt. maxmax) then 
          maxmax = tmax
          rmax = j
       endif

       dbuff = maxmax   ! dbuff should be linked to maxmax
       minsize = dbuff  ! minsize bisogna identificarlo meglio. 
    enddo

    ! Define central block 
    rn = rmax - maxmax/2 - dbuff 

    if(rn-dbuff.ge.0) then 

       blks(1) = rn-1  ! fine del blocco precedente

       nbl = 1

       do 
          
          do j = rn, dbuff, -1

             rnold = rn
             i1 = mat%rowpnt(j-dbuff+1)
             i2 = mat%rowpnt(j+1) - 1
             k = maxval(mat%colind(i1:i2))
       
             if(k.lt.rn) then
                rn = j       
                nbl = nbl + 1
                blks(nbl) = j-1 ! fine del blocco precedente
                exit
             endif
          enddo

          if(rn.le.minsize .or. rnold.eq.rn) then
             exit
          endif
          
       enddo

       rn = rmax - maxmax/2 - dbuff

    else
       nbl= 0
       rn = 1

    endif

    do 

       do j = rn, nrow-dbuff+1, 1

          rnold = rn
          i1 = mat%rowpnt(j)
          i2 = mat%rowpnt(j+dbuff) - 1
          k = minval(mat%colind(i1:i2)) 

          if(k.gt.rn) then  
             rn = j
             nbl = nbl + 1 
             blks(nbl) = j-1 ! fine del blocco 
             exit
          endif
       enddo

       if(nrow-rn.le.minsize .or. rnold.eq.rn) then
          exit
       endif
    enddo

    nbl = nbl + 1
    blks(nbl) = nrow

    print*, 'done. blks:',blks(1:nbl)


  end subroutine partition2

  !----------------------------------------------------------


  Subroutine sort(blks, Ipt)
    ! *
    ! ***********************************
    ! * Sort Array X(:) in ascendent order.
! * If present Ipt, a pointer with the 
! * changes is returned in Ipt. 
! ***********************************
 
    Integer, Intent (inout) :: blks(:)
    Integer, Intent (out), Optional :: Ipt(:)
 
    Integer :: Rtmp
    Integer :: i, j
 
    If (Present(Ipt)) Then
       Forall (i=1:Size(blks)) Ipt(i) = i
 
       Do i = 2, Size(blks)
          Rtmp = blks(i)
          Do j = i-1, 1, -1
             If (Rtmp < blks(j)) Then
                blks(j+1) = blks(j)
                call Swap(blks, j, j+1)
             Else
                Exit
             End If
          End Do
          blks(j+1) = Rtmp
       End Do
    Else
       Do i = 2, Size(blks)
          Rtmp = blks(i)
          Do j = i-1, 1, -1
             If (Rtmp < blks(j)) Then
                blks(j+1) = blks(j)
             Else
                Exit
             End If
          End Do
          blks(j+1) = Rtmp
       End Do
    End If
 
    Return
  End Subroutine sort
 
! ***********************************
! *
  Subroutine Swap(X, i, j)
! *
! ***********************************
! * Swaps elements I and J of array X(:). 
! ***********************************
 
    Integer, Intent (inout) :: X(:)
    Integer, Intent (in) :: i, j
 
    Integer :: Itmp
 
    Itmp = X(I)
    X(I) = X(J)
    X(J) = Itmp
 
    Return
  End Subroutine Swap
         
!-----------------------------------------------------       



end module libnegf
