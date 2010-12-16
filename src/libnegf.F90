module libnegf

 use precision
 use fermi_dist
 use lib_param
 use constants
 use mpi_globals, only : id, numprocs, id0
 use input_output
 use allocation
 use structure
 use sparsekit_drv
 use inversions
 use iterative
 use iterative_dns
 use mat_def
 use rcm_module
 use clock
 use extract
 use contselfenergy

 implicit none
 private

 public :: init_negf, destroy_negf
 public :: compute_dos, contour_int_n, contour_int_p, compute_contacts
 public :: reorder, partition, partition2
 public :: sort, swap

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

       call read_H(101,102,negf%S,negf%form)

       close(101)
       close(102)

    else
       ! create an Id matrix for S
       call create_id(negf%S,negf%H%nrow) 
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
    read(101,*) negf%Temp
    read(101,*) negf%Np_n(1)
    read(101,*) negf%Np_n(2)
    read(101,*) negf%Np_p(1)
    read(101,*) negf%Np_p(2)
    read(101,*) negf%n_kt
    read(101,*) negf%n_poles
    read(101,*) negf%spin
    read(101,*) negf%delta

    close(101)

    print*, 'before create T'
    call create_Tstruct(ncont, nbl, PL_end, cont_end, surf_end, cblk, negf%str)
    print*, 'done'

    call log_deallocate(PL_end)
    call log_deallocate(cblk)
    call log_deallocate(cont_end)
    call log_deallocate(surf_end)

  end subroutine init_negf

!--------------------------------------------------------------------

  subroutine destroy_negf(negf)
    type(Tnegf), pointer :: negf   
    integer :: i

    if (allocated(negf%H%nzval)) call destroy(negf%H) 
    if (allocated(negf%S%nzval)) call destroy(negf%S)

    if (allocated(negf%HM%nzval)) call destroy(negf%HM)
    if (allocated(negf%SM%nzval)) call destroy(negf%SM)

    do i=1,negf%str%num_conts
       if (allocated(negf%HC(i)%val)) call destroy(negf%HC(i))
       if (allocated(negf%SC(i)%val)) call destroy(negf%SC(i))
       if (allocated(negf%HMC(i)%nzval)) call destroy(negf%HMC(i))
       if (allocated(negf%SMC(i)%nzval)) call destroy(negf%SMC(i))
    enddo

    if (allocated(negf%rho%nzval)) call destroy(negf%rho)
    if (allocated(negf%rho_eps%nzval)) call destroy(negf%rho_eps)

    call kill_Tstruct(negf%str)   

  end subroutine destroy_negf

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
       
       if (id0.and.pnegf%verbose.gt.60) call message_clock('Compute Green`s funct ')
       call SelfEnergies(Ec,ncont,GS,Tlc,Tcl,SelfEneR)

     end subroutine compute_contacts

!-------------------------------------------------------------------------------

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

    open(101, file='dos.dat', form='formatted')

    it = negf%iteration
    ncont = negf%str%num_conts
    nbl = negf%str%num_PLs

    muref = negf%mu_n

    kbT = Kb * negf%Temp * negf%eneconv

    Omega = negf%n_kt * kbT
    Lambda = 2.d0* negf%n_poles * KbT * pi
    NumPoles = negf%n_poles

    outer = 1 !no contacts no outer

    call create(TmpMt,negf%H%nrow,negf%H%ncol,negf%H%nrow)
    call init(TmpMt)

    ! *******************************************************************************
    ! 1. INTEGRATION OVER THE SEGMENT [Ec - dEc , Ec - dEc + j*Lambda]
    ! *******************************************************************************
    ! NEW INTEGRATION FOR COMPLEX DENSITY (T>0):
    !----------------------------------------------------
    !   1  [ /                ]      1   [ /                       ] 
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

       if (negf%verbose.gt.80) then
          write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 1: point #',i,'/',iend,'  CPU=&
               &', id
       endif

       Ec = z1 + pnts(i) * z_diff

       ff = fermi_fc(Ec,muref,KbT)

       zt = z_diff * ff * wght(i) / (2.d0 *pi)

       call compute_contacts(Ec,negf,i,ncyc,Tlc,Tcl,SelfEneR,GS)

       call calls_eq_mem_dns(negf%HM,negf%SM,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer)

       call concat(TmpMt,zt,GreenR,1,1) 

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
    !   g  [ /                ]     g    [ /                       ] 
    !- --- [ |  Gr(z)*f(z) dz ] = - ---  [ | Gr(z)*(z2-z1)*f(z)*dt ]  
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

       if (negf%verbose.gt.80) then
          write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 2: point #',i,'/',iend,'  CPU=&
               &', id
       endif

       Ec = z1 + pnts(i) * z_diff

       ff = fermi_fc(Ec,muref,KbT)

       zt = z_diff * ff * wght(i) / (2.d0 *pi)

       call compute_contacts(Ec,negf,negf%Np_n(1)+i,ncyc,Tlc,Tcl,SelfEneR,GS)

       call calls_eq_mem_dns(negf%HM,negf%SM,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer) 

       call concat(TmpMt,zt,GreenR,1,1)  !TmpMt=TmpMt+GreenR

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
    !             [ 1                  ]          
    !  2 pi j* Res[ -- *Gr(z_k)f(z_k)  ] =  j*(-kb T)* Gr(z_k)     
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

       if (negf%verbose.gt.80) then
          write(6,'(a17,i3,a1,i3,a6,i3)') 'POLES: point #',i,'/',iend,'  CPU=', id
       endif

       Ec = muref + j * KbT *pi* (2.d0*real(i,dp) - 1.d0)   

       zt= -j*KbT*(1.d0,0.d0) 

       call compute_contacts(Ec,negf,negf%Np_n(1)+negf%Np_n(2)+i,ncyc,Tlc,Tcl,SelfEneR,GS)

       call calls_eq_mem_dns(negf%HM,negf%SM,Ec,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer)

       call concat(TmpMt,zt,GreenR,1,1) 

       call destroy(GreenR)   

       do i1=1,ncont
          call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
       enddo

    enddo

    call zspectral(TmpMt,TmpMt,0,negf%rho)

    call destroy(TmpMt)

  end subroutine contour_int_n


!--------------------------------------------------------------------------------

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

    open(101, file='dos.dat', form='formatted')

    it = negf%iteration
    ncont = negf%str%num_conts
    nbl = negf%str%num_PLs

    muref = negf%mu_p

    kbT = Kb * negf%Temp * negf%eneconv

    Omega = negf%n_kt * kbT
    Lambda = 2.d0* negf%n_poles * KbT * pi
    NumPoles = negf%n_poles
    
    outer = 1 !no contacts no outer

    call create(TmpMt,negf%H%nrow,negf%H%ncol,negf%H%nrow)
    call init(TmpMt)

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

       if (negf%verbose.gt.80) then
          write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 1: point #',i,'/',iend,'  CPU=&
               &', id
       endif

       Ev = z1 + pnts(i) * z_diff

       ff = (1.d0,0.d0) - fermi_fc(Ev,muref,KbT)

       zt = z_diff * ff * wght(i) / (2.d0 *pi)

       call compute_contacts(Ev,negf,i,ncyc,Tlc,Tcl,SelfEneR,GS)

       call calls_eq_mem_dns(negf%HM,negf%SM,Ev,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer)

       call concat(TmpMt,zt,GreenR,1,1) 

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

       if (negf%verbose.gt.80) then
          write(6,'(a17,i3,a1,i3,a6,i3)') 'INTEGRAL 2: point #',i,'/',iend,'  CPU=&
               &', id
       endif

       Ev = z1 + pnts(i) * z_diff

       ff = (1.d0,0.d0) - fermi_fc(Ev,muref,KbT)

       zt = z_diff * ff * wght(i) / (2.d0 *pi)

       call compute_contacts(Ev,negf,negf%Np_p(1)+i,ncyc,Tlc,Tcl,SelfEneR,GS)

       call calls_eq_mem_dns(negf%HM,negf%SM,Ev,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer) 

       call concat(TmpMt,zt,GreenR,1,1)  !TmpMt=TmpMt+GreenR

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

       if (negf%verbose.gt.80) then
          write(6,'(a17,i3,a1,i3,a6,i3)') 'POLES: point #',i,'/',iend,'  CPU=', id
       endif

       Ev =  muref + j * KbT *pi* (2.d0*real(i,dp) - 1.d0)   

       zt= j*KbT*(1.d0,0.d0) 

       call compute_contacts(Ev,negf,negf%Np_p(1)+negf%Np_p(2)+i,ncyc,Tlc,Tcl,SelfEneR,GS)

       call calls_eq_mem_dns(negf%HM,negf%SM,Ev,SelfEneR,Tlc,Tcl,GS,GreenR,negf%str,outer)

       call concat(TmpMt,zt,GreenR,1,1) 

       call destroy(GreenR)  

        do i1=1,ncont
          call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
       enddo

    enddo

    call zspectral(TmpMt,TmpMt,0,negf%rho)

    call destroy(TmpMt)

  end subroutine contour_int_p

!------------------------------------------------------------------------------------

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
