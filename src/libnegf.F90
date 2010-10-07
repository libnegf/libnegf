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
 use mat_def


 implicit none
 private

 public :: init_negf, destroy_negf, compute_dos, contour_int_n, contour_int_p


contains
  
  subroutine init_negf(negf)
    type(Tnegf), pointer :: negf
    Integer :: ncont, nbl
    Integer, dimension(:), allocatable :: PL_end, cont_end, surf_end, cblk

    open(101, file=negf%file_re_H, form='formatted')
    open(102, file=negf%file_im_H, form='formatted')   !open imaginary part of H

    call read_H(101,102,negf%H)

    close(101)
    close(102)

    if(.not.negf%isSid) then

       open(101, file=negf%file_re_S, form='formatted')
       open(102, file=negf%file_im_S, form='formatted')   !open imaginary part of H

       call read_H(101,102,negf%S)

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
    read(101,*) cblk(1:ncont)
    read(101,*) cont_end(1:ncont)
    read(101,*) surf_end(1:ncont)

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


  subroutine destroy_negf(negf)
    type(Tnegf), pointer :: negf   

    if (allocated(negf%H%nzval)) call destroy(negf%H) 
    if (allocated(negf%S%nzval)) call destroy(negf%S)

    if (allocated(negf%rho%nzval)) call destroy(negf%rho)
    if (allocated(negf%rho_eps%nzval)) call destroy(negf%rho_eps)

    call kill_Tstruct(negf%str)   

  end subroutine destroy_negf


  subroutine compute_dos(negf) 
    type(Tnegf), pointer :: negf

    integer :: outer, N, k, i
    complex(dp) :: Ec
    Type(z_CSR) :: SelfEneR(1), Tlc(1), Tcl(1), gsurfR(1), Gr
    Type(z_DNS) :: Inv   


    outer = 0

    N = nint((negf%Emax-negf%Emin)/negf%Estep)

    print*,'N=',N
    print*,'delta=',negf%delta


    call create(SelfEneR(1),0,0,0)
    call create(Tlc(1),0,0,0)
    call create(Tcl(1),0,0,0)
    call create(gsurfR(1),0,0,0)
    !call create(Inv,negf%H%nrow,negf%H%nrow)

    open(101,file='dos.dat')

    do i=1,N
  
       Ec=(negf%Emin+i*negf%Estep)+negf%delta*(0.d0,1.d0)

       call calls_eq_mem(negf%H,negf%S,Ec,SelfEneR,Tlc,Tcl,gsurfR,Gr,negf%str,outer)
            
       !call csr2dns(Gr,Inv)
       !negf%dos=0.d0
       !do k = 1,negf%H%nrow
       !   negf%dos = negf%dos - aimag( Inv%val(k,k) )/pi
       !end do
                       
       negf%dos = -aimag( trace(Gr) )/pi

       write(101,*) real(Ec), negf%dos
       print*, real(Ec), negf%dos

       call destroy(Gr)

    enddo

    call destroy(SelfEneR(1))
    call destroy(Tcl(1))
    call destroy(Tlc(1))
    call destroy(gsurfR(1))

    !call destroy(Inv)

    call writememinfo(6)

    close(101)   

  end subroutine compute_dos


  subroutine contour_int_n(negf)

    type(Tnegf), pointer :: negf 
    Type(z_CSR) :: SelfEneR(1), Tlc(1), Tcl(1), gsurfR(1)
    type(z_CSR) :: GreenR, TmpMt 

    integer :: npid, istart, iend, NumPoles
    integer :: i, outer

    real(dp), DIMENSION(:), ALLOCATABLE :: wght,pnts   ! Gauss-quadrature points
    real(dp) :: Omega, Lambda
    real(dp) :: muref, kbT

    complex(dp) :: z1,z2,z_diff, zt
    complex(dp) :: Ec, ff

    open(101, file='dos.dat', form='formatted')

    muref = negf%mu_n

    kbT = Kb * negf%Temp * negf%eneconv

    Omega = negf%n_kt * kbT
    Lambda = 2.d0* negf%n_poles * KbT * pi
    NumPoles = negf%n_poles

    call create(SelfEneR(1),0,0,0)
    call create(Tlc(1),0,0,0)
    call create(Tcl(1),0,0,0)
    call create(gsurfR(1),0,0,0)
    outer = 0 !no contacts no outer

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

       call calls_eq_mem(negf%H,negf%S,Ec,SelfEneR,Tlc,Tcl,gsurfR,GreenR,negf%str,outer)

       call concat(TmpMt,zt,GreenR,1,1) 

       call destroy(GreenR)      

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

       call calls_eq_mem(negf%H,negf%S,Ec,SelfEneR,Tlc,Tcl,gsurfR,GreenR,negf%str,outer) 

       call concat(TmpMt,zt,GreenR,1,1)  !TmpMt=TmpMt+GreenR

       call destroy(GreenR)

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

       call calls_eq_mem(negf%H,negf%S,Ec,SelfEneR,Tlc,Tcl,gsurfR,GreenR,negf%str,outer)

       call concat(TmpMt,zt,GreenR,1,1) 

       call destroy(GreenR)   

    enddo

    call zspectral(TmpMt,TmpMt,0,negf%rho)

    call destroy(TmpMt)
    call destroy(SelfEneR(1))
    call destroy(Tcl(1))
    call destroy(Tlc(1))
    call destroy(gsurfR(1)) 

  end subroutine contour_int_n




  subroutine contour_int_p(negf)

    type(Tnegf), pointer :: negf 
    Type(z_CSR) :: SelfEneR(1), Tlc(1), Tcl(1), gsurfR(1)
    type(z_CSR) :: GreenR, TmpMt 

    integer :: npid, istart, iend, NumPoles
    integer :: i, outer

    real(dp), DIMENSION(:), ALLOCATABLE :: wght,pnts   ! Gauss-quadrature points
    real(dp) :: Omega, Lambda
    real(dp) :: muref, kbT

    complex(dp) :: z1,z2,z_diff, zt
    complex(dp) :: Ev, ff

    open(101, file='dos.dat', form='formatted')

    muref = negf%mu_p

    kbT = Kb * negf%Temp * negf%eneconv

    Omega = negf%n_kt * kbT
    Lambda = 2.d0* negf%n_poles * KbT * pi
    NumPoles = negf%n_poles

    call create(SelfEneR(1),0,0,0)
    call create(Tlc(1),0,0,0)
    call create(Tcl(1),0,0,0)
    call create(gsurfR(1),0,0,0)
    
    outer = 0 !no contacts no outer

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

       call calls_eq_mem(negf%H,negf%S,Ev,SelfEneR,Tlc,Tcl,gsurfR,GreenR,negf%str,outer)

       call concat(TmpMt,zt,GreenR,1,1) 

       call destroy(GreenR)      

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

       call calls_eq_mem(negf%H,negf%S,Ev,SelfEneR,Tlc,Tcl,gsurfR,GreenR,negf%str,outer) 

       call concat(TmpMt,zt,GreenR,1,1)  !TmpMt=TmpMt+GreenR

       call destroy(GreenR)

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

       call calls_eq_mem(negf%H,negf%S,Ev,SelfEneR,Tlc,Tcl,gsurfR,GreenR,negf%str,outer)

       call concat(TmpMt,zt,GreenR,1,1) 

       call destroy(GreenR)   

    enddo

    call zspectral(TmpMt,TmpMt,0,negf%rho)

    call destroy(TmpMt)
    call destroy(SelfEneR(1))
    call destroy(Tcl(1))
    call destroy(Tlc(1))
    call destroy(gsurfR(1))

  end subroutine contour_int_p


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



end module libnegf
