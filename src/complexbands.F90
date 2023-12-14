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


module complexbands
  use ln_precision
  use ln_constants
  use inversions
  implicit none
  private

  public :: TComplexBandPar, compute_bands, complex_k
  public :: states_at_k, band_velocity, band_velocities
  public :: sort_and_normalize, TStatesSummary

  private :: diagonalize

  real(dp), parameter :: MAXK = 12.0_dp
  real(dp), parameter :: NULK = 15.0_dp
 
  type TComplexBandPar
    integer :: at_start,at_end
    integer :: mat_start,mat_end
    real(dp) :: L
    real(dp) :: emin, emax, estep
    integer :: GW
 end type TComplexBandPar

 type TStatesSummary
    integer :: prop_in, prop_out
    integer :: evan_in, evan_out
    integer :: null_in, null_out
 end type TStatesSummary

contains


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !% Subroutine ComplexBands ver. 1.0
  !% The DFTB hamiltonian is used to calculate the complex bandstructure
  !% Use Boykin algorithm of generalized eigenvalue problem
  !% March 2005
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine compute_bands(HH,SS,par)
    
    complex(dp) :: HH(:,:)    ! hamiltonian
    complex(dp) :: SS(:,:)    ! overlap
    type(TComplexBandPar) :: par    ! parameters 
    
    complex(dp), ALLOCATABLE, DIMENSION(:,:) :: HM,SM  ! block HAM and OVR
    complex(dp), ALLOCATABLE, DIMENSION(:,:) :: TM,ST  ! block-block HAM and OVR
    
    complex(dp), ALLOCATABLE, DIMENSION(:,:) :: SelfEneGW! GW Self-energy
    complex(dp), ALLOCATABLE, DIMENSION(:,:,:) :: Dummy! Dummy matrix for GWself  
    complex(dp), ALLOCATABLE, DIMENSION(:,:) :: Z0,Z1
    integer :: err,Nstep
    integer :: i,k,PLdim,Sdim
    real(dp) :: E
    
    complex(dp), ALLOCATABLE, DIMENSION(:) :: kz
    
    write(*,'(2(a,I5))') 'AT start  = ',par%at_start, &
         '  MAT start  = ',par%mat_start    
    write(*,'(2(a,I5))') 'AT end    = ',par%at_end, &
         '  MAT end    = ',par%mat_end
    
    Sdim = (par%mat_end - par%mat_start + 1)
    PLdim= Sdim/2 
    
    write(*,'(a,I5,a,F8.4)') 'PL dim    = ',PLdim,'  L= ',par%L*ATU
    
    ALLOCATE(HM(PLdim,PLdim),stat=err)
    ALLOCATE(SM(PLdim,PLdim),stat=err)
    ALLOCATE(TM(PLdim,PLdim),stat=err)
    ALLOCATE(ST(PLdim,PLdim),stat=err)
    ALLOCATE(Z0(PLdim,PLdim),stat=err)
    ALLOCATE(Z1(PLdim,PLdim),stat=err)
    if(err.ne.0) stop 'function complex bands: no space for allocation'
  
    write(*,*) 'extract PL from',size(HH,1),'x',size(HH,2)
    !extract PL Hamiltonian and PL-PL interaction 
    ! TM is the upper diagonal block linking i to i+1    
    HM(1:PLdim,1:PLdim)=HH(par%mat_start:par%mat_start+PLdim-1,&
         par%mat_start:par%mat_start+PLdim-1)
    SM(1:PLdim,1:PLdim)=SS(par%mat_start:par%mat_start+PLdim-1,&
         par%mat_start:par%mat_start+PLdim-1)
    TM(1:PLdim,1:PLdim)=HH(par%mat_start:par%mat_start+PLdim-1,&
         par%mat_start+PLdim:par%mat_end)
    ST(1:PLdim,1:PLdim)=SS(par%mat_start:par%mat_start+PLdim-1,&
         par%mat_start+PLdim:par%mat_end)
    ! -----------------------------------------------------------------  

    ALLOCATE(kz(Sdim),stat=err) 
    if(err.ne.0) stop 'function complex bands: no space for allocation' 
    ! -----------------------------------------------------------------      
       
    if (par%GW.gt.0) then
       ALLOCATE(SelfEneGW(Sdim,Sdim),stat=err)
       if(err.ne.0) stop 'function complex bands: no space for allocation'

       select case(par%GW)
       case(1)
          write(*,*) 'Read QPC.DAT for GW renormalization'
          allocate(Dummy(Sdim,1,1),stat=err)
          if(err.ne.0) stop 'function complex bands: no space for allocation'          
          !call selfGW_re(SelfEneGW,Dummy,ndim,mat_start,Sdim,1,1,0)
          deallocate(Dummy)
       case(2)
          write(*,*) 'Read Energy-dependent QP/QPC.DAT for GW renormalization'  
       end select
       
       HM = HM + SelfEneGW
       TM(1:PLdim,1:PLdim) = TM(1:PLdim,1:PLdim) + SelfEneGW(1:PLdim,PLdim+1:Sdim)  
       
    endif
    
    ! START ENERGY CYCLE:
    Nstep=nint((par%emax-par%emin)/par%estep) 
    write(*,*) 'Complex Band: ',par%emin*HAR,par%emax*HAR
    write(*,*) 'Number of steps: ',Nstep

    open(106,file="cmplxband.dat")
    open(107,file="cmplxband_re.dat")
    open(108,file="cmplxband_im.dat")
    
    write(*,*) 'Energy loop'

    do k=0,Nstep
       
       E=par%emin + k*par%estep
       
       if (par%GW.eq.2) then
          allocate(Dummy(Sdim,1,1),stat=err)
          !call selfGW_re(SelfEneGW,Dummy,ndim,mat_start,Sdim,1,1,k+1)    
          deallocate(Dummy)
          HM(1:PLdim,1:PLdim)=HH(par%mat_start:par%mat_start+PLdim-1, &
               par%mat_start:par%mat_start+PLdim-1)+ SelfEneGW(1:PLdim,1:PLdim)  
          TM(1:PLdim,1:PLdim)=HH(par%mat_start:par%mat_start+PLdim-1,&
               par%mat_start+PLdim:par%mat_end)+ SelfEneGW(1:PLdim,PLdim+1:Sdim)  
       endif
       
       ! --------------------------------------------------------------
       ! WRITE ON 3 FILES 
       ! The files contains both the real and imaginary bands
       ! one file contains the pure real and pure imaginary parts
       !
       write(106,'(F15.8)',ADVANCE='NO') E*HAR
       write(107,'(F15.8)',ADVANCE='NO') E*HAR
       write(108,'(F15.8)',ADVANCE='NO') E*HAR

       !rite(*,*) "E=",E
       Z0 = E*SM - HM
       Z1 = E*ST - TM

       call complex_k(E,PLdim,Z0,Z1,kz)
       
       do i=1,Sdim
          write(106,'(F15.8)',ADVANCE='NO') -abs(real(kz(i)))
          write(106,'(F15.8)',ADVANCE='NO') abs(aimag(kz(i)))            
       enddo
       !  Real and imaginary bands are on different files
       do i=1,Sdim
          write(107,'(F15.8)',ADVANCE='NO') real(kz(i))
          write(108,'(F15.8)',ADVANCE='NO') aimag(kz(i))
       enddo
       
       write(106,'(F15.8)',ADVANCE='YES') 
       write(107,'(F15.8)',ADVANCE='YES') 
       write(108,'(F15.8)',ADVANCE='YES') 
       
       !------------------------------------------------------------------
       
    end do
    
    close(106)
    close(107)
    close(108)
    
    deallocate(kz)
    deallocate(HM,SM)
    deallocate(TM,ST)
    
  end subroutine compute_bands
  
  ! ----------------------------------------------------------------------
  ! Calculates for given E the complex band structure,
  ! i.e. k_othagonal (complex) and the corresponding tight binding basis
  ! vectors.
  !
  ! the kz vectors are in units of Pi/L, the cell length (kz*L)
  ! kz <- [-Pi..Pi)
  subroutine complex_k(E,PLdim,Z11,Z12,kz,Z21,Cl,Cr,vf)
        
    real(dp), intent(in) :: E
    integer, intent(in) :: PLdim
    complex(dp), DIMENSION(PLdim,PLdim) :: Z11,Z12 !ES-H
    complex(dp), DIMENSION(2*PLdim) :: kz
    complex(dp), DIMENSION(PLdim,PLdim), optional :: Z21 
    complex(dp), DIMENSION(2*PLdim,2*PLdim), optional :: Cl
    complex(dp), DIMENSION(2*PLdim,2*PLdim), optional :: Cr    
    real(dp), DIMENSION(2*PLdim), optional :: vf

    ! locals ---
    complex(dp), DIMENSION(2*PLdim) :: Ad,Bd        ! 
    complex(dp), DIMENSION(2*PLdim,2*PLdim) :: TA,TB,Vl,Vr 
    complex(dp), DIMENSION(:), allocatable :: WORK, RWORK  
    integer :: i, Sdim, info, LWORK, err
    character(1) :: JOBVL,JOBVR

    JOBVL='N'
    if(present(Cl)) JOBVL='V'
    JOBVR='N'
    if(present(Cr)) JOBVR='V'

    Sdim = 2*PLdim
    ! -----------------------------------------------------------------
    !open(200,file="Z1")
    !call printmat_fr(200,"",Z1,PLdim,PLdim)
    !close(200)
    
    !write(*,*) "Z0 and Z1 defined"
    
    !------------------------------------------------------------------ 
    ! Diagonalize final matrix
    ! [ Z0   Z1 ] [C0]  =  exp(-ik) [ -Z1^H  0 ] [C0]
    ! [  I    0 ] [C1]              [  0     I ] [C1]

    TA=0.0_dp
    TB=0.0_dp
    
    TA(1:PLdim,1:PLdim) = Z11
    TA(1:PLdim,PLdim+1:Sdim)  = Z12
    do i=1,PLdim
       TA(PLdim+i,i)  = 1.0_dp
    enddo
    
    if(present(Z21)) then
       TB(1:PLdim,1:PLdim) = -Z21
    else
       TB(1:PLdim,1:PLdim) = -conjg(transpose(Z12))
    endif
    do i=1,PLdim
       TB(PLdim+i,PLdim+i) = 1.0_dp
    enddo
    
    !write(*,*) "TA,TB defined"
    !------------------------------------------------------------------ 
    Ad=0.0_dp
    Bd=0.0_dp
    
    ! compute generalized eigenproblem:
    ! TA * c = lambda * TB * c
    ! 
    ! lambda(j) = [Rd(j) + I * Id(j)] / Bd(j)
    !
    allocate(WORK(2))
    allocate(RWORK(8*Sdim), stat=err)
    if(err.ne.0) STOP 'complex_k: allocation error (RWORK)'
    
    LWORK = -1 ! Query LWORK size
    call zggev(JOBVL,JOBVR,Sdim, TA, Sdim, TB, Sdim, Ad, Bd, &
            Vl,Sdim,Vr,Sdim, WORK, LWORK, RWORK, info) 

    LWORK = nint(real(WORK(1)))
    deallocate(WORK)
    allocate(WORK(LWORK),stat=err)
    if(err.ne.0) STOP 'complex_k: allocation error (LWORK)'

    call zggev(JOBVL,JOBVR,Sdim, TA, Sdim, TB, Sdim, Ad, Bd, &
            Vl,Sdim,Vr,Sdim, WORK, LWORK, RWORK, info)     

    if(info.gt.0) STOP 'function complex band: dggev did not converge' 
    
    !write(*,*) "T diagonalized"
    !------------------------------------------------------------------ 
    ! find complex wavevectors by computing log(l)  
    
    do i=1,Sdim
       
       ! kz = j*log(z) = -Arg(z) + j*log|z|  : -Pi < Arg(z) <= Pi 
       ! ==> -Pi <= real(kz) < Pi 
       ! Numerical noise is set to im(kz) == -15.0 
       if (abs(Ad(i)).gt.EPS12 .and. abs(Bd(i)).gt.EPS12) then
          kz(i)=j*log(Ad(i)/Bd(i))
          ! Extremely fast-decaying or fast growing states are
          ! set to a finite value
          if (aimag(kz(i)).lt.-MAXK) then 
             kz(i)=-j*MAXK 
          endif
          if (aimag(kz(i)).gt.+MAXK) then 
             kz(i)=+j*MAXK
          endif

       else

          if(abs(Bd(i)).le.EPS12) then
             kz(i)=+j*NULK     ! this belongs to the null-space
          endif
          if(abs(Ad(i)).le.EPS12) then
             kz(i)=-j*NULK     ! this belongs to the null-space
          endif

       endif
                  
    enddo
    
    if(present(Cr)) then
       Cr = Vr
    end if
    if(present(vf)) then
       do i=1,Sdim
          if( abs(aimag(kz(i))) .lt. EPS12 ) then
             vf(i) = band_velocity(real(kz(i)),PLdim,Z12,Vr(:,i))
          else
             vf(i) = 0.0_dp
          endif
       enddo
    end if
    if(present(Cl)) then
       Cl = Vl       
    end if

  end subroutine complex_k

  !---------------------------------------------------
  ! Computes the real band-structure of a Layer
  subroutine states_at_k(kk,PLdim,HM,SM,TM,ST,Ev,Cr)
    complex(dp), intent(in) :: kk
    integer, intent(in) :: PLdim
    complex(dp), dimension(PLdim,PLdim) :: HM,SM,TM,ST
    real(dp), dimension(PLdim) :: Ev
    complex(dp), dimension(PLdim,PLdim), optional :: Cr

    ! locals:
    complex(dp), dimension(PLdim,PLdim) :: SS, HH


    HH = HM + TM*exp(j*kk) + transpose(TM)*exp(-j*conjg(kk))
    SS = SM + ST*exp(j*kk) + transpose(ST)*exp(-j*conjg(kk))
    
    if(present(Cr)) then
       call diagonalize(PLdim,HH,SS,Ev,Cr)
    else
       call diagonalize(PLdim,HH,SS,Ev)
    endif

  end subroutine states_at_k

  !--------------------------------------------------------------
  ! Computation of the band velocity using Hellmann-Feynamn th.
  !
  !
  !                [  ->+[          i k d ] ->  ]
  !  v =  - 2  * Im|  C  | |H      e      | C   |
  !   n            [   n [   n,n+1        ]  n  ]
  !  
  ! hbar == 1, L simplifies with the normalization of Cn
  ! Assumes that the Cn are orthonormal matrices
  !-------------------------------------------------------------- 
  subroutine band_velocities(kk,PLdim,TM,ST,Ev0,Cn,vf1)
    real(dp), intent(in) :: kk
    integer, intent(in) :: PLdim
    complex(dp), dimension(PLdim,PLdim), intent(in) :: TM,ST
    real(dp), dimension(PLdim), intent(in) :: Ev0
    complex(dp), dimension(PLdim,PLdim), intent(in) :: Cn    
    real(dp), dimension(PLdim), intent(out) :: vf1    
    ! locals:
    complex(dp), dimension(PLdim,PLdim) :: tmp, tmp2
    integer :: n

    ! Using Hellmann-Feynman we get velocity from dH/dk mat-el:
    tmp=matmul(transpose(conjg(Cn)),matmul(TM*exp(j*kk),Cn))
    tmp2=matmul(transpose(conjg(Cn)),matmul(ST*exp(j*kk),Cn))

    do n=1,PLdim
       vf1(n)=-2.0_dp*(aimag(tmp(n,n))-Ev0(n)*aimag(tmp2(n,n)))
    enddo

  end subroutine band_velocities
  !
  ! Like band_velocities, but one at the time
  ! Cn is NOT assumed normalized. 
  !-------------------------------------------------------------- 
  function band_velocity(kk,PLdim,Z12,Cn) result(vf)
    real(dp), intent(in) :: kk
    integer, intent(in) :: PLdim
    complex(dp), dimension(PLdim,PLdim), intent(in) :: Z12
    complex(dp), dimension(PLdim), intent(in) :: Cn    
    real(dp) :: vf    
    ! locals:
    complex(dp) :: tmp
    real(dp) :: norm

    ! Using Hellmann-Feynman we get velocity from dH/dk mat-el:
    ! ( In fact it is generalized for S as <Cn|(ES-H)exp(ik)|Cn>
    tmp = dot_product(Cn, matmul(Z12*exp(j*kk),Cn) )
    norm = real(dot_product(Cn,Cn))
    
    vf=2.0_dp*aimag(tmp)/norm

  end function band_velocity


  !------------------------------------------------------------
  ! Wrapper to LAPACK call for diagonalization
  !------------------------------------------------------------ 
  subroutine diagonalize(PLdim,HH,SS,Ev,Cr)
    integer :: PLdim
    complex(dp), dimension(PLdim,PLdim) :: HH, SS
    real(dp), dimension(PLdim) :: Ev
    complex(dp), dimension(PLdim,PLdim), optional :: Cr
    
    integer :: ier
    complex(dp), dimension(2*PLdim) :: zAux
    real(dp), dimension(3*PLdim) :: Aux
    character(1) :: JOBV 

    JOBV='N'
    if(present(Cr)) JOBV='V'
    

    call ZHEGV(1,JOBV,'L',PLdim,HH,PLdim,SS,PLdim,&
         Ev,Aux,3*PLdim,zAux,ier)

    if(ier.ne.0) write(*,*) 'ZHEGV ERROR n.',ier

    if(present(Cr)) then
       Cr=HH
    endif

  end subroutine diagonalize

  !------------------------------------------------------------
  ! Sort the PL eigenstates as in Vogl's PRB 62, 7289 (2000)
  !                               Ting's PRB 45, 3583 (1992)
  ! NB:
  ! -- -- -- -- ------ -- -- -- -- 
  ! ..|C2|C1|C0|Device|C0|C1|C2|...
  ! -- -- -- -- ------ -- -- -- -- 
  ! The contact layers are always going out from the Device.
  ! ==> C1 = exp(+i|k|L)*C0
  ! k = Re[k]+iIm[k] ==> C1=exp{iRe[k]}*exp{-Im[k]}*C0
  !
  ! Propagating waves (Im[k]=0)
  ! propagating towards the device vf[k]<0   (v-)
  ! propagating away have          vf[k]>0   (v+)
  ! 
  ! Evanescent waves towards the device:  Im[k]<0   (d-)
  ! Evanescent waves away from device:    Im[k]>0   (d+)
  !
  ! n- and n+ are the null spaces
  !
  ! The eigenstates are arranged such that:
  !
  ! [v-, d-, n-, v+, d+, n+]
  !
  ! Degenerate states are orthogonalized and normalized
  ! Note: The matrix is not orthogonal, but non-singular
  !------------------------------------------------------------     
  subroutine sort_and_normalize(Sdim,kz,vf,Cr,summary)
    integer, intent(in) :: Sdim
    complex(dp) :: kz(Sdim)
    real(dp) :: vf(Sdim)
    complex(dp) :: Cr(Sdim,Sdim)
    type(TStatesSummary), optional :: summary

    integer :: i,n
    integer :: n_prop_p, n_prop_n  ! _p == +  out from D
    integer :: n_dec_p, n_dec_n    ! _n == -  in the D
    integer :: n_null_p, n_null_n
    complex(dp) :: D(Sdim,Sdim)
    complex(dp) :: kt(Sdim)
    real(dp) :: vt(Sdim)
    real(dp) :: tmp, norm

    n = Sdim/2
    n_prop_p=0; n_prop_n=0; 
    n_dec_p=0; n_dec_n=0;
    n_null_p=0; n_null_n=0;

    ! PROPAGATING STATES
    do i = 1,Sdim
       if(vf(i).lt.0) then
          n_prop_n=n_prop_n+1
          D(:,n_prop_n) = Cr(:,i)
          kt(n_prop_n) = kz(i)
          vt(n_prop_n) = vf(i)
       endif
       if(vf(i).gt.0) then
          n_prop_p=n_prop_p+1
          D(:,n+n_prop_p) = Cr(:,i)
          kt(n+n_prop_p) = kz(i)
          vt(n+n_prop_p) = vf(i)
       endif  
    end do
    if (n_prop_p.ne.n_prop_n) then
       write(*,*) 'PROBLEM WITH P:',n_prop_p,'!=',n_prop_n
    endif

    ! EVANESCENT STATES
    do i = 1,Sdim
       if(vf(i).eq.0) then
          tmp = aimag(kz(i))
          if(tmp.ge.-MAXK .and. tmp.lt.0) then
             n_dec_n=n_dec_n+1
             D(:,n_prop_n+n_dec_n) = Cr(:,i)
             kt(n_prop_n+n_dec_n) = kz(i)
             vt(n_prop_n+n_dec_n) = vf(i)
          endif
          if(tmp.gt.0 .and. tmp.le.MAXK) then
             n_dec_p=n_dec_p+1
             D(:,n+n_prop_p+n_dec_p) = Cr(:,i)
             kt(n+n_prop_p+n_dec_p) = kz(i)
             vt(n+n_prop_p+n_dec_p) = vf(i)
          endif
       endif
    end do
    if (n_dec_p.ne.n_dec_n) then
       write(*,*) 'PROBLEM with D:',n_dec_p,'!=',n_dec_n
    endif    


    ! STATES IN THE NULL SPACES 
    do i = 1,Sdim
       if(vf(i).eq.0) then
          tmp = aimag(kz(i))
          if(tmp.eq.-NULK) then
             n_null_n=n_null_n+1
             D(:,n_prop_n+n_dec_n+n_null_n) = Cr(:,i)
             kt(n_prop_n+n_dec_n+n_null_n) = kz(i)
             vt(n_prop_n+n_dec_n+n_null_n) = vf(i)
          endif
          if(tmp.eq.NULK) then
             n_null_p=n_null_p+1
             D(:,n+n_prop_p+n_dec_p+n_null_p) = Cr(:,i)
             kt(n+n_prop_p+n_dec_p+n_null_p) = kz(i)
             vt(n+n_prop_p+n_dec_p+n_null_p) = vf(i)
          endif
       endif
    end do
    if (n_null_p.ne.n_null_n) then
       write(*,*) 'PROBLEM with N:',n_null_p,'!=',n_null_n
    endif    

    if (n_prop_p+n_dec_p+n_null_p.ne.n) then
       write(*,*) 'PROBLEM with + sum'
       write(*,*) 'Expected',n,'states'
       write(*,*) 'Got',n_prop_p+n_dec_p+n_null_p  
       write(*,*) '(',n_prop_p,n_dec_p,n_null_p,')'
    endif
    if (n_prop_n+n_dec_n+n_null_n.ne.n) then
       write(*,*) 'PROBLEM with - sum'
       write(*,*) 'Expected',n,'states'
       write(*,*) 'Got',n_prop_n+n_dec_n+n_null_n 
       write(*,*) '(',n_prop_n,n_dec_n,n_null_n,')'
    endif
    
    ! NORMALIZE ALL STATES
    do i=1,Sdim
       norm = sqrt(real(dot_product(D(:,i),D(:,i))));
       D(:,i)  = D(:,i) / norm;       
    enddo

    Cr = D
    kz = kt
    vf = vt

    if(present(summary)) then
       summary%prop_in=n_prop_n
       summary%prop_out=n_prop_p
       summary%evan_in=n_dec_n
       summary%evan_out=n_dec_p      
       summary%null_in=n_null_n
       summary%null_out=n_null_p
    endif

    ! CHECK INVERSION
    !call inverse(invD(1:n,1:n),D(1:n,1:n),n)

    !invD(1:n,1:n)=matmul(conjg(transpose(invD(1:n,1:n))),D(n+1:2*n,1:n))

    !write(*,*)
    !do i=1,n
    !   do m=1,n
    !      write(*,'(f15.8)',advance='NO') abs(D(m,i))
    !   enddo
    !   write(*,*)
    !   do m=1,n
    !      write(*,'(f15.8)',advance='NO') abs(D(n+m,i))
    !   enddo
    !   write(*,*)
    !enddo
          
    ! CHECK ORTHOGONALITY (among prop. states)
    !do i=2,n_prop_n
    !   do m=1,i-1
    !      if( abs(dot_product(D(:,i),D(:,m))).gt.1d-6) then
    !         write(*,*) 'WARNING: States not orth',i,m
    !      endif
    !   enddo
    !enddo
    !do i=n+1,n+n_prop_p
    !   do m=1,n_prop_n
    !      if( abs(dot_product(D(:,i),D(:,m))).gt.1d-6) then
    !         write(*,*) 'WARNING: States not orth',i,m
    !      endif
    !   enddo
    !   do m=n+1,i-1
    !      if( abs(dot_product(D(:,i),D(:,m))).gt.1d-6) then
    !         write(*,*) 'WARNING: States not orth',i,m
    !      endif
    !   enddo
    !enddo
  end subroutine sort_and_normalize

end module complexbands
