module complexbands
  use ln_precision
  use ln_constants
  use inversions
  use lowdin
  use outmatrix
  implicit none
  private

  public :: TComplexBandPar, compute_bands, complex_k
  public :: states_at_k, band_velocity, band_velocities
  public :: sort_and_normalize, TStatesSummary

  private :: diagonalize

  real(dp), parameter, public :: EPS18 = 1d-18
  real(dp), parameter, public :: EPS12 = 1d-12  
  real(dp), parameter :: MAXK = 18.0_dp
  real(dp), parameter :: NULK = 30.0_dp
 
  type TComplexBandPar
    integer :: at_start,at_end
    integer :: mat_start,mat_end
    real(dp) :: L
    real(dp) :: emin, emax, estep
    integer :: GW
    integer :: lowdin_order
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
    
    complex(dp) :: HH(:,:) 	    ! hamiltonian
    complex(dp) :: SS(:,:)	    ! overlap
    type(TComplexBandPar) :: par    ! parameters 
    
    complex(dp), ALLOCATABLE, DIMENSION(:,:) :: HM,SM  ! block HAM and OVR
    complex(dp), ALLOCATABLE, DIMENSION(:,:) :: TM,ST  ! block-block HAM and OVR
    
    complex(dp), ALLOCATABLE, DIMENSION(:,:) :: SelfEneGW! GW Self-energy
    complex(dp), ALLOCATABLE, DIMENSION(:,:,:) :: Dummy! Dummy matrix for GWself  
    complex(dp), ALLOCATABLE, DIMENSION(:,:) :: Z0,Z1,ZL,Z2
    integer :: err,Nstep, ndim
    logical :: exis
    integer :: i,k,PLdim,Sdim 
    real(dp) :: E, folding
    
    complex(dp), ALLOCATABLE, DIMENSION(:) :: kz
    
    call set_contactHS(HH,SS,par,HM,TM,SM,ST,PLdim)

    Sdim = 2*PLdim
    folding = 1.d0*(1+par%lowdin_order)
    IF(par%lowdin_order.le.-1) THEN
       folding=1.d0;
    ELSE IF(par%lowdin_order.eq.11) THEN
       folding=2.d0;
    END IF

    write(*,*) 'starting complex bandstructure'
    ! -----------------------------------------------------------------  
    ! START WITH COMPLEX BAND COMPUTATION
    ! -----------------------------------------------------------------
    ALLOCATE(Z0(PLdim,PLdim),stat=err)
    ALLOCATE(Z1(PLdim,PLdim),stat=err)
    ALLOCATE(Z2(PLdim,PLdim),stat=err)
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
    write(*,*) 'Complex Band: ',par%emin,par%emax
    write(*,*) 'Number of steps: ',Nstep+1

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
       write(106,'(F15.8)',ADVANCE='NO') E
       write(107,'(F15.8)',ADVANCE='NO') E
       write(108,'(F15.8)',ADVANCE='NO') E

       !write(*,*) "E=",E
       Z0 = E*SM - HM
       Z1 = E*ST - TM
       
       !Z0 = E*SM - HM
       !Z1 = - TM      
       !Z2 = - HH(par%mat_start+PLdim:par%mat_end,&
       !        par%mat_start:par%mat_start+PLdim-1)

       !if (ANY(conjg(transpose(Z2)).ne.Z1)) then
       !   print*, 'Z1 Z2+ dont agree'
       !end if
       !if (ANY(conjg(transpose(Z0)).ne.Z0)) then
       !   print*, 'Z0 Z0+ dont agree'
       !end if

       call complex_k(E,PLdim,Z0,Z1,kz)
       
       do i=1,Sdim
          write(106,'(F15.8)',ADVANCE='NO') -abs(real(kz(i)))/folding
          write(106,'(F15.8)',ADVANCE='NO') abs(aimag(kz(i)))/folding            
       enddo
       !  Real and imaginary bands are on different files
       do i=1,Sdim
          write(107,'(F15.8)',ADVANCE='NO') real(kz(i))/folding
          write(108,'(F15.8)',ADVANCE='NO') aimag(kz(i))/folding
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
    integer :: i,k, Sdim, info, LWORK, err
    complex(dp) :: zz
    real(dp) :: norm
    character(1) :: JOBVL,JOBVR

    JOBVL='N'
    if(present(Cl)) JOBVL='V'
    JOBVR='N'
    if(present(Cr)) JOBVR='V'

    Sdim = 2*PLdim
    ! -----------------------------------------------------------------
    !write(*,*) "Z0 and Z1 defined"
    
    !------------------------------------------------------------------ 
    ! Diagonalize final matrix
    ! [ Z0   Z1 ] [C0]  =  exp(-ik) [ -Z1^H  0 ] [C0]
    ! [  I    0 ] [C1]              [  0     I ] [C1]

    TA=0.d0
    TB=0.d0
    
    TA(1:PLdim,1:PLdim) = Z11
    TA(1:PLdim,PLdim+1:Sdim)  = Z12
    do i=1,PLdim
       TA(PLdim+i,i)  = 1.d0
    enddo
    
    if(present(Z21)) then
       TB(1:PLdim,1:PLdim) = -Z21
    else
       TB(1:PLdim,1:PLdim) = -conjg(transpose(Z12))
    endif
    do i=1,PLdim
       TB(PLdim+i,PLdim+i) = 1.d0
    enddo
    
    !write(*,*) "TA,TB defined"
    !------------------------------------------------------------------ 
    Ad=0.d0
    Bd=0.d0
    
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

    LWORK = WORK(1)
    deallocate(WORK)
    allocate(WORK(LWORK),stat=err)
    if(err.ne.0) STOP 'complex_k: allocation error (LWORK)'

    call zggev(JOBVL,JOBVR,Sdim, TA, Sdim, TB, Sdim, Ad, Bd, &
            Vl,Sdim,Vr,Sdim, WORK, LWORK, RWORK, info)     

    if(info.gt.0) STOP 'complex band: QZ failed in zggev' 
    if(info.gt.Sdim) STOP 'complex band: zggev failed'
    
    !write(*,*) "T diagonalized"
    !------------------------------------------------------------------ 
    ! find complex wavevectors by computing log(l)  
    
    do i=1,Sdim
       
       ! kz = j*log(z) = -Arg(z) + j*log|z|  : -Pi < Arg(z) <=t Pi 
       ! ==> -Pi <= real(kz) < Pi 
       ! Numerical noise is set to im(kz) == +/-30.0 
       if (abs(Ad(i)).gt.EPS18 .and. abs(Bd(i)).gt.EPS18) then
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

          if(abs(Bd(i)).le.EPS18) then
             kz(i)=+j*NULK     ! this belongs to the null-space
          endif
          if(abs(Ad(i)).le.EPS18) then
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
    complex(dp) :: tmp, tmp2
    real(dp) :: norm
    integer :: n

    ! Using Hellmann-Feynman we get velocity from dH/dk mat-el:
    ! ( In fact it is generalized for S as <Cn|(ES-H)exp(ik)|Cn>
    tmp = dot_product(Cn, matmul(Z12*exp(j*kk),Cn) )
    norm = dot_product(Cn,Cn)
    
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
    real(dp), dimension(PLdim) :: EW
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

    integer :: i,m,n, cnt
    integer :: n_prop_p, n_prop_n  ! _p == +  out from D
    integer :: n_dec_p, n_dec_n    ! _n == -  in the D
    integer :: n_null_p, n_null_n
    complex(dp) :: D(Sdim,Sdim),invD(Sdim,Sdim)
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
       write(*,*) 'Decaying -:'
       write(*,*) kt(n_prop_n+1:n_prop_n+n_dec_n), vt(n_prop_n+1:n_prop_n+n_dec_n)
       write(*,*) 'Decaying +:'
       write(*,*) kt(n+n_prop_p+1:n+n_prop_p+n_dec_p), vt(n+n_prop_p+1:n+n_prop_p+n_dec_p)
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
    n=Sdim/2
    do i=1,Sdim
       norm = sqrt(real(dot_product(D(:,i),D(:,i))));
       D(:,i)=D(:,i) / norm
    end do

    !do i=1,Sdim
    !  print*,i,aimag(kt(i)), vt(i)
    !  print*,'C0:',real(dot_product(D(1:n,i),D(1:n,i)))
    !  print*,'C1:',real(dot_product(D(n+1:2*n,i),D(n+1:2*n,i)))
    !enddo
    ! For incoming states |C1|^2 = 1 
    ! For outcoming states  |C0|^2 = 1
    ! transfer_mat relays on the fact that Im(kt(i)) = -Im(kt(n+i)) !!
    ! In this way column i and n+i are multiplied exactly by the same factor
    do i=1,n_prop_n+n_dec_n
       norm = sqrt( 1 + exp(2*aimag(kt(i))) )  
       D(:,i)  = D(:,i) * norm
    enddo
    
    do i=n+1,n+n_prop_p+n_dec_p
       norm = sqrt( 1 + exp(-2*aimag(kt(i))) )  
       D(:,i)  = D(:,i) * norm 
    enddo
   
    !do i=1,Sdim
    !  print*,i,aimag(kt(i)), vt(i)
    !  print*,'C0:',real(dot_product(D(1:n,i),D(1:n,i)))
    !  print*,'C1:',real(dot_product(D(n+1:2*n,i),D(n+1:2*n,i)))
    !enddo

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
    !      if( abs(dot_product(D(1:n,i),D(1:n,m))).gt.1d-6) then
    !         write(*,*) 'WARNING: States not orth',i,m
    !      endif
    !   enddo
    !enddo
    !do i=n+1,n+n_prop_p
    !   do m=1,n_prop_n
    !      if( abs(dot_product(D(1:n,i),D(1:n,m))).gt.1d-6) then
    !         write(*,*) 'WARNING: States not orth',i,m
    !      endif
    !   enddo
    !   do m=n+1,i-1
    !      if( abs(dot_product(D(1:n,i),D(1:n,m))).gt.1d-6) then
    !         write(*,*) 'WARNING: States not orth',i,m
    !      endif
    !   enddo
    !enddo

 end subroutine sort_and_normalize


  ! EXTRACT CONTACTS AND MAKE LOWDIN TRANSF 
  ! depending on par%lowdin_order
  ! par%lowdin_order = -1  => Exact
  ! par%lowdin_order = -2  => Neglect overlap
  ! par%lowdin_order = 11  => First order 
  ! par%lowdin_order =  n  => Exact and truncate to order n
  subroutine set_contactHS(HH,SS,par,HM,TM,SM,ST,PLdim)
    complex(dp) :: HH(:,:) 	    ! hamiltonian
    complex(dp) :: SS(:,:)	    ! overlap
    type(TComplexBandPar) :: par    ! parameters 
    complex(dp), ALLOCATABLE, DIMENSION(:,:) :: HM,SM  ! block HAM and OVR
    complex(dp), ALLOCATABLE, DIMENSION(:,:) :: TM,ST  ! block-block HAM and OVR
    integer :: PLdim

    ! locals:
    complex(dp), ALLOCATABLE, DIMENSION(:,:) :: HM0,SM0  ! block HAM and OVR
    complex(dp), ALLOCATABLE, DIMENSION(:,:) :: TM0,ST0  ! block-block HAM and OVR
    complex(dp), ALLOCATABLE, DIMENSION(:,:) :: Z0,Z1,ZL
    integer :: i,k,PLdim0,Sdim0,Sdim,ms,nPL, err 


    write(*,'((a,I5,I5))') 'MAT start/end  = ',par%mat_start, par%mat_end
     
    Sdim0 = (par%mat_end - par%mat_start + 1)
    PLdim0= Sdim0/2 
    ms = par%mat_start

    write(*,'(a,I5,a,F8.4)') 'PL dim    = ',PLdim0,'  L= ',par%L
    ! -----------------------------------------------------------------      
    ! EXTRACT PL BLOCKS
    ! -----------------------------------------------------------------      
    ALLOCATE(HM0(PLdim0,PLdim0),stat=err)
    ALLOCATE(SM0(PLdim0,PLdim0),stat=err)
    ALLOCATE(TM0(PLdim0,PLdim0),stat=err)
    ALLOCATE(ST0(PLdim0,PLdim0),stat=err)
    if(err.ne.0) stop 'function complex bands: no space for allocation'
  
    write(*,*) 'extract PL from',size(HH,1),'x',size(HH,2)
    !extract PL Hamiltonian and PL-PL interaction 
    ! TM is the upper diagonal block linking i to i+1    
    HM0(1:PLdim0,1:PLdim0) = HH(ms:ms+PLdim0-1,ms:ms+PLdim0-1)
    SM0(1:PLdim0,1:PLdim0) = SS(ms:ms+PLdim0-1,ms:ms+PLdim0-1)
    TM0(1:PLdim0,1:PLdim0) = HH(ms:ms+PLdim0-1,ms+PLdim0:par%mat_end)
    ST0(1:PLdim0,1:PLdim0) = SS(ms:ms+PLdim0-1,ms+PLdim0:par%mat_end)

    ! -----------------------------------------------------------------  
    ! APPROXIMATE LOWDIN TRUNCATED TO ORDER lowdin_order
    ! -----------------------------------------------------------------
    PLdim = (1+par%lowdin_order)*PLdim0
    Sdim = 2*PLdim

    IF(par%lowdin_order.le.-1) THEN
       PLdim = PLdim0; Sdim=2*PLdim
    ELSE IF(par%lowdin_order.eq.11) THEN
       PLdim = 2*PLdim0; Sdim=2*PLdim
    END IF

    ALLOCATE(HM(PLdim,PLdim),stat=err)
    ALLOCATE(SM(PLdim,PLdim),stat=err)
    ALLOCATE(TM(PLdim,PLdim),stat=err)
    ALLOCATE(ST(PLdim,PLdim),stat=err)
    if(err.ne.0) stop 'function complex bands: no space for allocation'

    IF(par%lowdin_order.eq.-1) THEN
       write(*,*) 'Exact '
       HM=HM0; SM=SM0; 
       TM=TM0; ST=ST0;

    ELSE IF(par%lowdin_order.eq.-2) THEN
       write(*,*) 'Neglecting Overlap '
       HM=HM0;  
       TM=TM0; 
       ! set overlap to Id
       SM = (0.d0,0.d0)
       ST = (0.d0,0.d0)
       do i = 1, PLdim
          SM(i,i) = 1.d0
       enddo
       
    ELSE IF(par%lowdin_order.eq.11) THEN
       write(*,*) 'first order lowdin'

       do i = 1, PLdim0
          SM0(i,i) = SM0(i,i) - 1.d0
       enddo
       HM(1:PLdim0,1:PLdim0) =  (HM0 - matmul(HM0,SM0)- &
                                 matmul(TM0,transpose(ST0))- &
                                 matmul(transpose(TM0),ST0)) * 0.5d0       
       HM(1:PLdim0,PLdim0+1:2*PLdim0) =  (TM0 - matmul(HM0,ST0) - &
                                         matmul(TM0,SM0) ) * 0.5d0
       HM(PLdim0+1:2*PLdim0,1:PLdim0) = (transpose(TM0)-matmul(HM0,transpose(ST0))-&
                                         matmul(transpose(TM0),SM0)) * 0.5d0
       HM(PLdim0+1:2*PLdim0,PLdim0+1:2*PLdim0) = HM(1:PLdim0,1:PLdim0)
       HM = HM + transpose(HM)

       TM(1:PLdim0,1:PLdim0) = (-matmul(TM0,ST0)-matmul(ST0,TM0)) * 0.5d0

       TM(1:PLdim0,PLdim0+1:2*PLdim0) = (0.d0,0.d0)

       TM(PLdim0+1:2*PLdim0,1:PLdim0) = TM0 + (matmul(HM0,ST0) +&
                                        matmul(ST0,HM0) + matmul(TM0,SM0) +&
                                        matmul(SM0,TM0) ) * (-0.5d0)

       TM(PLdim0+1:2*PLdim0,PLdim0+1:2*PLdim0) = TM(1:PLdim0,1:PLdim0)
       
       ! set overlap to Id
       SM = (0.d0,0.d0)
       ST = (0.d0,0.d0)
       do i = 1, PLdim
          SM(i,i) = 1.d0
       enddo

    ELSE
       write(*,*) 'exact lowdin and truncation to order',par%lowdin_order
       
       nPL = (1+par%lowdin_order) + 3 + 3
       Sdim = nPL * PLdim0

       ALLOCATE(Z0(Sdim,Sdim),stat=err)
       ALLOCATE(Z1(Sdim,Sdim),stat=err)      
       ALLOCATE(ZL(Sdim,Sdim),stat=err) 
       if(err.ne.0) stop 'function complex bands: no space for allocation' 

       Z0 = (0.d0,0.d0)
       Z1 = (0.d0,0.d0)
       ZL = (0.d0,0.d0)

       do i = 1, nPL - 1  
          Z0( (i-1)*PLdim0+1:i*PLdim0, (i-1)*PLdim0+1:i*PLdim0) = HM0 
          Z0( (i-1)*PLdim0+1:i*PLdim0, i*PLdim0+1:(i+1)*PLdim0) = TM0
          Z0( i*PLdim0+1:(i+1)*PLdim0, (i-1)*PLdim0+1:i*PLdim0) = transpose(TM0)
          
          Z1( (i-1)*PLdim0+1:i*PLdim0, (i-1)*PLdim0+1:i*PLdim0 ) = SM0 
          Z1( (i-1)*PLdim0+1:i*PLdim0, i*PLdim0+1:(i+1)*PLdim0 ) = ST0
          Z1( i*PLdim0+1:(i+1)*PLdim0, (i-1)*PLdim0+1:i*PLdim0 ) = transpose(ST0)    
       enddo

       i =  nPL
       Z0( (i-1)*PLdim0+1:i*PLdim0, (i-1)*PLdim0+1:i*PLdim0) = HM0        
       Z1( (i-1)*PLdim0+1:i*PLdim0, (i-1)*PLdim0+1:i*PLdim0) = SM0 

       write(*,*) 'Z0 and Z1 done, goto lowdin'

       call lowdin_trans('-',Z0,Z1,ZL,Sdim)

       write(*,*) 'Lowdin done'
       ! Extract HM from central part of ZL

       ms = 2*PLdim0
       nPL = (1+par%lowdin_order)
       HM(1:PLdim,1:PLdim) = ZL(ms+1:ms+PLdim,ms+1:ms+PLdim)
       TM(1:PLdim,1:PLdim) = ZL(ms+1:ms+PLdim,ms+PLdim+1:ms+2*PLdim)

       !do k = 1, nPL-1 
       !   do i = k+1, nPL
       !      TM((k-1)*PLdim0+1:k*PLdim0, (i-1)*PLdim0+1:i*PLdim0) = (0.d0,0.d0)
       !   end do
       !end do
       write(*,*) 'free mem'
              
       deallocate(Z0,Z1,ZL)

       write(*,*) 'define identiy'
       ! set overlap to Id
       SM = (0.d0,0.d0)
       ST = (0.d0,0.d0)
       do i = 1, PLdim
          SM(i,i) = 1.d0
       enddo

       Sdim = 2*PLdim

    END IF

    deallocate(HM0,SM0,TM0,ST0)

  
  end subroutine set_contactHS
    


end module complexbands
