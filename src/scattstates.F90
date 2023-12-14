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


module scattstates
  use ln_precision
  use ln_constants
  use ln_allocation
  use mat_def
  use sparsekit_drv
  use contselfenergy, only : surface_green,selfenergy
  use complexbands
  use inversions, only : inverse
  implicit none
  private

  public :: TStruct_info
  public :: contact_states, scatt_states, transfer_mat, transfer_mat2

  !----------------------------------------------------------------           
  type TStruct_Info
     integer, dimension(:), Pointer :: PL_start      !iatm(1)
     integer, dimension(:), Pointer :: PL_end        !iatm(2)
     integer, dimension(:), Pointer :: mat_PL_start  !ind(..)
     integer, dimension(:), Pointer :: mat_PL_end    !ind(..)    
     integer, dimension(:), Pointer :: cont_start    !iatc(3,cont)
     integer, dimension(:), Pointer :: cont_end      !iatc(2,cont)
     integer, dimension(:), Pointer :: mat_C_start   !ind(cont_start)+1 
     integer, dimension(:), Pointer :: mat_C_end     !ind(cont_end+1)
     integer, dimension(:), Pointer :: cblk          !contact int block
     integer :: num_PLs
     integer :: num_conts
     integer :: active_cont
  end type TStruct_Info

contains

  ! This subroutine diagonalizes the PL Hamiltonian
  ! and find the eigenenergies. 
  subroutine contact_states(cont,PLdim,HM,SM,TM,ST,Nk,Erange)
    integer, intent(in) :: cont
    integer, intent(in) :: PLdim
    complex(dp), dimension(PLdim,PLdim) :: HM,SM,TM,ST
    integer, dimension(3), intent(in) :: Nk
    real(dp), intent(in) :: Erange(2)

    ! local variables
    integer :: n, i, ik, err, numEn
    real(dp) :: kk
    complex(dp) :: kk_c
    real(dp), dimension(PLdim) :: Ev, vf
    logical, dimension(PLdim) :: maskEn
    integer, dimension(PLdim) :: ones
    complex(dp), dimension(PLdim,PLdim) :: Vr
    character(1) :: idc

    write(idc,'(i1)') cont    
    open(112,file='contact'//idc//'_bands.dat', &
         access='DIRECT',recl=dp+4+2*PLdim*dp)
    open(222,file='contact'//idc//'_states.dat')

    ones = 1

    do ik = 1, Nk(1)

       ! 0 and Pi are excluded
       kk = ik * PI/(Nk(1)+1) 
       kk_c = kk+j*0.0_dp

       call states_at_k(kk_c,PLdim,HM,SM,TM,ST,Ev,Cr=Vr)

       call band_velocities(kk,PLdim,TM,ST,Ev,Vr,vf)
       
       ! check E-range and save on files
       where(Ev.ge.Erange(1) .and. Ev.le.Erange(2)) 
          maskEn=.true.
       elsewhere 
          maskEn=.false.
       end where

       ! write down the eigenvectors
       do i=1,PLdim
          do n=1,PLdim
             if(maskEn(n)) then
                ! Contacts are always directed away from the 
                ! central region. Incident waves have negative vf              
                ! If v(k)>0 ==> v(-k)<0 and C(-k)=C*(k)
                if(vf(n).le.0) then
                   write(222,'(F15.8,F15.8)', advance='NO') Vr(i,n)
                else
                   write(222,'(F15.8,F15.8)', advance='NO') conjg(Vr(i,n))     
                endif
             endif
          enddo
          write(222,*) 
       enddo
       write(222,*) repeat('-',40)

       !numEn=PLdim
       numEn=sum(ones,maskEn)
       Ev=pack(Ev,maskEn)
       vf=pack(vf,maskEn)

       write(112,rec=ik) kk,numEn,Ev,vf

    end do    
    close(112)
    close(222)

    call direct2seq(trim('contact'//idc//'_bands.dat'),Nk(1),PLdim)

  end subroutine contact_states


  !-----------------------------------------------------------


  subroutine scatt_states(emcont,PLdim1,HM1,SM1,TM1,ST1,&
                              PLdim2,HM2,SM2,TM2,ST2,Nk)
    integer,intent(in)  :: emcont
    integer, intent(in) :: PLdim1
    complex(dp), dimension(PLdim1,PLdim1) :: HM1,SM1,TM1,ST1
    integer, intent(in) :: PLdim2
    complex(dp), dimension(PLdim2,PLdim2) :: HM2,SM2,TM2,ST2
    integer, dimension(3), intent(in) :: Nk

    integer :: i, n, ik, iE, NumEn, Sdim, numk1, numk2
    integer :: kstep, ik2, m, mb
    real(dp) :: kk, E, norm
    real(dp), dimension(PLdim1) :: Ev, vf, Ev1, vf1, Ev1b
    real(dp), dimension(PLdim2) :: Ev2, vf2, Ev2b
    complex(dp), dimension(PLdim1) :: q11
    complex(dp), dimension(PLdim2) :: q12 
    complex(dp), DIMENSION(2*PLdim1) :: kz1
    complex(dp), DIMENSION(2*PLdim2) :: kz2
    complex(dp), DIMENSION(2*PLdim1,2*PLdim1) :: Vr1
    complex(dp), DIMENSION(2*PLdim2,2*PLdim2) :: Vr2 
    complex(dp), DIMENSION(PLdim1,PLdim1) :: Vr1b, Z11, Z12
    complex(dp), DIMENSION(PLdim2,PLdim2) :: Vr2b, Z21, Z22
    logical, DIMENSION(2*PLdim1) :: mask1
    integer, dimension(2*PLdim1) :: ones1
    logical, DIMENSION(2*PLdim2) :: mask2
    integer, dimension(2*PLdim2) :: ones2
    character(1) :: idc, idc2

    ones1 = 1
    ones2 = 1
    Vr1=(0.0_dp,0.0_dp)
    Vr2=(0.0_dp,0.0_dp)

    if(emcont.eq.1) then
       open(111,file='contact1_bands.dat', &
            access='DIRECT',recl=dp+4+2*PLdim1*dp)
       open(221,file='q11.dat') !REFLECTED STATES
       open(222,file='q12.dat') !TRANSMITTED STATES  
    elseif(emcont.eq.2) then
       open(111,file='contact2_bands.dat', &
            access='DIRECT',recl=dp+4+2*PLdim2*dp)
       open(221,file='q22.dat') !REFLECTED
       open(222,file='q21.dat') !TRANSMITTED
    else
       write(*,*) 'Invalid Emitter'
       return
    endif

    kstep = (Nk(1)+1)/(Nk(2)+1)

    do ik = kstep, Nk(1), kstep

       ! k states incoming from contact 1
       read(111,rec=ik) kk, numEn, Ev, vf
       write(*,*) kk,numEn

       do iE = 1, numEn
          
          E = Ev(iE)
          write(*,*) iE,E

          ! STATES REFLECTED IN CONTACT 1
          ! ===================================================
          ! -------------------------------
          ! ..|C2|C1|C0|Device|C0|C1|C2|...
          ! -------------------------------
          ! C1 = exp(ikL)*C0
          ! k = Re[k]+iIm[k] ==> C1=exp{iRe[k]}*exp{-Im[k]}*C0
          ! Therefore, for evanescent waves:  Im[k]>0
          ! Propagating waves:  Im[k]=0  ==>  vf[k]>0
          Z11=E*SM1 - HM1
          Z12=E*ST1 - TM1

          call complex_k(E,PLdim1,Z11,Z12,kz1,Cr=Vr1)

          mask1=.false.

          do n = 1,2*PLdim1

             if(aimag(kz1(n)).gt.0.0_dp .and. &
                  abs(real(kz1(n))).lt.EPS12) then
                mask1(n)=.true.
                ! Build the bloch state: (diagonalize again:)
                kz1(n)=0.0_dp+j*aimag(kz1(n))              
                call states_at_k(kz1(n),PLdim1,HM1,SM1,TM1,ST1,Ev1b,Vr1b)
                ! Find the closest state
                mb = minloc(abs(Ev1b-E),1)
                ! Get corresponding Bloch state
                Vr1(1:PLdim1,n) = Vr1b(1:PLdim1,mb)
             endif
             ! for each real k-state check velocity.
             ! real kz are in (-Pi..Pi)
             ! but stored kz are only positive. 
             if(abs(aimag(kz1(n))).lt.EPS12 .and. &
                  abs(real(kz1(n))).gt.0.0_dp) then
                !interpolate the k-point (find closest value)
                ik2 = interp_k(real(kz1(n)),Nk(1))
                ! read record from direct file 
                read(111,rec=ik2) kk, numEn, Ev1, vf1
                ! search the Ev1 closest to input E 
                m = minloc(abs(Ev1-E),1)
                ! reflected states have vf(k) > 0
                ! if k<0 ==> vf(-k)<0 || k>0 ==> vf(k)>0
                if(vf1(m)*real(kz1(n)).gt.0) then
                   write(*,'(10x,i6,3(F15.8))') ik2,real(kz1(n)),Ev1(m),vf1(m) 
                   mask1(n)=.true. 
                   ! Build the bloch state: (diagonalize again:)
                   kz1(n)=real(kz1(n))+j*0.0_dp
                   call states_at_k(kz1(n),PLdim1,HM1,SM1,TM1,ST1,Ev1b,Vr1b)
                   ! Find the closest state
                   mb = minloc(abs(Ev1b-Ev1(m)),1)
                   ! Get corresponding Bloch state
                   Vr1(1:PLdim1,n) = Vr1b(1:PLdim1,mb)
                   !
                endif
             endif

          enddo

          ! WRITE DOWN THE STATES ON FILE           
          numk1=sum(ones1,mask1)
          
          if(numk1.gt.PLdim1) then
             write(*,*) 'ERROR: more than PLdim1 states',numk1
             stop
          endif
    
          ! store eigenvales
          q11= pack(kz1,mask1)

          write(221,'(i4)', advance='NO') numk1 
          do i=1,numk1
             write(221,'(F15.8,F15.8)', advance='NO') q11(i)          
          enddo
          write(221,*)

          ! write down the eigenvectors (Bloch-states)
          do i=1,PLdim1
             do n=1,2*PLdim1
                if(mask1(n)) then
                   write(221,'(F15.8,F15.8)', advance='NO') Vr1(i,n)
                endif
             enddo
             write(221,*)
          enddo

          ! STATES TRANSMITTED IN CONTACT 2          
          ! ===================================================
          Z21=E*SM2 - HM2
          Z22=E*ST2 - TM2

          call complex_k(E,PLdim1,Z21,Z22,kz2,Cr=Vr2)

          mask2=.false.

          do n = 1,2*PLdim2

             if(aimag(kz2(n)).gt.0.0_dp .and. &
                  abs(real(kz2(n))).lt.EPS12) then
                mask2(n)=.true.
                ! Build the bloch state: (diagonalize again:)
                kz2(n)=0.0_dp+j*aimag(kz2(n))
                call states_at_k(kz2(n),PLdim2,HM2,SM2,TM2,ST2,Ev2b,Vr2b)
                ! Find the closest state
                mb = minloc(abs(Ev2b-E),1)
                ! Get corresponding Bloch state
                Vr2(1:PLdim2,n) = Vr2b(1:PLdim2,mb)
             endif
             ! for each real k-state check velocity.
             ! real kz are in (-Pi..Pi)
             ! but stored kz are only positive. 
             if(abs(aimag(kz2(n))).lt.EPS12 .and. &
                  abs(real(kz2(n))).gt.0.0_dp) then
                !interpolates the k-point
                ik2 = interp_k(real(kz2(n)),Nk(1))
                ! read record from direct file 
                read(111,rec=ik2) kk, numEn, Ev2, vf2
                ! search the Ev1 closest to input E 
                m = minloc(abs(Ev2-E),1)
                ! reflected states have vf(k) > 0
                ! if k<0 ==> vf(-k)<0 || k>0 ==> vf(k)>0
                if(vf2(m)*real(kz2(n)).gt.0) then
                   write(*,'(10x,i6,3(F15.8))') ik2,real(kz2(n)),Ev2(m),vf2(m)
                   mask2(n)=.true.  
                   ! Build the bloch state: (diagonalize again:)
                   kz2(n)=real(kz2(n))+j*0.0_dp
                   call states_at_k(kz2(n),PLdim2,HM2,SM2,TM2,ST2,Ev2b,Vr2b)
                   ! Find the closest state
                   mb = minloc(abs(Ev2b-Ev2(m)),1)
                   ! Get corresponding Bloch state
                   Vr2(1:PLdim2,n) = Vr2b(1:PLdim2,mb)
                endif
             endif
          enddo
 
          ! WRITE DOWN THE STATES ON FILE          
          numk2=sum(ones2,mask2)
          
          if(numk2.gt.PLdim2) then
             write(*,*) 'ERROR: more than PLdim2 states'
             stop
          endif
          
          ! store eigenvales
          q12= pack(kz2,mask2)
          
          write(222,'(i4)', advance='NO') numk2 
          do i=1,numk2
             write(222,'(F15.8,F15.8)', advance='NO') q12(i)          
          enddo
          write(222,*)
          
          ! write down the eigenvector
          do i=1,PLdim2
             do n=1,2*PLdim2
                if(mask2(n)) then
                   write(222,'(F15.8,F15.8)', advance='NO') Vr2(i,n)
                endif
             enddo
             write(222,*)
          enddo       

       enddo

    end do
   
    close(111)
    close(221)
    close(222)

  end subroutine scatt_states

  !-----------------------------------------------------------
  function interp_k(kk,Nk) result(ik)
    real(dp), intent(in) :: kk
    integer, intent(in) :: Nk
    integer :: ik

    !find k (ik-index) closest to kk
    ! k <- Pi/(Nk+1),...,ik*Pi/(Nk+1), ... Nk*Pi/(Nk+1)
    ik = nint((Nk+1)*abs(kk)/Pi)
    
  end function interp_k


  !-----------------------------------------------------------
  subroutine direct2seq(filein,Nk,PLdim)
    character(50) :: filein
    integer, intent(in) :: Nk,PLdim
    
    integer :: ik, n, numEn
    real(dp) :: kk
    real(dp), dimension(PLdim) :: Ev,vf

    open(112,file=trim(filein),access='DIRECT',recl=dp+4+2*PLdim*dp)
    open(111,file='seq_'//trim(filein))


    do ik = 1,Nk

       read(112,rec=ik) kk,numEn,Ev,vf
       
       write(111,'(F15.8,i4)', advance='NO') kk, numEn 
       do n=1,numEn
          write(111,'(F15.8)', advance='NO') Ev(n)*HAR        
       enddo
       do n=1,numEn
          write(111,'(F15.8)', advance='NO') vf(n)*HAR
       enddo     
       do n=2*numEn+1,2*PLdim
          write(111,'(F15.8)', advance='NO') 0.0_dp          
       enddo   
       write(111,*)

    enddo

    close(112)
    close(111)

  end subroutine direct2seq

  !-----------------------------------------------------------

  subroutine transfer_mat(ni,nf,HH,SS,str_info,Erange)
    type array2
       real(dp), dimension(:), allocatable :: v
       complex(dp), dimension(:), allocatable :: c
    end type array2
    integer, intent(in) :: ni,nf
    type(z_CSR), intent(in) :: HH,SS
    type(TStruct_info) :: str_info
    real(dp), dimension(3), intent(in) :: Erange

    ! -------------------------------------------------------
    integer :: PLdim(2), Nk
    type(z_DNS) :: H11(2),S11(2),H12(2),S12(2),GS_d(2),Vr(2) 
    type(z_DNS) :: Z0(2),Z1(2),invD12(2),RHS,TT1,TT2,TT(2),TR,RR
    type(z_DNS) :: Self
    type(z_CSR) :: TM(2), ST(2), GS(2),Tlc(2),Tcl(2),SelfEne(2)
    type(z_CSR) :: HM, SM, ZM
    type(z_CSC) :: A_csc
    type(array2) :: vf(2), kz(2)
 
    integer :: i, i1, j1, j2, p, nrhs, nnz, m, iE
    integer :: n_prop_states(2), cinds(2),cinde(2), ierr, mem
    type(TStatesSummary) :: summ(2)
    real(dp) :: E, T, R, s, tmp
    complex(dp) :: Ec
    complex(dp), parameter :: minusONE = (-1.0_dp, 0.0_dp)

    character(4) :: Energy
    character(1) :: LR
    ! ---------------------------------------------------------------
    ! EXTRACT DEVICE AND CONTACT BLOCKS
    ! ---------------------------------------------------------------
    call extract_deviceHS(HH,SS,str_info,HM,SM)
    do i = 1, 2
       str_info%active_cont = i
       call extract_contactHS(HH,SS,str_info,&
            H11(i),S11(i),H12(i),S12(i),PLdim(i))
       
       write(*,*) 'PLdim',i,'=',PLdim(i)

       ! Extract Device/Contact interaction
       ! Just for the first contact PL
       cinds(i) = str_info%mat_PL_start(str_info%cblk(i))
       cinde(i) = str_info%mat_PL_end(str_info%cblk(i)) 
       j1=str_info%mat_C_start(i) 
       j2=j1+PLdim(i)-1
       call zextract(HH,cinds(i),cinde(i),j1,j2,TM(i))         
       call zextract(SS,cinds(i),cinde(i),j1,j2,ST(i))
    enddo
  
    do i = 1, 2 
       call create(Vr(i),2*PLdim(i),2*PLdim(i))
       call create(Z0(i),PLdim(i),PLdim(i))
       call create(Z1(i),PLdim(i),PLdim(i))
       allocate( vf(i)%v(2*PLdim(i)) )
       allocate( kz(i)%c(2*PLdim(i)) )
    end do

    LR = numtolabel(ni)

    open(1001,file='dos_'//LR//'.dat')
    open(1002,file='qq_'//LR//'.dat')
    open(1003,file='Trans_Refl_'//LR//'.dat')
    open(1004,file='Bloch_states_'//LR//'.dat')
    open(1005,file='Psi_C_'//LR//'.dat')
    ! ---------------------------------------------------------------
    ! BEGIN ENERGY LOOP 
    ! ---------------------------------------------------------------
    E = Erange(1)
    s = 0.0_dp
    iE = 0

    do while (E .le. Erange(2)) 
       
       iE = iE + 1
       ! ---------------------------------------------------------------
       ! Compute Complex bands for the two contacts
       ! Construct D matrices (Vr(i))
       !----------------------------------------------------------------
       do i = 1, 2

          Z0(i)%val = E*S11(i)%val - H11(i)%val
          Z1(i)%val = E*S12(i)%val - H12(i)%val
          call complex_k(E,PLdim(i),Z0(i)%val,Z1(i)%val,kz(i)%c,&
                                           Cr=Vr(i)%val,vf=vf(i)%v)
    
          call sort_and_normalize(2*PLdim(i),kz(i)%c,&
                                         vf(i)%v,Vr(i)%val,summ(i))

          !write(*,'(a,4(i3))') 'states of',i, &
          !     summ(i)%prop_in,summ(i)%evan_in,summ(i)%null_in

          if(summ(i)%prop_in.ne.summ(i)%prop_out .or. &
               summ(i)%evan_in.ne.summ(i)%evan_out ) &
               STOP 'ERROR: Asymmetry found betw. IN/OUT states'

          n_prop_states(i) = summ(i)%prop_out

       enddo

       ! Write down all k-values for propagating and evanescent states:
       ! E n_in Re(k_in) v(in) ...  n_out Re(k_out) v(out) ... n_ev_out Re(k_out) Im(k_out) ... 
       !       n_out_2 Re(k_out) v(out) ...

       write(1002,'(F15.8,i4)', advance='NO') E*HAR,summ(ni)%prop_in
       do i = 1, summ(ni)%prop_in
          write(1002,'(F15.8,F15.8)', advance='NO') real(kz(ni)%c(i)), vf(ni)%v(i)
       enddo

       write(1002,'(i4)', advance='NO') summ(ni)%prop_out
       do i = PLdim(ni)+1, PLdim(ni)+summ(ni)%prop_out
          write(1002,'(F15.8,F15.8)', advance='NO') real(kz(ni)%c(i)), vf(ni)%v(i)
       enddo 

       write(1002,'(i4)', advance='NO') summ(ni)%evan_out      
       do i = PLdim(ni)+summ(ni)%prop_out+1, PLdim(ni)+summ(ni)%prop_out+summ(ni)%evan_out
          write(1002,'(a2,F15.8,a1,F15.8,a2)', advance='NO') &
               ' (',real(kz(ni)%c(i)),',',aimag(kz(ni)%c(i)),') '
       enddo         

       write(1002,'(i4)', advance='NO') summ(nf)%prop_out
       do i = PLdim(nf)+1, PLdim(nf)+summ(nf)%prop_out
          write(1002,'(F15.8,F15.8)', advance='NO') real(kz(ni)%c(i)), vf(nf)%v(i)
       enddo       
       
       write(1002,'(i4)', advance='NO') summ(nf)%evan_out      
       do i = PLdim(nf)+summ(nf)%prop_out+1, PLdim(nf)+summ(nf)%prop_out+summ(nf)%evan_out
          write(1002,'(a2,F15.8,a1,F15.8,a2)', advance='NO') &
               ' (',real(kz(ni)%c(i)),',',aimag(kz(ni)%c(i)),') '
       enddo 
       write(1002,*) 

       !! Write down all Bloch states for propagating and evanescent states:
       !! C(in) ...  C(out) ... C(evan) ...
       
       write(1004,'(F15.8)') E*HAR 
       
       do j1 = 1, PLdim(ni)
          ! prop_in
          do i = 1, summ(ni)%prop_in      
             write(1004,'(F15.8,F15.8)', advance='NO') Vr(ni)%val(j1,i)
          enddo
          ! prop_out (riflessi)
          do i = PLdim(ni)+1, PLdim(ni)+summ(ni)%prop_out
             write(1004,'(F15.8,F15.8)', advance='NO') Vr(ni)%val(j1,i)
          enddo
          ! evan_out (riflessi)
          do i = PLdim(ni)+summ(ni)%prop_out+1, PLdim(ni)+summ(ni)%prop_out+summ(ni)%evan_out
             write(1004,'(F15.8,F15.8)', advance='NO') Vr(ni)%val(j1,i)
          enddo
          ! prop_out (trasm)
          do i = PLdim(nf)+1, PLdim(nf)+summ(nf)%prop_out
             write(1004,'(F15.8,F15.8)', advance='NO') Vr(nf)%val(j1,i)
          enddo
          ! evan_out (trasm)
          do i = PLdim(nf)+summ(nf)%prop_out+1, PLdim(nf)+summ(nf)%prop_out+summ(nf)%evan_out
             write(1004,'(F15.8,F15.8)', advance='NO') Vr(nf)%val(j1,i)
          enddo
          write(1004,*)           
       enddo
         

       if (n_prop_states(ni).eq.0) then
          !write(1000,'(3(F15.8),i4)') E*HAR,0.0_dp,0.0_dp,0
          write(1001,'(3(F15.8))') E*HAR,0.0_dp,s
          E = E + Erange(3)
          cycle
       else if(n_prop_states(nf).eq.0) then
          i1 = n_prop_states(ni)
          !write(1000,'(3(F15.8),i4)') E*HAR,0.0_dp,i1*0.0_dp,i1
          write(1001,'(3(F15.8))') E*HAR,0.0_dp,s
          E = E + Erange(3)
          cycle         
       endif

       ! ---------------------------------------------------------------
       ! Compute Surface green's functions
       !         and self-energies for a=L,R
       !----------------------------------------------------------------
       do i=1,2

          Ec = cmplx(E)   
          call prealloc_sum(TM(i),ST(i),minusONE,Ec,Tlc(i)) 
          call create(TT(i),Tlc(i)%nrow,Tlc(i)%ncol)
          call csr2dns(Tlc(i),TT(i))
          call destroy(Tlc(i))
         
          call surf_gf_(Ec,Z0(i),Z1(i),Vr(i),PLdim(i),invD12(i),GS_d(i))

          call create(Self,TT(i)%nrow,TT(i)%nrow)

          call create(TT1,GS_d(i)%nrow,TT(i)%nrow)
          TT1%val = matmul(GS_d(i)%val,conjg(transpose(TT(i)%val)))
          Self%val= matmul(TT(i)%val,TT1%val)
          call destroy(TT1)
          call create(SelfEne(i),Self%nrow,Self%ncol,nzdrop(Self,EPS))
          call dns2csr(Self,SelfEne(i))
          call destroy(Self)

          !call create(GS(i),GS_d(i)%nrow,GS_d(i)%ncol,nzdrop(GS_d(i),EPS))
          !call dns2csr(GS_d(i),GS(i))
          !call zdagacsr(Tlc(i),Tcl(i))
          !call SelfEnergy(GS(i),Tlc(i),Tcl(i),SelfEne(i))
          !call destroy(GS(i))
          !call destroy(Tcl(i))
       enddo

       ! ---------------------------------------------------
       ! Set R.H.S.
       ! ---------------------------------------------------
       nrhs = n_prop_states(ni)
       call set_rhs(PLdim(ni),nrhs,Z0(ni),Z1(ni),Vr(ni),&
                    GS_d(ni),TT(ni),RHS,HM%nrow,cinds(ni))

       ! ---------------------------------------------------
       ! ASSEMBLE MATRIX: ZM = ES-H-SelfEne(1)-SelfEne(2)
       ! ---------------------------------------------------
       call prealloc_sum(HM,SM,minusONE,Ec,ZM)
       do i=1,2
          CALL concat(ZM,minusONE,SelfEne(i),cinds(i),cinds(i))
          call destroy(SelfEne(i))
       enddo       

       ! ---------------------------------------------------       
       ! Solve linear systems of equations  (call superLU())
       ! ---------------------------------------------------
       !Trasposizione da csr a csc (needed by superLU)
       call create(A_csc,ZM%nrow,ZM%ncol,ZM%nnz)
       call csr2csc(ZM,A_csc)
       call destroy(ZM)

       call c_bridge_zgssv(A_csc%nrow, A_csc%nnz, nrhs, &
            A_csc%nzval, A_csc%rowind, A_csc%colpnt, &
            RHS%val, A_csc%nrow, ierr, mem)      

       call destroy(A_csc)

       ! ----------------- ->   R^-1 ->      R^-1    + ->  -----------
       ! Find transmission t = D     C  = - D    g  T  C
       ! -----------------      12    N      12   R  R  M  -----------
       
       if(TT(nf)%ncol.ne.PLdim(nf)) stop 'Error size!'
       
       call create(TT2,PLdim(nf),nrhs)
       TT2%val=matmul(conjg(transpose(TT(nf)%val)), &
                            RHS%val(cinds(nf):cinde(nf),1:nrhs))      

       call create(TT1,PLdim(nf),nrhs)
       TT1%val=matmul(GS_d(nf)%val,TT2%val)
       call destroy(TT2)

       call create(TR,PLdim(nf),nrhs)
       TR%val=matmul(invD12(nf)%val,TT1%val)
       call destroy(TT1)

       !------------------------------------------------------------
       ! For every real bloch state that travels towards the device 
       ! calculate the tranmission probability T:
       !      n_prop(ni)
       ! T  =    Sum     T     
       !          i       i
       !
       !      n_prop(nf)              2      |vj(E,kII,R)|
       ! T  =    Sum     | t (E,kII) |   *  ---------------
       !  i      j=1        j                |vi(E,kII,L)|
       !
       ! vi,j are the group velocities for each propagating 
       !      bloch wave 
       !------------------------------------------------------------ 
       i1 = n_prop_states(ni)
       j1 = 1
       j2 = n_prop_states(nf)
       m  = PLdim(nf)
       T=0.0_dp
       do i = 1, i1  ! For every possible incoming state
                  
          TR%val(j1:j2,i)=TR%val(j1:j2,i)*sqrt(vf(nf)%v(m+j1:m+j2))

          T=T+DOT_PRODUCT(TR%val(j1:j2,i),TR%val(j1:j2,i))/(-vf(ni)%v(i))

       end do
     
       write(1005,'(F15.8)') E*HAR 
       do i = 1, n_prop_states(ni)    
          do j2 = 1, HM%nrow     
             write(1005,'(F15.8,F15.8)',advance='NO') RHS%val(j2,i)
          enddo
          write(1005,*)
       enddo
       !------------------------------------------------------------ 
       !------------------------------------------------------------ 
       ! COMPUTE Reflection coefficient (cross check)

       call create(TT2,PLdim(ni),nrhs)
       TT2%val=matmul(conjg(transpose(TT(ni)%val)),RHS%val(cinds(ni):cinde(ni),1:nrhs))      

       call create(TT1,PLdim(ni),nrhs)
       TT1%val=matmul(GS_d(ni)%val,TT2%val)
       call destroy(TT2)

       call create(RR,PLdim(ni),nrhs)
       RR%val=matmul(invD12(ni)%val,TT1%val)
       call destroy(TT1)

       call destroy(RHS) 

       call set_rhs(PLdim(ni),nrhs,Z0(ni),Z1(ni),Vr(ni),&
                    GS_d(ni),TT(ni),RHS)

       call create(TT1,PLdim(ni),nrhs)
       TT1%val=matmul(invD12(ni)%val,RHS%val)

       RR%val = RR%val + TT1%val
       call destroy(TT1)

       i1 = n_prop_states(ni)
       j1 = 1
       j2 = n_prop_states(ni)
       m  = PLdim(ni)       

       R=0.0_dp
       do i = 1, i1  ! For every possible incoming state
                  
          RR%val(j1:j2,i)=RR%val(j1:j2,i)*sqrt(vf(ni)%v(m+j1:m+j2))

          R=R+DOT_PRODUCT(RR%val(j1:j2,i),RR%val(j1:j2,i))/(-vf(ni)%v(i))

       end do

       !------------------------------------------------------------ 
       write(*,'(a2,2(f15.8),3(i4))') &
            'E=',E*HAR,T,i1,summ(ni)%prop_in,summ(ni)%null_in 

       !write(1000,'(3(F15.8),i4)') E*HAR,T,R,j2
       !------------------------------------------------------------ 
       ! Compute D.O.S.
       !
       RHS%val(:,1)=abs(RHS%val(:,1))**2/(-vf(ni)%v(1))
       do i=2,nrhs
          RHS%val(:,1)=RHS%val(:,1)+abs(RHS%val(:,i))**2/(-vf(ni)%v(i))
       end do

       tmp = abs( sum(RHS%val(:,1),1))
       s = s + tmp * 2*Erange(3)/PI

       write(1001,'(3(F15.8))') E*HAR,tmp,s
       
       !------------------------------------------------------------ 
       ! Write r_ij, t_ij 
       write(1003,'(F15.8)') E*HAR 
       do i = 1,  n_prop_states(ni)
          do j1 = 1, summ(ni)%prop_out + summ(ni)%evan_out
             !write(1003,'(F15.8,F15.8)', advance='NO') -RR%val(j1,i)
             write(1003,'(a2,F15.8,a1,F15.8,a2)', advance='NO') &
                  ' (',real(-RR%val(j1,i)),',',aimag(-RR%val(j1,i)),') '
          enddo          
          do j1 = 1, summ(nf)%prop_out + summ(nf)%evan_out
             !write(1003,'(F15.8,F15.8)', advance='NO') -TR%val(j1,i)
             write(1003,'(a2,F15.8,a1,F15.8,a2)', advance='NO') &
                  ' (',real(-TR%val(j1,i)),',',aimag(-TR%val(j1,i)),') '

          enddo
          write(1003,*)
       enddo

       !------------------------------------------------------------ 
       ! release mem
       call destroy(RHS,TR,RR)
       do i=1,2
          call destroy(TT(i),invD12(i),GS_d(i))
       end do  

       ! update next energy and close loop 
       E = E + Erange(3)
    enddo

    do i=1,2
       call destroy(Vr(i),Z0(i),Z1(i),H11(i),H12(i),S11(i),S12(i))
       call destroy(TM(i),ST(i))
       deallocate(vf(i)%v, kz(i)%c)
    enddo

    call destroy(HM,SM)
    close(1001)
    close(1002)
    close(1003)
    close(1004)
    close(1005)

  end subroutine transfer_mat



  !----------------------------------------------------------------  
  ! Create Surface G.F. using complex band bloch-vectors
  ! Only for the first contact PL 
  !    ^-1    a    a   a   a ^-1
  !   g    = Z  + Z   D   D   
  !    a      11   12  22  12
  !                   +
  !     Self = T  g  T
  !         a   a  a  a
  ! ---------------------------------------------------------------
  subroutine surf_gf_(E,Z0,Z1,Vr,PLdim,invD12,GS_d)
    complex(dp), intent(in) :: E
    type(z_DNS), intent(in) :: Z0,Z1,Vr
    integer, intent(in) :: PLdim
    type(z_DNS), intent(out) :: invD12
    type(z_DNS), intent(out), optional :: GS_d

    complex(dp), dimension(:,:), allocatable :: TT1,TT2
    integer :: i1,i2,j1,j2,nnz

    ! Extract and invert D12
    allocate(TT1(PLdim,PLdim))
    allocate(TT2(PLdim,PLdim))
    i1 = 1; i2 = PLdim
    j1 = PLdim+1; j2 = 2*PLdim
    TT2=Vr%val(i1:i2,j1:j2)
    ! This inverse may not exist. 
    ! Probably one should invert removing the null-subspace
    call create(invD12,PLdim,PLdim)
    call inverse(invD12%val,TT2,PLdim)

    if(present(GS_d)) then
       ! Compute D22*D12^-1
       i1 = PLdim+1; i2 = 2*PLdim
       j1 = PLdim+1; j2 = 2*PLdim
       TT2 = matmul(Vr%val(i1:i2,j1:j2),invD12%val)
       ! Compute inverse (like G.F.)
       TT1 = Z0%val+matmul(Z1%val,TT2)
       
       call create(GS_d,PLdim,PLdim)   
       call inverse(GS_d%val,TT1,PLdim)
    end if

    deallocate(TT1,TT2)
        
  end subroutine surf_gf_


  ! ---------------------------------------------------------------
  ! ---------        [  L   L     L   L  ]         L
  ! set r.h.s.  T  g [ Z   D   + Z   D   ]In - T  D   In
  ! ---------    L  L[  12  21    11  11 ]      L  11
 
  subroutine set_rhs(PLdim,nrhs,Z0,Z1,Vr,GS_d,Tlc_d,RHS,nrow,cind)
    integer, intent(in) :: PLdim
    integer, intent(in) :: nrhs
    type(z_DNS), intent(in)  :: Z0,Z1,Vr,GS_d
    type(z_DNS), intent(in)  :: Tlc_d
    type(z_DNS), intent(out) :: RHS
    integer, intent(in), OPTIONAL :: nrow,cind
    
    integer :: i,i1,i2,j1,j2, nnz
    complex(dp), dimension(:,:), allocatable :: TT1,TT2,In 
    type(z_DNS) :: TT3,trhs
   
    !call log_allocate(In,nrhs,nrhs)

    !In=(0.0_dp,0.0_dp)
    ! Setup injection matrix
    !In(1:nrhs,1:nrhs) = &
    !     reshape((/(1.0_dp, (0.0_dp,i1=1,nrhs),&
    !     j1=1,nrhs-1), 1.0_dp/),(/nrhs,nrhs/));
    
    call log_allocate(TT1,PLdim,nrhs)
    call log_allocate(TT2,PLdim,nrhs)
    call create(TT3,PLdim,nrhs)
    
    ! Z12 * D21 * In:
    i1 = PLdim+1; i2 = 2*PLdim
    j1 = 1; j2 = nrhs 
    TT1 = Vr%val(i1:i2,j1:j2)        ! D21 (* In)
    TT2 = matmul(Z1%val,TT1)         ! Z12 * D21 * In
    
    ! Z11 * D11 * In:
    i1 = 1; i2 = PLdim
    j1 = 1; j2 = nrhs   
    TT1 = Vr%val(i1:i2,j1:j2)        ! D11 (*In)
    TT2 = TT2 + matmul(Z0%val,TT1)   ! Z12*D21*In + Z11*D11*In 
    
    !call log_deallocate(In)
    
    if(.not.present(nrow)) then

       ! gL * [Z11*D11*In + Z12*D21*In]        
       TT3%val = matmul(GS_d%val,TT2) 
       call log_deallocate(TT1)
       call log_deallocate(TT2)

       call create(RHS,TT3%nrow,TT3%ncol)    
       RHS%val = TT3%val

       call destroy(TT3)

       return
       
    else

       ! gL * [Z11*D11*In + Z12*D21*In] - D11*In       
       TT3%val = matmul(GS_d%val,TT2) - TT1 
       
       call log_deallocate(TT1)
       call log_deallocate(TT2)
       
       !Convert TL to a dense matrix 
       !and take the product TL * TT3 

       call create(trhs,Tlc_d%nrow,nrhs)
       trhs%val = matmul(Tlc_d%val,TT3%val)       
       call destroy(TT3)

       call create(RHS,nrow,nrhs)    
       RHS%val=(0.0_dp,0.0_dp)
       
       i2 = cind + trhs%nrow - 1
       
       !write(*,'(i4,i4,i4)')  cind, i2, trhs%nrow 
       
       RHS%val(cind:i2,1:nrhs)=trhs%val(1:i2-cind+1,1:nrhs)
       
       call destroy(trhs)
  
    endif

  end subroutine set_rhs



  !----------------------------------------------------------------   

  subroutine extract_contactHS(HH,SS,cont_info,HM,SM,TM,ST,PLdim)
    type(z_CSR) :: HH,SS
    type(TStruct_info), intent(in) :: cont_info
    type(z_DNS) :: HM,SM,TM,ST
    integer, intent(out) :: PLdim

    integer :: cont, mat_start, mat_end, err
    type(z_CSR) :: Hc,Sc
    type(z_DNS) :: Hc_d, Sc_d
 
    cont = cont_info%active_cont
    mat_start= cont_info%mat_c_start(cont)
    mat_end  = cont_info%mat_c_end(cont)
    PLdim = (mat_end - mat_start + 1)/2
    
    call zextract(HH,mat_start,mat_end,mat_start,mat_end,Hc)
    call zextract(SS,mat_start,mat_end,mat_start,mat_end,Sc)  
    
    call create(Hc_d,Hc%nrow,Hc%ncol)
    call create(Sc_d,Sc%nrow,Sc%ncol)
    
    call csr2dns(Hc,Hc_d)
    call csr2dns(Sc,Sc_d)
    call destroy(Hc)
    call destroy(Sc)

    call create(HM,PLdim,PLdim)
    call create(SM,PLdim,PLdim)
    call create(TM,PLdim,PLdim)
    call create(ST,PLdim,PLdim) 

    HM%val(1:PLdim,1:PLdim) = Hc_d%val(1:PLdim,1:PLdim)
    TM%val(1:PLdim,1:PLdim) = Hc_d%val(1:PLdim,PLdim+1:2*PLdim)
    SM%val(1:PLdim,1:PLdim) = Sc_d%val(1:PLdim,1:PLdim)
    ST%val(1:PLdim,1:PLdim) = Sc_d%val(1:PLdim,PLdim+1:2*PLdim)   
    
    call destroy(Hc_d)
    call destroy(Sc_d)  
        
  end subroutine extract_contactHS

  !----------------------------------------------------------------   

  subroutine extract_deviceHS(HH,SS,str_info,HM,SM)
    type(z_CSR) :: HH,SS
    type(TStruct_info), intent(in) :: str_info
    type(z_CSR) :: HM,SM
    
    integer :: mat_start,mat_end

    mat_start= str_info%mat_PL_start(1)
    mat_end  = str_info%mat_PL_end(str_info%num_PLs)   
    
    call zextract(HH,mat_start,mat_end,mat_start,mat_end,HM)
    call zextract(SS,mat_start,mat_end,mat_start,mat_end,SM)    

  end subroutine extract_deviceHS


  !-----------------------------------------------------------
  ! Implementation II with Klimeck trick.
  !-----------------------------------------------------------
  subroutine transfer_mat2(ni,nf,HH,SS,str_info,Erange)
    type array2
       real(dp), dimension(:), allocatable :: v
       complex(dp), dimension(:), allocatable :: c
    end type array2
    integer, intent(in) :: ni,nf
    type(z_CSR), intent(in) :: HH,SS
    type(TStruct_info) :: str_info
    real(dp), dimension(3), intent(in) :: Erange

    ! -------------------------------------------------------
    integer :: PLdim(2), Nk
    type(z_DNS) :: H11(2),S11(2),H12(2),S12(2),GS_d(2),Vr(2) 
    type(z_DNS) :: Z0(2),Z1(2),invD12(2),RHS,TT1,TT2,TT(2),TR
    type(z_DNS) :: D11,D21,D12,D22,D12h,L11,LM1(2),L1M(2),Self
    type(z_CSR) :: TM(2), ST(2), GS(2),Tlc(2),Tcl(2),SelfEne(2)
    type(z_CSR) :: HM, SM, ZM
    type(z_CSC) :: A_csc
    complex(dp), dimension(:,:), allocatable :: In
    type(array2) :: vf(2), kz(2)
 
    integer :: i, i1, j1, j2, p, nrhs, nnz, m, iE
    integer :: n_prop_states(2), cinds(2),cinde(2), ierr, mem
    type(TStatesSummary) :: summ(2)
    real(dp) :: E, T, s, tmp
    complex(dp) :: Ec

    character(4) :: Energy
    ! ---------------------------------------------------------------
    ! EXTRACT DEVICE AND CONTACT BLOCKS
    ! ---------------------------------------------------------------
    call extract_deviceHS(HH,SS,str_info,HM,SM)
    do i = 1, 2
       str_info%active_cont = i
       call extract_contactHS(HH,SS,str_info,&
            H11(i),S11(i),H12(i),S12(i),PLdim(i))
       
       ! Extract Device/Contact interaction
       ! Just for the first contact PL
       cinds(i) = str_info%mat_PL_start(str_info%cblk(i))
       cinde(i) = str_info%mat_PL_end(str_info%cblk(i)) 
       j1=str_info%mat_C_start(i) 
       j2=j1+PLdim(i)-1
       call zextract(HH,cinds(i),cinde(i),j1,j2,TM(i))         
       call zextract(SS,cinds(i),cinde(i),j1,j2,ST(i))
    enddo
  
    do i = 1, 2 
       call create(Vr(i),2*PLdim(i),2*PLdim(i))
       call create(Z0(i),PLdim(i),PLdim(i))
       call create(Z1(i),PLdim(i),PLdim(i))
       allocate( vf(i)%v(2*PLdim(i)) )
       allocate( kz(i)%c(2*PLdim(i)) )
    end do


    ! ---------------------------------------------------------------
    ! BEGIN ENERGY LOOP 
    ! ---------------------------------------------------------------
    E = Erange(1)
    iE = 1
    s = 0.0_dp
    do while (E .le. Erange(2)) 

       iE = iE + 1
       ! ---------------------------------------------------------------
       ! Compute Complex bands for the two contacts
       ! Construct D matrices (Vr(i))
       !----------------------------------------------------------------
       do i = 1, 2

          Z0(i)%val = E*S11(i)%val - H11(i)%val
          Z1(i)%val = E*S12(i)%val - H12(i)%val
          call complex_k(E,PLdim(i),Z0(i)%val,Z1(i)%val,kz(i)%c,&
                                           Cr=Vr(i)%val,vf=vf(i)%v)
    
          call sort_and_normalize(2*PLdim(i),kz(i)%c,&
                                         vf(i)%v,Vr(i)%val,summ(i))

          !write(*,'(a,4(i3))') 'states of',i, &
          !     summ(i)%prop_in,summ(i)%evan_in,summ(i)%null_in

          if(summ(i)%prop_in.ne.summ(i)%prop_out .or. &
               summ(i)%evan_in.ne.summ(i)%evan_out ) &
               STOP 'ERROR: Asymmetry found betw. IN/OUT states'

          n_prop_states(i) = summ(i)%prop_out

       enddo
       

       if (n_prop_states(ni).eq.0) then
          !write(1000,'(3(F15.8),i4)') E*HAR,0.0_dp,0.0_dp,0
          !write(1001,'(3(F15.8))') E*HAR,0.0_dp,0.0_dp
          E = E + Erange(3)
          cycle
       else if(n_prop_states(nf).eq.0) then
          i1 = n_prop_states(ni)
          !write(1000,'(3(F15.8),i4)') E*HAR,0.0_dp,i1*0.0_dp,i1
          !write(1001,'(3(F15.8))') E*HAR,0.0_dp,0.0_dp
          E = E + Erange(3)
          cycle         
       endif

       ! ---------------------------------------------------------------
       ! Compute Surface green's functions
       !         and self-energies for a=L,R
       !----------------------------------------------------------------
       p=+1
       if (nf.gt.ni) p=-1;

       do i=nf,ni,p   ! This trick is to compute ni last

          Ec = cmplx(E)   
          call prealloc_sum(TM(i),ST(i),minusONE,Ec,Tlc(i))
          call create(TT(i),Tlc(i)%nrow,Tlc(i)%ncol)
          call csr2dns(Tlc(i),TT(i))
          call destroy(Tlc(i))

          ! Klimeck implementation (trial):
          ! KLIMECK: L11 = D12h * (Z11*D12 + Z12*D22)
          !          LM1 = Tlc * D12;  L1M = LM1+
          !          Self= LM1 * L11^-1 * L1M

          i1 = PLdim(i)
          m = PLdim(i)
          call create(D12,m,i1)
          call create(D22,m,i1)
          call create(D12h,i1,m)
          call create(L11,i1,i1)   
          call create(GS_d(i),i1,i1)
          call create(LM1(i),TT(i)%nrow,m)
          call create(L1M(i),m,TT(i)%nrow)
          call create(D11,m,i1)
 
          D12%val = Vr(i)%val(1:m,m+1:m+i1)       !D12
          D22%val = Vr(i)%val(m+1:2*m,m+1:m+i1)   !D22      
          D12h%val = conjg(transpose(D12%val))    !D12+

          LM1(i)%val = matmul(TT(i)%val,D12%val)
          L1M(i)%val = conjg(transpose(LM1(i)%val))
    
          D11%val=matmul(Z0(i)%val,D12%val)+matmul(Z1(i)%val,D22%val)
          L11%val=matmul(D12h%val,D11%val)

          call inverse(GS_d(i)%val,L11%val,i1)
          
          call create(Self,TT(i)%nrow,TT(i)%nrow)
          Self%val = matmul(LM1(i)%val,matmul(GS_d(i)%val,L1M(i)%val))
          call create(SelfEne(i),Self%nrow,Self%ncol,nzdrop(Self,EPS))
          call dns2csr(Self,SelfEne(i))
          call destroy(Self)

          call destroy(D11,D12,D22,L11)

          ! this matrix is needed for r.h.s.
          if (i.ne.ni) then
             call destroy(D12h)                  
          endif

       enddo

       ! ---------------------------------------------------
       ! Set R.H.S.
       ! ---------------------------------------------------
       ! KLIMECK: I1 = D12h * (Z11*D11 + Z12*D21) * In
       !          I  = LM1*L11^-1 *I1  - Tlc*D11*In
       nrhs = n_prop_states(ni)
       m = PLdim(ni)

       call create(D11,m,nrhs)
       call create(D21,m,nrhs)
       call create(D22,m,nrhs)
       call create(D12,m,nrhs)
       call create(TT1,LM1(ni)%nrow,nrhs)
       call create(TT2,LM1(ni)%nrow,nrhs)

       D11%val = Vr(ni)%val(1:m,1:nrhs)         
       D21%val = Vr(ni)%val(m+1:2*m,1:nrhs)     

       ! TL * D11 * In  
       TT1%val = matmul(TT(ni)%val,D11%val)
       ! Z11* D11 * In +  Z12 * D21 * In
       D22%val=matmul(Z0(ni)%val,D11%val) + matmul(Z1(ni)%val,D21%val)
       ! I1 =  D12h * (Z11*D11 + Z12*D21) * In
       D12%val = matmul(D12h%val,D22%val)
       ! LM1* L11^-1 * I1
       D21%val = matmul(GS_d(ni)%val,D12%val)
       TT2%val = matmul(LM1(ni)%val,D21%val)

       call create(RHS,HM%nrow,nrhs)

       if(TT1%nrow.ne.cinde(ni)-cinds(ni)+1) stop 'wrong size'

       RHS%val = (0.0_dp,0.0_dp)

       RHS%val(cinds(ni):cinde(ni),1:nrhs)= TT2%val-TT1%val 

       !call log_deallocate(In)
       call destroy(D12,D11,D21,D22,TT1,TT2,D12h)
       ! ---------------------------------------------------
       ! ASSEMBLE MATRIX: ZM = ES-H-SelfEne(1)-SelfEne(2)
       ! ---------------------------------------------------
       call prealloc_sum(HM,SM,minusONE,Ec,ZM)
       do i=1,2
          CALL concat(ZM,minusONE,SelfEne(i),cinds(i),cinds(i))
          call destroy(SelfEne(i))
       enddo       

       ! ---------------------------------------------------       
       ! Solve linear systems of equations  (call superLU())
       ! ---------------------------------------------------
       !Trasposizione da csr a csc (needed by superLU)
       call create(A_csc,ZM%nrow,ZM%ncol,ZM%nnz)
       call csr2csc(ZM,A_csc)
       call destroy(ZM)

       call c_bridge_zgssv(A_csc%nrow, A_csc%nnz, nrhs, &
            A_csc%nzval, A_csc%rowind, A_csc%colpnt, &
            RHS%val, A_csc%nrow, ierr, mem)      

       call destroy(A_csc)

       ! ----------------- ->   R^-1 ->    R^-1    + ->  -----------
       ! Find transmission t = D     C  = D    g  T  C
       ! -----------------      12    N    12   R  R  M  -----------
       
       ! KLIMECK: t = - L11^-1 * L1M * CM
       !                (j2,j2) * (j2,nrow) * (nrow,nrhs) 

       call create(TT1,L1M(nf)%nrow,nrhs)
       call create(TR,PLdim(nf),nrhs)

       TT1%val = matmul(L1M(nf)%val,RHS%val(cinds(nf):cinde(nf),1:nrhs))
       TR%val = matmul(GS_d(nf)%val, TT1%val)

       call destroy(TT1)
       !------------------------------------------------------------
       ! For every real bloch state that travels towards the device 
       ! calculate the tranmission probability T:
       !      n_prop(ni)
       ! T  =    Sum     T     
       !          i       i
       !
       !      n_prop(nf)              2      |vj(E,kII,R)|
       ! T  =    Sum     | t (E,kII) |   *  ---------------
       !  i      j=1        j                |vi(E,kII,L)|
       !
       ! vi,j are the group velocities for each propagating 
       !      bloch wave 
       !------------------------------------------------------------ 
       i1 = n_prop_states(ni)
       j1 = 1
       j2 = n_prop_states(nf)
       m  = PLdim(nf)
       T=0.0_dp
       do i = 1, i1  ! For every possible incoming state
                  
          TR%val(j1:j2,i)=TR%val(j1:j2,i)*sqrt(vf(nf)%v(m+j1:m+j2))

          T=T+DOT_PRODUCT(TR%val(j1:j2,i),TR%val(j1:j2,i))/(-vf(ni)%v(i))

       end do

       !------------------------------------------------------------ 
       write(*,'(a2,2(f15.8),3(i4))') &
            'E=',E*HAR,T,i1,summ(ni)%prop_in,summ(ni)%null_in 

       !write(1000,'(2(F15.8),i4)') E*HAR,T,j2
       !------------------------------------------------------------ 
       ! Compute D.O.S.
       !
       RHS%val(:,1)=abs(RHS%val(:,1))**2/(-vf(ni)%v(1))
       do i=2,nrhs
          RHS%val(:,1)=RHS%val(:,1)+abs(RHS%val(:,i))**2/(-vf(ni)%v(i))
       end do

       tmp = abs( sum(RHS%val(:,1),1))
       s = s + tmp * 2*Erange(3)/PI

       !write(1001,'(3(F15.8))') E*HAR,tmp,s


       !------------------------------------------------------------ 
       ! release mem
       call destroy(RHS,TR)
       do i=1,2
          call destroy(LM1(i),L1M(i),TT(i),GS_d(i))
       end do  

       ! update next energy and close loop 
       E = E + Erange(3)
    enddo

    do i=1,2
       call destroy(Vr(i),Z0(i),Z1(i),H11(i),H12(i),S11(i),S12(i))
       call destroy(TM(i),ST(i))
       deallocate(vf(i)%v,kz(i)%c)
    enddo

    call destroy(HM,SM)

  end subroutine transfer_mat2

  !-----------------------------------------------------------
  ! Search of bound states
  !-----------------------------------------------------------
  
  function numtolabel(ni) result(LorR)
    integer :: ni
    character(1) :: LorR

    if(ni.eq.1) LorR='L'
    if(ni.eq.2) LorR='R'

  end function numtolabel


end module scattstates
