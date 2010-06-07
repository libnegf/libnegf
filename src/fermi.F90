module fermi_dist

  use precision
  
  implicit none
  private

  public :: fermi_f, fermi_fc

contains
  
  real(kind=dp) function fermi_f(E,Ef,kT) 

    implicit none

    real(kind=dp), intent(in) :: E, Ef, kT

    ! the check over 0 is important otherwise the next fails
    if (kT.eq.0.d0) then
      if(E.gt.Ef) then
        fermi_f = 0.D0
      else
        fermi_f = 1.D0
      end if
      return
    endif

    if (abs((E-Ef)/kT).gt.30) then
      if(E.gt.Ef) then
        fermi_f = 0.D0
      else
        fermi_f = 1.D0
      end if
      return
    else        
      fermi_f = 1.d0/(dexp((E-Ef)/kT)+1.d0);
      return
    endif

  end function fermi_f


  complex(kind=dp) function fermi_fc(Ec,Ef,kT)

    implicit none
    
    complex(kind=dp), intent(in) :: Ec
    real(kind=dp), intent(in) :: Ef, kT

    complex(kind=dp) :: Efc,kTc,ONE=(1.d0,0.d0),j=(0.d0,1.d0)

    Efc=Ef*ONE
    kTc=kT*ONE

    if (kT.eq.0.d0) then
      if(dreal(Ec).gt.Ef) then
        fermi_fc = (0.D0,0.D0)
      else
        fermi_fc = (1.D0,0.D0)
      end if
      return
    endif

    if (abs( (real(Ec)-Ef)/kT ).gt.30) then
      if(dreal(Ec).gt.Ef) then
        fermi_fc = (0.D0,0.D0)
      else
        fermi_fc = (1.D0,0.D0)
      end if
      return
    else        

      fermi_fc = ONE/(exp( (Ec-Efc)/kTc ) + ONE);
      return
    endif

  end function fermi_fc


  subroutine FERMI(nel,telec,ndim,ev,dacc,occ,efermi)

    implicit REAL*8 (A-H,O-Z)
    
    !
    integer ndim
    real*8 nel,telec,efermi,dacc,ev(*),occ(*)
    !
    ! Boltzmann constant in H/K
    !
    integer, parameter :: ckbol = 3.16679d-6
    !
    ! degeneracy tolerance for T = 0 
    !
    integer, parameter :: degtol = 1.0d-4 
    integer, parameter :: maxifermi = 1000


    real*8 beta,etol,occdg,chleft,ef1,ef2,ceps,eeps,charge,fac,racc
    integer i,nef1,nef2,nup,ndown,nocc2,ndeg,istart,iend,ifermi
    logical tzero

    ! machine accuracy
    racc=get_machineacc()
    ! start defining occupation numbers and their derivatives 
    !
    do i = 1,ndim 
      occ(i) = 0.d0
    enddo
    if (nel.lt.1.0d-5) goto 10
    if (nel.gt.2*ndim) then
      print *,'too many electrons'
      stop
    endif

    ! find energy interval for Fermi energy
    !
    if (telec.gt.5.0) then
      beta = 1.d0/(ckbol*telec) 
      etol = ckbol*telec*(log(beta)-log(racc)); tzero = .false.
    else 
      etol = degtol; tzero = .true.
    endif

    ! find energy range for Fermi energy
    !
    if(nel.gt.int(nel)) then
      nef1 = (nel+2)/2; nef2 = (nel+2)/2  
    else 
      nef1 = (nel+1)/2; nef2 = (nel+2)/2
    endif

    efermi = 0.5d0*(ev(nef1)+ev(nef2))
    nup = nef1; ndown = nef1;
    do while (nup.lt.ndim) 
      if (abs(ev(nup+1)-efermi) > etol) exit
      nup = nup+1
    enddo
    do while (ndown.gt.0)
      if (abs(ev(ndown)-efermi) > etol) exit
      ndown = ndown-1
    enddo
    ndeg = nup-ndown; nocc2 = ndown;
    do i = 1,nocc2
      occ(i) = 2.d0
    enddo
    if (ndeg == 0) goto 10

    ! for T == 0 occupy orbitals as usual
    !
    if (tzero) then
      occdg = ndeg; occdg = (nel-2*nocc2)/occdg
      do i = nocc2+1,nocc2+ndeg 
        occ(i) = occdg
      enddo
      ! for finite T, use Fermi distribution
      !
    else 
      chleft = nel-2*nocc2
      istart = nocc2+1; iend = istart+ndeg-1
      if (ndeg == 1) then
        occ(istart) = chleft
        goto 10
      endif

      ! bracket and iterate Fermi energy by bisection
      !
      ef1 = efermi-etol-degtol; ef2 = efermi+etol+degtol; 
      ceps = dacc*chleft; eeps = dacc*max(abs(ef1),abs(ef2))
      ifermi = 0
      charge = chleft+3*ceps

      do while ((abs(charge-chleft).gt.ceps).and.(abs(ef1-ef2).gt.eeps))

        efermi = 0.5*(ef1+ef2); charge = 0.d0;
        do i = istart,iend 
!!          occ(i) = 2.d0/(1.d0+exp(beta*(ev(i)-efermi)))
          charge = charge+occ(i)
        enddo
        if (charge.gt.chleft) then 
          ef2 = efermi 
        else 
          ef1 = efermi
        endif

        ifermi = ifermi + 1
        if (ifermi.gt.maxifermi) then
          write(*,*) 'fermi(): iteration stopped at charge, efermi:',charge, efermi
          exit
        endif

      enddo
      ! rescale occ(i) if the accuracy was not good enough 
      !
      if (abs(charge-chleft).gt.ceps) then
        fac = chleft/charge;
        do i = istart,iend 
          occ(i) = occ(i)*fac
        enddo
      endif

    endif

10  continue

  end subroutine FERMI
  
end module fermi_dist
