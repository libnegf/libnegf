program testdosdns

  use precision
  use constants
  use mat_def
  use allocation
  use structure
  use input_output
  use sparsekit_drv
  use inversions
  use iterative_dns
  use libnegf
  use lib_param
  use clock
  use extract
use iterative
  implicit none

  Type(Tnegf), target :: negf
  Type(Tnegf), pointer :: pnegf
  Integer :: outer, N, k, i
  Integer :: it, ncont, i1, nbl
  Real(dp) :: ncyc
  Complex(dp) :: Ec
  Type(z_CSR) :: Gr
  Type(z_CSR), Dimension(MAXNCONT) :: SelfEneR, Tlc, Tcl, GS


  pnegf => negf

  print*,'(main) testdosdns'

  negf%file_re_H='H_real2.dat'
  negf%file_im_H='H_imm2.dat'
  !negf%file_re_S='S_real.dat'
  !negf%file_im_S='S_imm.dat'
  negf%file_struct='driver'
  negf%verbose = 10
  negf%eneconv = 1.d0  !HAR  ! to convert Kb 
  negf%isSid = .true.  
  negf%form%formatted = .true.
  negf%form%type = 'PETSc'
  negf%form%fmt = 'F'
  negf%ReadOldSGF = 2

  print*,'(main) init'

  call init_negf(pnegf)

  print*, '(main) extract device'
  
  call extract_device(pnegf)

  print*, '(main) extract contacts'

  call extract_cont(pnegf)

  print*,'(main) computing dos'
  
  outer = 1
  it = negf%iteration
  nbl = negf%str%num_PLs
  ncont = negf%str%num_conts


  N = nint((negf%Emax-negf%Emin)/negf%Estep)

print*,'N=',N
print*,'delta=',negf%delta

print*,'(main) open dos.dat'

  open(101,file='dos.dat')

  do i=1,N

     Ec=(negf%Emin+i*negf%Estep)+negf%delta*(0.d0,1.d0)

     call compute_contacts(Ec,pnegf,it,ncyc,Tlc,Tcl,SelfEneR,GS)

     call calls_eq_mem_dns(negf%HM,negf%SM,Ec,SelfEneR,Tlc,Tcl,GS,Gr,negf%str,outer)

       do i1=1,ncont
          call destroy(Tlc(i1),Tcl(i1),SelfEneR(i1),GS(i1))
       enddo

       negf%dos = -aimag( trace(Gr) )/pi

       write(101,*) real(Ec), negf%dos
!       print*, real(Ec), negf%dos, Gr%nnz

       call destroy(Gr)

    enddo

    close(101)   


  print*,'(main) destroy negf'

  call destroy_negf(pnegf)

  call writepeakinfo(6)
  call writememinfo(6)


end program testdosdns


! negf:0 XXXXXXXXXXXX  negf:END
! pointer 
