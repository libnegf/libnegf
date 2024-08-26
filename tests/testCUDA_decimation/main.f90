program test1
  use ln_precision
  use libnegf
  use random
  use contselfenergy, only : decimation_gpu 
  implicit none
  complex(dp), dimension(:,:), allocatable :: H00, H01
  complex(dp), dimension(:,:), allocatable :: S00, S01
  complex(dp), dimension(:,:), allocatable, target :: A0, B0, C0, G0
  character(10) :: arg
  integer :: N, ncyc, fu, t1, t2, cr, cm, stat
  real(dp) :: delta, En
  complex(dp) :: z
  logical :: readfiles

  type(TNegf) :: negf

  readfiles=.true.
  if (command_argument_count() > 0) then
    call get_command_argument(1,arg)
    read(arg,*) N
    readfiles = .false.
  end if

  if (readfiles) then
    open(newunit=fu,file='size.dat',form='formatted', action='read')
    read(fu,*) N
    read(fu,*) En
    read(fu,*) delta
    close(fu)
    write(*,*) 'size:',N,'x',N
    write(*,*) 'delta:',delta
  end if

  allocate(H00(N,N))
  allocate(H01(N,N))
  allocate(S00(N,N))
  allocate(S01(N,N))

  if (readfiles) then
    open(newunit=fu,file='HC11.1.dat',access='stream',form='unformatted', action='read')
    read(fu) H00
    close(fu)
    open(newunit=fu,file='HC12.1.dat',access='stream',form='unformatted', action='read')
    read(fu) H01
    close(fu)
    open(newunit=fu,file='SC11.1.dat',access='stream',form='unformatted', action='read')
    read(fu) S00
    close(fu)
    open(newunit=fu,file='SC12.1.dat',access='stream',form='unformatted', action='read')
    read(fu) S01
    close(fu)
  else
    call zinit_rnd_matrix(H00)
    call zinit_rnd_matrix(H01)
    call zinit_rnd_matrix(S00)
    call zinit_rnd_matrix(S01)
    En = -0.1_dp
    delta = 1e-4_dp
  end if

  z = cmplx(En, delta, dp)

  allocate(A0(N,N))
  allocate(B0(N,N))
  allocate(C0(N,N))
  allocate(G0(N,N))

  A0 = z*S00 - H00
  B0 = z*S01 - H01
  C0 = z*conjg(transpose(S01))-conjg(transpose(H01))

  ! create cuBLAS instance
  write(*,*) 'Initialize negf (cuBLAS, cuSOLVER)'
  call init_negf(negf)

  call system_clock(t1, cr, cm)
  call decimation_gpu(negf, G0, A0, B0, C0, 2, .true., ncyc)
  call system_clock(t2, cr, cm)

  write(*,*) "decimation_cpu time: ",(t2-t1)*1.0/cr,"sec"
  print*,'converged in',ncyc,' iterations'
  
  deallocate(A0,B0,C0,G0)
  deallocate(H00,S00,H01,S01)

end program test1
