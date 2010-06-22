program test1

  use precision
  use parameters
  use allocation
  use structure
  use mat_def
  use load
  use greendftb

  implicit none

  type(Tstruct_info) :: str
  type(Tparam) :: param
  type(z_CSR) :: H, S, DensMat, EnMat

  integer :: ncont, nbl
  integer, dimension(:), allocatable :: PL_end, cont_end, surf_end, cblk
  
  integer :: ii, ka, jj, kb, jcol, nrow
  complex(dp) :: dd
  real(dp) :: qtot
  real(dp), dimension(:), allocatable :: qmulli

  call init_defaults(param)
  write(*,*) '(test) loading H,S...'

  call load_HS(H,S,.false.)

  write(*,*) '(test) H%nrow=',H%nrow
  write(*,*) '(test) H%ncol=',H%ncol
  write(*,*) '(test) H%nnz=',H%nnz
  write(*,*) '(test) S%nrow=',S%nrow
  write(*,*) '(test) S%ncol=',S%ncol
  write(*,*) '(test) S%nnz=',S%nnz

  read(*,*) ncont
  read(*,*) nbl
  
  call log_allocate(PL_end,nbl)
  call log_allocate(cblk,ncont)
  call log_allocate(cont_end,ncont)
  call log_allocate(surf_end,ncont)

  read(*,*) PL_end(1:nbl)
  read(*,*) cblk(1:ncont)
  read(*,*) cont_end(1:ncont)
  read(*,*) surf_end(1:ncont)
  if (ncont.eq.0) then
     read(*,*) param%Efermi(1)
  else
     read(*,*) param%Efermi(1:ncont)
  endif

  param%Efermi = param%Efermi/param%hartree

  read(*,*) param%Np(1:3)
  read(*,*) param%Temp
  read(*,*) param%nPoles
 
  !write(*,*) PL_end
  !write(*,*) cblk
  !write(*,*) cont_end
  !write(*,*) surf_end

  param%verbose=91

  write(*,*) '(test) create struct...'
  
  call create_Tstruct(ncont, nbl, PL_end, cont_end, surf_end, cblk, str)

  call print_Tstruct(str)

  print *, 'create DensMat, dim',str%total_dim

  call create(DensMat,str%total_dim,str%total_dim,str%total_dim)
  DensMat%nzval=0.d0

  write(*,*) '(test) start integration...'  

  call contour_int(H,S,param,str,DensMat,EnMat)

  ! -------------------------------------------------------------
  nrow = DensMat%nrow
  call log_allocate(qmulli,nrow)
  qmulli=0.d0

  do ii=1, nrow 
     do ka=DensMat%rowpnt(ii), DensMat%rowpnt(ii+1)-1 
        dd = DensMat%nzval(ka)
        jj = DensMat%colind(ka)
        
        do kb=S%rowpnt(jj),S%rowpnt(jj+1)-1
           jcol = S%colind(kb)
           if (jcol .eq. ii) then
              qmulli(jcol) = qmulli(jcol) + real(dd*S%nzval(kb))
           endif
        enddo
     enddo
  enddo
  
  open(11,file='qmulli.dat')
  qtot = 0.d0
  do ii = 1, str%central_dim
     write(11,*) ii,qmulli(ii)
     qtot = qtot+qmulli(ii)
  enddo
  close(11)

  write(*,*) 'qtot=',qtot
  write(*,*) 'should be',1.d0*str%central_dim 

  call destroy(H,S,DensMat)
  call log_deallocate(qmulli)
  call log_deallocate(PL_end)
  call log_deallocate(cblk)
  call log_deallocate(cont_end)
  call log_deallocate(surf_end)

  call writeMemInfo(6)
  call writePeakInfo(6)

end program test1


