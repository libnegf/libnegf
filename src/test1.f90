program test1

  use precision
  use parameters
  use allocation
  use structure
  use mat_def
  use load
  use population
  use greendftb

  implicit none

  type(Tstruct_info) :: str
  type(Tparam) :: param
  type(z_CSR) :: H, S, DensMat, EnMat

  integer :: ncont, nbl
  integer, dimension(:), allocatable :: PL_end, cont_end, surf_end, cblk
  
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

  param%verbose=51

  write(*,*) '(test) create struct...'
  
  call create_Tstruct(ncont, nbl, PL_end, cont_end, surf_end, cblk, str)

  call print_Tstruct(str)

  print *, 'create DensMat, dim',str%total_dim

  call create(DensMat,str%total_dim,str%total_dim,str%total_dim)
  DensMat%nzval=0.d0

  write(*,*) '(test) start integration...'  

  call contour_int(H,S,param,str,DensMat,EnMat)

  call log_allocate(qmulli,str%central_dim)

  call mulliken(S,DensMat,qmulli)

  call destroy(H,S,DensMat)
  call log_deallocate(qmulli)
  call log_deallocate(PL_end)
  call log_deallocate(cblk)
  call log_deallocate(cont_end)
  call log_deallocate(surf_end)

  call writeMemInfo(6)

end program test1


