!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Fortran 77 style subroutines for communication with UPTIGHT LIB
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! The order you should invoke the subroutines during initialisation is the
!! following:
!!
!!   negf_initsession
!!   negf_fillbasicparameters (after upt_initsession)
!!   negf_inituptight         (after upt_fillbasicparameters)
!!
!! After initialisation you can call the following methods:
!!
!!   negf_createhamiltonian    
!!   negf_lanczos              (after upt_createhamiltonian)
!!   negf_getzcsrhamiltonian   (after upt_createhamiltonian)
!!
!! In order to destroy a UPT instance call
!!
!!   negf_destructuptight
!!   negf_destructsession      (after upt_destructuptight)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!* Returns the size of the handler array
!!* @param  handlerSize  Contains the size of the handler array on exit.
subroutine negf_gethandlersize(handlerSize)
  use uptightAPICommon  ! if:mod:use
  implicit none
  integer :: handlerSize  ! if:var:out

  handlerSize = DAC_handlerSize

end subroutine negf_gethandlersize

!!* Initialises a new UPTIGHT instance
!!* @param  handler  Contains the handler for the new instance on return
subroutine negf_initsession(handler)
  use uptightAPICommon  ! if:mod:use
  use libnegf, only : negf_nullify_all
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:out
  
  type(TNEGF), pointer :: pNEGF

  IF( size(transfer(pNEGF, handler)) > size(handler) ) stop 'Handler size mismatch'
 
  !! Allocate resources
  !NULLIFY(pUPTs%pUPTIn)
  !ALLOCATE(pUPTs%pUPTIn)

  NULLIFY(pNEGF)
  ALLOCATE(pNEGF)

  handler(:) = 0
  handler = transfer(pNEGF, handler, size(handler))

  call negf_nullify_all(pNEGF)

end subroutine negf_initsession

!!* Destroys a certain UPTIGHT instance
!!* @param  handler  Handler for the instance to destroy
subroutine negf_destructsession(handler)
  use uptightAPICommon  ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:inout

  type(TNEGF), pointer :: pNEGF
  integer :: err

  pNEGF = transfer(handler, pNEGF)

  deallocate(pNEGF, stat= err)

  if (err.ne.0) write(*,*) '(negf_destructsession) Deallocation error'

end subroutine negf_destructsession


!!* Initialises a new UPTIGHT instance
!!* @param  handler  Contains the handler for the new instance on return
!subroutine negf_getversion(handler)
!  use uptightAPICommon                 !if:mod:use
!  use uptight, only : negf_version      !if:mod:use
!  implicit none
!  integer :: handler(DAC_handlerSize)  ! if:var:in
!
!  type(UPTPointers) :: pUPTs
!
!  pUPTS = transfer(handler, pUPTs)
!
!  call negf_version(pUPTS%pUPT)
!
!end subroutine negf_getversion

!!* Intialises a specific container for UPTIGHT input parameters.
!!* @param handler                : Handler of the instance.
!!* @param verbose_lev            : verbose level 
!!* @param database_path (negf_LC) : where to find the TB material database
!!* @param work_path     (negf_LC) : current working directory
!!* @param gen_filename  (negf_MC) : filename of input structure to load
!!* @param gen_out       (negf_MC) : output filename
!!* @param sparse_fmt    (negf_SC) : sparse format
!!* @param max_n_n                : to be initialized (usually 1 for n.n. tb)
!!* @param harrison_flag          : harrison scaling of TB parameters 1 = yes 
!!* @param relat_flag             : relativistic calculations         1 = yes
!!* @param potential_flag         : include potential shift           1 = yes
!!* @param optmat_flag            : computes optical matrix elements  1 = yes
!!* @param poldir                 : field polarization direction
!!* @param check_bondmap          : check bondmap internally
subroutine negf_fillbasicparameters(handler, verbose_lev, databasePath, workPath, outPath, &
                                  gen_filename, gen_outname, sparse_fmt, max_n_n, &
                                  harrison_flag, relat_flag, potential_flag, optmat_flag, & 
                                  poldir, c_axis_x, c_axis_y, c_axis_z, check_bondmap, &
                                  dg_coupl_scale, dg_onsite, hybrid_passivation)
  use precision         ! if:mod:use
  use globals           ! if:mod:use
  use uptightAPICommon  ! if:mod:use
  use negf_param         ! if:mod:use
  use uptight           ! if:mod:use
  implicit none
  integer :: handler(DAC_handlerSize)   ! if:var:in
  integer :: verbose_lev                ! if:var:in 
  character(LST) :: databasePath(1)     ! if:var:in
  character(LST) :: workPath(1)         ! if:var:in
  character(LST) :: outPath(1)          ! if:var:in
  character(MST) :: gen_filename(1)     ! if:var:in
  character(MST) :: gen_outname(1)      ! if:var:in
  character(MST) :: sparse_fmt(1)       ! if:var:in
  integer :: max_n_n                    ! if:var:in
  integer :: harrison_flag              ! if:var:in
  integer :: relat_flag                 ! if:var:in
  integer :: potential_flag             ! if:var:in
  integer :: optmat_flag                ! if:var:in
  integer :: poldir                     ! if:var:in
  integer :: check_bondmap              ! if:var:in
  integer :: hybrid_passivation         ! if:var:in
  real(dp) :: c_axis_x                  ! if:var:in  
  real(dp) :: c_axis_y                  ! if:var:in
  real(dp) :: c_axis_z                  ! if:var:in
  real(dp) :: dg_coupl_scale            ! if:var:in  
  real(dp) :: dg_onsite                 ! if:var:in  

  type(UPTPointers) :: pUPTs
  !type(TUPTIn), pointer :: pUPTIn
  type(OUPT), pointer :: pUPT
    

  pUPTs = transfer(handler, pUPTs)
  !pUPTIn => pUPTs%pUPTIn
  pUPT => pUPTs%pUPT

 
  !! Fill in parameters
  pUPT%database_path = trim( databasePath(1) ) // '/'
  pUPT%work_path = trim( workPath(1) ) // '/' 
  pUPT%gen_filename = trim(pUPT%work_path) // trim(gen_filename(1))
  pUPT%gen_out =  trim(gen_outname(1))   
  pUPT%state_file  = "states.info"
  select case(trim(sparse_fmt(1)))
  case('upper')
     pUPT%sparse_format = 'U'
  case('lower')
     pUPT%sparse_format = 'L'
  case('full')
     pUPT%sparse_format = 'F'
  end select
  pUPT%verbose = verbose_lev

  !! Initialize globals

  !write(*,*) '(upt-debug) db_path: ', trim(database_path)
  !write(*,*) '(upt-debug) wk_path: ', trim(work_path)
  !write(*,*) '(upt-debug) gen: ', trim(pUPTIn%gen_filename)

  !! Fill in UPT parameters

  pUPT%out_path = trim( outPath(1) ) // '/'
  !pUPT%nn_list%max_near = 2
  pUPT%potential_flag = .false.
  pUPT%relat = .false.
  pUPT%scaling = .false.
  pUPT%optmat = .false.
  pUPT%check_bondmap = .false.
  pUPT%hybrid_passivation = .false.
  if(relat_flag.eq.1) pUPT%relat = .true.
  if(harrison_flag.eq.1) pUPT%scaling = .true.
  if(optmat_flag.eq.1)  pUPT%optmat = .true.
  if(potential_flag.eq.1)  pUPT%potential_flag = .true.
  if(check_bondmap.eq.1)  pUPT%check_bondmap = .true.
  if(hybrid_passivation.eq.1)  pUPT%hybrid_passivation = .true.

  !write(*,*) '(upt-debug) Harrison scaling: ', pUPT%fuzzy
  !write(*,*) '(upt-debug) Optical Matrix: ', pUPT%optmat

  pUPT%d_onsite_shift_flag = .false.
  pUPT%syst_rotated = .false.
  pUPT%ioutput_flag = .false.
  pUPT%poldir = poldir    
  pUPT%nn_list%max_near = max_n_n

  pUPT%c_axis = (/ c_axis_x, c_axis_y, c_axis_z /)

  !write(*,*) '(upt-debug) c-axis: ', pUPT%c_axis

  pUPT%d_H = 1.d0/sqrt(dg_coupl_scale)
  pUPT%E_H = dg_onsite

  ! nullify some pointers

  NULLIFY(pUPT%materials)
  NULLIFY(pUPT%pot_data)
  NULLIFY(pUPT%eigen_values)
  NULLIFY(pUPT%eigen_vectors) 

  ! initialize n_basis for following checks
  pUPT%basis%n_basis = 0

  !set defaults for output that can be overridden wid setoutput
  pUPT%grid_step = 0.5d0
  pUPT%out_format = 'jvxl'


end subroutine negf_fillbasicparameters

!!* Set parameters for representation of output states
!!* @param handler
!!* @param format
!!* @param step
subroutine negf_setoutput(handler, out_format, step)
  use uptightAPICommon                   ! if:mod_use
  use globals                            ! if:mod:use
  use precision                          ! if:mod:use
  use negf_param                          ! if:mod:use
  implicit none  
  integer :: handler(DAC_handlerSize)     ! if:var:in
  character(SST) :: out_format(1)         ! if:var:in
  real(dp) :: step                        ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT

  pUPT%grid_step = step
  pUPT%out_format = trim(out_format(1))

end subroutine negf_setoutput

!!* Set the library to assemble the optical matrix
!!* @param handler                : Handler of the instance.
!!* @param optmat_flag            : computes optical matrix elements  1 = yes
!!* @param poldir                 : field polarization direction
subroutine negf_setpmatrix(handler, optmat_flag, poldir) 
  use uptightAPICommon  ! if:mod:use
  use negf_param         ! if:mod:use
  use uptight           ! if:mod:use
  implicit none 
  integer :: handler(DAC_handlerSize)   ! if:var:in  
  integer :: optmat_flag                ! if:var:in
  integer :: poldir                     ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT
    

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT

  pUPT%optmat = .false.
  if(optmat_flag.eq.1)  pUPT%optmat = .true.

  pUPT%poldir = poldir    

end subroutine negf_setpmatrix



!!* Initializes a given UPTIGHT instance. 
!!* Read structure file and database materials and initializes all of the
!!* UPTIGHT data containers.
!!* @param handler Number for the UPTIGHT instance to initialise and the number
!!* for the input data container, which should be used for the initialisation.
subroutine negf_inituptight(handler)
  use uptightAPICommon  ! if:mod:use
  use uptight, only : negf_init
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in

  type(UPTPointers) :: pUPTs
 
  pUPTs = transfer(handler, pUPTs)
  call negf_init(pUPTs%pUPT)

end subroutine negf_inituptight



!!* Destructs a given UPTIGHT instance.
!!* @param handler Number for the UPTIGHT instance to destroy.
subroutine negf_destructuptight(handler)
  use uptightAPICommon  ! if:mod:use
  use uptight, only : negf_destruct
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in

  type(UPTPointers) :: pUPTs

  pUPTs = transfer(handler, pUPTs)
  call negf_destruct(pUPTs%pUPT)

end subroutine negf_destructuptight



!!* Add an external potential on the atoms (external, piezo, piro, etc.) .
!!* @param handler       : handler for the UPTIGHT instance.
!!* @param nAtoms        : number of atoms in the structure
!!* @param potential     : atom-projected potential vector 
subroutine negf_addpotential(handler,nAtoms,potential)
  use uptightAPICommon !if:mod:use
  use precision, only : dp
  use negf_param        
  implicit none
  integer  :: handler(DAC_handlerSize)  ! if:var:in
  integer  :: nAtoms                    ! if:var:in
  real(dp) :: potential(nAtoms)         ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT
 
  if ( pUPT%basis%n_basis .eq. 0 ) then
     write(*,*) 'ERROR: negf_inituptight must be called first'
     return
  end if

  if ( pUPT%basis%n_basis .gt. nAtoms ) then
     write(*,*) 'ERROR: number of atoms mismatch'
     return
  end if

  if (.not.associated(pUPT%pot_data)) then
     allocate(pUPT%pot_data(nAtoms))
     pUPT%pot_data = 0.d0
  end if
  ! this will be deallocated in negf_destructuptight

  pUPT%pot_data(1:natoms) =  pUPT%pot_data(1:natoms) &
                            + potential(1:natoms)

  pUPT%potential_flag = .true. 

end subroutine negf_addpotential


subroutine negf_erasepotential(handler)
  use uptightAPICommon !if:mod:use
  use negf_param       
  implicit none
  integer  :: handler(DAC_handlerSize)  ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT
 
  if(associated(pUPT%pot_data)) then
     pUPT%pot_data=0.d0
  end if

end subroutine negf_erasepotential


!!* Computes the system hamiltonian
!!* @param handler Number for the UPTIGHT instance.
subroutine negf_createhamiltonian(handler)
  use uptightAPICommon  ! if:mod:use
  use uptight, only : negf_hamiltonian
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in

  type(UPTPointers) :: pUPTs

  pUPTs = transfer(handler, pUPTs)
  call negf_hamiltonian(pUPTs%pUPT)  


end subroutine negf_createhamiltonian

!!* Add k-points 
!!* @param handler       : handler for the UPTIGHT instance.
!!* @param k_vec(3)      : k-vector in fractional coord.  
subroutine negf_setkpoint(handler,k_vec)
  use uptightAPICommon !if:mod:use
  use precision, only : dp
  use negf_param        
  implicit none
  integer  :: handler(DAC_handlerSize)  ! if:var:in
  real(dp) :: k_vec(3)                  ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT

  ! this will be deallocated in negf_destructuptight
  pUPT%k_point(:) = k_vec(:)

end subroutine negf_setkpoint


!!* Diagonalize using lanczos iterative methods.
!!* @param handler                : Handler of the instance.
!!* @param n_vb                   : Number of valence eigenstates
!!* @param n_cb                   : Number of conduction eigenstates       
!!* @param min_iter               : minimum number of iterations (fast) (~2)
!!* @param long_iter              : number of lanczos iterations (~30) 
!!* @param max_iter               : maximum number of iterations (~3000)
!!* @param guess_vb               : valence bottom guess 
!!* @param guess_cb               : conduction bottom
!!* @param fast_tol               : fast tolerance  (1e-7) 
!!* @param long_tol               : long tolerance  (1e-10)
!!* @param ort_tol                : orthogonality tolerance (1e-6)
subroutine negf_lanczosdiag(handler, n_vb, n_cb, guess_vb, guess_cb, min_iter, &
                           long_iter, max_iter, fast_tol, long_tol, ort_tol)
  use precision, only : dp
  use uptightAPICommon  ! if:mod:use
  use negf_param
  use uptight, only : negf_lanczos
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in  
  integer :: n_vb                      ! if:var:in 
  integer :: n_cb                      ! if:var:in
  integer :: min_iter                  ! if:var:in
  integer :: long_iter                 ! if:var:in
  integer :: max_iter                  ! if:var:in
  real(dp) :: guess_vb                 ! if:var:in
  real(dp) :: guess_cb                 ! if:var:in 
  real(dp) :: fast_tol                 ! if:var:in
  real(dp) :: long_tol                 ! if:var:in
  real(dp) :: ort_tol                  ! if:var:in


  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: upt

  pUPTs = transfer(handler, pUPTs)

  upt => pUPTs%pUPT

  upt%num_vb = n_vb  
  upt%num_cb = n_cb

  upt%lambda_vb = guess_vb
  upt%lambda_cb = guess_cb

  upt%min_iter = min_iter
  upt%long_iter = long_iter
  upt%max_iter = max_iter

  upt%fast_tol = fast_tol
  upt%long_tol = long_tol
  upt%ort_tol  = ort_tol
     
  upt%fast_lanc_flag = .true.
  upt%res_flag       = .false. !not used
  upt%seed_flag      = .false.

  !..................................................................
  ! Checks consistence on the number of iterations
  !..................................................................
  if (min_iter.gt.long_iter) upt%long_iter = upt%min_iter

  if (min_iter.gt.max_iter .OR. long_iter.gt.max_iter) upt%max_iter = upt%long_iter
  
  call negf_lanczos(upt)

end subroutine negf_lanczosdiag



SUBROUTINE negf_feastsolver(handler, emin, emax, m0)
  use precision, only : dp
  use uptightAPICommon  ! if:mod:use
  use negf_param
  use uptight, only : negf_feast
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in  
  integer :: m0                      ! if:var:in 
  real(dp) :: emin                 ! if:var:in
  real(dp) :: emax                 ! if:var:in 



  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: upt

  pUPTs = transfer(handler, pUPTs)

  upt => pUPTs%pUPT

 upt%num_vb = m0
 upt%num_cb = m0
 upt%lambda_vb = emin
 upt%lambda_cb = emax

 call negf_feast(upt)

END SUBROUTINE negf_feastsolver



!!* Store eigenvectors on files
!!* @param handler Number for the UPTIGHT instance.
subroutine negf_write_states(handler)
  use uptightAPICommon  ! if:mod:use
  use uptight, only : negf_write_eigenvectors
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in

  type(UPTPointers) :: pUPTs

  pUPTs = transfer(handler, pUPTs)
  call negf_write_eigenvectors(pUPTs%pUPT)

end subroutine negf_write_states




!!* Set the number of states (used for initialize reading)
!!* @param handler Number for the UPTIGHT instance.
subroutine negf_set_num_states(handler, n_vb, n_cb)
  use uptightAPICommon                 ! if:mod:use
  implicit none 
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: n_vb                      ! if:var:in
  integer :: n_cb                      ! if:var:in
  
  type(UPTPointers) :: pUPTs
  
  pUPTs = transfer(handler, pUPTs)  
  
  pUPTs%pUPT%num_vb = n_vb  
  pUPTs%pUPT%num_cb = n_cb  
  
end subroutine negf_set_num_states




!!* Read eigenvectors from files
!!* @param handler Number for the UPTIGHT instance.
subroutine negf_read_states(handler)
  use uptightAPICommon  ! if:mod:use
  use globals, only : LST
  use savemofile, only : read_eigenstates, read_eigenvalues
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in

  type(UPTPointers) :: pUPTs
  character(LST) :: filename

  pUPTs = transfer(handler, pUPTs)


  filename = "eigv.dat"

  call read_eigenvalues(pUPTs%pUPT,filename) 

  filename = "eigvec"

  call read_eigenstates(pUPTs%pUPT,filename)


end subroutine negf_read_states


!!* Get the Hamiltonian dimension
!!* @param handler     : Number for the UPTIGHT instance.
!!* @param hdim        : Hamiltonian number of rows
subroutine negf_get_hamildim(handler, hdim)
  use uptightAPICommon  ! if:mod:use
  use sparse_matrix
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: hdim                      ! if:var:out

  type(UPTPointers) :: pUPTs

  pUPTs = transfer(handler, pUPTs)

  call get_nrow(pUPTs%pUPT%ham, hdim)

end subroutine negf_get_hamildim

!!* Get the Hamiltonian non zero values
!!* @param handler     : Number for the UPTIGHT instance.
!!* @param hdim        : Number of non zero values
subroutine negf_get_hamilnnz(handler, nnz)
  use uptightAPICommon  ! if:mod:use
  use sparse_matrix
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: nnz                       ! if:var:out

  type(UPTPointers) :: pUPTs

  pUPTs = transfer(handler, pUPTs)

  call get_nnz(pUPTs%pUPT%ham, nnz)

end subroutine negf_get_hamilnnz

!!* Get the Hamiltonian in sparse format
subroutine negf_get_csr_hamiltonian(handler, nrow, fmt, A, JA, IA)
  use precision, only : dp
  use uptightAPICommon  ! if:mod:use
  use uptight, only : negf_get_hamil
  use sparse_matrix, only : sprs_to_csr
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: nrow     ! if:var:out
  character(1) :: fmt    ! if:var:out
  complex(dp) :: A(1) ! if:var:out 
  integer :: JA(1)    ! if:var:out
  integer :: IA(1)    ! if:var:out

  type(UPTPointers) :: pUPTs
  COMPLEX ( dp ), DIMENSION( : ), POINTER :: M
  INTEGER,        DIMENSION( : ), POINTER :: Mij
  integer :: ncol

  pUPTs = transfer(handler, pUPTs)

  call negf_get_hamil(pUPTs%pUPT,nrow,ncol,fmt,M,Mij)
  
  call sprs_to_csr(M, Mij, A, JA, IA)

end subroutine negf_get_csr_hamiltonian


!!* Store eigenvectors on files
!!* @param handler    : number for the UPTIGHT instance.
!!* @param num_ev     : where to find the TB material database
!!* @param hdim       : dimensione di 
!!* @param eigenvals(num_ev)
!!* @param eigenstates(hdim,num_ev) 
subroutine negf_get_states(handler, num_ev, hdim, eigenvals, eigenstates)
  use uptightAPICommon  ! if:mod:use
  use negf_param
  use precision, only : dp
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: num_ev                        ! if:var:in
  integer :: hdim                          ! if:var:in
  real(dp) :: eigenvals(num_ev)             ! if:var:out
  complex(dp) :: eigenstates(hdim,num_ev)   ! if:var:out

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT

  eigenvals(1:num_ev) = pUPT%eigen_values(1:num_ev)

  if (hdim .ne. size(pUPT%eigen_vectors, 1)) then 
          write(*,*) 'ERROR: eigenvector dimension mismatch'
          return
  end if

  eigenstates = pUPT%eigen_vectors
  !eigenstates_im(1:hdim,1:num_ev) = imag(pUPT%eigen_vectors(1:hdim,1:num_ev))  

end subroutine negf_get_states

!!* Store eigenvectors on files
!!* @param handler    : number for the UPTIGHT instance.
!!* @param i          : left state index 
!!* @param j          : right state index
!!* @param matel      : complex matrix element
subroutine negf_get_matel(handler, i, j, matel)
  use uptightAPICommon  ! if:mod:use
  use negf_param
  use precision, only : dp
  use uptight, only : negf_get_mat_el
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: i                        ! if:var:in
  integer :: j                          ! if:var:in  
  complex(dp) :: matel                 ! if:var:out  


  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT
  
  call negf_get_mat_el(pUPT, i, j, matel)

end subroutine negf_get_matel
  

subroutine negf_project_pot(handler, i, potential, average)
  use uptightAPICommon  ! if:mod:use  use negf_param
  use precision, only : dp
  use negf_param
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: i                        ! if:var:in
  real(dp) :: potential(1)            ! if:var:in
  real(dp) :: average                 ! if:var:out

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT 
  integer :: j

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT

  if (.not.associated(pUPT%eigen_vectors)) then
     average = 1d+38
     return 
  endif
 
  do j = 1,pUPT%Ham%nrow
     ! to be completed
  end do
  
end subroutine negf_project_pot


subroutine negf_get_ion_numorbitals(handler,ion_block_vector)
  use uptightAPICommon  ! if:mod:use  use negf_param
  use precision, only : dp
  use TB_ham, only : get_ion_block_size
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: ion_block_vector(*)      ! if:var:inout
  
  type(UPTPointers) :: pUPTs
  
  pUPTs = transfer(handler, pUPTs)  
  
  call get_ion_block_size(pUPTs%pUPT, ion_block_vector)
  
end subroutine negf_get_ion_numorbitals


subroutine negf_set_verbosity(handler,verbose_lev)
  use uptightAPICommon  ! if:mod:use  use negf_param  
  implicit none
  integer :: handler(DAC_handlerSize)  ! if:var:in
  integer :: verbose_lev               ! if:var:in

  type(UPTPointers) :: pUPTs
  
  pUPTs = transfer(handler, pUPTs) 
   
  pUPTs%pUPT%verbose = verbose_lev

end subroutine negf_set_verbosity

subroutine real_test(re)
  use precision, only : dp
  real(dp) :: re     ! if:var:out  

  re = 3.14159265358979323844_dp
end subroutine real_test


subroutine complex_test(re,im,zz)
  use precision, only : dp
  implicit none
  real(dp) :: re     ! if:var:out
  real(dp) :: im     ! if:var:out
  complex(dp) :: zz  ! if:var:out

  re=4.d0; im=1.d0
  zz=re+(0.d0,1.d0)*im
 
end subroutine complex_test

!!* Add an external potential on the atoms (external, piezo, piro, etc.) .
!!* @param handler       : handler for the UPTIGHT instance.
!!* @param nAtoms        : number of atoms in the structure
!!* @param strain_xx     : atom-projected potential vector 
!!* @param strain_yy     : atom-projected potential vector 
!!* @param strain_zz     : atom-projected potential vector 
subroutine negf_setstrain(handler,nAtoms,strain_xx,strain_yy,strain_zz)
  use uptightAPICommon              !if:mod:use
  use precision, only : dp
  use negf_param        
  implicit none
  integer  :: handler(DAC_handlerSize)  ! if:var:in
  integer  :: nAtoms                    ! if:var:in
  real(dp) :: strain_xx(nAtoms)         ! if:var:in
  real(dp) :: strain_yy(nAtoms)         ! if:var:in
  real(dp) :: strain_zz(nAtoms)         ! if:var:in

  type(UPTPointers) :: pUPTs
  type(OUPT), pointer :: pUPT

  pUPTs = transfer(handler, pUPTs)
  pUPT => pUPTs%pUPT
 
  if ( pUPT%basis%n_basis .eq. 0 ) then
     write(*,*) 'ERROR: negf_inituptight must be called first'
     return
  end if

  if ( pUPT%basis%n_basis .gt. nAtoms ) then
     write(*,*) 'ERROR: number of atoms mismatch'
     return
  end if

  if (.not.associated(pUPT%basis%strain)) then
     allocate(pUPT%basis%strain(3,nAtoms))
  end if

  pUPT%basis%strain(1,1:nAtoms) = strain_xx(1:nAtoms)
  pUPT%basis%strain(2,1:nAtoms) = strain_yy(1:nAtoms)
  pUPT%basis%strain(3,1:nAtoms) = strain_zz(1:nAtoms)


  pUPT%d_onsite_shift_flag = .true.
  

end subroutine negf_setstrain
