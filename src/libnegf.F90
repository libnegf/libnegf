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


module libnegf

 use ln_precision
 use ln_constants
 use ln_allocation
 use lib_param
 use mpi_globals, only : id, numprocs, id0
 use input_output
 use ln_structure
 use rcm_module
 use mat_def
 use ln_extract
 use sparsekit_drv
 use integrations

 implicit none
 private

 public :: init_negf, destroy_negf, init_ldos
 public :: set_H, set_S, set_S_id, read_HS
 public :: read_negf_in
 public :: negf_version 
 public :: destroy_matrices ! cleanup matrices in Tnegf container (H,S,rho,rhoE)
 private :: block_partition ! chop structure into PLs (CAREFUL!!!)
                            ! H need to be already ordered properly 
 public :: negf_partition_info  !write down partition info
 private :: find_cblocks        ! Find interacting contact block
 public :: set_ref_cont

 public :: compute_density_dft      ! high-level wrapping
                                    ! Extract HM and SM
                                    ! run DM calculation
 public :: compute_density_efa      ! high-level wrapping
                                    ! Extract HM and SM
                                    ! run DM calculation
 public :: compute_current          ! high-level wrapping routines
                                    ! Extract HM and SM
                                    ! run total current calculation

 public :: compute_ldos                                   
 
 
 public :: reorder, sort, swap            ! not used 
 public :: printcsr   ! debugging routines
 public :: printcsrij   ! debugging routines
 public :: getel   ! debugging routines

 integer, PARAMETER :: VBT=70
 integer, PARAMETER :: MAXNUMPLs = 10000

contains
  
  !--------------------------------------------------------------------
  ! Init libNEGF
  ! General initializations of libNEGF are currently done via files.
  ! "negf.in" contains all parameters information
  ! all needed info about the structure and matrix format.
  ! H and S are also read-in from files
  ! Some parameters are still hard-coded and need to be set from api  
  !--------------------------------------------------------------------

  subroutine init_negf(negf)
    type(Tnegf) :: negf

    call set_defaults(negf)

    !print*, '(init_negf) Initializing libnegf...'

    negf%file_struct = "negf.in"
    negf%form%formatted = .true.
    negf%isSid = .false.
    negf%form%type = "PETSc" 
    negf%form%fmt = "F" 

   end subroutine init_negf

  !--------------------------------------------------------------------
  subroutine read_HS(negf)
    type(Tnegf) :: negf
    character(11) :: fmtstring
    logical :: doesexist  
    
    open(101, file=negf%file_struct, form='formatted')  

    read(101,*) negf%file_re_H 
    read(101,*) negf%file_im_H

    read(101,*) negf%file_re_S
    read(101,*) negf%file_im_S

    if (trim(negf%file_re_S).eq.'identity') then
         negf%isSid = .true.     
    endif
    
    !print*, '(init_negf) files: ', trim(negf%file_re_H)
    !print*, '(init_negf) files: ', trim(negf%file_im_H)
    !print*, '(init_negf) files: ', trim(negf%file_re_S) 
    !print*, '(init_negf) files: ', trim(negf%file_im_S) 
 
    close(101)

    if(negf%form%formatted) then
       fmtstring = 'formatted'
    else
       fmtstring = 'unformatted'
    endif
  
    inquire(file=trim(negf%file_re_H), exist= doesexist)  
    inquire(file=trim(negf%file_im_H), exist= doesexist)  
    if (.not.doesexist) then
       write(*,*) "libNEGF error. Hamiltonian files not found"
       stop  
    endif

    open(401, file=negf%file_re_H, form=trim(fmtstring))
    open(402, file=negf%file_im_H, form=trim(fmtstring))   !open imaginary part of H

    !print*, '(init_negf) Reading H...'
    allocate(negf%H)
    call read_H(401,402,negf%H,negf%form)
    close(401)
    close(402)

    if (trim(negf%file_re_S).eq.'') negf%isSid = .true.
         
    if(.not.negf%isSid) then
       inquire(file=trim(negf%file_re_H), exist= doesexist)  
       inquire(file=trim(negf%file_im_H), exist= doesexist)  
       if (.not.doesexist) then
          write(*,*) "libNEGF error. overlap files not found"
          stop  
       endif

       open(401, file=negf%file_re_S, form=trim(fmtstring))
       open(402, file=negf%file_im_S, form=trim(fmtstring))   !open imaginary part of S

       !print*, '(init_negf) Reading S...'
       allocate(negf%S)
       call read_H(401,402, negf%S,negf%form)
       
       close(401)
       close(402)

    else
       ! create an Id matrix for S
       !print*, '(init_negf) set identity S...'
       !call create_id( negf%S,negf%H%nrow) 
       allocate(negf%S)
       call create_id( negf%S, negf%H) 

    endif
    
    negf%intHS=.true.

  end subroutine read_HS

  !--------------------------------------------------------------------
  subroutine set_H(negf, nrow, nzval, colind, rowpnt)
    type(Tnegf) :: negf
    integer :: nrow
    complex(dp) :: nzval(*)
    integer :: colind(*)
    integer :: rowpnt(*)

    integer :: nnz, i, base

    base = 0
    if (rowpnt(1) == 0) base = 1

    nnz = rowpnt(nrow+1)-rowpnt(1)

    allocate(negf%H)
    call create(negf%H,nrow,nrow,nnz)

    do i = 1, nnz
      negf%H%nzval(i) = nzval(i)
      negf%H%colind(i) = colind(i) + base
    enddo
    do i = 1,nrow+1
      negf%H%rowpnt(i) = rowpnt(i) + base 
    enddo  
    negf%intHS=.true.

    !print*,'nrow = ',negf%H%nrow
    !print*,'nnz = ',negf%H%rowpnt(nrow+1)-1, negf%H%nnz

  end subroutine set_H

  !--------------------------------------------------------------------
  subroutine set_S(negf, nrow, nzval, colind, rowpnt)
    type(Tnegf) :: negf
    integer :: nrow
    complex(dp) :: nzval(*)
    integer :: colind(*)
    integer :: rowpnt(*)

    integer :: nnz, i, base

    base = 0
    if (rowpnt(1) == 0) base = 1

    nnz = rowpnt(nrow+1)-rowpnt(1)

    allocate(negf%S)
    call create(negf%S,nrow,nrow,nnz)

    do i = 1, nnz
      negf%S%nzval(i) = nzval(i)
      negf%S%colind(i) = colind(i) + base
    enddo
    do i = 1,nrow+1
      negf%S%rowpnt(i) = rowpnt(i) + base
    enddo  
    negf%intHS=.true.
    
    !print*,'nrow = ',negf%H%nrow
    !print*,'nnz = ',negf%H%rowpnt(nrow+1)-1, negf%H%nnz


  end subroutine set_S
  !--------------------------------------------------------------------
  subroutine set_S_id(negf, nrow)
    type(Tnegf) :: negf
    integer :: nrow
    !integer :: nnz

    allocate(negf%S)
    call create_id(negf%S, nrow) 

  end subroutine set_S_id
  
  !--------------------------------------------------------------------
  subroutine read_negf_in(negf)
    type(Tnegf) :: negf
    Integer :: ncont, nbl
    Integer, dimension(:), allocatable :: PL_end, cont_end, surf_end, cblk
    character(32) :: tmp

    open(101, file=trim(negf%file_struct), form='formatted')  
  
    read(101,*) negf%file_re_H 
    read(101,*) negf%file_im_H

    read(101,*) negf%file_re_S
    read(101,*) negf%file_im_S

    if (trim(negf%file_re_S).eq.'identity') then
         negf%isSid = .true.     
    endif

    read(101,*) tmp, ncont

    call log_allocate(cblk,ncont)
    call log_allocate(cont_end,ncont)
    call log_allocate(surf_end,ncont)

    read(101,*) tmp, nbl

    if (nbl .gt. 0) then
       call log_allocate(PL_end,nbl)
   print*,'read blocks ',nbl 
       read(101,*) PL_end(1:nbl)
    end if

    read(101,*) tmp,  cont_end(1:ncont)
    read(101,*) tmp,  surf_end(1:ncont)

    if (nbl .eq. 0) then
       call log_allocate(PL_end, MAXNUMPLs)  
       call block_partition(negf%H, surf_end(1), cont_end, surf_end, ncont, nbl, PL_end)   
    endif
    
    call find_cblocks(negf%H ,ncont, nbl, PL_end, cont_end, surf_end, cblk)

    !if (negf%verbose .gt. 50) then
    !   write(*,*) '(init NEGF) nbl: ', nbl
    !   write(*,*) "(init NEGF) blocks interactions:",cblk(1:ncont)     
    !endif

    call init_structure(negf, ncont, nbl, PL_end, cont_end, surf_end, cblk)
       
    call log_deallocate(PL_end)
    call log_deallocate(cblk)
    call log_deallocate(cont_end)
    call log_deallocate(surf_end)

    read(101,*) tmp,  negf%mu_n, negf%mu_p
    read(101,*) tmp,  negf%Ec, negf%Ev
    read(101,*) tmp,  negf%DeltaEc, negf%DeltaEv
    read(101,*) tmp,  negf%Emin, negf%Emax, negf%Estep
    read(101,*) tmp,  negf%kbT
    read(101,*) tmp,  negf%wght
    read(101,*) tmp,  negf%Np_n(1:2)
    read(101,*) tmp,  negf%Np_p(1:2)
    read(101,*) tmp,  negf%Np_real(1)
    read(101,*) tmp,  negf%n_kt
    read(101,*) tmp,  negf%n_poles
    read(101,*) tmp,  negf%g_spin
    read(101,*) tmp,  negf%delta
    read(101,*) tmp,  negf%nLDOS
    !call log_allocatep(negf%LDOS,2,negf%nLDOS)
    !read(101,*) tmp,  negf%LDOS
    read(101,*) tmp,  negf%Efermi(1:ncont)  ! Will be 0 from TC
    read(101,*) tmp,  negf%mu(1:ncont)      ! Will be the Electrochemical potential

    close(101)

    !print*, '(init NEGF) done'

  end subroutine read_negf_in
!--------------------------------------------------------------------

  subroutine negf_version(negf)
    type(Tnegf) :: negf
    !character(3), parameter :: SVNVER= __SVNREVISION 
    !character(3),parameter :: MODIF= __MODIFIED 
    character(3), parameter :: GITVER= __GITREVISION 
    character(10),parameter :: DATE= __COMPDATE 
 
    write(*,'(a21,a20,2x,a10)') '(libnegf) version: 1.',TRIM(GITVER), & 
                                         TRIM(DATE) 

  end subroutine negf_version

!--------------------------------------------------------------------
   subroutine negf_partition_info(negf)
      type(Tnegf) :: negf
       
      integer :: i

      write(*,*) "(LibNEGF) Partitioning:"
      write(*,*) "Number of blocks: ",negf%str%num_Pls
      !write(*,*) negf%str%mat_PL_end(:)
      write(*,*) "Contact interactions:",negf%str%cblk(:)     

      open(1001,file='blocks.dat')
        write(1001,*) 1
        do i = 1, negf%str%num_Pls       
           write(1001,*)  negf%str%mat_PL_end(i)
        enddo
      close(1001)

 end subroutine negf_partition_info

!--------------------------------------------------------------------
  subroutine init_structure(negf,ncont,nbl,PL_end,cont_end,surf_end,cblk)
    type(Tnegf) :: negf
    Integer :: ncont, nbl
    Integer, dimension(:) :: PL_end, cont_end, surf_end, cblk

    call create_Tstruct(ncont, nbl, PL_end, cont_end, surf_end, cblk, negf%str)

  end subroutine init_structure

!--------------------------------------------------------------------
  subroutine destroy_negf(negf)
    type(Tnegf) :: negf   

    call destroy_matrices(negf)

    call kill_Tstruct(negf%str) 

    if (allocated(negf%LDOS)) call destroy_ldos(negf%LDOS)
    if (associated(negf%tunn_mat)) call log_deallocatep(negf%tunn_mat)
    if (associated(negf%ldos_mat)) call log_deallocatep(negf%ldos_mat)    
    if (associated(negf%currents)) call log_deallocatep(negf%currents)    

    !call destroy_emesh(negf)

  end subroutine destroy_negf
  
  !--------------------------------------------------------------------
  subroutine init_ldos(negf,nregs,sizes)
    type(Tnegf) :: negf
    integer :: nregs
    integer, dimension(:) :: sizes

    integer :: err, i
    
    allocate(negf%ldos(nregs),stat=err)
    if (err/=0) stop 'allocation error of ldos'

    do i=1, nregs
      call log_allocate(negf%ldos(i)%indexes,sizes(i))
    end do
    
  end subroutine init_ldos
  !--------------------------------------------------------------------
  subroutine destroy_ldos(ldos)
    type(intarray), dimension(:), allocatable :: ldos 
      
    integer :: err, i

    do i=1, size(ldos)
      call log_deallocate(ldos(i)%indexes)
    end do

    deallocate(ldos)
    
  end subroutine destroy_ldos
  
  !-------------------------------------------------------------------- 
  subroutine create_DM(negf)
    type(Tnegf) :: negf   
  
    if (negf%intDM) then
       if (.not.associated(negf%rho)) allocate(negf%rho)
       if (.not.associated(negf%rho_eps)) allocate(negf%rho_eps)
    endif  

  end subroutine create_DM


  subroutine destroy_matrices(negf)
    type(Tnegf) :: negf   
    integer :: i

    do i=1,negf%str%num_conts
       if (allocated(negf%HC(i)%val)) call destroy(negf%HC(i))
       if (allocated(negf%SC(i)%val)) call destroy(negf%SC(i))
       if (allocated(negf%HMC(i)%val)) call destroy(negf%HMC(i))
       if (allocated(negf%SMC(i)%val)) call destroy(negf%SMC(i))
    enddo

    call destroy_HS(negf)
    call destroy_DM(negf)

  end subroutine destroy_matrices
 
!--------------------------------------------------------------------
  subroutine destroy_HS(negf)
    type(Tnegf) :: negf

    if (negf%intHS) then
      if (associated(negf%H)) then
        if (allocated(negf%H%nzval)) then
           !print*,'(destroy) deallocate negf%H',%LOC(negf%H%nzval)
           call destroy(negf%H)
        end if
        deallocate(negf%H)
        nullify(negf%H)
      endif

      if (associated(negf%S)) then
        if (allocated(negf%S%nzval)) then
           !print*,'(destroy) deallocate negf%S',%LOC(negf%S%nzval)
           call destroy(negf%S) 
        end if
        deallocate(negf%S)
        nullify(negf%S) 
      endif

    endif

  end subroutine destroy_HS   

!--------------------------------------------------------------------
  subroutine destroy_DM(negf)
    type(Tnegf) :: negf   

    if (negf%intDM) then
      
      if (associated(negf%rho)) then
        if (allocated(negf%rho%nzval)) then
           !print*,'(destroy) deallocate negf%rho',%LOC(negf%rho%nzval)
           call destroy(negf%rho) 
        end if
        deallocate(negf%rho)
        nullify(negf%rho)
      endif  
      
      if (associated(negf%rho_eps)) then
        if (allocated(negf%rho_eps%nzval)) then
           !print*,'(destroy) deallocate negf%rho_eps',%LOC(negf%rho_eps%nzval)
           call destroy(negf%rho_eps) 
        end if
        deallocate(negf%rho_eps)
        nullify(negf%rho_eps) 
      endif

    endif  

  end subroutine destroy_DM  
  
  !-------------------------------------------------------------------------------
  ! Compact collection of calls to extract device/contact H and S 
  ! and compute density matrix using contour + real axis integration
  ! Should be used for dftt calculations
  !  
  ! NOTE: the returned DensMat and EnMat are masked with S
  ! Matrix Structure:  CSR 
  !                    %nrow=%ncol=(Full squared Hamiltonian size)
  !                    %nnz = Only non-zero elements of the blocks
  !
  !                    +-----+--+--+--+
  !                    !  D  !C1!C2!C3!  masked with the S matrix 
  !                    !     !  !  !  !
  !                    +-----+--+--+--+
  !                    ! C1  !0 !0 !0 !  The lower part of DensMat
  !                    +-----+--+--+--+  is filled with 0.
  !                    ! C2  !0 !0 !0 !
  !                    +-----+--+--+--+  negf%outer=0,1,2 is used
  !                    ! C3  !0 !0 !0 !  in order to compute Ci
  !                    +-----+--+--+--+
  !-------------------------------------------------------------------------------
  subroutine compute_density_dft(negf)
    type(Tnegf) :: negf


    call extract_device(negf)

    call extract_cont(negf)

    ! Reference contact for contour/real axis separation
    call set_ref_cont(negf)
    
    !Decide what to do with surface GFs.
    !sets readOldSGF: if it is 0 or 1 it is left so 
    if (negf%readOldSGF.eq.2) then
      if(negf%iteration.eq.1) then        
        negf%readOldSGF=2  ! compute and save SGF on files
      else
        negf%readOldSGF=0  ! read from files
      endif
    endif

    if (negf%Np_n(1)+negf%Np_n(2)+negf%n_poles.gt.0) then
      call contour_int_def(negf)
      call contour_int(negf)
    endif

    if (negf%Np_real(1).gt.0) then
      call real_axis_int_def(negf)
      call real_axis_int(negf)
    endif

    call destroy_matrices(negf)

  end subroutine compute_density_dft


  !-------------------------------------------------------------------------------
  ! Compact collection of calls to extract device/contact H and S 
  ! and compute density matrix
  !
  ! It has been used to interface libnegf to TiberCAD
  ! Computes density for CB semiconductor 
  !-------------------------------------------------------------------------------
  subroutine compute_density_efa(negf, q)

    type(Tnegf) :: negf
    real(dp), dimension(:) :: q
    complex(dp), dimension(:), allocatable :: q_tmp

    integer :: k

    call extract_device(negf)

    call extract_cont(negf)

    call set_ref_cont(negf)

    call create_DM(negf)

    if (negf%Np_n(1)+negf%Np_n(2)+negf%n_poles.gt.0) then
       call contour_int_n_def(negf)
       call contour_int(negf)
    else 
       ! HACKING: THIS WAY COMPUTES DM FOR ALL CONTACTS
       negf%refcont = negf%str%num_conts+1  
    endif

    if (negf%Np_real(1).gt.0) then
       call real_axis_int_def(negf)
       call real_axis_int(negf)
    endif

    ! We need not to include S !!!!
    if (negf%rho%nrow.gt.0) then
       call log_allocate(q_tmp, negf%rho%nrow)

       call getdiag(negf%rho, q_tmp)

       do k = 1, size(q)
          q(k) = real(q_tmp(k))
       enddo

       call log_deallocate(q_tmp)
    else
       q = 0.d0
    endif

    call destroy_matrices(negf)

  end subroutine compute_density_efa

  !-------------------------------------------------------------------------------
  subroutine compute_ldos(negf)
    type(Tnegf) :: negf
    
    call extract_device(negf)
    call extract_cont(negf)
    call tunneling_int_def(negf)
    call ldos_int(negf)
    call destroy_matrices(negf)

  end subroutine compute_ldos  
  
  !-------------------------------------------------------------------------------
  subroutine compute_current(negf)

    type(Tnegf) :: negf
    
    integer :: flagbkup

    call extract_device(negf)
    
    call extract_cont(negf)
    
    flagbkup = negf%readOldSGF
    if (negf%readOldSGF.ne.1) then
       negf%readOldSGF = 1
    end if

    call tunneling_int_def(negf)

    call tunneling_and_current(negf)
   
    !!GP Locally writing energy dependent data is not meaningful in the MPI
    !implementation, because the gathering is done externally.
    ! An implementation node by node is still active, for debugging purposes 
    !call write_tunneling_and_dos(negf)
    
    call destroy_matrices(negf)
 
    negf%readOldSGF = flagbkup
  
  end subroutine compute_current


  ! --------------------------------------------------------------------------------
  ! GP Left in MPI version for debug purpose only. This will write a separate
  ! file for every ID, which is not possible on all architectures 
  subroutine write_current(negf)

    type(Tnegf), pointer :: negf

    integer :: i1
    logical :: lex
    character(6) :: idstr

    write(idstr,'(i6.6)') id
    inquire(file=trim(negf%out_path)//'current_'//idstr//'.dat',EXIST=lex)
    
    if (lex) then
       open(101,file=trim(negf%out_path)//'current_'//idstr//'.dat',position='APPEND')
    else
       open(101,file=trim(negf%out_path)//'current_'//idstr//'.dat')
    endif

    do i1=1,size(negf%currents)

       write(101,'(1x,a,i3,i3,a,i3,a,ES14.5,a,ES14.5,a)') 'contacts:',negf%ni(i1),negf%nf(i1), &
            '; k-point:',negf%kpoint,'; current:', negf%currents(i1),' A'

    end do

    close(101)

  end subroutine write_current  
  !-------------------------------------------------------------------------------
  
  !---- SAVE TUNNELING AND DOS ON FILES -----------------------------------------------
  ! GP Left in MPI version for debug purpose only. This will write a separate
  ! file for every ID, which is not possible on all architectures 
  subroutine write_tunneling_and_dos(negf)

    type(Tnegf) :: negf

    integer :: Nstep, i, i1, iLDOS, size_ni
    character(6) :: ofKP, idstr
    real(dp) :: E

    if (associated(negf%tunn_mat) .and. negf%writeTunn) then
        
        Nstep = size(negf%tunn_mat,1) 
        size_ni = size(negf%tunn_mat,2)
        
        write(ofKP,'(i6.6)') negf%kpoint
        write(idstr,'(i6.6)') id
        
        open(1021,file=trim(negf%out_path)//'tunneling_'//ofKP//'_'//idstr//'.dat')
        
        !print*,'ENE CONV=',negf%eneconv
        negf%eneconv=1.d0
        
        do i = 1,Nstep
       
          E=(negf%Emin+negf%Estep*(i-1))
          
          WRITE(1021,'(E17.8,20(E17.8))') E*negf%eneconv, &
              (negf%tunn_mat(i,i1), i1=1,size_ni)
          
        enddo
        
        close(1021)
        
    endif

    if (associated(negf%ldos_mat) .and. negf%writeLDOS .and. negf%nLDOS.gt.0) then
        
        Nstep = size(negf%ldos_mat,1)
        
        write(ofKP,'(i6.6)') negf%kpoint
        write(idstr,'(i6.6)') id
        
        open(1021,file=trim(negf%out_path)//'LEDOS_'//ofKP//'_'//idstr//'.dat')
        
        do i = 1,Nstep
          
          E=(negf%Emin+negf%Estep*(i-1))
          
          WRITE(1021,'(E17.8,10(E17.8))') E*negf%eneconv, & 
              ((negf%ldos_mat(i,iLDOS)/negf%eneconv), iLDOS=1,negf%nLDOS)        
          
        end do
        
        close(1021)
       
    endif
    
  end subroutine write_tunneling_and_dos
  !---------------------------------------------------------------------------
  ! Sets the Reference contact for non-eq calculations
  ! 
  ! The behaviour depends on how negf%minmax has been set.
  !
  ! minmax = 0 : refcont is chosen at the minimum   mu
  ! minmax = 1 : refcont is chosen at the maximum   mu 
  ! 
  subroutine set_ref_cont(negf)

    type(TNegf) :: negf

    integer :: nc_vec(1), ncont

    ncont = negf%str%num_conts

    if (ncont > 0) then
      if (negf%minmax .eq. 0) then
         negf%muref = minval(negf%mu(1:ncont))
         nc_vec = minloc(negf%mu(1:ncont))  
      else
         negf%muref = maxval(negf%mu(1:ncont))
         nc_vec = maxloc(negf%mu(1:ncont))
      endif
      negf%refcont = nc_vec(1)
    else
      negf%muref = negf%mu(1)
      negf%refcont = 1  
    endif  
     
  end subroutine set_ref_cont

  !////////////////////////////////////////////////////////////////////////
  ! RCM algorithm for reordering.
  ! 
  ! Actually not used because it is not suitable for contacted structures
  !////////////////////////////////////////////////////////////////////////
  subroutine reorder(mat)
    type(z_CSR) :: mat
    

    type(z_CSR) :: P, Tmp

    integer, dimension(:), allocatable :: perm
    integer :: i, nrow

    nrow=mat%nrow

    call log_allocate(perm,nrow)

    call genrcm(nrow, mat%nnz, mat%rowpnt, mat%colind, perm)

    call create(P,nrow,nrow,nrow)

    do i=1,nrow
       P%nzval(i)=1
       P%colind(i)=perm(i)
       P%rowpnt(i)=i
    enddo
    P%rowpnt(nrow+1)=nrow+1

    
    call create(Tmp,nrow,nrow,mat%nnz)

    call zamub_st(P,mat,Tmp)

    call ztransp_st(P)

    call zamub_st(Tmp,P,mat)   

    call destroy(P,Tmp)

    call log_deallocate(perm)
 
  end subroutine reorder

  !----------------------------------------------------------------------
  ! Authomatic Block partitioning. The Hamiltonian must be already sorted
  !----------------------------------------------------------------------
  subroutine block_partition(mat,nrow,cont_end,surf_end,ncont,nbl,blks)
    type(z_CSR), intent(in) :: mat
    integer, intent(in) :: nrow
    integer, dimension(:), intent(in) :: cont_end
    integer, dimension(:), intent(in) :: surf_end
   integer, intent(in) :: ncont
    integer, intent(out) :: nbl
    integer, dimension(:), intent(inout) :: blks 

    integer :: j, k, i
    integer :: i1, i2

    integer :: rn, rnold, tmax, rmax, maxmax
    integer :: dbuff, minsize, minv, maxv

    !nrow = mat%nrow
    
     minsize = 0
     do i1 = 1, ncont
        maxv = 0
        minv = 400000000
        do k = surf_end(i1)+1, cont_end(i1)
           do i = mat%rowpnt(k), mat%rowpnt(k+1)-1
            if (mat%colind(i).le.nrow .and.  mat%colind(i).lt.minv) minv = mat%colind(i)
            if (mat%colind(i).le.nrow .and.  mat%colind(i).gt.maxv) maxv = mat%colind(i)
           end do
        end do
        if (maxv-minv+1 .gt. minsize) minsize = maxv - minv + 1
    end do
                                                                                
    ! Find maximal stancil of the matrix and on which row
    !  ( Xx     )
    !  ( xXxx   )  
    !  (  xXxxx )  <- maxmax = 3 ; rmax = 3
    !  (   xXx  )
    maxmax = 0
    do j=1,nrow

       i1 = mat%rowpnt(j)
       i2 = mat%rowpnt(j+1) - 1

       tmax = 0
       do i = i1, i2
           if ( mat%colind(i).le.nrow .and. (mat%colind(i)-j) .gt. tmax) then
                tmax = mat%colind(i)-j
           endif
       enddo  
 
       if(tmax .gt. maxmax) then 
          maxmax = tmax
          rmax = j
       endif

       dbuff = maxmax        ! dbuff should be linked to maxmax
       minsize = max((dbuff+1)/2,minsize)  
    enddo

    !write(*,*) 'minsize=',minsize
    !write(*,*) 'maxrow=',rmax

    ! Define central block 
    rn = rmax - maxmax/2 - dbuff 

    !write(*,*) 'rn=',rn


    if(rn-dbuff.ge.0) then 

       blks(1) = rn-1  ! fine del blocco precedente

       nbl = 1

       do 
          
          do j = rn, minsize, -1

             rnold = rn
             i1 = mat%rowpnt(j-minsize+1)
             i2 = mat%rowpnt(j+1) - 1
 
             !k = maxval(mat%colind(i1:i2))
             k = 0
             do i = i1, i2
                if ( mat%colind(i).le.nrow .and. (mat%colind(i)) .gt. k) then
                   k = mat%colind(i)
                endif
             enddo
       
             if(k.lt.rn) then
                rn = j       
                nbl = nbl + 1
                if (nbl.gt.MAXNUMPLs) call errormsg()
                blks(nbl) = j-1 ! fine del blocco precedente
                exit
             endif
          enddo

          if(rn.le.minsize .or. rnold.eq.rn) then
             exit
          endif
          
       enddo

       rn = rmax - maxmax/2 - dbuff

    else
       nbl= 0
       rn = 1

    endif

    do 

       do j = rn, nrow-minsize+1, 1

          rnold = rn
          i1 = mat%rowpnt(j)
          i2 = mat%rowpnt(j+minsize) - 1

          !k = minval(mat%colind(i1:i2))
          k = nrow 
          do i = i1, i2
             if ( mat%colind(i).le.nrow .and. (mat%colind(i)) .lt. k) then
                k = mat%colind(i)
             endif
          enddo
         

          if(k.gt.rn) then  
             rn = j
             nbl = nbl + 1 
             if (nbl.gt.MAXNUMPLs) call errormsg()
             blks(nbl) = j-1 ! fine del blocco 
             exit
          endif
       enddo

       if(nrow-rn.le.minsize .or. rnold.eq.rn) then
          exit
       endif
    enddo

    nbl = nbl + 1
    if (nbl.gt.MAXNUMPLs) call errormsg()
    blks(nbl) = nrow

    ! Sorting blocks
    
    do i = 1, nbl
       do j = i+1, nbl 
          if(blks(j).lt.blks(i)) then
             k = blks(i)
             blks(i) = blks(j)
             blks(j) = k
          endif
       enddo
    enddo


  end subroutine block_partition

!----------------------------------------------------------------------------
  subroutine errormsg()

     write(*,*) "ERROR: Maximum number of PLs exceeded"
     write(*,*) "increase the value of MAXNUMPLS in libnegf.F90"
     write(*,*) "and recompile the library"

     STOP !should rise an exception

  end subroutine errormsg
 
!----------------------------------------------------------------------------
  subroutine find_cblocks(mat ,ncont, nbl, PL_end, cont_end, surf_end, cblk)
    type(z_CSR), intent(in) :: mat
    integer, intent(in) :: ncont
    integer, intent(in) :: nbl
    integer, dimension(:), intent(in) :: PL_end 
    integer, dimension(:), intent(in) :: cont_end
    integer, dimension(:), intent(in) :: surf_end
    integer, dimension(:), intent(inout) :: cblk

    integer :: j1,k,i,min,max
    integer, dimension(:), allocatable :: PL_start

    call log_allocate(PL_start,nbl)

    PL_start(1) = 1

    do i = 2, nbl
       PL_start(i) = PL_end(i-1) + 1
    enddo


    do j1 = 1, ncont

       max = 0
       min = 400000000

       do k = surf_end(j1)+1, cont_end(j1)

          do i = mat%rowpnt(k), mat%rowpnt(k+1)-1

             if (mat%colind(i).le.PL_end(nbl) .and.  mat%colind(i).lt.min) min = mat%colind(i)
             if (mat%colind(i).le.PL_end(nbl) .and.  mat%colind(i).gt.max) max = mat%colind(i)

          end do

       end do

       do k = 1, nbl

          if( max .le. PL_end(k) ) then
             cblk(j1) = k

             if( min .ge. PL_start(k) ) then                
                exit
             else
                write(*,*) "(LibNEGF) Partitioning:"
                write(*,*) "Number of blocks: ",nbl
                write(*,*) "PL_end: ",PL_end(1:nbl)
                write(*,*) "Contact interaction: ",cblk(j1)      
                write(*,'(a,i3,a)') " ERROR: contact",j1," interacting with more than one block"
                write(*,*) "min ",min,"max ",max
                stop
             end if

          end if

       end do

    end do

    call log_deallocate(PL_start)

  end subroutine find_cblocks

!----------------------------------------------------------------------------

  Subroutine sort(blks, Ipt)
    ! *
    ! ***********************************
    ! * Sort Array X(:) in ascendent order.
    ! * If present Ipt, a pointer with the 
    ! * changes is returned in Ipt. 
    ! ***********************************
 
    Integer, Intent (inout) :: blks(:)
    Integer, Intent (out), Optional :: Ipt(:)
 
    Integer :: Rtmp
    Integer :: i, j
 
    If (Present(Ipt)) Then
       Forall (i=1:Size(blks)) Ipt(i) = i
 
       Do i = 2, Size(blks)
          Rtmp = blks(i)
          Do j = i-1, 1, -1
             If (Rtmp < blks(j)) Then
                blks(j+1) = blks(j)
                call Swap(blks, j, j+1)
             Else
                Exit
             End If
          End Do
          blks(j+1) = Rtmp
       End Do
    Else
       Do i = 2, Size(blks)
          Rtmp = blks(i)
          Do j = i-1, 1, -1
             If (Rtmp < blks(j)) Then
                blks(j+1) = blks(j)
             Else
                Exit
             End If
          End Do
          blks(j+1) = Rtmp
       End Do
    End If
 
    Return
  End Subroutine sort
 
  ! ***********************************
  ! * Swaps elements I and J of array X(:). 
  ! ***********************************
  Subroutine Swap(X, i, j)
 
    Integer, Intent (inout) :: X(:)
    Integer, Intent (in) :: i, j
 
    Integer :: Itmp
 
    Itmp = X(I)
    X(I) = X(J)
    X(J) = Itmp
 
    Return
  End Subroutine Swap

  !------------------------------------------
  subroutine printcsr(id, mat)
    integer :: id
    type(z_csr) :: mat
    
    call zprint_csrcoo(id, mat, 'c')
  end subroutine printcsr

  subroutine printcsrij(id, mat, i, j)
    integer :: id
    type(z_csr) :: mat
    integer :: i, j

    write(id,*) i,j,getelement(i,j,mat)

  end subroutine printcsrij

  function getel(mat,i,j) result(el)
    type(z_csr) :: mat
    integer :: i, j
    complex(dp) :: el
    
    el = getelement(i,j,mat)

  end function getel

end module libnegf
