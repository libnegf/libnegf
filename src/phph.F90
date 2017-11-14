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


module phph

  use ln_precision, only : dp
  use globals
  use ln_allocation
  use mat_def

  implicit none
  private

  public :: Tphph
  public :: init_phph, destroy_phph
  public :: load_phph_couplings
  private :: load_cubic, load_quartic

  type Tphph
    integer :: total_dim
    type(r_DNS), dimension(:,:,:), allocatable :: T3    ! T[i,j,k]
    type(r_DNS), dimension(:,:,:), allocatable :: T4  ! T[i,j,k,k]
    integer, dimension(:), allocatable :: PL_start, PL_end

    integer :: scba_iterations
    integer :: scba_iter

    logical :: include_phph
    logical :: diagonal
    logical :: cubic
    logical :: quartic
    logical :: Selfene_Gr
    logical :: Selfene_Gless
    logical :: Selfene_Hilb
    logical :: memory
    logical :: check
  end type Tphph

  integer, parameter :: GAUSSIAN = 1

contains

  subroutine init_phph(phph,total_dim,order,PL_start,PL_end)
    Type(Tphph) :: phph
    integer, intent(in) :: total_dim 
    integer, intent(in) :: order
    integer, dimension(:), intent(in) :: PL_start,PL_end

    integer :: error, i,j,k,l, nbl, size_i, size_j, size_k, size_l

    phph%total_dim = total_dim 
    phph%scba_iterations = 0 ! starts from 0  
    phph%scba_iter = 0       ! initialize at 0
    phph%include_phph =.false.
    phph%cubic = .false.
    phph%quartic = .false.
   
    
    if (order == 3 .or. order == 34) then 
       phph%include_phph =.true.
       phph%cubic = .true.
       nbl = size(PL_end,1)
       allocate(phph%T3(nbl,nbl,total_dim),stat=error)
       if (error .ne. 0) STOP 'ALLOCATION ERROR of T3'
       do i = 1, nbl
         size_i = PL_end(i)-PL_start(i)+1
         do j = max(1,i-1), i
            size_j = PL_end(j)-PL_start(j)+1
            do k = 1, total_dim
               print*,'create',i,j,k,size_i, size_j
               call create(phph%T3(i,j,k), size_i, size_j)
            end do
          end do
       end do   
    end if 
    
    if (order == 4 .or. order == 34) then 
       phph%include_phph =.true.
       phph%quartic = .true.
       allocate(phph%T4(nbl,nbl,total_dim), stat=error)
       if (error /= 0) then
         write(*,*) "ALLOCATION ERROR"; STOP 
       end if
       do i = 1, nbl
         size_i = PL_end(i)-PL_start(i)
         do j = i-1, i
           size_j = PL_end(j)-PL_start(j)
           do k = 1, total_dim 
              call create(phph%T4(i,j,k), size_i, size_j)
           enddo 
         enddo
       enddo
    end if

    phph%diagonal = .false.
    phph%Selfene_Gr = .true.
    phph%Selfene_Gless = .true.
    phph%Selfene_Hilb = .true.
    
    phph%memory = .true.
    phph%check = .false. 

    call log_allocate(phph%PL_start, size(PL_start))
    phph%PL_start = PL_start
    call log_allocate(phph%PL_end, size(PL_end))
    phph%PL_end = PL_end

  end subroutine init_phph

  subroutine destroy_phph(phph)
    Type(Tphph) :: phph

    integer :: i, j, k ,l, nbl

    nbl = size(phph%T3,1)
    do i = 1, nbl
       do j = 1, i
          do k = 1, size(phph%T3,3) 
             if (allocated(phph%T3(i,j,k)%val)) then
               call destroy(phph%T3(i,j,k))
             endif
          end do
       end do
    end do

    if (allocated(phph%T3)) deallocate(phph%T3)
    
    nbl = size(phph%T4,1)
    do i = 1, nbl
       do j = 1, i
          do k = 1, size(phph%T4,3) 
             if (allocated(phph%T4(i,j,k)%val)) then
                 call destroy(phph%T4(i,j,k))
             endif
          end do
       end do
    end do

    if (allocated(phph%T4)) deallocate(phph%T4)
    call log_deallocate(phph%PL_start)
    call log_deallocate(phph%PL_end)

  end subroutine destroy_phph

  subroutine load_phph_couplings(phph, filename)
    Type(Tphph) :: phph 
    character(*) :: filename

    if (phph%cubic) call load_cubic(phph, filename)
    if (phph%quartic) call load_quartic(phph, filename)

  end subroutine load_phph_couplings  

  subroutine load_cubic(phph, filename)
    Type(Tphph) :: phph 
    character(*) :: filename

    character(LST) :: tmpline
    character(3) :: tmp
    logical :: foundline
    integer :: ii, jj, kk, ind_i, ind_j, ind_k, bl_i, bl_j, bl_k
    real(dp), dimension(:), allocatable :: rtmp

    allocate(rtmp(phph%total_dim))

    bl_i = 1
    open(105, file = trim(filename))

    do kk = 1, size(phph%T3,3)
      read(105,*) ! read line k = ....
      
      !bl_j = 1
      !bl_k = 1
      do ii = 1, kk     
        !call find_block(ii, phph%PL_start, phph%PL_end, bl_i, ind_i)
        !do jj = 1, ii
          !call find_block(jj, phph%PL_start, phph%PL_end, bl_j, ind_j)
          !read(105,*) (rtmp(ll), ll=1, jj)
          print*,kk,ii  
          read(105,'(*(E20.10))') phph%T3(1, 1, kk)%val(ii,1:kk) 
          !phph%T3(bl_i, bl_j, kk)%val(ind_j, ind_i) = rtmp(kk) 
        !end do 
      end do
    end do 
  
    close(105)


    call Permute_Cubic(phph%T3) 

  end subroutine load_cubic


  subroutine load_quartic(phph, filename)
    Type(Tphph) :: phph 
    character(*) :: filename

  end subroutine load_quartic

  !subroutine find_block(ii, PL_start, PL_end, bb, kk)
  !  integer, intent(in) :: ii
  !  integer, intent(inout) :: bb 
  !  integer, intent(out) :: kk
  !  integer, dimension(:), intent(in) :: PL_start, PL_end
  !
  !  integer :: j
  !
  !  do j = 1, size(PL_end)
  !    if (ii > PL_start(j) .and. ii < PL_end(j)) then
  !      bb = j
  !      kk = ii - PL_start(j) + 1  
  !      exit
  !    end if  
  !  end do
  !
  !end subroutine find_block


  !  phph%T3(1,1,kk)%MAT(ii,jj)
  subroutine Permute_Cubic(IMatrix)
    !There are two possible permutations for the Third order derivatives T[i,j,k]
    !       T[i,j,k]  -->   T[k,j,i]
    !       T[i,j,k]  -->   T[j,i,k]
    !
    implicit none
    type(r_DNS), dimension(:,:,:), intent(inout):: IMatrix
    !real(dp),dimension(:,:,:),allocatable,intent(out):: FMatrix
    
    real(dp) ::  vtmp
    integer :: ii, jj, kk, n
 
    !allocate(FMatrix(3*n,3*n,3*n))
 
    n = size(IMatrix,3)

    do kk =  1,  n
        do ii = 1,  n
            do jj  = 1,  ii
    
              ! vtmp=IMatrix(ii,jj,kk)
              vtmp=IMatrix(1,1,kk)%val(ii,jj)
 
                if  (dabs(vtmp).le.1.0d-50)  then
                    !vtmp  = IMatrix(kk,jj,ii)
                    vtmp  = IMatrix(1,1,ii)%val(kk,jj)
                    if  (dabs(vtmp).le.1.0d-50)  then
                        !vtmp  = IMatrix(jj,kk,ii)
                        vtmp  = IMatrix(1,1,ii)%val(jj,kk)
                        if  (dabs(vtmp).le.1.0d-50)  then
                            print*,'NO value'
                            stop
                        else
                            continue
                        end if
                    else
                        continue
                    end if
                else
                    continue
                end if
 
              IMatrix(1,1,kk)%val(ii,jj) = vtmp
              IMatrix(1,1,kk)%val(jj,ii) = vtmp
 
            end do
        end do
    end do

  end subroutine Permute_Cubic

end module phph


