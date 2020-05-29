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


 module negf_energy_mesh

  use negf_ln_precision
  implicit none
  private

  public :: elem, mesh
  public :: create_mesh
  public :: refine
  public :: destroy_mesh
  public :: test_mesh 
  public :: get_error
  
  type elem
    integer :: lev
    logical :: active
    integer :: ind
    real(dp) :: pnt(3)
    integer :: map(3)
    real(dp) :: error
    type(elem), pointer :: parent =>null()
    type(elem), pointer :: child1 =>null()
    type(elem), pointer :: child2 =>null()
  end type elem

  type TelemPointer
    type(elem), pointer :: pelem => null()
  end type TelemPointer   

  type mesh     
     type(elem), dimension(:), pointer :: el0 
     type(TelemPointer), dimension(:), allocatable :: pactive 
     integer :: maxreflev
     integer :: maxind 
     integer :: maxpnt
     integer :: iactive 
  end type mesh        

  integer :: meshmem

contains
 
  subroutine create_mesh(emesh,reflevel,x1,x2,n,ioff)
    type(mesh) :: emesh
    real(dp) :: x1,x2
    integer :: reflevel
    integer :: n
    integer, optional :: ioff


    real(dp) :: x(n), w(n)
    integer :: nelem, ierr, k, i
    type(elem), pointer :: el
    integer :: ioffset

    if (.not.present(ioff)) then
            ioffset = 0
    else 
            ioffset = ioff
    endif
    if (mod(n,2).eq.0) n=n+1 

    call trapez(x1,x2,x,w,n)
      
    nelem = (n-1)/2
    emesh%maxreflev = reflevel

    allocate(emesh%el0(nelem), stat=ierr)
    if(ierr.ne.0) stop 'ERROR: Cannot allocate mesh'

    allocate(emesh%pactive(2**(emesh%maxreflev-1)*nelem), stat=ierr)
    if(ierr.ne.0) stop 'ERROR: Cannot allocate mesh'

    do i=1,nelem
      emesh%pactive(i)%pelem => emesh%el0(i)
    enddo

    meshmem = 0
    k=1
    do i=2,n-1,2
      el => emesh%el0(k)
      el%pnt = (/x(i-1),x(i),x(i+1)/)
      el%map = (/ioffset+i-1,ioffset+i,ioffset+i+1/)
      el%lev = 1 
      el%ind = k
      el%active = .true.
      k = k + 1
      meshmem = meshmem + 10 
    enddo

    emesh%maxind = k-1
    emesh%maxpnt = n 

  end subroutine create_mesh 

  !--------------------------------------------------------------------  
  subroutine destroy_mesh(emesh)
    type(mesh) :: emesh
    
    integer :: nelem,i,k
    type(elem), pointer :: el

    nelem = size(emesh%el0)

    do i=1,nelem
      
      el=>emesh%el0(i)
      call traverse(el)

    enddo

    deallocate(emesh%el0)
    deallocate(emesh%pactive)

    meshmem = meshmem - 10 * nelem

  end subroutine destroy_mesh
  !--------------------------------------------------------------------  
  
  recursive subroutine traverse(el)
     type(elem), pointer :: el
   
     if (el%active) then 
          call destroy_el(el)
     else        
          call traverse(el%child1)
          call traverse(el%child2)
     endif

  end subroutine traverse 
  !--------------------------------------------------------------------  

  subroutine refine(emesh,el)
    type(mesh) :: emesh
    type(elem), pointer, intent(in) :: el

    integer :: maxind,ind
    integer :: i

    ind = el%ind
    !print*,'create children of',ind
 
    call create_children(el,emesh%maxpnt)

    !print*,'add new element in the list '

    maxind = emesh%maxind

    !print*,'update active elements: now', maxind,size(emesh%pactive)

    do i = maxind+1,ind+2,-1
       emesh%pactive(i)%pelem => emesh%pactive(i-1)%pelem
       emesh%pactive(i)%pelem%ind = i  
    enddo

    el%child2%ind = ind+1
    el%child2%active = .true.
    emesh%pactive(ind+1)%pelem => el%child2


    el%child1%ind = ind
    el%child1%active = .true.
    emesh%pactive(ind)%pelem => el%child1 

    el%active = .false.
    emesh%maxind = maxind + 1 
    emesh%maxpnt = emesh%maxpnt + 2 

  endsubroutine refine
  !--------------------------------------------------------------------  
  subroutine create_children(el,npoints)
    type(elem), pointer :: el
    integer, intent(in) :: npoints

    type(elem), pointer :: el1,el2

    call create_el(el1)
    call create_el(el2)

    el1%pnt(1) = el%pnt(1)
    el1%pnt(3) = el%pnt(2)
    el1%pnt(2) = (el%pnt(1) + el%pnt(2))/2.0_dp    

    el2%pnt(1) = el%pnt(2)
    el2%pnt(3) = el%pnt(3)
    el2%pnt(2) = (el%pnt(2) + el%pnt(3))/2.0_dp    

    el1%map(1) = el%map(1)
    el1%map(2) = npoints + 1    
    el1%map(3) = el%map(2)

    el2%map(1) = el%map(2)
    el2%map(2) = npoints + 2    
    el2%map(3) = el%map(3)

    el1%lev = el%lev + 1
    el2%lev = el%lev + 1
 
    el1%parent => el
    el2%parent => el

    el%child1 => el1
    el%child2 => el2

  end subroutine create_children
  !--------------------------------------------------------------------  
  
  subroutine create_el(el)
    type(elem), pointer :: el

    integer :: ierr

    allocate(el, stat=ierr)
    if(ierr.ne.0) stop 'ERROR: Cannot allocate el1'

    meshmem = meshmem + 10

  end subroutine create_el

  !--------------------------------------------------------------------  
  subroutine destroy_el(el)
    type(elem), pointer :: el

    deallocate(el)

    meshmem = meshmem - 10

  end subroutine destroy_el

  !--------------------------------------------------------------------  
  subroutine trapez(x1,x2,x,w,n)

    real(kind=dp), PARAMETER :: ACC = 1d-15

    INTEGER n
    real(kind=dp) :: x1,x2,x(n),w(n)

    real(dp) :: d
    integer :: i

    d = (x2-x1)/(n-1)

    w = d * 1.0_dp
    w(1) = d * 0.5_dp
    w(n) = d * 0.5_dp

    do i = 1, n
       x(i) = ( x1*(n-i) + x2*(i-1) ) / (n-1) 
    enddo
 
  end subroutine trapez
  !--------------------------------------------------------------------  

  subroutine test_mesh()
     type(mesh) :: emesh

     integer :: i,npoints
     type(elem), pointer :: el

     call create_mesh(emesh,5,-2.77_dp,0.0_dp,51)

     print*,'meshmem=', meshmem

     npoints = 51

     do i =1, size(emesh%el0)

       el => emesh%el0(i)
       call create_children(el,npoints)
       
     enddo 

     print*,'meshmem=', meshmem

     call destroy_mesh(emesh)
     
     print*,'meshmem=', meshmem

  end subroutine test_mesh

  function get_error(emesh,flag) result(error)
     type(mesh) :: emesh
     integer :: flag ! 0 = average
                     ! 1 = maximum

     integer :: i
     real(dp) :: error 

     error = 0.0_dp
 
     if (flag.eq.1) then

       do i = 1, emesh%maxind
         if (emesh%pactive(i)%pelem%error.gt.error) error = emesh%pactive(i)%pelem%error
       enddo

     else   

       do i = 1, emesh%maxind
         error = error + emesh%pactive(i)%pelem%error
       enddo

       error = error / emesh%maxind

     endif

  end function get_error   


end  module negf_energy_mesh
