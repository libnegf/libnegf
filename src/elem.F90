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


module energy_mesh

  use ln_precision
  implicit none

  public :: elem


  type elem
    integer :: lev
    logical :: active
    real(dp) :: pnt(3)
    integer :: map(3)
    type(elem), pointer :: parent =>null()
    type(elem), pointer :: child(2) =>null()
  end type elem

  type mesh     
     type(elem), dimension(:), pointer :: el0 
  end type mesh        

  integer :: meshmem

contains
 
  subroutine create(emesh,x1,x2,n)
    type(mesh) :: emesh
    real(dp) :: x1,x2
    integer :: n

    real(dp) :: x(n), w(n)
    integer :: nelem, ierr
    type(elem), pointer :: el

    if (mod(n,2).eq.0) n=n+1 

    call trapez(x1,x2,x,w,n)
      
    nelem = (n-1)/2

    allocate(emesh%el0(nelem), stat=ierr)
    if(ierr.ne.0) stop 'ERROR: Cannot allocate mesh'

    k = 1
    do i=1,n,2
      el => emesh%el0(i)
      el%pnt = /(x(i),x(i+1),x(i+2))/
      el%map = /(i,i+1,i+2)/
      el%lev = 0
      el%active = .true.
     
      k = k + 1     
    enddo

     meshmem = 10 * nelem

  end subroutine create 

  !--------------------------------------------------------------------  
  subroutine destroy(emesh)
    type(mesh) :: emesh
    
    integer :: nelem,i,k

    nelem = size(emesh%el0)

    do i=1,nelem
      
      call traverse(emesh%el0(i))

    enddo

    deallocate(emesh%el0)

    meshmem = meshmem - 10 * nelem

  end subroutine destroy
  !--------------------------------------------------------------------  
  
  recursive subroutine traverse(el)
     type(elem), pointer :: el
   
     if (el%active) then 
          call destroy_el(el)
     else        
          call traverse(el%child(1))
          call traverse(el%child(2))
     endif

  end subroutine traverse 
  !--------------------------------------------------------------------  

  subroutine create_children(el,npoints)
    type(elem), pointer :: el
    integer, intent(out) :: npoints

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
    el1%map(3) = el%map(2)
    el1%map(2) = npoints + 1    

    el2%map(1) = el%map(2)
    el2%map(3) = el%map(3)
    el2%map(2) = npoints + 2    
    
    el1%lev = el%lev + 1
    el2%lev = el%lev + 1
 
    el1%active = .true.
    el2%active = .true.
    el%active = .false.

    el1%parent => el
    el2%parent => el

    el%child(1) => el1
    el%child(2) => el2

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

  subroutine test_mesh()
     type(mesh) :: emesh

     integer :: i,npoints

     call create(emesh,-2.77,0,51)

     npoints = 51

     do i =1, size(emesh)

       call create_children(emesh%el(i),npoints)
       
     enddo 

     print*,'meshmem=', meshmem

     call destroy(emesh)

  end subroutine test_mesh


end module energy_mesh
