module random
  use ln_precision    
  implicit none
  private
  public :: init_seed
  public :: zinit_rnd_matrix

  contains
  ! ---------------------------------------------------------
  subroutine init_seed(seed_in)
    integer, intent(in), optional :: seed_in

    integer, dimension(:),allocatable :: seed
    integer :: j

    call random_seed(size=j)
    allocate(seed(j))

    if (present(seed_in)) then
      seed=seed_in
    else
      call system_clock(j)
      seed=j
    end if
    call random_seed(put=seed)
    deallocate(seed)

  end subroutine init_seed

  ! ---------------------------------------------------------
  subroutine zinit_rnd_matrix(Mat)
    complex(dp) :: Mat(:,:)

    integer :: i,j,M,N
    real(dp) :: rnd_r, rnd_i

    M = size(Mat,1) 
    N = size(Mat,2)   
    do j = 1, N
      do i = 1, M 
        call random_number(rnd_r)
        call random_number(rnd_i)
        Mat(i,j) = cmplx(rnd_r, rnd_i)
      end do 
    end do
  
  end subroutine zinit_rnd_matrix

end module random  
