module metis_interface

  use mat_def
  
  implicit none
  private

  public ::  METISpartition


contains

  subroutine METISpartition(mat,nbl,part)
    type(z_CSR) :: mat
    
    integer :: nbl
    integer :: n
    
    integer :: volume
    integer :: numflag
    integer :: wghts
    integer :: wgtflag

    integer, dimension(:) :: part


    integer, dimension(:), allocatable :: options
    integer, dimension(:), allocatable :: vwgt  
    integer, dimension(:), allocatable :: vsize
    

    external METIS_PartGraphVKway

    numflag = 1
    wghts = 0
    wgtflag = 0
    n = mat%nrow

    call log_allocate(vwgt, 0)
    call log_allocate(vsize, 0)
    call log_allocate(options, 5)
    options(1) = 0


    !call METIS_PartGraphVKway(n, mat%rowpnt, mat%colind, vwgt, vsize, wgtflag, &
    !                          numflag, nbl, options, volume, part)

    call METIS_PartGraphKway(n, mat%rowpnt, mat%colind, vwgt, vsize, wgtflag, &
                              numflag, nbl, options, volume, part)    

    call log_deallocate(vwgt)
    call log_deallocate(vsize)
    call log_deallocate(options)
    

  end subroutine METISpartition


  !----------------------------------------------------------

end module metis_interface
