!*****************************************************************************
!
! GENRCM finds the reverse Cuthill-Mckee ordering for a general graph.
!
!  Discussion:
!
!    For each connected component in the graph, the routine obtains
!    an ordering by calling RCM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 January 2003
!
!  Author:
!
!    Original FORTRAN77 version by Alan George, Joseph Liu.
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Alan George, Joseph Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, integer ( kind = 4 ) ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ( kind = 4 ) ADJ_ROW(NODE_NUM+1).  Information about 
!    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ( kind = 4 ) ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Output, integer ( kind = 4 ) PERM(NODE_NUM), the RCM ordering.
!
!  Local Parameters:
!
!    Local, integer LEVEL_ROW(NODE_NUM+1), the index vector for a level
!    structure.  The level structure is stored in the currently unused 
!    spaces in the permutation vector PERM.
!
!    Local, integer MASK(NODE_NUM), marks variables that have been numbered.
!
!  Additional notes:
!  The original collection of soubroutines can be found in
!  http://people.sc.fsu.edu/~jburkardt/f_src/rcm/rcm.f90
!  They have been readapted in this module by A. Pecchia.
!
!  Note: find_root crash for very large matrices because of stack overflow
!        part of the problem has been solved but not all. work in progress       
!
!*****************************************************************************

module rcm_module

  implicit none
  private

!!$  private ::  adj_bandwidth
!!$  private ::  adj_contains_ij
!!$  private ::  adj_insert_ij
!!$  private ::  adj_perm_bandwidth
!!$  private ::  adj_perm_show
!!$  private ::  adj_print
!!$  private ::  adj_print_some
!!$  private ::  adj_set
!!$  private ::  adj_show
!!$  private ::  d_is_int
!!$  private ::  d_uniform_01
!!$  private ::  d2vec_permute
!!$  private ::  degree
!!$  private ::  dmat_print_some
!!$  private ::  dmat_transpose_print_some 
     public  ::  genrcm
!!$  private ::  graph_01_adj
!!$  private ::  graph_01_size
!!$  private ::  graph_01_label    
!!$  private ::  i_swap
!!$  private ::  i_uniform
!!$  private ::  imat_print_some
!!$  private ::  imat_transpose_print
!!$  private ::  imat_transpose_print_some
!!$  private ::  irow_compare
!!$  private ::  irow_sort_a
!!$  private ::  irow_swap
!!$  private ::  ivec_indicator
!!$  private ::  ivec_reverse  
!!$  private ::  ivec_print
!!$  private ::  level_set
!!$  private ::  level_set_print
!!$  private ::  perm_check
!!$  private ::  perm_inverse
!!$  private ::  perm_random
!!$  private ::  rcm
!!$  private ::  root_find
!!$  private ::  sort_heap_external
!!$  private ::  timestamp
!!$  private ::  timestring
!!$  private ::  triangulation_order3_adj_count
!!$  private ::  triangulation_order3_adj_set
!!$  private ::  triangulation_order3_example
!!$  private ::  triangulation_order3_example_size
!!$  private ::  triangulation_order3_neighbor_triangles
!!$  private ::  triangulation_order6_adj_count 
!!$  private ::  triangulation_order6_adj_set
!!$  private ::  triangulation_order6_example
!!$  private ::  triangulation_order6_example_size
!!$  private ::  triangulation_order6_neighbor_triangles

contains
  
function adj_bandwidth ( node_num, adj_num, adj_row, adj )  

!*******************************************************************************
!
!! ADJ_BANDWIDTH computes the bandwidth of an adjacency matrix.
!
!  Modified:
!
!    11 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Output, integer ADJ_BANDWIDTH, the bandwidth of the adjacency
!    matrix.
!
  implicit none

  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_bandwidth
  integer adj_row(node_num+1)
  integer band_hi
  integer band_lo
  integer col
  integer i
  integer j

  band_lo = 0
  band_hi = 0

  do i = 1, node_num

    do j = adj_row(i), adj_row(i+1)-1
      col = adj(j)
      band_lo = max ( band_lo, i - col )
      band_hi = max ( band_hi, col - i )
    end do

  end do

  adj_bandwidth = band_lo + 1 + band_hi

  return
end function adj_bandwidth

logical function adj_contains_ij ( node_num, adj_num, adj_row, adj, i, j )

!*******************************************************************************
!
!! ADJ_CONTAINS_IJ determines if (I,J) is in an adjacency structure.
!
!  Modified:
!
!    23 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!
!    Input, integer I, J, the two nodes, for which we want to know
!    whether I is adjacent to J.
!
!    Output, logical ADJ_CONTAINS_IJ, is TRUE if I = J, or the adjacency
!    structure contains the information that I is adjacent to J.
!
  implicit none

  integer adj_num
  integer node_num

  integer adj(adj_num)
  !logical adj_contains_ij
  integer adj_row(node_num+1)
  integer i
  integer j
  integer k
  integer khi
  integer klo
!
!  Symmetric entries are not stored.
!
  if ( i == j ) then
    adj_contains_ij = .true.
    return
  end if
!
!  Illegal I, J entries.
!
  if ( node_num < i ) then
    adj_contains_ij = .false.
    return
  else if ( i < 1 ) then
    adj_contains_ij = .false.
    return
  else if ( node_num < j ) then
    adj_contains_ij = .false.
    return
  else if ( j < 1 ) then
    adj_contains_ij = .false.
    return
  end if
!
!  Search the adjacency entries already stored for row I,
!  to see if J has already been stored.
!
  klo = adj_row(i)
  khi = adj_row(i+1)-1

  do k = klo, khi

    if ( adj(k) == j ) then
      adj_contains_ij = .true.
      return
    end if

  end do

  adj_contains_ij = .false.

  return
end function adj_contains_ij

subroutine adj_insert_ij ( node_num, adj_max, adj_num, adj_row, adj, i, j )

!*******************************************************************************
!
!! ADJ_INSERT_IJ inserts (I,J) into an adjacency structure.
!
!  Modified:
!
!    23 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ADJ_MAX, the maximum number of adjacency entries.
!
!    Input/output, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!
!    Input, integer I, J, the two nodes which are adjacent.
!
  implicit none

  integer adj_max
  integer node_num

  integer adj(adj_max)
  integer adj_num
  integer adj_row(node_num+1)
  integer i
  integer j
  integer j_spot
  integer k
!
!  A new adjacency entry must be made.
!  Check that we're not exceeding the storage allocation for ADJ.
!
  if ( adj_max < adj_num + 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADJ_INSERT_IJ - Fatal error!'
    write ( *, '(a)' ) '  All available storage has been used.'
    write ( *, '(a)' ) '  No more information can be stored!'
    write ( *, '(a)' ) '  This error occurred for '
    write ( *, '(a,i6)' ) '  Row I = ', i
    write ( *, '(a,i6)' ) '  Column J = ', j
    stop
  end if
!
!  The action is going to occur between ADJ_ROW(I) and ADJ_ROW(I+1)-1:
!
  j_spot = adj_row(i)

  do k = adj_row(i), adj_row(i+1)-1
    if ( adj(k) == j ) then
      return
    else if ( adj(k) < j ) then
      j_spot = k + 1
    else 
      exit
    end if
  end do

  adj(j_spot+1:adj_num+1) = adj(j_spot:adj_num)
  adj(j_spot) = j

  adj_row(i+1:node_num+1) = adj_row(i+1:node_num+1) + 1

  adj_num = adj_num + 1

  return
end subroutine adj_insert_ij

function adj_perm_bandwidth ( node_num, adj_num, adj_row, adj, perm, perm_inv )

!*******************************************************************************
!
!! ADJ_PERM_BANDWIDTH computes the bandwidth of a permuted adjacency matrix.
!
!  Discussion:
!
!    The matrix is defined by the adjacency information and a permutation.  
!
!    The routine also computes the bandwidth and the size of the envelope.
!
!  Modified:
!
!    11 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input, integer PERM(NODE_NUM), integer PERM_INV(NODE_NUM), the permutation
!    and inverse permutation.
!
!    Output, integer ADJ_PERM_BANDWIDTH, the bandwidth of the permuted adjacency
!    matrix.
!
  implicit none

  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_perm_bandwidth
  integer adj_row(node_num+1)
  integer band_hi
  integer band_lo
  integer col
  integer i
  integer j
  integer perm(node_num)
  integer perm_inv(node_num)

  band_lo = 0
  band_hi = 0

  do i = 1, node_num

    do j = adj_row(perm(i)), adj_row(perm(i)+1)-1
      col = perm_inv(adj(j))
      band_lo = max ( band_lo, i - col )
      band_hi = max ( band_hi, col - i )
    end do

  end do

  adj_perm_bandwidth = band_lo + 1 + band_hi

  return
end function adj_perm_bandwidth

subroutine adj_perm_show ( node_num, adj_num, adj_row, adj, perm, perm_inv )

!*******************************************************************************
!
!! ADJ_PERM_SHOW displays a symbolic picture of a permuted adjacency matrix.
!
!  Discussion:
!
!    The matrix is defined by the adjacency information and a permutation.  
!
!    The routine also computes the bandwidth and the size of the envelope.
!
!    If no permutation has been done, you must set PERM(I) = PERM_INV(I) = I
!    before calling this routine.  
!
!  Modified:
!
!    28 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input, integer PERM(NODE_NUM), integer PERM_INV(NODE_NUM), the permutation
!    and inverse permutation.
!
  implicit none

  integer, parameter :: n_max = 100

  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_row(node_num+1)
  character band(n_max)
  integer band_lo
  integer col
  integer i
  integer j
  integer k
  integer nonzero_num
  integer perm(node_num)
  integer perm_inv(node_num)

  band_lo = 0
  nonzero_num = 0

  if ( n_max < node_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADJ_PERM_SHOW - Fatal error!'
    write ( *, '(a)' ) '  NODE_NUM is too large!'
    write ( *, '(a,i6)' ) '  Maximum legal value is ', n_max
    write ( *, '(a,i6)' ) '  Your input value was ', node_num
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Nonzero structure of matrix:'
  write ( *, '(a)' ) ' '

  do i = 1, node_num

    do k = 1, node_num
      band(k) = '.'
    end do

    band(i) = 'D'

    do j = adj_row(perm(i)), adj_row(perm(i)+1)-1
      col = perm_inv(adj(j))
      if ( col < i ) then
        nonzero_num = nonzero_num + 1
      end if
      band_lo = max ( band_lo, i-col )
      if ( col /= i ) then
        band(col) = 'X'
      end if
    end do

    write ( *, '(2x,i6,1x,100a1)' ) i, band(1:node_num)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Lower bandwidth = ', band_lo
  write ( *, '(a,i6,a)' ) '  Lower envelope contains ', &
    nonzero_num, ' nonzeros.'

  return
end subroutine adj_perm_show

subroutine adj_print ( node_num, adj_num, adj_row, adj, title )

!*******************************************************************************
!
!! ADJ_PRINT prints adjacency information.
!
!  Discussion:
!
!    The list has the form:
!
!    Row   Nonzeros
!
!    1       2   5   9
!    2       7   8   9   15   78   79   81  86  91  99
!          100 103
!    3      48  49  53
!
!  Modified:
!
!    18 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1), organizes the adjacency entries
!    into rows.  The entries for row I are in entries ADJ_ROW(I)
!    through ADJ_ROW(I+1)-1.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure, which contains,
!    for each row, the column indices of the nonzero entries.
!
!    Input, character ( len = * ) TITLE, a title to be printed.
!
  implicit none

  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_row(node_num+1)
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  call adj_print_some ( node_num, 1, node_num, adj_num, adj_row, adj )

  return
end subroutine adj_print

subroutine adj_print_some ( node_num, node_lo, node_hi, adj_num, adj_row, adj )

!*******************************************************************************
!
!! ADJ_PRINT_SOME prints some adjacency information.
!
!  Discussion:
!
!    The list has the form:
!
!    Row   Nonzeros
!
!    1       2   5   9
!    2       7   8   9   15   78   79   81  86  91  99
!          100 103
!    3      48  49  53
!
!  Modified:
!
!    18 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer NODE_LO, NODE_HI, the first and last nodes for
!    which the adjacency information is to be printed.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1), organizes the adjacency entries
!    into rows.  The entries for row I are in entries ADJ_ROW(I)
!    through ADJ_ROW(I+1)-1.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure, which contains,
!    for each row, the column indices of the nonzero entries.
!
  implicit none

  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_row(node_num+1)
  integer i
  integer jhi
  integer jlo
  integer jmax
  integer jmin
  integer node_hi
  integer node_lo

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Sparse adjacency structure:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of nodes       = ', node_num
  write ( *, '(a,i6)' ) '  Number of adjacencies = ', adj_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node Min Max      Nonzeros '
  write ( *, '(a)' ) ' '

  do i = node_lo, node_hi

    jmin = adj_row(i)
    jmax = adj_row(i+1) - 1

    if ( jmax < jmin ) then

      write ( *, '(2x,3i4,6x,10i6)' ) i, jmin, jmax

    else

      do jlo = jmin, jmax, 10

        jhi = min ( jlo + 9, jmax )

        if ( jlo == jmin ) then
          write ( *, '(2x,3i4,3x,10i6)' ) i, jmin, jmax, adj(jlo:jhi)
        else
          write ( *, '(2x,12x,3x,10i6)' )                adj(jlo:jhi)
        end if

      end do

    end if

  end do

  return
end subroutine adj_print_some

subroutine adj_set ( node_num, adj_max, adj_num, adj_row, adj, irow, jcol )

!*******************************************************************************
!
!! ADJ_SET sets up the adjacency information.
!
!  Discussion:
!
!    The routine records the locations of each nonzero element,
!    one at a time.
!
!    The first call for a given problem should be with IROW or ICOL
!    negative.  This is a signal indicating the data structure should
!    be initialized.
!
!    Then, for each case in which A(IROW,JCOL) is nonzero, or
!    in which IROW is adjacent to JCOL, call this routine once
!    to record that fact.
!
!    Diagonal entries are not to be stored.
!
!    The matrix is assumed to be symmetric, so setting I adjacent to J
!    will also set J adjacent to I.
!
!    Repeated calls with the same values of IROW and JCOL do not
!    actually hurt.  No extra storage will be allocated.
!
!  Modified:
!
!    23 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ADJ_MAX, the maximum dimension of the adjacency array.
!
!    Input/output, integer ADJ_NUM, the number of adjaceny entries.
!
!    Input/output, integer ADJ_ROW(NODE_NUM+1).  Information about 
!    row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input/output, integer ADJ(ADJ_NUM), the adjacency structure.
!
!    Input, integer IROW, JCOL, the row and column indices of a nonzero
!    entry of the matrix.
!
  implicit none

  integer adj_max
  integer node_num

  integer adj(adj_max)
  !logical adj_contains_ij
  integer adj_num
  integer adj_row(node_num+1)
  integer irow
  integer jcol
!
!  Negative IROW or JCOL indicates the data structure should be initialized.
!
  if ( irow < 0 .or. jcol < 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADJ_SET - Note:'
    write ( *, '(a)') '  Initializing adjacency information.'
    write ( *, '(a,i6)' ) '  Number of nodes NODE_NUM = ', node_num
    write ( *, '(a,i6)' ) '  Maximum adjacency ADJ_MAX = ', adj_max

    adj_num = 0
    adj_row(1:node_num+1) = 1
    adj(1:adj_max) = 0

    return

  end if
!
!  Diagonal entries are not stored.
!
  if ( irow == jcol ) then
    return
  end if

  if ( node_num < irow ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADJ_SET - Fatal error!'
    write ( *, '(a)' ) '  NODE_NUM < IROW.'
    write ( *, '(a,i6)' ) '  IROW = ', irow
    write ( *, '(a,i6)' ) '  NODE_NUM = ', node_num
    stop
  else if ( irow < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADJ_SET - Fatal error!'
    write ( *, '(a)' ) '  IROW < 1.'
    write ( *, '(a,i6)' ) '  IROW = ', irow
    stop
  else if ( node_num < jcol ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADJ_SET - Fatal error!'
    write ( *, '(a)' ) '  NODE_NUM < JCOL.'
    write ( *, '(a,i6)' ) '  JCOL = ', jcol
    write ( *, '(a,i6)' ) '  NODE_NUM = ', node_num
    stop
  else if ( jcol < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADJ_SET - Fatal error!'
    write ( *, '(a)' ) '  JCOL < 1.'
    write ( *, '(a,i6)' ) '  JCOL = ', jcol
    stop
  end if

  if ( .not. &
    adj_contains_ij ( node_num, adj_num, adj_row, adj, irow, jcol ) ) then
    call adj_insert_ij ( node_num, adj_max, adj_num, adj_row, adj, irow, jcol )
  end if

  if ( .not. &
    adj_contains_ij ( node_num, adj_num, adj_row, adj, jcol, irow ) ) then
    call adj_insert_ij ( node_num, adj_max, adj_num, adj_row, adj, jcol, irow )
  end if

  return
end subroutine adj_set

subroutine adj_show ( node_num, adj_num, adj_row, adj )

!*******************************************************************************
!
!! ADJ_SHOW displays a symbolic picture of an adjacency matrix.
!
!  Discussion:
!
!    The matrix is defined by the adjacency information and a permutation.  
!
!    The routine also computes the bandwidth and the size of the envelope.
!
!  Modified:
!
!    11 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
  implicit none

  integer, parameter :: n_max = 100

  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_row(node_num+1)
  character band(n_max)
  integer band_lo
  integer col
  integer i
  integer j
  integer k
  integer nonzero_num

  band_lo = 0
  nonzero_num = 0

  if ( n_max < node_num ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ADJ_SHOW - Fatal error!'
    write ( *, '(a)' ) '  NODE_NUM is too large!'
    write ( *, '(a,i6)' ) '  Maximum legal value is ', n_max
    write ( *, '(a,i6)' ) '  Your input value was ', node_num
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Nonzero structure of matrix:'
  write ( *, '(a)' ) ' '

  do i = 1, node_num

    do k = 1, node_num
      band(k) = '.'
    end do

    band(i) = 'D'

    do j = adj_row(i), adj_row(i+1)-1
      col = adj(j)
      if ( col < i ) then
        nonzero_num = nonzero_num + 1
      end if
      band_lo = max ( band_lo, i-col )
      band(col) = 'X'
    end do

    write ( *, '(2x,i6,1x,100a1)' ) i, band(1:node_num)

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Lower bandwidth = ', band_lo
  write ( *, '(a,i6,a)' ) '  Lower envelope contains ', &
    nonzero_num, ' nonzeros.'

  return
end subroutine adj_show

logical function d_is_int ( r )

!*******************************************************************************
!
!! D_IS_INT determines if a double precision number represents an integer value.
!
!  Modified:
!
!    07 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the number to be checked.
!
!    Output, logical D_IS_INT, is TRUE if R is an integer value.
!
  implicit none

  integer i
  real ( kind = 8 ) r

  if ( real ( huge ( i ), kind = 8 ) < r ) then
    d_is_int = .false.
  else if ( r < - real ( huge ( i ), kind = 8 ) ) then
    d_is_int = .false.
  else if ( r == real ( int ( r ), kind = 8 ) ) then
    d_is_int = .true.
  else
    d_is_int = .false.
  end if

  return
end function d_is_int

  real ( kind = 8 ) function d_uniform_01 ( seed )

!*******************************************************************************
!
!! D_UNIFORM_01 is a portable pseudorandom number generator.
!
!  Discussion:
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      unif = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!  Modified:
!
!    11 August 2004
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, L E Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!  Parameters:
!
!    Input/output, integer SEED, the "seed" value, which should NOT be 0.
!    (Otherwise, the output values of SEED and UNIFORM will be zero.)
!    On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) D_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer k
  integer seed
  !real ( kind = 8 ) d_uniform_01

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  d_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end function d_uniform_01

subroutine d2vec_permute ( n, a, p )

!*******************************************************************************
!
!! D2VEC_PERMUTE permutes a D2 vector in place.
!
!  Discussion:
!
!    This routine permutes an array of real "objects", but the same
!    logic can be used to permute an array of objects of any arithmetic
!    type, or an array of objects of any complexity.  The only temporary
!    storage required is enough to store a single object.  The number
!    of data movements made is N + the number of cycles of order 2 or more,
!    which is never more than N + N/2.
!
!  Example:
!
!    Input:
!
!      N = 5
!      P = (   2,    4,    5,    1,    3 )
!      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
!          (11.0, 22.0, 33.0, 44.0, 55.0 )
!
!    Output:
!
!      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
!             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
!
!  Modified:
!
!    13 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of objects.
!
!    Input/output, real ( kind = 8 ) A(2,N), the array to be permuted.
!
!    Input, integer P(N), the permutation.  P(I) = J means
!    that the I-th element of the output array should be the J-th
!    element of the input array.  P must be a legal permutation
!    of the integers from 1 to N, otherwise the algorithm will
!    fail catastrophically.
!
  implicit none

  integer n
  integer, parameter :: ndim = 2

  real ( kind = 8 ) a(ndim,n)
  real ( kind = 8 ) a_temp(ndim)
  integer ierror
  integer iget
  integer iput
  integer istart
  integer p(n)

  call perm_check ( n, p, ierror )

  if ( ierror /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'D2VEC_PERMUTE - Fatal error!'
    write ( *, '(a)' ) '  The input array does not represent'
    write ( *, '(a)' ) '  a proper permutation.  In particular, the'
    write ( *, '(a,i6)' ) '  array is missing the value ', ierror
    stop
  end if
!
!  Search for the next element of the permutation that has not been used.
!
  do istart = 1, n

    if ( p(istart) < 0 ) then

      cycle

    else if ( p(istart) == istart ) then

      p(istart) = -p(istart)
      cycle

    else

      a_temp(1:ndim) = a(1:ndim,istart)
      iget = istart
!
!  Copy the new value into the vacated entry.
!
      do

        iput = iget
        iget = p(iget)

        p(iput) = -p(iput)

        if ( iget < 1 .or. n < iget ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'D2VEC_PERMUTE - Fatal error!'
          write ( *, '(a)' ) '  A permutation index is out of range.'
          write ( *, '(a,i6,a,i6)' ) '  P(', iput, ') = ', iget
          stop
        end if

        if ( iget == istart ) then
          a(1:ndim,iput) = a_temp(1:ndim)
          exit
        end if

        a(1:ndim,iput) = a(1:ndim,iget)

      end do

    end if

  end do
!
!  Restore the signs of the entries.
!
  p(1:n) = -p(1:n)

  return
end subroutine d2vec_permute
  

subroutine degree(root, adj_num, adj_row, adj, mask, deg, iccsze, ls, node_num)

!*******************************************************************************
!
!! DEGREE computes the degrees of the nodes in the connected component.
!
!  Discussion:
!
!    The connected component is specified by MASK and ROOT.
!    Nodes for which MASK is zero are ignored.
!
!  Modified:
!
!   05 January 2003
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ROOT, the node that defines the connected component.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input, integer MASK(NODE_NUM), is nonzero for those nodes which are
!    to be considered.
!
!    Output, integer DEG(NODE_NUM), contains, for each  node in the connected
!    component, its degree.
!
!    Output, integer ICCSIZE, the number of nodes in the connected component.
!
!    Output, integer LS(NODE_NUM), stores in entries 1 through ICCSIZE the nodes
!    in the connected component, starting with ROOT, and proceeding 
!    by levels.
!
!    Input, integer NODE_NUM, the number of nodes.
!
  implicit none

  integer root
  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_row(node_num+1)
  integer deg(node_num)
  integer mask(node_num)
  integer ls(node_num)

  integer i
  integer iccsze
  integer ideg
  integer j
  integer jstop
  integer jstrt
  integer lbegin
  integer lvlend
  integer lvsize
  integer nbr
  integer node
!
!  The sign of ADJ_ROW(I) is used to indicate if node I has been considered.
!
  ls(1) = root
  adj_row(root) = -adj_row(root)
  lvlend = 0
  iccsze = 1
!
!  LBEGIN is the pointer to the beginning of the current level, and
!  LVLEND points to the end of this level.
!
  do

    lbegin = lvlend + 1
    lvlend = iccsze
!
!  Find the degrees of nodes in the current level,
!  and at the same time, generate the next level.
!
    do i = lbegin, lvlend

      node = ls(i)
      jstrt = -adj_row(node)
      jstop = abs ( adj_row(node+1) ) - 1
      ideg = 0

      do j = jstrt, jstop

        nbr = adj(j)

        if ( mask(nbr) /= 0 ) then

          ideg = ideg + 1

          if ( 0 <= adj_row(nbr) ) then
            adj_row(nbr) = -adj_row(nbr)
            iccsze = iccsze + 1
            ls(iccsze) = nbr
          end if

        end if

      end do

      deg(node) = ideg

    end do
!
!  Compute the current level width.
!
    lvsize = iccsze - lvlend
!
!  If the current level width is nonzero, generate another level.
!
    if ( lvsize == 0 ) then
      exit
    end if

  end do
!
!  Reset ADJ_ROW to its correct sign and return.
!
  do i = 1, iccsze
    node = ls(i)
    adj_row(node) = -adj_row(node)
  end do

  return
end subroutine degree

subroutine dmat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*******************************************************************************
!
!! DMAT_PRINT_SOME prints some of a double precision matrix.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer, parameter :: incx = 5
  integer m
  integer n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  !logical, external :: d_is_int
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7,7x)') j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( d_is_int ( a(i,j) ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end subroutine dmat_print_some

subroutine dmat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*******************************************************************************
!
!! DMAT_TRANSPOSE_PRINT_SOME prints some of a DMAT, transposed.
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer, parameter :: incx = 5
  integer m
  integer n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i7,7x)') i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,1x,5a14)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end subroutine dmat_transpose_print_some

subroutine genrcm ( node_num, adj_num, adj_row, adj, perm )

!*******************************************************************************
!
!! GENRCM finds the reverse Cuthill-Mckee ordering for a general graph.
!
!  Discussion:
!
!    For each connected component in the graph, the routine obtains
!    an ordering by calling RCM.
!
!  Modified:
!
!    04 January 2003
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Output, integer PERM(NODE_NUM), the RCM ordering.
!
!  Local Parameters:
!
!    Local, integer LEVEL_ROW(NODE_NUM+1), the index vector for a level
!    structure.  The level structure is stored in the currently unused 
!    spaces in the permutation vector PERM.
!
!    Local, integer MASK(NODE_NUM), marks variables that have been numbered.
!
  implicit none

  integer :: adj_num
  integer :: node_num

  integer, dimension(:) ::  adj
  integer, dimension(:) ::  adj_row
  integer ::  perm(node_num)
  ! locals
  integer ::  i
  integer ::  iccsze
  integer, dimension(:), allocatable ::  mask, level_row
  integer :: level_num, level, num, root, err

  allocate(mask(node_num),stat=err)
  allocate(level_row(node_num+1),stat=err)
  if (err.ne.0) STOP 'ERROR: RCM allocation error'

  mask(1:node_num) = 1

  num = 1

  do i = 1, node_num
!
!  For each masked connected component...
!
    if ( mask(i) /= 0 ) then

      root = i
!
!  Find a pseudo-peripheral node ROOT.  The level structure found by
!  ROOT_FIND is stored starting at PERM(NUM).
!  root_find ( root, adj_num, adj_row, adj, mask, level_num, &
!                   level_row, level, node_num )
   call root_find(root, adj_num, adj_row, adj, mask, level_num, &
                  level_row, perm(num:), node_num )

!
!  RCM orders the component using ROOT as the starting node.
!
   call rcm(root, adj_num, adj_row, adj, mask, perm(num:), iccsze, &
            node_num )
   
   num = num + iccsze
!
!  We can stop once every node is in one of the connected components.
!
      if ( node_num < num ) then
          return
         end if
   
    end if

  end do

  deallocate(mask)
  deallocate(level_row)

end subroutine genrcm

subroutine graph_01_adj ( node_num, adj_num, adj_row, adj )

!*******************************************************************************
!
!! GRAPH_01_ADJ returns the adjacency vector for graph 1.
!
!  Modified:
!
!    22 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer ADJ_NUM, the number of adjacencies.
!
!    Output, integer ADJ_ROW(NODE_NUM+1), node pointers into ADJ.
!
!    Output, integer ADJ(ADJ_NUM), the adjacency information.
!
  implicit none

  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_row(node_num+1)

  adj(1:adj_num) = (/ &
    4, 6, &
    3, 5, 7, 10, &
    2, 4, 5, &
    1, 3, 6, 9, &
    2, 3, 7, &
    1, 4, 7, 8, &
    2, 5, 6, 8, &
    6, 7, &
    4, &
    2 /)

  adj_row(1:node_num+1) = (/ 1, 3, 7, 10, 14, 17, 21, 25, 27, 28, 29 /)

  return
end subroutine graph_01_adj

subroutine graph_01_size ( node_num, adj_num )

!*******************************************************************************
!
!! GRAPH_01_ADJ_NUM returns the number of adjacencies for graph 1.
!
!  Modified:
!
!    22 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Output, integer NODE_NUM, the number of items that can be adjacent.
!
!    Output, integer ADJ_NUM, the number of adjacencies.
!
  implicit none

  integer adj_num
  integer node_num

  node_num = 10
  adj_num = 28

  return
end subroutine graph_01_size

subroutine graph_01_label ( node_num, label )

!*******************************************************************************
!
!! GRAPH_01_LABEL returns the labels for graph 1.
!
!  Modified:
!
!    22 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Output, character LABEL(NODE_NUM), the labels of the nodes.
!
  implicit none

  integer node_num

  character label(node_num)

  label(1:node_num) = (/ 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J' /)

  return
end subroutine graph_01_label

subroutine i_swap ( i, j )

!*******************************************************************************
!
!! I_SWAP swaps two integer values.
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer I, J.  On output, the values of I and
!    J have been interchanged.
!
  implicit none

  integer i
  integer j
  integer k

  k = i
  i = j
  j = k

  return
end subroutine i_swap

  integer function i_uniform ( a, b, seed )

!*******************************************************************************
!
!! I_UNIFORM returns a pseudorandom integer.
!
!  Discussion:
!
!    The pseudorandom number should be uniformly distributed
!    between A and B.
!
!  Modified:
!
!    21 November 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, B, the limits of the interval.
!
!    Input/output, integer SEED, the "seed" value, which should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, integer I_UNIFORM, a number between A and B.
!
  implicit none

  integer a
  integer b
  real ( kind = 8 ) d
  !real ( kind = 8 ) d_uniform_01
  !integer i_uniform
  integer seed

  d = real ( a, kind = 8 ) - 0.5D+00 &
    + real ( 1 + b - a, kind = 8 ) * d_uniform_01 ( seed )

  i_uniform = nint ( d )

  i_uniform = max ( i_uniform, a )
  i_uniform = min ( i_uniform, b )

  return
end function i_uniform

subroutine imat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*******************************************************************************
!
!! IMAT_PRINT_SOME prints some of an integer matrix.
!
!  Modified:
!
!    04 November 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer, parameter :: incx = 10
  integer m
  integer n

  integer a(m,n)
  character ( len = 7 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i7)') j
    end do

    write ( *, '(''  Col '',10a7)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(i7)' ) a(i,j)

      end do

      write ( *, '(i5,1x,10a7)' ) i, ( ctemp(j), j = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end subroutine imat_print_some

subroutine imat_transpose_print ( m, n, a, title )

!*******************************************************************************
!
!! IMAT_TRANSPOSE_PRINT prints an IMAT, transposed.
!
!  Modified:
!
!    28 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  character ( len = * ) title

  call imat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end subroutine imat_transpose_print

subroutine imat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*******************************************************************************
!
!! IMAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an integer matrix.
!
!  Modified:
!
!    09 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer, parameter :: incx = 10
  integer m
  integer n

  integer a(m,n)
  character ( len = 7 ) ctemp(incx)
  integer i
  integer i2
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i7)') i
    end do

    write ( *, '(''  Row '',10a7)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc

        i = i2lo - 1 + i2

        write ( ctemp(i2), '(i7)' ) a(i,j)

      end do

      write ( *, '(i5,1x,10a7)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  write ( *, '(a)' ) ' '

  return
end subroutine imat_transpose_print_some

subroutine irow_compare ( m, n, a, i, j, isgn )

!*******************************************************************************
!
!! IROW_COMPARE compares two rows of a integer array.
!
!  Example:
!
!    Input:
!
!  M = 3, N = 4, I = 2, J = 3
!
!  A = (
!    1  2  3  4
!    5  6  7  8
!    9 10 11 12 )
!
!    Output:
!
!  ISGN = -1
!
!  Modified:
!
!    14 November 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, integer A(M,N), an array of M rows of vectors of length N.
!
!    Input, integer I, J, the rows to be compared.
!    I and J must be between 1 and M.
!
!    Output, integer ISGN, the results of the comparison:
!    -1, row I < row J,
!     0, row I = row J,
!    +1, row J < row I.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer i
  integer isgn
  integer j
  integer k
!
!  Check that I and J are legal.
!
  if ( i < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index I is less than 1.'
    write ( *, '(a,i6)' ) '  I = ', i
    stop
  else if ( m < i ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index I is out of bounds.'
    write ( *, '(a,i6)' ) '  I = ', i
    write ( *, '(a,i6)' ) '  Maximum legal value is M = ', m
    stop
  end if

  if ( j < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index J is less than 1.'
    write ( *, '(a,i6)' ) '  J = ', j
    stop
  else if ( m < j ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IROW_COMPARE - Fatal error!'
    write ( *, '(a)' ) '  Row index J is out of bounds.'
    write ( *, '(a,i6)' ) '  J = ', j
    write ( *, '(a,i6)' ) '  Maximum legal value is M = ', m
    stop
  end if

  isgn = 0

  if ( i == j ) then
    return
  end if

  k = 1

  do while ( k <= n )

    if ( a(i,k) < a(j,k) ) then
      isgn = -1
      return
    else if ( a(j,k) < a(i,k) ) then
      isgn = +1
      return
    end if

    k = k + 1

  end do

  return
end subroutine irow_compare

subroutine irow_sort_a ( m, n, a )

!*******************************************************************************
!
!! IROW_SORT_A ascending sorts the rows of an integer array.
!
!  Discussion:
!
!    In lexicographic order, the statement "X < Y", applied to two
!    vectors X and Y of length M, means that there is some index I, with
!    1 <= I <= M, with the property that
!
!      X(J) = Y(J) for J < I,
!    and
!      X(I) < Y(I).
!
!    In other words, X is less than Y if, at the first index where they
!    differ, the X value is less than the Y value.
!
!  Example:
!
!    Input:
!
!      M = 5, N = 3
!
!      A =
!        3  2  1
!        2  4  3
!        3  1  8
!        2  4  2
!        1  9  9
!
!    Output:
!
!      A =
!        1  9  9
!        2  4  2
!        2  4  3
!        3  1  8
!        3  2  1
!
!  Modified:
!
!    25 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of rows of A.
!
!    Input, integer N, the number of columns of A.
!
!    Input/output, integer A(M,N).
!    On input, the array of M rows of N-vectors.
!    On output, the rows of A have been sorted in ascending
!    lexicographic order.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer i
  integer indx
  integer isgn
  integer j

  if ( m <= 1 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if
!
!  Initialize.
!
  i = 0
  indx = 0
  isgn = 0
  j = 0
!
!  Call the external heap sorter.
!
  do

    call sort_heap_external ( m, indx, i, j, isgn )
!
!  Interchange the I and J objects.
!
    if ( 0 < indx ) then

      call irow_swap ( m, n, a, i, j )
!
!  Compare the I and J objects.
!
    else if ( indx < 0 ) then

      call irow_compare ( m, n, a, i, j, isgn )

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end subroutine irow_sort_a

subroutine irow_swap ( m, n, a, irow1, irow2 )

!*******************************************************************************
!
!! IROW_SWAP swaps two rows of an integer array.
!
!  Modified:
!
!    04 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input/output, integer A(M,N), an array of data.
!
!    Input, integer IROW1, IROW2, the two rows to swap.
!
  implicit none

  integer m
  integer n

  integer a(m,n)
  integer irow1
  integer irow2
  integer row(n)
!
!  Check.
!
  if ( irow1 < 1 .or. m < irow1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  IROW1 is out of range.'
    stop
  end if

  if ( irow2 < 1 .or. m < irow2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'IROW_SWAP - Fatal error!'
    write ( *, '(a)' ) '  IROW2 is out of range.'
    stop
  end if

  if ( irow1 == irow2 ) then
    return
  end if

  row(1:n) = a(irow1,1:n)
  a(irow1,1:n) = a(irow2,1:n)
  a(irow2,1:n) = row(1:n)

  return
end subroutine irow_swap

subroutine ivec_indicator ( n, a )

!*******************************************************************************
!
!! IVEC_INDICATOR sets an integer vector to the vector A(I)=I.
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of elements of A.
!
!    Output, integer A(N), the array to be initialized.
!
  implicit none

  integer n

  integer a(n)
  integer i

  do i = 1, n
    a(i) = i
  end do

  return
end subroutine ivec_indicator

subroutine ivec_print ( n, a, title )

!*******************************************************************************
!
!! IVEC_PRINT prints an integer vector.
!
!  Modified:
!
!    28 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, integer A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title to be printed first.
!    TITLE may be blank.
!
  implicit none

  integer n

  integer a(n)
  integer big
  integer i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  big = maxval ( abs ( a(1:n) ) )

  write ( *, '(a)' ) ' '
  if ( big < 1000 ) then
    do i = 1, n
      write ( *, '(2x,i6,2x,i4)' ) i, a(i)
    end do
  else if ( big < 1000000 ) then
    do i = 1, n
      write ( *, '(2x,i6,2x,i7)' ) i, a(i)
    end do
  else
    do i = 1, n
      write ( *, '(2x,i6,2x,i12)' ) i, a(i)
    end do
  end if

  return
end subroutine ivec_print

subroutine ivec_reverse ( n, a )

!*******************************************************************************
!
!! IVEC_REVERSE reverses the elements of an integer vector.
!
!  Example:
!
!    Input:
!
!      N = 5,
!      A = ( 11, 12, 13, 14, 15 ).
!
!    Output:
!
!      A = ( 15, 14, 13, 12, 11 ).
!
!  Modified:
!
!    26 July 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the array.
!
!    Input/output, integer A(N), the array to be reversed.
!
  implicit none

  integer n

  integer a(n)
  integer i

  do i = 1, n/2
    call i_swap ( a(i), a(n+1-i) )
  end do

  return
end subroutine ivec_reverse

subroutine level_set ( root, adj_num, adj_row, adj, mask, level_num, &
  level_row, level, node_num )

!*******************************************************************************
!
!! LEVEL_SET generates the connected level structure rooted at a given node.
!
!  Discussion:
!
!    Only nodes for which MASK is nonzero will be considered.
!
!    The root node chosen by the user is assigned level 1, and masked.
!    All (unmasked) nodes reachable from a node in level 1 are
!    assigned level 2 and masked.  The process continues until there
!    are no unmasked nodes adjacent to any node in the current level.
!    The number of levels may vary between 2 and NODE_NUM.
!
!  Modified:
!
!    28 October 2003
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ROOT, the node at which the level structure
!    is to be rooted.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input/output, integer MASK(NODE_NUM).  On input, only nodes with nonzero
!    MASK are to be processed.  On output, those nodes which were included
!    in the level set have MASK set to 1.
!
!    Output, integer LEVEL_NUM, the number of levels in the level
!    structure.  ROOT is in level 1.  The neighbors of ROOT
!    are in level 2, and so on.
!
!    Output, integer LEVEL_ROW(NODE_NUM+1), LEVEL(NODE_NUM), the rooted 
!    level structure.
!
!    Input, integer NODE_NUM, the number of nodes.
!
  implicit none

  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_row(node_num+1)
  integer i
  integer iccsze
  integer j
  integer jstop
  integer jstrt
  integer lbegin
  integer level_num
  integer level_row(node_num+1)
  integer level(node_num)
  integer lvlend
  integer lvsize
  integer mask(node_num)
  integer nbr
  integer node
  integer root
  
  mask(root) = 0
  level(1) = root
  level_num = 0
  lvlend = 0
  iccsze = 1
!
!  LBEGIN is the pointer to the beginning of the current level, and
!  LVLEND points to the end of this level.
!

  do

    lbegin = lvlend + 1
    lvlend = iccsze
    level_num = level_num + 1
    level_row(level_num) = lbegin
!
!  Generate the next level by finding all the masked neighbors of nodes
!  in the current level.
!
    do i = lbegin, lvlend

      node = level(i)
      jstrt = adj_row(node)
      jstop = adj_row(node+1)-1
      
      do j = jstrt, jstop

        nbr = adj(j)
        
        if ( mask(nbr).ne.0 ) then
          if(iccsze.lt.node_num) then 
             iccsze = iccsze + 1
             level(iccsze) = nbr
             mask(nbr) = 0
          endif
        end if
      end do

    end do
!
!  Compute the current level width (the number of nodes encountered.)
!  If it is positive, generate the next level.
!
    lvsize = iccsze - lvlend

    if ( lvsize <= 0 ) then
      exit
    end if

  end do

! --------------------------------------------------------
  level_row(level_num+1) = lvlend + 1
! Print data
!
!  Reset MASK to 1 for the nodes in the level structure.
!
  mask(level(1:iccsze)) = 1
       
  return
end subroutine level_set

subroutine level_set_print ( node_num, level_num, level_row, level )

!*******************************************************************************
!
!! LEVEL_SET_PRINT prints level set information.
!
!  Modified:
!
!    26 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer LEVEL_NUM, the number of levels.
!
!    Input, integer LEVEL_ROW(LEVEL_NUM+1), organizes the entries of LEVEL.
!    The entries for level I are in entries LEVEL_ROW(I)
!    through LEVEL_ROW(I+1)-1.
!
!    Input, integer LEVEL(NODE_NUM), is simply a list of the nodes in an
!    order induced by the levels.
!
  implicit none

  integer level_num
  integer node_num

  integer level(node_num)
  integer level_row(level_num+1)
  integer i
  integer jhi
  integer jlo
  integer jmax
  integer jmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEVEL_PRINT'
  write ( *, '(a)' ) '  Show the level set structure of a rooted graph.'
  write ( *, '(a,i6)' ) '  The number of nodes is ', node_num
  write ( *, '(a,i6)' ) '  The number of levels is ', level_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Level Min Max      Nonzeros '
  write ( *, '(a)' ) ' '

  do i = 1, level_num

    jmin = level_row(i)
    jmax = level_row(i+1) - 1

    if ( jmax < jmin ) then

      write ( *, '(2x,3i4,6x,10i6)' ) i, jmin, jmax

    else

      do jlo = jmin, jmax, 10

        jhi = min ( jlo + 9, jmax )

        if ( jlo == jmin ) then
          write ( *, '(2x,3i4,3x,10i6)' ) i, jmin, jmax, level(jlo:jhi)
        else
          write ( *, '(2x,12x,3x,10i6)' )                level(jlo:jhi)
        end if

      end do

    end if

  end do

  return
end subroutine level_set_print

subroutine perm_check ( n, p, ierror )

!*******************************************************************************
!
!! PERM_CHECK checks that a vector represents a permutation.
!
!  Discussion:
!
!    The routine verifies that each of the integers from 1
!    to N occurs among the N entries of the permutation.
!
!  Modified:
!
!    01 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries.
!
!    Input, integer P(N), the array to check.
!
!    Output, integer IERROR, error flag.
!    0, the array represents a permutation.
!    nonzero, the array does not represent a permutation.  The smallest
!    missing value is equal to IERROR.
!
  implicit none

  integer n

  integer ierror
  integer ifind
  integer iseek
  integer p(n)

  ierror = 0

  do iseek = 1, n

    ierror = iseek

    do ifind = 1, n
      if ( p(ifind) == iseek ) then
        ierror = 0
        exit
      end if
    end do

    if ( ierror /= 0 ) then
      return
    end if

  end do

  return
end subroutine perm_check

subroutine perm_inverse ( n, perm, perm_inv )

!*******************************************************************************
!
!! PERM_INVERSE produces the inverse of a given permutation.
!
!  Modified:
!
!    28 October 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of items permuted.
!
!    Input, integer PERM(N), a permutation.
!
!    Output, integer PERM_INV(N), the inverse permutation.
!
  implicit none

  integer n

  integer i
  integer perm(n)
  integer perm_inv(n)

  do i = 1, n
    perm_inv(perm(i)) = i
  end do

  return
end subroutine perm_inverse

subroutine perm_random ( n, seed, p )

!*******************************************************************************
!
!! PERM_RANDOM selects a random permutation of N objects.
!
!  Discussion:
!
!    The routine assumes the objects are labeled 1, 2, ... N.
!
!  Modified:
!
!    12 May 2002
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the number of objects to be permuted.
!
!    Input/output, integer SEED, a seed for the random number generator.
!
!    Output, integer P(N), a permutation of ( 1, 2, ..., N ), in standard
!    index form.
!
  implicit none

  integer n

  integer i
  !integer i_uniform
  integer j
  integer p(n)
  integer seed

  call ivec_indicator ( n, p )

  do i = 1, n
    j = i_uniform ( i, n, seed )
    call i_swap ( p(i), p(j) )
  end do

  return
end subroutine perm_random

subroutine rcm ( root, adj_num, adj_row, adj, mask, perm, iccsze, node_num )

!*******************************************************************************
!
!! RCM renumbers a connected component by the reverse Cuthill McKee algorithm.
!
!  Discussion:
!
!    The connected component is specified by a node ROOT and a mask.
!    The numbering starts at the root node.
!
!    An outline of the algorithm is as follows:
!
!    X(1) = ROOT.
!
!    for ( I = 1 to N-1)
!      Find all unlabeled neighbors of X(I),
!      assign them the next available labels, in order of increasing degree.
!
!    When done, reverse the ordering.
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!  Parameters:
!
!    Input, integer ROOT, the node that defines the connected component.
!    It is used as the starting point for the RCM ordering.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input/output, integer MASK(NODE_NUM), a mask for the nodes.  Only 
!    those nodes with nonzero input mask values are considered by the 
!    routine.  The nodes numbered by RCM will have their mask values 
!    set to zero.
!
!    Output, integer PERM(NODE_NUM), the RCM ordering.
!
!    Output, integer ICCSZE, the size of the connected component
!    that has been numbered.
!
!    Input, integer NODE_NUM, the number of nodes.
!
!  Local parameters:
!
!    Workspace, integer DEG(NODE_NUM), a temporary vector used to hold 
!    the degree of the nodes in the section graph specified by mask and root.
!
  implicit none

  integer root
  integer adj_num
  integer node_num

  integer adj(adj_num)
  integer adj_row(node_num+1)
  integer mask(node_num)
  integer perm(node_num)


  integer, dimension(:), allocatable :: deg
  integer fnbr
  integer i
  integer iccsze
  integer j
  integer jstop
  integer jstrt
  integer k
  integer l
  integer lbegin
  integer lnbr
  integer lperm
  integer lvlend
  integer nbr
  integer node
  integer err

  allocate(deg(node_num),stat=err)
  if (err.ne.0) STOP 'ERROR allocation in rcm'
!
!  Find the degrees of the nodes in the component specified by MASK and ROOT.
!
  call degree(root, adj_num, adj_row, adj, mask, deg, iccsze, perm, node_num)

  mask(root) = 0

  if ( iccsze <= 1 ) then
    return
  end if

  lvlend = 0
  lnbr = 1
!
!  LBEGIN and LVLEND point to the beginning and
!  the end of the current level respectively.
!
  do while ( lvlend < lnbr )

    lbegin = lvlend + 1
    lvlend = lnbr

    do i = lbegin, lvlend
!
!  For each node in the current level...
!
      node = perm(i)
      jstrt = adj_row(node)
      jstop = adj_row(node+1) - 1
!
!  Find the unnumbered neighbors of NODE.
!
!  FNBR and LNBR point to the first and last neighbors
!  of the current node in PERM.
!
      fnbr = lnbr + 1

      do j = jstrt, jstop

        nbr = adj(j)

        if ( mask(nbr) /= 0 ) then
          lnbr = lnbr + 1
          mask(nbr) = 0
          perm(lnbr) = nbr
        end if

      end do
!
!  If no neighbors, skip to next node in this level.
!
      if ( lnbr <= fnbr ) then
        cycle
      end if
!
!  Sort the neighbors of NODE in increasing order by degree.
!  Linear insertion is used.
!
      k = fnbr

      do while ( k < lnbr )

        l = k
        k = k + 1
        nbr = perm(k)

        do while ( fnbr < l )

          lperm = perm(l)

          if ( deg(lperm) <= deg(nbr) ) then
            exit
          end if

          perm(l+1) = lperm
          l = l-1

        end do

        perm(l+1) = nbr

      end do

    end do

  end do

  deallocate(deg)
!
!  We now have the Cuthill-McKee ordering.  Reverse it.
!
  call ivec_reverse( iccsze, perm )

end subroutine rcm

subroutine root_find ( root, adj_num, adj_row, adj, mask, level_num, &
  level_row, level, node_num )

!*******************************************************************************
!
!! ROOT_FIND finds a pseudo-peripheral node.
!
!  Discussion:
!
!    The diameter of a graph is the maximum distance (number of edges)
!    between any two nodes of the graph.
!
!    The eccentricity of a node is the maximum distance between that
!    node and any other node of the graph.
!
!    A peripheral node is a node whose eccentricity equals the
!    diameter of the graph.
!
!    A pseudo-peripheral node is an approximation to a peripheral node;
!    it may be a peripheral node, but all we know is that we tried our
!    best.
!
!    The routine is given a graph, and seeks pseudo-peripheral nodes,
!    using a modified version of the scheme of Gibbs, Poole and
!    Stockmeyer.  It determines such a node for the section subgraph
!    specified by MASK and ROOT.
!
!    The routine also determines the level structure associated with
!    the given pseudo-peripheral node; that is, how far each node
!    is from the pseudo-peripheral node.  The level structure is
!    returned as a list of nodes LS, and pointers to the beginning
!    of the list of nodes that are at a distance of 0, 1, 2, ...,
!    NODE_NUM-1 from the pseudo-peripheral node.
!
!  Modified:
!
!    28 October 2003
!
!  Reference:
!
!    Alan George and J W Liu,
!    Computer Solution of Large Sparse Positive Definite Systems,
!    Prentice Hall, 1981.
!
!    Gibbs, Poole, Stockmeyer,
!    An Algorithm for Reducing the Bandwidth and Profile of a Sparse Matrix,
!    SIAM Journal on Numerical Analysis,
!    Volume 13, pages 236-250, 1976.
!
!    Gibbs,
!    Algorithm 509: A Hybrid Profile Reduction Algorithm,
!    ACM Transactions on Mathematical Software,
!    Volume 2, pages 378-387, 1976.
!
!  Parameters:
!
!    Input/output, integer ROOT.  On input, ROOT is a node in the
!    the component of the graph for which a pseudo-peripheral node is
!    sought.  On output, ROOT is the pseudo-peripheral node obtained.
!
!    Input, integer ADJ_NUM, the number of adjacency entries.
!
!    Input, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Input, integer ADJ(ADJ_NUM), the adjacency structure.
!    For each row, it contains the column indices of the nonzero entries.
!
!    Input, integer MASK(NODE_NUM), specifies a section subgraph.  Nodes 
!    for which MASK is zero are ignored by FNROOT.
!
!    Output, integer LEVEL_NUM, is the number of levels in the level structure
!    rooted at the node ROOT.
!
!    Output, integer LEVEL_ROW(NODE_NUM+1), integer LEVEL(NODE_NUM), the 
!    level structure array pair containing the level structure found.
!
!    Input, integer NODE_NUM, the number of nodes.
!
  implicit none

  integer adj_num
  integer node_num
  integer level_num
  integer adj(adj_num)
  integer adj_row(node_num+1)
  integer mask(node_num)
  integer level(node_num)
  integer level_row(node_num+1)
  ! locals
  integer iccsze
  integer j
  integer jstrt
  integer k
  integer kstop
  integer kstrt
  integer level_num2
  integer mindeg
  integer nabor
  integer ndeg
  integer node
  integer root
!
!  Determine the level structure rooted at ROOT.
!
  call level_set ( root, adj_num, adj_row, adj, mask, level_num, &
    level_row, level, node_num )
!
!  Count the number of nodes in this level structure.
!
  iccsze = level_row(level_num+1) - 1


!
!  Extreme case:
!    A complete graph has a level set of only a single level.
!    Every node is equally good (or bad).
!
  if ( level_num == 1 ) then
    return
  end if
!
!  Extreme case:
!    A "line graph" 0--0--0--0--0 has every node in its only level.
!    By chance, we've stumbled on the ideal root.
!
  if ( level_num == iccsze ) then
    return
  end if
!
!  Pick any node from the last level that has minimum degree
!  as the starting point to generate a new level set.
!
  do

    mindeg = iccsze

    jstrt = level_row(level_num)
    root = level(jstrt)

    if ( jstrt < iccsze ) then

      do j = jstrt, iccsze

        node = level(j)
        ndeg = 0

        kstrt = adj_row(node)
        kstop = adj_row(node+1)-1

        do k = kstrt, kstop
          nabor = adj(k)
          if ( 0 < mask(nabor) ) then
            ndeg = ndeg+1
          end if
        end do

        if ( ndeg < mindeg ) then
          root = node
          mindeg = ndeg
        end if

      end do
 
    end if
!
!  Generate the rooted level structure associated with this node.
!
    call level_set ( root, adj_num, adj_row, adj, mask, level_num2, &
      level_row, level, node_num )

!
!  If the number of levels did not increase, accept the new ROOT.
!
    if ( level_num2 <= level_num ) then
      exit
    end if

    level_num = level_num2
!
!  In the unlikely case that ROOT is one endpoint of a line graph,
!  we can exit now.
!
    if ( iccsze <= level_num ) then
      exit
    end if

  end do

  return
end subroutine root_find

subroutine sort_heap_external ( n, indx, i, j, isgn )

!*******************************************************************************
!
!! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
!
!  Discussion:
!
!    The actual list of data is not passed to the routine.  Hence this
!    routine may be used to sort integers, reals, numbers, names,
!    dates, shoe sizes, and so on.  After each call, the routine asks
!    the user to compare or interchange two items, until a special
!    return value signals that the sorting is completed.
!
!  Modified:
!
!    05 February 2004
!
!  Reference:
!
!    A Nijenhuis and H Wilf,
!    Combinatorial Algorithms,
!    Academic Press, 1978, second edition,
!    ISBN 0-12-519260-6.
!
!  Parameters:
!
!    Input, integer N, the number of items to be sorted.
!
!    Input/output, integer INDX, the main communication signal.
!
!    The user must set INDX to 0 before the first call.
!    Thereafter, the user should not change the value of INDX until
!    the sorting is done.
!
!    On return, if INDX is
!
!      greater than 0,
!      * interchange items I and J;
!      * call again.
!
!      less than 0,
!      * compare items I and J;
!      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
!      * call again.
!
!      equal to 0, the sorting is done.
!
!    Output, integer I, J, the indices of two items.
!    On return with INDX positive, elements I and J should be interchanged.
!    On return with INDX negative, elements I and J should be compared, and
!    the result reported in ISGN on the next call.
!
!    Input, integer ISGN, results of comparison of elements I and J.
!    (Used only when the previous call returned INDX less than 0).
!    ISGN <= 0 means I is less than or equal to J;
!    0 <= ISGN means I is greater than or equal to J.
!
  implicit none

  integer i
  integer, save :: i_save = 0
  integer indx
  integer isgn
  integer j
  integer, save :: j_save = 0
  integer, save :: k = 0
  integer, save :: k1 = 0
  integer n
  integer, save :: n1 = 0
!
!  INDX = 0: This is the first call.
!
  if ( indx == 0 ) then

    i_save = 0
    j_save = 0
    k = n / 2
    k1 = k
    n1 = n
!
!  INDX < 0: The user is returning the results of a comparison.
!
  else if ( indx < 0 ) then

    if ( indx == -2 ) then

      if ( isgn < 0 ) then
        i_save = i_save + 1
      end if

      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return

    end if

    if ( 0 < isgn ) then
      indx = 2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then

      if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
      end if

      i = i_save
      j = j_save
      return

    end if

    k = k - 1
    k1 = k
!
!  0 < INDX, the user was asked to make an interchange.
!
  else if ( indx == 1 ) then

    k1 = k

  end if

  do

    i_save = 2 * k1

    if ( i_save == n1 ) then
      j_save = k1
      k1 = i_save
      indx = -1
      i = i_save
      j = j_save
      return
    else if ( i_save <= n1 ) then
      j_save = i_save + 1
      indx = -2
      i = i_save
      j = j_save
      return
    end if

    if ( k <= 1 ) then
      exit
    end if

    k = k - 1
    k1 = k

  end do

  if ( n1 == 1 ) then
    i_save = 0
    j_save = 0
    indx = 0
    i = i_save
    j = j_save
  else
    i_save = n1
    n1 = n1 - 1
    j_save = 1
    indx = 1
    i = i_save
    j = j_save
  end if

  return
end subroutine sort_heap_external


subroutine timestamp ( )

!*******************************************************************************
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    15 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 40 ) string

  call timestring ( string )

  write ( *, '(a)' ) trim ( string )

  return
end subroutine timestamp

subroutine timestring ( string )

!*******************************************************************************
!
!! TIMESTRING writes the current YMDHMS date into a string.
!
!  Example:
!
!    STRING = 'May 31 2001   9:45:54.872 AM'
!
!  Modified:
!
!    15 March 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) STRING, contains the date information.
!    A character length of 40 should always be sufficient.
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = * ) string
  character ( len = 10 ) time
  integer values(8)
  integer y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( string, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end subroutine timestring

subroutine triangulation_order3_adj_count ( node_num, tri_num, triangle_node, &
  triangle_neighbor, adj_num, adj_row )

!*******************************************************************************
!
!! TRIANGULATION_ORDER3_ADJ_COUNT counts adjacencies in a triangulation.
!
!  Discussion:
!
!    This routine is called to count the adjacencies, so that the
!    appropriate amount of memory can be set aside for storage when
!    the adjacency structure is created.
!
!    The adjacency relation is assumed to be symmetric.  Only the
!    strict lower half of the adjacency matrix is to be considered.
!
!    The triangulation is assumed to involve 3-node triangles.
!
!    Two nodes are "adjacent" if they are both nodes in some triangle.
!
!  Diagram:
!
!       3
!    s  |\
!    i  | \
!    d  |  \
!    e  |   \  side 2
!       |    \
!    3  |     \
!       |      \
!       1-------2
!
!         side 1
!
!    The local node numbering
!
!
!   21-22-23-24-25
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   16-17-18-19-20
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   11-12-13-14-15
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    6--7--8--9-10
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    1--2--3--4--5
!
!    A sample grid, with the following adjacencies:
!
!    1:  *  2 6
!    2:  1  *  3  6  7
!    3:  2  *  4  7  8
!    4:  3  *  5  8  9
!    5:  4  *  9 10
!    6:  1  2  *  7 11
!    7:  2  3  6  *  8 11 12
!    8:  3  4  7  *  9 12 13
!    9:  4  5  8  * 10 13 14
!   10:  5  9  * 14 15
!   11:  6  7  * 12 16
!   12:  7  8 11  * 13 16 17
!   13:  8  9 12  * 14 17 18
!   14:  9 10 13  * 15 18 19
!   15: 10 14  * 19 20
!   16: 11 12  * 17 21
!   17: 12 13 16  * 18 21 22
!   18: 13 14 17  * 19 22 23
!   19: 14 15 18  * 20 23 24
!   20: 15 19  * 24 25
!   21: 16 17  * 22
!   22: 17 18 21  * 23
!   23: 18 19 22  * 24
!   24: 19 20 23  * 25
!   25: 20 24  *
!
!    Counting the adjacencies, we get:
!
!      2, 4, 4, 4, 3,
!      4, 6, 6, 6, 4,
!      4, 6, 6, 6, 4,
!      4, 6, 6, 6, 4,
!      3, 4, 4, 4, 2.
!
!    yielding ADJ_ROW:
!
!      1,   3,   7,  11,  15,
!     18,  22,  28,  34,  40,
!     44,  48,  54,  60,  66,
!     70,  74,  80,  86,  92,
!     96,  99, 103, 107, 111,
!    113
!      
!  Modified:
!
!    11 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Input, integer TRIANGLE_NODE(3,TRI_NUM), lists the nodes that
!    make up each triangle, in counterclockwise order. 
!
!    Input, integer TRIANGLE_NEIGHBOR(3,TRI_NUM), for each side of
!    a triangle, lists the neighboring triangle, or -1 if there is
!    no neighbor.
!
!    Output, integer ADJ_NUM, the number of adjacencies.
!
!    Output, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
  implicit none

  integer node_num
  integer tri_num

  integer adj_num
  integer adj_row(node_num+1)
  integer i
  integer n1
  integer n2
  integer n3
  integer triangle
  integer triangle2
  integer triangle_neighbor(3,tri_num)
  integer triangle_node(3,tri_num)

  adj_num = 0
  adj_row(1:node_num) = 0
!
!  Examine each triangle.
!
  do triangle = 1, tri_num

    n1 = triangle_node(1,triangle)
    n2 = triangle_node(2,triangle)
    n3 = triangle_node(3,triangle)
!
!  Add edge (1,2) if this is the first occurrence,
!  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
!  or if this triangle is the first of the pair in which the edge
!  occurs (TRIANGLE < TRIANGLE2).
!
    triangle2 = triangle_neighbor(1,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_row(n1) = adj_row(n1) + 1
      adj_row(n2) = adj_row(n2) + 1
    end if
!
!  Add edge (2,3).
!
    triangle2 = triangle_neighbor(2,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_row(n2) = adj_row(n2) + 1
      adj_row(n3) = adj_row(n3) + 1
    end if
!
!  Add edge (3,1).
!
    triangle2 = triangle_neighbor(3,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_row(n1) = adj_row(n1) + 1
      adj_row(n3) = adj_row(n3) + 1
    end if
      
  end do
!
!  We used ADJ_ROW to count the number of entries in each row.
!  Convert it to pointers into the ADJ array.
!
  adj_row(2:node_num+1) = adj_row(1:node_num)

  adj_row(1) = 1
  do i = 2, node_num+1
    adj_row(i) = adj_row(i-1) + adj_row(i)
  end do

  adj_num = adj_row(node_num+1) - 1

  return
end subroutine triangulation_order3_adj_count

subroutine triangulation_order3_adj_set ( node_num, tri_num, triangle_node, &
  triangle_neighbor, adj_num, adj_row, adj )

!*******************************************************************************
!
!! TRIANGULATION_ORDER3_ADJ_SET sets adjacencies in a triangulation.
!
!  Discussion:
!
!    This routine is called to count the adjacencies, so that the
!    appropriate amount of memory can be set aside for storage when
!    the adjacency structure is created.
!
!    The adjacency relation is assumed to be symmetric.  Only the
!    strict lower half of the adjacency matrix is to be considered.
!
!    The triangulation is assumed to involve 3-node triangles.
!
!    Two nodes are "adjacent" if they are both nodes in some triangle.
!
!  Diagram:
!
!       3
!    s  |\
!    i  | \
!    d  |  \
!    e  |   \  side 2
!       |    \
!    3  |     \
!       |      \
!       1-------2
!
!         side 1
!
!    The local node numbering
!
!
!   21-22-23-24-25
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   16-17-18-19-20
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   11-12-13-14-15
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    6--7--8--9-10
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    1--2--3--4--5
!
!    A sample grid, with the following adjacencies:
!
!    1:  *  2 6
!    2:  1  *  3  6
!    3:  2  *  4  7  8
!    4:  3  *  5  8  9
!    5:  4  *  9 10
!    6:  1  2  *  7 11
!    7:  2  3  6  *  8 11 12
!    8:  3  4  7  *  9 12 13
!    9:  4  5  8  * 10 13 14
!   10:  5  9  * 14 15
!   11:  6  7  * 12 16
!   12:  7  8 11  * 13 16 17
!   13:  8  9 12  * 14 17 18
!   14:  9 10 13  * 15 18 19
!   15: 10 14  * 19 20
!   16: 11 12  * 17 21
!   17: 12 13 16  * 18 21 22
!   18: 13 14 17  * 19 22 23
!   19: 14 15 18  * 20 23 24
!   20: 15 19  * 24 25
!   21: 16 17  * 22
!   22: 17 18 21  * 23
!   23: 18 19 22  * 24
!   24: 19 20 23  * 25
!   25: 20 24  *
!
!    Counting the adjacencies, we get:
!
!      2, 4, 4, 4, 3,
!      4, 6, 6, 6, 4,
!      4, 6, 6, 6, 4,
!      4, 6, 6, 6, 4,
!      3, 4, 4, 4, 2.
!
!    yielding ADJ_ROW:
!
!      1,   3,   7,  11,  15,
!     18,  22,  28,  34,  40,
!     44,  48,  54,  60,  66,
!     70,  74,  80,  86,  92,
!     96,  99, 103, 107, 111,
!    113
!         
!  Modified:
!
!    11 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Input, integer TRIANGLE_NODE(3,TRI_NUM), lists the nodes that
!    make up each triangle in counterclockwise order.
!
!    Input, integer TRIANGLE_NEIGHBOR(3,TRI_NUM), for each side of
!    a triangle, lists the neighboring triangle, or -1 if there is
!    no neighbor.
!
!    Input, integer ADJ_NUM, the number of adjacencies.
!
!    Output, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Output, integer ADJ(ADJ_NUM), the adjacency information.
!
  implicit none

  integer adj_num
  integer node_num
  integer tri_num

  integer adj(adj_num)
  integer adj_copy(node_num)
  integer adj_row(node_num+1)
  integer n1
  integer n2
  integer n3
  integer triangle
  integer triangle2
  integer triangle_neighbor(3,tri_num)
  integer triangle_node(3,tri_num)

  adj(1:adj_num) = -1
  adj_copy(1:node_num) = adj_row(1:node_num)
!
!  Examine each triangle.
!
  do triangle = 1, tri_num

    n1 = triangle_node(1,triangle)
    n2 = triangle_node(2,triangle)
    n3 = triangle_node(3,triangle)
!
!  Add edge (1,2) if this is the first occurrence,
!  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
!  or if this triangle is the first of the pair in which the edge
!  occurs (TRIANGLE < TRIANGLE2).
!
    triangle2 = triangle_neighbor(1,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj(adj_copy(n1)) = n2
      adj_copy(n1) = adj_copy(n1) + 1
      adj(adj_copy(n2)) = n1
      adj_copy(n2) = adj_copy(n2) + 1
    end if
!
!  Add edge (2,3).
!
    triangle2 = triangle_neighbor(2,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj(adj_copy(n2)) = n3
      adj_copy(n2) = adj_copy(n2) + 1
      adj(adj_copy(n3)) = n2
      adj_copy(n3) = adj_copy(n3) + 1
    end if
!
!  Add edge (3,1).
!
    triangle2 = triangle_neighbor(3,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj(adj_copy(n1)) = n3
      adj_copy(n1) = adj_copy(n1) + 1
      adj(adj_copy(n3)) = n1
      adj_copy(n3) = adj_copy(n3) + 1
    end if
      
  end do

  return
end subroutine triangulation_order3_adj_set

subroutine triangulation_order3_example ( node_num, tri_num, triangle_node, &
  triangle_neighbor )

!*******************************************************************************
!
!! TRIANGULATION_ORDER3_EXAMPLE returns an example triangulation.
!
!  Diagram:
!
!
!   21-22-23-24-25
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   16-17-18-19-20
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   11-12-13-14-15
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    6--7--8--9-10
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    1--2--3--4--5
!
!  Modified:
!
!    12 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Output, integer TRIANGLE_NODE(3,TRI_NUM), lists the nodes that
!    make up each triangle, in counterclockwise order.
!
!    Output, integer TRIANGLE_NEIGHBOR(3,TRI_NUM), for each side of
!    a triangle, lists the neighboring triangle, or -1 if there is
!    no neighbor.
!
  implicit none

  integer node_num
  integer tri_num

  integer triangle_neighbor(3,tri_num)
  integer triangle_node(3,tri_num)

  triangle_node(1:3,1:tri_num) = reshape ( (/ &
     1,  2,  6, &
     7,  6,  2, &
     2,  3,  7, &
     8,  7,  3, &
     3,  4,  8, &
     9,  8,  4, &
     4,  5,  9, &
    10,  9,  5, &
     6,  7, 11, &
    12, 11,  7, &
     7,  8, 12, &
    13, 12,  8, &
     8,  9, 13, &
    14, 13,  9, &
     9, 10, 14, &
    15, 14, 10, &
    11, 12, 16, &
    17, 16, 12, &
    12, 13, 17, &
    18, 17, 13, &
    13, 14, 18, &
    19, 18, 14, &
    14, 15, 19, &
    20, 19, 15, &
    16, 17, 21, &
    22, 21, 17, &
    17, 18, 22, &
    23, 22, 18, &
    18, 19, 23, &
    24, 23, 19, &
    19, 20, 24, &
    25, 24, 20 /), (/ 3, tri_num /) )

  triangle_neighbor(1:3,1:tri_num) = reshape ( (/ &
    -1,  2, -1, &
     9,  1,  3, &
    -1,  4,  2, &
    11,  3,  5, &
    -1,  6,  4, &
    13,  5,  7, &
    -1,  8,  6, &
    15,  7, -1, & 
     2, 10, -1, &
    17,  9, 11, &
     4, 12, 10, &
    19, 11, 13, &
     6, 14, 12, &
    21, 13, 15, &
     8, 16, 14, &
    23, 15, -1, & 
    10, 18, -1, &
    25, 17, 19, &
    12, 20, 18, &
    27, 19, 21, &
    14, 22, 20, &
    29, 21, 23, &
    16, 24, 22, &
    31, 23, -1, & 
    18, 26, -1, &
    -1, 25, 27, &
    20, 28, 26, &
    -1, 27, 29, &
    22, 30, 28, &
    -1, 29, 31, &
    24, 32, 30, &
    -1, 31, -1 /), (/ 3, tri_num /) )

  return
end subroutine triangulation_order3_example

subroutine triangulation_order3_example_size ( node_num, tri_num )

!*******************************************************************************
!
!! TRIANGULATION_ORDER3_EXAMPLE_SIZE returns the size of an example.
!
!  Diagram:
!
!
!   21-22-23-24-25
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   16-17-18-19-20
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!   11-12-13-14-15
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    6--7--8--9-10
!    |\ |\ |\ |\ |
!    | \| \| \| \|
!    1--2--3--4--5
!
!  Modified:
!
!    12 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer TRI_NUM, the number of triangles.
!
  implicit none

  integer node_num
  integer tri_num

  node_num = 25
  tri_num = 32

  return
end subroutine triangulation_order3_example_size

subroutine triangulation_order3_neighbor_triangles ( tri_num, tri_vert, &
  tri_nabe )

!*******************************************************************************
!
!! TRIANGULATION_ORDER3_NEIGHBOR_TRIANGLES determines triangle neighbors.
!
!  Discussion:
!
!    A triangulation of a set of nodes can be completely described by
!    the coordinates of the nodes, and the list of nodes that make up
!    each triangle.  However, in some cases, it is necessary to know
!    triangle adjacency information, that is, which triangle, if any,
!    is adjacent to a given triangle on a particular side.
!
!    This routine creates a data structure recording this information.
!
!    The primary amount of work occurs in sorting a list of 3 * TRI_NUM
!    data items.
!
!  Example:
!
!    The input information from TRI_VERT:
!
!    Triangle   Nodes
!    --------   ---------------
!     1         3      4      1
!     2         3      1      2
!     3         3      2      8
!     4         2      1      5
!     5         8      2     13
!     6         8     13      9
!     7         3      8      9
!     8        13      2      5
!     9         9     13      7
!    10         7     13      5
!    11         6      7      5
!    12         9      7      6
!    13        10      9      6
!    14         6      5     12
!    15        11      6     12
!    16        10      6     11
!
!    The output information in TRI_NABE:
!
!    Triangle  Neighboring Triangles
!    --------  ---------------------
!
!     1        -1     -1      2
!     2         1      4      3
!     3         2      5      7
!     4         2     -1      8
!     5         3      8      6
!     6         5      9      7
!     7         3      6     -1
!     8         5      4     10
!     9         6     10     12
!    10         9      8     11
!    11        12     10     14
!    12         9     11     13
!    13        -1     12     16
!    14        11     -1     15
!    15        16     14     -1
!    16        13     15     -1
!
!  Modified:
!
!    17 July 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Input, integer TRI_VERT(3,TRI_NUM), the nodes that make up each triangle.
!
!    Output, integer TRI_NABE(3,TRI_NUM), the three triangles that are direct
!    neighbors of a given triangle.  TRI_NABE(1,I) is the index of the triangle
!    which touches side 1, defined by nodes 2 and 3, and so on.  TRI_NABE(1,I)
!    is negative if there is no neighbor on that side.  In this case, that
!    side of the triangle lies on the boundary of the triangulation.
!
  implicit none

  integer tri_num

  integer i
  integer irow
  integer j
  integer k
  integer row(3*tri_num,4)
  integer side1
  integer side2
  integer tri_nabe(3,tri_num)
  integer tri
  integer tri_vert(3,tri_num)
  integer tri1
  integer tri2
!
!  Step 1.
!  From the list of nodes for triangle T, of the form: (I,J,K)
!  construct the three neighbor relations:
!
!    (I,J,1,T) or (J,I,1,T),
!    (J,K,2,T) or (K,J,2,T),
!    (K,I,3,T) or (I,K,3,T)
!
!  where we choose (I,J,1,T) if I < J, or else (J,I,1,T)
!
  do tri = 1, tri_num

    i = tri_vert(1,tri)
    j = tri_vert(2,tri)
    k = tri_vert(3,tri)

    if ( i < j ) then
      row(3*(tri-1)+1,1:4) = (/ i, j, 1, tri /)
    else
      row(3*(tri-1)+1,1:4) = (/ j, i, 1, tri /)
    end if

    if ( j < k ) then
      row(3*(tri-1)+2,1:4) = (/ j, k, 2, tri /)
    else
      row(3*(tri-1)+2,1:4) = (/ k, j, 2, tri /)
    end if

    if ( k < i ) then
      row(3*(tri-1)+3,1:4) = (/ k, i, 3, tri /)
    else
      row(3*(tri-1)+3,1:4) = (/ i, k, 3, tri /)
    end if

  end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!  We only intend to sort on columns 1 and 2; the routine we call here
!  sorts on columns 1 through 4 but that won't hurt us.
!
!  What we need is to find cases where two triangles share an edge.
!  Say they share an edge defined by the nodes I and J.  Then there are
!  two rows of ROW that start out ( I, J, ?, ? ).  By sorting ROW,
!  we make sure that these two rows occur consecutively.  That will
!  make it easy to notice that the triangles are neighbors.
!
  call irow_sort_a ( 3*tri_num, 4, row )
!
!  Step 3. Neighboring triangles show up as consecutive rows with
!  identical first two entries.  Whenever you spot this happening,
!  make the appropriate entries in TRI_NABE.
!
  tri_nabe(1:3,1:tri_num) = -1

  irow = 1

  do

    if ( 3 * tri_num <= irow ) then
      exit
    end if

    if ( row(irow,1) /= row(irow+1,1) .or. row(irow,2) /= row(irow+1,2) ) then
      irow = irow + 1
      cycle
    end if

    side1 = row(irow,3)
    tri1 = row(irow,4)
    side2 = row(irow+1,3)
    tri2 = row(irow+1,4)

    tri_nabe(side1,tri1) = tri2
    tri_nabe(side2,tri2) = tri1

    irow = irow + 2

  end do

  return
end subroutine triangulation_order3_neighbor_triangles

subroutine triangulation_order6_adj_count ( node_num, tri_num, triangle_node, &
  triangle_neighbor, adj_num, adj_row )

!*******************************************************************************
!
!! TRIANGULATION_ORDER6_ADJ_COUNT counts adjacencies in a triangulation.
!
!  Discussion:
!
!    This routine is called to count the adjacencies, so that the
!    appropriate amount of memory can be set aside for storage when
!    the adjacency structure is created.
!
!    The adjacency relation is assumed to be symmetric.  Only the
!    strict lower half of the adjacency matrix is to be considered.
!
!    The triangulation is assumed to involve 6-node triangles.
!
!    Two nodes are "adjacent" if they are both nodes in some triangle.
!
!  Diagram:
!
!       3
!    s  |\
!    i  | \
!    d  |  \
!    e  6   5  side 2
!       |    \
!    3  |     \
!       |      \
!       1---4---2
!
!         side 1
!
!    The local node numbering
!
!
!   21-22-23-24-25
!    |\    |\    |
!    | \   | \   |
!   16 17 18 19 20
!    |   \ |   \ |
!    |    \|    \|
!   11-12-13-14-15
!    |\    |\    |
!    | \   | \   |
!    6  7  8  9 10
!    |   \ |   \ |
!    |    \|    \|
!    1--2--3--4--5
!
!    A sample grid, with the following adjacencies:
!
!    1:  *  2  3  6  7 11
!    2:  1  *  3  6  7 11
!    3:  1  2  *  4  5  6  7  8  9 11 12 13
!    4:  3  *  5  8  9 13
!    5:  3  4  *  8  9 10 13 14 15
!    6:  1  2  3  *  7 11
!    7:  1  2  3  6  *  8 11 12 13
!    8:  3  4  5  7  *  9 11 12 13
!    9:  3  4  5  8  * 10 13 14 15
!   10:  5  9  * 13 14 15
!   11:  1  2  3  6  7  8  * 12 13 16 17 21
!   12:  3  7  8 11  * 13 16 17 21
!   13:  3  4  5  7  8  9 10 11 12  * 14 15 16 17 18 19 21 22 23
!   14:  5  9 10 13  * 15 18 19 23
!   15:  5  9 10 13 14  * 18 19 20 23 24 25
!   16: 11 12 13  * 17 21
!   17: 11 12 13 16  * 18 21 22 23
!   18: 13 14 15 17  * 19 21 22 23
!   19: 13 14 15 18  * 20 23 24 25
!   20: 15 19  * 23 24 25
!   21: 11 12 13 16 17 18  * 22 23
!   22: 13 17 18 21  * 23
!   23: 13 14 15 17 18 19 20 21 22  * 24 25
!   24: 15 19 20 23  * 25
!   25: 15 19 20 23 24  *
!
!    Counting only the strict lower triangle of adjacency, we get
!
!      0, 1, 2, 1, 2, 
!      3, 4, 4, 4, 2,
!      6, 4, 9, 4, 5,
!      3, 4, 4, 4, 2,
!      6, 4, 9, 4, 5
!
!    Yielding ADJ_ROW:
!
!      1,  1,  2,  4,  5,
!      7, 10, 14, 18, 22,
!     24, 30, 34, 43, 47,
!     52, 55, 59, 63, 67,
!     69, 75, 79, 88, 92,
!     97.
!      
!  Modified:
!
!    11 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Input, integer TRIANGLE_NODE(6,TRI_NUM), lists the nodes that
!    make up each triangle.  The first three nodes are the vertices,
!    in counterclockwise order.  The fourth value is the midside
!    node between nodes 1 and 2; the fifth and sixth values are
!    the other midside nodes in the logical order.
!
!    Input, integer TRIANGLE_NEIGHBOR(3,TRI_NUM), for each side of
!    a triangle, lists the neighboring triangle, or -1 if there is
!    no neighbor.
!
!    Output, integer ADJ_NUM, the number of adjacencies.
!
!    Output, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
  implicit none

  integer node_num
  integer tri_num

  integer adj_num
  integer adj_row(node_num+1)
  integer i
  integer n1
  integer n2
  integer n3
  integer n4
  integer n5
  integer n6
  integer triangle
  integer triangle2
  integer triangle_neighbor(3,tri_num)
  integer triangle_node(6,tri_num)

  adj_num = 0
  adj_row(1:node_num) = 0
!
!  Examine each triangle.
!
  do triangle = 1, tri_num

    n1 = triangle_node(1,triangle)
    n2 = triangle_node(2,triangle)
    n3 = triangle_node(3,triangle)
    n4 = triangle_node(4,triangle)
    n5 = triangle_node(5,triangle)
    n6 = triangle_node(6,triangle)
!
!  For sure, we add the adjacencies:
!    43 / (34)
!    51 / (15)
!    54 / (45)
!    62 / (26)
!    64 / (46)
!    65 / (56)
!
    adj_row(n3) = adj_row(n3) + 1
    adj_row(n4) = adj_row(n4) + 1
    adj_row(n1) = adj_row(n1) + 1
    adj_row(n5) = adj_row(n5) + 1
    adj_row(n4) = adj_row(n4) + 1
    adj_row(n5) = adj_row(n5) + 1
    adj_row(n2) = adj_row(n2) + 1
    adj_row(n6) = adj_row(n6) + 1
    adj_row(n4) = adj_row(n4) + 1
    adj_row(n6) = adj_row(n6) + 1
    adj_row(n5) = adj_row(n5) + 1
    adj_row(n6) = adj_row(n6) + 1
!
!  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
!  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
!  or if this triangle is the first of the pair in which the edge
!  occurs (TRIANGLE < TRIANGLE2).
!
!  Maybe add
!    21 / 12
!    41 / 14
!    42 / 24
!
    triangle2 = triangle_neighbor(1,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_row(n1) = adj_row(n1) + 1
      adj_row(n2) = adj_row(n2) + 1
      adj_row(n1) = adj_row(n1) + 1
      adj_row(n4) = adj_row(n4) + 1
      adj_row(n2) = adj_row(n2) + 1
      adj_row(n4) = adj_row(n4) + 1
    end if
!
!  Maybe add
!    32 / 23
!    52 / 25
!    53 / 35
!
    triangle2 = triangle_neighbor(2,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_row(n2) = adj_row(n2) + 1
      adj_row(n3) = adj_row(n3) + 1
      adj_row(n2) = adj_row(n2) + 1
      adj_row(n5) = adj_row(n5) + 1
      adj_row(n3) = adj_row(n3) + 1
      adj_row(n5) = adj_row(n5) + 1
    end if
!
!  Maybe add
!    31 / 13
!    61 / 16
!    63 / 36
!
    triangle2 = triangle_neighbor(3,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj_row(n1) = adj_row(n1) + 1
      adj_row(n3) = adj_row(n3) + 1
      adj_row(n1) = adj_row(n1) + 1
      adj_row(n6) = adj_row(n6) + 1
      adj_row(n3) = adj_row(n3) + 1
      adj_row(n6) = adj_row(n6) + 1
    end if
      
  end do
!
!  We used ADJ_ROW to count the number of entries in each row.
!  Convert it to pointers into the ADJ array.
!
  adj_row(2:node_num+1) = adj_row(1:node_num)

  adj_row(1) = 1
  do i = 2, node_num+1
    adj_row(i) = adj_row(i-1) + adj_row(i)
  end do

  adj_num = adj_row(node_num+1) - 1

  return
end subroutine triangulation_order6_adj_count

subroutine triangulation_order6_adj_set ( node_num, tri_num, triangle_node, &
  triangle_neighbor, adj_num, adj_row, adj )

!*******************************************************************************
!
!! TRIANGULATION_ORDER6_ADJ_SET sets adjacencies in a triangulation.
!
!  Discussion:
!
!    This routine is called to count the adjacencies, so that the
!    appropriate amount of memory can be set aside for storage when
!    the adjacency structure is created.
!
!    The adjacency relation is assumed to be symmetric.  Only the
!    strict lower half of the adjacency matrix is to be considered.
!
!    The triangulation is assumed to involve 6-node triangles.
!
!    Two nodes are "adjacent" if they are both nodes in some triangle.
!
!  Diagram:
!
!       3
!    s  |\
!    i  | \
!    d  |  \
!    e  6   5  side 2
!       |    \
!    3  |     \
!       |      \
!       1---4---2
!
!         side 1
!
!    The local node numbering
!
!
!   21-22-23-24-25
!    |\    |\    |
!    | \   | \   |
!   16 17 18 19 20
!    |   \ |   \ |
!    |    \|    \|
!   11-12-13-14-15
!    |\    |\    |
!    | \   | \   |
!    6  7  8  9 10
!    |   \ |   \ |
!    |    \|    \|
!    1--2--3--4--5
!
!    A sample grid, with the following adjacencies:
!
!    1:  *  2  3  6  7 11
!    2:  1  *  3  6  7 11
!    3:  1  2  *  4  5  6  7  8  9 11 12 13
!    4:  3  *  5  8  9 13
!    5:  3  4  *  8  9 10 13 14 15
!    6:  1  2  3  *  7 11
!    7:  1  2  3  6  *  8 11 12 13
!    8:  3  4  5  7  *  9 11 12 13
!    9:  3  4  5  8  * 10 13 14 15
!   10:  5  9  * 13 14 15
!   11:  1  2  3  6  7  8  * 12 13 16 17 21
!   12:  3  7  8 11  * 13 16 17 21
!   13:  3  4  5  7  8  9 10 11 12  * 14 15 16 17 18 19 21 22 23
!   14:  5  9 10 13  * 15 18 19 23
!   15:  5  9 10 13 14  * 18 19 20 23 24 25
!   16: 11 12 13  * 17 21
!   17: 11 12 13 16  * 18 21 22 23
!   18: 13 14 15 17  * 19 21 22 23
!   19: 13 14 15 18  * 20 23 24 25
!   20: 15 19  * 23 24 25
!   21: 11 12 13 16 17 18  * 22 23
!   22: 13 17 18 21  * 23
!   23: 13 14 15 17 18 19 20 21 22  * 24 25
!   24: 15 19 20 23  * 25
!   25: 15 19 20 23 24  *
!
!    Counting only the strict lower triangle of adjacency, we get
!
!      0, 1, 2, 1, 2, 
!      3, 4, 4, 4, 2,
!      6, 4, 9, 4, 5,
!      3, 4, 4, 4, 2,
!      6, 4, 9, 4, 5
!
!    Yielding ADJ_ROW:
!
!      1,  1,  2,  4,  5,
!      7, 10, 14, 18, 22,
!     24, 30, 34, 43, 47,
!     52, 55, 59, 63, 67,
!     69, 75, 79, 88, 92,
!     97.
!      
!  Modified:
!
!    11 November 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Input, integer TRIANGLE_NODE(6,TRI_NUM), lists the nodes that
!    make up each triangle.  The first three nodes are the vertices,
!    in counterclockwise order.  The fourth value is the midside
!    node between nodes 1 and 2; the fifth and sixth values are
!    the other midside nodes in the logical order.
!
!    Input, integer TRIANGLE_NEIGHBOR(3,TRI_NUM), for each side of
!    a triangle, lists the neighboring triangle, or -1 if there is
!    no neighbor.
!
!    Input, integer ADJ_NUM, the number of adjacencies.
!
!    Output, integer ADJ_ROW(NODE_NUM+1).  Information about row I is stored
!    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
!
!    Output, integer ADJ(ADJ_NUM), the adjacency information.
!
  implicit none

  integer adj_num
  integer node_num
  integer tri_num

  integer adj(adj_num)
  integer adj_copy(node_num)
  integer adj_row(node_num+1)
  integer n1
  integer n2
  integer n3
  integer n4
  integer n5
  integer n6
  integer triangle
  integer triangle2
  integer triangle_neighbor(3,tri_num)
  integer triangle_node(6,tri_num)

  adj(1:adj_num) = -1
  adj_copy(1:node_num) = adj_row(1:node_num)
!
!  Examine each triangle.
!
  do triangle = 1, tri_num

    n1 = triangle_node(1,triangle)
    n2 = triangle_node(2,triangle)
    n3 = triangle_node(3,triangle)
    n4 = triangle_node(4,triangle)
    n5 = triangle_node(5,triangle)
    n6 = triangle_node(6,triangle)
!
!  For sure, we add the adjacencies:
!    43 / (34)
!    51 / (15)
!    54 / (45)
!    62 / (26)
!    64 / (46)
!    65 / (56)
!
    adj(adj_copy(n3)) = n4
    adj_copy(n3) = adj_copy(n3) + 1
    adj(adj_copy(n4)) = n3
    adj_copy(n4) = adj_copy(n4) + 1

    adj(adj_copy(n1)) = n5
    adj_copy(n1) = adj_copy(n1) + 1
    adj(adj_copy(n5)) = n1
    adj_copy(n5) = adj_copy(n5) + 1

    adj(adj_copy(n4)) = n5
    adj_copy(n4) = adj_copy(n4) + 1
    adj(adj_copy(n5)) = n4
    adj_copy(n5) = adj_copy(n5) + 1

    adj(adj_copy(n2)) = n6
    adj_copy(n2) = adj_copy(n2) + 1
    adj(adj_copy(n6)) = n2
    adj_copy(n6) = adj_copy(n6) + 1

    adj(adj_copy(n4)) = n6
    adj_copy(n4) = adj_copy(n4) + 1
    adj(adj_copy(n6)) = n4
    adj_copy(n6) = adj_copy(n6) + 1

    adj(adj_copy(n5)) = n6
    adj_copy(n5) = adj_copy(n5) + 1
    adj(adj_copy(n6)) = n5
    adj_copy(n6) = adj_copy(n6) + 1
!
!  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
!  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
!  or if this triangle is the first of the pair in which the edge
!  occurs (TRIANGLE < TRIANGLE2).
!
!  Maybe add
!    21 / 12
!    41 / 14
!    42 / 24
!
    triangle2 = triangle_neighbor(1,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj(adj_copy(n1)) = n2
      adj_copy(n1) = adj_copy(n1) + 1
      adj(adj_copy(n2)) = n1
      adj_copy(n2) = adj_copy(n2) + 1
      adj(adj_copy(n1)) = n4
      adj_copy(n1) = adj_copy(n1) + 1
      adj(adj_copy(n4)) = n1
      adj_copy(n4) = adj_copy(n4) + 1
      adj(adj_copy(n2)) = n4
      adj_copy(n2) = adj_copy(n2) + 1
      adj(adj_copy(n4)) = n2
      adj_copy(n4) = adj_copy(n4) + 1
    end if
!
!  Maybe add
!    32 / 23
!    52 / 25
!    53 / 35
!
    triangle2 = triangle_neighbor(2,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj(adj_copy(n2)) = n3
      adj_copy(n2) = adj_copy(n2) + 1
      adj(adj_copy(n3)) = n2
      adj_copy(n3) = adj_copy(n3) + 1
      adj(adj_copy(n2)) = n5
      adj_copy(n2) = adj_copy(n2) + 1
      adj(adj_copy(n5)) = n2
      adj_copy(n5) = adj_copy(n5) + 1
      adj(adj_copy(n3)) = n5
      adj_copy(n3) = adj_copy(n3) + 1
      adj(adj_copy(n5)) = n3
      adj_copy(n5) = adj_copy(n5) + 1
    end if
!
!  Maybe add
!    31 / 13
!    61 / 16
!    63 / 36
!
    triangle2 = triangle_neighbor(3,triangle)

    if ( triangle2 < 0 .or. triangle < triangle2 ) then
      adj(adj_copy(n1)) = n3
      adj_copy(n1) = adj_copy(n1) + 1
      adj(adj_copy(n3)) = n1
      adj_copy(n3) = adj_copy(n3) + 1
      adj(adj_copy(n1)) = n6
      adj_copy(n1) = adj_copy(n1) + 1
      adj(adj_copy(n6)) = n1
      adj_copy(n6) = adj_copy(n6) + 1
      adj(adj_copy(n3)) = n6
      adj_copy(n3) = adj_copy(n3) + 1
      adj(adj_copy(n6)) = n3
      adj_copy(n6) = adj_copy(n6) + 1
    end if
      
  end do

  return
end subroutine triangulation_order6_adj_set

subroutine triangulation_order6_example ( node_num, tri_num, triangle_node, &
  triangle_neighbor )

!*******************************************************************************
!
!! TRIANGULATION_ORDER6_EXAMPLE returns an example triangulation.
!
!  Diagram:
!
!
!   21-22-23-24-25
!    |\  6 |\  8 |
!    | \   | \   |
!   16 17 18 19 20
!    |   \ |   \ |
!    | 5  \| 7  \|
!   11-12-13-14-15
!    |\  2 |\  4 |
!    | \   | \   |
!    6  7  8  9 10
!    | 1 \ | 3 \ |
!    |    \|    \|
!    1--2--3--4--5
!
!  Modified:
!
!    12 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Output, integer TRIANGLE_NODE(6,TRI_NUM), lists the nodes that
!    make up each triangle.  The first three nodes are the vertices,
!    in counterclockwise order.  The fourth value is the midside
!    node between nodes 1 and 2; the fifth and sixth values are
!    the other midside nodes in the logical order.
!
!    Output, integer TRIANGLE_NEIGHBOR(3,TRI_NUM), for each side of
!    a triangle, lists the neighboring triangle, or -1 if there is
!    no neighbor.
!
  implicit none

  integer node_num
  integer tri_num

  integer triangle_neighbor(3,tri_num)
  integer triangle_node(6,tri_num)

  triangle_node(1:6,1:tri_num) = reshape ( (/ &
     1,  3, 11,  2,  7,  6, &
    13, 11,  3, 12,  7,  8, &
     3,  5, 13,  4,  9,  8, &
    15, 13,  5, 14,  9, 10, &
    11, 13, 21, 12, 17, 16, &
    23, 21, 13, 22, 17, 18, &
    13, 15, 23, 14, 19, 18, &
    25, 23, 15, 24, 19, 20 /), (/ 6, tri_num /) )

  triangle_neighbor(1:3,1:tri_num) = reshape ( (/ &
    -1,  2, -1, &
     5,  1,  3, &
    -1,  4,  2, &
     7,  3, -1, &
     2,  6, -1, &
    -1,  5,  7, &
     4,  8,  6, &
    -1,  7, -1 /), (/ 3, tri_num /) )

  return
end subroutine triangulation_order6_example

subroutine triangulation_order6_example_size ( node_num, tri_num )

!*******************************************************************************
!
!! TRIANGULATION_ORDER6_EXAMPLE_SIZE returns the size of an example.
!
!  Diagram:
!
!
!   21-22-23-24-25
!    |\  6 |\  8 |
!    | \   | \   |
!   16 17 18 19 20
!    |   \ |   \ |
!    | 5  \| 7  \|
!   11-12-13-14-15
!    |\  2 |\  4 |
!    | \   | \   |
!    6  7  8  9 10
!    | 1 \ | 3 \ |
!    |    \|    \|
!    1--2--3--4--5
!
!  Modified:
!
!    12 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, integer TRI_NUM, the number of triangles.
!
  implicit none

  integer node_num
  integer tri_num

  node_num = 25
  tri_num = 8

  return
end subroutine triangulation_order6_example_size



subroutine triangulation_order6_neighbor_triangles ( tri_num, tri_vert, &
  tri_nabe )

!*******************************************************************************
!
!! TRIANGULATION_ORDER6_NEIGHBOR_TRIANGLES determines triangle neighbors.
!
!  Discussion:
!
!    A triangulation of a set of nodes can be completely described by
!    the coordinates of the nodes, and the list of nodes that make up
!    each triangle.  However, in some cases, it is necessary to know
!    triangle adjacency information, that is, which triangle, if any,
!    is adjacent to a given triangle on a particular side.
!
!    This routine creates a data structure recording this information.
!
!    The primary amount of work occurs in sorting a list of 3 * TRI_NUM
!    data items.
!
!  Modified:
!
!    12 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer TRI_NUM, the number of triangles.
!
!    Input, integer TRI_VERT(6,TRI_NUM), the nodes that make up each triangle.
!
!    Output, integer TRI_NABE(3,TRI_NUM), the three triangles that are direct
!    neighbors of a given triangle.  TRI_NABE(1,I) is the index of the triangle
!    which touches side 1, defined by nodes 2 and 3, and so on.  TRI_NABE(1,I)
!    is negative if there is no neighbor on that side.  In this case, that
!    side of the triangle lies on the boundary of the triangulation.
!
  implicit none

  integer tri_num

  integer i
  integer irow
  integer j
  integer k
  integer row(3*tri_num,4)
  integer side1
  integer side2
  integer tri_nabe(3,tri_num)
  integer tri
  integer tri_vert(6,tri_num)
  integer tri1
  integer tri2
!
!  Step 1.
!  From the list of vertices for triangle T, of the form: (I,J,K)
!  construct the three neighbor relations:
!
!    (I,J,1,T) or (J,I,1,T),
!    (J,K,2,T) or (K,J,2,T),
!    (K,I,3,T) or (I,K,3,T)
!
!  where we choose (I,J,1,T) if I < J, or else (J,I,1,T)
!
  do tri = 1, tri_num

    i = tri_vert(1,tri)
    j = tri_vert(2,tri)
    k = tri_vert(3,tri)

    if ( i < j ) then
      row(3*(tri-1)+1,1:4) = (/ i, j, 1, tri /)
    else
      row(3*(tri-1)+1,1:4) = (/ j, i, 1, tri /)
    end if

    if ( j < k ) then
      row(3*(tri-1)+2,1:4) = (/ j, k, 2, tri /)
    else
      row(3*(tri-1)+2,1:4) = (/ k, j, 2, tri /)
    end if

    if ( k < i ) then
      row(3*(tri-1)+3,1:4) = (/ k, i, 3, tri /)
    else
      row(3*(tri-1)+3,1:4) = (/ i, k, 3, tri /)
    end if

  end do
!
!  Step 2. Perform an ascending dictionary sort on the neighbor relations.
!  We only intend to sort on columns 1 and 2; the routine we call here
!  sorts on columns 1 through 4 but that won't hurt us.
!
!  What we need is to find cases where two triangles share an edge.
!  Say they share an edge defined by the nodes I and J.  Then there are
!  two rows of ROW that start out ( I, J, ?, ? ).  By sorting ROW,
!  we make sure that these two rows occur consecutively.  That will
!  make it easy to notice that the triangles are neighbors.
!
  call irow_sort_a ( 3*tri_num, 4, row )
!
!  Step 3. Neighboring triangles show up as consecutive rows with
!  identical first two entries.  Whenever you spot this happening,
!  make the appropriate entries in TRI_NABE.
!
  tri_nabe(1:3,1:tri_num) = -1

  irow = 1

  do

    if ( 3 * tri_num <= irow ) then
      exit
    end if

    if ( row(irow,1) /= row(irow+1,1) .or. row(irow,2) /= row(irow+1,2) ) then
      irow = irow + 1
      cycle
    end if

    side1 = row(irow,3)
    tri1 = row(irow,4)
    side2 = row(irow+1,3)
    tri2 = row(irow+1,4)

    tri_nabe(side1,tri1) = tri2
    tri_nabe(side2,tri2) = tri1

    irow = irow + 2

  end do

  return
end subroutine triangulation_order6_neighbor_triangles

end module rcm_module
