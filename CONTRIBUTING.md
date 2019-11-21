Contributing to libNEGF
=======================

How to contribute
=================

The preferred method is to fork the project on and then create a
pull request. Your changes should be based on the master branch. Before you
start, please read the style convention below.

libNEGF has been contributed by several authors, and a lot of the code
is not homogeneous in style. Nevertheless, future development should follow
this guideline in order to increase the readability and homogeneity of the project.


Style guide
===============


Line length and indentation
---------------------------

* Maximal **line length** is **100** characters. For lines longer than that, use
  continuation lines.

* **Nested blocks** are indented by **2** white spaces::

      write(*, *) "Nested block follows"
      do ii = 1, 100
        write(*, *) "This is the nested block"
        if (ii == 50) then
        write(*, *) "Next nested block"
        end if
      end do

* **Continuation lines** are indented by **4** white spaces. Make sure to
  place continuation characters (`&`) both at the end of the line as well as at
  the beginning of the continuation line::

      call someRoutineWithManyParameters(param1, param2, param3, param4,&
          & param5)

  Try to break lines at natural places (e.g. at white space characters) and
  include one white space character after the opening ampersand in the
  continuation line.

* **Single line preprocessor directives** are indented as normal code::

      @:ASSERT(someCondition)
      call someRoutine(...)

* **Preprocessor block directives** (directives with starting and ending
  constructs) are outdented by **2** characters with respect of the code they
  enclose. The enclosed code must be aligned as if the preprocessor directives
  were not present::

      call doSomething()
      #:if WITH_SCALAPACK
      call someRoutineScalapackVersion(...)
      #:else
      call someRoutineSerialVersion(...)
      #:endif

      do iKS = 1, nKS
        #:if WITH_SCALAPACK
        call someRoutineScalapackVersion(iKS, ...)
        #:else
        call someRoutineSerialVersion(iKS, ...)
        #:endif
      end do



Naming
------

The naming conventions are similar to those in the [Google Style Guide for
C++ naming convention](https://google.github.io/styleguide/cppguide.html#Naming) with some modifications.

* Use **lowercase** for all Fortran constructs and prefer the `end do`, `end if` to `enddo` and `endif`::

      do i = 1, 2
        j = i + 1
      end do

* **Variable** names follow the **snake_case** convention::

      logical :: has_component

    Variable names should be descriptive, excessive use of abbreviations is discouraged.
    However full lowercase names without lowercase can be used when the name is short and
    clear. For example::

      integer :: numofprms ! Don't do this, multiple abbreviations without separator.
      ... :: kpoint        ! This is ok
      ... :: meshgrid      ! This is also ok
      ... :: save_to_mem   ! Ok, it is obvious that mem stands for memory

* **Constants** (parameters) can be declared using the **snake_case** convention similar to variables (especially for parameters with limited scope), or can be **UPPER_CASE** (especially for parameters defined at module level and externally accessible) ::

      integer, parameter :: max_array_size = 100
      integer, parameter :: MAX_ARRAY_SIZE = 100

  with the exception of the constants used to define the kind parameter for
  intrinsic types, which should be all lowercase (and short)::

      integer, parameter :: dp = kind(1.0d0)
      real(dp) :: val


* **Subroutine** and **function** names follow also the **snake_case**
  notation::

      subroutine test_some_functionality()
      myValue = get_some_value(...)


* **Type** (object) names are written **UpperCamelCase**::

      type :: RealList
      type(RealList) :: my_list

  Type names can be prefixed with a capital 'T' when they need to be clearly distinguished from other components, but it is not demanded necessary::

      type :: TNegf
      :
      end type TNegf
      :


* **Module** names follow **snake_case** convention::

      use dftb_common_accuracy

  Underscores are used for name-spacing only, so the module above would be
  typically found at the path `dftb/common/accuracy.f90`. The individual
  component names (``dftb``, ``common``, ``accuracy``) may not contain any
  underscores and must be shorter than 15 characters.


* **Preprocessor** variables and macros follow **UPPER_CASE_WITH_UNDERSCORE**
  convention::

    #:if WITH_MPI
      with_mpi = ${FORTRAN_LOGICAL(WITH_MPI)}$
    #:endif


White spaces
------------

Please use white spaces to make the code readable.

**Avoid trailing whitespaces**. They slow down code navigation and most modern editors trim trailing whitespaces by default. Please make sure that your editor does it too.
If a file comes with trailing whitespaces, feel free to clean it in a separate commit (see [above](#Refactoring-existing-code)).

In general, you **must use** white spaces in following situations:

* Around arithmetic operators::

      2 + 2

* Around assignment and pointer assignment operators::

      aa = 3 + 2
      window => array(1:3)

* Around the ``::`` separator in declarations::

      integer :: ind

* After commas (``,``) in general and especially in declarations, calls and
  lists::

      real(wp), allocatable :: array(:)
      type, extends(TBaseType) :: TDerivedType
      subroutine my_routine(par1, par2)
      call my_routine(val1, val2)
      print *, 'My value:', val
      do ii = 1, 3
      array(1:3) = [1, 2, 3]

* When separating array indices, when the actual index value for an index
  contains an expression::

      my_array(ii + 2, jj) = 12

You **may omit** white space in following cases:

* When separating array indices and the actual index values are simple and
  short (typically two letters) variable names, one or two digit integers or the
  range operator ``:``::

      my_array(:,1) = vector
      lat_vecs(1,1) = 1.0_wp
      my_array(ii,jj) = my_array(jj,ii)

You **must omit** white spaces in following cases:

* Around opening and closing braces of any kind::

      call my_subroutine(aa, bb)  ! and NOT call my_subroutine( aa, bb )
      my_vector(:) = [1, 2, 3]    ! instead of my_vector(:) = [ 1, 2, 3 ]
      tmp = 2 * (aa + bb)        ! instead of 2 * ( aa + bb )

* Around the equal (``=``) sign, when passing named arguments to a function or
  subroutine::

      call my_subroutine(aa, optional_argument=.true.)

* Around the power operator::

      val = base**power   (instead of val = base ** power)

**Avoid** white spaces for **visual aligning** of code, use::

      integer, intent(in) :: num_neighbors
      real(wp), intent(out) :: interaction

instead of::

      integer, intent(in)   :: num_eighbors
      real(wp), intent(out) :: energy

Although latter may look more readable, it makes rather difficult to track real
changes in the code with the revision control system. For example when a new
line is added to the block making the realignment of previous (but otherwise
unchanged) lines necessary ::

      integer, intent(in)             :: num_neighbors
      real(wp), intent(out)           :: energy
      real(wp), intent(out), optional :: forces(:)

the version control system will indicate all of those lines having been
modified, although only the alignment (but not the actual instructions) were
changed.


Comments
--------

* **Module**, **Subroutine** and **function** comments should be consistent with
  `doxygen <http://doxygen.org/>`_ / `FORD
  <https://github.com/cmacmackin/ford>`_ literate comments for publicly visible
  interfaces and variables.

* Comments are indented to the same position as the code they document::

      ! Take spin degeneracy into account
      energy = 2.0_wp * energy

* Generally, write the comment *before* the code snippet it documents::

      ! Loop over all neighbours
      do i_neigh = 1, num_neighbours
        :
      end do

* Try to avoid mixing code and comments within one line as this is often hard to
  read::

      bb = 2 * aa   ! this comment should be before the line.

* Never use multi-line suffix comments, as an indenting editor would mess up the
  indentation of subsequent lines::

      bb = 2 * aa  ! This comment goes over multiple lines, therefore, it
                   ! should stay ALWAYS before the code snippet and NOT HERE.

* Specifically comment any workarounds, include the compiler name and the
  version number for which the workaround had to be made. Always use the
  following pattern, so that searching for workarounds which can be possibly
  removed is easy::

      ! Workaround: gfortran 4.8
      ! Finalisation not working, we have to deallocate explicitly
      deallocate(myPointer)


* Comments should always start with one bang only. Comments with two bangs are
  reserved for source code documentation systems::

      ! This block needs a documentation
      do ii = 1, 2
        :
      end do

* If you need a comment for a longer block of code, consider instead packaging
  that block of code into a properly named function (if the additional function
  call would be performance critical, write it as an internal procedure)::

      some_previous_statement
      ind = get_first_non_zero(array)
      some_statement_after

  instead of ::

      some_previous_statement

      ! Look for the first nonzero element
      found = .false.
      do ind = 1, size(array)
        if (array(ind) > 0) then
	  found = .true.
	  exit
	    end if
      end do
      if (.not. found) then
        ind = 0
      end if

      some_statement_after

Allocation status
=================

At several places, the allocation status of a variable is used to signal choices
about logical flow in the code::

      !> SCC module internal variables
      type(TScc), allocatable :: scc_calc
      .
      .
      .
      if (allocated(scc_calc)) then

      end if

This is to be preferred to the use of additional logical variables if possible.

Part of the reason for this choice is that from Fortran 2008 onwards, optional
arguments to subroutines and functions are treated as not-present if not
allocated.


Refactoring existing code
=========================

Everyone is welcome to refactor existing code and namings, as long as this is harmless
from a user perspective (e.g. fixing indentation, whitespaces etc.). However, it is suggested to separate commits related to clean-up and refactoring from the commits
containing the actual features.

Impactful refactoring should be done in specific pull requests, to avoid compatibility issues.
