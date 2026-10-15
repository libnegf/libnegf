!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Greens functions.         !
!! Copyright (C) 2012 - 2026                                                !
!!                                                                          !
!! This file is part of libNEGF: a library for                              !
!! Non Equilibrium Green's Functions calculations                           !
!!                                                                          !
!! Developers: Alessandro Pecchia, Daniele Soccodato                        !
!! Former Contributors: Gabriele Penazzi, Luca Latessa, Aldo Di Carlo       !
!!                                                                          !
!! libNEGF is free software: you can redistribute and/or modify it          !
!! under the terms of the GNU Lesser General Public License as published    !
!! by the Free Software Foundation, either version 3 of the License, or     !
!! (at your option) any later version.                                      !
!!                                                                          !
!!  You should have received a copy of the GNU Lesser General Public        !
!!  License along with libNEGF.  If not, see                                !
!!  <http://www.gnu.org/licenses/>.                                         !
!!--------------------------------------------------------------------------!

!> Labelled timers accumulating wall time of the main calculation steps.
!>
!> Timers are identified by a fixed enumeration (timers%Gr, timers%SGF, ...), so that every
!> process holds the same list in the same order and results can be reduced over MPI.
!> Every process accumulates its own times. gather_timers() reduces min/avg/max over a
!> communicator and write_timers() prints the summary on the print node.
!>
!> Optionally set_timer/stop_timer also print a progress line on the print node, as the old
!> message_clock/write_clock, when the passed verbosity is larger than VBT.
!>
!> Note: timers are not thread-safe; call them outside OpenMP parallel regions.
module ln_timers

  use ln_precision, only : dp
  use iso_fortran_env, only : int64
  use mpi_globals, only : id0
#:if defined("MPI")
  use libmpifx_module, only : mpifx_comm, mpifx_reduce
  use mpi_f08, only : MPI_SUM, MPI_MIN, MPI_MAX
#:endif

  implicit none
  private

  public :: timers, NTIMERS, timer_names
  public :: set_timer, stop_timer, reset_timers
  public :: gather_timers, write_timers
  public :: get_timer_elapsed, get_timer_calls

  !> Enumeration of the timers
  type :: TTimerEnum
    integer :: ContactSE = 1
    integer :: SGF = 2
    integer :: Gr = 3
    integer :: Gn = 4
    integer :: ElasticSCBA = 5
    integer :: SigmaElph = 6
    integer :: SigmaElpt = 7
    integer :: Transmission = 8
    integer :: MeirWingreen = 9
    integer :: LayerCurrent = 10
    integer :: MPIGather = 11
  end type TTimerEnum

  type(TTimerEnum), parameter :: timers = TTimerEnum()

  integer, parameter :: NTIMERS = 11

  !> Labels used in the summary and as default progress message
  character(len=*), parameter :: timer_names(NTIMERS) = [character(len=24) :: &
      & 'Contact self-energy', &
      & 'Surface Green function', &
      & 'Green function Gr', &
      & 'Correlation function Gn', &
      & 'Elastic SCBA', &
      & 'El-phonon self-energy', &
      & 'El-photon self-energy', &
      & 'Transmission', &
      & 'Meir-Wingreen current', &
      & 'Layer currents', &
      & 'MPI gather']

  !> Verbosity above which progress lines are printed
  integer, parameter :: VBT = 70

  !> Width of the progress message field
  integer, parameter :: MSGWIDTH = 54

  !> Accumulated state of one timer
  type :: TTimer
    integer(int64) :: tstart = 0
    integer(int64) :: counts = 0
    integer :: ncalls = 0
    integer :: depth = 0
    logical :: printing = .false.
    logical :: suspended = .false.
  end type TTimer

  type(TTimer), save :: timerList(NTIMERS)

  !> Timer whose progress line is currently open (0 if none)
  integer, save :: openLine = 0

  !> Results of gather_timers (meaningful on the root of the communicator)
  real(dp), save :: tmin(NTIMERS) = 0.0_dp
  real(dp), save :: tavg(NTIMERS) = 0.0_dp
  real(dp), save :: tmax(NTIMERS) = 0.0_dp
  integer, save :: calls_max(NTIMERS) = 0
  integer, save :: nprocs_gathered = 1

contains

  !> Start a timer. If verbose > VBT a progress line is printed on the print node,
  !> using message if present or the timer label otherwise.
  subroutine set_timer(itimer, verbose, message)
    integer, intent(in) :: itimer
    integer, intent(in), optional :: verbose
    character(*), intent(in), optional :: message

    if (.not. valid(itimer)) return

    associate(tt => timerList(itimer))
      tt%depth = tt%depth + 1
      ! Nested calls of the same timer are accumulated only once
      if (tt%depth > 1) return

      tt%printing = .false.
      if (present(verbose)) then
        tt%printing = (id0 .and. verbose > VBT)
      end if
      if (tt%printing) then
        if (present(message)) then
          call open_progress_line(itimer, message)
        else
          call open_progress_line(itimer, trim(timer_names(itimer)))
        end if
      end if

      call system_clock(tt%tstart)
    end associate

  end subroutine set_timer

  !> Stop a timer, accumulate the elapsed time and close the progress line, if any
  subroutine stop_timer(itimer)
    integer, intent(in) :: itimer

    integer(int64) :: tstop, rate
    real(dp) :: elapsed

    if (.not. valid(itimer)) return

    associate(tt => timerList(itimer))
      if (tt%depth == 0) then
        if (id0) write(*,*) 'WARNING: stop_timer called for inactive timer ', &
            & trim(timer_names(itimer))
        return
      end if
      tt%depth = tt%depth - 1
      if (tt%depth > 0) return

      call system_clock(tstop, rate)
      tt%counts = tt%counts + (tstop - tt%tstart)
      tt%ncalls = tt%ncalls + 1

      if (tt%printing) then
        elapsed = real(tstop - tt%tstart, dp) / real(rate, dp)
        call close_progress_line(itimer, elapsed)
        tt%printing = .false.
      end if
    end associate

  end subroutine stop_timer

  !> Reset all timers
  subroutine reset_timers()
    integer :: ii

    do ii = 1, NTIMERS
      timerList(ii)%counts = 0
      timerList(ii)%ncalls = 0
      timerList(ii)%depth = 0
      timerList(ii)%printing = .false.
      timerList(ii)%suspended = .false.
    end do
    openLine = 0
    tmin = 0.0_dp
    tavg = 0.0_dp
    tmax = 0.0_dp
    calls_max = 0
    nprocs_gathered = 1

  end subroutine reset_timers

  !> Local accumulated time of a timer in seconds
  function get_timer_elapsed(itimer) result(elapsed)
    integer, intent(in) :: itimer
    real(dp) :: elapsed

    integer(int64) :: tdummy, rate

    elapsed = 0.0_dp
    if (.not. valid(itimer)) return
    call system_clock(tdummy, rate)
    elapsed = real(timerList(itimer)%counts, dp) / real(rate, dp)

  end function get_timer_elapsed

  !> Local number of calls of a timer
  function get_timer_calls(itimer) result(ncalls)
    integer, intent(in) :: itimer
    integer :: ncalls

    ncalls = 0
    if (.not. valid(itimer)) return
    ncalls = timerList(itimer)%ncalls

  end function get_timer_calls

#:if defined("MPI")
  !> Reduce min/avg/max of the timers over the processes of comm.
  !> Results are meaningful on the root of comm.
  subroutine gather_timers(comm)
    type(mpifx_comm), intent(in) :: comm

    real(dp) :: tloc(NTIMERS)
    integer :: cloc(NTIMERS)

    call local_times(tloc, cloc)

    call mpifx_reduce(comm, tloc, tmin, MPI_MIN)
    call mpifx_reduce(comm, tloc, tmax, MPI_MAX)
    call mpifx_reduce(comm, tloc, tavg, MPI_SUM)
    call mpifx_reduce(comm, cloc, calls_max, MPI_MAX)
    nprocs_gathered = comm%size
    tavg = tavg / real(nprocs_gathered, dp)

  end subroutine gather_timers
#:else
  !> Serial version: min = avg = max = local times
  subroutine gather_timers()

    real(dp) :: tloc(NTIMERS)
    integer :: cloc(NTIMERS)

    call local_times(tloc, cloc)
    tmin = tloc
    tavg = tloc
    tmax = tloc
    calls_max = cloc
    nprocs_gathered = 1

  end subroutine gather_timers
#:endif

  !> Write the summary of the last gather_timers() on the print node.
  !> Timers never called on any process are omitted.
  subroutine write_timers(unit)
    integer, intent(in), optional :: unit

    integer :: iu, ii
    real(dp) :: imbalance

    if (.not. id0) return
    iu = 6
    if (present(unit)) iu = unit

    write(iu, '(/,1x,a,i0,a)') 'libNEGF timers (wall time over ', nprocs_gathered, &
        & ' processes, seconds)'
    write(iu, '(1x,a24,a10,3a12,a10)') 'Label', 'Calls', 'Min', 'Avg', 'Max', 'Max/Avg'
    write(iu, '(1x,78("-"))')
    do ii = 1, NTIMERS
      if (calls_max(ii) == 0) cycle
      if (tavg(ii) > 0.0_dp) then
        imbalance = tmax(ii) / tavg(ii)
      else
        imbalance = 1.0_dp
      end if
      write(iu, '(1x,a24,i10,3f12.3,f10.2)') timer_names(ii), calls_max(ii), tmin(ii), &
          & tavg(ii), tmax(ii), imbalance
    end do
    write(iu, '(1x,78("-"))')
    write(iu, '(1x,a)') 'Calls: maximum over processes. Nested timers are included in their parents.'
    flush(iu)

  end subroutine write_timers

  ! ----------------------------------------------------------------------------------
  ! Private helpers
  ! ----------------------------------------------------------------------------------

  logical function valid(itimer)
    integer, intent(in) :: itimer
    valid = (itimer >= 1 .and. itimer <= NTIMERS)
  end function valid

  subroutine local_times(tloc, cloc)
    real(dp), intent(out) :: tloc(:)
    integer, intent(out) :: cloc(:)

    integer :: ii

    do ii = 1, NTIMERS
      tloc(ii) = get_timer_elapsed(ii)
      cloc(ii) = timerList(ii)%ncalls
    end do

  end subroutine local_times

  !> Print "message ......" without advancing. A line left open by another timer
  !> is terminated first and marked as suspended (as the old message_clock).
  subroutine open_progress_line(itimer, message)
    integer, intent(in) :: itimer
    character(*), intent(in) :: message

    integer :: ndots

    if (openLine > 0) then
      write(6,*)
      timerList(openLine)%suspended = .true.
    end if

    ndots = max(MSGWIDTH - len(message), 1)
    write(6, '(a)', advance='no') message // repeat('.', ndots)
    flush(6)
    openLine = itimer

  end subroutine open_progress_line

  !> Complete the progress line with the elapsed time
  subroutine close_progress_line(itimer, elapsed)
    integer, intent(in) :: itimer
    real(dp), intent(in) :: elapsed

    if (timerList(itimer)%suspended) then
      write(6, '(a)', advance='no') repeat('.', MSGWIDTH)
      timerList(itimer)%suspended = .false.
    end if
    write(6, '(f12.4,a)') elapsed, ' sec'
    flush(6)
    openLine = 0

  end subroutine close_progress_line

end module ln_timers
