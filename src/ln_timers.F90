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
!> The timer state is held in a TTimerSet (one per Tnegf instance). Timers are identified
!> by a fixed enumeration (timers%Gr, timers%SGF, ...), so that every process holds the
!> same list in the same order and results can be reduced over MPI.
!>
!> Statistics are accumulated separately for each calculation phase (density matrix,
!> current, ldos, ...), selected with set_timer_phase() by the driver routines, so that
!> times of different calculations are not mixed.
!>
!> For every timer and phase two kinds of statistics are collected:
!>  - per-process totals, whose min/avg/max over processes measure the load balance;
!>  - per-call min, max and mode. Min and max carry the (iE_path, iK) point where they
!>    occurred; the mode is obtained from a log-spaced histogram of the call times.
!>
!> Optionally timer_start/timer_stop also print a progress line on the print node, as the
!> old message_clock/write_clock, when the passed verbosity is larger than VBT.
!>
!> Note: timers are not thread-safe; call them outside OpenMP parallel regions.
module ln_timers

  use ln_precision, only : dp
  use iso_fortran_env, only : int64
  use mpi_globals, only : id0
#:if defined("MPI")
  use libmpifx_module, only : mpifx_comm, mpifx_reduce, mpifx_gather
  use mpi_f08, only : MPI_SUM, MPI_MIN, MPI_MAX
#:endif

  implicit none
  private

  public :: TTimerSet, TTimerSummary
  public :: timers, NTIMERS, timer_names
  public :: timer_phases, NPHASES, phase_names
  public :: timer_start, timer_stop, set_timer_phase, reset_timers
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

  !> Whether a call of the timer refers to a single (E, k) point. The others run over
  !> all the local points (self-energies, MPI reductions) and carry no point.
  logical, parameter :: point_resolved(NTIMERS) = [ &
      & .true., .true., .true., .true., .true., &
      & .false., .false., &
      & .true., .true., .true., &
      & .false.]

  !> Enumeration of the calculation phases
  type :: TPhaseEnum
    integer :: Density = 1
    integer :: Current = 2
    integer :: Ldos = 3
    integer :: PhononCurrent = 4
  end type TPhaseEnum

  type(TPhaseEnum), parameter :: timer_phases = TPhaseEnum()

  integer, parameter :: NPHASES = 4

  character(len=*), parameter :: phase_names(NPHASES) = [character(len=16) :: &
      & 'Density matrix', &
      & 'Current', &
      & 'Local DOS', &
      & 'Phonon current']

  !> Verbosity above which progress lines are printed
  integer, parameter :: VBT = 70

  !> Width of the progress message field
  integer, parameter :: MSGWIDTH = 54

  !> Histogram of the single-call times: log-spaced bins starting at 10**TLOG_MIN sec
  integer, parameter :: NBINS = 100
  integer, parameter :: BINS_PER_DECADE = 10
  real(dp), parameter :: TLOG_MIN = -6.0_dp

  !> State of a running timer
  type :: TTimerRun
    integer(int64) :: tstart = 0
    integer :: depth = 0
    integer :: phase = 1
    integer :: point(2) = 0
    logical :: printing = .false.
    logical :: suspended = .false.
  end type TTimerRun

  !> Accumulated statistics of one timer in one phase
  type :: TTimerStat
    integer(int64) :: counts = 0
    integer :: ncalls = 0
    real(dp) :: tcallMin = huge(1.0_dp)
    real(dp) :: tcallMax = 0.0_dp
    integer :: pointMin(2) = 0
    integer :: pointMax(2) = 0
    integer :: hist(NBINS) = 0
  end type TTimerStat

  !> Set of timers of a Tnegf instance
  type :: TTimerSet
    type(TTimerRun) :: run(NTIMERS)
    type(TTimerStat) :: stat(NTIMERS, NPHASES)
    !> Phase the next started timers are accumulated into
    integer :: phase = 1
    !> Timer whose progress line is currently open (0 if none)
    integer :: openLine = 0
  end type TTimerSet

  !> Statistics reduced over the processes, as produced by gather_timers().
  !> Meaningful on the root of the communicator.
  type :: TTimerSummary
    integer :: nprocs = 1
    !> Per-process totals: min, average and max over the processes
    real(dp) :: tmin(NTIMERS, NPHASES) = 0.0_dp
    real(dp) :: tavg(NTIMERS, NPHASES) = 0.0_dp
    real(dp) :: tmax(NTIMERS, NPHASES) = 0.0_dp
    !> Maximum number of calls over the processes
    integer :: callsMax(NTIMERS, NPHASES) = 0
    !> Per-call min, mode and max over all calls of all processes
    real(dp) :: tcallMin(NTIMERS, NPHASES) = 0.0_dp
    real(dp) :: tcallMode(NTIMERS, NPHASES) = 0.0_dp
    real(dp) :: tcallMax(NTIMERS, NPHASES) = 0.0_dp
    !> (iE_path, iK) of the fastest and slowest call
    integer :: pointMin(2, NTIMERS, NPHASES) = 0
    integer :: pointMax(2, NTIMERS, NPHASES) = 0
  end type TTimerSummary

contains

  !> Start a timer. point = (iE_path, iK) of the processed point, used for timers that
  !> are point_resolved. If verbose > VBT a progress line is printed on the print node,
  !> using message if present or the timer label otherwise.
  subroutine timer_start(ts, itimer, point, verbose, message)
    type(TTimerSet), intent(inout) :: ts
    integer, intent(in) :: itimer
    integer, intent(in) :: point(2)
    integer, intent(in), optional :: verbose
    character(*), intent(in), optional :: message

    if (.not. valid(itimer)) return

    associate(tr => ts%run(itimer))
      tr%depth = tr%depth + 1
      ! Nested calls of the same timer are accumulated only once
      if (tr%depth > 1) return

      tr%phase = ts%phase
      if (point_resolved(itimer)) then
        tr%point = point
      else
        tr%point = 0
      end if

      tr%printing = .false.
      if (present(verbose)) then
        tr%printing = (id0 .and. verbose > VBT)
      end if
      if (tr%printing) then
        if (present(message)) then
          call open_progress_line(ts, itimer, message)
        else
          call open_progress_line(ts, itimer, trim(timer_names(itimer)))
        end if
      end if

      call system_clock(tr%tstart)
    end associate

  end subroutine timer_start

  !> Stop a timer, accumulate the elapsed time and close the progress line, if any
  subroutine timer_stop(ts, itimer)
    type(TTimerSet), intent(inout) :: ts
    integer, intent(in) :: itimer

    integer(int64) :: tstop, rate
    real(dp) :: elapsed
    integer :: ib

    if (.not. valid(itimer)) return

    associate(tr => ts%run(itimer))
      if (tr%depth == 0) then
        if (id0) write(*,*) 'WARNING: timer_stop called for inactive timer ', &
            & trim(timer_names(itimer))
        return
      end if
      tr%depth = tr%depth - 1
      if (tr%depth > 0) return

      call system_clock(tstop, rate)
      elapsed = real(tstop - tr%tstart, dp) / real(rate, dp)

      associate(st => ts%stat(itimer, tr%phase))
        st%counts = st%counts + (tstop - tr%tstart)
        st%ncalls = st%ncalls + 1
        if (elapsed < st%tcallMin) then
          st%tcallMin = elapsed
          st%pointMin = tr%point
        end if
        if (elapsed > st%tcallMax) then
          st%tcallMax = elapsed
          st%pointMax = tr%point
        end if
        ib = time_bin(elapsed)
        st%hist(ib) = st%hist(ib) + 1
      end associate

      if (tr%printing) then
        call close_progress_line(ts, itimer, elapsed)
        tr%printing = .false.
      end if
    end associate

  end subroutine timer_stop

  !> Select the phase the timers started from now on are accumulated into
  subroutine set_timer_phase(ts, phase)
    type(TTimerSet), intent(inout) :: ts
    integer, intent(in) :: phase

    if (phase < 1 .or. phase > NPHASES) return
    ts%phase = phase

  end subroutine set_timer_phase

  !> Reset all timers
  subroutine reset_timers(ts)
    type(TTimerSet), intent(inout) :: ts

    ts = TTimerSet()

  end subroutine reset_timers

  !> Local accumulated time of a timer in seconds (in phase, or summed over all phases)
  function get_timer_elapsed(ts, itimer, phase) result(elapsed)
    type(TTimerSet), intent(in) :: ts
    integer, intent(in) :: itimer
    integer, intent(in), optional :: phase
    real(dp) :: elapsed

    integer(int64) :: tdummy, rate

    elapsed = 0.0_dp
    if (.not. valid(itimer)) return
    call system_clock(tdummy, rate)
    if (present(phase)) then
      elapsed = real(ts%stat(itimer, phase)%counts, dp) / real(rate, dp)
    else
      elapsed = real(sum(ts%stat(itimer, :)%counts), dp) / real(rate, dp)
    end if

  end function get_timer_elapsed

  !> Local number of calls of a timer (in phase, or summed over all phases)
  function get_timer_calls(ts, itimer, phase) result(ncalls)
    type(TTimerSet), intent(in) :: ts
    integer, intent(in) :: itimer
    integer, intent(in), optional :: phase
    integer :: ncalls

    ncalls = 0
    if (.not. valid(itimer)) return
    if (present(phase)) then
      ncalls = ts%stat(itimer, phase)%ncalls
    else
      ncalls = sum(ts%stat(itimer, :)%ncalls)
    end if

  end function get_timer_calls

#:if defined("MPI")
  !> Reduce the timer statistics over the processes of comm.
  !> All the processes of comm must call it. Results are meaningful on the root of comm.
  subroutine gather_timers(ts, summary, comm)
    type(TTimerSet), intent(in) :: ts
    type(TTimerSummary), intent(out) :: summary
    type(mpifx_comm), intent(in) :: comm

    real(dp) :: tloc(NTIMERS, NPHASES), tsum(NTIMERS, NPHASES)
    integer :: cloc(NTIMERS, NPHASES)
    integer :: hloc(NBINS, NTIMERS, NPHASES), hsum(NBINS, NTIMERS, NPHASES)
    ! Per-process extremes: tcallMin, iEmin, iKmin, tcallMax, iEmax, iKmax.
    ! The indices are stored as reals (exact) to need a single gather.
    real(dp) :: ext(6, NTIMERS, NPHASES)
    real(dp), allocatable :: extAll(:,:,:,:)

    call local_stats(ts, tloc, cloc, ext, hloc)
    tsum = 0.0_dp
    hsum = 0

    call mpifx_reduce(comm, tloc, summary%tmin, MPI_MIN)
    call mpifx_reduce(comm, tloc, summary%tmax, MPI_MAX)
    call mpifx_reduce(comm, tloc, tsum, MPI_SUM)
    call mpifx_reduce(comm, cloc, summary%callsMax, MPI_MAX)
    call mpifx_reduce(comm, hloc, hsum, MPI_SUM)
    summary%nprocs = comm%size
    summary%tavg = tsum / real(summary%nprocs, dp)

    if (comm%lead) then
      allocate(extAll(6, NTIMERS, NPHASES, comm%size))
    else
      allocate(extAll(6, NTIMERS, NPHASES, 1))
    end if
    call mpifx_gather(comm, ext, extAll)

    if (comm%lead) then
      call select_extremes(extAll, hsum, summary)
    end if

  end subroutine gather_timers
#:else
  !> Serial version: statistics of the local process
  subroutine gather_timers(ts, summary)
    type(TTimerSet), intent(in) :: ts
    type(TTimerSummary), intent(out) :: summary

    real(dp) :: tloc(NTIMERS, NPHASES)
    integer :: cloc(NTIMERS, NPHASES)
    integer :: hloc(NBINS, NTIMERS, NPHASES)
    real(dp) :: ext(6, NTIMERS, NPHASES)
    real(dp) :: extAll(6, NTIMERS, NPHASES, 1)

    call local_stats(ts, tloc, cloc, ext, hloc)
    summary%nprocs = 1
    summary%tmin = tloc
    summary%tavg = tloc
    summary%tmax = tloc
    summary%callsMax = cloc
    extAll(:,:,:,1) = ext
    call select_extremes(extAll, hloc, summary)

  end subroutine gather_timers
#:endif

  !> Write the summary on the print node, one block per phase with calls.
  !> Timers never called on any process are omitted.
  subroutine write_timers(summary, unit)
    type(TTimerSummary), intent(in) :: summary
    integer, intent(in), optional :: unit

    integer :: iu, ii, ip
    real(dp) :: imbalance

    if (.not. id0) return
    iu = 6
    if (present(unit)) iu = unit

    write(iu, '(/,1x,a,i0,a)') 'libNEGF timers (wall time in seconds, ', summary%nprocs, &
        & ' processes)'

    do ip = 1, NPHASES
      if (all(summary%callsMax(:, ip) == 0)) cycle

      write(iu, '(/,1x,a)') trim(phase_names(ip))
      write(iu, '(1x,a)') 'Totals per process:'
      write(iu, '(1x,a24,a10,3a12,a10)') 'Label', 'Calls', 'Min', 'Avg', 'Max', 'Max/Avg'
      write(iu, '(1x,78("-"))')
      do ii = 1, NTIMERS
        if (summary%callsMax(ii, ip) == 0) cycle
        if (summary%tavg(ii, ip) > 0.0_dp) then
          imbalance = summary%tmax(ii, ip) / summary%tavg(ii, ip)
        else
          imbalance = 1.0_dp
        end if
        write(iu, '(1x,a24,i10,3f12.3,f10.2)') timer_names(ii), summary%callsMax(ii, ip), &
            & summary%tmin(ii, ip), summary%tavg(ii, ip), summary%tmax(ii, ip), imbalance
      end do

      write(iu, '(1x,a)') 'Single call; (iE,iK) = energy point on path and k-point:'
      write(iu, '(1x,a24,a10,2a6,a10,a10,2a6)') 'Label', 'Min', 'iE', 'iK', 'Mode', 'Max', &
          & 'iE', 'iK'
      write(iu, '(1x,78("-"))')
      do ii = 1, NTIMERS
        if (summary%callsMax(ii, ip) == 0) cycle
        if (point_resolved(ii)) then
          write(iu, '(1x,a24,es10.3,2i6,es10.3,es10.3,2i6)') timer_names(ii), &
              & summary%tcallMin(ii, ip), summary%pointMin(:, ii, ip), &
              & summary%tcallMode(ii, ip), summary%tcallMax(ii, ip), &
              & summary%pointMax(:, ii, ip)
        else
          write(iu, '(1x,a24,es10.3,2a6,es10.3,es10.3,2a6)') timer_names(ii), &
              & summary%tcallMin(ii, ip), '-', '-', summary%tcallMode(ii, ip), &
              & summary%tcallMax(ii, ip), '-', '-'
        end if
      end do
      write(iu, '(1x,78("-"))')
    end do

    write(iu, '(1x,a)') 'Calls: maximum over processes. Nested timers are included in their parents.'
    write(iu, '(1x,a,i0,a)') 'Mode: centre of the most populated bin (', BINS_PER_DECADE, &
        & ' bins per decade).'
    flush(iu)

  end subroutine write_timers

  ! ----------------------------------------------------------------------------------
  ! Private helpers
  ! ----------------------------------------------------------------------------------

  logical function valid(itimer)
    integer, intent(in) :: itimer
    valid = (itimer >= 1 .and. itimer <= NTIMERS)
  end function valid

  !> Local totals, calls, extremes and histograms
  subroutine local_stats(ts, tloc, cloc, ext, hloc)
    type(TTimerSet), intent(in) :: ts
    real(dp), intent(out) :: tloc(:,:)
    integer, intent(out) :: cloc(:,:)
    real(dp), intent(out) :: ext(:,:,:)
    integer, intent(out) :: hloc(:,:,:)

    integer :: ii, ip

    do ip = 1, NPHASES
      do ii = 1, NTIMERS
        associate(st => ts%stat(ii, ip))
          tloc(ii, ip) = get_timer_elapsed(ts, ii, ip)
          cloc(ii, ip) = st%ncalls
          ext(1, ii, ip) = st%tcallMin
          ext(2:3, ii, ip) = real(st%pointMin, dp)
          ext(4, ii, ip) = st%tcallMax
          ext(5:6, ii, ip) = real(st%pointMax, dp)
          hloc(:, ii, ip) = st%hist
        end associate
      end do
    end do

  end subroutine local_stats

  !> Select the overall per-call extremes among the processes and the mode
  subroutine select_extremes(extAll, hist, summary)
    real(dp), intent(in) :: extAll(:,:,:,:)
    integer, intent(in) :: hist(:,:,:)
    type(TTimerSummary), intent(inout) :: summary

    integer :: ii, ip, ir

    do ip = 1, NPHASES
      do ii = 1, NTIMERS
        if (summary%callsMax(ii, ip) == 0) cycle
        ! Processes without calls hold tcallMin = huge and tcallMax = 0
        ir = minloc(extAll(1, ii, ip, :), dim=1)
        summary%tcallMin(ii, ip) = extAll(1, ii, ip, ir)
        summary%pointMin(:, ii, ip) = nint(extAll(2:3, ii, ip, ir))
        ir = maxloc(extAll(4, ii, ip, :), dim=1)
        summary%tcallMax(ii, ip) = extAll(4, ii, ip, ir)
        summary%pointMax(:, ii, ip) = nint(extAll(5:6, ii, ip, ir))
        summary%tcallMode(ii, ip) = mode_from_hist(hist(:, ii, ip))
      end do
    end do

  end subroutine select_extremes

  !> Histogram bin of a call time
  integer function time_bin(time)
    real(dp), intent(in) :: time

    time_bin = 1 + floor((log10(max(time, tiny(time))) - TLOG_MIN) * BINS_PER_DECADE)
    time_bin = min(max(time_bin, 1), NBINS)

  end function time_bin

  !> Geometric centre of the most populated bin
  real(dp) function mode_from_hist(hist)
    integer, intent(in) :: hist(:)

    integer :: ib

    mode_from_hist = 0.0_dp
    if (sum(hist) == 0) return
    ib = maxloc(hist, dim=1)
    mode_from_hist = 10.0_dp**(TLOG_MIN + (real(ib, dp) - 0.5_dp) / BINS_PER_DECADE)

  end function mode_from_hist

  !> Print "message ......" without advancing. A line left open by another timer
  !> is terminated first and marked as suspended (as the old message_clock).
  subroutine open_progress_line(ts, itimer, message)
    type(TTimerSet), intent(inout) :: ts
    integer, intent(in) :: itimer
    character(*), intent(in) :: message

    integer :: ndots

    if (ts%openLine > 0) then
      write(6,*)
      ts%run(ts%openLine)%suspended = .true.
    end if

    ndots = max(MSGWIDTH - len(message), 1)
    write(6, '(a)', advance='no') message // repeat('.', ndots)
    flush(6)
    ts%openLine = itimer

  end subroutine open_progress_line

  !> Complete the progress line with the elapsed time
  subroutine close_progress_line(ts, itimer, elapsed)
    type(TTimerSet), intent(inout) :: ts
    integer, intent(in) :: itimer
    real(dp), intent(in) :: elapsed

    if (ts%run(itimer)%suspended) then
      write(6, '(a)', advance='no') repeat('.', MSGWIDTH)
      ts%run(itimer)%suspended = .false.
    end if
    write(6, '(f12.4,a)') elapsed, ' sec'
    flush(6)
    ts%openLine = 0

  end subroutine close_progress_line

end module ln_timers
