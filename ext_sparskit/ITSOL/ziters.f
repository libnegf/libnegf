c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c         Basic Iterative Solvers with Reverse Communication           c
c----------------------------------------------------------------------c
c     This file currently has several basic iterative linear system    c
c     solvers. They are:                                               c
c     CG       -- Conjugate Gradient Method                            c
c     CGNR     -- Conjugate Gradient Method (Normal Residual equation) c
c     BCG      -- Bi-Conjugate Gradient Method                         c
c     DBCG     -- BCG with partial pivoting                            c
c     BCGSTAB  -- BCG stabilized                                       c
c     TFQMR    -- Transpose-Free Quasi-Minimum Residual method         c
c     FOM      -- Full Orthogonalization Method                        c
c     GMRES    -- Generalized Minimum RESidual method                  c
c     FGMRES   -- Flexible version of Generalized Minimum              c
c                 RESidual method                                      c
c     DQGMRES  -- Direct versions of Quasi Generalize Minimum          c
c                 Residual method                                      c
c----------------------------------------------------------------------c
c     They all have the following calling sequence:
c      subroutine solver(n, rhs, sol, ipar, fpar, w)
c      integer n, ipar(16)
c      real*8 rhs(n), sol(n), fpar(16), w(*)
c     Where
c     (1) 'n' is the size of the linear system,
c     (2) 'rhs' is the right-hand side of the linear system,
c     (3) 'sol' is the solution to the linear system,
c     (4) 'ipar' is an integer parameter array for the reverse
c     communication protocol,
c     (5) 'fpar' is an floating-point parameter array storing
c     information to and from the iterative solvers.
c     (6) 'w' is the work space (size is specified in ipar)
c
c     They are preconditioned iterative solvers with reverse
c     communication. The preconditioners can be applied from either
c     from left or right or both (specified by ipar(2), see below).
c
c     Author: Kesheng John Wu (kewu@mail.cs.umn.edu) 1993
c
c     NOTES:
c
c     (1) Work space required by each of the iterative solver
c     routines is as follows:
c       CG      == 5 * n
c       CGNR    == 5 * n
c       BCG     == 7 * n
c       DBCG    == 11 * n
c       BCGSTAB == 8 * n
c       TFQMR   == 11 * n
c       FOM     == (n+3)*(m+2) + (m+1)*m/2 (m = ipar(5), default m=15)
c       GMRES   == (n+3)*(m+2) + (m+1)*m/2 (m = ipar(5), default m=15)
c       FGMRES  == 2*n*(m+1) + (m+1)*m/2 + 3*m + 2 (m = ipar(5),
c                  default m=15)
c       DQGMRES == n + lb * (2*n+4) (lb=ipar(5)+1, default lb = 16)
c
c     (2) ALL iterative solvers require a user-supplied DOT-product
c     routine named DISTDOT. The prototype of DISTDOT is
c
c     real*8 function distdot(n,x,ix,y,iy)
c     integer n, ix, iy
c     real*8 x(1+(n-1)*ix), y(1+(n-1)*iy)
c
c     This interface of DISTDOT is exactly the same as that of
c     DDOT (or SDOT if real == real*8) from BLAS-1. It should have
c     same functionality as DDOT on a single processor machine. On a
c     parallel/distributed environment, each processor can perform
c     DDOT on the data it has, then perform a summation on all the
c     partial results.
c
c     (3) To use this set of routines under SPMD/MIMD program paradigm,
c     several things are to be noted: (a) 'n' should be the number of
c     vector elements of 'rhs' that is present on the local processor.
c     (b) if RHS(i) is on processor j, it is expected that SOL(i)
c     will be on the same processor, i.e. the vectors are distributed
c     to each processor in the same way. (c) the preconditioning and
c     stopping criteria specifications have to be the same on all
c     processor involved, ipar and fpar have to be the same on each
c     processor. (d) DISTDOT should be replaced by a distributed
c     dot-product function.
c
c     ..................................................................
c     Reverse Communication Protocols
c
c     When a reverse-communication routine returns, it could be either
c     that the routine has terminated or it simply requires the caller
c     to perform one matrix-vector multiplication. The possible matrices
c     that involve in the matrix-vector multiplications are:
c     A       (the matrix of the linear system),
c     A^T     (A transposed),
c     Ml^{-1} (inverse of the left preconditioner),
c     Ml^{-T} (inverse of the left preconditioner transposed),
c     Mr^{-1} (inverse of the right preconditioner),
c     Mr^{-T} (inverse of the right preconditioner transposed).
c     For all the matrix vector multiplication, v = A u. The input and
c     output vectors are supposed to be part of the work space 'w', and
c     the starting positions of them are stored in ipar(8:9), see below.
c
c     The array 'ipar' is used to store the information about the solver.
c     Here is the list of what each element represents:
c
c     ipar(1) -- status of the call/return.
c     A call to the solver with ipar(1) == 0 will initialize the
c     iterative solver. On return from the iterative solver, ipar(1)
c     carries the status flag which indicates the condition of the
c     return. The status information is divided into two categories,
c     (1) a positive value indicates the solver requires a matrix-vector
c     multiplication,
c     (2) a non-positive value indicates termination of the solver.
c     Here is the current definition:
c       1 == request a matvec with A,
c       2 == request a matvec with A^T,
c       3 == request a left preconditioner solve (Ml^{-1}),
c       4 == request a left preconditioner transposed solve (Ml^{-T}),
c       5 == request a right preconditioner solve (Mr^{-1}),
c       6 == request a right preconditioner transposed solve (Mr^{-T}),
c      10 == request the caller to perform stopping test,
c       0 == normal termination of the solver, satisfied the stopping
c            criteria,
c      -1 == termination because iteration number is greater than the
c            preset limit,
c      -2 == return due to insufficient work space,
c      -3 == return due to anticipated break-down / divide by zero,
c            in the case where Arnoldi procedure is used, additional
c            error code can be found in ipar(12), where ipar(12) is
c            the error code of orthogonalization procedure MGSRO:
c               -1: zero input vector
c               -2: input vector contains abnormal numbers
c               -3: input vector is a linear combination of others
c               -4: trianguler system in GMRES/FOM/etc. has nul rank
c      -4 == the values of fpar(1) and fpar(2) are both <= 0, the valid
c            ranges are 0 <= fpar(1) < 1, 0 <= fpar(2), and they can
c            not be zero at the same time
c      -9 == while trying to detect a break-down, an abnormal number is
c            detected.
c     -10 == return due to some non-numerical reasons, e.g. invalid
c            floating-point numbers etc.
c
c     ipar(2) -- status of the preconditioning:
c       0 == no preconditioning
c       1 == left preconditioning only
c       2 == right preconditioning only
c       3 == both left and right preconditioning
c
c     ipar(3) -- stopping criteria (details of this will be
c     discussed later).
c
c     ipar(4) -- number of elements in the array 'w'. if this is less
c     than the desired size, it will be over-written with the minimum
c     requirement. In which case the status flag ipar(1) = -2.
c
c     ipar(5) -- size of the Krylov subspace (used by GMRES and its
c     variants), e.g. GMRES(ipar(5)), FGMRES(ipar(5)),
c     DQGMRES(ipar(5)).
c
c     ipar(6) -- maximum number of matrix-vector multiplies, if not a
c     positive number the iterative solver will run till convergence
c     test is satisfied.
c
c     ipar(7) -- current number of matrix-vector multiplies. It is
c     incremented after each matrix-vector multiplication. If there
c     is preconditioning, the counter is incremented after the
c     preconditioning associated with each matrix-vector multiplication.
c
c     ipar(8) -- pointer to the input vector to the requested matrix-
c     vector multiplication.
c
c     ipar(9) -- pointer to the output vector of the requested matrix-
c     vector multiplication.
c
c     To perform v = A * u, it is assumed that u is w(ipar(8):ipar(8)+n-1)
c     and v is stored as w(ipar(9):ipar(9)+n-1).
c
c     ipar(10) -- the return address (used to determine where to go to
c     inside the iterative solvers after the caller has performed the
c     requested services).
c
c     ipar(11) -- the result of the external convergence test
c     On final return from the iterative solvers, this value
c     will be reflected by ipar(1) = 0 (details discussed later)
c
c     ipar(12) -- error code of MGSRO, it is
c                  1 if the input vector to MGSRO is linear combination
c                    of others,
c                  0 if MGSRO was successful,
c                 -1 if the input vector to MGSRO is zero,
c                 -2 if the input vector contains invalid number.
c
c     ipar(13) -- number of initializations. During each initilization
c                 residual norm is computed directly from M_l(b - A x).
c
c     ipar(14) to ipar(16) are NOT defined, they are NOT USED by
c     any iterative solver at this time.
c
c     Information about the error and tolerance are stored in the array
c     FPAR. So are some internal variables that need to be saved from
c     one iteration to the next one. Since the internal variables are
c     not the same for each routine, we only define the common ones.
c
c     The first two are input parameters:
c     fpar(1) -- the relative tolerance,
c     fpar(2) -- the absolute tolerance (details discussed later),
c
c     When the iterative solver terminates,
c     fpar(3) -- initial residual/error norm,
c     fpar(4) -- target residual/error norm,
c     fpar(5) -- current residual norm (if available),
c     fpar(6) -- current residual/error norm,
c     fpar(7) -- convergence rate,
c
c     fpar(8:10) are used by some of the iterative solvers to save some
c     internal information.
c
c     fpar(11) -- number of floating-point operations. The iterative
c     solvers will add the number of FLOPS they used to this variable,
c     but they do NOT initialize it, nor add the number of FLOPS due to
c     matrix-vector multiplications (since matvec is outside of the
c     iterative solvers). To insure the correct FLOPS count, the
c     caller should set fpar(11) = 0 before invoking the iterative
c     solvers and account for the number of FLOPS from matrix-vector
c     multiplications and preconditioners.
c
c     fpar(12:16) are not used in current implementation.
c
c     Whether the content of fpar(3), fpar(4) and fpar(6) are residual
c     norms or error norms depends on ipar(3). If the requested
c     convergence test is based on the residual norm, they will be
c     residual norms. If the caller want to test convergence based the
c     error norms (estimated by the norm of the modifications applied
c     to the approximate solution), they will be error norms.
c     Convergence rate is defined by (Fortran 77 statement)
c     fpar(7) = log10(fpar(3) / fpar(6)) / (ipar(7)-ipar(13))
c     If fpar(7) = 0.5, it means that approximately every 2 (= 1/0.5)
c     steps the residual/error norm decrease by a factor of 10.
c
c     ..................................................................
c     Stopping criteria,
c
c     An iterative solver may be terminated due to (1) satisfying
c     convergence test; (2) exceeding iteration limit; (3) insufficient
c     work space; (4) break-down. Checking of the work space is
c     only done in the initialization stage, i.e. when it is called with
c     ipar(1) == 0. A complete convergence test is done after each
c     update of the solutions. Other conditions are monitored
c     continuously.
c
c     With regard to the number of iteration, when ipar(6) is positive,
c     the current iteration number will be checked against it. If
c     current iteration number is greater the ipar(6) than the solver
c     will return with status -1. If ipar(6) is not positive, the
c     iteration will continue until convergence test is satisfied.
c
c     Two things may be used in the convergence tests, one is the
c     residual 2-norm, the other one is 2-norm of the change in the
c     approximate solution. The residual and the change in approximate
c     solution are from the preconditioned system (if preconditioning
c     is applied). The DQGMRES and TFQMR use two estimates for the
c     residual norms. The estimates are not accurate, but they are
c     acceptable in most of the cases. Generally speaking, the error
c     of the TFQMR's estimate is less accurate.
c
c     The convergence test type is indicated by ipar(3). There are four
c     type convergence tests: (1) tests based on the residual norm;
c     (2) tests based on change in approximate solution; (3) caller
c     does not care, the solver choose one from above two on its own;
c     (4) caller will perform the test, the solver should simply continue.
c     Here is the complete definition:
c      -2 == || dx(i) || <= rtol * || rhs || + atol
c      -1 == || dx(i) || <= rtol * || dx(1) || + atol
c       0 == solver will choose test 1 (next)
c       1 == || residual || <= rtol * || initial residual || + atol
c       2 == || residual || <= rtol * || rhs || + atol
c     999 == caller will perform the test
c     where dx(i) denote the change in the solution at the ith update.
c     ||.|| denotes 2-norm. rtol = fpar(1) and atol = fpar(2).
c
c     If the caller is to perform the convergence test, the outcome
c     should be stored in ipar(11).
c     ipar(11) = 0 -- failed the convergence test, iterative solver
c     should continue
c     ipar(11) = 1 -- satisfied convergence test, iterative solver
c     should perform the clean up job and stop.
c
c     Upon return with ipar(1) = 10,
c     ipar(8)  points to the starting position of the change in
c              solution Sx, where the actual solution of the step is
c              x_j = x_0 + M_r^{-1} Sx.
c              Exception: ipar(8) < 0, Sx = 0. It is mostly used by
c              GMRES and variants to indicate (1) Sx was not necessary,
c              (2) intermediate result of Sx is not computed.
c     ipar(9)  points to the starting position of a work vector that
c              can be used by the caller.
c
c     NOTE: the caller should allow the iterative solver to perform
c     clean up job after the external convergence test is satisfied,
c     since some of the iterative solvers do not directly
c     update the 'sol' array. A typical clean-up stage includes
c     performing the final update of the approximate solution and
c     computing the convergence information (e.g. values of fpar(3:7)).
c
c     NOTE: fpar(4) and fpar(6) are not set by the accelerators (the
c     routines implemented here) if ipar(3) = 999.
c
c     ..................................................................
c     Usage:
c
c     To start solving a linear system, the user needs to specify
c     first 6 elements of the ipar, and first 2 elements of fpar.
c     The user may optionally set fpar(11) = 0 if one wants to count
c     the number of floating-point operations. (Note: the iterative
c     solvers will only add the floating-point operations inside
c     themselves, the caller will have to add the FLOPS from the
c     matrix-vector multiplication routines and the preconditioning
c     routines in order to account for all the arithmetic operations.)
c
c     Here is an example:
c     ipar(1) = 0	! always 0 to start an iterative solver
c     ipar(2) = 2	! right preconditioning
c     ipar(3) = 1	! use convergence test scheme 1
c     ipar(4) = 10000	! the 'w' has 10,000 elements
c     ipar(5) = 10	! use *GMRES(10) (e.g. FGMRES(10))
c     ipar(6) = 100	! use at most 100 matvec's
c     fpar(1) = 1.0E-6	! relative tolerance 1.0E-6
c     fpar(2) = 1.0E-10 ! absolute tolerance 1.0E-10
c     fpar(11) = 0.0	! clearing the FLOPS counter
c
c     After the above specifications, one can start to call an iterative
c     solver, say BCG. Here is a piece of pseudo-code showing how it can
c     be done,
c
c 10   call zbcg(n,rhs,sol,ipar,fpar,w)
c      if (ipar(1).eq.1) then
c         call zamux(n,w(ipar(8)),w(ipar(9)),a,ja,ia)
c         goto 10
c      else if (ipar(1).eq.2) then
c         call zatmux(n,w(ipar(8)),w(ipar(9)),a,ja,ia)
c         goto 10
c      else if (ipar(1).eq.3) then
c         left preconditioner solver
c         goto 10
c      else if (ipar(1).eq.4) then
c         left preconditioner transposed solve
c         goto 10
c      else if (ipar(1).eq.5) then
c         right preconditioner solve
c         goto 10
c      else if (ipar(1).eq.6) then
c         right preconditioner transposed solve
c         goto 10
c      else if (ipar(1).eq.10) then
c         call my own stopping test routine
c         goto 10
c      else if (ipar(1).gt.0) then
c         ipar(1) is an unspecified code
c      else
c         the iterative solver terminated with code = ipar(1)
c      endif
c
c     This segment of pseudo-code assumes the matrix is in CSR format,
c     AMUX and ATMUX are two routines from the SPARSKIT MATVEC module.
c     They perform matrix-vector multiplications for CSR matrices,
c     where w(ipar(8)) is the first element of the input vectors to the
c     two routines, and w(ipar(9)) is the first element of the output
c     vectors from them. For simplicity, we did not show the name of
c     the routine that performs the preconditioning operations or the
c     convergence tests.
c$$$c-----------------------------------------------------------------------
c$$$      subroutine zcg(n, rhs, sol, ipar, fpar, w)
c$$$      implicit none
c$$$      integer n, ipar(16)
c$$$      complex*16 rhs(n), sol(n), w(n,*)
c$$$      real*8 fpar(16)
c$$$c-----------------------------------------------------------------------
c$$$c     This is a implementation of the Conjugate Gradient (CG) method
c$$$c     for solving linear system.
c$$$c
c$$$c     NOTE: This is not the PCG algorithm. It is a regular CG algorithm.
c$$$c     To be consistent with the other solvers, the preconditioners are
c$$$c     applied by performing Ml^{-1} A Mr^{-1} P in place of A P in the
c$$$c     CG algorithm.  PCG uses the preconditioner differently.
c$$$c
c$$$c     fpar(7) is used here internally to store <r, r>.
c$$$c     w(:,1) -- residual vector
c$$$c     w(:,2) -- P, the conjugate direction
c$$$c     w(:,3) -- A P, matrix multiply the conjugate direction
c$$$c     w(:,4) -- temporary storage for results of preconditioning
c$$$c     w(:,5) -- change in the solution (sol) is stored here until
c$$$c               termination of this solver
c$$$c-----------------------------------------------------------------------
c$$$c     external functions used
c$$$c
c$$$      complex*16 distdot
c$$$      logical stopbis, brkdn
c$$$      external distdot, stopbis, brkdn, bisinit
c$$$c
c$$$c     local variables
c$$$c
c$$$      integer i
c$$$      complex*16 alpha, one
c$$$      parameter(one=(1.0D0,0.D0))
c$$$      logical lp,rp
c$$$      save
c$$$c
c$$$c     check the status of the call
c$$$c
c$$$      if (ipar(1).le.0) ipar(10) = 0
c$$$      goto (10, 20, 40, 50, 60, 70, 80), ipar(10)
c$$$c
c$$$c     initialization
c$$$c
c$$$      call zbisinit(ipar,fpar,5*n,1,lp,rp,w)
c$$$      if (ipar(1).lt.0) return
c$$$c
c$$$c     request for matrix vector multiplication A*x in the initialization
c$$$c
c$$$      ipar(1) = 1
c$$$      ipar(8) = n+1
c$$$      ipar(9) = ipar(8) + n
c$$$      ipar(10) = 1
c$$$      do i = 1, n
c$$$         w(i,2) = sol(i)
c$$$      enddo
c$$$      return
c$$$ 10   ipar(7) = ipar(7) + 1
c$$$      ipar(13) = 1
c$$$      do i = 1, n
c$$$         w(i,2) = rhs(i) - w(i,3)
c$$$      enddo
c$$$      fpar(11) = fpar(11) + n
c$$$c
c$$$c     if left preconditioned
c$$$c
c$$$      if (lp) then
c$$$         ipar(1) = 3
c$$$         ipar(9) = 1
c$$$         ipar(10) = 2
c$$$         return
c$$$      endif
c$$$c
c$$$ 20   if (lp) then
c$$$         do i = 1, n
c$$$            w(i,2) = w(i,1)
c$$$         enddo
c$$$      else
c$$$         do i = 1, n
c$$$            w(i,1) = w(i,2)
c$$$         enddo
c$$$      endif
c$$$c
c$$$C      fpar(7) = distdot(n,w,1,w,1)
c$$$      fpar(3) = dznrm2(n,w,1)
c$$$      fpar(7) = fpar(3)**2
c$$$      fpar(11) = fpar(11) + 2 * n
c$$$      fpar(5) = fpar(3)
c$$$      if (abs(ipar(3)).eq.2) then
c$$$         fpar(4) = fpar(1) * dznrm2(n,rhs,1) + fpar(2)
c$$$         fpar(11) = fpar(11) + 2 * n
c$$$      else if (ipar(3).ne.999) then
c$$$         fpar(4) = fpar(1) * fpar(3) + fpar(2)
c$$$      endif
c$$$c
c$$$c     before iteration can continue, we need to compute A * p, which
c$$$c     includes the preconditioning operations
c$$$c
c$$$ 30   if (rp) then
c$$$         ipar(1) = 5
c$$$         ipar(8) = n + 1
c$$$         if (lp) then
c$$$            ipar(9) = ipar(8) + n
c$$$         else
c$$$            ipar(9) = 3*n + 1
c$$$         endif
c$$$         ipar(10) = 3
c$$$         return
c$$$      endif
c$$$c
c$$$ 40   ipar(1) = 1
c$$$      if (rp) then
c$$$         ipar(8) = ipar(9)
c$$$      else
c$$$         ipar(8) = n + 1
c$$$      endif
c$$$      if (lp) then
c$$$         ipar(9) = 3*n+1
c$$$      else
c$$$         ipar(9) = n+n+1
c$$$      endif
c$$$      ipar(10) = 4
c$$$      return
c$$$c
c$$$ 50   if (lp) then
c$$$         ipar(1) = 3
c$$$         ipar(8) = ipar(9)
c$$$         ipar(9) = n+n+1
c$$$         ipar(10) = 5
c$$$         return
c$$$      endif
c$$$c
c$$$c     continuing with the iterations
c$$$c
c$$$ 60   ipar(7) = ipar(7) + 1
c$$$c     alpha becomes complex and nothing works !!!  
c$$$      alpha = distdot(n,w(1,2),1,w(1,3),1)
c$$$c     ??????????????????????????????????????????? 
c$$$      fpar(11) = fpar(11) + 2*n
c$$$      if (brkdn(alpha,ipar)) goto 900
c$$$      alpha = fpar(7) * one / alpha
c$$$      do i = 1, n
c$$$         w(i,5) = w(i,5) + alpha * w(i,2)
c$$$         w(i,1) = w(i,1) - alpha * w(i,3)
c$$$      enddo
c$$$      fpar(11) = fpar(11) + 4*n
c$$$c
c$$$c     are we ready to terminate ?
c$$$c
c$$$      if (ipar(3).eq.999) then
c$$$         ipar(1) = 10
c$$$         ipar(8) = 4*n + 1
c$$$         ipar(9) = 3*n + 1
c$$$         ipar(10) = 6
c$$$         return
c$$$      endif
c$$$ 70   if (ipar(3).eq.999) then
c$$$         if (ipar(11).eq.1) goto 900
c$$$      else if (stopbis(n,ipar,1,fpar,w,w(1,2),alpha)) then
c$$$         goto 900
c$$$      endif
c$$$c
c$$$c     continue the iterations
c$$$c
c$$$      alpha = fpar(5)*fpar(5) / fpar(7)
c$$$      fpar(7) = fpar(5)*fpar(5)
c$$$      do i = 1, n
c$$$         w(i,2) = w(i,1) + alpha * w(i,2)
c$$$      enddo
c$$$      fpar(11) = fpar(11) + 2*n
c$$$      goto 30
c$$$c
c$$$c     clean up -- necessary to accommodate the right-preconditioning
c$$$c
c$$$ 900  if (rp) then
c$$$         if (ipar(1).lt.0) ipar(12) = ipar(1)
c$$$         ipar(1) = 5
c$$$         ipar(8) = 4*n + 1
c$$$         ipar(9) = ipar(8) - n
c$$$         ipar(10) = 7
c$$$         return
c$$$      endif
c$$$ 80   if (rp) then
c$$$         call ztidycg(n,ipar,fpar,sol,w(1,4))
c$$$      else
c$$$         call ztidycg(n,ipar,fpar,sol,w(1,5))
c$$$      endif
c$$$c
c$$$      return
c$$$      end
c$$$c-----end-of-cg
c$$$c-----------------------------------------------------------------------
c$$$      subroutine zcgnr(n,rhs,sol,ipar,fpar,wk)
c$$$      implicit none
c$$$      integer n, ipar(16)
c$$$      complex*16 rhs(n),sol(n),wk(n,*)
c$$$      real*8 fpar(16)
c$$$c-----------------------------------------------------------------------
c$$$c     CGNR -- Using CG algorithm solving A x = b by solving
c$$$c     Normal Residual equation: A^T A x = A^T b
c$$$c     As long as the matrix is not singular, A^T A is symmetric
c$$$c     positive definite, therefore CG (CGNR) will converge.
c$$$c
c$$$c     Usage of the work space:
c$$$c     wk(:,1) == residual vector R
c$$$c     wk(:,2) == the conjugate direction vector P
c$$$c     wk(:,3) == a scratch vector holds A P, or A^T R
c$$$c     wk(:,4) == a scratch vector holds intermediate results of the
c$$$c                preconditioning
c$$$c     wk(:,5) == a place to hold the modification to SOL
c$$$c
c$$$c     size of the work space WK is required = 5*n
c$$$c-----------------------------------------------------------------------
c$$$c     external functions used
c$$$c
c$$$      real*8 distdot
c$$$      logical stopbis, brkdn
c$$$      external distdot, stopbis, brkdn, bisinit
c$$$c
c$$$c     local variables
c$$$c
c$$$      integer i
c$$$      real*8 alpha, zz, zzm1
c$$$      logical lp, rp
c$$$      save
c$$$c
c$$$c     check the status of the call
c$$$c
c$$$      if (ipar(1).le.0) ipar(10) = 0
c$$$      goto (10, 20, 40, 50, 60, 70, 80, 90, 100, 110), ipar(10)
c$$$c
c$$$c     initialization
c$$$c
c$$$      call zbisinit(ipar,fpar,5*n,1,lp,rp,wk)
c$$$      if (ipar(1).lt.0) return
c$$$c
c$$$c     request for matrix vector multiplication A*x in the initialization
c$$$c
c$$$      ipar(1) = 1
c$$$      ipar(8) = 1
c$$$      ipar(9) = 1 + n
c$$$      ipar(10) = 1
c$$$      do i = 1, n
c$$$         wk(i,1) = sol(i)
c$$$      enddo
c$$$      return
c$$$ 10   ipar(7) = ipar(7) + 1
c$$$      ipar(13) = ipar(13) + 1
c$$$      do i = 1, n
c$$$         wk(i,1) = rhs(i) - wk(i,2)
c$$$      enddo
c$$$      fpar(11) = fpar(11) + n
c$$$c
c$$$c     if left preconditioned, precondition the initial residual
c$$$c
c$$$      if (lp) then
c$$$         ipar(1) = 3
c$$$         ipar(10) = 2
c$$$         return
c$$$      endif
c$$$c
c$$$ 20   if (lp) then
c$$$         do i = 1, n
c$$$            wk(i,1) = wk(i,2)
c$$$         enddo
c$$$      endif
c$$$c
c$$$      zz = distdot(n,wk,1,wk,1)
c$$$      fpar(11) = fpar(11) + 2 * n
c$$$      fpar(3) = sqrt(zz)
c$$$      fpar(5) = fpar(3)
c$$$      if (abs(ipar(3)).eq.2) then
c$$$         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
c$$$         fpar(11) = fpar(11) + 2 * n
c$$$      else if (ipar(3).ne.999) then
c$$$         fpar(4) = fpar(1) * fpar(3) + fpar(2)
c$$$      endif
c$$$c
c$$$c     normal iteration begins here, first half of the iteration
c$$$c     computes the conjugate direction
c$$$c
c$$$ 30   continue
c$$$c
c$$$c     request the caller to perform a A^T r --> wk(:,3)
c$$$c
c$$$      if (lp) then
c$$$         ipar(1) = 4
c$$$         ipar(8) = 1
c$$$         if (rp) then
c$$$            ipar(9) = n + n + 1
c$$$         else
c$$$            ipar(9) = 3*n + 1
c$$$         endif
c$$$         ipar(10) = 3
c$$$         return
c$$$      endif
c$$$c
c$$$ 40   ipar(1) = 2
c$$$      if (lp) then
c$$$         ipar(8) = ipar(9)
c$$$      else
c$$$         ipar(8) = 1
c$$$      endif
c$$$      if (rp) then
c$$$         ipar(9) = 3*n + 1
c$$$      else
c$$$         ipar(9) = n + n + 1
c$$$      endif
c$$$      ipar(10) = 4
c$$$      return
c$$$c
c$$$ 50   if (rp) then
c$$$         ipar(1) = 6
c$$$         ipar(8) = ipar(9)
c$$$         ipar(9) = n + n + 1
c$$$         ipar(10) = 5
c$$$         return
c$$$      endif
c$$$c
c$$$ 60   ipar(7) = ipar(7) + 1
c$$$      zzm1 = zz
c$$$      zz = distdot(n,wk(1,3),1,wk(1,3),1)
c$$$      fpar(11) = fpar(11) + 2 * n
c$$$      if (brkdn(zz,ipar)) goto 900
c$$$      if (ipar(7).gt.3) then
c$$$         alpha = zz / zzm1
c$$$         do i = 1, n
c$$$            wk(i,2) = wk(i,3) + alpha * wk(i,2)
c$$$         enddo
c$$$         fpar(11) = fpar(11) + 2 * n
c$$$      else
c$$$         do i = 1, n
c$$$            wk(i,2) = wk(i,3)
c$$$         enddo
c$$$      endif
c$$$c
c$$$c     before iteration can continue, we need to compute A * p
c$$$c
c$$$      if (rp) then
c$$$         ipar(1) = 5
c$$$         ipar(8) = n + 1
c$$$         if (lp) then
c$$$            ipar(9) = ipar(8) + n
c$$$         else
c$$$            ipar(9) = 3*n + 1
c$$$         endif
c$$$         ipar(10) = 6
c$$$         return
c$$$      endif
c$$$c
c$$$ 70   ipar(1) = 1
c$$$      if (rp) then
c$$$         ipar(8) = ipar(9)
c$$$      else
c$$$         ipar(8) = n + 1
c$$$      endif
c$$$      if (lp) then
c$$$        ipar(9) = 3*n+1
c$$$      else
c$$$         ipar(9) = n+n+1
c$$$      endif
c$$$      ipar(10) = 7
c$$$      return
c$$$c
c$$$ 80   if (lp) then
c$$$         ipar(1) = 3
c$$$         ipar(8) = ipar(9)
c$$$         ipar(9) = n+n+1
c$$$         ipar(10) = 8
c$$$         return
c$$$      endif
c$$$c
c$$$c     update the solution -- accumulate the changes in w(:,5)
c$$$c
c$$$ 90   ipar(7) = ipar(7) + 1
c$$$      alpha = distdot(n,wk(1,3),1,wk(1,3),1)
c$$$      fpar(11) = fpar(11) + 2 * n
c$$$      if (brkdn(alpha,ipar)) goto 900
c$$$      alpha = zz / alpha
c$$$      do i = 1, n
c$$$         wk(i,5) = wk(i,5) + alpha * wk(i,2)
c$$$         wk(i,1) = wk(i,1) - alpha * wk(i,3)
c$$$      enddo
c$$$      fpar(11) = fpar(11) + 4 * n
c$$$c
c$$$c     are we ready to terminate ?
c$$$c
c$$$      if (ipar(3).eq.999) then
c$$$         ipar(1) = 10
c$$$         ipar(8) = 4*n + 1
c$$$         ipar(9) = 3*n + 1
c$$$         ipar(10) = 9
c$$$         return
c$$$      endif
c$$$ 100  if (ipar(3).eq.999) then
c$$$         if (ipar(11).eq.1) goto 900
c$$$      else if (stopbis(n,ipar,1,fpar,wk,wk(1,2),alpha)) then
c$$$         goto 900
c$$$      endif
c$$$c
c$$$c     continue the iterations
c$$$c
c$$$      goto 30
c$$$c
c$$$c     clean up -- necessary to accommodate the right-preconditioning
c$$$c
c$$$ 900  if (rp) then
c$$$         if (ipar(1).lt.0) ipar(12) = ipar(1)
c$$$         ipar(1) = 5
c$$$         ipar(8) = 4*n + 1
c$$$         ipar(9) = ipar(8) - n
c$$$         ipar(10) = 10
c$$$         return
c$$$      endif
c$$$ 110  if (rp) then
c$$$         call ztidycg(n,ipar,fpar,sol,wk(1,4))
c$$$      else
c$$$         call ztidycg(n,ipar,fpar,sol,wk(1,5))
c$$$      endif
c$$$      return
c$$$      end
c$$$c-----end-of-cgnr
c$$$c-----------------------------------------------------------------------
c$$$      subroutine zbcg(n,rhs,sol,ipar,fpar,w)
c$$$      implicit none
c$$$      integer n, ipar(16)
c$$$      real*8  fpar(16)
c$$$      complex*16 rhs(n), sol(n), w(n,*)
c$$$c-----------------------------------------------------------------------
c$$$c     BCG: Bi Conjugate Gradient method. Programmed with reverse
c$$$c     communication, see the header for detailed specifications
c$$$c     of the protocol.
c$$$c
c$$$c     in this routine, before successful return, the fpar's are
c$$$c     fpar(3) == initial residual norm
c$$$c     fpar(4) == target residual norm
c$$$c     fpar(5) == current residual norm
c$$$c     fpar(7) == current rho (rhok = <r, s>)
c$$$c     fpar(8) == previous rho (rhokm1)
c$$$c
c$$$c     w(:,1) -- r, the residual
c$$$c     w(:,2) -- s, the dual of the 'r'
c$$$c     w(:,3) -- p, the projection direction
c$$$c     w(:,4) -- q, the dual of the 'p'
c$$$c     w(:,5) -- v, a scratch vector to store A*p, or A*q.
c$$$c     w(:,6) -- a scratch vector to store intermediate results
c$$$c     w(:,7) -- changes in the solution
c$$$c-----------------------------------------------------------------------
c$$$c     external routines used
c$$$c
c$$$      real*8 distdot
c$$$      logical stopbis,brkdn
c$$$      external distdot, stopbis, brkdn
c$$$c
c$$$      real*8 one
c$$$      parameter(one=1.0D0)
c$$$c
c$$$c     local variables
c$$$c
c$$$      integer i
c$$$      real*8 alpha
c$$$      logical rp, lp
c$$$      save
c$$$c
c$$$c     status of the program
c$$$c
c$$$      if (ipar(1).le.0) ipar(10) = 0
c$$$      goto (10, 20, 40, 50, 60, 70, 80, 90, 100, 110), ipar(10)
c$$$c
c$$$c     initialization, initial residual
c$$$c
c$$$      call zbisinit(ipar,fpar,7*n,1,lp,rp,w)
c$$$      if (ipar(1).lt.0) return
c$$$c
c$$$c     compute initial residual, request a matvecc
c$$$c
c$$$      ipar(1) = 1
c$$$      ipar(8) = 3*n+1
c$$$      ipar(9) = ipar(8) + n
c$$$      do i = 1, n
c$$$         w(i,4) = sol(i)
c$$$      enddo
c$$$      ipar(10) = 1
c$$$      return
c$$$ 10   ipar(7) = ipar(7) + 1
c$$$      ipar(13) = ipar(13) + 1
c$$$      do i = 1, n
c$$$         w(i,1) = rhs(i) - w(i,5)
c$$$      enddo
c$$$      fpar(11) = fpar(11) + n
c$$$      if (lp) then
c$$$         ipar(1) = 3
c$$$         ipar(8) = 1
c$$$         ipar(9) = n+1
c$$$         ipar(10) = 2
c$$$         return
c$$$      endif
c$$$c
c$$$ 20   if (lp) then
c$$$         do i = 1, n
c$$$            w(i,1) = w(i,2)
c$$$            w(i,3) = w(i,2)
c$$$            w(i,4) = w(i,2)
c$$$         enddo
c$$$      else
c$$$         do i = 1, n
c$$$            w(i,2) = w(i,1)
c$$$            w(i,3) = w(i,1)
c$$$            w(i,4) = w(i,1)
c$$$         enddo
c$$$      endif
c$$$c
c$$$      fpar(7) = distdot(n,w,1,w,1)
c$$$      fpar(11) = fpar(11) + 2 * n
c$$$      fpar(3) = sqrt(fpar(7))
c$$$      fpar(5) = fpar(3)
c$$$      fpar(8) = one
c$$$      if (abs(ipar(3)).eq.2) then
c$$$         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
c$$$         fpar(11) = fpar(11) + 2 * n
c$$$      else if (ipar(3).ne.999) then
c$$$         fpar(4) = fpar(1) * fpar(3) + fpar(2)
c$$$      endif
c$$$      if (ipar(3).ge.0.and.fpar(5).le.fpar(4)) then
c$$$         fpar(6) = fpar(5)
c$$$         goto 900
c$$$      endif
c$$$c
c$$$c     end of initialization, begin iteration, v = A p
c$$$c
c$$$ 30   if (rp) then
c$$$         ipar(1) = 5
c$$$         ipar(8) = n + n + 1
c$$$         if (lp) then
c$$$            ipar(9) = 4*n + 1
c$$$         else
c$$$            ipar(9) = 5*n + 1
c$$$         endif
c$$$         ipar(10) = 3
c$$$         return
c$$$      endif
c$$$c
c$$$ 40   ipar(1) = 1
c$$$      if (rp) then
c$$$         ipar(8) = ipar(9)
c$$$      else
c$$$         ipar(8) = n + n + 1
c$$$      endif
c$$$      if (lp) then
c$$$         ipar(9) = 5*n + 1
c$$$      else
c$$$         ipar(9) = 4*n + 1
c$$$      endif
c$$$      ipar(10) = 4
c$$$      return
c$$$c
c$$$ 50   if (lp) then
c$$$         ipar(1) = 3
c$$$         ipar(8) = ipar(9)
c$$$         ipar(9) = 4*n + 1
c$$$         ipar(10) = 5
c$$$         return
c$$$      endif
c$$$c
c$$$ 60   ipar(7) = ipar(7) + 1
c$$$      alpha = distdot(n,w(1,4),1,w(1,5),1)
c$$$      fpar(11) = fpar(11) + 2 * n
c$$$      if (brkdn(alpha,ipar)) goto 900
c$$$      alpha = fpar(7) / alpha
c$$$      do i = 1, n
c$$$         w(i,7) = w(i,7) + alpha * w(i,3)
c$$$         w(i,1) = w(i,1) - alpha * w(i,5)
c$$$      enddo
c$$$      fpar(11) = fpar(11) + 4 * n
c$$$      if (ipar(3).eq.999) then
c$$$         ipar(1) = 10
c$$$         ipar(8) = 6*n + 1
c$$$         ipar(9) = 5*n + 1
c$$$         ipar(10) = 6
c$$$         return
c$$$      endif
c$$$ 70   if (ipar(3).eq.999) then
c$$$         if (ipar(11).eq.1) goto 900
c$$$      else if (stopbis(n,ipar,1,fpar,w,w(1,3),alpha)) then
c$$$         goto 900
c$$$      endif
c$$$c
c$$$c     A^t * x
c$$$c
c$$$      if (lp) then
c$$$         ipar(1) = 4
c$$$         ipar(8) = 3*n + 1
c$$$         if (rp) then
c$$$            ipar(9) = 4*n + 1
c$$$         else
c$$$            ipar(9) = 5*n + 1
c$$$         endif
c$$$         ipar(10) = 7
c$$$         return
c$$$      endif
c$$$c
c$$$ 80   ipar(1) = 2
c$$$      if (lp) then
c$$$         ipar(8) = ipar(9)
c$$$      else
c$$$         ipar(8) = 3*n + 1
c$$$      endif
c$$$      if (rp) then
c$$$         ipar(9) = 5*n + 1
c$$$      else
c$$$         ipar(9) = 4*n + 1
c$$$      endif
c$$$      ipar(10) = 8
c$$$      return
c$$$c
c$$$ 90   if (rp) then
c$$$         ipar(1) = 6
c$$$         ipar(8) = ipar(9)
c$$$         ipar(9) = 4*n + 1
c$$$         ipar(10) = 9
c$$$         return
c$$$      endif
c$$$c
c$$$ 100  ipar(7) = ipar(7) + 1
c$$$      do i = 1, n
c$$$         w(i,2) = w(i,2) - alpha * w(i,5)
c$$$      enddo
c$$$      fpar(8) = fpar(7)
c$$$      fpar(7) = distdot(n,w,1,w(1,2),1)
c$$$      fpar(11) = fpar(11) + 4 * n
c$$$      if (brkdn(fpar(7), ipar)) return
c$$$      alpha = fpar(7) / fpar(8)
c$$$      do i = 1, n
c$$$         w(i,3) = w(i,1) + alpha * w(i,3)
c$$$         w(i,4) = w(i,2) + alpha * w(i,4)
c$$$      enddo
c$$$      fpar(11) = fpar(11) + 4 * n
c$$$c
c$$$c     end of the iterations
c$$$c
c$$$      goto 30
c$$$c
c$$$c     some clean up job to do
c$$$c
c$$$ 900  if (rp) then
c$$$         if (ipar(1).lt.0) ipar(12) = ipar(1)
c$$$         ipar(1) = 5
c$$$         ipar(8) = 6*n + 1
c$$$         ipar(9) = ipar(8) - n
c$$$         ipar(10) = 10
c$$$         return
c$$$      endif
c$$$ 110  if (rp) then
c$$$         call ztidycg(n,ipar,fpar,sol,w(1,6))
c$$$      else
c$$$         call ztidycg(n,ipar,fpar,sol,w(1,7))
c$$$      endif
c$$$      return
c$$$c-----end-of-bcg
c$$$      end
c$$$c-----------------------------------------------------------------------
c$$$      subroutine zbcgstab(n, rhs, sol, ipar, fpar, w)
c$$$      implicit none
c$$$      integer n, ipar(16)
c$$$      real*8 rhs(n), sol(n), w(n,8)
c$$$      real*8 fpar(16)
c$$$c-----------------------------------------------------------------------
c$$$c     BCGSTAB --- Bi Conjugate Gradient stabilized (BCGSTAB)
c$$$c     This is an improved BCG routine. (1) no matrix transpose is
c$$$c     involved. (2) the convergence is smoother.
c$$$c
c$$$c
c$$$c     Algorithm:
c$$$c     Initialization - r = b - A x, r0 = r, p = r, rho = (r0, r),
c$$$c     Iterate -
c$$$c     (1) v = A p
c$$$c     (2) alpha = rho / (r0, v)
c$$$c     (3) s = r - alpha v
c$$$c     (4) t = A s
c$$$c     (5) omega = (t, s) / (t, t)
c$$$c     (6) x = x + alpha * p + omega * s
c$$$c     (7) r = s - omega * t
c$$$c     convergence test goes here
c$$$c     (8) beta = rho, rho = (r0, r), beta = rho * alpha / (beta * omega)
c$$$c         p = r + beta * (p - omega * v)
c$$$c
c$$$c     in this routine, before successful return, the fpar's are
c$$$c     fpar(3) == initial (preconditionied-)residual norm
c$$$c     fpar(4) == target (preconditionied-)residual norm
c$$$c     fpar(5) == current (preconditionied-)residual norm
c$$$c     fpar(6) == current residual norm or error
c$$$c     fpar(7) == current rho (rhok = <r, r0>)
c$$$c     fpar(8) == alpha
c$$$c     fpar(9) == omega
c$$$c
c$$$c     Usage of the work space W
c$$$c     w(:, 1) = r0, the initial residual vector
c$$$c     w(:, 2) = r, current residual vector
c$$$c     w(:, 3) = s
c$$$c     w(:, 4) = t
c$$$c     w(:, 5) = v
c$$$c     w(:, 6) = p
c$$$c     w(:, 7) = tmp, used in preconditioning, etc.
c$$$c     w(:, 8) = delta x, the correction to the answer is accumulated
c$$$c               here, so that the right-preconditioning may be applied
c$$$c               at the end
c$$$c-----------------------------------------------------------------------
c$$$c     external routines used
c$$$c
c$$$      real*8 distdot
c$$$      logical stopbis, brkdn
c$$$      external distdot, stopbis, brkdn
c$$$c
c$$$      real*8 one
c$$$      parameter(one=1.0D0)
c$$$c
c$$$c     local variables
c$$$c
c$$$      integer i
c$$$      real*8 alpha,beta,rho,omega
c$$$      logical lp, rp
c$$$      save lp, rp
c$$$c
c$$$c     where to go
c$$$c
c$$$      if (ipar(1).gt.0) then
c$$$         goto (10, 20, 40, 50, 60, 70, 80, 90, 100, 110) ipar(10)
c$$$      else if (ipar(1).lt.0) then
c$$$         goto 900
c$$$      endif
c$$$c
c$$$c     call the initialization routine
c$$$c
c$$$      call zbisinit(ipar,fpar,8*n,1,lp,rp,w)
c$$$      if (ipar(1).lt.0) return
c$$$c
c$$$c     perform a matvec to compute the initial residual
c$$$c
c$$$      ipar(1) = 1
c$$$      ipar(8) = 1
c$$$      ipar(9) = 1 + n
c$$$      do i = 1, n
c$$$         w(i,1) = sol(i)
c$$$      enddo
c$$$      ipar(10) = 1
c$$$      return
c$$$ 10   ipar(7) = ipar(7) + 1
c$$$      ipar(13) = ipar(13) + 1
c$$$      do i = 1, n
c$$$         w(i,1) = rhs(i) - w(i,2)
c$$$      enddo
c$$$      fpar(11) = fpar(11) + n
c$$$      if (lp) then
c$$$         ipar(1) = 3
c$$$         ipar(10) = 2
c$$$         return
c$$$      endif
c$$$c
c$$$ 20   if (lp) then
c$$$         do i = 1, n
c$$$            w(i,1) = w(i,2)
c$$$            w(i,6) = w(i,2)
c$$$         enddo
c$$$      else
c$$$         do i = 1, n
c$$$            w(i,2) = w(i,1)
c$$$            w(i,6) = w(i,1)
c$$$         enddo
c$$$      endif
c$$$c
c$$$      fpar(7) = distdot(n,w,1,w,1)
c$$$      fpar(11) = fpar(11) + 2 * n
c$$$      fpar(5) = sqrt(fpar(7))
c$$$      fpar(3) = fpar(5)
c$$$      if (abs(ipar(3)).eq.2) then
c$$$         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
c$$$         fpar(11) = fpar(11) + 2 * n
c$$$      else if (ipar(3).ne.999) then
c$$$         fpar(4) = fpar(1) * fpar(3) + fpar(2)
c$$$      endif
c$$$      if (ipar(3).ge.0) fpar(6) = fpar(5)
c$$$      if (ipar(3).ge.0 .and. fpar(5).le.fpar(4) .and.
c$$$     +     ipar(3).ne.999) then
c$$$         goto 900
c$$$      endif
c$$$c
c$$$c     beginning of the iterations
c$$$c
c$$$c     Step (1), v = A p
c$$$ 30   if (rp) then
c$$$         ipar(1) = 5
c$$$         ipar(8) = 5*n+1
c$$$         if (lp) then
c$$$            ipar(9) = 4*n + 1
c$$$         else
c$$$            ipar(9) = 6*n + 1
c$$$         endif
c$$$         ipar(10) = 3
c$$$         return
c$$$      endif
c$$$c
c$$$ 40   ipar(1) = 1
c$$$      if (rp) then
c$$$         ipar(8) = ipar(9)
c$$$      else
c$$$         ipar(8) = 5*n+1
c$$$      endif
c$$$      if (lp) then
c$$$         ipar(9) = 6*n + 1
c$$$      else
c$$$         ipar(9) = 4*n + 1
c$$$      endif
c$$$      ipar(10) = 4
c$$$      return
c$$$ 50   if (lp) then
c$$$         ipar(1) = 3
c$$$         ipar(8) = ipar(9)
c$$$         ipar(9) = 4*n + 1
c$$$         ipar(10) = 5
c$$$         return
c$$$      endif
c$$$c
c$$$ 60   ipar(7) = ipar(7) + 1
c$$$c
c$$$c     step (2)
c$$$      alpha = distdot(n,w(1,1),1,w(1,5),1)
c$$$      fpar(11) = fpar(11) + 2 * n
c$$$      if (brkdn(alpha, ipar)) goto 900
c$$$      alpha = fpar(7) / alpha
c$$$      fpar(8) = alpha
c$$$c
c$$$c     step (3)
c$$$      do i = 1, n
c$$$         w(i,3) = w(i,2) - alpha * w(i,5)
c$$$      enddo
c$$$      fpar(11) = fpar(11) + 2 * n
c$$$c
c$$$c     Step (4): the second matvec -- t = A s
c$$$c
c$$$      if (rp) then
c$$$         ipar(1) = 5
c$$$         ipar(8) = n+n+1
c$$$         if (lp) then
c$$$            ipar(9) = ipar(8)+n
c$$$         else
c$$$            ipar(9) = 6*n + 1
c$$$         endif
c$$$         ipar(10) = 6
c$$$         return
c$$$      endif
c$$$c
c$$$ 70   ipar(1) = 1
c$$$      if (rp) then
c$$$         ipar(8) = ipar(9)
c$$$      else
c$$$         ipar(8) = n+n+1
c$$$      endif
c$$$      if (lp) then
c$$$         ipar(9) = 6*n + 1
c$$$      else
c$$$         ipar(9) = 3*n + 1
c$$$      endif
c$$$      ipar(10) = 7
c$$$      return
c$$$ 80   if (lp) then
c$$$         ipar(1) = 3
c$$$         ipar(8) = ipar(9)
c$$$         ipar(9) = 3*n + 1
c$$$         ipar(10) = 8
c$$$         return
c$$$      endif
c$$$ 90   ipar(7) = ipar(7) + 1
c$$$c
c$$$c     step (5)
c$$$      omega = distdot(n,w(1,4),1,w(1,4),1)
c$$$      fpar(11) = fpar(11) + n + n
c$$$      if (brkdn(omega,ipar)) goto 900
c$$$      omega = distdot(n,w(1,4),1,w(1,3),1) / omega
c$$$      fpar(11) = fpar(11) + n + n
c$$$      if (brkdn(omega,ipar)) goto 900
c$$$      fpar(9) = omega
c$$$      alpha = fpar(8)
c$$$c
c$$$c     step (6) and (7)
c$$$      do i = 1, n
c$$$         w(i,7) = alpha * w(i,6) + omega * w(i,3)
c$$$         w(i,8) = w(i,8) + w(i,7)
c$$$         w(i,2) = w(i,3) - omega * w(i,4)
c$$$      enddo
c$$$      fpar(11) = fpar(11) + 6 * n + 1
c$$$c
c$$$c     convergence test
c$$$      if (ipar(3).eq.999) then
c$$$         ipar(1) = 10
c$$$         ipar(8) = 7*n + 1
c$$$         ipar(9) = 6*n + 1
c$$$         ipar(10) = 9
c$$$         return
c$$$      endif
c$$$      if (stopbis(n,ipar,2,fpar,w(1,2),w(1,7),one))  goto 900
c$$$ 100  if (ipar(3).eq.999.and.ipar(11).eq.1) goto 900
c$$$c
c$$$c     step (8): computing new p and rho
c$$$      rho = fpar(7)
c$$$      fpar(7) = distdot(n,w(1,2),1,w(1,1),1)
c$$$      omega = fpar(9)
c$$$      beta = fpar(7) * fpar(8) / (fpar(9) * rho)
c$$$      do i = 1, n
c$$$         w(i,6) = w(i,2) + beta * (w(i,6) - omega * w(i,5))
c$$$      enddo
c$$$      fpar(11) = fpar(11) + 6 * n + 3
c$$$      if (brkdn(fpar(7),ipar)) goto 900
c$$$c
c$$$c     end of an iteration
c$$$c
c$$$      goto 30
c$$$c
c$$$c     some clean up job to do
c$$$c
c$$$ 900  if (rp) then
c$$$         if (ipar(1).lt.0) ipar(12) = ipar(1)
c$$$         ipar(1) = 5
c$$$         ipar(8) = 7*n + 1
c$$$         ipar(9) = ipar(8) - n
c$$$         ipar(10) = 10
c$$$         return
c$$$      endif
c$$$ 110  if (rp) then
c$$$         call ztidycg(n,ipar,fpar,sol,w(1,7))
c$$$      else
c$$$         call ztidycg(n,ipar,fpar,sol,w(1,8))
c$$$      endif
c$$$c
c$$$      return
c$$$c-----end-of-bcgstab
c$$$      end
c$$$c-----------------------------------------------------------------------
      subroutine ztfqmr(n, rhs, sol, ipar, fpar, w)
      implicit none
      integer n, ipar(16)
      complex*16 rhs(n), sol(n), w(n,*)
      real*8 fpar(16)
c-----------------------------------------------------------------------
c     TFQMR --- transpose-free Quasi-Minimum Residual method
c     This is developed from BCG based on the principle of Quasi-Minimum
c     Residual, and it is transpose-free.
c
c     It uses approximate residual norm.
c
c     Internally, the fpar's are used as following:
c     fpar(3) --- initial residual norm squared
c     fpar(4) --- target residual norm squared
c     fpar(5) --- current residual norm squared
c
c     w(:,1) -- R, residual
c     w(:,2) -- R0, the initial residual
c     w(:,3) -- W
c     w(:,4) -- Y
c     w(:,5) -- Z
c     w(:,6) -- A * Y
c     w(:,7) -- A * Z
c     w(:,8) -- V
c     w(:,9) -- D
c     w(:,10) -- intermediate results of preconditioning
c     w(:,11) -- changes in the solution
c-----------------------------------------------------------------------
c     external functions
c
      complex*16 zdistdot
      logical zbrkdn
      external zbrkdn, zdistdot
c
      real*8 one,zero
      parameter(one=1.0D0,zero=0.0D0)
c
c     local variables
c
      integer i
      logical lp, rp
      real*8 eta,sigma,theta,te,alpha,rho,tao
      save
c
c     status of the call (where to go)
c
      if (ipar(1).le.0) ipar(10) = 0
      goto (10,20,40,50,60,70,80,90,100,110), ipar(10)
c
c     initializations
c
      call zbisinit(ipar,fpar,11*n,2,lp,rp,w)
      if (ipar(1).lt.0) return
      ipar(1) = 1
      ipar(8) = 1
      ipar(9) = 1 + 6*n
      do i = 1, n
         w(i,1) = sol(i)
      enddo
      ipar(10) = 1
      return
 10   ipar(7) = ipar(7) + 1
      ipar(13) = ipar(13) + 1
      do i = 1, n
         w(i,1) = rhs(i) - w(i,7)
         w(i,9) = zero
      enddo
      fpar(11) = fpar(11) + n
c
      if (lp) then
         ipar(1) = 3
         ipar(9) = n+1
         ipar(10) = 2
         return
      endif
 20   continue
      if (lp) then
         do i = 1, n
            w(i,1) = w(i,2)
            w(i,3) = w(i,2)
         enddo
      else
         do i = 1, n
            w(i,2) = w(i,1)
            w(i,3) = w(i,1)
         enddo
      endif
c
      fpar(5) = dsqrt(dble(zdistdot(n,w,1,w,1)))
      fpar(3) = fpar(5)
      tao = fpar(5)
      fpar(11) = fpar(11) + n + n
      if (abs(ipar(3)).eq.2) then
         fpar(4)= fpar(1)*dsqrt(dble(zdistdot(n,rhs,1,rhs,1))) +fpar(2)
         fpar(11) = fpar(11) + n + n
      else if (ipar(3).ne.999) then
         fpar(4) = fpar(1) * tao + fpar(2)
      endif
      te = zero
      rho = zero
c
c     begin iteration
c
 30   sigma = rho
      rho = zdistdot(n,w(1,2),1,w(1,3),1)
      fpar(11) = fpar(11) + n + n
      if (zbrkdn(rho,ipar)) goto 900
      if (ipar(7).eq.1) then
         alpha = zero
      else
         alpha = rho / sigma
      endif
      do i = 1, n
         w(i,4) = w(i,3) + alpha * w(i,5)
      enddo
      fpar(11) = fpar(11) + n + n
c
c     A * x -- with preconditioning
c
      if (rp) then
         ipar(1) = 5
         ipar(8) = 3*n + 1
         if (lp) then
            ipar(9) = 5*n + 1
         else
            ipar(9) = 9*n + 1
         endif
         ipar(10) = 3
         return
      endif
c
 40   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = 3*n + 1
      endif
      if (lp) then
         ipar(9) = 9*n + 1
      else
         ipar(9) = 5*n + 1
      endif
      ipar(10) = 4
      return
c
 50   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = 5*n + 1
         ipar(10) = 5
         return
      endif
 60   ipar(7) = ipar(7) + 1
      do i = 1, n
         w(i,8) = w(i,6) + alpha * (w(i,7) + alpha * w(i,8))
      enddo
      sigma = zdistdot(n,w(1,2),1,w(1,8),1)
      fpar(11) = fpar(11) + 6 * n
      if (zbrkdn(sigma,ipar)) goto 900
      alpha = rho / sigma
      do i = 1, n
         w(i,5) = w(i,4) - alpha * w(i,8)
      enddo
      fpar(11) = fpar(11) + 2*n
c
c     the second A * x
c
      if (rp) then
         ipar(1) = 5
         ipar(8) = 4*n + 1
         if (lp) then
            ipar(9) = 6*n + 1
         else
            ipar(9) = 9*n + 1
         endif
         ipar(10) = 6
         return
      endif
c
 70   ipar(1) = 1
      if (rp) then
         ipar(8) = ipar(9)
      else
         ipar(8) = 4*n + 1
      endif
      if (lp) then
         ipar(9) = 9*n + 1
      else
         ipar(9) = 6*n + 1
      endif
      ipar(10) = 7
      return
c
 80   if (lp) then
         ipar(1) = 3
         ipar(8) = ipar(9)
         ipar(9) = 6*n + 1
         ipar(10) = 8
         return
      endif
 90   ipar(7) = ipar(7) + 1
      do i = 1, n
         w(i,3) = w(i,3) - alpha * w(i,6)
      enddo
c
c     update I
c
      theta = zdistdot(n,w(1,3),1,w(1,3),1) / (tao*tao)
      sigma = one / (one + theta)
      tao = tao * dsqrt(sigma * theta)
      fpar(11) = fpar(11) + 4*n + 6
      if (zbrkdn(tao,ipar)) goto 900
      eta = sigma * alpha
      sigma = te / alpha
      te = theta * eta
      do i = 1, n
         w(i,9) = w(i,4) + sigma * w(i,9)
         w(i,11) = w(i,11) + eta * w(i,9)
         w(i,3) = w(i,3) - alpha * w(i,7)
      enddo
      fpar(11) = fpar(11) + 6 * n + 6
      if (ipar(7).eq.1) then
         if (ipar(3).eq.-1) then
            fpar(3) = eta * dsqrt(dble(zdistdot(n,w(1,9),1,w(1,9),1)))
            fpar(4) = fpar(1)*fpar(3) + fpar(2)
            fpar(11) = fpar(11) + n + n + 4
         endif
      endif
c
c     update II
c
      theta = dble(zdistdot(n,w(1,3),1,w(1,3),1)) / (tao*tao)
      sigma = one / (one + theta)
      tao = tao * dsqrt(sigma * theta)
      fpar(11) = fpar(11) + 8 + 2*n
      if (zbrkdn(tao,ipar)) goto 900
      eta = sigma * alpha
      sigma = te / alpha
      te = theta * eta
      do i = 1, n
         w(i,9) = w(i,5) + sigma * w(i,9)
         w(i,11) = w(i,11) + eta * w(i,9)
      enddo
      fpar(11) = fpar(11) + 4*n + 3
c
c     this is the correct over-estimate
c      fpar(5) = sqrt(real(ipar(7)+1)) * tao
c     this is an approximation
      fpar(5) = tao
      if (ipar(3).eq.999) then
         ipar(1) = 10
         ipar(8) = 10*n + 1
         ipar(9) = 9*n + 1
         ipar(10) = 9
         return
      else if (ipar(3).lt.0) then
         fpar(6) = eta * dsqrt(dble(zdistdot(n,w(1,9),1,w(1,9),1)))
         fpar(11) = fpar(11) + n + n + 2
      else
         fpar(6) = fpar(5)
      endif
      if (fpar(6).gt.fpar(4) .and. (ipar(7).lt.ipar(6)
     +     .or. ipar(6).le.0)) goto 30
 100  if (ipar(3).eq.999.and.ipar(11).eq.0) goto 30
c
c     clean up
c
 900  if (rp) then
         if (ipar(1).lt.0) ipar(12) = ipar(1)
         ipar(1) = 5
         ipar(8) = 10*n + 1
         ipar(9) = ipar(8) - n
         ipar(10) = 10
         return
      endif
 110  if (rp) then
         call ztidycg(n,ipar,fpar,sol,w(1,10))
      else
         call ztidycg(n,ipar,fpar,sol,w(1,11))
      endif
c
      return
      end
c-----end-of-tfqmr
c$$$c-----------------------------------------------------------------------
c$$$      subroutine zfom(n, rhs, sol, ipar, fpar, w)
c$$$      implicit none
c$$$      integer n, ipar(16)
c$$$      real*8 rhs(n), sol(n), fpar(16), w(*)
c$$$c-----------------------------------------------------------------------
c$$$c     This a version of The Full Orthogonalization Method (FOM) 
c$$$c     implemented with reverse communication. It is a simple restart 
c$$$c     version of the FOM algorithm and is implemented with plane 
c$$$c     rotations similarly to GMRES.
c$$$c
c$$$c  parameters:
c$$$c  ----------- 
c$$$c     ipar(5) == the dimension of the Krylov subspace
c$$$c     after every ipar(5) iterations, the FOM will restart with
c$$$c     the updated solution and recomputed residual vector.
c$$$c
c$$$c     the work space in `w' is used as follows:
c$$$c     (1) the basis for the Krylov subspace, size n*(m+1);
c$$$c     (2) the Hessenberg matrix, only the upper triangular
c$$$c     portion of the matrix is stored, size (m+1)*m/2 + 1
c$$$c     (3) three vectors, all are of size m, they are
c$$$c     the cosine and sine of the Givens rotations, the third one holds
c$$$c     the residuals, it is of size m+1.
c$$$c
c$$$c     TOTAL SIZE REQUIRED == (n+3)*(m+2) + (m+1)*m/2
c$$$c     Note: m == ipar(5). The default value for this is 15 if
c$$$c     ipar(5) <= 1.
c$$$c-----------------------------------------------------------------------
c$$$c     external functions used
c$$$c
c$$$      real*8 distdot
c$$$      external distdot
c$$$c
c$$$      real*8 one, zero
c$$$      parameter(one=1.0D0, zero=0.0D0)
c$$$c
c$$$c     local variables, ptr and p2 are temporary pointers,
c$$$c     hes points to the Hessenberg matrix,
c$$$c     vc, vs point to the cosines and sines of the Givens rotations
c$$$c     vrn points to the vectors of residual norms, more precisely
c$$$c     the right hand side of the least square problem solved.
c$$$c
c$$$      integer i,ii,idx,k,m,ptr,p2,prs,hes,vc,vs,vrn
c$$$      real*8 alpha, c, s 
c$$$      logical lp, rp
c$$$      save
c$$$c
c$$$c     check the status of the call
c$$$c
c$$$      if (ipar(1).le.0) ipar(10) = 0
c$$$      goto (10, 20, 30, 40, 50, 60, 70) ipar(10)
c$$$c
c$$$c     initialization
c$$$c
c$$$      if (ipar(5).le.1) then
c$$$         m = 15
c$$$      else
c$$$         m = ipar(5)
c$$$      endif
c$$$      idx = n * (m+1)
c$$$      hes = idx + n
c$$$      vc = hes + (m+1) * m / 2 + 1
c$$$      vs = vc + m
c$$$      vrn = vs + m
c$$$      i = vrn + m + 1
c$$$      call zbisinit(ipar,fpar,i,1,lp,rp,w)
c$$$      if (ipar(1).lt.0) return
c$$$c
c$$$c     request for matrix vector multiplication A*x in the initialization
c$$$c
c$$$ 100  ipar(1) = 1
c$$$      ipar(8) = n+1
c$$$      ipar(9) = 1
c$$$      ipar(10) = 1
c$$$      k = 0
c$$$      do i = 1, n
c$$$         w(n+i) = sol(i)
c$$$      enddo
c$$$      return
c$$$ 10   ipar(7) = ipar(7) + 1
c$$$      ipar(13) = ipar(13) + 1
c$$$      if (lp) then
c$$$         do i = 1, n
c$$$            w(n+i) = rhs(i) - w(i)
c$$$         enddo
c$$$         ipar(1) = 3
c$$$         ipar(10) = 2
c$$$         return
c$$$      else
c$$$         do i = 1, n
c$$$            w(i) = rhs(i) - w(i)
c$$$         enddo
c$$$      endif
c$$$      fpar(11) = fpar(11) + n
c$$$c
c$$$ 20   alpha = sqrt(distdot(n,w,1,w,1))
c$$$      fpar(11) = fpar(11) + 2*n + 1
c$$$      if (ipar(7).eq.1 .and. ipar(3).ne.999) then
c$$$         if (abs(ipar(3)).eq.2) then
c$$$            fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
c$$$            fpar(11) = fpar(11) + 2*n
c$$$         else
c$$$            fpar(4) = fpar(1) * alpha + fpar(2)
c$$$         endif
c$$$         fpar(3) = alpha
c$$$      endif
c$$$      fpar(5) = alpha
c$$$      w(vrn+1) = alpha
c$$$      if (alpha.le.fpar(4) .and. ipar(3).ge.0 .and. ipar(3).ne.999) then
c$$$         ipar(1) = 0
c$$$         fpar(6) = alpha
c$$$         goto 300
c$$$      endif
c$$$      alpha = one / alpha
c$$$      do ii = 1, n
c$$$         w(ii) = alpha * w(ii)
c$$$      enddo
c$$$      fpar(11) = fpar(11) + n
c$$$c
c$$$c     request for (1) right preconditioning
c$$$c     (2) matrix vector multiplication
c$$$c     (3) left preconditioning
c$$$c
c$$$ 110  k = k + 1
c$$$      if (rp) then
c$$$         ipar(1) = 5
c$$$         ipar(8) = k*n - n + 1
c$$$         if (lp) then
c$$$            ipar(9) = k*n + 1
c$$$         else
c$$$            ipar(9) = idx + 1
c$$$         endif
c$$$         ipar(10) = 3
c$$$         return
c$$$      endif
c$$$c
c$$$ 30   ipar(1) = 1
c$$$      if (rp) then
c$$$         ipar(8) = ipar(9)
c$$$      else
c$$$         ipar(8) = (k-1)*n + 1
c$$$      endif
c$$$      if (lp) then
c$$$         ipar(9) = idx + 1
c$$$      else
c$$$         ipar(9) = 1 + k*n
c$$$      endif
c$$$      ipar(10) = 4
c$$$      return
c$$$c
c$$$ 40   if (lp) then
c$$$         ipar(1) = 3
c$$$         ipar(8) = ipar(9)
c$$$         ipar(9) = k*n + 1
c$$$         ipar(10) = 5
c$$$         return
c$$$      endif
c$$$c
c$$$c     Modified Gram-Schmidt orthogonalization procedure
c$$$c     temporary pointer 'ptr' is pointing to the current column of the
c$$$c     Hessenberg matrix. 'p2' points to the new basis vector
c$$$c
c$$$ 50   ipar(7) = ipar(7) + 1
c$$$      ptr = k * (k - 1) / 2 + hes
c$$$      p2 = ipar(9)
c$$$      call zmgsro(.false.,n,n,k+1,k+1,fpar(11),w,w(ptr+1),ipar(12))
c$$$      if (ipar(12).lt.0) goto 200
c$$$c
c$$$c     apply previous Givens rotations to column.
c$$$c
c$$$      p2 = ptr + 1
c$$$      do i = 1, k-1
c$$$         ptr = p2
c$$$         p2 = p2 + 1
c$$$         alpha = w(ptr)
c$$$         c = w(vc+i)
c$$$         s = w(vs+i)
c$$$         w(ptr) = c * alpha + s * w(p2)
c$$$         w(p2) = c * w(p2) - s * alpha
c$$$      enddo
c$$$c
c$$$c     end of one Arnoldi iteration, alpha will store the estimated
c$$$c     residual norm at current stage
c$$$c
c$$$      fpar(11) = fpar(11) + 6*k
c$$$
c$$$      prs = vrn+k
c$$$      alpha = fpar(5) 
c$$$      if (w(p2) .ne. zero) alpha = abs(w(p2+1)*w(prs)/w(p2)) 
c$$$      fpar(5) = alpha
c$$$c
c$$$      if (k.ge.m .or. (ipar(3).ge.0 .and. alpha.le.fpar(4))
c$$$     +     .or. (ipar(6).gt.0 .and. ipar(7).ge.ipar(6)))
c$$$     +     goto 200
c$$$c
c$$$      call zgivens(w(p2), w(p2+1), c, s)
c$$$      w(vc+k) = c
c$$$      w(vs+k) = s
c$$$      alpha = - s * w(prs)
c$$$      w(prs) = c * w(prs)
c$$$      w(prs+1) = alpha
c$$$c
c$$$      if (w(p2).ne.zero) goto 110
c$$$c
c$$$c     update the approximate solution, first solve the upper triangular
c$$$c     system, temporary pointer ptr points to the Hessenberg matrix,
c$$$c     prs points to the right-hand-side (also the solution) of the system.
c$$$c
c$$$ 200  ptr = hes + k * (k + 1) / 2
c$$$      prs = vrn + k
c$$$      if (w(ptr).eq.zero) then
c$$$c
c$$$c     if the diagonal elements of the last column is zero, reduce k by 1
c$$$c     so that a smaller trianguler system is solved
c$$$c
c$$$         k = k - 1
c$$$         if (k.gt.0) then
c$$$            goto 200
c$$$         else
c$$$            ipar(1) = -3
c$$$            ipar(12) = -4
c$$$            goto 300
c$$$         endif
c$$$      endif
c$$$      w(prs) = w(prs) / w(ptr)
c$$$      do i = k-1, 1, -1
c$$$         ptr = ptr - i - 1
c$$$         do ii = 1, i
c$$$            w(vrn+ii) = w(vrn+ii) - w(prs) * w(ptr+ii)
c$$$         enddo
c$$$         prs = prs - 1
c$$$         w(prs) = w(prs) / w(ptr)
c$$$      enddo
c$$$c
c$$$      do ii = 1, n
c$$$         w(ii) = w(ii) * w(prs)
c$$$      enddo
c$$$      do i = 1, k-1
c$$$         prs = prs + 1
c$$$         ptr = i*n
c$$$         do ii = 1, n
c$$$            w(ii) = w(ii) + w(prs) * w(ptr+ii)
c$$$         enddo
c$$$      enddo
c$$$      fpar(11) = fpar(11) + 2*(k-1)*n + n + k*(k+1)
c$$$c
c$$$      if (rp) then
c$$$         ipar(1) = 5
c$$$         ipar(8) = 1
c$$$         ipar(9) = idx + 1
c$$$         ipar(10) = 6
c$$$         return
c$$$      endif
c$$$c
c$$$ 60   if (rp) then
c$$$         do i = 1, n
c$$$            sol(i) = sol(i) + w(idx+i)
c$$$         enddo
c$$$      else
c$$$         do i = 1, n
c$$$            sol(i) = sol(i) + w(i)
c$$$         enddo
c$$$      endif
c$$$      fpar(11) = fpar(11) + n
c$$$c
c$$$c     process the complete stopping criteria
c$$$c
c$$$      if (ipar(3).eq.999) then
c$$$         ipar(1) = 10
c$$$         ipar(8) = -1
c$$$         ipar(9) = idx + 1
c$$$         ipar(10) = 7
c$$$         return
c$$$      else if (ipar(3).lt.0) then
c$$$         if (ipar(7).le.m+1) then
c$$$            fpar(3) = abs(w(vrn+1))
c$$$            if (ipar(3).eq.-1) fpar(4) = fpar(1)*fpar(3)+fpar(2)
c$$$         endif
c$$$         alpha = abs(w(vrn+k))
c$$$      endif
c$$$      fpar(6) = alpha
c$$$c
c$$$c     do we need to restart ?
c$$$c
c$$$ 70   if (ipar(12).ne.0) then
c$$$         ipar(1) = -3
c$$$         goto 300
c$$$      endif
c$$$      if (ipar(7).lt.ipar(6) .or. ipar(6).le.0) then
c$$$         if (ipar(3).ne.999) then
c$$$            if (fpar(6).gt.fpar(4)) goto 100
c$$$         else
c$$$            if (ipar(11).eq.0) goto 100
c$$$         endif
c$$$      endif
c$$$c
c$$$c     termination, set error code, compute convergence rate
c$$$c
c$$$      if (ipar(1).gt.0) then
c$$$         if (ipar(3).eq.999 .and. ipar(11).eq.1) then
c$$$            ipar(1) = 0
c$$$         else if (ipar(3).ne.999 .and. fpar(6).le.fpar(4)) then
c$$$            ipar(1) = 0
c$$$         else if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
c$$$            ipar(1) = -1
c$$$         else
c$$$            ipar(1) = -10
c$$$         endif
c$$$      endif
c$$$ 300  if (fpar(3).ne.zero .and. fpar(6).ne.zero .and.
c$$$     +     ipar(7).gt.ipar(13)) then
c$$$         fpar(7) = log10(fpar(3) / fpar(6)) / dble(ipar(7)-ipar(13))
c$$$      else
c$$$         fpar(7) = zero
c$$$      endif
c$$$      return
c$$$      end
c$$$c-----end-of-fom-------------------------------------------------------- 
c-----------------------------------------------------------------------
c$$$      subroutine zgmres(n, rhs, sol, ipar, fpar, w)
c$$$      implicit none
c$$$      integer n, ipar(16)
c$$$      real*8 rhs(n), sol(n), fpar(16), w(*)
c$$$c-----------------------------------------------------------------------
c$$$c     This a version of GMRES implemented with reverse communication.
c$$$c     It is a simple restart version of the GMRES algorithm.
c$$$c
c$$$c     ipar(5) == the dimension of the Krylov subspace
c$$$c     after every ipar(5) iterations, the GMRES will restart with
c$$$c     the updated solution and recomputed residual vector.
c$$$c
c$$$c     the space of the `w' is used as follows:
c$$$c     (1) the basis for the Krylov subspace, size n*(m+1);
c$$$c     (2) the Hessenberg matrix, only the upper triangular
c$$$c     portion of the matrix is stored, size (m+1)*m/2 + 1
c$$$c     (3) three vectors, all are of size m, they are
c$$$c     the cosine and sine of the Givens rotations, the third one holds
c$$$c     the residuals, it is of size m+1.
c$$$c
c$$$c     TOTAL SIZE REQUIRED == (n+3)*(m+2) + (m+1)*m/2
c$$$c     Note: m == ipar(5). The default value for this is 15 if
c$$$c     ipar(5) <= 1.
c$$$c-----------------------------------------------------------------------
c$$$c     external functions used
c$$$c
c$$$      real*8 distdot
c$$$      external distdot
c$$$c
c$$$      real*8 one, zero
c$$$      parameter(one=1.0D0, zero=0.0D0)
c$$$c
c$$$c     local variables, ptr and p2 are temporary pointers,
c$$$c     hess points to the Hessenberg matrix,
c$$$c     vc, vs point to the cosines and sines of the Givens rotations
c$$$c     vrn points to the vectors of residual norms, more precisely
c$$$c     the right hand side of the least square problem solved.
c$$$c
c$$$      integer i,ii,idx,k,m,ptr,p2,hess,vc,vs,vrn
c$$$      real*8 alpha, c, s
c$$$      logical lp, rp
c$$$      save
c$$$c
c$$$c     check the status of the call
c$$$c
c$$$      if (ipar(1).le.0) ipar(10) = 0
c$$$      goto (10, 20, 30, 40, 50, 60, 70) ipar(10)
c$$$c
c$$$c     initialization
c$$$c
c$$$      if (ipar(5).le.1) then
c$$$         m = 15
c$$$      else
c$$$         m = ipar(5)
c$$$      endif
c$$$      idx = n * (m+1)
c$$$      hess = idx + n
c$$$      vc = hess + (m+1) * m / 2 + 1
c$$$      vs = vc + m
c$$$      vrn = vs + m
c$$$      i = vrn + m + 1
c$$$      call zbisinit(ipar,fpar,i,1,lp,rp,w)
c$$$      if (ipar(1).lt.0) return
c$$$c
c$$$c     request for matrix vector multiplication A*x in the initialization
c$$$c
c$$$ 100  ipar(1) = 1
c$$$      ipar(8) = n+1
c$$$      ipar(9) = 1
c$$$      ipar(10) = 1
c$$$      k = 0
c$$$      do i = 1, n
c$$$         w(n+i) = sol(i)
c$$$      enddo
c$$$      return
c$$$ 10   ipar(7) = ipar(7) + 1
c$$$      ipar(13) = ipar(13) + 1
c$$$      if (lp) then
c$$$         do i = 1, n
c$$$            w(n+i) = rhs(i) - w(i)
c$$$         enddo
c$$$         ipar(1) = 3
c$$$         ipar(10) = 2
c$$$         return
c$$$      else
c$$$         do i = 1, n
c$$$            w(i) = rhs(i) - w(i)
c$$$         enddo
c$$$      endif
c$$$      fpar(11) = fpar(11) + n
c$$$c
c$$$ 20   alpha = sqrt(distdot(n,w,1,w,1))
c$$$      fpar(11) = fpar(11) + 2*n
c$$$      if (ipar(7).eq.1 .and. ipar(3).ne.999) then
c$$$         if (abs(ipar(3)).eq.2) then
c$$$            fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
c$$$            fpar(11) = fpar(11) + 2*n
c$$$         else
c$$$            fpar(4) = fpar(1) * alpha + fpar(2)
c$$$         endif
c$$$         fpar(3) = alpha
c$$$      endif
c$$$      fpar(5) = alpha
c$$$      w(vrn+1) = alpha
c$$$      if (alpha.le.fpar(4) .and. ipar(3).ge.0 .and. ipar(3).ne.999) then
c$$$         ipar(1) = 0
c$$$         fpar(6) = alpha
c$$$         goto 300
c$$$      endif
c$$$      alpha = one / alpha
c$$$      do ii = 1, n
c$$$         w(ii) = alpha * w(ii)
c$$$      enddo
c$$$      fpar(11) = fpar(11) + n
c$$$c
c$$$c     request for (1) right preconditioning
c$$$c     (2) matrix vector multiplication
c$$$c     (3) left preconditioning
c$$$c
c$$$ 110  k = k + 1
c$$$      if (rp) then
c$$$         ipar(1) = 5
c$$$         ipar(8) = k*n - n + 1
c$$$         if (lp) then
c$$$            ipar(9) = k*n + 1
c$$$         else
c$$$            ipar(9) = idx + 1
c$$$         endif
c$$$         ipar(10) = 3
c$$$         return
c$$$      endif
c$$$c
c$$$ 30   ipar(1) = 1
c$$$      if (rp) then
c$$$         ipar(8) = ipar(9)
c$$$      else
c$$$         ipar(8) = (k-1)*n + 1
c$$$      endif
c$$$      if (lp) then
c$$$         ipar(9) = idx + 1
c$$$      else
c$$$         ipar(9) = 1 + k*n
c$$$      endif
c$$$      ipar(10) = 4
c$$$      return
c$$$c
c$$$ 40   if (lp) then
c$$$         ipar(1) = 3
c$$$         ipar(8) = ipar(9)
c$$$         ipar(9) = k*n + 1
c$$$         ipar(10) = 5
c$$$         return
c$$$      endif
c$$$c
c$$$c     Modified Gram-Schmidt orthogonalization procedure
c$$$c     temporary pointer 'ptr' is pointing to the current column of the
c$$$c     Hessenberg matrix. 'p2' points to the new basis vector
c$$$c
c$$$ 50   ipar(7) = ipar(7) + 1
c$$$      ptr = k * (k - 1) / 2 + hess
c$$$      p2 = ipar(9)
c$$$      call zmgsro(.false.,n,n,k+1,k+1,fpar(11),w,w(ptr+1),ipar(12))
c$$$      if (ipar(12).lt.0) goto 200
c$$$c
c$$$c     apply previous Givens rotations and generate a new one to eliminate
c$$$c     the subdiagonal element.
c$$$c
c$$$      p2 = ptr + 1
c$$$      do i = 1, k-1
c$$$         ptr = p2
c$$$         p2 = p2 + 1
c$$$         alpha = w(ptr)
c$$$         c = w(vc+i)
c$$$         s = w(vs+i)
c$$$         w(ptr) = c * alpha + s * w(p2)
c$$$         w(p2) = c * w(p2) - s * alpha
c$$$      enddo
c$$$      call zgivens(w(p2), w(p2+1), c, s)
c$$$      w(vc+k) = c
c$$$      w(vs+k) = s
c$$$      p2 = vrn + k
c$$$      alpha = - s * w(p2)
c$$$      w(p2) = c * w(p2)
c$$$      w(p2+1) = alpha
c$$$c
c$$$c     end of one Arnoldi iteration, alpha will store the estimated
c$$$c     residual norm at current stage
c$$$c
c$$$      fpar(11) = fpar(11) + 6*k + 2
c$$$      alpha = abs(alpha)
c$$$      fpar(5) = alpha
c$$$      if (k.lt.m .and. .not.(ipar(3).ge.0 .and. alpha.le.fpar(4))
c$$$     +     .and. (ipar(6).le.0 .or. ipar(7).lt.ipar(6))) goto 110
c$$$c
c$$$c     update the approximate solution, first solve the upper triangular
c$$$c     system, temporary pointer ptr points to the Hessenberg matrix,
c$$$c     p2 points to the right-hand-side (also the solution) of the system.
c$$$c
c$$$ 200  ptr = hess + k * (k + 1) / 2
c$$$      p2 = vrn + k
c$$$      if (w(ptr).eq.zero) then
c$$$c
c$$$c     if the diagonal elements of the last column is zero, reduce k by 1
c$$$c     so that a smaller trianguler system is solved [It should only
c$$$c     happen when the matrix is singular, and at most once!]
c$$$c
c$$$         k = k - 1
c$$$         if (k.gt.0) then
c$$$            goto 200
c$$$         else
c$$$            ipar(1) = -3
c$$$            ipar(12) = -4
c$$$            goto 300
c$$$         endif
c$$$      endif
c$$$      w(p2) = w(p2) / w(ptr)
c$$$      do i = k-1, 1, -1
c$$$         ptr = ptr - i - 1
c$$$         do ii = 1, i
c$$$            w(vrn+ii) = w(vrn+ii) - w(p2) * w(ptr+ii)
c$$$         enddo
c$$$         p2 = p2 - 1
c$$$         w(p2) = w(p2) / w(ptr)
c$$$      enddo
c$$$c
c$$$      do ii = 1, n
c$$$         w(ii) = w(ii) * w(p2)
c$$$      enddo
c$$$      do i = 1, k-1
c$$$         ptr = i*n
c$$$         p2 = p2 + 1
c$$$         do ii = 1, n
c$$$            w(ii) = w(ii) + w(p2) * w(ptr+ii)
c$$$         enddo
c$$$      enddo
c$$$      fpar(11) = fpar(11) + 2*k*n - n + k*(k+1)
c$$$c
c$$$      if (rp) then
c$$$         ipar(1) = 5
c$$$         ipar(8) = 1
c$$$         ipar(9) = idx + 1
c$$$         ipar(10) = 6
c$$$         return
c$$$      endif
c$$$c
c$$$ 60   if (rp) then
c$$$         do i = 1, n
c$$$            sol(i) = sol(i) + w(idx+i)
c$$$         enddo
c$$$      else
c$$$         do i = 1, n
c$$$            sol(i) = sol(i) + w(i)
c$$$         enddo
c$$$      endif
c$$$      fpar(11) = fpar(11) + n
c$$$c
c$$$c     process the complete stopping criteria
c$$$c
c$$$      if (ipar(3).eq.999) then
c$$$         ipar(1) = 10
c$$$         ipar(8) = -1
c$$$         ipar(9) = idx + 1
c$$$         ipar(10) = 7
c$$$         return
c$$$      else if (ipar(3).lt.0) then
c$$$         if (ipar(7).le.m+1) then
c$$$            fpar(3) = abs(w(vrn+1))
c$$$            if (ipar(3).eq.-1) fpar(4) = fpar(1)*fpar(3)+fpar(2)
c$$$         endif
c$$$         fpar(6) = abs(w(vrn+k))
c$$$      else
c$$$         fpar(6) = fpar(5)
c$$$      endif
c$$$c
c$$$c     do we need to restart ?
c$$$c
c$$$ 70   if (ipar(12).ne.0) then
c$$$         ipar(1) = -3
c$$$         goto 300
c$$$      endif
c$$$      if ((ipar(7).lt.ipar(6) .or. ipar(6).le.0) .and.
c$$$     +     ((ipar(3).eq.999.and.ipar(11).eq.0) .or.
c$$$     +     (ipar(3).ne.999.and.fpar(6).gt.fpar(4)))) goto 100
c$$$c
c$$$c     termination, set error code, compute convergence rate
c$$$c
c$$$      if (ipar(1).gt.0) then
c$$$         if (ipar(3).eq.999 .and. ipar(11).eq.1) then
c$$$            ipar(1) = 0
c$$$         else if (ipar(3).ne.999 .and. fpar(6).le.fpar(4)) then
c$$$            ipar(1) = 0
c$$$         else if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
c$$$            ipar(1) = -1
c$$$         else
c$$$            ipar(1) = -10
c$$$         endif
c$$$      endif
c$$$ 300  if (fpar(3).ne.zero .and. fpar(6).ne.zero .and.
c$$$     +     ipar(7).gt.ipar(13)) then
c$$$         fpar(7) = log10(fpar(3) / fpar(6)) / dble(ipar(7)-ipar(13))
c$$$      else
c$$$         fpar(7) = zero
c$$$      endif
c$$$      return
c$$$      end
c$$$c-----end-of-gmres
c$$$c-----------------------------------------------------------------------
c$$$      subroutine zdqgmres(n, rhs, sol, ipar, fpar, w)
c$$$      implicit none
c$$$      integer n, ipar(16)
c$$$      complex*16 rhs(n), sol(n), w(*)
c$$$      real*8 fpar(16)
c$$$c-----------------------------------------------------------------------
c$$$c     DQGMRES -- Flexible Direct version of Quasi-General Minimum
c$$$c     Residual method. The right preconditioning can be varied from
c$$$c     step to step.
c$$$c
c$$$c     Work space used = n + lb * (2*n+4)
c$$$c     where lb = ipar(5) + 1 (default 16 if ipar(5) <= 1)
c$$$c-----------------------------------------------------------------------
c$$$c     local variables
c$$$c
c$$$      real*8 one,zero,deps
c$$$      parameter(one=1.0D0,zero=0.0D0)
c$$$      parameter(deps=1.0D-33)
c$$$c
c$$$      integer i,ii,j,jp1,j0,k,ptrw,ptrv,iv,iw,ic,is,ihm,ihd,lb,ptr
c$$$      real*8  alpha,beta,psi,distdot
c$$$      complex*16 c,s,alp,bet
c$$$      logical lp,rp,full
c$$$      external distdot,bisinit
c$$$      save
c$$$c
c$$$c     where to go
c$$$c
c$$$      if (ipar(1).le.0) ipar(10) = 0
c$$$      goto (10, 20, 40, 50, 60, 70) ipar(10)
c$$$c
c$$$c     locations of the work arrays. The arrangement is as follows:
c$$$c     w(1:n) -- temporary storage for the results of the preconditioning
c$$$c     w(iv+1:iw) -- the V's
c$$$c     w(iw+1:ic) -- the W's
c$$$c     w(ic+1:is) -- the COSINEs of the Givens rotations
c$$$c     w(is+1:ihm) -- the SINEs of the Givens rotations
c$$$c     w(ihm+1:ihd) -- the last column of the Hessenberg matrix
c$$$c     w(ihd+1:i) -- the inverse of the diagonals of the Hessenberg matrix
c$$$c
c$$$      if (ipar(5).le.1) then
c$$$         lb = 16
c$$$      else
c$$$         lb = ipar(5) + 1
c$$$      endif
c$$$      iv = n
c$$$      iw = iv + lb * n
c$$$      ic = iw + lb * n
c$$$      is = ic + lb
c$$$      ihm = is + lb
c$$$      ihd = ihm + lb
c$$$      i = ihd + lb
c$$$c
c$$$c     parameter check, initializations
c$$$c
c$$$      full = .false.
c$$$      call zbisinit(ipar,fpar,i,1,lp,rp,w)
c$$$      if (ipar(1).lt.0) return
c$$$      ipar(1) = 1
c$$$      if (lp) then
c$$$         do ii = 1, n
c$$$            w(iv+ii) = sol(ii)
c$$$         enddo
c$$$         ipar(8) = iv+1
c$$$         ipar(9) = 1
c$$$      else
c$$$         do ii = 1, n
c$$$            w(ii) = sol(ii)
c$$$         enddo
c$$$         ipar(8) = 1
c$$$         ipar(9) = iv+1
c$$$      endif
c$$$      ipar(10) = 1
c$$$      return
c$$$c
c$$$ 10   ipar(7) = ipar(7) + 1
c$$$      ipar(13) = ipar(13) + 1
c$$$      if (lp) then
c$$$         do i = 1, n
c$$$            w(i) = rhs(i) - w(i)
c$$$         enddo
c$$$         ipar(1) = 3
c$$$         ipar(8) = 1
c$$$         ipar(9) = iv+1
c$$$         ipar(10) = 2
c$$$         return
c$$$      else
c$$$         do i = 1, n
c$$$            w(iv+i) = rhs(i) - w(iv+i)
c$$$         enddo
c$$$      endif
c$$$      fpar(11) = fpar(11) + n
c$$$c
c$$$ 20   alpha = sqrt(distdot(n, w(iv+1), 1, w(iv+1), 1))
c$$$      fpar(11) = fpar(11) + (n + n)
c$$$      if (abs(ipar(3)).eq.2) then
c$$$         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
c$$$         fpar(11) = fpar(11) + 2*n
c$$$      else if (ipar(3).ne.999) then
c$$$         fpar(4) = fpar(1) * alpha + fpar(2)
c$$$      endif
c$$$      fpar(3) = alpha
c$$$      fpar(5) = alpha
c$$$      psi = alpha
c$$$      if (alpha.le.fpar(4)) then
c$$$         ipar(1) = 0
c$$$         fpar(6) = alpha
c$$$         goto 80
c$$$      endif
c$$$      alpha = one / alpha
c$$$      do i = 1, n
c$$$         w(iv+i) = w(iv+i) * alpha
c$$$      enddo
c$$$      fpar(11) = fpar(11) + n
c$$$      j = 0
c$$$c
c$$$c     iterations start here
c$$$c
c$$$ 30   j = j + 1
c$$$      if (j.gt.lb) j = j - lb
c$$$      jp1 = j + 1
c$$$      if (jp1.gt.lb) jp1 = jp1 - lb
c$$$      ptrv = iv + (j-1)*n + 1
c$$$      ptrw = iv + (jp1-1)*n + 1
c$$$      if (.not.full) then
c$$$         if (j.gt.jp1) full = .true.
c$$$      endif
c$$$      if (full) then
c$$$         j0 = jp1+1
c$$$         if (j0.gt.lb) j0 = j0 - lb
c$$$      else
c$$$         j0 = 1
c$$$      endif
c$$$c
c$$$c     request the caller to perform matrix-vector multiplication and
c$$$c     preconditioning
c$$$c
c$$$      if (rp) then
c$$$         ipar(1) = 5
c$$$         ipar(8) = ptrv
c$$$         ipar(9) = ptrv + iw - iv
c$$$         ipar(10) = 3
c$$$         return
c$$$      else
c$$$         do i = 0, n-1
c$$$            w(ptrv+iw-iv+i) = w(ptrv+i)
c$$$         enddo
c$$$      endif
c$$$c
c$$$ 40   ipar(1) = 1
c$$$      if (rp) then
c$$$         ipar(8) = ipar(9)
c$$$      else
c$$$         ipar(8) = ptrv
c$$$      endif
c$$$      if (lp) then
c$$$         ipar(9) = 1
c$$$      else
c$$$         ipar(9) = ptrw
c$$$      endif
c$$$      ipar(10) = 4
c$$$      return
c$$$c
c$$$ 50   if (lp) then
c$$$         ipar(1) = 3
c$$$         ipar(8) = ipar(9)
c$$$         ipar(9) = ptrw
c$$$         ipar(10) = 5
c$$$         return
c$$$      endif
c$$$c
c$$$c     compute the last column of the Hessenberg matrix
c$$$c     modified Gram-schmidt procedure, orthogonalize against (lb-1)
c$$$c     previous vectors
c$$$c
c$$$ 60   continue
c$$$      call zmgsro(full,n,n,lb,jp1,fpar(11),w(iv+1),w(ihm+1),ipar(12))
c$$$      if (ipar(12).lt.0) then
c$$$         ipar(1) = -3
c$$$         goto 80
c$$$      endif
c$$$      beta = w(ihm+jp1)
c$$$c
c$$$c     incomplete factorization (QR factorization through Givens rotations)
c$$$c     (1) apply previous rotations [(lb-1) of them]
c$$$c     (2) generate a new rotation
c$$$c
c$$$      if (full) then
c$$$         w(ihm+jp1) = w(ihm+j0) * w(is+jp1)
c$$$         w(ihm+j0) = w(ihm+j0) * w(ic+jp1)
c$$$      endif
c$$$      i = j0
c$$$      do while (i.ne.j)
c$$$         k = i+1
c$$$         if (k.gt.lb) k = k - lb
c$$$         c = w(ic+i)
c$$$         s = w(is+i)
c$$$         alp = w(ihm+i)
c$$$         w(ihm+i) = c * alp + s * w(ihm+k)
c$$$         w(ihm+k) = c * w(ihm+k) - s * alp
c$$$         i = k
c$$$      enddo
c$$$      call zgivens(w(ihm+j), beta, c, s)
c$$$      if (full) then
c$$$         fpar(11) = fpar(11) + 6 * lb
c$$$      else
c$$$         fpar(11) = fpar(11) + 6 * j
c$$$      endif
c$$$c
c$$$c     detect whether diagonal element of this column is zero
c$$$c
c$$$      if (abs(w(ihm+j)).lt.deps) then
c$$$         ipar(1) = -3
c$$$         goto 80
c$$$      endif
c$$$      w(ihd+j) = one / w(ihm+j)
c$$$      w(ic+j) = c
c$$$      w(is+j) = s
c$$$c
c$$$c     update the W's (the conjugate directions) -- essentially this is one
c$$$c     step of triangular solve.
c$$$c
c$$$      ptrw = iw+(j-1)*n + 1
c$$$      if (full) then
c$$$         do i = j+1, lb
c$$$            alp = -w(ihm+i)*w(ihd+i)
c$$$            ptr = iw+(i-1)*n+1
c$$$            do ii = 0, n-1
c$$$               w(ptrw+ii) = w(ptrw+ii) + alp * w(ptr+ii)
c$$$            enddo
c$$$         enddo
c$$$      endif
c$$$      do i = 1, j-1
c$$$         alp = -w(ihm+i)*w(ihd+i)
c$$$         ptr = iw+(i-1)*n+1
c$$$         do ii = 0, n-1
c$$$            w(ptrw+ii) = w(ptrw+ii) + alp * w(ptr+ii)
c$$$         enddo
c$$$      enddo
c$$$c
c$$$c     update the solution to the linear system
c$$$c
c$$$      alp = psi * c * w(ihd+j)
c$$$      psi = - s * psi
c$$$      do i = 1, n
c$$$         sol(i) = sol(i) + alp * w(ptrw-1+i)
c$$$      enddo
c$$$      if (full) then
c$$$         fpar(11) = fpar(11) + lb * (n+n)
c$$$      else
c$$$         fpar(11) = fpar(11) + j * (n+n)
c$$$      endif
c$$$c
c$$$c     determine whether to continue,
c$$$c     compute the desired error/residual norm
c$$$c
c$$$      ipar(7) = ipar(7) + 1
c$$$      fpar(5) = abs(psi)
c$$$      if (ipar(3).eq.999) then
c$$$         ipar(1) = 10
c$$$         ipar(8) = -1
c$$$         ipar(9) = 1
c$$$         ipar(10) = 6
c$$$         return
c$$$      endif
c$$$      if (ipar(3).lt.0) then
c$$$         alpha = abs(alpha)
c$$$         if (ipar(7).eq.2 .and. ipar(3).eq.-1) then
c$$$            fpar(3) = alpha*sqrt(distdot(n, w(ptrw), 1, w(ptrw), 1))
c$$$            fpar(4) = fpar(1) * fpar(3) + fpar(2)
c$$$            fpar(6) = fpar(3)
c$$$         else
c$$$            fpar(6) = alpha*sqrt(distdot(n, w(ptrw), 1, w(ptrw), 1))
c$$$         endif
c$$$         fpar(11) = fpar(11) + 2 * n
c$$$      else
c$$$         fpar(6) = fpar(5)
c$$$      endif
c$$$      if (ipar(1).ge.0 .and. fpar(6).gt.fpar(4) .and. (ipar(6).le.0
c$$$     +     .or. ipar(7).lt.ipar(6))) goto 30
c$$$ 70   if (ipar(3).eq.999 .and. ipar(11).eq.0) goto 30
c$$$c
c$$$c     clean up the iterative solver
c$$$c
c$$$ 80   fpar(7) = zero
c$$$      if (fpar(3).ne.zero .and. fpar(6).ne.zero .and.
c$$$     +     ipar(7).gt.ipar(13))
c$$$     +     fpar(7) = log10(fpar(3) / fpar(6)) / dble(ipar(7)-ipar(13))
c$$$      if (ipar(1).gt.0) then
c$$$         if (ipar(3).eq.999 .and. ipar(11).ne.0) then
c$$$            ipar(1) = 0
c$$$         else if (fpar(6).le.fpar(4)) then
c$$$            ipar(1) = 0
c$$$         else if (ipar(6).gt.0 .and. ipar(7).ge.ipar(6)) then
c$$$            ipar(1) = -1
c$$$         else
c$$$            ipar(1) = -10
c$$$         endif
c$$$      endif
c$$$      return
c$$$      end
c$$$c-----end-of-dqgmres
c-----------------------------------------------------------------------
c$$$      subroutine zfgmres(n, rhs, sol, ipar, fpar, w)
c$$$      implicit none
c$$$      integer n, ipar(16)
c$$$      real*8 rhs(n), sol(n), fpar(16), w(*)
c$$$c-----------------------------------------------------------------------
c$$$c     This a version of FGMRES implemented with reverse communication.
c$$$c
c$$$c     ipar(5) == the dimension of the Krylov subspace
c$$$c
c$$$c     the space of the `w' is used as follows:
c$$$c     >> V: the bases for the Krylov subspace, size n*(m+1);
c$$$c     >> W: the above bases after (left-)multiplying with the
c$$$c     right-preconditioner inverse, size m*n;
c$$$c     >> a temporary vector of size n;
c$$$c     >> the Hessenberg matrix, only the upper triangular portion
c$$$c     of the matrix is stored, size (m+1)*m/2 + 1
c$$$c     >> three vectors, first two are of size m, they are the cosine
c$$$c     and sine of the Givens rotations, the third one holds the
c$$$c     residuals, it is of size m+1.
c$$$c
c$$$c     TOTAL SIZE REQUIRED == n*(2m+1) + (m+1)*m/2 + 3*m + 2
c$$$c     Note: m == ipar(5). The default value for this is 15 if
c$$$c     ipar(5) <= 1.
c$$$c-----------------------------------------------------------------------
c$$$c     external functions used
c$$$c
c$$$      real*8 zdistdot
c$$$      external zdistdot
c$$$c
c$$$      real*8 one, zero
c$$$      parameter(one=1.0D0, zero=0.0D0)
c$$$c
c$$$c     local variables, ptr and p2 are temporary pointers,
c$$$c     hess points to the Hessenberg matrix,
c$$$c     vc, vs point to the cosines and sines of the Givens rotations
c$$$c     vrn points to the vectors of residual norms, more precisely
c$$$c     the right hand side of the least square problem solved.
c$$$c
c$$$      integer i,ii,idx,iz,k,m,ptr,p2,hess,vc,vs,vrn
c$$$      real*8 alpha, c, s
c$$$      logical lp, rp
c$$$      save
c$$$c
c$$$c     check the status of the call
c$$$c
c$$$      if (ipar(1).le.0) ipar(10) = 0
c$$$      goto (10, 20, 30, 40, 50, 60) ipar(10)
c$$$c
c$$$c     initialization
c$$$c
c$$$      if (ipar(5).le.1) then
c$$$         m = 15
c$$$      else
c$$$         m = ipar(5)
c$$$      endif
c$$$      idx = n * (m+1)
c$$$      iz = idx + n
c$$$      hess = iz + n*m
c$$$      vc = hess + (m+1) * m / 2 + 1
c$$$      vs = vc + m
c$$$      vrn = vs + m
c$$$      i = vrn + m + 1
c$$$      call zbisinit(ipar,fpar,i,1,lp,rp,w)
c$$$      if (ipar(1).lt.0) return
c$$$c
c$$$c     request for matrix vector multiplication A*x in the initialization
c$$$c
c$$$ 100  ipar(1) = 1
c$$$      ipar(8) = n+1
c$$$      ipar(9) = 1
c$$$      ipar(10) = 1
c$$$      k = 0
c$$$      do ii = 1, n
c$$$         w(ii+n) = sol(ii)
c$$$      enddo
c$$$      return
c$$$ 10   ipar(7) = ipar(7) + 1
c$$$      ipar(13) = ipar(13) + 1
c$$$      fpar(11) = fpar(11) + n
c$$$      if (lp) then
c$$$         do i = 1, n
c$$$            w(n+i) = rhs(i) - w(i)
c$$$         enddo
c$$$         ipar(1) = 3
c$$$         ipar(10) = 2
c$$$         return
c$$$      else
c$$$         do i = 1, n
c$$$            w(i) = rhs(i) - w(i)
c$$$         enddo
c$$$      endif
c$$$c
c$$$ 20   alpha = sqrt(zdistdot(n,w,1,w,1))
c$$$      fpar(11) = fpar(11) + n + n
c$$$      if (ipar(7).eq.1 .and. ipar(3).ne.999) then
c$$$         if (abs(ipar(3)).eq.2) then
c$$$            fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
c$$$            fpar(11) = fpar(11) + 2*n
c$$$         else
c$$$            fpar(4) = fpar(1) * alpha + fpar(2)
c$$$         endif
c$$$         fpar(3) = alpha
c$$$      endif
c$$$      fpar(5) = alpha
c$$$      w(vrn+1) = alpha
c$$$      if (alpha.le.fpar(4) .and. ipar(3).ge.0 .and. ipar(3).ne.999) then
c$$$         ipar(1) = 0
c$$$         fpar(6) = alpha
c$$$         goto 300
c$$$      endif
c$$$      alpha = one / alpha
c$$$      do ii = 1, n
c$$$         w(ii) = w(ii) * alpha
c$$$      enddo
c$$$      fpar(11) = fpar(11) + n
c$$$c
c$$$c     request for (1) right preconditioning
c$$$c     (2) matrix vector multiplication
c$$$c     (3) left preconditioning
c$$$c
c$$$ 110  k = k + 1
c$$$      if (rp) then
c$$$         ipar(1) = 5
c$$$         ipar(8) = k*n - n + 1
c$$$         ipar(9) = iz + ipar(8)
c$$$         ipar(10) = 3
c$$$         return
c$$$      else
c$$$         do ii = 0, n-1
c$$$            w(iz+k*n-ii) = w(k*n-ii)
c$$$         enddo
c$$$      endif
c$$$c
c$$$ 30   ipar(1) = 1
c$$$      if (rp) then
c$$$         ipar(8) = ipar(9)
c$$$      else
c$$$         ipar(8) = (k-1)*n + 1
c$$$      endif
c$$$      if (lp) then
c$$$         ipar(9) = idx + 1
c$$$      else
c$$$         ipar(9) = 1 + k*n
c$$$      endif
c$$$      ipar(10) = 4
c$$$      return
c$$$c
c$$$ 40   if (lp) then
c$$$         ipar(1) = 3
c$$$         ipar(8) = ipar(9)
c$$$         ipar(9) = k*n + 1
c$$$         ipar(10) = 5
c$$$         return
c$$$      endif
c$$$c
c$$$c     Modified Gram-Schmidt orthogonalization procedure
c$$$c     temporary pointer 'ptr' is pointing to the current column of the
c$$$c     Hessenberg matrix. 'p2' points to the new basis vector
c$$$c
c$$$ 50   ptr = k * (k - 1) / 2 + hess
c$$$      p2 = ipar(9)
c$$$      ipar(7) = ipar(7) + 1
c$$$      call zmgsro(.false.,n,n,k+1,k+1,fpar(11),w,w(ptr+1),ipar(12))
c$$$      if (ipar(12).lt.0) goto 200
c$$$c
c$$$c     apply previous Givens rotations and generate a new one to eliminate
c$$$c     the subdiagonal element.
c$$$c
c$$$      p2 = ptr + 1
c$$$      do i = 1, k-1
c$$$         ptr = p2
c$$$         p2 = p2 + 1
c$$$         alpha = w(ptr)
c$$$         c = w(vc+i)
c$$$         s = w(vs+i)
c$$$         w(ptr) = c * alpha + s * w(p2)
c$$$         w(p2) = c * w(p2) - s * alpha
c$$$      enddo
c$$$      call zgivens(w(p2), w(p2+1), c, s)
c$$$      w(vc+k) = c
c$$$      w(vs+k) = s
c$$$      p2 = vrn + k
c$$$      alpha = - s * w(p2)
c$$$      w(p2) = c * w(p2)
c$$$      w(p2+1) = alpha
c$$$      fpar(11) = fpar(11) + 6 * k
c$$$c
c$$$c     end of one Arnoldi iteration, alpha will store the estimated
c$$$c     residual norm at current stage
c$$$c
c$$$      alpha = abs(alpha)
c$$$      fpar(5) = alpha
c$$$      if (k.lt.m .and. .not.(ipar(3).ge.0 .and. alpha.le.fpar(4))
c$$$     +      .and. (ipar(6).le.0 .or. ipar(7).lt.ipar(6))) goto 110
c$$$c
c$$$c     update the approximate solution, first solve the upper triangular
c$$$c     system, temporary pointer ptr points to the Hessenberg matrix,
c$$$c     p2 points to the right-hand-side (also the solution) of the system.
c$$$c
c$$$ 200  ptr = hess + k * (k + 1 ) / 2
c$$$      p2 = vrn + k
c$$$      if (w(ptr).eq.zero) then
c$$$c
c$$$c     if the diagonal elements of the last column is zero, reduce k by 1
c$$$c     so that a smaller trianguler system is solved [It should only
c$$$c     happen when the matrix is singular!]
c$$$c
c$$$         k = k - 1
c$$$         if (k.gt.0) then
c$$$            goto 200
c$$$         else
c$$$            ipar(1) = -3
c$$$            ipar(12) = -4
c$$$            goto 300
c$$$         endif
c$$$      endif
c$$$      w(p2) = w(p2) / w(ptr)
c$$$      do i = k-1, 1, -1
c$$$         ptr = ptr - i - 1
c$$$         do ii = 1, i
c$$$            w(vrn+ii) = w(vrn+ii) - w(p2) * w(ptr+ii)
c$$$         enddo
c$$$         p2 = p2 - 1
c$$$         w(p2) = w(p2) / w(ptr)
c$$$      enddo
c$$$c
c$$$      do i = 0, k-1
c$$$         ptr = iz+i*n
c$$$         do ii = 1, n
c$$$            sol(ii) = sol(ii) + w(p2)*w(ptr+ii)
c$$$         enddo
c$$$         p2 = p2 + 1
c$$$      enddo
c$$$      fpar(11) = fpar(11) + 2*k*n + k*(k+1)
c$$$c
c$$$c     process the complete stopping criteria
c$$$c
c$$$      if (ipar(3).eq.999) then
c$$$         ipar(1) = 10
c$$$         ipar(8) = -1
c$$$         ipar(9) = idx + 1
c$$$         ipar(10) = 6
c$$$         return
c$$$      else if (ipar(3).lt.0) then
c$$$         if (ipar(7).le.m+1) then
c$$$            fpar(3) = abs(w(vrn+1))
c$$$            if (ipar(3).eq.-1) fpar(4) = fpar(1)*fpar(3)+fpar(2)
c$$$         endif
c$$$         fpar(6) = abs(w(vrn+k))
c$$$      else if (ipar(3).ne.999) then
c$$$         fpar(6) = fpar(5)
c$$$      endif
c$$$c
c$$$c     do we need to restart ?
c$$$c
c$$$ 60   if (ipar(12).ne.0) then
c$$$         ipar(1) = -3
c$$$         goto 300
c$$$      endif
c$$$      if ((ipar(7).lt.ipar(6) .or. ipar(6).le.0).and.
c$$$     +     ((ipar(3).eq.999.and.ipar(11).eq.0) .or.
c$$$     +     (ipar(3).ne.999.and.fpar(6).gt.fpar(4)))) goto 100
c$$$c
c$$$c     termination, set error code, compute convergence rate
c$$$c
c$$$      if (ipar(1).gt.0) then
c$$$         if (ipar(3).eq.999 .and. ipar(11).eq.1) then
c$$$            ipar(1) = 0
c$$$         else if (ipar(3).ne.999 .and. fpar(6).le.fpar(4)) then
c$$$            ipar(1) = 0
c$$$         else if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
c$$$            ipar(1) = -1
c$$$         else
c$$$            ipar(1) = -10
c$$$         endif
c$$$      endif
c$$$ 300  if (fpar(3).ne.zero .and. fpar(6).ne.zero .and.
c$$$     $     ipar(7).gt.ipar(13)) then
c$$$         fpar(7) = log10(fpar(3) / fpar(6)) / dble(ipar(7)-ipar(13))
c$$$      else
c$$$         fpar(7) = zero
c$$$      endif
c$$$      return
c$$$      end
c$$$c-----end-of-fgmres
c$$$c-----------------------------------------------------------------------
c$$$      subroutine zdbcg (n,rhs,sol,ipar,fpar,w)
c$$$      implicit none
c$$$      integer n,ipar(16)
c$$$      real*8 rhs(n), sol(n), fpar(16), w(n,*)
c$$$c-----------------------------------------------------------------------
c$$$c Quasi GMRES method for solving a linear
c$$$c system of equations a * sol = y.  double precision version.
c$$$c this version is without restarting and without preconditioning.
c$$$c parameters :
c$$$c -----------
c$$$c n     = dimension of the problem
c$$$c
c$$$c y     = w(:,1) a temporary storage used for various operations
c$$$c z     = w(:,2) a work vector of length n.
c$$$c v     = w(:,3:4) size n x 2
c$$$c w     = w(:,5:6) size n x 2
c$$$c p     = w(:,7:9) work array of dimension n x 3
c$$$c del x = w(:,10)  accumulation of the changes in solution
c$$$c tmp   = w(:,11)  a temporary vector used to hold intermediate result of
c$$$c                  preconditioning, etc.
c$$$c
c$$$c sol   = the solution of the problem . at input sol must contain an
c$$$c         initial guess to the solution.
c$$$c    ***  note:   y is destroyed on return.
c$$$c
c$$$c-----------------------------------------------------------------------
c$$$c subroutines and functions called:
c$$$c 1) matrix vector multiplication and preconditioning through reverse
c$$$c     communication
c$$$c
c$$$c 2) implu, uppdir, distdot (blas)
c$$$c-----------------------------------------------------------------------
c$$$c aug. 1983  version.    author youcef saad. yale university computer
c$$$c science dept. some  changes made july 3, 1986.
c$$$c references: siam j. sci. stat. comp., vol. 5, pp. 203-228 (1984)
c$$$c-----------------------------------------------------------------------
c$$$c     local variables
c$$$c
c$$$      real*8 one,zero
c$$$      parameter(one=1.0D0,zero=0.0D0)
c$$$c
c$$$      real*8 t,sqrt,distdot,ss,res,beta,ss1,delta,x,zeta,umm
c$$$      integer k,j,i,i2,ip2,ju,lb,lbm1,np,indp
c$$$      logical lp,rp,full, perm(3)
c$$$      real*8 ypiv(3),u(3),usav(3)
c$$$      external tidycg
c$$$      save
c$$$c
c$$$c     where to go
c$$$c
c$$$      if (ipar(1).le.0) ipar(10) = 0
c$$$      goto (110, 120, 130, 140, 150, 160, 170, 180, 190, 200) ipar(10)
c$$$c
c$$$c     initialization, parameter checking, clear the work arrays
c$$$c
c$$$      call zbisinit(ipar,fpar,11*n,1,lp,rp,w)
c$$$      if (ipar(1).lt.0) return
c$$$      perm(1) = .false.
c$$$      perm(2) = .false.
c$$$      perm(3) = .false.
c$$$      usav(1) = zero
c$$$      usav(2) = zero
c$$$      usav(3) = zero
c$$$      ypiv(1) = zero
c$$$      ypiv(2) = zero
c$$$      ypiv(3) = zero
c$$$c-----------------------------------------------------------------------
c$$$c     initialize constants for outer loop :
c$$$c-----------------------------------------------------------------------
c$$$      lb = 3
c$$$      lbm1 = 2
c$$$c
c$$$c     get initial residual vector and norm
c$$$c
c$$$      ipar(1) = 1
c$$$      ipar(8) = 1
c$$$      ipar(9) = 1 + n
c$$$      do i = 1, n
c$$$         w(i,1) = sol(i)
c$$$      enddo
c$$$      ipar(10) = 1
c$$$      return
c$$$ 110  ipar(7) = ipar(7) + 1
c$$$      ipar(13) = ipar(13) + 1
c$$$      if (lp) then
c$$$         do i = 1, n
c$$$            w(i,1) = rhs(i) - w(i,2)
c$$$         enddo
c$$$         ipar(1) = 3
c$$$         ipar(8) = 1
c$$$         ipar(9) = n+n+1
c$$$         ipar(10) = 2
c$$$         return
c$$$      else
c$$$         do i = 1, n
c$$$            w(i,3) = rhs(i) - w(i,2)
c$$$         enddo
c$$$      endif
c$$$      fpar(11) = fpar(11) + n
c$$$c
c$$$ 120  fpar(3) = sqrt(distdot(n,w(1,3),1,w(1,3),1))
c$$$      fpar(11) = fpar(11) + n + n
c$$$      fpar(5) = fpar(3)
c$$$      fpar(7) = fpar(3)
c$$$      zeta = fpar(3)
c$$$      if (abs(ipar(3)).eq.2) then
c$$$         fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
c$$$         fpar(11) = fpar(11) + 2*n
c$$$      else if (ipar(3).ne.999) then
c$$$         fpar(4) = fpar(1) * zeta + fpar(2)
c$$$      endif
c$$$      if (ipar(3).ge.0.and.fpar(5).le.fpar(4)) then
c$$$         fpar(6) = fpar(5)
c$$$         goto 900
c$$$      endif
c$$$c
c$$$c     normalize first arnoldi vector
c$$$c
c$$$      t = one/zeta
c$$$      do 22 k=1,n
c$$$         w(k,3) = w(k,3)*t
c$$$         w(k,5) = w(k,3)
c$$$ 22   continue
c$$$      fpar(11) = fpar(11) + n
c$$$c
c$$$c     initialize constants for main loop
c$$$c
c$$$      beta = zero
c$$$      delta = zero
c$$$      i2 = 1
c$$$      indp = 0
c$$$      i = 0
c$$$c
c$$$c     main loop: i = index of the loop.
c$$$c
c$$$c-----------------------------------------------------------------------
c$$$ 30   i = i + 1
c$$$c
c$$$      if (rp) then
c$$$         ipar(1) = 5
c$$$         ipar(8) = (1+i2)*n+1
c$$$         if (lp) then
c$$$            ipar(9) = 1
c$$$         else
c$$$            ipar(9) = 10*n + 1
c$$$         endif
c$$$         ipar(10) = 3
c$$$         return
c$$$      endif
c$$$c
c$$$ 130  ipar(1) = 1
c$$$      if (rp) then
c$$$         ipar(8) = ipar(9)
c$$$      else
c$$$         ipar(8) = (1+i2)*n + 1
c$$$      endif
c$$$      if (lp) then
c$$$         ipar(9) = 10*n + 1
c$$$      else
c$$$         ipar(9) = 1
c$$$      endif
c$$$      ipar(10) = 4
c$$$      return
c$$$c
c$$$ 140  if (lp) then
c$$$         ipar(1) = 3
c$$$         ipar(8) = ipar(9)
c$$$         ipar(9) = 1
c$$$         ipar(10) = 5
c$$$         return
c$$$      endif
c$$$c
c$$$c     A^t * x
c$$$c
c$$$ 150  ipar(7) = ipar(7) + 1
c$$$      if (lp) then
c$$$         ipar(1) = 4
c$$$         ipar(8) = (3+i2)*n + 1
c$$$         if (rp) then
c$$$            ipar(9) = n + 1
c$$$         else
c$$$            ipar(9) = 10*n + 1
c$$$         endif
c$$$         ipar(10) = 6
c$$$         return
c$$$      endif
c$$$c
c$$$ 160  ipar(1) = 2
c$$$      if (lp) then
c$$$         ipar(8) = ipar(9)
c$$$      else
c$$$         ipar(8) = (3+i2)*n + 1
c$$$      endif
c$$$      if (rp) then
c$$$         ipar(9) = 10*n + 1
c$$$      else
c$$$         ipar(9) = n + 1
c$$$      endif
c$$$      ipar(10) = 7
c$$$      return
c$$$c
c$$$ 170  if (rp) then
c$$$         ipar(1) = 6
c$$$         ipar(8) = ipar(9)
c$$$         ipar(9) = n + 1
c$$$         ipar(10) = 8
c$$$         return
c$$$      endif
c$$$c-----------------------------------------------------------------------
c$$$c     orthogonalize current v against previous v's and
c$$$c     determine relevant part of i-th column of u(.,.) the
c$$$c     upper triangular matrix --
c$$$c-----------------------------------------------------------------------
c$$$ 180  ipar(7) = ipar(7) + 1
c$$$      u(1) = zero
c$$$      ju = 1
c$$$      k = i2
c$$$      if (i .le. lbm1) ju = 0
c$$$      if (i .lt. lb) k = 0
c$$$ 31   if (k .eq. lbm1) k=0
c$$$      k=k+1
c$$$c
c$$$      if (k .ne. i2) then
c$$$         ss  = delta
c$$$         ss1 = beta
c$$$         ju = ju + 1
c$$$         u(ju) = ss
c$$$      else
c$$$         ss = distdot(n,w(1,1),1,w(1,4+k),1)
c$$$         fpar(11) = fpar(11) + 2*n
c$$$         ss1= ss
c$$$         ju = ju + 1
c$$$         u(ju) = ss
c$$$      endif
c$$$c
c$$$      do 32  j=1,n
c$$$         w(j,1) = w(j,1) - ss*w(j,k+2)
c$$$         w(j,2) = w(j,2) - ss1*w(j,k+4)
c$$$ 32   continue
c$$$      fpar(11) = fpar(11) + 4*n
c$$$c
c$$$      if (k .ne. i2) goto 31
c$$$c
c$$$c     end of Mod. Gram. Schmidt loop
c$$$c
c$$$      t = distdot(n,w(1,2),1,w(1,1),1)
c$$$c
c$$$      beta   = sqrt(abs(t))
c$$$      delta  = t/beta
c$$$c
c$$$      ss = one/beta
c$$$      ss1 = one/ delta
c$$$c
c$$$c     normalize and insert new vectors
c$$$c
c$$$      ip2 = i2
c$$$      if (i2 .eq. lbm1) i2=0
c$$$      i2=i2+1
c$$$c
c$$$      do 315 j=1,n
c$$$         w(j,i2+2)=w(j,1)*ss
c$$$         w(j,i2+4)=w(j,2)*ss1
c$$$ 315  continue
c$$$      fpar(11) = fpar(11) + 4*n
c$$$c-----------------------------------------------------------------------
c$$$c     end of orthogonalization.
c$$$c     now compute the coefficients u(k) of the last
c$$$c     column of the  l . u  factorization of h .
c$$$c-----------------------------------------------------------------------
c$$$      np = min0(i,lb)
c$$$      full = (i .ge. lb)
c$$$      call zimplu(np, umm, beta, ypiv, u, perm, full)
c$$$c-----------------------------------------------------------------------
c$$$c     update conjugate directions and solution
c$$$c-----------------------------------------------------------------------
c$$$      do 33 k=1,n
c$$$         w(k,1) = w(k,ip2+2)
c$$$ 33   continue
c$$$      call zuppdir(n, w(1,7), np, lb, indp, w, u, usav, fpar(11))
c$$$c-----------------------------------------------------------------------
c$$$      if (i .eq. 1) goto 34
c$$$      j = np - 1
c$$$      if (full) j = j-1
c$$$      if (.not.perm(j)) zeta = -zeta*ypiv(j)
c$$$ 34   x = zeta/u(np)
c$$$      if (perm(np))goto 36
c$$$      do 35 k=1,n
c$$$         w(k,10) = w(k,10) + x*w(k,1)
c$$$ 35   continue
c$$$      fpar(11) = fpar(11) + 2 * n
c$$$c-----------------------------------------------------------------------
c$$$ 36   if (ipar(3).eq.999) then
c$$$         ipar(1) = 10
c$$$         ipar(8) = 9*n + 1
c$$$         ipar(9) = 10*n + 1
c$$$         ipar(10) = 9
c$$$         return
c$$$      endif
c$$$      res = abs(beta*zeta/umm)
c$$$      fpar(5) = res * sqrt(distdot(n, w(1,i2+2), 1, w(1,i2+2), 1))
c$$$      fpar(11) = fpar(11) + 2 * n
c$$$      if (ipar(3).lt.0) then
c$$$         fpar(6) = x * sqrt(distdot(n,w,1,w,1))
c$$$         fpar(11) = fpar(11) + 2 * n
c$$$         if (ipar(7).le.3) then
c$$$            fpar(3) = fpar(6)
c$$$            if (ipar(3).eq.-1) then
c$$$               fpar(4) = fpar(1) * sqrt(fpar(3)) + fpar(2)
c$$$            endif
c$$$         endif
c$$$      else
c$$$         fpar(6) = fpar(5)
c$$$      endif
c$$$c---- convergence test -----------------------------------------------
c$$$ 190  if (ipar(3).eq.999.and.ipar(11).eq.0) then
c$$$         goto 30
c$$$      else if (fpar(6).gt.fpar(4) .and. (ipar(6).gt.ipar(7) .or.
c$$$     +        ipar(6).le.0)) then
c$$$         goto 30
c$$$      endif
c$$$c-----------------------------------------------------------------------
c$$$c     here the fact that the last step is different is accounted for.
c$$$c-----------------------------------------------------------------------
c$$$      if (.not. perm(np)) goto 900
c$$$      x = zeta/umm
c$$$      do 40 k = 1,n
c$$$         w(k,10) = w(k,10) + x*w(k,1)
c$$$ 40   continue
c$$$      fpar(11) = fpar(11) + 2 * n
c$$$c
c$$$c     right preconditioning and clean-up jobs
c$$$c
c$$$ 900  if (rp) then
c$$$         if (ipar(1).lt.0) ipar(12) = ipar(1)
c$$$         ipar(1) = 5
c$$$         ipar(8) = 9*n + 1
c$$$         ipar(9) = ipar(8) + n
c$$$         ipar(10) = 10
c$$$         return
c$$$      endif
c$$$ 200  if (rp) then
c$$$         call ztidycg(n,ipar,fpar,sol,w(1,11))
c$$$      else
c$$$         call ztidycg(n,ipar,fpar,sol,w(1,10))
c$$$      endif
c$$$      return
c$$$      end
c$$$c-----end-of-dbcg-------------------------------------------------------
c$$$c-----------------------------------------------------------------------
c$$$      subroutine zimplu(np,umm,beta,ypiv,u,permut,full)
c$$$      real*8 umm,beta,ypiv(*),u(*),x, xpiv
c$$$      logical full, perm, permut(*)
c$$$      integer np,k,npm1
c$$$c-----------------------------------------------------------------------
c$$$c     performs implicitly one step of the lu factorization of a
c$$$c     banded hessenberg matrix.
c$$$c-----------------------------------------------------------------------
c$$$      if (np .le. 1) goto 12
c$$$      npm1 = np - 1
c$$$c
c$$$c     -- perform  previous step of the factorization-
c$$$c
c$$$      do 6 k=1,npm1
c$$$         if (.not. permut(k)) goto 5
c$$$         x=u(k)
c$$$         u(k) = u(k+1)
c$$$         u(k+1) = x
c$$$ 5       u(k+1) = u(k+1) - ypiv(k)*u(k)
c$$$ 6    continue
c$$$c-----------------------------------------------------------------------
c$$$c     now determine pivotal information to be used in the next call
c$$$c-----------------------------------------------------------------------
c$$$ 12   umm = u(np)
c$$$      perm = (beta .gt. abs(umm))
c$$$      if (.not. perm) goto 4
c$$$      xpiv = umm / beta
c$$$      u(np) = beta
c$$$      goto 8
c$$$ 4    xpiv = beta/umm
c$$$ 8    permut(np) = perm
c$$$      ypiv(np) = xpiv
c$$$      if (.not. full) return
c$$$c     shift everything up if full...
c$$$      do 7 k=1,npm1
c$$$         ypiv(k) = ypiv(k+1)
c$$$         permut(k) = permut(k+1)
c$$$ 7    continue
c$$$      return
c$$$c-----end-of-implu
c$$$      end
c$$$c-----------------------------------------------------------------------
c$$$      subroutine zuppdir(n,p,np,lbp,indp,y,u,usav,flops)
c$$$      real*8 p(n,lbp), y(*), u(*), usav(*), x, flops
c$$$      integer k,np,n,npm1,j,ju,indp,lbp
c$$$c-----------------------------------------------------------------------
c$$$c     updates the conjugate directions p given the upper part of the
c$$$c     banded upper triangular matrix u.  u contains the non zero
c$$$c     elements of the column of the triangular matrix..
c$$$c-----------------------------------------------------------------------
c$$$      real*8 zero
c$$$      parameter(zero=0.0D0)
c$$$c
c$$$      npm1=np-1
c$$$      if (np .le. 1) goto 12
c$$$      j=indp
c$$$      ju = npm1
c$$$ 10   if (j .le. 0) j=lbp
c$$$      x = u(ju) /usav(j)
c$$$      if (x .eq. zero) goto 115
c$$$      do 11 k=1,n
c$$$         y(k) = y(k) - x*p(k,j)
c$$$ 11   continue
c$$$      flops = flops + 2*n
c$$$ 115  j = j-1
c$$$      ju = ju -1
c$$$      if (ju .ge. 1) goto 10
c$$$ 12   indp = indp + 1
c$$$      if (indp .gt. lbp) indp = 1
c$$$      usav(indp) = u(np)
c$$$      do 13 k=1,n
c$$$         p(k,indp) = y(k)
c$$$ 13   continue
c$$$ 208  return
c$$$c-----------------------------------------------------------------------
c$$$c-------end-of-uppdir---------------------------------------------------
c$$$      end
c$$$      subroutine zgivens(x,y,c,s)
c$$$      real*8 x,y,c,s
c$$$c-----------------------------------------------------------------------
c$$$c     Given x and y, this subroutine generates a Givens' rotation c, s.
c$$$c     And apply the rotation on (x,y) ==> (sqrt(x**2 + y**2), 0).
c$$$c     (See P 202 of "matrix computation" by Golub and van Loan.)
c$$$c-----------------------------------------------------------------------
c$$$      real*8 t,one,zero
c$$$      parameter (zero=0.0D0,one=1.0D0)
c$$$c
c$$$      if (x.eq.zero .and. y.eq.zero) then
c$$$         c = one
c$$$         s = zero
c$$$      else if (abs(y).gt.abs(x)) then
c$$$         t = x / y
c$$$         x = sqrt(one+t*t)
c$$$         s = sign(one / x, y)
c$$$         c = t*s
c$$$      else if (abs(y).le.abs(x)) then
c$$$         t = y / x
c$$$         y = sqrt(one+t*t)
c$$$         c = sign(one / y, x)
c$$$         s = t*c
c$$$      else
c$$$c
c$$$c     X or Y must be an invalid floating-point number, set both to zero
c$$$c
c$$$         x = zero
c$$$         y = zero
c$$$         c = one
c$$$         s = zero
c$$$      endif
c$$$      x = abs(x*y)
c$$$c
c$$$c     end of givens
c$$$c
c$$$      return
c$$$      end
c$$$c-----end-of-givens
c$$$c-----------------------------------------------------------------------
c$$$      logical function zstopbis(n,ipar,mvpi,fpar,r,delx,sx)
c$$$      implicit none
c$$$      integer n,mvpi,ipar(16)
c$$$      complex*16 fpar(16), r(n), delx(n), sx, distdot
c$$$      external distdot
c$$$c-----------------------------------------------------------------------
c$$$c     function for determining the stopping criteria. return value of
c$$$c     true if the stopbis criteria is satisfied.
c$$$c-----------------------------------------------------------------------
c$$$      if (ipar(11) .eq. 1) then
c$$$         stopbis = .true.
c$$$      else
c$$$         stopbis = .false.
c$$$      endif
c$$$      if (ipar(6).gt.0 .and. ipar(7).ge.ipar(6)) then
c$$$         ipar(1) = -1
c$$$         stopbis = .true.
c$$$      endif
c$$$      if (stopbis) return
c$$$c
c$$$c     computes errors
c$$$c
c$$$      fpar(5) = dznrm2(n,r,1) 
c$$$      fpar(11) = fpar(11) + 2 * n
c$$$      if (ipar(3).lt.0) then
c$$$c
c$$$c     compute the change in the solution vector
c$$$c
c$$$         fpar(6) = sx * dznrm2(n,delx,1)
c$$$         fpar(11) = fpar(11) + 2 * n
c$$$         if (ipar(7).lt.mvpi+mvpi+1) then
c$$$c
c$$$c     if this is the end of the first iteration, set fpar(3:4)
c$$$c
c$$$            fpar(3) = fpar(6)
c$$$            if (ipar(3).eq.-1) then
c$$$               fpar(4) = fpar(1) * fpar(3) + fpar(2)
c$$$            endif
c$$$         endif
c$$$      else
c$$$         fpar(6) = fpar(5)
c$$$      endif
c$$$c
c$$$c     .. the test is struct this way so that when the value in fpar(6)
c$$$c       is not a valid number, STOPBIS is set to .true.
c$$$c
c$$$      if (fpar(6).gt.fpar(4)) then
c$$$         stopbis = .false.
c$$$         ipar(11) = 0
c$$$      else
c$$$         stopbis = .true.
c$$$         ipar(11) = 1
c$$$      endif
c$$$c
c$$$      return
c$$$      end
c$$$c-----end-of-stopbis
c$$$c-----------------------------------------------------------------------
      subroutine ztidycg(n,ipar,fpar,sol,delx)
      implicit none
      integer i,n,ipar(16)
      real*8 fpar(16)
      complex*16 sol(n),delx(n)
c-----------------------------------------------------------------------
c     Some common operations required before terminating the CG routines
c-----------------------------------------------------------------------
      real*8 zero
      parameter(zero=0.0D0)
c
      if (ipar(12).ne.0) then
         ipar(1) = ipar(12) 
      else if (ipar(1).gt.0) then
         if ((ipar(3).eq.999 .and. ipar(11).eq.1) .or.
     +        fpar(6).le.fpar(4)) then
            ipar(1) = 0
         else if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
            ipar(1) = -1
         else
            ipar(1) = -10
         endif
      endif
      if (fpar(3).gt.zero .and. fpar(6).gt.zero .and.
     +     ipar(7).gt.ipar(13)) then
         fpar(7) = log10(fpar(3) / fpar(6)) / dble(ipar(7)-ipar(13))
      else
         fpar(7) = zero
      endif
      do i = 1, n
         sol(i) = sol(i) + delx(i)
      enddo
      return
      end
c-----end-of-tidycg
c-----------------------------------------------------------------------
      logical function zbrkdn(alpha, ipar)
      implicit none
      integer ipar(16)
      real*8 alpha
      real*8 zero, one, beta
      parameter (zero=0.0D0, one=1.0D0)
c-----------------------------------------------------------------------
c     test whether alpha is zero or an abnormal number, if yes,
c     this routine will return .true.
c
c     If alpha == 0, ipar(1) = -3,
c     if alpha is an abnormal number, ipar(1) = -9.
c-----------------------------------------------------------------------
      zbrkdn = .false.
      if (abs(alpha).gt.zero) then
         beta = one / abs(alpha)
         if (.not. beta.gt.zero) then
            zbrkdn = .true.
            ipar(1) = -9
         endif
      else if (abs(alpha).lt.zero) then
         beta = one / abs(alpha)
         if (.not. beta.lt.zero) then
            zbrkdn = .true.
            ipar(1) = -9
         endif
      else if (abs(alpha).eq.zero) then
         zbrkdn = .true.
         ipar(1) = -3
      else
         zbrkdn = .true.
         ipar(1) = -9
      endif
      return
      end
c-----end-of-zbrkdn
c-----------------------------------------------------------------------
      subroutine zbisinit(ipar,fpar,wksize,dsc,lp,rp,wk)
      implicit none
      integer i,ipar(16),wksize,dsc
      logical lp,rp
      real*8  fpar(16)
      complex*16 wk(*)
c-----------------------------------------------------------------------
c     some common initializations for the iterative solvers
c-----------------------------------------------------------------------
      real*8 zero, one
      parameter(zero=0.0D0, one=1.0D0)
c
c     ipar(1) = -2 inidcate that there are not enough space in the work
c     array
c
      if (ipar(4).lt.wksize) then
         ipar(1) = -2
         ipar(4) = wksize
         return
      endif
c
      if (ipar(2).gt.2) then
         lp = .true.
         rp = .true.
      else if (ipar(2).eq.2) then
         lp = .false.
         rp = .true.
      else if (ipar(2).eq.1) then
         lp = .true.
         rp = .false.
      else
         lp = .false.
         rp = .false.
      endif
      if (ipar(3).eq.0) ipar(3) = dsc
c     .. clear the ipar elements used
      ipar(7) = 0
      ipar(8) = 0
      ipar(9) = 0
      ipar(10) = 0
      ipar(11) = 0
      ipar(12) = 0
      ipar(13) = 0
c
c     fpar(1) must be between (0, 1), fpar(2) must be positive,
c     fpar(1) and fpar(2) can NOT both be zero
c     Normally return ipar(1) = -4 to indicate any of above error
c
      if (fpar(1).lt.zero .or. fpar(1).ge.one .or. fpar(2).lt.zero .or.
     &     (fpar(1).eq.zero .and. fpar(2).eq.zero)) then
         if (ipar(1).eq.0) then
            ipar(1) = -4
            return
         else
            fpar(1) = 1.0D-6
            fpar(2) = 1.0D-16
         endif
      endif
c     .. clear the fpar elements
      do i = 3, 10
         fpar(i) = zero
      enddo
      if (fpar(11).lt.zero) fpar(11) = zero
c     .. clear the used portion of the work array to zero
      do i = 1, wksize
         wk(i) = (zero,zero)
      enddo
c
      return
c-----end-of-bisinit
      end
c$$$c-----------------------------------------------------------------------
c$$$      subroutine zmgsro(full,lda,n,m,ind,ops,vec,hh,ierr)
c$$$      implicit none
c$$$      logical full
c$$$      integer lda,m,n,ind,ierr
c$$$      real*8  ops,hh(m),vec(lda,m)
c$$$c-----------------------------------------------------------------------
c$$$c     MGSRO  -- Modified Gram-Schmidt procedure with Selective Re-
c$$$c               Orthogonalization
c$$$c     The ind'th vector of VEC is orthogonalized against the rest of
c$$$c     the vectors.
c$$$c
c$$$c     The test for performing re-orthogonalization is performed for
c$$$c     each indivadual vectors. If the cosine between the two vectors
c$$$c     is greater than 0.99 (REORTH = 0.99**2), re-orthogonalization is
c$$$c     performed. The norm of the 'new' vector is kept in variable NRM0,
c$$$c     and updated after operating with each vector.
c$$$c
c$$$c     full   -- .ture. if it is necessary to orthogonalize the ind'th
c$$$c               against all the vectors vec(:,1:ind-1), vec(:,ind+2:m)
c$$$c               .false. only orthogonalize againt vec(:,1:ind-1)
c$$$c     lda    -- the leading dimension of VEC
c$$$c     n      -- length of the vector in VEC
c$$$c     m      -- number of vectors can be stored in VEC
c$$$c     ind    -- index to the vector to be changed
c$$$c     ops    -- operation counts
c$$$c     vec    -- vector of LDA X M storing the vectors
c$$$c     hh     -- coefficient of the orthogonalization
c$$$c     ierr   -- error code
c$$$c               0 : successful return
c$$$c               -1: zero input vector
c$$$c               -2: input vector contains abnormal numbers
c$$$c               -3: input vector is a linear combination of others
c$$$c
c$$$c     External routines used: real*8 distdot
c$$$c-----------------------------------------------------------------------
c$$$      integer i,k
c$$$      real*8  nrm0, nrm1, fct, thr, distdot, zero, one, reorth
c$$$      parameter (zero=0.0D0, one=1.0D0, reorth=0.98D0)
c$$$      external distdot
c$$$c
c$$$c     compute the norm of the input vector
c$$$c
c$$$      nrm0 = distdot(n,vec(1,ind),1,vec(1,ind),1)
c$$$      ops = ops + n + n
c$$$      thr = nrm0 * reorth
c$$$      if (nrm0.le.zero) then
c$$$         ierr = - 1
c$$$         return
c$$$      else if (nrm0.gt.zero .and. one/nrm0.gt.zero) then
c$$$         ierr = 0
c$$$      else
c$$$         ierr = -2
c$$$         return
c$$$      endif
c$$$c
c$$$c     Modified Gram-Schmidt loop
c$$$c
c$$$      if (full) then
c$$$         do 40 i = ind+1, m
c$$$            fct = distdot(n,vec(1,ind),1,vec(1,i),1)
c$$$            hh(i) = fct
c$$$            do 20 k = 1, n
c$$$               vec(k,ind) = vec(k,ind) - fct * vec(k,i)
c$$$ 20         continue
c$$$            ops = ops + 4 * n + 2
c$$$            if (fct*fct.gt.thr) then
c$$$               fct = distdot(n,vec(1,ind),1,vec(1,i),1)
c$$$               hh(i) = hh(i) + fct
c$$$               do 30 k = 1, n
c$$$                  vec(k,ind) = vec(k,ind) - fct * vec(k,i)
c$$$ 30            continue
c$$$               ops = ops + 4*n + 1
c$$$            endif
c$$$            nrm0 = nrm0 - hh(i) * hh(i)
c$$$            if (nrm0.lt.zero) nrm0 = zero
c$$$            thr = nrm0 * reorth
c$$$ 40      continue
c$$$      endif
c$$$c
c$$$      do 70 i = 1, ind-1
c$$$         fct = distdot(n,vec(1,ind),1,vec(1,i),1)
c$$$         hh(i) = fct
c$$$         do 50 k = 1, n
c$$$            vec(k,ind) = vec(k,ind) - fct * vec(k,i)
c$$$ 50      continue
c$$$         ops = ops + 4 * n + 2
c$$$         if (fct*fct.gt.thr) then
c$$$            fct = distdot(n,vec(1,ind),1,vec(1,i),1)
c$$$            hh(i) = hh(i) + fct
c$$$            do 60 k = 1, n
c$$$               vec(k,ind) = vec(k,ind) - fct * vec(k,i)
c$$$ 60         continue
c$$$            ops = ops + 4*n + 1
c$$$         endif
c$$$         nrm0 = nrm0 - hh(i) * hh(i)
c$$$         if (nrm0.lt.zero) nrm0 = zero
c$$$         thr = nrm0 * reorth
c$$$ 70   continue
c$$$c
c$$$c     test the resulting vector
c$$$c
c$$$      nrm1 = sqrt(distdot(n,vec(1,ind),1,vec(1,ind),1))
c$$$      ops = ops + n + n
c$$$ 75   hh(ind) = nrm1
c$$$      if (nrm1.le.zero) then
c$$$         ierr = -3
c$$$         return
c$$$      endif
c$$$c
c$$$c     scale the resulting vector
c$$$c
c$$$      fct = one / nrm1
c$$$      do 80 k = 1, n
c$$$         vec(k,ind) = vec(k,ind) * fct
c$$$ 80   continue
c$$$      ops = ops + n + 1
c$$$c
c$$$c     normal return
c$$$c
c$$$      ierr = 0
c$$$      return
c$$$c     end surbotine mgsro
c$$$      end
c$$$c-----------------------------------------------------------------------
      function zdistdot(n,x,ix,y,iy)
      integer n, ix, iy
      complex*16 zdistdot, x(*), y(*), zdotc
      external zdotc
      zdistdot = zdotc(n,x,ix,y,iy)
      return
      end
c-----end-of-distdot
