!----------------------------------------------------------------------
! libNEGF  
!
! Source codes for sparse matrix operations and algebra 
!----------------------------------------------------------------------
! Source codes originally coded in Sparsekit
!           Y. Saad, Sep. 21 1989                                      
!----------------------------------------------------------------------
! 
! Re-written in modern Fortran in 2022
! A. Pecchia 
! D. Soccodato
!
!----------------------------------------------------------------------

module skit_blassm
  use ln_precision
  implicit none
  private

  contains


end module skit_blassm
! BLAS OPERATIONS
! ----------------------
!    subroutine amub(nrow, ncol, job, aa, ja, ia, bb, jb, ib, cc, jc, ic, nnz, iw, ierr)
!    
!    subroutine zamub(nrow, ncol, job, aa, ja, ia, bb, jb, ib, cc, jc, ic, nnz, iw, ierr)
!          
!    subroutine zamubs(nrow, ncol, job, aa, ja, ia, s, bb, jb, ib, cc, jc, ic, nnz, iw, ierr)
!    
!    subroutine aplb(nrow, ncol, job, aa, ja, ia, bb, jb, ib, cc, jc, ic, nnz, iw, ierr)
!    
!    subroutine aplsb(nrow, ncol, job, aa, ja, ia, s, bb, jb, ib, cc, jc, ic, nnz, iw, ierr)
!
!    subroutine zaplb(nrow, ncol, job, aa, ja, ia, bb, jb, ib, cc, jc, ic, nnz, iw, ierr)
!    
!    subroutine zaplsb(nrow, ncol, job, aa, ja, ia, s, bb, jb, ib, cc, jc, ic, nnz, iw, ierr)
!
!    subroutine zaplb1(nrow, ncol, job, aa, ja, ia, bb, jb, ib, cc, jc, ic, nnz, ierr)
!    
!    subroutine zcplsamub(nrow, ncol, job, aa, ja, ia, s, bb, jb, ib, cc, jc, ic, nnz, iw, ierr)
!
!    subroutine zcopmat(nrow,aa,ja,ia,bb,jb,ib,job,ierr)
!  
!    subroutine zas1pls2b(nrow, ncol, job, aa, ja, ia, s1, s2, bb, jb, ib, &
!                 & cc, jc, ic, nnz, iw, ierr)
