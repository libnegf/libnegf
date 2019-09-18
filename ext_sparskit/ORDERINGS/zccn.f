c----------------------------------------------------------------------c
c                          S P A R S K I T                             c
c----------------------------------------------------------------------c
c          REORDERING ROUTINES -- STRONGLY CONNECTED COMPONENTS        c 
c----------------------------------------------------------------------c
c     Contributed by:
C     Laura C. Dutto - email: dutto@cerca.umontreal.ca
c                      July 1992 - Update: March 1994
C-----------------------------------------------------------------------
c     CONTENTS:
c     --------
c     blccnx : Driver routine to reduce the structure of a  matrix 
c              to its strongly connected components.
c     cconex : Main routine to compute the strongly connected components
c              of a (block diagonal) matrix.
c     anccnx : We put in ICCNEX the vertices marked in the component MCCNEX.
c     newcnx : We put in ICCNEX the vertices marked in the component
c              MCCNEX. We modify also the vector KPW.
c     blccn1 : Parallel computation of the connected components of a
c              matrix. The parallel loop is performed only if the matrix
c              has a block diagonal structure.
c     ccnicopy:We copy an integer vector into anothoer.
c     compos : We calculate the composition between two permutation
c              vectors.
c     invlpw : We calculate the inverse of a permutation vector.
c     numini : We initialize a vector to the identity.
c     tbzero : We initialize to ZERO an integer vector.
c     iplusa : Given two integers IALPHA and IBETA, for an integer vector 
c              IA we calculate IA(i) = ialpha + ibeta * ia(i)
C
c----------------------------------------------------------------------c
      subroutine ZBLCCNX(n, nbloc, nblcmx, nsbloc, job, lpw, amat, ja,
     *                  ia, iout, ier, izs, nw)
C-----------------------------------------------------------------------
c     
c     This routine determines if the matrix given by the structure
c     IA et JA is irreductible. If not, it orders the unknowns such
c     that all the consecutive unknowns in KPW between NSBLOC(i-1)+1
c     and NSBLOC(i) belong to the ith component of the matrix.
c     The numerical values of the matrix are in AMAT. They are modified
c     only if JOB = 1 and if we have more than one connected component.
c
c     On entry:
c     --------
c     n      = row and column dimension of the matrix
c     nblcmx = maximum number of connected components allowed. The size
c              of NSBLOC is nblcmx + 1 (in fact, it starts at 0).
c     job    = integer indicating the work to be done:
c              job = 1  if the permutation LPW is modified, we
c                       permute not only the structure of the matrix
c                       but also its numerical values.
c              job.ne.1 if the permutation LPW is modified, we permute 
c                       the structure of the matrix ignoring real values.
c     iout   = impression parameter. If 0 < iout < 100, we print
c              comments and error messages on unit IOUT.
c     nw     = length of the work vector IZS.
c
c     Input / output:
c     --------------
c     nbloc  = number of connected components of the matrix. If the
c              matrix is not irreductible, nbloc > 1. We allow
c              nbloc > 1 on entry; in this case we calculate the
c              number of connected components in each previous one.
c     nsbloc = integer array of length NBLOC + 1 containing the pointers
c              to the first node of each component on the old (input)
c              and on the new (output) ordering.
c     lpw    = integer array of length N corresponding to the
c              permutation of the unknowns. We allow LPW to be a vector 
c              different from the identity on input.
c     amat   = complex(kind(1.0d0)) values of the matrix given by the structure IA, JA.
c     ja     = integer array of length NNZERO (= IA(N+1)-IA(1)) corresponding
c              to the column indices of nonzero elements of the matrix, stored
c              rowwise. It is modified only if the matrix has more
c              than one connected component.
c     ia     = integer array of length N+1 corresponding to the
c              pointer to the beginning of each row in JA (compressed
c              sparse row storage). It is modified only if
c              the matrix has more than one connected component.
c
c     On return:
c     ----------
c     ier    = integer. Error message. Normal return ier = 0.
c
c     Work space:
c     ----------
c     izs    = integer vector of length NW
c
C-----------------------------------------------------------------------
C     Laura C. Dutto - email: dutto@cerca.umontreal.ca
c                      July 1992 - Update: March 1994
C-----------------------------------------------------------------------
      integer izs(nw), lpw(n), nsbloc(0:nblcmx), ia(n+1), ja(*)
      complex(kind(1.0d0)) amat(*), zizs(nw)
      logical impr
      character(6) chsubr
C-----------------------------------------------------------------------
      ier    = 0
      impr   = iout.gt.0.and.iout.le.99
      ntb    = ia(n+1) - 1
      mxccex = max(nblcmx,20)
c.....The matrix AMAT is a complex(kind(1.0d0)) vector
      ireal  = 2
c
c.....MXPTBL: maximal number of vertices by block
      mxptbl = 0
      do ibloc = 1, nbloc
         mxptbl = max( mxptbl, nsbloc(ibloc) - nsbloc(ibloc-1))
      enddo
c
      long1 = nbloc * mxptbl
      long2 = nbloc * (mxccex+1)
c.....Dynamic allocation of memory
      iend   = 1
      iiend  = iend
      ilpw   = iiend 
      ikpw   = ilpw   + n
      ilccnx = ikpw   + long1
      imark  = ilccnx + long2
      iend   = imark  + n
      if(iend .gt. nw) go to 220
c
      nbloc0 = nbloc
      chsubr = 'BLCCN1'
c.....We determine if the matrix has more than NBLOC0 connected components.
      call BLCCN1(n, nbloc, nblcmx, nsbloc, izs(ilpw), izs(ikpw), ia,
     *            ja, izs(imark), mxccex, izs(ilccnx), mxptbl, iout,
     *            ier)
      if(ier.ne.0) go to 210
c
      if(nbloc .gt. nbloc0) then
c..........The matrix has more than NBLOC0 conneted components. So, we
c..........modify the vectors IA and JA to take account of the new permutation.
           nfree = iend - ikpw
           call tbzero(izs(ikpw), nfree)
           iiat  = ikpw
           ijat  = iiat + n + 1
           iamat = ijat + ntb
           iend  = iamat
           if(job .eq. 1) iend = iamat + ireal * ntb
           if(iend .gt. nw) go to 220
c
c..........We copy IA and JA on IAT and JAT respectively
           call ccnicopy(n+1, ia, izs(iiat))
           call ccnicopy(ntb, ja, izs(ijat))
           if(job .eq. 1) call zcopy(ntb, amat, 1, izs(iamat), 1)
           zizs(1:nw)=cmplx(izs(1:nw))
           call zdperm(n, zizs(iamat), izs(ijat), izs(iiat), amat,
     *                ja, ia, izs(ilpw), izs(ilpw), job)
           ipos = 1
c..........We sort columns inside JA.
           zizs(1:nw)=cmplx(izs(1:nw))
           call zcsrcsc(n, job, ipos, amat, ja, ia, zizs(iamat),
     *                 izs(ijat), izs(iiat))
           zizs(1:nw)=cmplx(izs(1:nw))
           call zcsrcsc(n, job, ipos, zizs(iamat), izs(ijat), izs(iiat),
     *                 amat, ja, ia)
      endif
c.....We modify the ordering of unknowns in LPW
      call compos(n, lpw, izs(ilpw))
c
 120  nfree = iend - iiend
      call tbzero(izs(iiend), nfree)
      iend = iiend
      return
c
 210  IF(IMPR) WRITE(IOUT,310) chsubr,ier
      go to 120
 220  IF(IMPR) WRITE(IOUT,320) nw, iend
      if(ier.eq.0) ier = -1
      go to 120
c
 310  FORMAT(' ***BLCCNX*** ERROR IN ',a6,'. IER = ',i8)
 320  FORMAT(' ***BLCCNX*** THERE IS NOT ENOUGH MEMORY IN THE WORK',
     1       ' VECTOR.'/13X,' ALLOWED MEMORY = ',I10,'  - NEEDED',
     2       ' MEMORY = ',I10)
      end 
