#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#ifdef __INTEL_COMPILER
  #define MKL_Complex16 double complex
  #include "mkl.h"
#endif
#include <omp.h>

//#include "global_parameters.h"
//#include "simulation_parameters.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define KK(im, iK, iQ) ( KK[ im + Ndz * iK + Ndz * NK * iQ ] )
#define GG(iE, iQ)  ( GG[ iQ*NEloc + iE ] )
#define Sigma(iE, iQ)  ( Sigma[ iQ*NEloc + iE ] )
#define pSigma(mu, nu) ( pSigma[ nu*Np + mu ] )
#define pbuff1(mu, nu)  ( pbuff1[ nu*Np + mu ] )
#define pbuff2(mu, nu)  ( pbuff2[ nu*Np + mu ] )
#define sbuffH(mu, nu)  ( sbuffH[ nu*Np + mu ] )
#define rbuffH(mu, nu)  ( rbuffH[ nu*Np + mu ] )

// ======================================================================================
// Self-energy for polar optical electron-phonon interaction
//
// Sigma_mn(k,E) = Sum_q,G K(|zm - zn|, k, q+G) [f(-)*Gmn(q,E-wq) + f(+)*Gmn(q,E+wq)]
// ======================================================================================

int self_energy(
  // Fortran communicator of the 2d cartesian grid
  MPI_Fint *fcomm2d,
  // number of rows of the matrix block (Np*Mp)
  int Np,
  // number of cols of the matrix block (Np*Mp)
  int Mp,
  // global number of K-points
  int NK,
  // global number of E-points
  int NE,
  // local number of K-points
  int NKloc,
  // local number of E-points
  int NEloc,
  // hbar wq in energy grid units
  int iEhbaromega,
  // array of pointers to G (double complex *)
  void **GG,
  // array of pointers to Sigma (double complex *)
  void **Sigma,
  // a series of buffers size Np*Np (can they be reused ?)
  double complex *rbuff1,
  double complex *rbuff2,
  double complex *sbuffH,
  double complex *rbuffH,
  // fac_minus * G(E-wq)
  double complex fac_minus,
  // fac_plus * G(E+wq)
  double complex fac_plus,
  // Array of z coordinates on the grid indexed as the matrix rows
  int *izr,
  // Array of z coordinates on the grid indexed as the matrix cols
  int *izc,
  // Pretubulated array KK(|zi-zj|, iQ, iK)
  double *KK,
  // Leading dimension of array KK
  int Ndz
)
{
  int myid, size;
  int iK,iE, iQ;
  int iKglo,iEglo,iQglo, iQglo2;
  int mu, nu, im, im0, im1, im2, im3, im4, im5, im6, im7, in;
  int coords[2], coordsH[2], dims[2];
  int iEminus, iEplus;

  int msource, mdest;
  int psource, pdest;
  int hsource, hdest;
  int ndiff;
  int ndims;

  double complex *pbuff1, *pbuff2, *pGG, *pSigma;

  clock_t clock_start, clock_end;

  MPI_Comm comm2d;

  // dimensions of the cartesian grid
  ndims = 2;
  dims[0] = NK/NKloc;
  dims[1] = NE/NEloc;

  MPI_Request rqE[4];
  MPI_Status statusE[4];

  MPI_Request rqH[2];
  MPI_Status statusH[2];

  comm2d = MPI_Comm_f2c(*fcomm2d);

  MPI_Barrier(comm2d);
  MPI_Comm_size(comm2d, &size);
  MPI_Comm_rank(comm2d, &myid);
  // get process coordinates
  MPI_Cart_coords(comm2d,myid,ndims,coords);
  //DEBUG PRINTING
  //printf("Check Cart coords: %d %d \n", coords[0], coords[1]);
  //printf("check Np %d  Mp %d\n", Np, Mp);
  //printf("check NK %d  NKloc %d\n", NK, NKloc);
  //printf("check NE %d  NEloc %d\n", NE, NEloc);
  //printf("check iShift %d \n", iEhbaromega);
  //printf("check fac_minus %f %f \n", creal(fac1), cimag(fac1));
  //printf("check fac_plus %f %f \n", creal(fac2), cimag(fac2));
  //printf("check Kmat: \n");
  //for ( mu = 0; mu<Np; mu++){
  //  im = abs(izr[mu]-izc[0]);
  //  printf("KK(%d,0,0)=%f \n",im,KK(im,0,0));
  //}
  //
  //
  // Sigma_ij(iK, iE) = Sum_iQ   KK(|z_i-z_j|, iK, iQ) *
  //                        * (fac_minus * GG_ij(iQ, E-wq) + fac_plus * GG_ij(iQ, E+wq))
  //
  MPI_Barrier(comm2d);

  for( iQ=0; iQ<NKloc; iQ++ )
  {
    iQglo=iQ+coords[0]*NKloc;
    for( iE=0; iE<NEloc; iE++ )
    {
      iEglo=iE+coords[1]*NEloc;
      //printf("iQloc: %d iQglo: %d iEloc: %d iEglo: %d\n",iQ, iQglo, iE,iEglo);

      ////////////////////////////////////////////////////////////////////////////////////
      // Communications of G(k, E-hwq)
      //  012345 012345 012345 012345    iE < NEloc = 6
      //  -------------------------------------------------------------------------------
      //  000000 000011 111111 112222
      //  012345 678901 234567 890123    iEglo = 1 + 6 = 7
      // |oooooo|oooooo|oooooo|oooooo|
      //  ^-----|-E----|--^              ihbarOmega = 7; ndiff = 1 iMinus = 0
      //     ^--|----E-|-----^                           ndiff = 1 iMinus = 3
      ////////////////////////////////////////////////////////////////////////////////////

      // pbuff1 points to G(k,E-wq)
      // checks if iE-iEhbaromega is on the same processor => no communication
      // if iEglo<iEhbaromega => the processor is on the lower end of the energy grid
      //                      => communication is local or forced to be by truncation of G
      //                         such that G(E<0) = G(0), e.g., G(E-hw) = G(Elow)
      iEminus=iE-iEhbaromega;
      if (  iEminus >=0 || iEglo < iEhbaromega )
      {
        if (iEglo < iEhbaromega) {iEminus = 0;}
        pbuff1 = (double complex *) GG(iEminus, iQ);
        //printf("G(%d)=%g %g \n",iEminus+1,creal(pbuff1(0,0)),cimag(pbuff1(0,0)));
        mdest=MPI_PROC_NULL;
        msource=MPI_PROC_NULL;
      }

      // MPI Communication for G(E-wq)
      // get global iEminus and originating process that has to process
      // current working on iE
      if(dims[1] > 1)
      {
        if( iEhbaromega >= NEloc )
          ndiff = iEhbaromega / NEloc;
        else
          ndiff = 1;

        MPI_Cart_shift( comm2d, 1, ndiff, &msource, &mdest );
        // iE < iEhbaromega ensure MPI communication.
        if(mdest != MPI_PROC_NULL && iE < iEhbaromega)
        {
          // gets the local point to be sent to current iE such that pbuff1 => G(E-wq)
          // Original formula iEminus = (iE+NEloc-iEhbaromega) % NEloc;  Assumes iEhbaromega <= NEloc
          iEminus = (iE + (ndiff+1)*NEloc - iEhbaromega) % NEloc;
          //if (iEminus < 0 || iEminus >= NEloc){ printf("ERROR\n");}
          //printf("CPU %d iEminus=%d mdest=%d\n",coords[1], iEminus, mdest);
          pGG = (double complex *) GG(iEminus, iQ);
          //printf("CPU %d GG(%d,%d)=%x \n",coords[1], iEminus, iQ, pGG);
          MPI_Isend(pGG, Mp*Np, MPI_DOUBLE_COMPLEX, mdest, 41, comm2d, &rqE[0]);
        }

        if(msource != MPI_PROC_NULL && iE < iEhbaromega)
        {
          //printf("CPU %d iE=%d msource=%d\n",coords[1], iE, msource);
          MPI_Irecv(rbuff1, Mp*Np, MPI_DOUBLE_COMPLEX, msource, 41, comm2d, &rqE[1]);
          pbuff1 = rbuff1;
        }
      }

      ////////////////////////////////////////////////////////////////////////////////////
      // Communications of G(k, E+hwq)
      ////////////////////////////////////////////////////////////////////////////////////

      // checks if iE+iEhbaromega is on the same processor => no communication
      // pbuff2 points to G(k,E+wq)
      iEplus = iE+iEhbaromega;
      if( iEplus < NEloc || iEglo+iEhbaromega > NE )
      {
        if(iEglo+iEhbaromega>NE) {iEplus = NEloc-1;}
        pbuff2 = (double complex *) GG(iEplus, iQ);
        //printf("G(%d)=%g %g \n",iEplus+1,creal(pbuff2(0,0)),cimag(pbuff2(0,0)));
        pdest=MPI_PROC_NULL;
        psource=MPI_PROC_NULL;
      }

      // MPI Communication for G(E+wq)
      // get global iEminus and originating process that has to process
      // current working on iE
      if(dims[1] > 1)
      {
        if( iEhbaromega >= NEloc )
          ndiff = iEhbaromega / NEloc;
        else
          ndiff = (NEloc+iEhbaromega)/ NEloc;

        MPI_Cart_shift( comm2d, 1, -ndiff, &psource, &pdest );

        if(pdest != MPI_PROC_NULL && iE >= NEloc-iEhbaromega)
        {
          // gets the local point to be sent to current iE such that pbuff2 => G(E+wq)
          // Original formula iEplus = (iE-NEloc+iEhbaromega)   Assumes iEhbaromega <= NEloc
          iEplus = (iEglo + iEhbaromega) % NEloc;
          // printf("CPU# %d iE = %d iEglo = %d iEplus= %d send to: %d \n",coords[1], iE, iEglo, iEplus, pdest);
          //printf("CPU %d iEplus=%d  pdest=%d\n",coords[1], iEplus, pdest);
          pGG = (double complex *) GG(iEplus, iQ);
          //printf("CPU %d GG(%d,%d)=%x \n",coords[1], iEplus, iQ, pGG);
          MPI_Isend(pGG, Mp*Np, MPI_DOUBLE_COMPLEX, pdest, 42, comm2d, &rqE[2]);
        }
        if(psource != MPI_PROC_NULL && iE >= NEloc-iEhbaromega)
        {
          //printf("CPU# %d iE = %d iEglo = %d recv from: %d \n",coords[1], iE, iEglo, psource);
          //printf("CPU %d iE=%d  psource=%d\n",coords[1], iEplus, psource);
          MPI_Irecv(rbuff2, Mp*Np,MPI_DOUBLE_COMPLEX,psource,42,comm2d,&rqE[3]);
          pbuff2 = rbuff2;
        }
      }

      // Wait nodes for completed communications
      if(dims[1] > 1)
      {
        if(mdest != MPI_PROC_NULL  && iE<iEhbaromega) MPI_Wait(&rqE[0],&statusE[0]);
        if(msource != MPI_PROC_NULL  && iE<iEhbaromega) MPI_Wait(&rqE[1],&statusE[1]);

        if(pdest != MPI_PROC_NULL && iE >= NEloc-iEhbaromega) MPI_Wait(&rqE[2],&statusE[2]);
        if(psource != MPI_PROC_NULL && iE >= NEloc-iEhbaromega) MPI_Wait(&rqE[3],&statusE[3]);
      }

      ////////////////////////////////////////////////////////////////////////////////////
      // Communications of k-points
      ////////////////////////////////////////////////////////////////////////////////////
      // Update local iQ
      // KK(im, iK, iQ) ( KK[ im + Ndz * iK + Ndz * NK * iQ ] )
      // Sigma_ij(iQ, iE) = Sum_iK   KK(|z_i-z_j|, iK, iQ) *
      //                        * (fac_minus * GG_ij(iK, E-wq) + fac_plus * GG_ij(iK, E+wq))

      for( iK=0; iK<NKloc; iK++ )
      {
        iKglo=iK+coords[0]*NKloc;
        pSigma = (double complex *) Sigma(iE, iK);
        #pragma omp parallel for private(mu,nu,im)
        for( nu=0; nu<Mp; nu++ )
        {
          for( mu=0; mu<Np; mu++ )
          {
            im0 = abs(izr[mu]-izc[nu]);
            pSigma(mu,nu) += KK(im0, iKglo, iQglo) *
                             (fac_minus * pbuff1(mu,nu) + fac_plus * pbuff2(mu,nu);
          }
        }  
      }

      // MPI Communication over the k-grid
      if(dims[0] > 1)
      {
        #pragma omp parallel for private(mu,nu)
        for( nu=0; nu<Mp; nu++ )
          for( mu=0; mu<Np; mu++ )
            sbuffH(mu,nu) = fac_minus * pbuff1(mu,nu) + fac_plus * pbuff2(mu,nu);

        for( in=1; in<dims[0]; in++ )
        {
          MPI_Cart_shift(comm2d, 0, in, &hsource, &hdest);

          MPI_Cart_coords(comm2d, hsource, ndims, coordsH);

          iQglo2 = iQ + coordsH[0]*NKloc;

          // send to in+1%dims[0] = hdest
          MPI_Isend(sbuffH, Mp*Np, MPI_DOUBLE_COMPLEX, hdest, 43, comm2d, &rqH[0]);
          // recv from in+dims[0]%dims[0] = hsource
          MPI_Irecv(rbuffH, Mp*Np, MPI_DOUBLE_COMPLEX, hsource, 43, comm2d, &rqH[1]);

          MPI_Wait(&rqH[1], &statusH[1]);

          for( iK=0; iK<NKloc; iK++ )
          {
            iKglo = iK + coords[0]*NKloc;
            pSigma = (double complex *) Sigma(iE, iK);
            #pragma omp parallel for private(mu,nu,im)
            for( nu=0; nu<Mp; nu++ )
            {
              for( mu=0; mu<Np; mu++ )
              {
                im = abs(izr[mu]-izc[nu]);
                pSigma(mu,nu) += KK(im, iKglo, iQglo2) * rbuffH(mu,nu);
              }
            }
          }
          MPI_Wait(&rqH[0],&statusH[0]);

        }//k loop
      }//if dim[0]>1

    } //end of E-loop
  } //end of k-loop

  MPI_Barrier(comm2d);

  return 0;

}
