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
#include <time.h>
#include <omp.h>

//#include "global_parameters.h"
//#include "simulation_parameters.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define KK(im, iK, iQ) ( KK[ im + Ndz * iK + Ndz * NK * iQ ] )


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
  double complex *sbuff1,
  double complex *sbuff2,
  double complex *rbuff1,
  double complex *rbuff2,
  double complex *rbuffH,
  double complex *sbuffH,
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
  int mu, nu, im, in;
  int coords[2], coordsH[2], dims[2];
  int iEminus, iEplus;

  int msource, mdest;
  int psource, pdest;
  int hsource, hdest;
  int ndiff;
  int ndims;

  double complex *pbuff1, *pbuff2, *pGG, *pSigma;

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
      
      // checks if iE-iEhbaromega is on the same processor => no communication
      // pbuff1 points to G(k,E-wq)
      if( (((int)iE - iEhbaromega) >=0 && (iE - iEhbaromega) < NEloc) || iEglo < NEloc)
      {
        iEminus=iE-iEhbaromega;
        if(iEminus<=0) {iEminus = 0;}
        pbuff1 = (double complex *) GG[iQ*NEloc+iEminus];
        mdest=MPI_PROC_NULL;
        msource=MPI_PROC_NULL;
      }

      // MPI Communication for G(E-wq) 
      // get global iEminus and originating process that has to process
      // current working on iE 
      if(dims[1] > 1)
      {
        iEminus = iEglo-iEhbaromega;
        if( iEhbaromega >= NEloc )
          ndiff = iEhbaromega / NEloc;
        else
          ndiff = (NEloc+iEhbaromega)/ NEloc;

        MPI_Cart_shift( comm2d, 1, ndiff, &msource, &mdest );

        // iE < iEhbaromega ensuhre MPI communication.
        if(mdest != MPI_PROC_NULL && iE < iEhbaromega)
        {
          // iEminus sends to current iE so pbuff1 => G(E-wq) 
          iEminus = (NEloc+iE-iEhbaromega) % NEloc;
          pGG = (double complex *) GG[iQ*NEloc+iEminus];
          MPI_Isend(pGG, Mp*Np, MPI_DOUBLE_COMPLEX, mdest, 41, comm2d, &rqE[0]);
        }

        if(msource != MPI_PROC_NULL && iE < iEhbaromega)
        {
          MPI_Irecv(rbuff1, Mp*Np, MPI_DOUBLE_COMPLEX, msource, 41, comm2d, &rqE[1]);
          pbuff1 = rbuff1;
        }
      }

      ////////////////////////////////////////////////////////////////////////////////////
      // Communications of G(k, E+hwq)
      ////////////////////////////////////////////////////////////////////////////////////

      // checks if iE+iEhbaromega is on the same processor => no communication
      // pbuff2 points to G(k,E+wq)
      if( (((int)iE + iEhbaromega) >=0 && (iE + iEhbaromega) < NEloc) || iEglo >= NE-NEloc)
      {
        iEplus = iE+iEhbaromega;
        if(iEglo+iEhbaromega>=NE) {iEplus = NEloc-1;}
        pbuff2 = (double complex *) GG[iQ*NEloc+iEplus];
        pdest=MPI_PROC_NULL;
        psource=MPI_PROC_NULL;
      }

      // MPI Communication for G(E+wq) 
      // get global iEminus and originating process that has to process
      // current working on iE 
      if(dims[1] > 1)
      {
        iEplus = iEglo + iEhbaromega;
        if( iEhbaromega >= NEloc )
          ndiff = iEhbaromega / NEloc;
        else
          ndiff = (NEloc+iEhbaromega)/ NEloc;

        MPI_Cart_shift( comm2d, 1, -ndiff, &psource, &pdest );
        // iEplus sends to current processing iE so pbuff2 => G(E+wq) 
        if(pdest != MPI_PROC_NULL && iE >= NEloc-iEhbaromega)
        {
          iEplus = iE - (NEloc - iEhbaromega);
          pGG = (double complex *) GG[iQ*NEloc+iEplus];
          MPI_Isend(pGG, Mp*Np, MPI_DOUBLE_COMPLEX, pdest, 42, comm2d, &rqE[2]);
        }
        if(psource != MPI_PROC_NULL && iE >= NEloc-iEhbaromega)
        {
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
     #pragma omp parallel for private(iK,iKglo,mu,nu,im) collapse(2)
      for( iK=0; iK<NKloc; iK++ )
      {
        iKglo=iK+coords[0]*NKloc;
        pSigma = (double complex *) Sigma[iK*NEloc + iE];
        for( nu=0; nu<Mp; nu++ )
        {
          for( mu=0; mu<Np; mu++ )
          {
            im = abs(izr[mu]-izc[nu]);
            pSigma[nu*Np + mu] += KK(im, iKglo, iQglo) *
                                 (fac_minus * pbuff1[nu*Np+mu] + fac_plus * pbuff2[nu*Np+mu]);

            //printf("Sigma(mu,nu): %f \n",pSigma[nu*Np+mu]);
            //im = abs(iz-jz);
            // F(iQ, iK, |iz-jz|) = F[Np*NK*iQglo + Np*iKglo + im]
            //sigless[(iQ*NEloc+iE)*Np*Np+jz*Np+iz] += dK * qe * hbaromegaLO/(4.0*pow(M_PI,2.0)) *
            //  (1.0/(epsinfr*eps0) - 1.0/(eps0r*eps0)) * K[iKglo] * F[ im+ Np*iKglo + Np*NK*iQglo ] *
            //  Mtilde * (fac_minus * pbuff1[jz*Np+iz] + fac_plus * pbuff2[jz*Np+iz]);
          }
        }
      }

      // MPI Communication over the k-grid
      //printf("global k communications if %d > 1 \n",dims[0]);
      if(dims[0] > 1)
      {
        //printf("Communication for Q integration\n");
        #pragma omp parallel for private(mu,nu)
        for( nu=0; nu<Mp; nu++ )
          for( mu=0; mu<Np; mu++ )
            sbuffH[nu*Np+mu] = fac_minus * pbuff1[nu*Np+mu] + fac_plus * pbuff2[nu*Np+mu];


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

          #pragma omp parallel for private(iK,iKglo,mu,nu,im) collapse(2)
          for( iK=0; iK<NKloc; iK++ )
          {
            iKglo = iK + coords[0]*NKloc;
            pSigma = (double complex *) Sigma[iK*NEloc + iE];
            for( nu=0; nu<Mp; nu++ )
            {
              for( mu=0; mu<Np; mu++ )
              {
                im = abs(izr[mu]-izc[nu]);
                pSigma[nu*Np + mu] += KK(im, iKglo, iQglo2) * rbuffH[nu*Np+mu];
                //im = abs(iz-jz);
                //sigless[(iQ*NEloc+iE)*Np*Np+jz*Np+iz] += dK * qe * hbaromegaLO/(4.0*pow(M_PI,2.0)) *
                //  (1.0/(epsinfr*eps0) - 1.0/(eps0r*eps0)) * K[iKglo2] * F[ im+ Np*iKglo + Np*NK*iQglo2 ] *
                //  Mtilde * rbuffH[jz*Np+iz];
              }
            }
          }

          MPI_Wait(&rqH[0],&statusH[0]);
        }
      }

    } //end of E-loop
  } //end of k-loop

  MPI_Barrier(comm2d);

  return 0;

}
