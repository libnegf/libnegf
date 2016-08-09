#ifndef _LNPARAMS_H
#define _LNPARAMS_H

#define MAXNCONT 10


struct lnparams {
  int verbose;
  int readoldsgf; 
  int spin;
  int kpoint;
  int iteration;
  double g_spin;
  double delta;
  double dos_delta;
  double eneconv;
  double wght;
  double ec;
  double ev;
  double emin;
  double emax;
  double estep;
  double mu_n[MAXNCONT];
  double mu_p[MAXNCONT];
  double mu[MAXNCONT];
  double contact_dos[MAXNCONT];
  int fictcont[MAXNCONT];
  double kbt[MAXNCONT];
  int np_n[2];
  int np_p[2];
  int np_real[11];
  int n_kt;
  int n_poles;
  int ni[MAXNCONT];
  int nf[MAXNCONT];
  char dore[1];
  };

#endif
