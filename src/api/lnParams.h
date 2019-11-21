#ifndef _LNPARAMS_H
#define _LNPARAMS_H

#define MAXNCONT 10


struct lnparams {
  int verbose;
  int readold_t_sgf;
  int readold_dm_sgf;
  int spin;
  int kpoint;
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
  _Bool fictcont[MAXNCONT];
  double kbt_dm[MAXNCONT];
  double kbt_t[MAXNCONT];
  int np_n[2];
  int np_p[2];
  int np_real[11];
  int n_kt;
  int n_poles;
  int ni[MAXNCONT];
  int nf[MAXNCONT];
  char dore[1];
  int min_or_max;
  _Bool is_s_is;
  };

#endif
