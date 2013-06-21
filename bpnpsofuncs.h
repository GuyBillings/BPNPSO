#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include <gsl/gsl_matrix.h>
#include "swarm.h"
#include <time.h>
#include <string.h>
#include <mpi.h>

//#define __GSL_RANGE_CHECK_OFF__

extern gsl_rng *rgen;
extern const int dimensions;
extern int my_id,root_process,num_procs,particles_to_compute;
extern double c, p, lfitness, xi, yi, delta;
extern int m, g, d, t, i;
extern MPI_Status status;

void binomial(mpz_t out,int n,int k);
void inpdf(mpfr_t *out,mpz_t n,int psize,double pini,double pinc,unsigned long mu,mpfr_t bcs[],mpfr_prec_t prec);
void get_hypergeometricpdf(mpfr_t *rval,unsigned int mi,int sclass,int conns, int phi,mpz_t key,HFNARGS table_prms,NEXT_ILIST_UNIT* keytable,NEXT_DLIST_UNIT* datatable,int prec);
void pphi(mpfr_t *pout,mpz_t mu,mpz_t dmp,int m,int d, int phi,mpz_t bc_b_ai[], mpz_t mu_bc[], mpfr_prec_t prec);
void uniquecodes(mpfr_t *ucodes, int psiz, mpz_t n, unsigned long int mu, mpz_t ncodes ,mpfr_t *pdf,mpfr_prec_t prec);
void scjoint(mpfr_t *distptr, int mu, int gamma, int psiz, int conns, int phi, mpfr_t *bcs,mpfr_t *binpdf,mpfr_t *pp,mpfr_prec_t prec);
void tabulate_bins_fr(mpfr_t *bcs, int Nini, int Nend, int mini, int mend, mpfr_prec_t prec);
void tabulate_bins_z(mpz_t bcs[], int Nini, int Nend, int mini, int mend);
void eff_point(mpfr_t *ent, mpfr_t *ipdf, mpfr_t *jdist, mpfr_t *pp, mpfr_t *gamma_bcs, mpfr_t *p0, int conns, int phi, int mu, int gamma, double pin, double xi, double yi, double delta, mpfr_prec_t prec);
unsigned long long fact(unsigned long long x);
void binomialpdf(double *bcs, int n,double p);
double network_eff(double bfrac, int mu, int gamma, int d, int phi, double pin,double xi, double yi, double delta, int prec);
