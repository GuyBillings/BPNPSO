#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//#define __GSL_RANGE_CHECK_OFF__

void update_posn(gsl_matrix *posns,gsl_matrix *tprobs,int limits[],int dmin,int particles,int init);
void update_fitness(double fitness[],gsl_matrix *sw_current_posn,int particles,int mumin,int gammamin,int dmin,
        int pmin,double pinc,double capacity[],double xi,double yi,double delta,int prec,int send_data_tag,
        int return_data_tag,int debug);
double update_global_pos(double maximum_fitness,int particles,int dimensions,gsl_matrix *posns,int global_best_position[],double global_best_fitness,double fitness[]);
void update_amps(int particles,int dimesnions,gsl_matrix *sw_tamps,gsl_matrix *sw_current_posn,gsl_matrix *sw_bestpos,
        int limits[],int global_best_pos[],double omega,double rp,double rg,double stat_amp);
void amps_to_probs(int particles,gsl_matrix *tprobs, gsl_matrix *amps,double stationary_prob);
void update_capacity(gsl_matrix *posn,int tsteps,int particles,double capacity[]);