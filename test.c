/* 
 * File:   test.c
 * Author: guybillings
 *
 * Test program for functions in the PSO algorithm
 */

#include "bpnpsofuncs.h"
gsl_rng *rgen;
const int dimensions=5;

#define send_data_tag 2001
#define return_data_tag 2002

int my_id,root_process,num_procs;
double c;
int m;
int g;
int d;
int t;
int i;
double p;
double lfitness;
double xi;
double yi;
double delta;
int particles_to_compute=1;
MPI_Status status;

int main(int argc, char** argv) {

    int prec=200;
    int delta=1;
    int xi=1;
    int yi=1;
    int size;
    int maxsize=140;
    int minsize=40;
    int size_inc=20;
    double c=0.5;
    int d=4;
    int t=3;
    double p=0.1;
    double teff;

    for(size=minsize;size<=maxsize;size=size+size_inc)
    {
       teff=network_eff(c,size,size,d,t,p,xi,yi,delta,prec);
       printf("%s %i %s %f\n","Network size",size,":",teff);
    }
    
}

