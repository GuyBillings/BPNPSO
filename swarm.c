/* Functions required by the particle swarm
   Guy Billings UCL 2011                                */

#include "bpnpsofuncs.h"

void update_posn(gsl_matrix *posns,gsl_matrix *tprobs,int limits[],int dmin,int particles,int init)
{
    int particle;
    int dim;
    int lim;
    int x;
    double coin=0;
    double probe=0;
    int newpos;

    //Format of positions is {mu, gamma, d, phi, pin}
    if(init==1)//Then intialise positions at random
    {
        for(particle=0;particle<particles;particle++)
        {
            for(dim=0;dim<dimensions;dim++)
            {
              if(dim!=3)
                x=gsl_rng_uniform_int(rgen,limits[dim])+1;
              else if(dim==3)
                x=gsl_rng_uniform_int(rgen,(int)gsl_matrix_get(posns,particle,dim-1)+dmin-1)+1;
              
              gsl_matrix_set(posns,particle,dim,(double)x);
            }

        }
    }
    else//Update positions based on the transition probabilities
    {
        for(particle=0;particle<particles;particle++)
        {
            for(dim=0;dim<dimensions;dim++)
            {
              //Flip a coin to decide whether to
              //test for an up move or a down move
              coin=gsl_rng_uniform(rgen);
              if(coin<0.5)//Then probe for moving up
              {
                  probe=gsl_rng_uniform(rgen);
                  if(probe<=gsl_matrix_get(tprobs,particle,2*dim))
                  {

                      if(dim==3)
                          lim=(int)gsl_matrix_get(posns,particle,2)+dmin-1;
                      else
                          lim=limits[dim];

                      newpos=(int)gsl_matrix_get(posns,particle,dim)+1;
                      if(newpos>lim)
                          newpos=lim;
                      gsl_matrix_set(posns,particle,dim,(double)newpos);
                  }
              }
              else if(coin>0.5)//Then probe for moving down
              {
                  probe=gsl_rng_uniform(rgen);
                  if(probe<=gsl_matrix_get(tprobs,particle,2*dim+1))
                  {
                      newpos=(int)gsl_matrix_get(posns,particle,dim)-1;
                      if(newpos<1)
                          newpos=1;
                      gsl_matrix_set(posns,particle,dim,(double)newpos);
                  }
              }

              if(gsl_matrix_get(posns,particle,3)>gsl_matrix_get(posns,particle,2))
                  gsl_matrix_set(posns,particle,3,gsl_matrix_get(posns,particle,2));
              
            }
        }
    }


}

void update_fitness(double fitness[],gsl_matrix *sw_current_posn,int particles,int mumin,int gammamin,int dmin,
        int pmin,double pinc,double capacity[],double xi,double yi,double delta,int prec,int send_data_tag,
        int return_data_tag,int debug)
{

     int block=num_procs;
     int particle=0;
     int rank;
 
     while(particles_to_compute==1)
     {

       for(rank=1;rank<block;rank++)
       {
          c=capacity[particle];
          MPI_Send(&particles_to_compute,1,MPI_INT,rank,send_data_tag,MPI_COMM_WORLD);
          MPI_Send(&c,1,MPI_DOUBLE,rank,send_data_tag,MPI_COMM_WORLD);
          m=(int)gsl_matrix_get(sw_current_posn,particle,0)+mumin-1;
          MPI_Send(&m,1,MPI_INT,rank,send_data_tag,MPI_COMM_WORLD);
          g=(int)gsl_matrix_get(sw_current_posn,particle,1)+gammamin-1;
          MPI_Send(&g,1,MPI_INT,rank,send_data_tag,MPI_COMM_WORLD);
          d=(int)gsl_matrix_get(sw_current_posn,particle,2)+dmin-1;
          MPI_Send(&d,1,MPI_INT,rank,send_data_tag,MPI_COMM_WORLD);
          t=(int)gsl_matrix_get(sw_current_posn,particle,3);
          MPI_Send(&t,1,MPI_INT,rank,send_data_tag,MPI_COMM_WORLD);
          p=pinc*((int)gsl_matrix_get(sw_current_posn,particle,4)-1)+pmin;
          MPI_Send(&p,1,MPI_DOUBLE,rank,send_data_tag,MPI_COMM_WORLD);
          MPI_Send(&xi,1,MPI_DOUBLE,rank,send_data_tag,MPI_COMM_WORLD);
          MPI_Send(&yi,1,MPI_DOUBLE,rank,send_data_tag,MPI_COMM_WORLD);
          MPI_Send(&delta,1,MPI_DOUBLE,rank,send_data_tag,MPI_COMM_WORLD);
          MPI_Send(&prec,1,MPI_INT,rank,send_data_tag,MPI_COMM_WORLD);
          if(debug==1)
          {
            printf("%s %i %s %i %s %i\n","Rank: ",my_id," Pariticle: ",particle," Sent to rank: ",rank);
            fflush(stdout);
          }
          particle=particle+1;
       }

       if(particle<particles)
       {
           if(debug==1)
           {
             printf("%s %i\n","Root process computing particle: ",particle);
             fflush(stdout);
           }

         fitness[particle]=network_eff(capacity[particle],(int)gsl_matrix_get(sw_current_posn,particle,0)+mumin-1,(int)gsl_matrix_get(sw_current_posn,particle,1)+gammamin-1,
                (int)gsl_matrix_get(sw_current_posn,particle,2)+dmin-1,(int)gsl_matrix_get(sw_current_posn,particle,3),pinc*((int)gsl_matrix_get(sw_current_posn,particle,4)-1)+pmin,
                xi,yi,delta,prec);
       particle=particle+1;
       }

       for(rank=1;rank<block;rank++)
       {
          MPI_Recv(&lfitness,1,MPI_DOUBLE,rank,return_data_tag,MPI_COMM_WORLD,&status);
          //sender=status.MPI_SOURCE;
          fitness[rank]=lfitness;
          if(debug==1)
          {
             printf("%s %i\n","Data received from rank: ",rank);
             fflush(stdout);
          }
       }

       if((particles-particle)<block)
           block=particles-particle;

       if(particle>=particles)
           particles_to_compute=0;

     }

}

double update_global_pos(double maximum_fitness,int particles,int dimensions,gsl_matrix *posns,int global_best_position[],double global_best_fitness,double fitness[])
{
    int particle, dim;
    int fittest_particle;
    int mod=0;
    for(particle=0;particle<particles;particle++)
    {
     if(fitness[particle]>maximum_fitness)
     {
         maximum_fitness=fitness[particle];
         fittest_particle=particle;
         mod=1;
     }
    }
    if(mod!=0)
    {
      for(dim=0;dim<dimensions;dim++)
          global_best_position[dim]=(int)gsl_matrix_get(posns,fittest_particle,dim);
    }

    return maximum_fitness;
}

void update_amps(int particles,int dimesnions,gsl_matrix *sw_tamps,gsl_matrix *sw_current_posn,gsl_matrix *sw_bestpos,
        int limits[],int global_best_pos[],double omega,double rp,double rg, double stat_amp)
{
    int particle;
    int degree_f;
    int degrees_of_f=2*dimensions;
    double amp;
    int l;

    for(particle=0;particle<particles;particle++)
    {
        for(degree_f=1;degree_f<=degrees_of_f;degree_f++)
        {

            if(degree_f%2!=0)
            {

            if(((degree_f+1)/2)==4)
                l=limits[((degree_f+1)/2)-2];
            else
                l=limits[((degree_f+1)/2)-1];

            amp=omega*gsl_matrix_get(sw_tamps,particle,degree_f-1)+rp*(gsl_matrix_get(sw_bestpos,particle,((degree_f+1)/2)-1)-
                    gsl_matrix_get(sw_current_posn,particle,((degree_f+1)/2)-1))/l+rg*(global_best_pos[((degree_f+1)/2)-1]-
                    gsl_matrix_get(sw_current_posn,particle,((degree_f+1)/2)-1))/l;

            if(amp<0)
                amp=0;
            else if(amp>1)
                amp=1;

            gsl_matrix_set(sw_tamps,particle,degree_f-1,amp);
            }
            else if(degree_f%2==0)
            {
            amp=1-gsl_matrix_get(sw_tamps,particle,degree_f-2)-stat_amp;
            if(amp<0)
                amp=0;
            else if(amp>1)
                amp=1;
            gsl_matrix_set(sw_tamps,particle,degree_f-1,amp);
            }

        }
    }
}
void amps_to_probs(int particles,gsl_matrix *tprobs, gsl_matrix *amps,double stationary_prob)
{
    int particle;
    int degrees_of_f=2*dimensions;
    int degree_f;
    double total=0;

    for(particle=0;particle<particles;particle++)
    {
        for(degree_f=1;degree_f<=degrees_of_f;degree_f++)
        {
            total=total+gsl_matrix_get(amps,particle,degree_f-1);
            if(degree_f%2==0)
            {
                //printf("%f\n",gsl_matrix_get(amps,particle,degree_f-2)/(total+stationary_prob));
                gsl_matrix_set(tprobs,particle,degree_f-2,gsl_matrix_get(amps,particle,degree_f-2)/(total+stationary_prob));
                gsl_matrix_set(tprobs,particle,degree_f-1,gsl_matrix_get(amps,particle,degree_f-1)/(total+stationary_prob));
                total=0;
            }
        }
    }
}

void update_capacity(gsl_matrix *posn,int tsteps,int particles,double capacity[])
{
    int particle;
    for(particle=0;particle<particles;particle++)
    {
        capacity[particle]=(log(tsteps)/log(2))/gsl_matrix_get(posn,particle,0);
        if(capacity[particle]>1)
            capacity[particle]=1;
        //THIS could be a fudge. Should probably calculate in terms of the
        //source entropy (rather than raw timestep entropy)
    }
}
