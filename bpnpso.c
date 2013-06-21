/* BPNPSO: Particle Swarm Optimisation of
   Binary Partition Networks
   Guy Billings UCL 2011*/
   
#include "bpnpsofuncs.h"
gsl_rng *rgen;
const int dimensions=5;

#define send_data_tag 2001
#define return_data_tag 2002

NEXT_DLIST_UNIT* mpfr_bcs_dattab;
NEXT_ILIST_UNIT* mpfr_bcs_keytab;
HFNARGS mpfr_bcs_prms;

NEXT_DLIST_UNIT* mpfr_ipdf_dattab;
NEXT_ILIST_UNIT* mpfr_ipdf_keytab;
HFNARGS mpfr_ipdf_prms;

NEXT_ILIST_UNIT* mpfr_p0_keytab;
NEXT_DLIST_UNIT* mpfr_p0_dattab;
HFNARGS mpfr_p0_prms;

NEXT_ILIST_UNIT* mpz_bcs_dattab;
NEXT_ILIST_UNIT* mpz_bcs_keytab;
HFNARGS mpz_bcs_prms;

NEXT_ILIST_UNIT* hyp_keytab;
NEXT_DLIST_UNIT* hyp_dattab;
HFNARGS hyp_table_prms;

NEXT_ILIST_UNIT* mpfr_joint_keytab;
NEXT_DLIST_UNIT* mpfr_joint_dattab;
HFNARGS mpfr_joint_prms;

NEXT_ILIST_UNIT* mpfr_pp_keytab;
NEXT_DLIST_UNIT* mpfr_pp_dattab;
HFNARGS mpfr_pp_prms;

NEXT_ILIST_UNIT* mpfr_ent_keytab;
NEXT_DLIST_UNIT* mpfr_ent_dattab;
HFNARGS mpfr_ent_prms;

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

int main(int argc, char* argv[])
{
   
   //SET UP MPI ENVIRONMENT:
   MPI_Init(&argc, &argv);
   root_process=0;
   MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   int debug=0;
   int rank;
   int prec;

   //CODE FOR ROOT PROCESS TO EXECUTE:
   if(my_id == root_process)
   {
     int particles;
     int iterations;
     int tsteps;
     int mumin;
     int mumax;
     int gammamin;
     int gammamax;
     int dmin;
     int dmax;
     double pmin;
     double pinc;
     double pmax;
     double omega;
     double rp;
     double rg;
     double stationary_prob;

     printf("\n");
     printf("%s\n","BPNPSO: The Binary Partition Network Particle Swarm Optimiser");
     printf("%s\n","Guy Billings UCL 2011");
     printf("%s\n","---------------------");

     double args[20];
     int conf_file=1;
     FILE *parameters;
     char instring[150];

     char *fullppath;
     char *fullbpath;
     char *fullfpath;

     char pfname[]={'p','o','s','i','t','i','o','n','s','_','\0'};
     char bfname[]={'b','e','s','t','p','o','s','_','\0'};
     char ffname[]={'f','i','t','n','e','s','s','_','\0'};
     char ext[]={'.','t','x','t','\0'};

     time_t starttime;
     time_t endtime;
     struct tm *timeinfo;
     time(&starttime);
     timeinfo=localtime(&starttime);
     char *stamp;
     stamp=asctime(timeinfo);
     time(&starttime);
     char* strendtime;

     char fstamp[]={stamp[4],stamp[5],stamp[6],'_',stamp[8],stamp[9],'_',stamp[20],stamp[21],stamp[22],stamp[23],'_',
                    stamp[11],stamp[12],'-',stamp[14],stamp[15],'-',stamp[17],stamp[18],'\0'};

     if(argc==1)
     {
       conf_file=0;
       prec=200;
       tsteps=1000000000;
       particles=20;
       mumin=20;
       mumax=100;
       gammamin=20;
       gammamax=100;
       dmin=1;
       dmax=20;
       pmin=0.1;
       pinc=0.1;
       pmax=0.9;
       xi=50;
       yi=50;
       delta=1;
       omega=1;
       rp=0;
       rg=0;
       stationary_prob=0;
       iterations=2;
       fullppath=(char*)malloc(39);
       fullbpath=(char*)malloc(35);
       fullfpath=(char*)malloc(37);
       char path[]={'.','/','\0'};
       strcpy(fullppath,path);
       strcpy(fullbpath,path);
       strcpy(fullfpath,path);
       strcat(fullppath,pfname);
       strcat(fullbpath,bfname);
       strcat(fullfpath,ffname);
       // TODO : space is being inserted in the file name after the month
       strcat(fullppath,fstamp);
       strcat(fullppath,ext);
       strcat(fullbpath,fstamp);
       strcat(fullbpath,ext);
       strcat(fullfpath,fstamp);
       strcat(fullfpath,ext);
     }
     else
     {
       parameters=fopen(argv[1],"r");
       if(parameters==NULL)
       {
          printf("%s\n","Error reading configuration file.");
          printf("%s\n","Exiting.");
          exit(EXIT_FAILURE);
       }
       for(i=0;i<21;i++)
       {
          fgets(instring,50,parameters);
          if(i<20)
          {
            args[i]=atof(instring);
          }
          //if(instring==NULL)
          //{
          //    instring[0]='.';
          //    instring[1]='/';
          //    instring[2]='\0';
          //}

       }

       char *path;
       int plen=strlen(instring);
       path=(char*)malloc(plen);
       for(i=0;i<plen;i++)
         path[i]=instring[i];
       //if(path[i-1]=='\n')
       path[i-1]='\0';

       fgets(instring,150,parameters);
       fclose(parameters);

       char *fileid;
       int fdlen=strlen(instring);
       fileid=(char*)malloc(fdlen);
       for(i=0;i<fdlen;i++)
         fileid[i]=instring[i];
       //if(fileid[i-1]=='\n')
       fileid[i-1]='\0';

       //File contents (in separate lines)
       //prec,bfrac,particles,mumin,mumax,gammamin,gammamax,dmin,dmax,pmin,pinc,pmax,xi,yi,delta,
       //omega,rp,rg,stationary_prob,iterations,output path
       prec=(int)args[0];
       tsteps=(int)args[1];
       particles=(int)args[2];
       mumin=(int)args[3];
       mumax=(int)args[4];
       gammamin=(int)args[5];
       gammamax=(int)args[6];
       dmin=(int)args[7];
       dmax=(int)args[8];
       pmin=args[9];
       pinc=args[10];
       pmax=args[11];
       xi=args[12];
       yi=args[13];
       delta=args[14];
       omega=args[15];
       rp=args[16];
       rg=args[17];
       stationary_prob=args[18];
       iterations=(int)args[19];

       //if(path[i]!='\n')
       //{
       //    printf("%s\n","Warning: No newline at end of parameters file.");
       //}
       //if(path[i-1]!='/')
       //{
       //    printf("%s\n","Warning: No trailing / in path.");
       //}
                    
       fullppath=(char*)malloc(plen+16+fdlen);
       fullbpath=(char*)malloc(plen+14+fdlen);
       fullfpath=(char*)malloc(plen+14+fdlen);
       strcpy(fullppath,path);
       strcat(fullppath,pfname);
       strcpy(fullbpath,path);
       strcat(fullbpath,bfname);
       strcpy(fullfpath,path);
       strcat(fullfpath,ffname);
       strcat(fullppath,fileid);
       strcat(fullppath,ext);
       strcat(fullbpath,fileid);
       strcat(fullbpath,ext);
       strcat(fullfpath,fileid);
       strcat(fullfpath,ext);
       
     }
     
     //passes=(int)floor(particles/num_procs);
     //remainder=particles-(passes*num_procs);

     //for(i=1;i<num_procs;i++)
     //{
     //  MPI_Send(&passes,1,MPI_INT,i,send_data_tag,MPI_COMM_WORLD);
     //  MPI_Send(&remainder,1,MPI_INT,i,send_data_tag,MPI_COMM_WORLD);
     //  MPI_Send(&iterations,1,MPI_INT,i,send_data_tag,MPI_COMM_WORLD);
     //}

     if(conf_file==0)
      printf("%s\n","No path to configuration supplied, running test.");

     rgen=gsl_rng_alloc(gsl_rng_taus);
     int seed= 1;
     gsl_rng_set(rgen,seed);
     
     double fitness[particles];
     double best_fitness[particles];
     double global_best_fitness=0;
     int global_best_position[dimensions];
     int iteration;
     int j;
     int dim;
     int particle;
     int block;
     double capacity[particles];
 
     FILE* positions;
     FILE* bestpositions;
     FILE* fitfile;
     positions=fopen(fullppath,"a+");
     bestpositions=fopen(fullbpath,"a+");
     fitfile=fopen(fullfpath,"a+");
     
     if(num_procs>particles)
     {
       printf("%s\n","The number of processors exceeds the number of particles!");
       printf("%s\n","Exiting.");
       exit(EXIT_FAILURE);
     }

     //Format of positions is {mu, gamma, d, phi, pin}
     // TODO : limits[4]=round does not work on cluster (differing compiler version
     // need a fix that works more robustly
     int limits[dimensions];
     limits[0]=mumax-mumin+1;
     limits[1]=gammamax-gammamin+1;
     limits[2]=dmax-dmin+1;
     limits[3]=0;
     limits[4]=round((pmax-pmin)/pinc)+1;
     gsl_matrix *sw_current_posn;
     sw_current_posn=gsl_matrix_calloc(particles,dimensions);
     gsl_matrix *sw_tprobs;
     sw_tprobs=gsl_matrix_calloc(particles,2*dimensions);
     gsl_matrix *sw_amps;
     sw_amps=gsl_matrix_calloc(particles,2*dimensions);
     for(i=0;i<particles;i++)
       for(j=0;j<2*dimensions;j++)
           gsl_matrix_set(sw_amps,i,j,gsl_rng_uniform(rgen));
     gsl_matrix *sw_bestpos;
     sw_bestpos=gsl_matrix_calloc(particles,dimensions);
     printf("%s %s","Running at ",stamp);
     printf("%s %i %s\n","Using ",num_procs," threads.");
     
     //Initialise the swarm:
     update_posn(sw_current_posn,sw_tprobs,limits,dmin,particles,1);
     //Initialise capacity
     update_capacity(sw_current_posn,tsteps,particles,capacity);
     //Calculate fitness
     update_fitness(fitness,sw_current_posn,particles,mumin,gammamin,dmin,pmin,pinc,capacity,xi,yi,delta,prec, send_data_tag,
                    return_data_tag,debug);
     //Set bestpos = currentpos

     for(i=0;i<particles;i++)
     {
         for(j=0;j<dimensions;j++)
         {
          gsl_matrix_set(sw_bestpos,i,j,gsl_matrix_get(sw_current_posn,i,j));
         }
     }
     //Set bestfitness = fitness
     for(i=0;i<particles;i++)
     {
         best_fitness[i]=fitness[i];
     }
     //Set global best fitness and position
     global_best_fitness=update_global_pos(global_best_fitness,particles,dimensions,sw_current_posn,global_best_position,global_best_fitness,fitness);
     //Initialise transition amplitudes
     update_amps(particles,dimensions,sw_amps,sw_current_posn,sw_bestpos,limits,global_best_position,omega,rp,rg,stationary_prob);
     //Nomalise to produce transition probabilities
     amps_to_probs(particles,sw_tprobs,sw_amps,stationary_prob);
     //Save positions to output file
     for(i=0;i<particles;i++)
     {
         for(dim=0;dim<dimensions;dim++)
         {
           fprintf(positions,"%i ",(int)gsl_matrix_get(sw_current_posn,i,dim));
         }
         fprintf(positions,"\n");
     }
     fflush(positions);
     //Save bestposition to output file
     for(dim=0;dim<dimensions;dim++)
     {
         fprintf(bestpositions,"%i ",global_best_position[dim]);
     }
     fprintf(bestpositions,"\n");
     fflush(bestpositions);
     //Save fitness to output file
     for(i=0;i<particles;i++)
     {
         fprintf(fitfile,"%f ",fitness[i]);
     }
     fprintf(fitfile,"\n");
     fflush(fitfile);

     //Iterate:
     if(debug==0)
     {
        printf("%s","Iteration: ");
        fflush(stdout);
     }
     
     for(iteration=0;iteration<iterations;iteration++)
     {
     particles_to_compute=1;
     if(debug==0)
     {
        printf("\r%s %i","Iteration: ",iteration+1);
        fflush(stdout);
     }
     else
         printf("%s %i\n","Iteration: ",iteration+1);
     //Update current position
     update_posn(sw_current_posn,sw_tprobs,limits,dmin,particles,0);
     //Update capacity
     update_capacity(sw_current_posn,tsteps,particles,capacity);
     //Update fitness
     update_fitness(fitness,sw_current_posn,particles,mumin,gammamin,dmin,pmin,pinc,capacity,xi,yi,delta,prec, send_data_tag,
                    return_data_tag,debug);
     //Update bestfitness and bestposition
     for(i=0;i<particles;i++)
     {
         if(fitness[i]>best_fitness[i])
         {
           best_fitness[i]=fitness[i];
           for(j=0;j<dimensions;j++)
           {
            gsl_matrix_set(sw_bestpos,i,j,gsl_matrix_get(sw_current_posn,i,j));
           }
         }
     }
     //Update global best position and fitness
     global_best_fitness=update_global_pos(global_best_fitness,particles,dimensions,sw_current_posn,global_best_position,global_best_fitness,fitness);
     //Update transition amplitudes
     update_amps(particles,dimensions,sw_amps,sw_current_posn,sw_bestpos,limits,global_best_position,omega,rp,rg,stationary_prob);
     //Normalise to produce transition probabilities
     amps_to_probs(particles,sw_tprobs,sw_amps,stationary_prob);
     //Save positions to output file
     for(i=0;i<particles;i++)
     {
         for(dim=0;dim<dimensions;dim++)
         {
           fprintf(positions,"%i ",(int)gsl_matrix_get(sw_current_posn,i,dim));
         }
         fprintf(positions,"\n");
     }
     fflush(positions);
     //Save bestposition to output file
     for(dim=0;dim<dimensions;dim++)
     {
         fprintf(bestpositions,"%i ",global_best_position[dim]);
     }
     fprintf(bestpositions,"\n");
     fflush(bestpositions);
     //Save fitness to output file
     for(i=0;i<particles;i++)
     {
         fprintf(fitfile,"%f ",fitness[i]);
     }
     fprintf(fitfile,"\n");
     fflush(fitfile);
     }

     for(rank=1;rank<num_procs;rank++)
     {
       MPI_Send(&particles_to_compute,1,MPI_INT,rank,send_data_tag,MPI_COMM_WORLD);
     }
     
     printf("\n");
     fclose(positions);
     fclose(bestpositions);
     fclose(fitfile);

     time(&endtime);
     timeinfo=localtime(&endtime);
     strendtime=asctime(timeinfo);
     printf("%s %s","Finished at ",strendtime);
     printf("%s %i %s\n","Runtime: ",(int)endtime-(int)starttime," seconds.");
   }
   //CODE FOR SLAVE PROCESS TO EXECUTE
   else
   {

     particles_to_compute=1;

     while(particles_to_compute==1)
     {

         MPI_Recv(&particles_to_compute,1,MPI_INT,root_process,send_data_tag,MPI_COMM_WORLD,&status);
         if(debug==1)
         {
                 printf("%s %i %s %i\n","Rank: ",my_id," particles_to_compute =",particles_to_compute);
                 fflush(stdout);
         }
         if(particles_to_compute==0)
             break;
         
         MPI_Recv(&c,1,MPI_DOUBLE,root_process,send_data_tag,MPI_COMM_WORLD,&status);
         MPI_Recv(&m,1,MPI_INT,root_process,send_data_tag,MPI_COMM_WORLD,&status);
         MPI_Recv(&g,1,MPI_INT,root_process,send_data_tag,MPI_COMM_WORLD,&status);
         MPI_Recv(&d,1,MPI_INT,root_process,send_data_tag,MPI_COMM_WORLD,&status);
         MPI_Recv(&t,1,MPI_INT,root_process,send_data_tag,MPI_COMM_WORLD,&status);
         MPI_Recv(&p,1,MPI_DOUBLE,root_process,send_data_tag,MPI_COMM_WORLD,&status);
         MPI_Recv(&xi,1,MPI_DOUBLE,root_process,send_data_tag,MPI_COMM_WORLD,&status);
         MPI_Recv(&yi,1,MPI_DOUBLE,root_process,send_data_tag,MPI_COMM_WORLD,&status);
         MPI_Recv(&delta,1,MPI_DOUBLE,root_process,send_data_tag,MPI_COMM_WORLD,&status);
         MPI_Recv(&prec,1,MPI_INT,root_process,send_data_tag,MPI_COMM_WORLD,&status);
         if(debug==1)
         {
                 printf("%s %i %s\n","Rank: ",my_id," MPI_Recv");
                 fflush(stdout);
         }
         lfitness=network_eff(c,m,g,d,t,p,xi,yi,delta,prec);
         MPI_Send(&lfitness,1,MPI_DOUBLE,root_process,return_data_tag,MPI_COMM_WORLD);
         if(debug==1)
         {
                 printf("%s %i %s\n","Rank: ",my_id," MPI_Send");
                 fflush(stdout);
         }
       }
       
     }
   

 MPI_Finalize();
}
