/* Functions required by BPN
   Guy Billings UCL 2011                                */

#include "bpnpsofuncs.h"

//////////////////////////////////////////////////////////
/*Function to precompute all binomial coefficients     
   
 inpdf:
 mpz_bin_ui(bcnum,n,i);

 hypergeometricpdf:
 mpz_bin_ui(bc_beta_active_inputs,beta,ai_cnt);
 mpz_bin_ui(bc_sub,mu_sub_beta,di-ai_cnt);
 mpz_bin_ui(bc_mu_d,mu,di);

 uniquecodes:
 mpz_bin_ui(bcnum,n,i);

 scjoint:
 mpz_bin_ui(bc_gamma_sy,gamp,(unsigned long int)sy);

 1) mu,sx;
 2) gamma,sy;

 3) mu,d; 				   (always subset of 1)
 4) mu_sub_beta,di-ai_cnt; (for all mu, find all coeffs for 
                            less than or equal to d draws)
 5) beta,ai_cnt;		   (always a subset of 2)       
 
 Function finds binomial coefficient for all (N,m) between Nini, 
 Nend and mini mend. Wrties to mpfr variables.            */
void tabulate_bins_fr(mpfr_t *bcs, int Nini, int Nend, int mini, int mend, mpfr_prec_t prec)
{
 int N;
 int m;
 int done=mend-mini+1;
 mpz_t bc;
 mpz_init(bc);
 mpz_t Nz;
 mpz_init(Nz);

 for(N=Nini;N<=Nend;N++)
 {
  mpz_set_ui(Nz,N);
  for(m=mini;m<=mend;m++)
  {
   mpz_bin_ui(bc,Nz,m);
   mpfr_set_z(*(bcs+(N-Nini)*done+(m-mini)),bc,MPFR_RNDN);
  }
 }
 mpz_clear(bc);
 mpz_clear(Nz);
}   
//////////////////////////////////////////////////////////  
/* Function finds binomial coefficient for all (N,m) between Nini, 
   Nend and mini mend. Wrties to mpz variables.        */
void tabulate_bins_z(mpz_t bcs[], int Nini, int Nend, int mini, int mend)
{
 int N;
 int m;
 int done=mend-mini+1;
 mpz_t bc;
 mpz_init(bc);
 mpz_t Nz;
 mpz_init(Nz);

 for(N=Nini;N<=Nend;N++)
 {
  mpz_set_ui(Nz,N);
  for(m=mini;m<=mend;m++)
  {
   if(N<1)
   {
    mpz_set_ui(bcs[(N-Nini)*done+(m-mini)],0);
   }
   else
   {
    mpz_bin_ui(bc,Nz,m);
    mpz_set(bcs[(N-Nini)*done+(m-mini)],bc);
   } 
  }
 }
 mpz_clear(bc);
 mpz_clear(Nz);
}      
//////////////////////////////////////////////////////////
/* Function to generate the binomial pdf                */
void inpdf(mpfr_t *out,mpz_t n,int psize,double pini,double pinc,unsigned long mu,mpfr_t bcs[],mpfr_prec_t prec)
{
 mpfr_t s;
 mpfr_init2(s,prec);
 mpfr_t f;
 mpfr_init2(f,prec);
 mpfr_t unity;
 mpfr_init2(unity,prec);
 mpfr_set_ui(unity,(unsigned long int) 1,MPFR_RNDN);
 mpfr_t Z;
 mpfr_init2(Z,prec);
 mpfr_set_d(Z,0,MPFR_RNDN);
 mpfr_t pfr;
 mpfr_init2(pfr,prec);
 int pdex; 
 int odex;
 int odn;
 int i;
  
 for(pdex=0;pdex<=psize-1;pdex=pdex++)  
 {
  mpfr_set_d(pfr,pini+pdex*pinc,MPFR_RNDN);
  for(i=0;i<=mu;i++)
  {  
     odex=pdex*(mu+1)+i;
     mpfr_pow_ui(s, pfr,(unsigned long int)i,MPFR_RNDN);
     mpfr_sub(f, unity, pfr,MPFR_RNDN);
     mpfr_pow_ui(f,f,(unsigned long int)mu-i,MPFR_RNDN);
     mpfr_mul(*(out+odex),s,f,MPFR_RNDN);
     mpfr_mul(*(out+odex),bcs[i],*(out+odex),MPFR_RNDN);
     mpfr_add(Z,Z,*(out+odex),MPFR_RNDN);  
   }
   for(i=0;i<=mu;i++)
     odn=pdex*(mu+1)+i;  
     mpfr_div(*(out+odn), *(out+odn),Z,MPFR_RNDN); 
 }    
 
mpfr_clear(pfr); 
mpfr_clear(s);
mpfr_clear(f);
mpfr_clear(unity);
mpfr_clear(Z);
 
}
//////////////////////////////////////////////////////////
/* Function to generate the hypergeometric pdf          */
void hypergeometricpdf(mpfr_t *out,mpz_t mu, mpz_t d,unsigned int mi,unsigned int conns,unsigned int phi, mpz_t bc_b_ai[], mpz_t mu_bc[])
{
 unsigned int dini=(unsigned int)conns;
 unsigned int dmax=(unsigned int)conns;
 
 mpz_t bcm_b_ai;
 mpz_init(bcm_b_ai);
 mpq_t quotient;
 mpq_init(quotient);
 unsigned int odex=0;
 unsigned int di;
// unsigned int phi;
 unsigned int bi;
 unsigned int ai_cnt;
 //loop over d, phi and beta 
 for(di=dini;di<=dmax;di=di++)
 {
  //for(phi=1;phi<=di;phi++)
  //{
   for(bi=0;bi<=mi;bi++)
   {
    for(ai_cnt=phi;ai_cnt<=di;ai_cnt++)
    { 
      if((mi-bi)<=0)           //Then the number of successful picks is the number of draws with probability 1
      {
        if(ai_cnt==di)              
         mpq_set_d(quotient,1);
        else
         mpq_set_d(quotient,0);
      }
      else
      {
       mpz_mul(bcm_b_ai,bc_b_ai[bi*(dmax+1)+(ai_cnt)],bc_b_ai[(mi-bi)*(dmax+1)+(di-ai_cnt)]);
       mpq_set_num(quotient, bcm_b_ai);
       mpq_set_den(quotient, mu_bc[di-1]);
       mpq_canonicalize(quotient);
      }
      mpfr_set_q(*(out+odex),quotient,MPFR_RNDN);
      odex++; //Address array with simple counter increment;
    }
   }
  //}
 }
 mpz_clear(bcm_b_ai);
 mpq_clear(quotient);
} 
//////////////////////////////////////////////////////////
/* Function to calculate the probability that a cell is on
   given some number of active inputs and the threshold */
void pphi(mpfr_t *pout,mpz_t mu,mpz_t dmp,int m,int d, int phi,mpz_t bc_b_ai[], mpz_t mu_bc[], mpfr_prec_t prec)
{
 int dini=d;
 int dmax=d;

 //Required size of hp array is (1/6)*mi*(dmax*(dmax + 1)*(dmax + 2) - (dini - 1)*((dini - 1) + 1)*((dini - 1) + 2))
 //int hpsize=((m+1)*(dmax*(dmax+1)*(dmax+2)-(dini-1)*((dini-1)+1)*((dini-1)+2)))/6;
 int hpsize=(m+1)*(d-phi+1);
// unsigned int index1=0;
 unsigned int index2=0;
 unsigned int di;
 //unsigned int phi;
 unsigned int bi;
 unsigned int ai;
 unsigned int i;

 mpfr_t *hp;
 hp=(mpfr_t *) malloc(hpsize*sizeof(mpfr_t)); 
 for(i=0;i<=hpsize-1;i++)
  mpfr_init2(*(hp+i),prec);
 
 hypergeometricpdf(hp,mu,dmp,(unsigned long int)m,(unsigned long int)d,(unsigned long int)phi,bc_b_ai,mu_bc);

 for(di=dini;di<=dmax;di++)
 {
//  for(phi=1;phi<=di;phi++)
//  {
   for(bi=0;bi<=m;bi++)
   {
    for(ai=0;ai<=di-phi;ai++)
    {
     //mpfr_set(*(pout+index2),*(hp+index2),MPFR_RNDN);
     mpfr_add(*(pout+bi),*(pout+bi),*(hp+index2),MPFR_RNDN);
     index2++;
    }
//    index1++;
   }  
//  }   
 }
 for(i=0;i<=hpsize-1;i++)
   mpfr_clear(*(hp+i));

}
//////////////////////////////////////////////////////////
/* Function to calculate the expected number of unique codes 
   within each input success class                      */
void uniquecodes(mpfr_t *ucodes, int psiz, mpz_t n, unsigned long int mu, mpz_t ncodes ,mpfr_t *pdf,mpfr_prec_t prec)
{
//200 bits precision gives log[2^200] around 60 digits decimal precision
 mpfr_t unity;
 mpfr_init2(unity,prec);
 mpfr_set_ui(unity,(unsigned long int) 1,MPFR_RNDN);
 mpfr_t nunity;
 mpfr_init2(nunity,prec);
 mpfr_set_si(nunity,(signed long int) -1,MPFR_RNDN);
 mpz_t bcnum;
 mpz_init(bcnum);
 mpfr_t bcnumf;
 mpfr_init2(bcnumf,prec);
 mpq_t cfrac;
 mpq_init(cfrac);
 mpfr_t cfracf;
 mpfr_init2(cfracf,prec);
 mpfr_t expargr;
 mpfr_init2(expargr,prec); 			
 mpfr_t r1;
 mpfr_init2(r1,prec);
 mpfr_t r2;
 mpfr_init2(r2,prec);
// mpfr_t pdf[mu];
 unsigned long int i;
 int pdex;
 int addr;
// for(i=0;i<=mu;i++)
//  mpfr_init2(pdf[i],prec);
 
// inpdf(pdf,n,pin,mu,bcs,prec);

for(pdex=0;pdex<=psiz-1;pdex++)
{
 for(i=0;i<=mu;i++)
 {
  addr=pdex*(mu+1)+i;
  mpz_bin_ui(bcnum,n,i);
  mpq_set_num(cfrac,ncodes);
  mpq_set_den(cfrac,bcnum);
  mpq_canonicalize(cfrac);
  mpfr_set_q(cfracf,cfrac,MPFR_RNDN);
  mpfr_mul(expargr,cfracf,*(pdf+addr),MPFR_RNDN);   
  mpfr_mul(expargr,expargr,nunity,MPFR_RNDN);
  mpfr_exp(r1,expargr,MPFR_RNDN);      
  mpfr_sub(r1,unity,r1,MPFR_RNDN);
  mpfr_set_z(bcnumf,bcnum,MPFR_RNDN); 			
  mpfr_mul(r2,bcnumf,r1,MPFR_RNDN);
  mpfr_set((*ucodes+addr),r2,MPFR_RNDN);          
  mpfr_round((*ucodes+addr),(*ucodes+addr));
 } 
} 
 mpfr_clear(unity);
 mpfr_clear(nunity);
 mpz_clear(bcnum);
 mpfr_clear(bcnumf);
 mpq_clear(cfrac);
 mpfr_clear(cfracf);
 mpfr_clear(expargr);
 mpfr_clear(r1);
 mpfr_clear(r2);
// for(i=0;i<=mu;i++)
//  mpfr_clear(pdf[i]);
} 
////////////////////////////////////////////////////////// 
/* Function to calculate the probability of joint 
   occurance of input class sx and output class sy      */
void scjoint(mpfr_t *distptr, int mu, int gamma, int psiz, int conns, int phi, mpfr_t *bcs,mpfr_t *binpdf,mpfr_t *pp,mpfr_prec_t prec)
{

 int dini=conns;
 int dmax=conns;

 mpfr_t poff;
 mpfr_init2(poff,prec);
 mpfr_t bcf;
 mpfr_init2(bcf,prec);
 mpfr_t f1;
 mpfr_init2(f1,prec);
 mpfr_t f2;
 mpfr_init2(f2,prec);
 mpfr_t f3;
 mpfr_init2(f3,prec);
 unsigned int i;
 unsigned int d;
// unsigned int phi; 
 unsigned int sy;
 unsigned int sx;
 unsigned int pdex;
 int done=gamma+1;
 int dist_dex=0;
 int ppdex=0;
 int bdex;
 
 for(pdex=0;pdex<=psiz-1;pdex++)
 {
  for(d=dini;d<=dmax;d++)
  {
//   for(phi=1;phi<=d;phi++)
//   {   
    for(sx=0;sx<=mu;sx++)
    {
     mpfr_ui_sub(poff,(unsigned long int)1,*(pp+ppdex),MPFR_RNDN); 
     bdex=pdex*(mu+1)+sx;
     for(sy=0;sy<=gamma;sy++)
     {
      mpfr_mul(f1,*(binpdf+bdex),*(bcs+sy),MPFR_RNDN); 
      mpfr_pow_ui(f2,*(pp+ppdex),(unsigned long int)sy,MPFR_RNDN);
      mpfr_pow_ui(f3,poff,(unsigned long int)(gamma-sy),MPFR_RNDN);
      mpfr_mul(f1,f1,f2,MPFR_RNDN);
      mpfr_mul(*(distptr+dist_dex),f1,f3,MPFR_RNDN);
      dist_dex++;
     }
     ppdex++;
    }    
//   }
  }
  ppdex=0;
 } 
 
 mpfr_clear(poff);
 mpfr_clear(bcf);
 mpfr_clear(f1);
 mpfr_clear(f2);
 mpfr_clear(f3);
 
}
//////////////////////////////////////////////////////////
/* Function to compute the factorial                    */
unsigned long long fact(unsigned long long x)
{
 register int i;
 unsigned long long fact=1;
 for(i=1;i<=x;i++)
   fact=i*fact;
 return fact;  
}
//////////////////////////////////////////////////////////
/* Function to compute binomial pdf for output cells    */
void binomialpdf(double *bcs, int n,double p)
{
 int i;
 double tot=0;
 for(i=0;i<=n;i++)
 {
   *(bcs+i)=1;//;(fact(n)/(fact(n-i)*fact(i)))*pow(p,i)*pow(1-p,n-i);
   tot=tot+*(bcs+i);
 } 
 for(i=0;i<=n;i++)
 {
   *(bcs+i)=*(bcs+i)/tot;
 } 
}
////////////////////////////////////////////////////////// 
/* Function to calculate the entropy	                */
void eff_point(mpfr_t *ent, mpfr_t *ipdf, mpfr_t *jdist, mpfr_t *pp, mpfr_t *gamma_bcs, mpfr_t *p0, int conns, int phi, int mu, int gamma, double pin, double xi, double yi, double delta, mpfr_prec_t prec)
{
 unsigned int z,sx,sy,i;
 int ppdex=0,p0dex=0,thr;
 double po,cost;
 double *bpdf;
 bpdf=(double *) malloc((conns+1)*sizeof(double));
 binomialpdf(bpdf,conns,pin); //CONSIDER USING GSL version here?

 mpfr_t rp;
 mpfr_init2(rp,prec);
 mpfr_set_d(rp,0,MPFR_RNDN);
 mpfr_t nunity;
 mpfr_init2(nunity,prec);
 mpfr_set_d(nunity,-1,MPFR_RNDN);
 mpfr_t tbin;
 mpfr_init2(tbin,prec);
 mpfr_t t1;
 mpfr_init2(t1,prec);
 mpfr_t gate;
 mpfr_init2(gate,prec);
 mpfr_t qa;
 mpfr_init2(qa,prec);
 mpfr_t s0;
 mpfr_init2(s0,prec);
 mpfr_t jdist_tot;
 mpfr_init2(jdist_tot,prec);
 
 for(thr=phi;thr<=conns;thr++)
   po=po+*(bpdf+thr);
 cost=mu*(pin*xi+(1-pin))+delta*gamma*(po*yi+(1-po));

 for(sx=0;sx<=mu;sx++)
 {
   ppdex=sx;
   p0dex=sx;
        for(sy=1;sy<=gamma;sy++)
        {
           mpfr_set_d(qa,0,MPFR_RNDN);
           for(z=0;z<=gamma-sy;z++)
           {      
              mpfr_pow_ui(rp,nunity,(unsigned long int)z,MPFR_RNDN);
              mpfr_mul(rp,rp,*(gamma_bcs+(gamma-sy)*(gamma+1)+z),MPFR_RNDN);
       
                  if(mpfr_cmp_d(*(pp+ppdex),0)>0)
                  {
        
                   if(mpfr_cmp(*(p0+p0dex),*(pp+ppdex))>0)
                    mpfr_set_d(gate,1,MPFR_RNDN);
                   else
                   {
                     mpfr_sub(gate,*(pp+ppdex),*(p0+p0dex),MPFR_RNDN);
                     mpfr_div(gate,gate,*(pp+ppdex),MPFR_RNDN);
                   }
                   mpfr_pow_ui(t1,*(pp+ppdex),(unsigned long int)z+sy,MPFR_RNDN);
                   mpfr_mul(rp,rp,gate,MPFR_RNDN); 
                  }  
                  else
                    mpfr_set_d(t1,0,MPFR_RNDN);
        
                  mpfr_mul(rp,rp,t1,MPFR_RNDN);
                  mpfr_add(qa,qa,rp,MPFR_RNDN);
           }
          
           mpfr_add(qa,qa,*(p0+p0dex),MPFR_RNDN);
           mpfr_set_d(jdist_tot,0,MPFR_RNDN);
           for(i=0;i<=mu;i++)
           { 
             if(i!=sx)
               mpfr_add(jdist_tot, jdist_tot,*(jdist+i*(gamma+1)+sy),MPFR_RNDN);
           }
      
           mpfr_div(jdist_tot,jdist_tot,*(gamma_bcs+gamma*(gamma+1)+sy),MPFR_RNDN);
           mpfr_mul(qa,qa,*(ipdf+sx),MPFR_RNDN);
           mpfr_add(qa,qa,jdist_tot,MPFR_RNDN);
      
           if(mpfr_cmp_d(qa,0)>0)
           {
             mpfr_log2(qa,qa,MPFR_RNDN);
             mpfr_mul(qa,qa,*(jdist+sx*(gamma+1)+sy),MPFR_RNDN);
           }
           else
             mpfr_set_d(qa,0,MPFR_RNDN);
           mpfr_sub(*ent,*ent,qa,MPFR_RNDN);

         }
 }
 mpfr_set_d(jdist_tot,0,MPFR_RNDN);  
 for(i=0;i<=mu;i++) 
   mpfr_add(jdist_tot, jdist_tot,*(jdist+i*(gamma+1)),MPFR_RNDN);
 if(mpfr_cmp_d(jdist_tot,0)>0)  
 {
    mpfr_log2(s0,jdist_tot,MPFR_RNDN);  
    mpfr_mul(s0,jdist_tot,s0,MPFR_RNDN);  
    mpfr_sub(*ent,*ent,s0,MPFR_RNDN);  
    mpfr_div_d(*ent,*ent,cost,MPFR_RNDN);
 }  

 mpfr_clear(s0);
 mpfr_clear(rp);
 mpfr_clear(nunity);
 mpfr_clear(tbin);
 mpfr_clear(t1);
 free(bpdf);
 
}
//////////////////////////////////////////////////////////
/*Standalone function to calculate the efficiency of a 
  network                                               */
// TODO : Refactor function so that memory allocation and deallocation
//calls are minimized.
double network_eff(double bfrac, int mu, int gamma, int d, int phi, double pin,double xi, double yi, double delta, int prec)
{
 int psiz=1,i,dmax=d;
 double pini=pin,pmax=pin, pinc=0.1,effout;

 mpz_t mump;
 mpz_init_set_ui(mump,(unsigned long int)mu);
 mpz_t dmp;
 mpz_init_set_ui(dmp,(unsigned long int)dmax);
 mpz_t gamp;
 mpz_init(gamp);
 mpz_set_ui(gamp,gamma);
 mpz_t ncodes;
 mpz_init(ncodes);
 mpfr_t ncodesfr;
 mpfr_t nc_pow;
 mpfr_init2(nc_pow,prec);
 mpfr_init2(ncodesfr,prec);
 mpfr_set_d(ncodesfr,2,MPFR_RNDN);
 mpfr_set_d(nc_pow,bfrac*mu,MPFR_RNDN);
 mpfr_pow(ncodesfr,ncodesfr,nc_pow,MPFR_RNDN);
 mpfr_get_z(ncodes,ncodesfr,MPFR_RNDN);

 int pdfsize=(mu+1)*psiz;
 mpfr_t * ipdfptr;
 mpfr_t * ucodesptr;
 mpfr_t * p0ptr;
 ipdfptr=(mpfr_t *) malloc(pdfsize*sizeof(mpfr_t));
 ucodesptr=(mpfr_t *) malloc(pdfsize*sizeof(mpfr_t));
 p0ptr=(mpfr_t *) malloc(pdfsize*sizeof(mpfr_t));
 for(i=0;i<=pdfsize-1;i++)
 {
  mpfr_init2(*(ipdfptr+i),prec);
  mpfr_set_d(*(ipdfptr+i),0,MPFR_RNDN);
  mpfr_init2(*(ucodesptr+i),prec);
  mpfr_init2(*(p0ptr+i),prec);
 }

 int rsize=(mu+1)*(gamma+1);
 mpfr_t * distptr;
 distptr=(mpfr_t *) malloc(rsize*sizeof(mpfr_t));
 for(i=0;i<=rsize-1;i++)
 {
   mpfr_init2(*(distptr+i),prec);
 }

 mpfr_t *eff;
 int esize=1;
 eff=(mpfr_t *) malloc(esize*sizeof(mpfr_t));
 for(i=0;i<=esize-1;i++)
 {
   mpfr_init2(*(eff+i),prec);
   mpfr_set_d(*(eff+i),0,MPFR_RNDN);
 }

 //Binomial coeffs: (z,x) with z:0->gamma, x:0->gamma
 mpfr_t * gamma_bcfr;
 int gbf_size=(gamma+1)*(gamma+1);
 gamma_bcfr=(mpfr_t *)malloc(gbf_size*sizeof(mpfr_t));
 for(i=0;i<=gbf_size-1;i++)
   mpfr_init2(*(gamma_bcfr+i),prec);

 //Binomial coeffs: (y,x) with y:0->mu, x:0->d
 int bacsize=(mu+1)*(dmax+1);
 mpz_t beta_ai_cnt[bacsize];
 for(i=0;i<=bacsize-1;i++)
   mpz_init(beta_ai_cnt[i]);

 //Binomial coeffs: (mu,x) with x:0->mu
 mpfr_t mu_bcfr[mu+1];
 for(i=0;i<=mu;i++)
   mpfr_init2(mu_bcfr[i],prec);

 //Binomial coeff: (mu,d) x:1->d (integer)
 mpz_t mu_bc[dmax];
 for(i=0;i<=dmax-1;i++)
   mpz_init(mu_bc[i]);

 int ppsize=(mu+1);//0.5*(dmax*(dmax+1)-(dini-1)*dini)*(mu+1);
 mpfr_t * ppptr;
 ppptr=(mpfr_t *) malloc(ppsize*sizeof(mpfr_t));
 for(i=0;i<=ppsize-1;i++)
 {
   mpfr_init2(*(ppptr+i),prec);
   mpfr_set_d(*(ppptr+i),0,MPFR_RNDN);
 }

 tabulate_bins_fr(gamma_bcfr,0,gamma,0,gamma,prec);
 tabulate_bins_fr(mu_bcfr,mu,mu,0,mu,prec);
 tabulate_bins_z(mu_bc,mu,mu,1,dmax);
 tabulate_bins_z(beta_ai_cnt,0,mu,0,dmax);
 inpdf(ipdfptr,mump,psiz,pini,pinc,mu,mu_bcfr,prec);
 uniquecodes(ucodesptr,psiz,mump,mu,ncodes,ipdfptr,prec);
 pphi(ppptr,mump,dmp,mu,d,phi,beta_ai_cnt,mu_bc,prec);
 scjoint(distptr,mu,gamma,psiz,d,phi,(gamma_bcfr+gamma*(gamma+1)),ipdfptr,ppptr,prec);

 for(i=0;i<=pdfsize-1;i++)
 {
    if(mpfr_cmp_si(*(ucodesptr+i),0)<=0)
      mpfr_set_ui(*(p0ptr+i),0,MPFR_RNDN);
    else
      mpfr_ui_div(*(p0ptr+i),1,*(ucodesptr+i),MPFR_RNDN);
 }

 eff_point(eff,ipdfptr,distptr,ppptr,gamma_bcfr,p0ptr,d,phi,mu,gamma,pin,xi,yi,delta,prec);
 effout=mpfr_get_d(*eff,MPFR_RNDN);

 mpz_clear(mump);
 mpz_clear(dmp);
 mpz_clear(gamp);
 mpz_clear(ncodes);
 mpfr_clear(ncodesfr);
 mpfr_clear(nc_pow);

 for(i=0;i<=pdfsize-1;i++)
 {
  mpfr_clear(*(ipdfptr+i));
  mpfr_clear(*(ucodesptr+i));
  mpfr_clear(*(p0ptr+i));
 }
 free(ipdfptr);
 free(ucodesptr);
 free(p0ptr);

 for(i=0;i<=rsize-1;i++)
   mpfr_clear(*(distptr+i));
 free(distptr);

 for(i=0;i<=esize-1;i++)
   mpfr_clear(*(eff+i));
 free(eff);

 for(i=0;i<=(gamma+1)*(gamma+1)-1;i++)
   mpfr_clear(*(gamma_bcfr+i));
 free(gamma_bcfr);

 for(i=0;i<=bacsize-1;i++)
   mpz_clear(beta_ai_cnt[i]);

 for(i=0;i<=mu;i++)
   mpfr_clear(mu_bcfr[i]);

 for(i=0;i<=dmax-1;i++)
   mpz_clear(mu_bc[i]);

 for(i=0;i<=ppsize-1;i++)
   mpfr_clear(*(ppptr+i));
 free(ppptr);

 return effout;
}

